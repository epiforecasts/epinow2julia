#' Conversion layer between R and Julia
#'
#' @description Internal functions for converting R objects (dist_spec,
#'   *_opts) to their Julia equivalents and converting Julia results back to
#'   R data.tables.
#'
#' @keywords internal

# ── Distribution conversion ──────────────────────────────────────────────

#' Convert an R dist_spec to a Julia distribution
#'
#' @param dist A `dist_spec` object
#' @return A JuliaCall proxy object representing the Julia distribution
#' @keywords internal
dist_spec_to_julia <- function(dist) {
  if (inherits(dist, "multi_dist_spec")) {
    # Composite distribution: convert each component and combine
    components <- lapply(seq_len(ndist(dist)), function(i) {
      single <- extract_single_dist(dist, i)
      dist_spec_to_julia(single)
    })
    return(Reduce(function(a, b) {
      juliaCall("+", a, b)
    }, components))
  }

  family <- get_distribution(dist)

  if (family == "nonparametric") {
    return(juliaCall("NonParametricDist", get_pmf(dist)))
  }

  if (family == "fixed") {
    value <- dist$parameters$value
    return(juliaCall("Dirac", as.numeric(value)))
  }

  # Check for uncertain parameters (dist_spec nested inside parameters)
  params <- dist$parameters
  uncertain <- vapply(params, function(p) inherits(p, "dist_spec"), logical(1))

  if (any(uncertain)) {
    return(dist_spec_uncertain_to_julia(dist, family, params, uncertain))
  }

  # Fixed parametric distribution
  dist_max <- max(dist)
  jdist <- switch(family,
    lognormal = juliaCall(
      "LogNormal", as.numeric(params$meanlog), as.numeric(params$sdlog)
    ),
    gamma = juliaCall(
      "Gamma", as.numeric(params$shape), 1.0 / as.numeric(params$rate)
    ),
    normal = juliaCall(
      "Normal", as.numeric(params$mean), as.numeric(params$sd)
    )
  )

  if (is.finite(dist_max)) {
    jdist <- juliaCall(
      "truncated", jdist,
      lower = 0.0, upper = as.numeric(dist_max)
    )
  }
  jdist
}

#' Convert an uncertain dist_spec to a Julia UncertainDistribution
#'
#' @param dist The dist_spec object
#' @param family Character, the distribution family
#' @param params Named list of parameters
#' @param uncertain Logical vector indicating which params are uncertain
#' @return A JuliaCall proxy object
#' @keywords internal
dist_spec_uncertain_to_julia <- function(dist, family, params, uncertain) {
  dist_max <- max(dist)
  if (!is.finite(dist_max)) dist_max <- 20.0

  # Build Julia Normal priors for each parameter
  param_priors <- lapply(params, function(p) {
    if (inherits(p, "dist_spec")) {
      m <- mean(p, ignore_uncertainty = TRUE)
      s <- sd(p, ignore_uncertainty = TRUE)
      juliaCall("Normal", as.numeric(m), as.numeric(s))
    } else {
      juliaCall("Normal", as.numeric(p), 0.001)
    }
  })

  # Build the constructor using the natural parameterisation.
  # The priors are truncated at lower bounds so values stay valid.
  # For Gamma: R uses (shape, rate) but Julia uses (shape, scale).
  # We pass scale directly to avoid 1/rate singularity during AD.
  constructor_str <- switch(family,
    lognormal = "(mu, sigma) -> LogNormal(mu, sigma)",
    gamma = "(alpha, scale) -> Gamma(alpha, scale)",
    normal = "(mu, sigma) -> Normal(mu, sigma)"
  )

  constructor <- juliaEval(constructor_str)

  # Build priors vector as Julia expression.
  # Use truncated Normal for params that must be positive.
  # For Gamma: convert rate prior to scale prior (1/rate) via delta method
  # to avoid the 1/x singularity during AD.
  lb <- lower_bounds(family)
  param_names <- natural_params(family)
  prior_strs <- vapply(seq_along(params), function(j) {
    p <- params[[j]]
    pname <- param_names[j]
    bound <- lb[pname]
    is_gamma_rate <- (family == "gamma" && pname == "rate")
    if (inherits(p, "dist_spec")) {
      m <- mean(p, ignore_uncertainty = TRUE)
      s <- sd(p, ignore_uncertainty = TRUE)
      if (is_gamma_rate) {
        # Convert rate prior to scale prior via delta method
        scale_m <- 1.0 / as.numeric(m)
        scale_s <- as.numeric(s) / as.numeric(m)^2
        sprintf(
          "truncated(Normal(%s, %s); lower=0.01)",
          scale_m, scale_s
        )
      } else if (!is.na(bound) && is.finite(bound) && bound >= 0) {
        sprintf(
          "truncated(Normal(%s, %s); lower=%s)",
          as.numeric(m), as.numeric(s), as.numeric(bound)
        )
      } else {
        sprintf("Normal(%s, %s)", as.numeric(m), as.numeric(s))
      }
    } else {
      if (is_gamma_rate) {
        sprintf("Normal(%s, 0.001)", 1.0 / as.numeric(p))
      } else {
        sprintf("Normal(%s, 0.001)", as.numeric(p))
      }
    }
  }, character(1))

  juliaEval(sprintf(
    'UncertainDistribution(%s, Distribution[%s], %.1f)',
    constructor_str,
    paste(prior_strs, collapse = ", "),
    as.numeric(dist_max)
  ))
}

#' Convert a Julia distribution back to an R dist_spec
#'
#' @param jdist A JuliaCall proxy for a Julia distribution
#' @return A `dist_spec` object
#' @keywords internal
julia_dist_to_r <- function(jdist) {
  type <- juliaCall("string", juliaCall("typeof", jdist))

  if (grepl("NonParametricDist", type)) {
    pmf <- juliaCall("getfield", jdist, juliaEval(":pmf"))
    return(NonParametric(pmf))
  }

  if (grepl("Dirac", type)) {
    val <- juliaCall("params", jdist)
    return(Fixed(value = val[[1]]))
  }

  if (grepl("LogNormal", type)) {
    p <- juliaCall("params", jdist)
    return(LogNormal(meanlog = p[[1]], sdlog = p[[2]]))
  }

  if (grepl("Gamma", type)) {
    p <- juliaCall("params", jdist)
    return(Gamma(shape = p[[1]], rate = 1.0 / p[[2]]))
  }

  if (grepl("Normal", type)) {
    p <- juliaCall("params", jdist)
    return(Normal(mean = p[[1]], sd = p[[2]]))
  }

  # Fallback
  Fixed(value = 0)
}


# ── Opts conversion ──────────────────────────────────────────────────────

#' Convert R gt_opts to Julia GTOpts
#' @param opts A `generation_time_opts` object
#' @return JuliaCall proxy for Julia GTOpts
#' @keywords internal
r_gt_opts_to_julia <- function(opts) {
  # gt_opts IS a dist_spec with weight_prior as attribute
  jdist <- dist_spec_to_julia(opts)
  weight_prior <- isTRUE(attr(opts, "weight_prior"))
  juliaCall("gt_opts", jdist, weight_prior = weight_prior)
}

#' Convert R delay_opts to Julia DelayOpts
#' @param opts A `delay_opts` object
#' @return JuliaCall proxy for Julia DelayOpts
#' @keywords internal
r_delay_opts_to_julia <- function(opts) {
  # delay_opts IS a dist_spec with weight_prior as attribute
  jdist <- dist_spec_to_julia(opts)
  weight_prior <- isTRUE(attr(opts, "weight_prior"))
  juliaCall("delay_opts", jdist, weight_prior = weight_prior)
}

#' Convert R trunc_opts to Julia TruncOpts
#' @param opts A `trunc_opts` object
#' @return JuliaCall proxy for Julia TruncOpts
#' @keywords internal
r_trunc_opts_to_julia <- function(opts) {
  # trunc_opts IS a dist_spec with weight_prior as attribute
  jdist <- dist_spec_to_julia(opts)
  weight_prior <- isTRUE(attr(opts, "weight_prior"))
  juliaCall("trunc_opts", jdist, weight_prior = weight_prior)
}

#' Convert R rt_opts to Julia RtOpts
#' @param opts An `rt_opts` object or NULL
#' @return JuliaCall proxy for Julia RtOpts, or NULL
#' @keywords internal
r_rt_opts_to_julia <- function(opts) {
  if (is.null(opts) || !isTRUE(opts$use_rt)) {
    return(juliaEval("rt_opts(use_rt = false)"))
  }

  prior <- dist_spec_to_julia(opts$prior)

  # Map R strings to Julia enum names
  gp_on <- switch(opts$gp_on,
    "R_t-1" = "gp_Rt", "R0" = "gp_R0", "gp_Rt"
  )
  future_str <- switch(opts$future,
    "latest" = "latest", "project" = "project",
    "estimate" = "estimate", "latest"
  )
  pop_period <- switch(opts$pop_period,
    "forecast" = "pop_forecast", "all" = "pop_all", "pop_forecast"
  )

  pop_val <- if (inherits(opts$pop, "dist_spec") &&
                 get_distribution(opts$pop) != "fixed") {
    mean(opts$pop, ignore_uncertainty = TRUE)
  } else if (inherits(opts$pop, "dist_spec")) {
    as.numeric(opts$pop$parameters$value)
  } else {
    0.0
  }

  # Use juliaLet to pass the prior as a variable while using enums inline
  juliaEval(sprintf(
    'function _make_rt_opts(p) rt_opts(prior = p, use_rt = true, rw = %d, future = %s, gp_on = %s, pop = %.1f, pop_period = %s, pop_floor = %.1f) end',
    as.integer(opts$rw), future_str, gp_on,
    as.numeric(pop_val), pop_period,
    as.numeric(opts$pop_floor %||% 1.0)
  ))
  juliaCall("_make_rt_opts", prior)
}

#' Convert R gp_opts to Julia GPOpts
#' @param opts A `gp_opts` object or NULL
#' @return JuliaCall proxy for Julia GPOpts, or NULL
#' @keywords internal
r_gp_opts_to_julia <- function(opts) {
  if (is.null(opts) || opts$basis_prop == 0) {
    return(juliaEval("gp_opts(basis_prop = 0.0)"))
  }

  ls_dist <- dist_spec_to_julia(opts$ls)
  alpha_dist <- dist_spec_to_julia(opts$alpha)

  kernel <- switch(opts$kernel,
    "matern" = "matern", "se" = "se", "ou" = "matern",
    "periodic" = "periodic", "matern"
  )

  juliaEval(sprintf(
    'function _make_gp_opts(ls, alpha) gp_opts(basis_prop = %.4f, boundary_scale = %.4f, ls = ls, alpha = alpha, kernel = %s, matern_order = %.1f, w0 = %.4f) end',
    as.numeric(opts$basis_prop), as.numeric(opts$boundary_scale),
    kernel, as.numeric(opts$matern_order), as.numeric(opts$w0)
  ))
  juliaCall("_make_gp_opts", ls_dist, alpha_dist)
}

#' Convert R obs_opts to Julia ObsOpts
#' @param opts An `obs_opts` object
#' @return JuliaCall proxy for Julia ObsOpts
#' @keywords internal
r_obs_opts_to_julia <- function(opts) {
  family <- switch(opts$family,
    "negbin" = "negbin", "poisson" = "poisson", "negbin"
  )

  dispersion <- dist_spec_to_julia(opts$dispersion)

  # Fixed scale should be passed as Float64, not Dirac
  scale <- if (inherits(opts$scale, "dist_spec") &&
               get_distribution(opts$scale) != "fixed") {
    dist_spec_to_julia(opts$scale)
  } else if (inherits(opts$scale, "dist_spec")) {
    as.numeric(opts$scale$parameters$value)
  } else {
    as.numeric(opts$scale)
  }

  juliaEval(sprintf(
    'function _make_obs_opts(disp, sc) obs_opts(family = %s, dispersion = disp, weight = %.4f, week_effect = %s, week_length = %d, scale = sc, likelihood = %s) end',
    family,
    as.numeric(opts$weight),
    tolower(as.character(isTRUE(opts$week_effect))),
    as.integer(opts$week_length),
    tolower(as.character(isTRUE(opts$likelihood)))
  ))
  juliaCall("_make_obs_opts", dispersion, scale)
}

#' Convert R backcalc_opts to Julia BackcalcOpts
#' @param opts A `backcalc_opts` object
#' @return JuliaCall proxy for Julia BackcalcOpts
#' @keywords internal
r_backcalc_opts_to_julia <- function(opts) {
  prior <- switch(opts$prior,
    "infections" = "bc_infections", "none" = "bc_none",
    "reports" = "bc_infections", "growth_rate" = "bc_growth_rate",
    "bc_infections"
  )
  juliaEval(sprintf(
    'backcalc_opts(prior = %s, prior_window = %d, rt_window = %d)',
    prior, as.integer(opts$prior_window), as.integer(opts$rt_window)
  ))
}

#' Convert R forecast_opts to Julia ForecastOpts
#' @param opts A `forecast_opts` object
#' @return JuliaCall proxy for Julia ForecastOpts
#' @keywords internal
r_forecast_opts_to_julia <- function(opts) {
  juliaCall(
    "forecast_opts",
    horizon = as.integer(opts$horizon)
  )
}

#' Convert R secondary_opts to Julia SecondaryOpts
#' @param opts A `secondary_opts` object
#' @return JuliaCall proxy for Julia SecondaryOpts
#' @keywords internal
r_secondary_opts_to_julia <- function(opts) {
  # R secondary_opts stores flags, not a type field.
  # Infer type from the flag pattern.
  if (isTRUE(opts$cumulative) || opts$cumulative == 1) {
    juliaEval("secondary_opts(prevalence)")
  } else {
    juliaEval("secondary_opts(incidence)")
  }
}

#' Convert stan_opts to Julia InferenceOpts
#'
#' @description Maps the user-facing `stan_opts()` parameters to Julia's
#'   `inference_opts()`. The `stan_opts()` signature is preserved for backward
#'   compatibility.
#'
#' @param stan A `stan_opts` object
#' @return JuliaCall proxy for Julia InferenceOpts
#' @keywords internal
stan_opts_to_inference_opts <- function(stan) {
  method <- stan$method %||% "sampling"
  if (method != "sampling") {
    cli::cli_warn(
      c(
        "!" = "Method {.val {method}} is not yet supported by the Julia backend.",
        "i" = "Falling back to NUTS sampling."
      )
    )
  }

  # Extract sampling parameters
  samples <- as.integer(stan$args$iter_sampling %||%
    stan$args$samples %||% 2000L)
  warmup <- as.integer(stan$args$iter_warmup %||%
    stan$args$warmup %||% 250L)
  chains <- as.integer(stan$args$chains %||% 4L)
  seed <- stan$args$seed

  # Extract control parameters
  control <- stan$args$control %||% list()
  target_acceptance <- as.numeric(control$adapt_delta %||% 0.9)
  max_treedepth <- as.integer(control$max_treedepth %||% 12L)

  seed_str <- if (!is.null(seed)) sprintf(", seed = %d", as.integer(seed)) else ""
  juliaEval(sprintf(
    'inference_opts(samples = %d, warmup = %d, chains = %d, target_acceptance = %.4f, max_treedepth = %d%s)',
    samples, warmup, chains, target_acceptance, max_treedepth, seed_str
  ))
}


# ── Data conversion ──────────────────────────────────────────────────────

#' Convert R data.frame to Julia DataFrame
#'
#' @param data A data.frame with columns `date` and `confirm`
#' @return JuliaCall proxy for a Julia DataFrame
#' @keywords internal
r_data_to_julia <- function(data) {
  dates <- as.character(data$date)
  confirm <- as.integer(data$confirm)

  juliaEval(sprintf(
    'DataFrame(date = Date.([%s]), confirm = [%s])',
    paste(sprintf('"%s"', dates), collapse = ", "),
    paste(confirm, collapse = ", ")
  ))
}

#' Convert R secondary data.frame to Julia DataFrame
#'
#' @param data A data.frame with columns `date`, `primary`, `secondary`
#' @return JuliaCall proxy for a Julia DataFrame
#' @keywords internal
r_secondary_data_to_julia <- function(data) {
  dates <- as.character(data$date)
  primary <- as.integer(data$primary)
  secondary <- as.integer(data$secondary)

  juliaEval(sprintf(
    'DataFrame(date = Date.([%s]), primary = [%s], secondary = [%s])',
    paste(sprintf('"%s"', dates), collapse = ", "),
    paste(primary, collapse = ", "),
    paste(secondary, collapse = ", ")
  ))
}


# ── Result conversion ────────────────────────────────────────────────────

#' Convert Julia EstimateInfectionsResult samples to R data.table
#'
#' @param julia_result JuliaCall proxy for the Julia result
#' @param observations The original R observations data.frame
#' @param horizon Integer, the forecast horizon
#' @param seeding_time Integer, the seeding time
#' @return A `data.table` with columns: variable, time, date, sample, value,
#'   strat, type
#' @importFrom data.table data.table rbindlist setcolorder fcase
#' @keywords internal
julia_samples_to_r <- function(julia_result, observations, horizon = 0,
                               seeding_time = 0) {
  # Get samples from Julia - use a helper to convert dates/symbols to strings
  # to avoid JuliaConnectoR serialisation issues
  convert_fn <- juliaEval('function _convert_samples_df(result)
    df = get_samples(result)
    df2 = copy(df)
    if hasproperty(df2, :date)
      df2.date = string.(df2.date)
    end
    if hasproperty(df2, :variable)
      df2.variable = string.(df2.variable)
    end
    df2
  end')
  julia_df <- juliaCall("_convert_samples_df", julia_result)
  samples <- data.table::as.data.table(julia_df)

  # Ensure correct column types
  if ("date" %in% names(samples)) {
    samples[, date := as.Date(date)]
  }
  if ("variable" %in% names(samples)) {
    samples[, variable := as.character(variable)]
  }
  if ("sample" %in% names(samples)) {
    samples[, sample := as.integer(sample)]
  }
  if ("value" %in% names(samples)) {
    samples[, value := as.numeric(value)]
  }

  # Map Julia variable names to R names
  samples[variable == "R", variable := "R"]
  samples[variable == "infections", variable := "infections"]
  samples[variable == "reports", variable := "reported_cases"]

  # Remove log_R (not part of the standard R output)
  samples <- samples[variable != "log_R"]

  # Compute growth rate from infections: r(t) = log(I(t+1)/I(t))
  inf_samples <- samples[variable == "infections"]
  if (nrow(inf_samples) > 0) {
    data.table::setorder(inf_samples, sample, date)
    growth_rate <- inf_samples[
      ,
      .(
        date = date[-1],
        value = diff(log(pmax(value, 1e-6)))
      ),
      by = sample
    ]
    growth_rate[, variable := "growth_rate"]
    samples <- data.table::rbindlist(
      list(samples, growth_rate),
      fill = TRUE
    )
  }

  # Add time index per variable
  samples[, time := as.integer(date - min(date, na.rm = TRUE)) + 1L,
    by = variable
  ]

  # Add strat column
  samples[, strat := NA_character_]

  # Add type column
  max_obs_date <- max(observations$date, na.rm = TRUE)
  samples[, type := data.table::fcase(
    date > max_obs_date, "forecast",
    date > (max_obs_date - seeding_time), "estimate based on partial data",
    is.na(date), NA_character_,
    default = "estimate"
  )]

  data.table::setcolorder(
    samples,
    c("variable", "time", "date", "sample", "value", "strat", "type")
  )

  samples[]
}

#' Convert Julia summary DataFrame to R data.table
#'
#' @param julia_df JuliaCall proxy for a Julia DataFrame with summary stats
#' @return A `data.table` with columns: date, mean, median, sd, and CrI columns
#' @importFrom data.table as.data.table
#' @keywords internal
julia_summary_to_r <- function(julia_df) {
  dt <- data.table::as.data.table(julia_df)
  if ("date" %in% names(dt)) {
    dt[, date := as.Date(date)]
  }
  dt[]
}

#' Convert Julia parameters to R dist_spec list
#'
#' @param julia_result JuliaCall proxy for the Julia result
#' @return A named list of `dist_spec` objects
#' @keywords internal
julia_parameters_to_r <- function(julia_result) {
  julia_params <- juliaCall("get_parameters", julia_result)
  param_names <- juliaCall("collect",
    juliaCall("keys", julia_params)
  )

  result <- list()
  for (name in param_names) {
    name_str <- as.character(name)
    vals <- juliaCall("getindex", julia_params, name)
    # Convert posterior samples to Normal(mean, sd)
    m <- mean(as.numeric(vals))
    s <- stats::sd(as.numeric(vals))
    result[[name_str]] <- Normal(mean = round(m, 3), sd = round(s, 3))
  }
  result
}
