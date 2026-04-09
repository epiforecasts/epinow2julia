#' Fit an Integer Adjusted Exponential, Gamma or Lognormal distributions
#'
#' @description `r lifecycle::badge("stable")`
#' Fits an integer adjusted exponential, gamma or lognormal distribution using
#' stan.
#' @param values Numeric vector of values
#'
#' @param samples Numeric, number of samples to take. Must be >= 1000.
#' Defaults to 1000.
#'
#' @param dist Character string, which distribution to fit. Defaults to
#' exponential (`"exp"`) but gamma (`"gamma"`) and lognormal (`"lognormal"`) are
#' also supported.
#'
#' @param cores Numeric, defaults to 1. Number of CPU cores to use (no effect
#' if greater than the number of chains).
#'
#' @param chains Numeric, defaults to 2. Number of MCMC chains to use. More is
#' better with the minimum being two.
#'
#' @param verbose Logical, defaults to FALSE. Should verbose progress messages
#' be printed.
#'
#' @return A stan fit of an interval censored distribution
#' @export
#' @inheritParams stan_opts
#' @importFrom cli cli_warn col_blue
#' @examples
#' \donttest{
#' # integer adjusted exponential model
#' dist_fit(rexp(1:100, 2),
#'   samples = 1000, dist = "exp",
#'   cores = ifelse(interactive(), 4, 1), verbose = TRUE
#' )
#'
#'
#' # integer adjusted gamma model
#' dist_fit(rgamma(1:100, 5, 5),
#'   samples = 1000, dist = "gamma",
#'   cores = ifelse(interactive(), 4, 1), verbose = TRUE
#' )
#'
#' # integer adjusted lognormal model
#' dist_fit(rlnorm(1:100, log(5), 0.2),
#'   samples = 1000, dist = "lognormal",
#'   cores = ifelse(interactive(), 4, 1), verbose = TRUE
#' )
#' }
dist_fit <- function(values = NULL, samples = 1000, cores = 1,
                     chains = 2, dist = "exp", verbose = FALSE,
                     backend = "rstan") {
  ensure_julia()

  # Convert values to Julia DataFrame
  julia_data <- juliaEval(sprintf(
    'DataFrame(delay = [%s])',
    paste(as.numeric(values), collapse = ", ")
  ))

  family <- switch(dist,
    exp = juliaEval(":exponential"),
    lognormal = juliaEval(":lognormal"),
    gamma = juliaEval(":gamma")
  )

  # Call Julia's estimate_dist
  julia_result <- juliaCall(
    "estimate_dist",
    julia_data,
    family = family
  )

  julia_result
}


#' Fit a Subsampled Bootstrap to Integer Values and Summarise Distribution
#' Parameters
#'
#' @description `r lifecycle::badge("stable")`
#' Fits an integer adjusted distribution to a subsampled bootstrap of data and
#' then integrates the posterior samples into a single set of summary
#' statistics. Can be used to generate a robust reporting delay that accounts
#' for the fact the underlying delay likely varies over time or that the size
#' of the available reporting delay sample may not be representative of the
#' current case load.
#'
#' @param values Integer vector of values.
#'
#' @param dist Character string, which distribution to fit. Defaults to
#' lognormal (`"lognormal"`) but gamma (`"gamma"`) is also supported.
#'
#' @param verbose Logical, defaults to `FALSE`. Should progress messages be
#' printed.
#'
#' @param samples Numeric, number of samples to take overall from the
#' bootstrapped posteriors.
#'
#' @param bootstraps Numeric, defaults to 1. The number of bootstrap samples
#' (with replacement) of the delay distribution to take. If `samples` is less
#' than `bootstraps`, `samples` takes the value of `bootstraps`.
#'
#' @param bootstrap_samples Numeric, defaults to 250. The number of samples to
#' take in each bootstrap if the sample size of the supplied delay
#' distribution is less than its value.
#'
#' @param max_value Numeric, defaults to the maximum value in the observed
#' data. Maximum delay to  allow (added to output but does impact fitting).
#'
#' @return A `<dist_spec>` object summarising the bootstrapped distribution
#' @importFrom purrr list_transpose
#' @importFrom data.table data.table rbindlist
#' @importFrom cli cli_abort col_blue
#' @export
#' @examples
#' \donttest{
#' # lognormal
#' # bootstraps and samples have been reduced for this example
#' # for real analyses, use more
#' delays <- rlnorm(500, log(5), 1)
#' out <- bootstrapped_dist_fit(delays,
#'   samples = 500, bootstraps = 2,
#'   dist = "lognormal"
#' )
#' out
#' }
bootstrapped_dist_fit <- function(values, dist = "lognormal",
                                  samples = 2000, bootstraps = 10,
                                  bootstrap_samples = 250, max_value,
                                  verbose = FALSE) {
  if (!dist %in% c("gamma", "lognormal")) {
    cli_abort(
      c(
        "x" = "Unsupported distribution.",
        "i" = "Only {col_blue(\"lognormal\")} and {col_blue(\"gamma\")}
      distributions are supported"
      )
    )
  }

  ensure_julia()

  ## Make values integer if not
  values <- as.integer(values)
  ## Remove NA values
  values <- values[!is.na(values)]
  ## Filter out negative values
  values <- values[values >= 0]

  if (!missing(max_value)) {
    dist_max <- max_value
  } else {
    dist_max <- max(values)
  }

  # Call Julia's bootstrapped_dist_fit
  julia_data <- juliaEval(sprintf(
    'DataFrame(delay = [%s])',
    paste(as.numeric(values), collapse = ", ")
  ))

  family <- switch(dist,
    lognormal = juliaEval(":lognormal"),
    gamma = juliaEval(":gamma")
  )

  julia_result <- juliaCall(
    "bootstrapped_dist_fit", julia_data,
    family = family,
    max_delay = as.integer(dist_max),
    n_bootstraps = as.integer(bootstraps)
  )

  # Extract the parameter priors from the UncertainDistribution
  param_priors <- juliaEval(
    'let ud = _bootstrapped_result
      [(string(typeof(p)), Distributions.mean(p), Distributions.std(p))
       for p in ud.param_priors]
    end'
  )

  # Build dist_spec from the bootstrap results
  if (dist == "lognormal") {
    # UncertainDistribution has priors on (meanlog, sdlog)
    meanlog_mean <- juliaCall("Distributions.mean",
      juliaCall("getindex",
        juliaCall("getfield", julia_result,
          juliaEval(":param_priors")), 1L))
    meanlog_sd <- juliaCall("Distributions.std",
      juliaCall("getindex",
        juliaCall("getfield", julia_result,
          juliaEval(":param_priors")), 1L))
    sdlog_mean <- juliaCall("Distributions.mean",
      juliaCall("getindex",
        juliaCall("getfield", julia_result,
          juliaEval(":param_priors")), 2L))
    sdlog_sd <- juliaCall("Distributions.std",
      juliaCall("getindex",
        juliaCall("getfield", julia_result,
          juliaEval(":param_priors")), 2L))

    params <- list(
      meanlog = Normal(mean = meanlog_mean, sd = meanlog_sd),
      sdlog = Normal(mean = sdlog_mean, sd = sdlog_sd)
    )
  } else {
    shape_mean <- juliaCall("Distributions.mean",
      juliaCall("getindex",
        juliaCall("getfield", julia_result,
          juliaEval(":param_priors")), 1L))
    shape_sd <- juliaCall("Distributions.std",
      juliaCall("getindex",
        juliaCall("getfield", julia_result,
          juliaEval(":param_priors")), 1L))
    scale_mean <- juliaCall("Distributions.mean",
      juliaCall("getindex",
        juliaCall("getfield", julia_result,
          juliaEval(":param_priors")), 2L))
    scale_sd <- juliaCall("Distributions.std",
      juliaCall("getindex",
        juliaCall("getfield", julia_result,
          juliaEval(":param_priors")), 2L))

    params <- list(
      shape = Normal(mean = shape_mean, sd = shape_sd),
      rate = Normal(mean = 1.0 / scale_mean, sd = scale_sd / scale_mean^2)
    )
  }

  new_dist_spec(params = params, max = dist_max, distribution = dist)
}

#' Estimate a Delay Distribution
#'
#' @description `r lifecycle::badge("maturing")`
#' Estimate a log normal delay distribution from a vector of integer delays.
#' Currently this function is a simple wrapper for [bootstrapped_dist_fit()].
#'
#' @param delays Integer vector of delays
#'
#' @param ... Arguments to pass to internal methods.
#'
#' @return A `<dist_spec>` summarising the bootstrapped distribution
#' @export
#' @seealso [bootstrapped_dist_fit()]
#' @examples
#' \donttest{
#' # bootstraps and samples have been reduced for this example
#' delays <- rlnorm(500, log(5), 1)
#' estimate_delay(delays, samples = 500, bootstraps = 2)
#' }
estimate_delay <- function(delays, ...) {
  bootstrapped_dist_fit(
    values = delays,
    dist = "lognormal", ...
  )
}
