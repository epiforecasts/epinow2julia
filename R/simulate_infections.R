#' Simulate infections using the renewal equation
#'
#' Simulations are done from given initial infections and, potentially
#' time-varying, reproduction numbers. Delays and parameters of the observation
#' model can be specified using the same options as in [estimate_infections()].
#'
#' In order to simulate, all parameters that are specified such as the mean and
#' standard deviation of delays or observation scaling, must be fixed.
#' Uncertain parameters are not allowed.
#'
#' @param R a data frame of reproduction numbers (column `R`) by date (column
#'   `date`). Column `R` must be numeric and `date` must be in date format. If
#'   not all days between the first and last day in the `date` are present,
#'   it will be assumed that R stays the same until the next given date.
#' @param initial_infections numeric; the initial number of infections (i.e.
#'   before `R` applies). Note that results returned start the day after, i.e.
#'   the initial number of infections is not reported again. See also
#'   `seeding_time`
#' @param day_of_week_effect either `NULL` (no day of the week effect) or a
#'   numerical vector of length specified in [obs_opts()] as `week_length`
#'   (default: 7) if `week_effect` is set to TRUE. Each element of the vector
#'   gives the weight given to reporting on this day (normalised to 1).
#'   The default is `NULL`.
#' @param seeding_time Integer; the number of days before the first time point
#'   of `R`; default is `NULL`, in which case it is set to the maximum of the
#'   generation time. The minimum is 1 , i.e. the first reproduction number
#'   given applies on the day after the index cases given by
#'   `initial_infections`. If the generation time is longer than 1 day on
#'   average, a seeding time of 1 will always lead to an initial decline (as
#'   there are no infections before the initial ones). Instead, if this is
#'   greater than 1, an initial part of the epidemic (before the first value of
#'   R given) of `seeding_time` days is assumed to have followed exponential
#'   growth roughly in line with the growth rate implied by the first value of
#'   R.
#' @inheritParams estimate_infections
#' @inheritParams calc_CrIs
#' @inheritParams rt_opts
#' @inheritParams stan_opts
#' @importFrom checkmate assert_data_frame assert_date assert_numeric
#'   assert_subset assert_integer
#' @importFrom data.table data.table merge.data.table nafill rbindlist
#' @importFrom cli cli_abort
#' @return A data.table of simulated infections (variable `infections`) and
#'   reported cases (variable `reported_cases`) by date.
#' @export
#' @examples
#' \donttest{
#' R <- data.frame(
#'   date = seq.Date(as.Date("2023-01-01"), length.out = 14, by = "day"),
#'   R = c(rep(1.2, 7), rep(0.8, 7))
#' )
#' sim <- simulate_infections(
#'   R = R,
#'   initial_infections = 100,
#'   generation_time = generation_time_opts(
#'     fix_parameters(example_generation_time)
#'   ),
#'   delays = delay_opts(fix_parameters(example_reporting_delay)),
#'   obs = obs_opts(family = "poisson")
#' )
#' }
simulate_infections <- function(R,
                                initial_infections,
                                day_of_week_effect = NULL,
                                generation_time = generation_time_opts(),
                                delays = delay_opts(),
                                truncation = trunc_opts(),
                                obs = obs_opts(),
                                CrIs = c(0.2, 0.5, 0.9),
                                backend = "rstan",
                                seeding_time = NULL,
                                pop = Fixed(0),
                                pop_period = c("forecast", "all"),
                                pop_floor = 1.0,
                                growth_method = c("infections",
                                                  "infectiousness")) {
  if (is.numeric(pop)) {
    lifecycle::deprecate_stop(
      "1.9.0",
      "simulate_infections(pop = 'must be a `<dist_spec>`')",
      details = paste(
        "Population size must now be specified as a distribution.",
        "For a fixed known population, wrap the value with `Fixed()`.",
        "For example: `simulate_infections(..., pop = Fixed(1000000))`."
      )
    )
  }
  assert_class(pop, "dist_spec")
  pop_period <- arg_match(pop_period)

  ## check inputs
  assert_data_frame(R, any.missing = FALSE)
  assert_subset(c("date", "R"), colnames(R))
  assert_date(R$date)
  assert_numeric(R$R, lower = 0)
  assert_numeric(initial_infections, lower = 0)
  assert_class(generation_time, "generation_time_opts")
  assert_class(delays, "delay_opts")
  assert_class(obs, "obs_opts")

  # Ensure Julia backend is ready
  ensure_julia()

  # Convert R trajectory to Julia DataFrame
  r_dates <- as.character(R$date)
  r_vals <- as.numeric(R$R)
  julia_r_traj <- juliaEval(sprintf(
    'DataFrame(date = Date.([%s]), R = [%s])',
    paste(sprintf('"%s"', r_dates), collapse = ", "),
    paste(r_vals, collapse = ", ")
  ))

  # Convert opts
  julia_gt <- r_gt_opts_to_julia(generation_time)
  julia_delays <- r_delay_opts_to_julia(delays)
  julia_obs <- r_obs_opts_to_julia(obs)

  pop_val <- if (inherits(pop, "dist_spec") && get_distribution(pop) != "fixed") {
    mean(pop, ignore_uncertainty = TRUE)
  } else if (inherits(pop, "dist_spec")) {
    as.numeric(pop$parameters$value)
  } else {
    0.0
  }

  # Call Julia's simulate_infections, converting dates to strings
  # to avoid JuliaConnectoR Date serialisation issue
  juliaEval('function _sim_inf_wrapper(rt_traj; kwargs...)
    df = simulate_infections(rt_traj; kwargs...)
    df2 = copy(df)
    df2.date = string.(df2.date)
    df2
  end')
  julia_result <- juliaCall(
    "_sim_inf_wrapper",
    julia_r_traj,
    generation_time = julia_gt,
    delays = julia_delays,
    obs = julia_obs,
    initial_infections = as.numeric(initial_infections),
    pop = as.numeric(pop_val)
  )

  # Convert Julia DataFrame to R data.table
  out <- data.table::as.data.table(julia_result)
  if ("date" %in% names(out)) out[, date := as.Date(date)]

  # Reshape to long format matching expected output
  infections_dt <- out[, .(date, variable = "infections", value = infections)]
  reports_dt <- out[, .(date, variable = "reported_cases", value = reports)]
  out <- rbindlist(list(infections_dt, reports_dt))

  out[]
}

#' Forecast infections from a given fit and trajectory of the time-varying
#' reproduction number
#'
#' @description `r lifecycle::badge("stable")`
#' This function simulates infections using an existing fit to observed cases
#' but with a modified time-varying reproduction number. This can be used to
#' explore forecast models or past counterfactuals. Simulations can be run in
#' parallel using [future::plan()].
#'
#' @param estimates The \code{estimates} element of an [epinow()] run that
#' has been done with output = "fit", or the result of
#' [estimate_infections()] with \code{return_fit} set to TRUE.
#'
#' @param model A compiled stan model as returned by [rstan::stan_model()].
#'
#' @param R A numeric vector of reproduction numbers; these will overwrite the
#' reproduction numbers contained in \code{estimates}, except elements set to
#' NA. Alternatively accepts a `<data.frame>` containing at least `date` and
#' `value` (integer) variables and optionally `sample`. More (or fewer) days
#' than in the original fit can be simulated.
#'
#' @param samples Numeric, number of posterior samples to simulate from. The
#' default is to use all samples in the `estimates` input.
#'
#' @param batch_size Numeric, defaults to 10. Size of batches in which to
#' simulate. May decrease run times due to reduced IO costs but this is still
#' being evaluated. If set to NULL then all simulations are done at once.
#'
#' @param verbose Logical defaults to [interactive()]. If the `progressr`
#' package is available, a progress bar will be shown.
#' @inheritParams stan_opts
#' @importFrom purrr list_transpose map safely compact
#' @importFrom data.table rbindlist as.data.table
#' @importFrom lubridate days
#' @importFrom checkmate assert_class assert_names test_numeric test_data_frame
#' assert_numeric assert_integerish assert_logical
#' @importFrom cli cli_abort
#' @return A `<forecast_infections>` object containing simulated infections and
#' cases from the specified scenario. The structure is similar to
#' [estimate_infections()] output but contains `samples` rather than `fit`.
#' @seealso [generation_time_opts()] [delay_opts()] [rt_opts()]
#' [estimate_infections()] [trunc_opts()] [stan_opts()] [obs_opts()]
#' [gp_opts()]
#' @export
#' @examples
#' \donttest{
#' # set number of cores to use
#' old_opts <- options()
#' options(mc.cores = ifelse(interactive(), 4, 1))
#'
#' # get example case counts
#' reported_cases <- example_confirmed[1:40]
#'
#' # fit model to data to recover Rt estimates
#' # samples and calculation time have been reduced for this example
#' # for real analyses, use at least samples = 2000
#' est <- estimate_infections(reported_cases,
#'   generation_time = generation_time_opts(example_generation_time),
#'   delays = delay_opts(example_incubation_period + example_reporting_delay),
#'   rt = rt_opts(prior = LogNormal(mean = 2, sd = 0.1), rw = 7),
#'   obs = obs_opts(scale = Normal(mean = 0.1, sd = 0.01)),
#'   gp = NULL,
#'   forecast = forecast_opts(horizon = 0),
#'   stan = stan_opts(samples = 100, warmup = 200)
#' )
#'
#' # update Rt trajectory and simulate new infections using it
#' # keeping the first 30 days' estimates and adding a 10-day forecast
#' R <- c(rep(NA_real_, 30), rep(0.8, 10))
#' sims <- forecast_infections(est, R)
#' plot(sims)
#'
#' options(old_opts)
#' }
forecast_infections <- function(estimates,
                                R = NULL,
                                model = NULL,
                                samples = NULL,
                                batch_size = 10,
                                backend = "rstan",
                                verbose = interactive()) {
  ## check inputs
  assert_class(estimates, "estimate_infections")
  if (!(test_numeric(R, lower = 0, null.ok = TRUE) ||
          test_data_frame(R, null.ok = TRUE))) {
    cli_abort(
      c(
        "!" = "R must either be a {.cls numeric} vector or
        a {.cls data.frame}."
      )
    )
  }
  if (test_data_frame(R)) {
    assert_names(names(R), must.include = c("date", "value"))
    assert_numeric(R$value, lower = 0)
  }
  assert_integerish(samples, lower = 1, null.ok = TRUE)
  assert_logical(verbose)

  ensure_julia()

  # Extract fitted R samples to fill NA values
  fitted_samples <- get_samples(estimates)
  fitted_r <- fitted_samples[variable == "R"]
  # Use median R per date as the baseline
  fitted_r_median <- fitted_r[, .(R = median(value)), by = date]
  data.table::setorder(fitted_r_median, date)

  if (is.numeric(R) && !inherits(R, "data.frame")) {
    # Numeric vector: build dates and fill NAs from fitted values
    n_fitted <- nrow(fitted_r_median)
    n_total <- length(R)
    max_date <- max(estimates$observations$date, na.rm = TRUE)

    if (n_total > n_fitted) {
      # Extend dates for forecast horizon
      extra_dates <- seq(max_date + 1, by = "day",
        length.out = n_total - n_fitted)
      r_dates <- c(fitted_r_median$date, extra_dates)
    } else {
      r_dates <- fitted_r_median$date[seq_len(n_total)]
    }
    r_vals <- R
    # Fill NAs with fitted median R
    for (i in seq_along(r_vals)) {
      if (is.na(r_vals[i]) && i <= n_fitted) {
        r_vals[i] <- fitted_r_median$R[i]
      } else if (is.na(r_vals[i])) {
        r_vals[i] <- fitted_r_median$R[n_fitted]
      }
    }
    r_df <- data.table::data.table(date = r_dates, R = r_vals)
  } else {
    r_df <- data.table::as.data.table(R)
    if (!"R" %in% names(r_df) && "value" %in% names(r_df)) {
      r_df[, R := value]
    }
  }

  # Build Julia R trajectory DataFrame (no NAs)
  r_dates_str <- as.character(r_df$date)
  r_vals <- as.numeric(r_df$R)
  julia_r_traj <- juliaEval(sprintf(
    'DataFrame(date = Date.([%s]), R = [%s])',
    paste(sprintf('"%s"', r_dates_str), collapse = ", "),
    paste(r_vals, collapse = ", ")
  ))

  # Call Julia's forecast_infections using the stored Julia result
  juliaEval('function _forecast_inf_wrapper(result, r_traj)
    df = forecast_infections(result, r_traj)
    df2 = copy(df)
    if hasproperty(df2, :date)
      df2.date = string.(df2.date)
    end
    df2
  end')
  julia_result <- juliaCall(
    "_forecast_inf_wrapper",
    estimates$fit,
    julia_r_traj
  )

  # Convert to R data.table
  forecast_dt <- data.table::as.data.table(julia_result)
  if ("date" %in% names(forecast_dt)) forecast_dt[, date := as.Date(date)]

  # Build output in expected format
  format_out <- list()
  format_out$summarised <- forecast_dt
  format_out$observations <- estimates$observations

  # Create samples placeholder
  format_out$samples <- data.table::data.table(
    variable = character(0), time = integer(0), date = as.Date(character(0)),
    sample = integer(0), value = numeric(0), strat = character(0),
    type = character(0)
  )

  class(format_out) <- c("forecast_infections", class(format_out))
  format_out
}
