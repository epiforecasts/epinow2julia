#' Get accumulation flags from data
#'
#' @description Returns the `accumulate` column if present, otherwise
#'   a vector of `FALSE` (no accumulation).
#' @param data A data.table that may contain an `accumulate` column.
#' @return A logical vector of length `nrow(data)`.
#' @keywords internal
get_accumulate <- function(data) {
  if ("accumulate" %in% colnames(data)) {
    data$accumulate
  } else {
    rep(FALSE, nrow(data))
  }
}

#' Create Delay Shifted Cases
#'
#' @description `r lifecycle::badge("stable")`
#'
#' This functions creates a data frame of reported cases that has been smoothed
#' using a centred partial rolling average (with a period set by
#' `smoothing_window`) and shifted back in time by some delay. It is used by
#' [estimate_infections()] to generate the mean shifted prior on which the back
#' calculation method (see [backcalc_opts()]) is based.
#'
#' @details
#' The function first shifts all the data back in time by `shift` days (thus
#' discarding the first `shift` days of data) and then applies a centred
#' rolling mean of length `smoothing_window` to the shifted data except for
#' the final period. The final period (the forecast horizon plus half the
#' smoothing window) is instead replaced by a log-linear model fit (with 1
#' added to the data for fitting to avoid zeroes and later subtracted again),
#' projected to the end of the forecast horizon. The initial part of the data
#' (corresponding to the length of the smoothing window) is then removed, and
#' any non-integer resulting values rounded up.
#'
#' @param smoothing_window Numeric, the rolling average smoothing window
#' to apply. Must be odd in order to be defined as a centred average.
#'
#' @param shift Numeric, mean delay shift to apply.
#'
#' @inheritParams estimate_infections
#' @importFrom data.table copy shift frollmean fifelse .N
#' @importFrom stats lm
#' @importFrom runner mean_run
#' @return A `<data.frame>` for shifted reported cases
#' @keywords internal
#' @examples
#' \dontrun{
#' shift <- 7
#' horizon <- 7
#' smoothing_window <- 14
#' ## add NAs for horizon
#' cases <- add_horizon(example_confirmed[1:30], horizon)
#' ## add zeroes initially
#' cases <- data.table::rbindlist(list(
#'   data.table::data.table(
#'     date = seq(
#'       min(cases$date) - 10,
#'       min(cases$date) - 1,
#'       by = "days"
#'     ),
#'     confirm = 0, breakpoint = 0
#'   ),
#'   cases
#' ))
#' create_shifted_cases(cases, shift, smoothing_window, horizon)
#' }
create_shifted_cases <- function(data, shift,
                                 smoothing_window, horizon) {
  shifted_reported_cases <- copy(data)
  ## turn initial NAs into zeroes
  shifted_reported_cases[cumsum(!is.na(confirm)) == 0L, confirm := 0.0]
  ## pad with additional zeroes
  shifted_reported_cases <- pad_reported_cases(data, smoothing_window, 0.0)

  if ("accumulate" %in% colnames(data)) {
    shifted_reported_cases[
      is.na(confirm) & accumulate,
      confirm := 0
    ]
  }
  shifted_reported_cases[
    ,
    confirm := data.table::shift(confirm,
      n = shift,
      type = "lead", fill = NA
    )
  ][
    ,
    confirm := runner::mean_run(
      confirm,
      k = smoothing_window, lag = -floor(smoothing_window / 2)
    )
  ]

  ## Forecast trend on reported cases using the last week of data
  final_period <- shifted_reported_cases[!is.na(confirm)][
    max(1, .N - smoothing_window):.N
  ][
    ,
    t := seq_len(.N)
  ]
  lm_model <- stats::lm(log(confirm + 1) ~ t, data = final_period)
  ## Estimate unreported future infections using a log linear model
  shifted_reported_cases <- shifted_reported_cases[
    date >= min(final_period$date), t := seq_len(.N)
  ][
    ,
    confirm := data.table::fifelse(
      !is.na(t) & t >= 0,
      exp(lm_model$coefficients[1] + lm_model$coefficients[2] * t) - 1,
      confirm
    )
  ][, t := NULL]

  ## Drop median generation interval initial values
  shifted_reported_cases <- shifted_reported_cases[
    ,
    confirm := ceiling(confirm)
  ]
  shifted_reported_cases <- shifted_reported_cases[-(1:smoothing_window)]
  if (anyNA(shifted_reported_cases$confirm)) {
    cli::cli_abort(
      c(
        "!" = "Some values are missing after prior smoothing. Consider
        increasing the smoothing using the {.var prior_window} argument in
        {.fn backcalc_opts}."
      )
    )
  }
  shifted_reported_cases
}

#' Construct the Required Future Rt assumption
#'
#' @description `r lifecycle::badge("stable")`
#' Converts the `future` argument from [rt_opts()] into arguments that can be
#' passed to stan.
#'
#' @param future A character string or integer. This argument indicates how to
#' set future Rt values. Supported options are to project using the Rt model
#' ("project"), to use the latest estimate based on partial data ("latest"),
#' to use the latest estimate based on data that is over 50% complete
#' ("estimate"). If an integer is supplied then the Rt estimate from this many
#' days into the future (or past if negative) past will be used forwards in
#' time.
#'
#' @param delay Numeric mean delay
#' @importFrom rlang arg_match
#' @keywords internal
#' @return A list containing a logical called fixed and an integer called from
create_future_rt <- function(future = c("latest", "project", "estimate"),
                             delay = 0) {
  out <- list(fixed = FALSE, from = 0)
  if (is.character(future)) {
    future <- arg_match(future)
    if (future != "project") {
      out$fixed <- TRUE
      out$from <- ifelse(future == "latest", 0, -delay)
    }
  } else if (is.numeric(future)) {
    out$fixed <- TRUE
    out$from <- as.integer(future)
  }
  out
}

#' Create summary output from infection estimation objects
#'
#' @description `r lifecycle::badge("stable")`
#'
#' This function creates summary output from infection estimation objects
#' (either `estimate_infections` or `forecast_infections`). It is used
#' internally by [summary.estimate_infections()] and
#' [summary.forecast_infections()] to provide a consistent summary interface.
#'
#' @param object An infection estimation object (either from
#'   [estimate_infections()] or [forecast_infections()]).
#'
#' @param type A character vector of data types to return. Defaults to
#'   "snapshot" but also supports "parameters". "snapshot" returns
#'   a summary at a given date (by default the latest date informed by data).
#'   "parameters" returns summarised parameter estimates that can be further
#'   filtered using `params` to show just the parameters of interest and date.
#'
#' @inheritParams summary.estimate_infections
#'
#' @param CrIs Numeric vector of credible intervals to calculate. Defaults
#'   to c(0.2, 0.5, 0.9).
#'
#' @param ... Additional arguments passed to [report_summary()].
#'
#' @return A `<data.frame>` of summary output, either a snapshot summary
#'   (via [report_summary()]) or parameter summaries (via
#'   [calc_summary_measures()]).
#'
#' @importFrom rlang arg_match
#' @seealso [summary.estimate_infections()] [summary.forecast_infections()]
#'   [report_summary()] [calc_summary_measures()]
#' @keywords internal
create_infection_summary <- function(object,
                                     type = c("snapshot", "parameters"),
                                     target_date = NULL, params = NULL,
                                     CrIs = c(0.2, 0.5, 0.9), ...) {
  type <- arg_match(type)

  samples <- get_samples(object)

  summarised <- calc_summary_measures(
    samples,
    summarise_by = c("date", "variable", "strat", "type"),
    order_by = c("variable", "date"),
    CrIs = CrIs
  )

  if (type == "snapshot") {
    if (is.null(target_date)) {
      target_date <- max(object$observations$date)
    } else {
      target_date <- as.Date(target_date)
    }
    out <- report_summary(
      summarised_estimates = summarised[date == target_date],
      rt_samples = samples[variable == "R"][
        date == target_date, .(sample, value)
      ],
      ...
    )
  } else if (type == "parameters") {
    out <- summarised
    if (!is.null(target_date)) {
      out <- out[date == as.Date(target_date)]
    }
    if (!is.null(params)) {
      out <- out[variable %in% params]
    }
  }
  out[]
}
