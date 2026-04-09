#' Format sample predictions
#'
#' Helper function to format posterior samples into the structure expected by
#' [scoringutils::as_forecast_sample()].
#'
#' @param samples Data.table with date, sample, and value columns
#' @param forecast_date Date when the forecast was made
#' @return Data.table with columns: forecast_date, date, horizon, sample,
#'   predicted
#' @keywords internal
format_sample_predictions <- function(samples, forecast_date) {
  predictions <- samples[, .(date, sample, predicted = value)]
  predictions[, forecast_date := forecast_date]
  predictions[, horizon := as.numeric(date - forecast_date)]
  data.table::setcolorder(
    predictions,
    c("forecast_date", "date", "horizon", "sample", "predicted")
  )
  predictions[]
}

#' Format quantile predictions
#'
#' Helper function to format posterior samples into quantiles in the structure
#' expected by [scoringutils::as_forecast_quantile()].
#'
#' @param samples Data.table with date and value columns
#' @param quantiles Numeric vector of quantile levels
#' @param forecast_date Date when the forecast was made
#' @return Data.table with columns: forecast_date, date, horizon,
#'   quantile_level, predicted
#' @keywords internal
format_quantile_predictions <- function(samples, quantiles, forecast_date) {
  predictions <- samples[
    ,
    .(predicted = quantile(value, probs = quantiles)),
    by = date
  ]
  predictions[, quantile_level := rep(quantiles, .N / length(quantiles))]
  predictions[, forecast_date := forecast_date]
  predictions[, horizon := as.numeric(date - forecast_date)]
  data.table::setcolorder(
    predictions,
    c("forecast_date", "date", "horizon", "quantile_level", "predicted")
  )
  predictions[]
}
