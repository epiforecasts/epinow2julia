#' Simulate secondary observations from primary observations
#'
#' Simulations are done from a given trajectory of primary observations by
#' applying any given delays and observation parameters.
#'
#' In order to simulate, all parameters that are specified such as the mean and
#' standard deviation of delays or observation scaling, must be fixed.
#' Uncertain parameters are not allowed.
#'
#' A function of the same name that was previously based on a reimplementation
#' of that model in R with potentially time-varying scalings and delays is
#' available as `convolve_and_scale()
#' @param primary a data frame of primary reports (column `primary`) by date
#'   (column `date`). Column `primary` must be numeric and `date` must be in
#'   date format.
#' @inheritParams simulate_infections
#' @inheritParams estimate_secondary
#' @importFrom checkmate assert_data_frame assert_date assert_numeric
#'   assert_subset
#' @importFrom cli cli_abort
#' @return A data.table of simulated secondary observations (column `secondary`)
#'   by date.
#' @export
#' @examples
#' \dontrun{
#' ## load data.table to manipulate `example_confirmed` below
#' library(data.table)
#' cases <- as.data.table(example_confirmed)[, primary := confirm]
#' sim <- simulate_secondary(
#'   cases,
#'   delays = delay_opts(fix_parameters(example_reporting_delay)),
#'   obs = obs_opts(family = "poisson")
#' )
#' }
simulate_secondary <- function(primary,
                               day_of_week_effect = NULL,
                               secondary = secondary_opts(),
                               delays = delay_opts(),
                               truncation = trunc_opts(),
                               obs = obs_opts(),
                               CrIs = c(0.2, 0.5, 0.9),
                               backend = "rstan") {
  assert_data_frame(primary, any.missing = FALSE)
  assert_subset(c("date", "primary"), colnames(primary))
  assert_date(primary$date)
  assert_numeric(primary$primary, lower = 0)
  assert_class(secondary, "secondary_opts")
  assert_class(delays, "delay_opts")
  assert_class(obs, "obs_opts")

  # Ensure Julia backend is ready
  ensure_julia()

  # Convert primary to Julia DataFrame
  dates_str <- as.character(primary$date)
  primary_vals <- as.integer(primary$primary)
  julia_primary <- juliaEval(sprintf(
    'DataFrame(date = Date.([%s]), primary = [%s])',
    paste(sprintf('"%s"', dates_str), collapse = ", "),
    paste(primary_vals, collapse = ", ")
  ))

  # Convert opts
  julia_delays <- r_delay_opts_to_julia(delays)
  julia_obs <- r_obs_opts_to_julia(obs)
  julia_secondary <- r_secondary_opts_to_julia(secondary)

  # Call Julia's simulate_secondary, converting dates to strings
  juliaEval('function _sim_sec_wrapper(primary; kwargs...)
    df = simulate_secondary(primary; kwargs...)
    df2 = copy(df)
    df2.date = string.(df2.date)
    df2
  end')
  julia_result <- juliaCall(
    "_sim_sec_wrapper",
    julia_primary,
    delays = julia_delays,
    obs = julia_obs,
    secondary = julia_secondary
  )

  # Convert result to R data.table
  out <- data.table::as.data.table(julia_result)
  if ("date" %in% names(out)) out[, date := as.Date(date)]
  out <- out[, .(date, secondary)]

  out[]
}
