library("data.table")
library("lifecycle")

if (requireNamespace("future", quietly = TRUE)) {
  withr::defer(future::plan("sequential"), teardown_env())
}

# Disable progressr output during tests
if (requireNamespace("progressr", quietly = TRUE)) {
  progressr::handlers("void")
}

# Test categorisation helpers -----------------------------------------------

#' Check if integration tests should be run
#'
#' Integration tests require Julia and run MCMC inference. By default, these
#' are skipped to speed up test runs. Set EPINOW2_SKIP_INTEGRATION=false and
#' ensure Julia + EpiNow2.jl are available to run them.
#'
#' @return Logical indicating whether to run integration tests
integration_test <- function() {
  skip_integration <- Sys.getenv("EPINOW2_SKIP_INTEGRATION", "true")
  !isTRUE(as.logical(skip_integration))
}

#' Skip test if not running integration tests
#'
#' Helper to skip integration tests with a consistent message.
#' Use at the start of test_that blocks for slow MCMC-based tests.
#'
#' @return Invisibly returns NULL, called for side effect of skipping test
skip_integration <- function() {
  testthat::skip_if_not(integration_test(), "Skipping integration test")
}

#' Skip test if Julia is not available
#'
#' @return Invisibly returns NULL, called for side effect of skipping test
skip_if_no_julia <- function() {
  testthat::skip_if_not(
    julia_available(),
    "Julia backend not available"
  )
}

#' Check if full test suite should be run
#'
#' Full tests include all integration tests and are typically run on a schedule
#' rather than on every commit. Set EPINOW2_FULL_TESTS=true to run them.
#'
#' @return Logical indicating whether to run full test suite
full_tests <- function() {
  isTRUE(as.logical(Sys.getenv("EPINOW2_FULL_TESTS", "false")))
}

# Shared test fixtures -----------------------------------------------------
# Run regional_epinow() once and reuse output for all downstream tests.
# This avoids running MCMC multiple times while testing the full pipeline.

#' Get shared test fixtures
#'
#' Runs regional_epinow() once (lazily) and caches the result.
#' Requires Julia to be available.
#'
#' @return List with regional_epinow output and extracted estimate_infections
#'   objects
get_test_fixtures <- local({
  fixtures <- NULL
  function() {
    if (is.null(fixtures)) {
      if (!julia_available()) {
        return(NULL)
      }
      futile.logger::flog.threshold("FATAL")

      # Create test data with 2 regions
      cases <- EpiNow2::example_confirmed[1:30]
      cases <- data.table::rbindlist(list(
        data.table::copy(cases)[, region := "testland"],
        data.table::copy(cases)[, region := "realland"]
      ))

      # Run regional_epinow once with Julia backend
      suppressWarnings(suppressMessages({
        regional_out <- regional_epinow(
          data = cases,
          generation_time = gt_opts(example_generation_time),
          delays = delay_opts(
            example_incubation_period + example_reporting_delay
          ),
          rt = rt_opts(prior = LogNormal(mean = 2, sd = 0.2)),
          stan = stan_opts(
            samples = 25, warmup = 25,
            chains = 2
          ),
          output = c(
            "regions", "summary", "samples", "plots", "timing",
            "estimate_infections"
          ),
          verbose = FALSE
        )
      }))

      fixtures <<- list(
        regional = regional_out,
        estimate_infections = regional_out$regional$testland,
        estimate_infections_alt = regional_out$regional$realland
      )
    }
    fixtures
  }
})
