# Snapshot tests for the R -> Julia opts translators.
#
# Each translator is fragile: it produces Julia code (often via sprintf) that
# the Julia side has to parse and accept. A rename in either the R *_opts()
# API or the Julia *_opts() function silently breaks the bridge. These tests
# build a non-default options object on the R side, send it through the
# translator, render its Julia string repr via `show`, and snapshot the
# result. Drift on either side flips the snapshot.

show_julia <- function(x) {
  JuliaConnectoR::juliaCall(
    "sprint",
    JuliaConnectoR::juliaEval("(io, x) -> show(io, MIME(\"text/plain\"), x)"),
    x
  )
}

test_that("r_gt_opts_to_julia round-trips a fixed lognormal", {
  skip_if_no_julia()
  d <- LogNormal(meanlog = 1.6, sdlog = 0.5, max = 14)
  julia_opts <- r_gt_opts_to_julia(gt_opts(fix_parameters(d)))
  expect_snapshot(cat(show_julia(julia_opts)))
})

test_that("r_delay_opts_to_julia round-trips a fixed lognormal", {
  skip_if_no_julia()
  d <- LogNormal(meanlog = 0.6, sdlog = 0.5, max = 10)
  julia_opts <- r_delay_opts_to_julia(delay_opts(fix_parameters(d)))
  expect_snapshot(cat(show_julia(julia_opts)))
})

test_that("r_rt_opts_to_julia round-trips a weekly random walk", {
  skip_if_no_julia()
  julia_opts <- r_rt_opts_to_julia(
    rt_opts(prior = LogNormal(mean = 2, sd = 0.2), rw = 7)
  )
  expect_snapshot(cat(show_julia(julia_opts)))
})

test_that("r_obs_opts_to_julia round-trips negbin + week effect", {
  skip_if_no_julia()
  julia_opts <- r_obs_opts_to_julia(
    obs_opts(family = "negbin", week_effect = TRUE)
  )
  expect_snapshot(cat(show_julia(julia_opts)))
})

test_that("r_forecast_opts_to_julia round-trips a 7-day horizon", {
  skip_if_no_julia()
  julia_opts <- r_forecast_opts_to_julia(forecast_opts(horizon = 7))
  expect_snapshot(cat(show_julia(julia_opts)))
})

test_that("r_inference_opts_to_julia maps sampling args", {
  skip_if_no_julia()
  julia_opts <- r_inference_opts_to_julia(
    inference_opts(samples = 250, warmup = 100, chains = 2)
  )
  expect_snapshot(cat(show_julia(julia_opts)))
})

test_that("r_inference_opts_to_julia forwards adtype", {
  skip_if_no_julia()
  julia_opts <- r_inference_opts_to_julia(
    inference_opts(samples = 50, warmup = 50, chains = 1,
                   adtype = "AutoForwardDiff()")
  )
  ad_str <- JuliaConnectoR::juliaCall(
    "string",
    JuliaConnectoR::juliaCall("getfield", julia_opts,
                              JuliaConnectoR::juliaEval(":adtype"))
  )
  expect_match(ad_str, "ForwardDiff")
})

test_that("stan_opts() still works with a deprecation warning", {
  skip_if_no_julia()
  expect_warning(
    julia_opts <- r_inference_opts_to_julia(
      suppressWarnings(stan_opts(samples = 100, warmup = 50, chains = 1))
    ),
    regexp = NA
  )
  # Check basic shape survives the old-style input
  samples_val <- JuliaConnectoR::juliaCall(
    "getfield", julia_opts, JuliaConnectoR::juliaEval(":samples")
  )
  expect_equal(as.integer(samples_val), 100L)
})
