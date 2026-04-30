#' Julia backend initialisation
#'
#' @description Internal functions for managing the Julia backend used by
#'   EpiNow2 for Bayesian inference via Turing.jl. Heavy lifting is
#'   delegated to [juliaready::julia_ready()].
#'
#' @keywords internal

# Package-level environment for Julia state
.julia_env <- new.env(parent = emptyenv())

#' Set up Julia for use with EpiNow2
#'
#' @description Initialises Julia and loads the EpiNow2.jl package. Called
#'   automatically on first use. The path to EpiNow2.jl can be configured via
#'   the `EpiNow2.julia_project` option or the `EPINOW2_JULIA_PROJECT`
#'   environment variable.
#'
#' @param project_path Character string giving the path to the EpiNow2.jl
#'   project. If `NULL` (default), uses the `EpiNow2.julia_project` option or
#'   the `EPINOW2_JULIA_PROJECT` environment variable.
#' @param threads Integer; number of threads Julia should start with.
#'   Controls chain-level parallelism via `MCMCThreads()`. If `NULL`
#'   (default), respects an existing `JULIA_NUM_THREADS` environment
#'   variable; otherwise uses one fewer than the number of physical cores
#'   (minimum 1). Must be set before Julia starts.
#'
#' @importFrom cli cli_inform cli_abort
#' @return Invisible `NULL`, called for side effects.
#' @keywords internal
setup_julia <- function(project_path = NULL, threads = NULL) {
  if (isTRUE(.julia_env$ready)) return(invisible(NULL))

  # Determine project path
  if (is.null(project_path)) {
    project_path <- getOption(
      "EpiNow2.julia_project",
      Sys.getenv("EPINOW2_JULIA_PROJECT", unset = "")
    )
    if (project_path == "") {
      cli_abort(
        c(
          "!" = "No EpiNow2.jl project path specified.",
          "i" = paste(
            "Set {.code options(EpiNow2.julia_project = \"/path/to/epinow2.jl\")}",
            "or the {.envvar EPINOW2_JULIA_PROJECT} environment variable."
          )
        )
      )
    }
  }

  # Determine thread count. Must be set before Julia is launched, since
  # JULIA_NUM_THREADS is read once at startup.
  existing_threads <- Sys.getenv("JULIA_NUM_THREADS", unset = "")
  if (is.null(threads)) {
    if (nzchar(existing_threads)) {
      threads <- existing_threads
    } else {
      n_cores <- tryCatch(
        parallel::detectCores(logical = FALSE),
        error = function(e) NA_integer_
      )
      threads <- if (is.na(n_cores)) 1L else max(1L, n_cores - 1L)
    }
  }
  Sys.setenv(JULIA_NUM_THREADS = as.character(threads))

  # Pin the project so JuliaConnectoR / juliaready use it.
  Sys.setenv(JULIA_PROJECT = project_path)

  cli_inform("Setting up Julia backend ({threads} thread{?s})...")

  # juliaready handles binary detection, idempotency, and the lazy guard.
  # We use the project-pinned environment via JULIA_PROJECT (set above)
  # so the requested packages are resolved within that project.
  juliaready::julia_ready(
    packages  = c("EpiNow2", "DataFrames", "Distributions", "Dates"),
    state_env = .julia_env,
    install   = FALSE,
    verbose   = FALSE
  )

  n_threads <- as.integer(juliaready::eval_julia("Threads.nthreads()"))
  .julia_env$threads <- n_threads
  cli_inform("Julia backend ready ({n_threads} thread{?s}).")

  invisible(NULL)
}

#' Ensure Julia is initialised
#'
#' @description Lazy guard: calls [setup_julia()] on first use.
#' @keywords internal
ensure_julia <- function() {
  juliaready::ensure_julia(.julia_env, setup_julia)
}

#' Check whether Julia is available
#'
#' @description Tests whether Julia and the EpiNow2.jl package are
#'   available. Useful for conditional test execution.
#' @return Logical; `TRUE` if Julia is ready, `FALSE` otherwise.
#' @export
julia_available <- function() {
  tryCatch(
    {
      ensure_julia()
      TRUE
    },
    error = function(e) FALSE
  )
}

#' Number of threads the Julia backend is using
#'
#' @description Reports how many threads the Julia backend was started
#'   with. Chain-level parallelism requires threads > 1.
#' @return Integer number of threads, or `NA` if Julia has not been
#'   initialised.
#' @export
julia_threads <- function() {
  if (!isTRUE(.julia_env$ready)) {
    return(NA_integer_)
  }
  .julia_env$threads
}
