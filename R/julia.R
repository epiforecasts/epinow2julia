#' Julia backend initialisation
#'
#' @description Internal functions for managing the Julia backend used by
#'   EpiNow2 for Bayesian inference via Turing.jl.
#'
#' @keywords internal

# Package-level environment for Julia state
.julia_env <- new.env(parent = emptyenv())
.julia_env$initialised <- FALSE

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
#'
#' @importFrom cli cli_inform cli_abort
#' @importFrom JuliaConnectoR juliaEval juliaCall juliaImport
#' @return Invisible `NULL`, called for side effects.
#' @keywords internal
setup_julia <- function(project_path = NULL) {
  if (.julia_env$initialised) {
    return(invisible(NULL))
  }

  if (!requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    cli_abort(
      c(
        "!" = "Package {.pkg JuliaConnectoR} is required but not installed.",
        "i" = "Install with {.code install.packages(\"JuliaConnectoR\")}."
      )
    )
  }

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

  cli_inform("Setting up Julia backend...")

  # Set JULIA_PROJECT so JuliaConnectoR starts Julia in the right env
  Sys.setenv(JULIA_PROJECT = project_path)

  # Load the EpiNow2 Julia module
  JuliaConnectoR::juliaEval("using EpiNow2")
  JuliaConnectoR::juliaEval("using DataFrames")
  JuliaConnectoR::juliaEval("using Distributions")
  JuliaConnectoR::juliaEval("using Dates")

  .julia_env$initialised <- TRUE
  cli_inform("Julia backend ready.")

  invisible(NULL)
}

#' Ensure Julia is initialised
#'
#' @description Checks whether Julia has been set up and calls
#'   [setup_julia()] if not. Should be called at the start of every
#'   user-facing function that needs the Julia backend.
#' @keywords internal
ensure_julia <- function() {
  if (!.julia_env$initialised) {
    setup_julia()
  }
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
