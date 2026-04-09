#' Stan model compatibility layer
#'
#' @description These functions are retained for backward compatibility.
#'   The Julia backend is used for all model fitting.
#' @keywords internal

#' Load and compile an EpiNow2 cmdstanr model
#'
#' @description `r lifecycle::badge("deprecated")`
#' This function is deprecated. The Julia backend is now used for inference.
#'
#' @param model Character string indicating the model to use.
#' @param dir Character string specifying the path to Stan files.
#' @param verbose Logical. Should verbose messages be shown.
#' @param ... Additional arguments.
#' @return Errors with deprecation message.
#' @export
epinow2_cmdstan_model <- function(model = "estimate_infections",
                                  dir = NULL,
                                  verbose = FALSE,
                                  ...) {
  lifecycle::deprecate_stop(
    "2.0.0",
    "epinow2_cmdstan_model()",
    details = "The Julia backend is now used for inference."
  )
}
