#' Extract all samples from a fitted model
#'
#' @description `r lifecycle::badge("stable")`
#' Extracts posterior samples from a fitted model object. For Julia-backed
#' fits, this calls Julia's `get_samples()` and converts the result.
#'
#' @param stan_fit A Julia result object (from the Julia backend) or a
#'   legacy stanfit object.
#' @param pars Any selection of parameters to extract
#' @param include whether the parameters specified in `pars` should be included
#' (`TRUE`, the default) or excluded (`FALSE`)
#' @return List of data.tables with samples
#' @export
extract_samples <- function(stan_fit, pars = NULL, include = TRUE) {
  ensure_julia()
  # Get samples from Julia result
  julia_df <- juliaCall("get_samples", stan_fit)
  samples_dt <- data.table::as.data.table(julia_df)

  if ("date" %in% names(samples_dt)) {
    samples_dt[, date := as.Date(date)]
  }
  if ("variable" %in% names(samples_dt)) {
    samples_dt[, variable := as.character(variable)]
  }

  # Filter by pars if specified
  if (!is.null(pars)) {
    if (include) {
      samples_dt <- samples_dt[variable %in% pars]
    } else {
      samples_dt <- samples_dt[!variable %in% pars]
    }
  }

  # Convert to list-of-matrices format for backward compatibility.
  # Reshape via dcast so that values are placed by (sample, date) rather
  # than relying on the row order of the long-format table — which is not
  # guaranteed by the Julia conversion and would silently mis-align the
  # matrix if it happened to be sorted (sample, date) instead of
  # (date, sample).
  vars <- unique(samples_dt$variable)
  result <- list()
  for (v in vars) {
    sub <- samples_dt[variable == v]
    if ("date" %in% names(sub) && !all(is.na(sub$date))) {
      wide <- data.table::dcast(
        sub, sample ~ date, value.var = "value"
      )
      sample_ids <- wide$sample
      wide[, sample := NULL]
      mat <- as.matrix(wide)
      rownames(mat) <- sample_ids
      result[[v]] <- mat
    } else {
      # Scalar parameter: vector ordered by sample
      data.table::setorder(sub, sample)
      result[[v]] <- sub$value
    }
  }

  result
}

#' Extract a parameter summary from a fitted model
#'
#' @description `r lifecycle::badge("stable")`
#' Extracts summarised parameter posteriors from a fitted model.
#'
#' @param fit A Julia result object.
#' @param params A character vector of parameters to extract.
#' @param var_names Logical defaults to `FALSE`. Should variables be named.
#' @inheritParams calc_summary_measures
#' @return A `<data.table>` summarising parameter posteriors.
#' @export
#' @importFrom data.table as.data.table :=
extract_stan_param <- function(fit, params = NULL,
                               CrIs = c(0.2, 0.5, 0.9), var_names = FALSE) {
  ensure_julia()

  samples <- extract_samples(fit, pars = params)

  # Build summary data.table from samples
  CrIs <- sort(CrIs)
  sym_CrIs <- c(0.5, 0.5 - CrIs / 2, 0.5 + CrIs / 2)
  sym_CrIs <- sort(sym_CrIs)
  CrI_names <- round(100 * CrIs, 0)
  CrI_labels <- c(
    paste0("lower_", rev(CrI_names)), "median", paste0("upper_", CrI_names)
  )

  results <- lapply(names(samples), function(v) {
    vals <- as.numeric(samples[[v]])
    if (length(vals) == 0) return(NULL)

    qs <- quantile(vals, probs = sym_CrIs, na.rm = TRUE)
    dt <- data.table::data.table(
      variable = v,
      mean = mean(vals, na.rm = TRUE),
      se_mean = stats::sd(vals, na.rm = TRUE) / sqrt(length(vals)),
      sd = stats::sd(vals, na.rm = TRUE)
    )
    for (i in seq_along(CrI_labels)) {
      dt[[CrI_labels[i]]] <- unname(qs[i])
    }
    dt
  })

  param_summary <- data.table::rbindlist(results)
  if (!var_names && "variable" %in% names(param_summary)) {
    param_summary[, variable := NULL]
  }
  param_summary
}

#' Generate initial conditions from a fit
#'
#' @description `r lifecycle::badge("experimental")`
#' Extracts posterior samples to use as initial conditions for a new fit.
#'
#' @param fit A Julia result object.
#' @param current_inits A function returning initial conditions.
#' @param exclude_list Parameters to exclude.
#' @param samples Number of posterior samples.
#' @return A function that returns initial conditions.
#' @export
extract_inits <- function(fit, current_inits,
                          exclude_list = NULL,
                          samples = 50) {
  cli::cli_warn(
    "extract_inits() is not fully supported with the Julia backend."
  )
  current_inits
}

#' Extract samples for a latent state from a model
#'
#' @param param Character string indicating the latent state to extract
#' @param samples Extracted model samples (list of arrays)
#' @param dates A vector of dates for the time dimension
#' @return A `<data.frame>` with time, date, sample, and value columns
#' @importFrom data.table melt as.data.table
#' @keywords internal
extract_latent_state <- function(param, samples, dates) {
  if (!(param %in% names(samples))) {
    return(NULL)
  }

  param_data <- samples[[param]]

  if (is.matrix(param_data)) {
    param_df <- data.table::as.data.table(t(param_data))
    param_df[, time := seq_len(.N)]
    param_df <- data.table::melt(param_df,
      id.vars = "time", variable.name = "var"
    )
    param_df[, var := NULL]
    param_df[, sample := seq_len(.N), by = .(time)]
    param_df[, date := dates[time], by = .(sample)]
  } else {
    # Vector (scalar per sample)
    param_df <- data.table::data.table(
      time = 1L,
      sample = seq_along(param_data),
      value = param_data,
      date = if (length(dates) >= 1) dates[1] else NA
    )
  }

  param_df[, .(time, date, sample, value)]
}

#' Extract samples from all parameters
#' @param samples Extracted model samples
#' @param args Model arguments list
#' @return A `<data.table>` with variable, sample, value columns
#' @keywords internal
extract_parameters <- function(samples, args) {
  NULL
}

#' Extract samples from all delay parameters
#' @param samples Extracted model samples
#' @param args Model arguments list
#' @return A `<data.table>` with variable, sample, value columns
#' @keywords internal
extract_delays <- function(samples, args) {
  NULL
}
