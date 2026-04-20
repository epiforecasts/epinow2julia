test_that("extract_samples reshape preserves (sample, date) alignment", {
  # Build a long-format samples table with values that encode their
  # (sample, date) coordinates: value = 100 * sample + date_idx. Reshape
  # via the same dcast logic used in extract_samples() and check that
  # mat[s, t] decodes back to (s, t). Catches the silent transposition
  # that the previous matrix(byrow = FALSE) reshape was vulnerable to.
  dates <- as.Date("2020-01-01") + 0:4
  n_samples <- 7L
  long <- data.table::CJ(
    sample = seq_len(n_samples),
    date = dates
  )
  long[, value := 100 * sample + as.integer(date - min(date)) + 1L]
  # Shuffle the rows so any reshape that relies on row order will fail
  set.seed(1)
  long <- long[sample(.N)]

  wide <- data.table::dcast(long, sample ~ date, value.var = "value")
  sample_ids <- wide$sample
  wide[, sample := NULL]
  mat <- as.matrix(wide)
  rownames(mat) <- sample_ids

  # mat must be n_samples x n_times with rows ordered by sample id
  expect_equal(dim(mat), c(n_samples, length(dates)))
  expect_equal(as.integer(rownames(mat)), seq_len(n_samples))

  for (s in seq_len(n_samples)) {
    for (t in seq_along(dates)) {
      expect_equal(mat[s, t], 100 * s + t,
        info = sprintf("mismatch at sample=%d, date_idx=%d", s, t)
      )
    }
  }
})
