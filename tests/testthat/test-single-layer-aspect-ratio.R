aspect_ratio_fixture <- function() {
  angles <- 2 * pi * (seq_len(64L) - 0.5) / 64L
  input <- cbind(matrix(rep(cos(angles), 80L), nrow = 64L),
                 matrix(rep(sin(angles), 16L), nrow = 64L))
  rownames(input) <- paste0("specimen_", seq_len(nrow(input)))
  input
}

train_for_aspect_ratio_test <- function(input, verbose = FALSE) {
  fit <- delimSOM::train.SOM(input, N.steps = 12L, N.replicates = 1L,
                             parallel = FALSE, verbose = verbose)
  grid <- fit$som_models[[1L]]$grid
  c(grid$xdim, grid$ydim)
}

test_that("a wide single layer uses the two leading axes to shape the grid", {
  input <- aspect_ratio_fixture()
  expect_gt(ncol(input), nrow(input))
  expect_equal(train_for_aspect_ratio_test(input), c(10, 4))
})

test_that("the aspect ratio agrees with covariance eigenvalues on full-rank data", {
  angles <- 2 * pi * (seq_len(81L) - 0.5) / 81L
  input <- cbind(x = cos(angles), y = 0.8 * cos(angles) + 0.6 * sin(angles))
  rownames(input) <- paste0("specimen_", seq_len(nrow(input)))
  eigenvalues <- eigen(stats::cov(input), symmetric = TRUE, only.values = TRUE)$values
  expect_equal(delimSOM:::estimate.single.layer.aspect.ratio.SOM(input),
               sqrt(eigenvalues[1L] / eigenvalues[2L]))
  expect_equal(train_for_aspect_ratio_test(input), c(11, 4))
})

test_that("single-layer grid sizing still falls back when a second axis is absent", {
  input <- aspect_ratio_fixture()
  input[] <- input[, 1L]
  expect_message(
    dimensions <- train_for_aspect_ratio_test(input, verbose = TRUE),
    "Using a square-like SOM grid because two axes of variation could not be estimated",
    fixed = TRUE
  )
  expect_equal(dimensions, c(7, 6))
})

test_that("single-layer aspect estimation still uses the mean-imputation path", {
  input <- aspect_ratio_fixture()
  rows <- seq_len(nrow(input))
  input[cbind(rows, rows)] <- NA_real_
  expect_false(any(complete.cases(input)))
  expect_true(all(colMeans(is.na(input)) < 0.25))
  expect_equal(train_for_aspect_ratio_test(input), c(10, 4))
})
