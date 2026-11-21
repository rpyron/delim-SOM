grid_size_fixture <- function(n) {
  angles <- 2 * pi * (seq_len(n) - 0.5) / n
  input <- cbind(x = cos(angles), y = sin(angles))
  rownames(input) <- paste0("specimen_", seq_len(n))
  input
}

train_for_grid_size_test <- function(input, ..., verbose = FALSE) {
  delimSOM::train.SOM(input, N.steps = 12L, N.replicates = 1L,
                      parallel = FALSE, verbose = verbose, ...)
}

fitted_grid_size <- function(trained) {
  grid <- trained$som_models[[1L]]$grid
  c(grid$xdim, grid$ydim)
}

test_that("oversized automatic grids use the retained specimen count", {
  input <- grid_size_fixture(16L)
  input[13:16, ] <- NA_real_
  expect_message(
    trained <- train_for_grid_size_test(input, verbose = TRUE),
    "Reduced automatic SOM grid to 4 x 3 (12 units) to fit 12 retained samples",
    fixed = TRUE
  )
  expect_equal(fitted_grid_size(trained), c(4L, 3L))
  expect_identical(rownames(trained$input_data[[1L]]), rownames(input)[1:12])
  explicit <- train_for_grid_size_test(input, grid.size = c(4, 3))
  expect_identical(trained$som_models, explicit$som_models)
})

test_that("rounding an automatic grid cannot exceed the specimen count", {
  # The target is 8 units, but ordinary rounding proposes 3 x 3.
  trained <- train_for_grid_size_test(grid_size_fixture(8L), grid.multiplier = 3)
  expect_true(all(fitted_grid_size(trained) >= 2L))
  expect_lte(prod(fitted_grid_size(trained)), 8L)
  expect_equal(fitted_grid_size(train_for_grid_size_test(grid_size_fixture(4L))),
               c(2L, 2L))
})

test_that("feasible automatic grids still use the original dimensions", {
  expect_equal(fitted_grid_size(train_for_grid_size_test(grid_size_fixture(64L))),
               c(7L, 6L))
  expect_equal(fitted_grid_size(train_for_grid_size_test(grid_size_fixture(16L),
                                                        grid.multiplier = 4)),
               c(4L, 4L))
  # A target of 22 rounds to a feasible 4 x 5 grid for 20 specimens.
  expect_equal(fitted_grid_size(train_for_grid_size_test(grid_size_fixture(20L))),
               c(4L, 5L))
})

test_that("infeasible explicit grids and cohorts below the minimum still fail", {
  expect_error(
    train_for_grid_size_test(grid_size_fixture(12L), grid.size = c(4L, 4L)),
    "requested SOM grid size 4 x 4 (16 units) exceeds number of samples (12)",
    fixed = TRUE
  )
  expect_error(
    train_for_grid_size_test(grid_size_fixture(3L)),
    "a minimum 2 x 2 grid requires at least 4 retained samples (found 3)",
    fixed = TRUE
  )
})
