kohonen_unit_distances <- function(model) {
  # Map to one unit at a time to get every distance through kohonen's public API.
  n_units <- nrow(model$codes[[1L]])
  n_samples <- nrow(model$data[[1L]])
  distances <- vapply(seq_len(n_units), function(unit) {
    one_unit <- model
    one_unit$codes <- lapply(model$codes, function(layer) layer[unit, , drop = FALSE])
    kohonen::map(one_unit)$distances
  }, numeric(n_samples))
  matrix(distances, nrow = n_samples, ncol = n_units)
}

for (metric in c("sumofsquares", "euclidean", "manhattan", "tanimoto")) {
  test_that(paste("sample-to-unit distances agree with kohonen for", metric), {
    samples <- rbind(c(0, 1, .5, 0), c(0, NA, NA, NA), c(1, NA, .5, 0), rep(NA, 4))
    codes <- rbind(c(1, 1, .5, 0), c(.5, 0, 1, .49), c(.49, .5, 0, 1))
    model <- structure(list(data = list(samples), codes = list(codes), dist.fcts = metric,
                            user.weights = 1, distance.weights = 1, whatmap = 1L,
                            maxNA.fraction = .75),
                       class = "kohonen")
    expected <- kohonen_unit_distances(model)
    actual <- delimSOM:::compute.layer.sample.to.unit.distance.SOM(
      samples, codes, metric, max_missing_fraction = model$maxNA.fraction
    )
    expect_equal(unname(actual), expected)

    complete <- delimSOM:::compute.layer.sample.to.unit.distance.SOM(
      samples[1L, , drop = FALSE], codes, metric, use_crossproduct = TRUE
    )
    expect_equal(unname(complete), expected[1L, , drop = FALSE])
  })
}

test_that("missing-data scaling identifies the correct runner-up on a two-layer map", {
  grid <- kohonen::somgrid(3, 2, topo = "hexagonal")
  samples <- list(matrix(c(0, NA, NA, NA), nrow = 1L), matrix(0, nrow = 1L))
  codes <- list(
    cbind(c(0, .8, .2, .95, .95, .95), matrix(1, nrow = 6L, ncol = 3L)),
    matrix(c(0, 0, 1, 1, 1, 1), ncol = 1L)
  )
  model <- structure(list(data = samples, codes = codes, grid = grid, whatmap = 1:2,
                          dist.fcts = rep("sumofsquares", 2), user.weights = c(.5, .5),
                          distance.weights = c(.4, .6), maxNA.fraction = .75),
                     class = "kohonen")
  model$unit.classif <- kohonen::map(model)$unit.classif
  ranked_units <- order(kohonen_unit_distances(model)[1L, ])
  expect_equal(c(model$unit.classif, ranked_units[2L]), c(1L, 3L))
  expect_gt(kohonen::unit.distances(grid)[1L, ranked_units[2L]], 1)
  expect_equal(delimSOM:::calculate.topographic.error.SOM(model), 1)
})
