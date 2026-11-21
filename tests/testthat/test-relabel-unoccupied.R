test_that("soft weights survive relabeling when a grid cluster has no specimens", {
  grid <- c(1L, 2L, 3L, 4L, 5L, 5L)
  specimens <- c(1L, 2L, 3L, 4L, 1L, 4L)
  reference <- c(2L, 4L, 1L, 3L, 2L, 3L)
  weights <- matrix(rep(c(.4, .2, .15, .15, .1), each = 6), nrow = 6,
                    dimnames = list(paste0("s", 1:6), as.character(1:5)))
  align <- function(hard, units, soft) {
    map <- delimSOM:::synchronize.cluster.labels.SOM(hard, units, reference)
    list(hard = unname(map[as.character(hard)]),
         soft = delimSOM:::relabel.soft.ancestry.SOM(soft, map, 5L))
  }

  original <- align(specimens, grid, weights)
  permutation <- c(5L, 1L, 3L, 4L, 2L)
  permuted <- weights[, order(permutation), drop = FALSE]
  colnames(permuted) <- as.character(1:5)
  expect_identical(align(permutation[specimens], permutation[grid], permuted), original)
  expect_identical(original$hard, reference)
  expect_equal(unname(original$soft[, 5]), rep(.1, 6))
  expect_equal(rowSums(original$soft), rowSums(weights))
})
