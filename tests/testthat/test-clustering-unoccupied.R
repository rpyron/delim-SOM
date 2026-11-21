test_that("clustering.SOM preserves soft weights for an unoccupied grid cluster", {
  # Use a fixed codebook so the empty cluster does not depend on SOM training.
  codebook <- matrix(c(0, .02, .4, .42, .9, .92), ncol = 1)
  specimens <- matrix(c(.004, .016, .405, .419), ncol = 1,
                      dimnames = list(paste0("s", 1:4), NULL))
  model <- list(codes = list(codebook), data = list(specimens),
                unit.classif = 1:4, user.weights = 1, distance.weights = 1)
  trained <- list(som_models = list(R1 = model), codebook_vectors = list(codebook),
                  replicate_ids = "R1", quantization_error = c(R1 = .01),
                  topographic_error = c(R1 = 0), layer.distance.functions = "sumofsquares")

  fit <- delimSOM::clustering.SOM(
    trained, max.k = 3L, set.k = 3L, clustering.method = "hierarchical+DB",
    quantization.error.quantile = NULL, topographic.error.quantile = NULL,
    calculate.soft.ancestry = TRUE, calculate.variable.importance = FALSE,
    save.SOM.results = FALSE, overwrite.SOM.results = TRUE, verbose = FALSE
  )

  expect_equal(unname(fit$optim_k_vals[1, "R1"]), 3L)
  expect_identical(unname(fit$som_clusters$R1), rep(1:3, each = 2L))
  expected_groups <- c(1L, 1L, 2L, 2L)
  hard <- unname(fit$cluster_assignment[, "R1"])
  expect_identical(outer(hard, hard, `==`), outer(expected_groups, expected_groups, `==`))

  # Compute the soft mass from squared distances to the six map units.
  distances <- outer(as.numeric(specimens), as.numeric(codebook), `-`)^2
  unit_weights <- exp(-distances / median(distances))
  expected_empty_weight <- rowSums(unit_weights[, 5:6]) / rowSums(unit_weights)
  soft <- fit$replicate_ancestry_matrices$R1
  expect_equal(colnames(soft), paste0("Cluster_", 1:3))
  expect_equal(unname(soft[, colnames(soft) == "Cluster_3"]), expected_empty_weight)
})
