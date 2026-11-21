occupied_group_model <- function(specimens, codes, winners) {
  list(data = list(specimens), codes = list(codes), unit.classif = winners,
       user.weights = 1, distance.weights = 1)
}

fit_occupied_group_models <- function(models, set.k = NULL, max.k = 3L, verbose = FALSE) {
  ids <- names(models)
  trained <- list(
    som_models = models, replicate_ids = ids,
    codebook_vectors = list(do.call(rbind, lapply(models, function(model) model$codes[[1L]]))),
    quantization_error = setNames(rep(0, length(ids)), ids),
    topographic_error = setNames(rep(0, length(ids)), ids),
    layer.distance.functions = "sumofsquares"
  )
  delimSOM::clustering.SOM(
    trained, set.k = set.k, max.k = max.k, clustering.method = "kmeans+BICthreshold",
    BIC.thresh = 3, quantization.error.quantile = NULL, topographic.error.quantile = NULL,
    calculate.variable.importance = FALSE, verbose = verbose
  )
}

occupied_group_fixture <- function() {
  profiles <- rep(c(0, .42), each = 3L)
  specimens <- cbind(x = profiles, y = profiles)
  rownames(specimens) <- paste0("s", seq_len(nrow(specimens)))
  first <- c(0, .02, .42, .9)
  second <- c(0, .02, .98, 1)
  list(
    R1 = occupied_group_model(specimens, cbind(x = first, y = first), rep(c(1L, 3L), each = 3L)),
    R2 = occupied_group_model(specimens, cbind(x = second, y = second), rep(c(1L, 2L), each = 3L))
  )
}

first_seen_labels <- function(labels) match(labels, unique(labels))

expected_grid_bic <- function(model, groups) {
  codes <- model$codes[[1L]]
  withinss <- sum(vapply(split(seq_len(nrow(codes)), groups), function(indices) {
    group <- codes[indices, , drop = FALSE]
    sum(sweep(group, 2L, colMeans(group))^2)
  }, numeric(1)))
  n_units <- nrow(codes)
  n_units * log((withinss + .Machine$double.eps) / n_units) +
    log(n_units) * length(unique(groups))
}

expected_soft_weight <- function(model, units) {
  samples <- model$data[[1L]]
  codes <- model$codes[[1L]]
  distances <- vapply(seq_len(nrow(codes)), function(i) {
    rowSums(sweep(samples, 2L, codes[i, ])^2)
  }, numeric(nrow(samples)))
  weights <- exp(-distances / median(distances))
  unname(rowSums(weights[, units, drop = FALSE]) / rowSums(weights))
}

test_that("a grid with three clusters can have only two occupied specimen groups", {
  models <- occupied_group_fixture()["R1"]
  expect_equal(nrow(unique(models$R1$data[[1L]])), 2L)
  fit <- fit_occupied_group_models(models, set.k = 3L)

  grid <- c(1L, 1L, 2L, 3L)
  expect_equal(unname(fit$optim_k_vals[1L, "R1"]), 3L)
  expect_identical(first_seen_labels(fit$som_clusters$R1), grid)
  expect_equal(unname(fit$BIC_values["k3", "R1"]), expected_grid_bic(models$R1, grid))
  expect_identical(first_seen_labels(fit$cluster_assignment[, "R1"]), rep(1:2, each = 3L))
  soft <- fit$replicate_ancestry_matrices$R1
  expect_equal(unname(soft[, "Cluster_3"]), expected_soft_weight(models$R1, 4L))
  expect_equal(unname(rowSums(soft)), rep(1, 6L))
})

test_that("all retained maps contribute to sample summaries in a mixed grid-K run", {
  models <- occupied_group_fixture()
  fit <- fit_occupied_group_models(models)

  expect_equal(as.integer(fit$optim_k_vals), c(3L, 2L))
  expect_identical(colnames(fit$BIC_values), names(models))
  expect_equal(unname(fit$BIC_values["k3", "R1"]),
               expected_grid_bic(models$R1, c(1, 1, 2, 3)))
  expect_equal(unname(fit$BIC_values["k2", "R2"]),
               expected_grid_bic(models$R2, c(1, 1, 2, 2)))
  expect_identical(colnames(fit$cluster_assignment), names(models))
  expect_identical(names(fit$replicate_ancestry_matrices), names(models))
  expect_identical(first_seen_labels(fit$cluster_assignment[, "R1"]), rep(1:2, each = 3L))
  expect_equal(unname(fit$cluster_assignment[, "R2"]), rep(1L, 6L))
  expect_equal(sort(unname(fit$ancestry_matrix[, "Cluster_2"])), rep(c(0, .5), each = 3L))
  expect_equal(unname(fit$replicate_ancestry_matrices$R2[, "Cluster_2"]),
               expected_soft_weight(models$R2, 3:4))
})

test_that("the first hard partition supplies a reference when profiles are insufficient", {
  specimens <- cbind(x = rep(.5, 6L), y = rep(c(0, 1), c(4L, 2L)))
  rownames(specimens) <- paste0("s", seq_len(nrow(specimens)))
  codes <- cbind(x = c(.25, .75, .5, .625, .875, .99),
                 y = c(0, 0, 1, 1, .625, .625))
  # The first four specimens are tied between units 1 and 2.
  first <- c(1L, 1L, 2L, 2L, 3L, 3L)
  second <- c(2L, 1L, 2L, 2L, 3L, 3L)
  models <- list(R1 = occupied_group_model(specimens, codes, first),
                 R2 = occupied_group_model(specimens, codes, second))
  expect_equal(nrow(unique(specimens)), 2L)
  expect_message(fit <- fit_occupied_group_models(models, set.k = 4L, max.k = 4L, verbose = TRUE),
                 "Using hard assignments from replicate R1 as the relabeling reference for 3 occupied groups",
                 fixed = TRUE)

  grid <- c(1L, 2L, 3L, 3L, 4L, 4L)
  expect_equal(as.integer(fit$optim_k_vals), c(4L, 4L))
  expect_equal(unname(fit$BIC_values["k4", ]), rep(expected_grid_bic(models$R1, grid), 2L))
  expect_identical(lapply(fit$som_clusters, first_seen_labels), list(R1 = grid, R2 = grid))
  expect_identical(unname(fit$cluster_assignment), unname(cbind(first, second)))
  expect_equal(unname(fit$replicate_ancestry_matrices$R1[, "Cluster_4"]),
               expected_soft_weight(models$R1, 5:6))
  expect_equal(unname(fit$replicate_ancestry_matrices$R2[, "Cluster_4"]),
               expected_soft_weight(models$R2, 5:6))
})
