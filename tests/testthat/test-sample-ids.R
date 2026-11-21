sample_id_fixture <- function() {
  matrix(sin(seq_len(48) * sqrt(2)), nrow = 12L,
         dimnames = list(sprintf("sample_%d", 1:12), paste0("v", 1:4)))
}

train_for_sample_id_test <- function(input, verbose = FALSE) {
  delimSOM::train.SOM(
    input, N.steps = 12L, N.replicates = 2L, grid.size = c(2L, 2L),
    parallel = FALSE, verbose = verbose
  )
}

cluster_for_sample_id_test <- function(trained, k, verbose = FALSE) {
  delimSOM::clustering.SOM(
    trained, max.k = 3L, set.k = k, clustering.method = "kmeans+BICelbow",
    quantization.error.quantile = NULL, topographic.error.quantile = NULL,
    calculate.variable.importance = FALSE, verbose = verbose
  )
}

test_that("malformed sample IDs are rejected before matching layers", {
  named <- sample_id_fixture()
  duplicate <- named
  rownames(duplicate)[1L] <- rownames(named)[2L]
  expect_error(train_for_sample_id_test(list(first = named, second = duplicate)),
               "Layer 2 (second) has duplicate sample IDs at rows 1, 2", fixed = TRUE)
  expect_error(train_for_sample_id_test(duplicate),
               "Single input layer has duplicate sample IDs at rows 1, 2", fixed = TRUE)

  for (id in c("", " \t", NA_character_)) {
    invalid <- named
    rownames(invalid)[1L] <- id
    expect_error(train_for_sample_id_test(list(first = named, second = invalid)),
                 "Layer 2 (second) has missing or blank sample IDs at rows 1", fixed = TRUE)
    expect_error(train_for_sample_id_test(invalid),
                 "Single input layer has missing or blank sample IDs at rows 1", fixed = TRUE)
  }

  unnamed <- named
  rownames(unnamed) <- NULL
  numeric_ids <- named
  rownames(numeric_ids) <- as.character(seq_len(nrow(numeric_ids)))
  for (layer in list(unnamed, as.data.frame(unnamed))) {
    expect_error(train_for_sample_id_test(list(first = numeric_ids, second = layer)),
                 "Layer 2 (second) must have explicit sample IDs as row names", fixed = TRUE)
  }
})

test_that("explicit IDs match across layer order and partial overlap", {
  first <- as.data.frame(sample_id_fixture())
  rownames(first) <- as.character(seq_len(nrow(first)))
  second <- first[12:3, , drop = FALSE]
  expected_ids <- as.character(3:12)
  trained <- train_for_sample_id_test(list(first = first, second = second))
  ordered <- train_for_sample_id_test(list(first = first[expected_ids, , drop = FALSE],
                                          second = second[expected_ids, , drop = FALSE]))

  expect_identical(lapply(trained$input_data, rownames),
                   list(first = expected_ids, second = expected_ids))
  expect_identical(trained$som_models, ordered$som_models)
})

test_that("single-layer IDs preserve input positions after filtering", {
  named <- sample_id_fixture()
  named[3L, ] <- NA_real_
  expected_ids <- rownames(named)[-3L]
  unnamed <- named
  rownames(unnamed) <- NULL
  named_models <- train_for_sample_id_test(named)$som_models
  for (input in list(unnamed, as.data.frame(unnamed), list(data = unnamed))) {
    trained <- train_for_sample_id_test(input)
    expect_identical(rownames(trained$input_data[[1L]]), expected_ids)
    for (id in names(trained$som_models)) {
      model <- trained$som_models[[id]]
      expect_identical(rownames(model$data[[1L]]), expected_ids)
      expect_identical(model$codes[[1L]], named_models[[id]]$codes[[1L]])
      expect_identical(model$unit.classif, named_models[[id]]$unit.classif)
    }
    for (k in 1:2) {
      clustered <- cluster_for_sample_id_test(trained, k)
      expect_identical(rownames(clustered$cluster_assignment), expected_ids)
      expect_identical(rownames(clustered$ancestry_matrix), expected_ids)
      for (id in clustered$retained_replicates) {
        expect_identical(rownames(clustered$cluster_gridcell_assignments[[id]]), expected_ids)
        soft <- clustered$replicate_ancestry_matrices[[id]]
        expect_identical(rownames(soft), expected_ids)
        expect_equal(unname(rowSums(soft)), rep(1, length(expected_ids)))
        if (k == 1L) expect_equal(as.numeric(soft), rep(1, length(expected_ids)))
      }
    }
  }
})

test_that("older unnamed fits receive IDs only when replicate rows align", {
  named <- train_for_sample_id_test(sample_id_fixture())
  older <- named
  rownames(older$input_data[[1L]]) <- NULL
  for (i in seq_along(older$som_models)) rownames(older$som_models[[i]]$data[[1L]]) <- NULL

  for (k in 1:2) {
    expected <- cluster_for_sample_id_test(named, k)
    clustered <- cluster_for_sample_id_test(older, k)
    expect_identical(clustered$som_models, expected$som_models)
    expect_identical(clustered$input_data, expected$input_data)
    expect_identical(clustered$cluster_assignment, expected$cluster_assignment)
    expect_identical(clustered$replicate_ancestry_matrices, expected$replicate_ancestry_matrices)
  }

  reordered <- older
  reorder <- rev(seq_len(nrow(older$input_data[[1L]])))
  reordered$som_models[[2L]]$data[[1L]] <- older$som_models[[2L]]$data[[1L]][reorder, , drop = FALSE]
  reordered$som_models[[2L]]$unit.classif <- older$som_models[[2L]]$unit.classif[reorder]
  expect_error(cluster_for_sample_id_test(reordered, 2L), "cannot align unnamed samples")
})

test_that("positional sample IDs are reported at default verbosity", {
  named <- sample_id_fixture()
  unnamed <- named
  rownames(unnamed) <- NULL
  expect_message(train_for_sample_id_test(unnamed, verbose = TRUE),
                 "Assigning IDs such as sample_1 from input row positions", fixed = TRUE)

  older <- train_for_sample_id_test(named)
  for (i in seq_along(older$som_models)) rownames(older$som_models[[i]]$data[[1L]]) <- NULL
  expect_message(cluster_for_sample_id_test(older, 1L, verbose = TRUE),
                 "Assigning IDs such as sample_1 from the stored row positions", fixed = TRUE)
})
