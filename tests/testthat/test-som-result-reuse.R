with_temporary_results_dir <- function(code) {
  path <- tempfile()
  dir.create(path)
  previous_dir <- setwd(path)
  on.exit({
    setwd(previous_dir)
    unlink(path, recursive = TRUE)
  })
  force(code)
}

result_test_specimens <- matrix(
  sin(seq_len(48) * sqrt(2)), nrow = 12,
  dimnames = list(paste0("s", 1:12), paste0("v", 1:4))
)

train_result_test_som <- function(specimens = result_test_specimens, ...) {
  delimSOM::train.SOM(
    specimens, N.steps = 12, N.replicates = 1, grid.size = c(2, 2),
    parallel = FALSE, verbose = FALSE, ...
  )
}

cluster_result_test_som <- function(trained, k, ...) {
  delimSOM::clustering.SOM(
    trained, set.k = k, max.k = 3, clustering.method = "kmeans+BICelbow",
    quantization.error.quantile = NULL, topographic.error.quantile = NULL,
    calculate.variable.importance = FALSE, verbose = FALSE, ...
  )
}

test_that("training uses current specimens when saving and reuse are disabled", {
  with_temporary_results_dir({
    saved <- train_result_test_som(save.SOM.results = TRUE)
    default_path <- saved$train.SOM.args$save.SOM.results.name
    explicit_path <- "previous-training.Rdata"
    stopifnot(file.copy(default_path, explicit_path))
    original_checksums <- tools::md5sum(c(default_path, explicit_path))
    specimens <- result_test_specimens[1:10, , drop = FALSE]

    from_default_path <- train_result_test_som(specimens)
    from_explicit_path <- train_result_test_som(specimens, save.SOM.results.name = explicit_path)

    expect_identical(rownames(from_default_path$som_models[[1]]$data[[1]]), rownames(specimens))
    expect_identical(rownames(from_explicit_path$som_models[[1]]$data[[1]]), rownames(specimens))
    expect_identical(tools::md5sum(names(original_checksums)), original_checksums)
  })
})

test_that("clustering uses current K when saving and reuse are disabled", {
  with_temporary_results_dir({
    trained <- train_result_test_som()
    saved <- cluster_result_test_som(trained, 2L, save.SOM.results = TRUE)
    default_path <- saved$clustering.SOM.args$save.SOM.results.name
    explicit_path <- "previous-clustering.Rdata"
    stopifnot(file.copy(default_path, explicit_path))
    original_checksums <- tools::md5sum(c(default_path, explicit_path))

    from_default_path <- cluster_result_test_som(trained, 1L)
    from_explicit_path <- cluster_result_test_som(trained, 1L, save.SOM.results.name = explicit_path)

    expect_equal(unname(from_default_path$optim_k_vals[1, 1]), 1L)
    expect_equal(unname(from_explicit_path$optim_k_vals[1, 1]), 1L)
    expect_identical(tools::md5sum(names(original_checksums)), original_checksums)
  })
})

test_that("saving training and clustering results protects existing files", {
  with_temporary_results_dir({
    trained <- train_result_test_som(save.SOM.results = TRUE)
    training_path <- trained$train.SOM.args$save.SOM.results.name
    clustering_path <- "results/clustering.Rdata"
    cluster_result_test_som(trained, 2L, save.SOM.results = TRUE,
                            save.SOM.results.name = clustering_path)
    original_checksums <- tools::md5sum(c(training_path, clustering_path))
    specimens <- result_test_specimens[1:10, , drop = FALSE]

    expect_error(train_result_test_som(specimens, save.SOM.results = TRUE),
                 "SOM results file already exists")
    expect_error(cluster_result_test_som(trained, 1L, save.SOM.results = TRUE,
                                        save.SOM.results.name = clustering_path),
                 "SOM results file already exists")
    expect_identical(tools::md5sum(names(original_checksums)), original_checksums)

    new_training <- train_result_test_som(specimens, save.SOM.results = TRUE,
                                          overwrite.SOM.results = TRUE)
    new_clustering <- cluster_result_test_som(trained, 1L, save.SOM.results = TRUE,
                                              save.SOM.results.name = clustering_path,
                                              overwrite.SOM.results = TRUE)
    expect_identical(rownames(new_training$som_models[[1]]$data[[1]]), rownames(specimens))
    expect_equal(unname(new_clustering$optim_k_vals[1, 1]), 1L)
    expect_identical(delimSOM::train.SOM(save.SOM.results.name = training_path,
                                        reuse.SOM.results = TRUE, verbose = FALSE),
                     new_training)
    expect_identical(delimSOM::clustering.SOM(save.SOM.results.name = clustering_path,
                                             reuse.SOM.results = TRUE, verbose = FALSE),
                     new_clustering)
  })
})

test_that("explicit reuse validates the path and saved results", {
  with_temporary_results_dir({
    expect_error(delimSOM::train.SOM(reuse.SOM.results = TRUE), "explicit .Rdata file path",
                 fixed = TRUE)
    expect_error(delimSOM::clustering.SOM(reuse.SOM.results = TRUE,
                                         save.SOM.results.name = "previous.rds"),
                 "must end with '.Rdata'", fixed = TRUE)
    expect_error(delimSOM::train.SOM(reuse.SOM.results = TRUE,
                                    save.SOM.results.name = "missing.Rdata"),
                 "does not exist")
    expect_error(delimSOM::clustering.SOM(reuse.SOM.results = TRUE,
                                         save.SOM.results = TRUE),
                 "requires save.SOM.results and overwrite.SOM.results to be FALSE",
                 fixed = TRUE)

    SOM_results <- list(som_models = list())
    save(SOM_results, file = "malformed.Rdata")
    expect_error(delimSOM::train.SOM(reuse.SOM.results = TRUE,
                                    save.SOM.results.name = "malformed.Rdata",
                                    verbose = FALSE),
                 "does not contain valid SOM training results")
    expect_error(delimSOM::clustering.SOM(reuse.SOM.results = TRUE,
                                         save.SOM.results.name = "malformed.Rdata",
                                         verbose = FALSE),
                 "does not contain valid SOM clustering results")
  })
})
