test_that("leave-one-layer-out margins require two grid clusters on both maps", {
  group <- rep(c(0, 1), each = 6L)
  offset <- rep(seq(0, .05, length.out = 6L), 2L)
  specimens <- cbind(x = group + offset, y = group - offset)
  rownames(specimens) <- paste0("s", seq_len(nrow(specimens)))
  trained <- delimSOM::train.SOM(
    list(first = specimens, second = specimens[, 2:1]),
    N.steps = 40L, N.replicates = 2L, grid.size = c(2L, 2L),
    parallel = FALSE, verbose = FALSE, set.seed.N = 20L
  )
  # A square codebook gives R1 a reproducible K = 1 while the real
  # omission reruns recover the two specimen groups in each input layer.
  square <- cbind(x = c(0, 0, 1, 1), y = c(0, 1, 0, 1))
  trained$som_models$R1$codes[] <- list(square)
  trained$som_models$R1$unit.classif <- ifelse(group == 0, 1L, 4L)
  fit <- delimSOM::clustering.SOM(
    trained, max.k = 2L, clustering.method = "kmeans+BICthreshold", BIC.thresh = 3,
    quantization.error.quantile = NULL, topographic.error.quantile = NULL,
    calculate.variable.importance = FALSE, verbose = FALSE
  )
  expect_equal(as.integer(fit$optim_k_vals), c(1L, 2L))
  expect_equal(ncol(fit$replicate_ancestry_matrices$R1), 2L)

  plot_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_path, width = 12, height = 8)
  device <- grDevices::dev.cur()
  results <- tryCatch(
    delimSOM::plot.layer.importance.leaveoneout.SOM(
      fit, save.leave.one.layer.out.results = FALSE, add.points = FALSE, verbose = FALSE
    ),
    finally = {
      grDevices::dev.off(device)
      unlink(plot_path)
    }
  )
  comparisons <- results$replicate.matched.results
  expect_true(all(is.na(comparisons$error)))
  expect_equal(comparisons$leave.one.layer.out.modal.k, rep(2L, 4L))
  single <- comparisons[comparisons$baseline.modal.k == 1L, ]
  expect_equal(single$baseline.mean.assignment.margin, rep(NA_real_, 2L))
  expect_equal(single$delta.mean.assignment.margin, rep(NA_real_, 2L))
  comparable <- comparisons[comparisons$baseline.modal.k == 2L, ]
  expect_true(all(is.finite(comparable$delta.mean.assignment.margin)))
})
