model_plot_fixture <- function(mixed = FALSE) {
  x <- c(.02, .16, .32, .60, .77, .98)
  specimens <- cbind(x = x, y = 0)
  rownames(specimens) <- paste0("s", seq_len(nrow(specimens)))
  model <- function(codes, winners) {
    structure(list(data = list(specimens), codes = list(codes), whatmap = 1L,
                   dist.fcts = "sumofsquares", user.weights = 1, distance.weights = 1,
                   maxNA.fraction = .5, unit.classif = winners,
                   grid = kohonen::somgrid(3, 2, topo = "hexagonal")),
              class = "kohonen")
  }
  models <- list(R1 = model(cbind(x = c(0, .05, .95, .43, .52, 1),
                                  y = c(0, 0, 0, 1, 1, 0)),
                            c(1L, 2L, 2L, 3L, 3L, 6L)))
  grid_clusters <- list(R1 = c(2L, 2L, 4L, 1L, 3L, 4L))
  hard <- cbind(R1 = rep(c(2L, 1L), each = 3L))
  if (mixed) {
    models$R2 <- model(specimens, seq_len(nrow(specimens)))
    grid_clusters$R2 <- c(3L, 5L, 5L, 1L, 2L, 4L)
    hard <- cbind(hard, R2 = c(1L, 2L, 2L, 3L, 4L, 5L))
  }
  rownames(hard) <- rownames(specimens)
  ancestry <- vapply(seq_len(max(hard)), function(k) rowMeans(hard == k), numeric(nrow(hard)))
  list(som_models = models, som_clusters = grid_clusters, cluster_assignment = hard,
       ancestry_matrix = ancestry)
}

render_model_plot_pdf <- function(fit, ...) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path, width = 7, height = 7, compress = FALSE, useKerning = FALSE)
  device <- grDevices::dev.cur()
  on.exit({
    if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
    unlink(path)
  })
  expect_null(delimSOM::plot.model.SOM(fit, verbose = FALSE, ...))
  grDevices::dev.off(device)
  expect_gt(file.info(path)$size, 100)
  readLines(path, warn = FALSE, skipNul = TRUE)
}

has_unoccupied_plot_key <- function(pdf_lines) {
  any(grepl("No hard specimen assignments", pdf_lines, fixed = TRUE, useBytes = TRUE))
}

test_that("default plotting shows unoccupied grid clusters", {
  fit <- model_plot_fixture()
  before <- serialize(fit, NULL)
  expect_true(has_unoccupied_plot_key(render_model_plot_pdf(fit)))
  expect_identical(serialize(fit, NULL), before)
})

test_that("fixed-K plotting selects by grid K with unoccupied clusters", {
  fit <- model_plot_fixture()
  expect_true(has_unoccupied_plot_key(render_model_plot_pdf(fit, set.k = 4L)))
})

test_that("mixed-run plotting matches the selected occupied groups to the common palette", {
  fit <- model_plot_fixture(mixed = TRUE)
  expect_true(has_unoccupied_plot_key(render_model_plot_pdf(fit, replicate.mode = "first")))
  expect_false(has_unoccupied_plot_key(render_model_plot_pdf(fit, set.k = 5L)))
})
