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

nested_split_model_plot_fixture <- function() {
  template <- model_plot_fixture()$som_models$R1
  specimens <- cbind(x = seq_len(17L) / 17, y = 0)
  rownames(specimens) <- sprintf("s%02d", seq_len(nrow(specimens)))
  model <- function(winners) {
    template$data <- list(specimens)
    template$unit.classif <- winners
    template
  }
  models <- list(R1 = model(rep(c(1L, 4L), c(12L, 5L))),
                 R2 = model(rep(c(1L, 3L, 5L), c(7L, 5L, 5L))),
                 R3 = model(rep(1:4, c(4L, 3L, 5L, 5L))))
  grid_clusters <- list(R1 = rep(1:2, each = 3L),
                        R2 = rep(c(1L, 3L, 2L), each = 2L),
                        R3 = c(1L, 4L, 3L, 2L, 2L, 3L))
  hard <- vapply(names(models), function(id) grid_clusters[[id]][models[[id]]$unit.classif],
                 integer(nrow(specimens)))
  rownames(hard) <- rownames(specimens)
  ancestry <- vapply(1:4, function(k) rowMeans(hard == k), numeric(nrow(hard)))
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

expect_rendered_model_cell_colors <- function(pdf_lines, colors) {
  # The lower panel draws one closed, filled hexagon per unit. Points and the
  # rectangular legend key use different PDF drawing commands.
  hexagons <- tail(which(grepl("^h B$", pdf_lines, useBytes = TRUE)), length(colors))
  fill_commands <- which(grepl(" (scn|rg)$", pdf_lines, useBytes = TRUE))
  cell_fills <- pdf_lines[fill_commands[findInterval(hexagons, fill_commands)]]
  rgb <- scan(text = sub(" (scn|rg)$", "", cell_fills), quiet = TRUE)
  expect_equal(matrix(rgb, nrow = 3L), round(unname(grDevices::col2rgb(colors)) / 255, 3))
}

test_that("default plotting shows unoccupied grid clusters", {
  fit <- model_plot_fixture()
  before <- serialize(fit, NULL)
  plot <- render_model_plot_pdf(fit)
  expect_true(has_unoccupied_plot_key(plot))
  palette <- viridis::viridis(2L)
  expect_rendered_model_cell_colors(plot, c(palette[c(2, 2, 1)], "#D9D9D9", "#D9D9D9", palette[1]))
  expect_identical(serialize(fit, NULL), before)
})

test_that("fixed-K plotting selects by grid K with unoccupied clusters", {
  fit <- model_plot_fixture()
  plot <- render_model_plot_pdf(fit, set.k = 4L)
  expect_true(has_unoccupied_plot_key(plot))
  palette <- viridis::viridis(4L)
  expect_rendered_model_cell_colors(plot, c(palette[c(2, 2, 1)], "#D9D9D9", "#D9D9D9", palette[1]))
})

test_that("mixed-run plotting matches the selected occupied groups to the common palette", {
  fit <- model_plot_fixture(mixed = TRUE)
  expect_true(has_unoccupied_plot_key(render_model_plot_pdf(fit, replicate.mode = "first")))
  plot <- render_model_plot_pdf(fit, set.k = 5L)
  expect_false(has_unoccupied_plot_key(plot))
  expect_rendered_model_cell_colors(plot, viridis::viridis(5L)[c(1, 2, 2, 3, 4, 5)])
})

test_that("a lower-K map displays stored labels through nested splits", {
  fit <- nested_split_model_plot_fixture()
  plot <- render_model_plot_pdf(fit, replicate.mode = "first")
  expect_rendered_model_cell_colors(plot, viridis::viridis(4L)[rep(1:2, each = 3L)])
})

test_that("unused palette columns in a saved fit do not prevent rendering", {
  fit <- nested_split_model_plot_fixture()
  fit$ancestry_matrix <- cbind(fit$ancestry_matrix, 0, 0)
  plot <- render_model_plot_pdf(fit, replicate.mode = "first")
  expect_rendered_model_cell_colors(plot, viridis::viridis(6L)[rep(1:2, each = 3L)])
})
