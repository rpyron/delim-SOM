#### Supplementary Figure XX ##################################################


## Set directories
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package/Figure_files")
figure.dir <- getwd()
intermediate_files_dir <- "intermediate_files"
dir.create(intermediate_files_dir, showWarnings = FALSE, recursive = TRUE)


## Load delim-SOM functions
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")


## Set benchmark settings
set.seed(1)
N.samples <- 50
core.values <- 1:10
overwrite.benchmark.results <- FALSE
benchmark.csv <- file.path(intermediate_files_dir, "train_SOM_runtime_summary.csv")
benchmark.rds <- file.path(intermediate_files_dir, "train_SOM_runtime_summary.rds")


## Create function to make two-layer test data
make.test.SOM.data <- function(N.samples, N.variables.total, seed = 1) {
  set.seed(seed)
  N.variables.layer1 <- floor(N.variables.total / 2)
  N.variables.layer2 <- N.variables.total - N.variables.layer1
  sample.names <- paste0("sample_", seq_len(N.samples))
  layer1 <- matrix(sample(0:2, N.samples * N.variables.layer1, replace = TRUE), nrow = N.samples, ncol = N.variables.layer1)
  layer2 <- matrix(stats::rnorm(N.samples * N.variables.layer2), nrow = N.samples, ncol = N.variables.layer2)
  rownames(layer1) <- sample.names
  rownames(layer2) <- sample.names
  colnames(layer1) <- paste0("SNP_", seq_len(N.variables.layer1))
  colnames(layer2) <- paste0("continuous_", seq_len(N.variables.layer2))
  input.data <- list(
    SNPs = layer1,
    Continuous = layer2
  )
  return(input.data)
}


## Create or load benchmark datasets
data.file.1k <- file.path(intermediate_files_dir, "test_data_50_samples_1000_variables_2_layers.Rdata")
data.file.25k <- file.path(intermediate_files_dir, "test_data_50_samples_25000_variables_2_layers.Rdata")
if (!file.exists(data.file.1k)) {
  input.data.1k <- make.test.SOM.data(N.samples = N.samples, N.variables.total = 1000, seed = 1)
  save(input.data.1k, file = data.file.1k)
} else {
  load(data.file.1k)
}
if (!file.exists(data.file.25k)) {
  input.data.25k <- make.test.SOM.data(N.samples = N.samples, N.variables.total = 25000, seed = 2)
  save(input.data.25k, file = data.file.25k)
} else {
  load(data.file.25k)
}


## Create benchmark table
benchmark.results <- data.frame(
  dataset = character(0),
  N.samples = integer(0),
  N.variables.total = integer(0),
  run.type = character(0),
  parallel = logical(0),
  N.cores = integer(0),
  N.steps = integer(0),
  N.replicates = integer(0),
  grid.x = integer(0),
  grid.y = integer(0),
  status = character(0),
  elapsed.seconds = numeric(0),
  user.seconds = numeric(0),
  system.seconds = numeric(0),
  output.file = character(0),
  error.message = character(0),
  stringsAsFactors = FALSE
)
if (file.exists(benchmark.rds)) benchmark.results <- readRDS(benchmark.rds)


## Create function to save benchmark table after every run
save.benchmark.results <- function(benchmark.results) {
  saveRDS(benchmark.results, benchmark.rds)
  utils::write.csv(benchmark.results, benchmark.csv, row.names = FALSE)
}


## Create function to remove previous benchmark row
remove.previous.benchmark.result <- function(dataset.name, run.type) {
  keep.rows <- !(benchmark.results$dataset == dataset.name & benchmark.results$run.type == run.type)
  benchmark.results <<- benchmark.results[keep.rows, , drop = FALSE]
  save.benchmark.results(benchmark.results)
}


## Create function to run one benchmark
run.train.SOM.benchmark <- function(input.data, dataset.name, N.variables.total, run.type, parallel.run, N.cores.run) {
  output.file <- file.path(
    intermediate_files_dir,
    paste0("SOM_", dataset.name, "_", run.type, ".Rdata")
  )
  completed.row <- benchmark.results$dataset == dataset.name &
    benchmark.results$run.type == run.type &
    benchmark.results$status == "completed"
  already.done <- any(completed.row) && file.exists(output.file)
  if (already.done && !overwrite.benchmark.results) {
    message("Skipping completed run: ", dataset.name, " / ", run.type)
    return(invisible(NULL))
  }
  if (overwrite.benchmark.results) {
    remove.previous.benchmark.result(dataset.name = dataset.name, run.type = run.type)
    if (file.exists(output.file)) unlink(output.file)
  }
  message("Running: ", dataset.name, " / ", run.type)
  gc()
  run.time <- NULL
  error.message <- NA_character_
  status <- "completed"
  run.time <- system.time({
    SOM.output <- tryCatch(
      train.SOM(
        input_data = input.data,
        parallel = parallel.run,
        N.cores = N.cores.run,
        learning.rate.tuning = FALSE,
        save.SOM.results = FALSE,
        overwrite.SOM.results = TRUE,
        verbose = TRUE,
        message.N.replicates = 20,
        set.seed.N = 1
      ),
      error = function(e) {
        status <<- "failed"
        error.message <<- conditionMessage(e)
        return(NULL)
      }
    )
  })
  if (!is.null(SOM.output)) {
    save(SOM.output, file = output.file)
    N.steps.used <- SOM.output$train.SOM.args$N.steps
    N.replicates.used <- SOM.output$train.SOM.args$N.replicates
    grid.x.used <- SOM.output$som_models[[1]]$grid$xdim
    grid.y.used <- SOM.output$som_models[[1]]$grid$ydim
  } else {
    N.steps.used <- NA_integer_
    N.replicates.used <- NA_integer_
    grid.x.used <- NA_integer_
    grid.y.used <- NA_integer_
  }
  new.result <- data.frame(
    dataset = dataset.name,
    N.samples = N.samples,
    N.variables.total = N.variables.total,
    run.type = run.type,
    parallel = parallel.run,
    N.cores = N.cores.run,
    N.steps = N.steps.used,
    N.replicates = N.replicates.used,
    grid.x = grid.x.used,
    grid.y = grid.y.used,
    status = status,
    elapsed.seconds = unname(run.time[["elapsed"]]),
    user.seconds = unname(run.time[["user.self"]]),
    system.seconds = unname(run.time[["sys.self"]]),
    output.file = ifelse(is.null(SOM.output), NA_character_, output.file),
    error.message = error.message,
    stringsAsFactors = FALSE
  )
  benchmark.results <<- rbind(benchmark.results, new.result)
  save.benchmark.results(benchmark.results)
  rm(SOM.output)
  gc()
}


## Run 1k-variable tests
run.train.SOM.benchmark(
  input.data = input.data.1k,
  dataset.name = "50_samples_1000_variables_2_layers",
  N.variables.total = 1000,
  run.type = "no_parallel",
  parallel.run = FALSE,
  N.cores.run = NA_integer_
)
for (core.i in core.values) {
  run.train.SOM.benchmark(
    input.data = input.data.1k,
    dataset.name = "50_samples_1000_variables_2_layers",
    N.variables.total = 1000,
    run.type = paste0("parallel_", sprintf("%02d", core.i), "_cores"),
    parallel.run = TRUE,
    N.cores.run = core.i
  )
}


## Run 25k-variable tests
run.train.SOM.benchmark(
  input.data = input.data.25k,
  dataset.name = "50_samples_25000_variables_2_layers",
  N.variables.total = 25000,
  run.type = "no_parallel",
  parallel.run = FALSE,
  N.cores.run = NA_integer_
)
for (core.i in core.values) {
  run.train.SOM.benchmark(
    input.data = input.data.25k,
    dataset.name = "50_samples_25000_variables_2_layers",
    N.variables.total = 25000,
    run.type = paste0("parallel_", sprintf("%02d", core.i), "_cores"),
    parallel.run = TRUE,
    N.cores.run = core.i
  )
}


## Load final benchmark summary
benchmark.results <- readRDS(benchmark.rds)
print(benchmark.results)


## Plot train.SOM benchmark results
benchmark.plot <- benchmark.results
benchmark.plot$core.plot <- ifelse(benchmark.plot$parallel, benchmark.plot$N.cores, 0)
benchmark.plot$core.label <- ifelse(benchmark.plot$parallel, as.character(benchmark.plot$N.cores), "no parallel")
benchmark.plot <- benchmark.plot[benchmark.plot$status == "completed", , drop = FALSE]
runtime.plot <- ggplot2::ggplot(benchmark.plot, ggplot2::aes(x = core.plot, y = elapsed.seconds, group = dataset)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~ dataset, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = 0:10, labels = c("no\nparallel", as.character(1:10))) +
  ggplot2::labs(x = "Number of cores", y = "Elapsed time (seconds)") +
  ggplot2::theme_classic()
ggplot2::ggsave(file.path(figure.dir, "Supplementary_figure_SXX.svg"), runtime.plot, width = 8, height = 5)
print(runtime.plot)

