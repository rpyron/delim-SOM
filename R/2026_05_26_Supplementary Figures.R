#### Supplementary Figure XX ##################################################


## Set directories
setwd("C:/Users/daniel/Desktop/SOM package/Figure_files")
figure.dir <- getwd()
intermediate_files_dir <- "intermediate_files"
dir.create(intermediate_files_dir, showWarnings = FALSE, recursive = TRUE)


## Load delim-SOM functions
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")


## Set benchmark settings
set.seed(1)
N.samples <- 50
core.values <- 1:10
overwrite.benchmark.results <- FALSE
benchmark.csv <- file.path(intermediate_files_dir, "train_SOM_runtime_summary.csv")
benchmark.rds <- file.path(intermediate_files_dir, "train_SOM_runtime_summary.rds")


## Set benchmark dataset settings
benchmark.datasets <- data.frame(
  dataset.name = c("50_samples_1000_variables_2_layers",
                   "50_samples_10000_variables_2_layers",
                   "50_samples_25000_variables_2_layers",
                   "50_samples_60000_variables_2_layers"),
  N.variables.total = c(1000, 10000, 25000, 60000),
  seed = c(1, 2, 3, 4),
  stringsAsFactors = FALSE
)


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
benchmark.input.data <- list()
for (dataset.i in seq_len(nrow(benchmark.datasets))) {
  dataset.name.i <- benchmark.datasets$dataset.name[dataset.i]
  N.variables.total.i <- benchmark.datasets$N.variables.total[dataset.i]
  seed.i <- benchmark.datasets$seed[dataset.i]
  data.file.i <- file.path(intermediate_files_dir, paste0("test_data_", dataset.name.i, ".Rdata"))
  input.object.name.i <- paste0("input.data.", N.variables.total.i)
  if (!file.exists(data.file.i)) {
    input.data.i <- make.test.SOM.data(N.samples = N.samples, N.variables.total = N.variables.total.i, seed = seed.i)
    save(input.data.i, file = data.file.i)
  } else {
    load(data.file.i)
  }
  benchmark.input.data[[dataset.name.i]] <- input.data.i
  assign(input.object.name.i, input.data.i)
  rm(input.data.i)
  gc()
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
  SOM.output <- NULL
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


## Run train.SOM benchmark tests
for (dataset.i in seq_len(nrow(benchmark.datasets))) {
  dataset.name.i <- benchmark.datasets$dataset.name[dataset.i]
  N.variables.total.i <- benchmark.datasets$N.variables.total[dataset.i]
  input.data.i <- benchmark.input.data[[dataset.name.i]]
  run.train.SOM.benchmark(
    input.data = input.data.i,
    dataset.name = dataset.name.i,
    N.variables.total = N.variables.total.i,
    run.type = "no_parallel",
    parallel.run = FALSE,
    N.cores.run = NA_integer_
  )
  for (core.i in core.values) {
    run.train.SOM.benchmark(
      input.data = input.data.i,
      dataset.name = dataset.name.i,
      N.variables.total = N.variables.total.i,
      run.type = paste0("parallel_", sprintf("%02d", core.i), "_cores"),
      parallel.run = TRUE,
      N.cores.run = core.i
    )
  }
  rm(input.data.i)
  gc()
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
ggplot2::ggsave(file.path(figure.dir, "Supplementary_figure_SXX_train_SOM.svg"), runtime.plot, width = 8, height = 6)
print(runtime.plot)


## Set clustering benchmark settings
clustering.method <- "kmeans+BICelbow"
clustering.core.values <- core.values
trained.SOM.run.type.for.clustering <- "parallel_03_cores"
calculate.soft.ancestry.clustering.benchmark <- FALSE
save.clustering.outputs <- FALSE
overwrite.clustering.benchmark.results <- FALSE
clustering.benchmark.csv <- file.path(intermediate_files_dir, "clustering_SOM_kmeansBICelbow_runtime_summary.csv")
clustering.benchmark.rds <- file.path(intermediate_files_dir, "clustering_SOM_kmeansBICelbow_runtime_summary.rds")


## Create clustering benchmark table
clustering.benchmark.results <- data.frame(
  dataset = character(0),
  N.samples = integer(0),
  N.variables.total = integer(0),
  trained.SOM.run.type = character(0),
  clustering.method = character(0),
  run.type = character(0),
  parallel = logical(0),
  N.cores = integer(0),
  N.replicates.input = integer(0),
  N.replicates.clustered = integer(0),
  status = character(0),
  elapsed.seconds = numeric(0),
  user.seconds = numeric(0),
  system.seconds = numeric(0),
  input.SOM.file = character(0),
  output.file = character(0),
  optim.k.summary = character(0),
  error.message = character(0),
  stringsAsFactors = FALSE
)
if (file.exists(clustering.benchmark.rds)) clustering.benchmark.results <- readRDS(clustering.benchmark.rds)


## Create function to save clustering benchmark table after every run
save.clustering.benchmark.results <- function(clustering.benchmark.results) {
  saveRDS(clustering.benchmark.results, clustering.benchmark.rds)
  utils::write.csv(clustering.benchmark.results, clustering.benchmark.csv, row.names = FALSE)
}


## Create function to remove previous clustering benchmark row
remove.previous.clustering.benchmark.result <- function(dataset.name, trained.SOM.run.type, run.type) {
  keep.rows <- !(clustering.benchmark.results$dataset == dataset.name &
                   clustering.benchmark.results$trained.SOM.run.type == trained.SOM.run.type &
                   clustering.benchmark.results$run.type == run.type)
  clustering.benchmark.results <<- clustering.benchmark.results[keep.rows, , drop = FALSE]
  save.clustering.benchmark.results(clustering.benchmark.results)
}


## Create function to summarize optimal K
summarize.optim.k <- function(clustering.output) {
  if (is.null(clustering.output$optim_k_summary)) return(NA_character_)
  optim.k.summary <- clustering.output$optim_k_summary
  if (is.data.frame(optim.k.summary)) {
    if (nrow(optim.k.summary) == 0) return(NA_character_)
    if (all(c("Count", "Proportion") %in% colnames(optim.k.summary))) {
      return(paste0(rownames(optim.k.summary), "=", optim.k.summary$Count, " (", round(optim.k.summary$Proportion, 3), ")", collapse = "; "))
    }
    return(paste(capture.output(print(optim.k.summary)), collapse = " | "))
  }
  if (is.matrix(optim.k.summary)) {
    if (nrow(optim.k.summary) == 0) return(NA_character_)
    if (all(c("Count", "Proportion") %in% colnames(optim.k.summary))) {
      return(paste0(rownames(optim.k.summary), "=", optim.k.summary[, "Count"], " (", round(optim.k.summary[, "Proportion"], 3), ")", collapse = "; "))
    }
    return(paste(capture.output(print(optim.k.summary)), collapse = " | "))
  }
  if (is.atomic(optim.k.summary)) {
    if (length(optim.k.summary) == 0) return(NA_character_)
    if (!is.null(names(optim.k.summary))) return(paste0(names(optim.k.summary), "=", optim.k.summary, collapse = "; "))
    return(paste(optim.k.summary, collapse = "; "))
  }
  return(paste(capture.output(print(optim.k.summary)), collapse = " | "))
}


## Create function to count clustered replicates
count.clustered.replicates <- function(clustering.output) {
  if (!is.null(clustering.output$som_models)) return(length(clustering.output$som_models))
  if (!is.null(clustering.output$replicate_ids)) return(length(clustering.output$replicate_ids))
  if (!is.null(clustering.output$sample_assignments) && length(clustering.output$sample_assignments) > 0) return(length(clustering.output$sample_assignments))
  return(NA_integer_)
}


## Create function to run one clustering benchmark
run.clustering.SOM.benchmark <- function(SOM.file, dataset.name, N.variables.total, trained.SOM.run.type, run.type, parallel.run, N.cores.run) {
  output.file <- file.path(
    intermediate_files_dir,
    paste0("clustering_", dataset.name, "_", trained.SOM.run.type, "_kmeansBICelbow_", run.type, ".Rdata")
  )
  completed.row <- clustering.benchmark.results$dataset == dataset.name &
    clustering.benchmark.results$trained.SOM.run.type == trained.SOM.run.type &
    clustering.benchmark.results$run.type == run.type &
    clustering.benchmark.results$status == "completed"
  already.done <- any(completed.row) && (!save.clustering.outputs || file.exists(output.file))
  if (already.done && !overwrite.clustering.benchmark.results) {
    message("Skipping completed clustering run: ", dataset.name, " / ", trained.SOM.run.type, " / ", run.type)
    return(invisible(NULL))
  }
  if (overwrite.clustering.benchmark.results) {
    remove.previous.clustering.benchmark.result(dataset.name = dataset.name, trained.SOM.run.type = trained.SOM.run.type, run.type = run.type)
    if (file.exists(output.file)) unlink(output.file)
  }
  message("Running clustering: ", dataset.name, " / ", trained.SOM.run.type, " / ", run.type)
  gc()
  load(SOM.file)
  if (!exists("SOM.output")) stop("Clustering benchmark aborted: SOM.output not found in ", SOM.file)
  run.time <- NULL
  error.message <- NA_character_
  status <- "completed"
  clustering.output <- NULL
  run.time <- system.time({
    clustering.output <- tryCatch(
      clustering.SOM(
        SOM.output,
        clustering.method = clustering.method,
        parallel = parallel.run,
        N.cores = N.cores.run,
        calculate.soft.ancestry = calculate.soft.ancestry.clustering.benchmark,
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
  if (!is.null(clustering.output)) {
    if (save.clustering.outputs) save(clustering.output, file = output.file)
    N.replicates.input <- length(SOM.output$som_models)
    N.replicates.clustered <- count.clustered.replicates(clustering.output)
    optim.k.summary <- summarize.optim.k(clustering.output)
  } else {
    N.replicates.input <- length(SOM.output$som_models)
    N.replicates.clustered <- NA_integer_
    optim.k.summary <- NA_character_
  }
  new.result <- data.frame(
    dataset = dataset.name,
    N.samples = N.samples,
    N.variables.total = N.variables.total,
    trained.SOM.run.type = trained.SOM.run.type,
    clustering.method = clustering.method,
    run.type = run.type,
    parallel = parallel.run,
    N.cores = N.cores.run,
    N.replicates.input = N.replicates.input,
    N.replicates.clustered = N.replicates.clustered,
    status = status,
    elapsed.seconds = unname(run.time[["elapsed"]]),
    user.seconds = unname(run.time[["user.self"]]),
    system.seconds = unname(run.time[["sys.self"]]),
    input.SOM.file = SOM.file,
    output.file = ifelse(is.null(clustering.output) || !save.clustering.outputs, NA_character_, output.file),
    optim.k.summary = optim.k.summary,
    error.message = error.message,
    stringsAsFactors = FALSE
  )
  clustering.benchmark.results <<- rbind(clustering.benchmark.results, new.result)
  save.clustering.benchmark.results(clustering.benchmark.results)
  rm(SOM.output, clustering.output)
  gc()
}


## Select trained SOM outputs for clustering benchmark
benchmark.results <- readRDS(benchmark.rds)
trained.SOM.files.for.clustering <- benchmark.results[
  benchmark.results$status == "completed" &
    benchmark.results$run.type == trained.SOM.run.type.for.clustering &
    !is.na(benchmark.results$output.file) &
    file.exists(benchmark.results$output.file),
  , drop = FALSE
]


## Run clustering benchmarks with trained SOM outputs
for (row.i in seq_len(nrow(trained.SOM.files.for.clustering))) {
  SOM.file.i <- trained.SOM.files.for.clustering$output.file[row.i]
  dataset.name.i <- trained.SOM.files.for.clustering$dataset[row.i]
  N.variables.total.i <- trained.SOM.files.for.clustering$N.variables.total[row.i]
  trained.SOM.run.type.i <- trained.SOM.files.for.clustering$run.type[row.i]
  run.clustering.SOM.benchmark(
    SOM.file = SOM.file.i,
    dataset.name = dataset.name.i,
    N.variables.total = N.variables.total.i,
    trained.SOM.run.type = trained.SOM.run.type.i,
    run.type = "no_parallel",
    parallel.run = FALSE,
    N.cores.run = NA_integer_
  )
  for (core.i in clustering.core.values) {
    run.clustering.SOM.benchmark(
      SOM.file = SOM.file.i,
      dataset.name = dataset.name.i,
      N.variables.total = N.variables.total.i,
      trained.SOM.run.type = trained.SOM.run.type.i,
      run.type = paste0("parallel_", sprintf("%02d", core.i), "_cores"),
      parallel.run = TRUE,
      N.cores.run = core.i
    )
  }
}


## Load final clustering benchmark summary
clustering.benchmark.results <- readRDS(clustering.benchmark.rds)
print(clustering.benchmark.results)


## Plot clustering.SOM benchmark results
clustering.benchmark.plot <- clustering.benchmark.results
clustering.benchmark.plot$core.plot <- ifelse(clustering.benchmark.plot$parallel, clustering.benchmark.plot$N.cores, 0)
clustering.benchmark.plot$core.label <- ifelse(clustering.benchmark.plot$parallel, as.character(clustering.benchmark.plot$N.cores), "no parallel")
clustering.benchmark.plot <- clustering.benchmark.plot[clustering.benchmark.plot$status == "completed", , drop = FALSE]
clustering.runtime.plot <- ggplot2::ggplot(clustering.benchmark.plot, ggplot2::aes(x = core.plot, y = elapsed.seconds, group = dataset)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~ dataset, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = 0:10, labels = c("no\nparallel", as.character(1:10))) +
  ggplot2::labs(x = "Number of cores", y = "Elapsed time (seconds)") +
  ggplot2::theme_classic()
ggplot2::ggsave(file.path(figure.dir, "Supplementary_figure_SXX_clustering_SOM_kmeansBICelbow.svg"), clustering.runtime.plot, width = 8, height = 6)
print(clustering.runtime.plot)
