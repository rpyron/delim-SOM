run_delimSOM_tests <-
function ()
{
    testthat::test_that("All exported delimSOM functions work", {
        test_failures <- character()
        tested_functions <- character()
        run_test <- function(test_name, test_expression) {
            function_name <- strsplit(test_name, " ", fixed = TRUE)[[1]][1]
            tested_functions <<- unique(c(tested_functions, function_name))
            tryCatch({
                test_result <- force(test_expression)
                message("[PASS] ", test_name)
                test_result
            }, error = function(error) {
                test_failures <<- c(test_failures, paste0(test_name, ": ", conditionMessage(error)))
                message("[FAIL] ", test_name, ": ", conditionMessage(error))
                NULL
            })
        }
        test_plot <- function(plot_expression) {
            plot_path <- tempfile(fileext = ".pdf")
            grDevices::pdf(plot_path, width = 12, height = 9)
            plot_device <- grDevices::dev.cur()
            tryCatch({
                force(plot_expression)
            }, finally = {
                if (plot_device %in% grDevices::dev.list()) {
                  grDevices::dev.off(which = plot_device)
                }
            })
            stopifnot(file.exists(plot_path))
            stopifnot(file.info(plot_path)$size > 100)
            invisible(TRUE)
        }
        set.seed(123)
        N_samples <- 60
        sample_groups <- rep(1:3, each = 20)
        sample_names <- paste0("Individual_", seq_len(N_samples))
        SNP_data <- replicate(24, rbinom(N_samples, 2, c(0.15, 0.5, 0.85)[sample_groups]))
        morphology_data <- cbind(Wing_length = rnorm(N_samples, sample_groups * 2, 0.5), Body_length = rnorm(N_samples, sample_groups, 0.5), Wing_width = rnorm(N_samples, sample_groups * 1.5, 0.5), Antenna_length = rnorm(N_samples, sample_groups * 0.8, 0.3))
        environment_data <- cbind(Temperature = rnorm(N_samples, sample_groups * 3, 1), Precipitation = rnorm(N_samples, sample_groups * 5, 2), Elevation = rnorm(N_samples, sample_groups * 50, 10))
        rownames(SNP_data) <- sample_names
        rownames(morphology_data) <- sample_names
        rownames(environment_data) <- sample_names
        colnames(SNP_data) <- paste0("SNP_", seq_len(ncol(SNP_data)))
        input_data_multi <- list(SNPs = SNP_data, Morphology = morphology_data, Environment = environment_data)
        Coordinates <- data.frame(Latitude = runif(N_samples, 36, 39), Longitude = runif(N_samples, -89, -82), row.names = sample_names)
        trait_data <- data.frame(Color = rep(c("Red", "Blue", "Green"), length.out = N_samples), Habitat = rep(c("Forest", "Grassland"), length.out = N_samples), Wing_length = morphology_data[, "Wing_length"], row.names = sample_names)
        run_test("make.cols.binary.SOM", {
            test_result <- delimSOM::make.cols.binary.SOM(dataframe = trait_data, make.binary.cols = c("Color", "Habitat"), append.to.original = TRUE, remove.original.cols = TRUE, verbose = FALSE)
            stopifnot(nrow(test_result) == N_samples)
            stopifnot(!any(c("Color", "Habitat") %in% colnames(test_result)))
            test_result
        })
        run_test("remove.lowCV.multicollinearity.SOM", {
            filter_input <- as.data.frame(morphology_data)
            filter_input$Constant_variable <- 1
            filter_input$Duplicate_wing_length <- filter_input$Wing_length
            test_result <- delimSOM::remove.lowCV.multicollinearity.SOM(input.dataframe = filter_input, CV.threshold = 0.05, cor.threshold = 0.9, prevalence.threshold = 0.05, verbose = FALSE)
            stopifnot(nrow(test_result) == N_samples)
            stopifnot(!("Constant_variable" %in% colnames(test_result)))
            test_result
        })
        run_test("process.SNP.data.SOM", {
            test_result <- delimSOM::process.SNP.data.SOM(snp.matrix.input = SNP_data, singleton.loci.filter = FALSE, verbose = FALSE)
            stopifnot(is.data.frame(test_result))
            stopifnot(nrow(test_result) == N_samples)
            stopifnot(ncol(test_result) > 0)
            stopifnot(all(as.matrix(test_result) %in% 0:2))
            test_result
        })
        SOM_single <- run_test("train.SOM single-layer", {
            test_result <- delimSOM::train.SOM(input_data = morphology_data, N.steps = 12, N.replicates = 4, grid.size = c(4, 4), parallel = FALSE, save.SOM.results = FALSE, overwrite.SOM.results = TRUE, verbose = FALSE, set.seed.N = 123)
            stopifnot(length(test_result$som_models) == 4)
            stopifnot(all(is.finite(test_result$quantization_error)))
            test_result
        })
        SOM_multi <- run_test("train.SOM multilayer", {
            test_result <- delimSOM::train.SOM(input_data = input_data_multi, N.steps = 12, N.replicates = 4, grid.size = c(4, 4), parallel = FALSE, save.SOM.results = FALSE, overwrite.SOM.results = TRUE, verbose = FALSE, set.seed.N = 123)
            stopifnot(length(test_result$som_models) == 4)
            stopifnot(all(is.finite(test_result$quantization_error)))
            test_result
        })
        SOM_parallel <- run_test("train.SOM parallel", {
            test_result <- delimSOM::train.SOM(input_data = morphology_data, N.steps = 12, N.replicates = 4, grid.size = c(4, 4), parallel = TRUE, N.cores = 2, save.SOM.results = FALSE, overwrite.SOM.results = TRUE, verbose = FALSE, set.seed.N = 123)
            stopifnot(length(test_result$som_models) == 4)
            stopifnot(all(is.finite(test_result$quantization_error)))
            test_result
        })
        run_test("train.SOM reproducibility", {
            stopifnot(!is.null(SOM_single))
            test_result <- delimSOM::train.SOM(input_data = morphology_data, N.steps = 12, N.replicates = 4, grid.size = c(4, 4), parallel = FALSE, overwrite.SOM.results = TRUE, verbose = FALSE, set.seed.N = 123)
            stopifnot(isTRUE(all.equal(SOM_single$quantization_error, test_result$quantization_error, tolerance = 1e-10)))
            TRUE
        })
        clustering_methods <- c("kmeans+BICelbow", "kmeans+BICthreshold", "GMM+BICthreshold", "hierarchical+DB", "HDBSCAN", "OPTICS+Silhouette")
        clustering_results <- list()
        for (clustering_method in clustering_methods) {
            clustering_results[[clustering_method]] <- run_test(paste("clustering.SOM", clustering_method), {
                stopifnot(!is.null(SOM_multi))
                test_result <- delimSOM::clustering.SOM(SOM.output = SOM_multi, max.k = 4, clustering.method = clustering_method, quantization.error.quantile = NULL, topographic.error.quantile = NULL, save.SOM.results = FALSE, overwrite.SOM.results = TRUE, verbose = FALSE, set.seed.N = 123)
                stopifnot(length(test_result$som_clusters) > 0)
                stopifnot(nrow(test_result$cluster_assignment) == N_samples)
                test_result
            })
        }
        SOM_clustered <- run_test("clustering.SOM K=2", {
            stopifnot(!is.null(SOM_multi))
            test_result <- delimSOM::clustering.SOM(SOM.output = SOM_multi, max.k = 4, set.k = 2, clustering.method = "kmeans+BICelbow", quantization.error.quantile = NULL, topographic.error.quantile = NULL, calculate.soft.ancestry = TRUE, calculate.variable.importance = TRUE, overwrite.SOM.results = TRUE, verbose = FALSE, set.seed.N = 123)
            stopifnot(all(test_result$optim_k_vals == 2))
            stopifnot(!is.null(test_result$ancestry_matrix))
            test_result
        })
        run_test("clustering.SOM K=1", {
            stopifnot(!is.null(SOM_multi))
            test_result <- delimSOM::clustering.SOM(SOM.output = SOM_multi, max.k = 4, set.k = 1, clustering.method = "kmeans+BICelbow", quantization.error.quantile = NULL, topographic.error.quantile = NULL, overwrite.SOM.results = TRUE, verbose = FALSE, set.seed.N = 123)
            stopifnot(all(test_result$optim_k_vals == 1))
            TRUE
        })
        run_test("plot.learning.SOM", {
            stopifnot(!is.null(SOM_multi))
            test_plot(delimSOM::plot.learning.SOM(SOM.output = SOM_multi))
        })
        run_test("plot.layer.distance.scale.SOM", {
            stopifnot(!is.null(SOM_multi))
            test_plot(delimSOM::plot.layer.distance.scale.SOM(SOM.output = SOM_multi))
        })
        run_test("plot.K.SOM", {
            stopifnot(!is.null(clustering_results[["kmeans+BICelbow"]]))
            test_plot(delimSOM::plot.K.SOM(SOM.output = clustering_results[["kmeans+BICelbow"]]))
        })
        run_test("plot.model.SOM", {
            stopifnot(!is.null(SOM_clustered))
            test_plot(delimSOM::plot.model.SOM(SOM.output = SOM_clustered, replicate.mode = "first"))
        })
        run_test("plot.structure.SOM", {
            stopifnot(!is.null(SOM_clustered))
            test_plot(delimSOM::plot.structure.SOM(SOM.output = SOM_clustered))
        })
        run_test("plot.map.SOM", {
            stopifnot(!is.null(SOM_clustered))
            test_plot(delimSOM::plot.map.SOM(SOM.output = SOM_clustered, Coordinates = Coordinates, USA.add.states = FALSE, USA.add.counties = FALSE))
        })
        for (importance_mode in c("Cluster.separation", "Map.variance")) {
            run_test(paste("plot.variable.importance.SOM", importance_mode), {
                stopifnot(!is.null(SOM_clustered))
                test_plot(delimSOM::plot.variable.importance.SOM(SOM.output = SOM_clustered, mode = importance_mode))
            })
        }
        run_test("plot.layer.importance.varimp.SOM", {
            stopifnot(!is.null(SOM_clustered))
            test_plot(delimSOM::plot.layer.importance.varimp.SOM(SOM.output = SOM_clustered))
        })
        run_test("plot.layer.importance.leaveoneout.SOM", {
            stopifnot(!is.null(SOM_clustered))
            test_plot(delimSOM::plot.layer.importance.leaveoneout.SOM(SOM_output = SOM_clustered, save.leave.one.layer.out.results = FALSE, overwrite.leave.one.layer.out.results = TRUE))
        })
        for (plot_format in c("svg", "png", "jpg")) {
            run_test(paste("plot.learning.SOM", plot_format), {
                stopifnot(!is.null(SOM_multi))
                plot_path <- tempfile(fileext = paste0(".", plot_format))
                delimSOM::plot.learning.SOM(SOM.output = SOM_multi, save = TRUE, overwrite = TRUE, plot.type = plot_format, file.name = plot_path, resolution = 100)
                stopifnot(file.exists(plot_path))
                stopifnot(file.info(plot_path)$size > 100)
                TRUE
            })
        }
        expected_functions <- c("clustering.SOM", "make.cols.binary.SOM", "plot.K.SOM", "plot.layer.distance.scale.SOM", "plot.layer.importance.leaveoneout.SOM", "plot.layer.importance.varimp.SOM", "plot.learning.SOM", "plot.map.SOM", "plot.model.SOM", "plot.structure.SOM", "plot.variable.importance.SOM", "process.SNP.data.SOM", "remove.lowCV.multicollinearity.SOM", "train.SOM")
        testthat::expect_setequal(getNamespaceExports("delimSOM"), expected_functions)
        testthat::expect_setequal(tested_functions, expected_functions)
        if (length(test_failures) > 0) {
            message("\nFAILED FUNCTION TESTS:")
            message(paste(test_failures, collapse = "\n"))
        }
        testthat::expect_length(test_failures, 0)
    })
}

run_delimSOM_tests()
