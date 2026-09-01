#### Figure 1 ##################################################################

## Set environment
rm(list = ls()) #clear environment
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")


## Install and load required packages
required_packages <- c("emmeans",
                       "glmmTMB",
                       "nlme",
                       "patchwork",
                       "svglite",
                       "tidyverse")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
  library(package_name, character.only = TRUE)
}


## Set main plotting parameters
plot_width_cm <- 16
plot_height_cm <- 15

plot_font_family <- "Arial"
plot_text_color <- "black"

plot_base_size <- 9
plot_axis_title_size <- 9
plot_axis_text_size <- 7.1
plot_column_title_size <- 9

plot_raw_point_alpha <- 0.35
plot_raw_point_size <- 1
plot_estimate_point_size <- 2.6
plot_estimate_linewidth <- 0.7
plot_errorbar_linewidth <- 1.3
plot_errorbar_width <- 0
plot_categorical_jitter_width <- 0.06
plot_SE_multiplier <- 1
X_text_angle_plot <- 90
minimum_column_width_units <- 3.5
panel_margins_mm <- 2.5

plot_color_clustering_methods <- "#56B4E9"
plot_color_neighborhood <- "#009E73"
plot_color_N_steps <- "#E69F00"
plot_color_NA <- "#D55E00"
plot_color_learning_rate_tuning <- "#CC80A7"


## Set input and output
simulation_k3_dir <- file.path("Simulations", "Simulation_set_1", "k3")
sim_results_clustering_methods_csv <- "Sim_results_clustering_methods.csv"
sim_results_neighborhoods_csv <- "Sim_results_neighborhoods.csv"
sim_results_N_steps_csv <- "Sim_results_N_steps.csv"
sim_results_NA_csv <- "Sim_results_NA.csv"
sim_results_learning_rate_tuning_csv <- "Sim_results_learning_rate_tuning.csv"

plot_dir <- "Figure_files"
plot_file_name <- "Figure_1.svg"


## Combine directories and file names
sim_results_clustering_methods_csv <- file.path(simulation_k3_dir, sim_results_clustering_methods_csv)
sim_results_neighborhoods_csv <- file.path(simulation_k3_dir, sim_results_neighborhoods_csv)
sim_results_N_steps_csv <- file.path(simulation_k3_dir, sim_results_N_steps_csv)
sim_results_NA_csv <- file.path(simulation_k3_dir, sim_results_NA_csv)
sim_results_learning_rate_tuning_csv <- file.path(simulation_k3_dir, sim_results_learning_rate_tuning_csv)
combined_svg <- file.path(plot_dir, plot_file_name)


## Create output directory if needed
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


## Check that all completed result files exist
required_result_files <- c(sim_results_clustering_methods_csv,
                           sim_results_neighborhoods_csv,
                           sim_results_N_steps_csv,
                           sim_results_NA_csv,
                           sim_results_learning_rate_tuning_csv)
missing_result_files <- required_result_files[!file.exists(required_result_files)]
if (length(missing_result_files) > 0) stop("The following completed result files are missing: ", paste(missing_result_files, collapse = ", "))


## Define tests
test_specs <- list(
  clustering_methods = list(
    title = "Clustering\nmethod",
    condition_variable = "clustering_method",
    condition_levels = c("kmeans+BICelbow",
                         "kmeans+BICthreshold",
                         "GMM+BICthreshold",
                         "hierarchical+DB",
                         "HDBSCAN",
                         "OPTICS+Silhouette"),
    estimate_color = plot_color_clustering_methods,
    connect_estimates = FALSE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  neighborhoods = list(
    title = "Neighborhood\nfunction",
    condition_variable = "neighborhood_function",
    condition_levels = c("gaussian", "bubble"),
    estimate_color = plot_color_neighborhood,
    connect_estimates = FALSE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  N_steps = list(
    title = "N training\niterations",
    condition_variable = "N_steps",
    condition_levels = c(20,
                         50,
                         100,
                         200,
                         500,
                         1000,
                         2000),
    estimate_color = plot_color_N_steps,
    connect_estimates = TRUE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  learning_rate_tuning = list(
    title = "Learning-rate\ntuning",
    condition_variable = "learning_rate_tuning",
    condition_levels = c(FALSE, TRUE),
    estimate_color = plot_color_learning_rate_tuning,
    connect_estimates = FALSE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  missing_data = list(
    title = "Missing data\nproportion",
    condition_variable = "missing_data_prop",
    condition_levels = c(0.0,
                         0.1,
                         0.2,
                         0.3,
                         0.4,
                         0.5,
                         0.6,
                         0.7,
                         0.8),
    estimate_color = plot_color_NA,
    connect_estimates = TRUE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  )
)


## Calculate category-based column widths
test_category_counts <- vapply(test_specs, function(test_spec) length(test_spec$condition_levels), numeric(1))
column_width_units <- pmax(test_category_counts, minimum_column_width_units)


## Create additional functions
safe.mean <- function(input_vector) {
  if (all(is.na(input_vector))) return(NA_real_)
  mean(input_vector, na.rm = TRUE)
}
convert.to.logical <- function(input_vector) {
  if (is.logical(input_vector)) return(input_vector)
  if (is.numeric(input_vector)) return(ifelse(is.na(input_vector), NA, input_vector != 0))
  input_vector <- tolower(trimws(as.character(input_vector)))
  output_vector <- rep(NA, length(input_vector))
  output_vector[input_vector %in% c("true", "t", "1", "yes", "y")] <- TRUE
  output_vector[input_vector %in% c("false", "f", "0", "no", "n")] <- FALSE
  output_vector
}
prepare.continuous.benchmark.data <- function(input_data,
                                              condition_variable,
                                              response_variable) {
  model_data <- input_data[ , c("sim", condition_variable, response_variable), drop = FALSE]
  colnames(model_data)[colnames(model_data) == condition_variable] <- "condition"
  colnames(model_data)[colnames(model_data) == response_variable] <- "response"
  model_data$response <- suppressWarnings(as.numeric(model_data$response))
  model_data <- model_data[!is.na(model_data$condition) & !is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data <- model_data %>%
    dplyr::group_by(sim, condition) %>%
    dplyr::summarise(response = safe.mean(response), .groups = "drop")
  model_data <- model_data[!is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data
}
prepare.binary.benchmark.data <- function(input_data,
                                          condition_variable,
                                          response_variable) {
  model_data <- input_data[ , c("sim", condition_variable, response_variable), drop = FALSE]
  colnames(model_data)[colnames(model_data) == condition_variable] <- "condition"
  colnames(model_data)[colnames(model_data) == response_variable] <- "response"
  model_data$response <- as.numeric(convert.to.logical(model_data$response))
  model_data <- model_data[!is.na(model_data$condition) & !is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data
}
extract.response.emmeans.columns <- function(emmeans_output_df) {
  response_column <- if ("prob" %in% colnames(emmeans_output_df)) {
    "prob"
  } else if ("response" %in% colnames(emmeans_output_df)) {
    "response"
  } else {
    stop("Could not find the response column in the emmeans output")
  }
  emmeans_output_df$estimated_response <- emmeans_output_df[[response_column]]
  emmeans_output_df$lower_SE <- pmax(0, emmeans_output_df$estimated_response - plot_SE_multiplier * emmeans_output_df$SE)
  emmeans_output_df$upper_SE <- pmin(1, emmeans_output_df$estimated_response + plot_SE_multiplier * emmeans_output_df$SE)
  emmeans_output_df
}


## Create plotting functions
create.benchmark.plot.theme <- function() {
  theme_classic(base_size = plot_base_size, base_family = plot_font_family) +
    theme(text = element_text(family = plot_font_family, colour = plot_text_color),
          axis.title = element_text(family = plot_font_family, size = plot_axis_title_size),
          axis.text = element_text(family = plot_font_family, size = plot_axis_text_size),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          axis.ticks.length = grid::unit(0.8, "mm"),
          axis.text.x.top = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.line.x.top = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          plot.title = element_text(family = plot_font_family,
                                    size = plot_column_title_size,
                                    face = "bold",
                                    hjust = 0.5,
                                    lineheight = 0.9,
                                    margin = margin(b = 1.2, unit = "mm")),
          plot.margin = margin(1, 1, 1, 1, unit = "mm"))
}
plot.categorical.binary.benchmark.model <- function(input_data,
                                                    condition_variable,
                                                    response_variable,
                                                    condition_levels,
                                                    estimate_color,
                                                    connect_estimates = FALSE) {
  model_data <- prepare.binary.benchmark.data(input_data = input_data, condition_variable = condition_variable, response_variable = response_variable)
  model_data$condition_factor <- factor(as.character(model_data$condition), levels = as.character(condition_levels))
  model_data <- model_data[!is.na(model_data$condition_factor), , drop = FALSE]
  model_data$sim_condition <- interaction(model_data$sim, model_data$condition_factor, drop = TRUE)
  model_data$sim_condition <- factor(model_data$sim_condition)
  if (length(unique(model_data$sim)) < 2) stop("Binary mixed model requires at least two simulation replicates")
  model_output <- glmmTMB::glmmTMB(response ~ condition_factor + (1 | sim) + (1 | sim_condition),family = binomial, data = model_data)
  emmeans_output <- emmeans::emmeans(model_output, ~ condition_factor, type = "response")
  emmeans_output_df <- as.data.frame(emmeans_output)
  emmeans_output_df <- extract.response.emmeans.columns(emmeans_output_df)
  raw_plot_data <- model_data %>%
    dplyr::group_by(sim, condition_factor) %>%
    dplyr::summarise(response = safe.mean(response), .groups = "drop")
  model_plot <- ggplot() +
    geom_point(data = raw_plot_data,
               aes(x = condition_factor, y = response),
               alpha = plot_raw_point_alpha,
               size = plot_raw_point_size,
               position = position_jitter(width = plot_categorical_jitter_width, height = 0)) +
    geom_errorbar(data = emmeans_output_df,
                  aes(x = condition_factor, ymin = lower_SE, ymax = upper_SE),
                  colour = estimate_color,
                  linewidth = plot_errorbar_linewidth,
                  width = plot_errorbar_width)
  if (connect_estimates) {
    model_plot <- model_plot +
      geom_line(data = emmeans_output_df,
                aes(x = condition_factor, y = estimated_response, group = 1),
                colour = estimate_color,
                linewidth = plot_estimate_linewidth)
  }
  model_plot <- model_plot +
    geom_point(data = emmeans_output_df,
               aes(x = condition_factor, y = estimated_response),
               colour = estimate_color,
               size = plot_estimate_point_size) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
    create.benchmark.plot.theme()
  list(model = model_output,
       emmeans = emmeans_output,
       emmeans_df = emmeans_output_df,
       plot = model_plot)
}
plot.categorical.gaussian.benchmark.model <- function(input_data,
                                                      condition_variable,
                                                      response_variable,
                                                      condition_levels,
                                                      estimate_color,
                                                      connect_estimates = FALSE) {
  model_data <- prepare.continuous.benchmark.data(input_data = input_data, condition_variable = condition_variable, response_variable = response_variable)
  model_data$condition_factor <- factor(as.character(model_data$condition), levels = as.character(condition_levels))
  model_data <- model_data[!is.na(model_data$condition_factor), , drop = FALSE]
  if (length(unique(model_data$sim)) < 2) stop("Gaussian mixed model requires at least two simulation replicates")
  model_output <- nlme::lme(response ~ condition_factor, random = ~ 1 | sim, data = model_data, na.action = na.omit, method = "REML")
  emmeans_output <- emmeans::emmeans(model_output, ~ condition_factor)
  emmeans_output_df <- as.data.frame(emmeans_output)
  emmeans_output_df$lower_SE <- emmeans_output_df$emmean - plot_SE_multiplier * emmeans_output_df$SE
  emmeans_output_df$upper_SE <- emmeans_output_df$emmean + plot_SE_multiplier * emmeans_output_df$SE
  model_plot <- ggplot() +
    geom_point(data = model_data,
               aes(x = condition_factor, y = response),
               alpha = plot_raw_point_alpha,
               size = plot_raw_point_size,
               position = position_jitter(width = plot_categorical_jitter_width, height = 0)) +
    geom_errorbar(data = emmeans_output_df,
                  aes(x = condition_factor, ymin = lower_SE, ymax = upper_SE),
                  colour = estimate_color,
                  linewidth = plot_errorbar_linewidth,
                  width = plot_errorbar_width)
  if (connect_estimates) {
    model_plot <- model_plot +
      geom_line(data = emmeans_output_df,
                aes(x = condition_factor, y = emmean, group = 1),
                colour = estimate_color,
                linewidth = plot_estimate_linewidth)
  }
  model_plot <- model_plot +
    geom_point(data = emmeans_output_df,
               aes(x = condition_factor, y = emmean),
               colour = estimate_color,
               size = plot_estimate_point_size) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
    create.benchmark.plot.theme()
  list(model = model_output,
       emmeans = emmeans_output,
       emmeans_df = emmeans_output_df,
       plot = model_plot)
}
format.benchmark.panel <- function(input_plot,
                                   test_title = NULL,
                                   metric_label = NULL,
                                   show_x_axis = FALSE,
                                   x_text_angle = 0,
                                   x_text_hjust = 1,
                                   x_text_vjust = 0.5,
                                   panel_tag = NULL) {
  output_plot <- input_plot +
    labs(title = test_title,
         x = NULL,
         y = metric_label,
         tag = panel_tag) +
    theme(axis.title.x = element_blank(),
          axis.text.y = element_text(family = plot_font_family, size = plot_axis_text_size),
          axis.ticks.y = element_line(linewidth = 0.3),
          axis.line.y = element_line(linewidth = 0.3))
  if (is.null(test_title)) output_plot <- output_plot + theme(plot.title = element_blank())
  if (is.null(metric_label)) {
    output_plot <- output_plot + theme(axis.title.y = element_blank())
  } else {
    output_plot <- output_plot + theme(axis.title.y = element_text(family = plot_font_family, size = plot_axis_title_size, margin = margin(r = 2, unit = "mm")))
  }
  if (show_x_axis) {
    output_plot <- output_plot +
      theme(axis.text.x = element_text(family = plot_font_family, size = plot_axis_text_size, angle = x_text_angle, hjust = x_text_hjust, vjust = x_text_vjust),
            axis.ticks.x = element_line(linewidth = 0.3),
            axis.line.x = element_line(linewidth = 0.3))
  } else {
    output_plot <- output_plot +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_line(linewidth = 0.3))
  }
  output_plot
}


## Load and check completed analyses
analysis_data <- list(
  clustering_methods = read.csv(sim_results_clustering_methods_csv, stringsAsFactors = FALSE),
  neighborhoods = read.csv(sim_results_neighborhoods_csv, stringsAsFactors = FALSE),
  N_steps = read.csv(sim_results_N_steps_csv, stringsAsFactors = FALSE),
  missing_data = read.csv(sim_results_NA_csv, stringsAsFactors = FALSE),
  learning_rate_tuning = read.csv(sim_results_learning_rate_tuning_csv, stringsAsFactors = FALSE)
)
for (test_name in names(test_specs)) {
  condition_variable <- test_specs[[test_name]]$condition_variable
  required_columns <- c("sim", condition_variable, "K_correct", "Acc")
  missing_columns <- setdiff(required_columns, colnames(analysis_data[[test_name]]))
  if (length(missing_columns) > 0) stop("The result table for ", test_specs[[test_name]]$title, " is missing the following columns: ", paste(missing_columns, collapse = ", "))
}


## Create top-row K-correct plots
K_correct_plots <- vector(mode = "list", length = length(test_specs))
K_correct_models <- vector(mode = "list", length = length(test_specs))
K_correct_emmeans <- vector(mode = "list", length = length(test_specs))
for (test_index in seq_along(test_specs)) {
  test_name <- names(test_specs)[test_index]
  test_spec <- test_specs[[test_name]]
  model_plot_output <- plot.categorical.binary.benchmark.model(
    input_data = analysis_data[[test_name]],
    condition_variable = test_spec$condition_variable,
    response_variable = "K_correct",
    condition_levels = test_spec$condition_levels,
    estimate_color = test_spec$estimate_color,
    connect_estimates = test_spec$connect_estimates
  )
  metric_label <- if (test_index == 1) "Proportion of correct K" else NULL
  K_correct_plots[[test_index]] <- format.benchmark.panel(input_plot = model_plot_output$plot,
                                                          test_title = test_spec$title,
                                                          metric_label = metric_label,
                                                          show_x_axis = FALSE,
                                                          panel_tag = NULL)
  K_correct_models[[test_index]] <- model_plot_output$model
  K_correct_emmeans[[test_index]] <- model_plot_output$emmeans
  names(K_correct_models)[test_index] <- test_name
  names(K_correct_emmeans)[test_index] <- test_name
}


## Create bottom-row assignment-accuracy plots
Acc_plots <- vector(mode = "list", length = length(test_specs))
Acc_models <- vector(mode = "list", length = length(test_specs))
Acc_emmeans <- vector(mode = "list", length = length(test_specs))
for (test_index in seq_along(test_specs)) {
  test_name <- names(test_specs)[test_index]
  test_spec <- test_specs[[test_name]]
  model_plot_output <- plot.categorical.gaussian.benchmark.model(input_data = analysis_data[[test_name]],
                                                                 condition_variable = test_spec$condition_variable,
                                                                 response_variable = "Acc",
                                                                 condition_levels = test_spec$condition_levels,
                                                                 estimate_color = test_spec$estimate_color,
                                                                 connect_estimates = test_spec$connect_estimates)
  metric_label <- if (test_index == 1) "Assignment accuracy" else NULL
  Acc_plots[[test_index]] <- format.benchmark.panel(input_plot = model_plot_output$plot,
                                                    metric_label = metric_label,
                                                    show_x_axis = TRUE,
                                                    x_text_angle = test_spec$x_text_angle,
                                                    x_text_hjust = test_spec$x_text_hjust,
                                                    x_text_vjust = test_spec$x_text_vjust,
                                                    panel_tag = NULL)
  Acc_models[[test_index]] <- model_plot_output$model
  Acc_emmeans[[test_index]] <- model_plot_output$emmeans
  names(Acc_models)[test_index] <- test_name
  names(Acc_emmeans)[test_index] <- test_name
}


## Combine plots using category-based column widths
combined_plot <- patchwork::wrap_plots(plotlist = c(K_correct_plots, Acc_plots),
                                       ncol = length(test_specs),
                                       nrow = 2,
                                       byrow = TRUE) +
  patchwork::plot_layout(widths = unname(column_width_units), heights = c(1, 1)) &
  theme(text = element_text(family = plot_font_family),
        plot.margin = margin(panel_margins_mm, panel_margins_mm, panel_margins_mm, panel_margins_mm, unit = "mm"))


## Show and save combined plot
combined_plot
ggsave(filename = combined_svg,
       plot = combined_plot,
       device = svglite::svglite,
       width = plot_width_cm,
       height = plot_height_cm,
       units = "cm",
       bg = "white",
       limitsize = FALSE)
message("Combined SVG saved to: ", combined_svg)




#### Figure 2 ##################################################################

## Set environment
rm(list = ls()) #clear environment
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
invisible(gc())


## Install and load required packages
required_packages <- c("svglite", "tidyverse", "viridisLite")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
  library(package_name, character.only = TRUE)
}


## Set main plotting parameters
plot_width_cm <- 15.36
plot_height_cm <- 19.24

plot_font_family <- "Arial"
plot_text_color <- "black"

plot_base_size <- 9
plot_axis_title_size <- 9
plot_axis_text_size <- 7.1
plot_column_title_size <- 9
plot_row_title_size <- 9
plot_legend_title_size <- 9
plot_legend_text_size <- 9

plot_axis_title_margin_mm <- 3
legend_box_margin_mm <- 1
panel_spacing_mm <- 2.5
plot_column_title_margin_mm <- 2.5
plot_x_axis_left_expansion <- 10000
plot_x_axis_right_expansion <- 10000

plot_point_size <- 1
plot_point_alpha <- 0.7
plot_fitted_line_width <- 1
plot_threshold_line_width <- 0.6
plot_threshold_line_type <- "dashed"
plot_K2_threshold_proportion <- 0.5
plot_reference_line_positions <- c(50000, 100000, 150000, 200000)
plot_reference_line_width <- 0.1
legend_box_line_width <- 0.3
migration_values <- c(0, 1e-6, 4e-6, 7e-6)
migration_labels <- c("0", "1e-6", "4e-6", "7e-6")
migration_colors <- stats::setNames(c("#2B0A3D", "#31668EFF", "#30B57BFF", "#F0DD1CFF"), migration_labels)


## Set input and output
simulation_dir <- file.path("Simulations", "Simulation_set_2")
results_root_dir <- file.path(simulation_dir, "fastsimcoal2_results")
balanced_results_dir <- file.path(results_root_dir, "symmetric_24")
unbalanced_results_dir <- file.path(results_root_dir, "asymmetric_24")
balanced_SOM_results_csv <- file.path(balanced_results_dir, "fastsimcoal2_SOM_results.csv")
balanced_SOM_optim_k_csv <- file.path(balanced_results_dir, "fastsimcoal2_SOM_optim_k_summary.csv")
balanced_STRUCTURE_csv <- file.path(balanced_results_dir, "STRUCTURE_results", "STRUCTURE_binomial_table_by_vcf.csv")
unbalanced_SOM_results_csv <- file.path(unbalanced_results_dir, "fastsimcoal2_SOM_results.csv")
unbalanced_SOM_optim_k_csv <- file.path(unbalanced_results_dir, "fastsimcoal2_SOM_optim_k_summary.csv")
unbalanced_STRUCTURE_csv <- file.path(unbalanced_results_dir, "STRUCTURE_results", "STRUCTURE_binomial_table_by_vcf.csv")
plot_dir <- "Figure_files"
plot_file_name <- "Figure_2.svg"
combined_svg <- file.path(plot_dir, plot_file_name)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


## Check that all completed result files exist
required_result_files <- c(balanced_SOM_results_csv,
                           balanced_SOM_optim_k_csv,
                           balanced_STRUCTURE_csv,
                           unbalanced_SOM_results_csv,
                           unbalanced_SOM_optim_k_csv,
                           unbalanced_STRUCTURE_csv)
missing_result_files <- required_result_files[!file.exists(required_result_files)]
if (length(missing_result_files) > 0) stop("The following completed result files are missing: ", paste(missing_result_files, collapse = ", "))


## Create additional functions
check.required.columns <- function(input_table,
                                   required_columns,
                                   table_name) {
  missing_columns <- setdiff(required_columns, colnames(input_table))
  if (length(missing_columns) > 0) stop(table_name, " is missing the following columns: ", paste(missing_columns, collapse = ", "))
}
standardize.migration.labels <- function(input_data) {
  migration_indices <- vapply(input_data$mig, function(current_migration) {
    migration_difference <- abs(migration_values - current_migration)
    matching_index <- which.min(migration_difference)
    if (length(matching_index) != 1 || !is.finite(migration_difference[matching_index]) || migration_difference[matching_index] > sqrt(.Machine$double.eps)) return(NA_integer_)
    matching_index
  }, integer(1))
  if (any(is.na(migration_indices))) stop("At least one migration rate could not be matched to the expected migration rates")
  input_data$mig.tag <- factor(migration_labels[migration_indices], levels = migration_labels)
  input_data
}
create.analysis.binomial.table <- function(SOM_results_csv,
                                           SOM_optim_k_csv,
                                           STRUCTURE_csv,
                                           design_label) {
  result_table <- read.csv(SOM_results_csv, stringsAsFactors = FALSE)
  optim_k_result_table <- read.csv(SOM_optim_k_csv, stringsAsFactors = FALSE)
  STRUCTURE_binomial_table <- read.csv(STRUCTURE_csv, stringsAsFactors = FALSE)
  check.required.columns(result_table, c("file", "status", "mig", "tdiv", "deNovo.kmeans.best.k", "sNMF.best.k"), paste0(design_label, " SOM result table"))
  check.required.columns(optim_k_result_table, c("file", "mig", "tdiv", "Count", "k.label"), paste0(design_label, " SOM optim-k table"))
  check.required.columns(STRUCTURE_binomial_table, c("file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2"), paste0(design_label, " STRUCTURE binomial table"))
  plot_result_table <- result_table[result_table$status == "ok", ]
  plot_result_table <- plot_result_table[order(plot_result_table$mig, plot_result_table$tdiv), ]
  if (nrow(plot_result_table) == 0) stop("No successful SOM result rows are available for ", design_label)
  SOM_optim_k_count_table <- optim_k_result_table[!is.na(optim_k_result_table$Count) & optim_k_result_table$file %in% plot_result_table$file, ]
  SOM_optim_k_count_table$Count <- as.numeric(as.character(SOM_optim_k_count_table$Count))
  SOM_total_count_table <- stats::aggregate(Count ~ file + mig + tdiv, data = SOM_optim_k_count_table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE))
  colnames(SOM_total_count_table)[colnames(SOM_total_count_table) == "Count"] <- "n.total"
  SOM_k2_count_table <- SOM_optim_k_count_table[SOM_optim_k_count_table$k.label == "k2", ]
  if (nrow(SOM_k2_count_table) > 0) {
    SOM_k2_count_table <- stats::aggregate(Count ~ file + mig + tdiv, data = SOM_k2_count_table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE))
    colnames(SOM_k2_count_table)[colnames(SOM_k2_count_table) == "Count"] <- "n.k2"
  } else {
    SOM_k2_count_table <- SOM_total_count_table[, c("file", "mig", "tdiv")]
    SOM_k2_count_table$n.k2 <- 0
  }
  SOM_binomial_table <- merge(SOM_total_count_table, SOM_k2_count_table, by = c("file", "mig", "tdiv"), all.x = TRUE)
  SOM_binomial_table$n.k2[is.na(SOM_binomial_table$n.k2)] <- 0
  SOM_binomial_table$n.not.k2 <- SOM_binomial_table$n.total - SOM_binomial_table$n.k2
  SOM_binomial_table$proportion.k2 <- SOM_binomial_table$n.k2 / SOM_binomial_table$n.total
  SOM_binomial_table$method <- "SOM"
  SOM_binomial_table <- SOM_binomial_table[, c("method", "file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")]
  DAPC_binomial_source_table <- result_table[result_table$status == "ok" & !is.na(result_table$deNovo.kmeans.best.k), ]
  if ("deNovo.kmeans.status" %in% colnames(DAPC_binomial_source_table)) DAPC_binomial_source_table <- DAPC_binomial_source_table[DAPC_binomial_source_table$deNovo.kmeans.status == "ok", ]
  DAPC_best_k <- as.integer(as.character(DAPC_binomial_source_table$deNovo.kmeans.best.k))
  DAPC_binomial_table <- data.frame(method = "DAPC",
                                    file = DAPC_binomial_source_table$file,
                                    mig = DAPC_binomial_source_table$mig,
                                    tdiv = DAPC_binomial_source_table$tdiv,
                                    n.k2 = as.integer(DAPC_best_k == 2L),
                                    n.not.k2 = as.integer(DAPC_best_k != 2L),
                                    proportion.k2 = as.numeric(DAPC_best_k == 2L),
                                    stringsAsFactors = FALSE)
  sNMF_binomial_source_table <- result_table[result_table$status == "ok" & !is.na(result_table$sNMF.best.k), ]
  if ("sNMF.status" %in% colnames(sNMF_binomial_source_table)) sNMF_binomial_source_table <- sNMF_binomial_source_table[sNMF_binomial_source_table$sNMF.status == "ok", ]
  sNMF_best_k <- as.integer(as.character(sNMF_binomial_source_table$sNMF.best.k))
  sNMF_binomial_table <- data.frame(method = "sNMF",
                                    file = sNMF_binomial_source_table$file,
                                    mig = sNMF_binomial_source_table$mig,
                                    tdiv = sNMF_binomial_source_table$tdiv,
                                    n.k2 = as.integer(sNMF_best_k == 2L),
                                    n.not.k2 = as.integer(sNMF_best_k != 2L),
                                    proportion.k2 = as.numeric(sNMF_best_k == 2L),
                                    stringsAsFactors = FALSE)
  STRUCTURE_binomial_table <- STRUCTURE_binomial_table[, c("file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")]
  STRUCTURE_binomial_table$method <- "STRUCTURE"
  STRUCTURE_binomial_table <- STRUCTURE_binomial_table[, c("method", "file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")]
  combined_binomial_table <- rbind(SOM_binomial_table, DAPC_binomial_table, sNMF_binomial_table, STRUCTURE_binomial_table)
  combined_binomial_table$design <- design_label
  combined_binomial_table <- standardize.migration.labels(combined_binomial_table)
  combined_binomial_table
}
fit.K2.binomial.models <- function(input_data) {
  prediction_list <- list()
  threshold_list <- list()
  prediction_index <- 1
  threshold_index <- 1
  for (current_design in levels(input_data$design)) {
    for (current_method in levels(input_data$method)) {
      for (current_migration in levels(input_data$mig.tag)) {
        current_binomial_table <- input_data[input_data$design == current_design & input_data$method == current_method & input_data$mig.tag == current_migration, ]
        current_binomial_table <- current_binomial_table[is.finite(current_binomial_table$tdiv), ]
        if (nrow(current_binomial_table) == 0) next
        current_binomial_table <- current_binomial_table[order(current_binomial_table$tdiv), ]
        prediction_tdiv <- seq(min(current_binomial_table$tdiv), max(current_binomial_table$tdiv), length.out = 1000)
        if (all(current_binomial_table$n.k2 == 0)) {
          current_prediction_table <- data.frame(design = current_design,
                                                 method = current_method,
                                                 mig.tag = current_migration,
                                                 tdiv = prediction_tdiv,
                                                 fitted.proportion.k2 = 0,
                                                 stringsAsFactors = FALSE)
        } else if (all(current_binomial_table$n.not.k2 == 0)) {
          current_prediction_table <- data.frame(design = current_design,
                                                 method = current_method,
                                                 mig.tag = current_migration,
                                                 tdiv = prediction_tdiv,
                                                 fitted.proportion.k2 = 1,
                                                 stringsAsFactors = FALSE)
        } else {
          current_fit <- stats::glm(cbind(n.k2, n.not.k2) ~ tdiv, data = current_binomial_table, family = stats::binomial(link = "logit"), control = stats::glm.control(maxit = 100))
          current_prediction_table <- data.frame(design = current_design,
                                                 method = current_method,
                                                 mig.tag = current_migration,
                                                 tdiv = prediction_tdiv,
                                                 fitted.proportion.k2 = stats::predict(current_fit, newdata = data.frame(tdiv = prediction_tdiv), type = "response"),
                                                 stringsAsFactors = FALSE)
          current_coefficients <- stats::coef(current_fit)
          current_intercept <- unname(current_coefficients["(Intercept)"])
          current_tdiv_slope <- unname(current_coefficients["tdiv"])
          if (is.finite(current_intercept) && is.finite(current_tdiv_slope) && current_tdiv_slope != 0) {
            current_threshold_tdiv <- (stats::qlogis(plot_K2_threshold_proportion) - current_intercept) / current_tdiv_slope
            threshold_within_range <- current_threshold_tdiv >= min(current_binomial_table$tdiv) && current_threshold_tdiv <= max(current_binomial_table$tdiv)
            if (threshold_within_range) {
              threshold_list[[threshold_index]] <- data.frame(design = current_design,
                                                              method = current_method,
                                                              mig.tag = current_migration,
                                                              tdiv.at.threshold.k2 = current_threshold_tdiv,
                                                              stringsAsFactors = FALSE)
              threshold_index <- threshold_index + 1
            }
          }
        }
        prediction_list[[prediction_index]] <- current_prediction_table
        prediction_index <- prediction_index + 1
      }
    }
  }
  prediction_table <- do.call(rbind, prediction_list)
  if (length(threshold_list) > 0) {
    threshold_table <- do.call(rbind, threshold_list)
  } else {
    threshold_table <- data.frame(design = character(0), method = character(0), mig.tag = character(0), tdiv.at.threshold.k2 = numeric(0))
  }
  prediction_table$design <- factor(prediction_table$design, levels = levels(input_data$design))
  prediction_table$method <- factor(prediction_table$method, levels = levels(input_data$method))
  prediction_table$mig.tag <- factor(prediction_table$mig.tag, levels = migration_labels)
  threshold_table$design <- factor(threshold_table$design, levels = levels(input_data$design))
  threshold_table$method <- factor(threshold_table$method, levels = levels(input_data$method))
  threshold_table$mig.tag <- factor(threshold_table$mig.tag, levels = migration_labels)
  list(prediction_table = prediction_table, threshold_table = threshold_table)
}


## Load balanced analysis
balanced_binomial_table <- create.analysis.binomial.table(SOM_results_csv = balanced_SOM_results_csv,
                                                          SOM_optim_k_csv = balanced_SOM_optim_k_csv,
                                                          STRUCTURE_csv = balanced_STRUCTURE_csv,
                                                          design_label = "Balanced")


## Load unbalanced analysis
unbalanced_binomial_table <- create.analysis.binomial.table(SOM_results_csv = unbalanced_SOM_results_csv,
                                                            SOM_optim_k_csv = unbalanced_SOM_optim_k_csv,
                                                            STRUCTURE_csv = unbalanced_STRUCTURE_csv,
                                                            design_label = "Unbalanced")


## Combine balanced and unbalanced analyses
combined_binomial_table <- rbind(balanced_binomial_table, unbalanced_binomial_table)
combined_binomial_table$design <- factor(combined_binomial_table$design, levels = c("Balanced", "Unbalanced"))
combined_binomial_table$method <- factor(combined_binomial_table$method, levels = c("SOM", "DAPC", "sNMF", "STRUCTURE"))


## Fit binomial models
fitted_K2_output <- fit.K2.binomial.models(combined_binomial_table)
fitted_prediction_table <- fitted_K2_output$prediction_table
threshold_table <- fitted_K2_output$threshold_table


## Create Figure 2
combined_plot <- ggplot() +
  geom_vline(xintercept = plot_reference_line_positions,
             linewidth = plot_reference_line_width,
             colour = plot_text_color) +
  geom_point(data = combined_binomial_table,
             aes(x = tdiv, y = proportion.k2, colour = mig.tag, group = mig.tag),
             size = plot_point_size,
             alpha = plot_point_alpha) +
  geom_line(data = fitted_prediction_table,
            aes(x = tdiv, y = fitted.proportion.k2, colour = mig.tag, group = mig.tag),
            linewidth = plot_fitted_line_width) +
  geom_vline(data = threshold_table,
             aes(xintercept = tdiv.at.threshold.k2, colour = mig.tag),
             linewidth = plot_threshold_line_width,
             linetype = plot_threshold_line_type,
             show.legend = FALSE) +
  facet_grid(rows = vars(method),
             cols = vars(design),
             labeller = labeller(design = c("Balanced" = "Even sampling", "Unbalanced" = "Uneven sampling")),
             axes = "all",
             axis.labels = "margins") +
  scale_color_manual(values = migration_colors,
                     breaks = migration_labels,
                     labels = migration_labels,
                     drop = FALSE) +
  scale_x_continuous(expand = expansion(add = c(plot_x_axis_left_expansion, plot_x_axis_right_expansion))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.25),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "Divergence time (generations)",
       y = "Proportion choosing k = 2",
       colour = "Migration rate") +
  theme_classic(base_size = plot_base_size,
                base_family = plot_font_family) +
  theme(text = element_text(family = plot_font_family,
                            colour = plot_text_color),
        axis.title.x = element_text(family = plot_font_family,
                                    size = plot_axis_title_size,
                                    colour = plot_text_color,
                                    margin = margin(t = plot_axis_title_margin_mm, unit = "mm")),
        axis.title.y = element_text(family = plot_font_family,
                                    size = plot_axis_title_size,
                                    colour = plot_text_color,
                                    margin = margin(r = plot_axis_title_margin_mm, unit = "mm")),
        axis.text = element_text(family = plot_font_family,
                                 size = plot_axis_text_size,
                                 colour = plot_text_color),
        strip.background = element_blank(),
        strip.text.x = element_text(family = plot_font_family,
                                    size = plot_column_title_size,
                                    face = "bold",
                                    colour = plot_text_color,
                                    margin = margin(b = plot_column_title_margin_mm, unit = "mm")),
        strip.text.y = element_text(family = plot_font_family,
                                    size = plot_row_title_size,
                                    face = "bold",
                                    colour = plot_text_color),
        legend.position = "bottom",
        legend.title = element_text(family = plot_font_family,
                                    size = plot_legend_title_size,
                                    face = "bold",
                                    colour = plot_text_color),
        legend.text = element_text(family = plot_font_family,
                                   size = plot_legend_text_size,
                                   colour = plot_text_color),
        legend.box.background = element_rect(colour = plot_text_color,
                                             fill = NA,
                                             linewidth = legend_box_line_width),
        legend.box.margin = margin(legend_box_margin_mm,
                                   legend_box_margin_mm,
                                   legend_box_margin_mm,
                                   legend_box_margin_mm,
                                   unit = "mm"),
        panel.spacing = grid::unit(panel_spacing_mm, "mm"),
        plot.margin = margin(1, 1, 1, 1, unit = "mm")) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))


## Show and save figure
combined_plot
ggsave(filename = combined_svg,
       plot = combined_plot,
       device = svglite::svglite,
       width = plot_width_cm,
       height = plot_height_cm,
       units = "cm",
       bg = "white",
       limitsize = FALSE)
message("Combined SVG saved to: ", combined_svg)




#### Figure 3 ##################################################################

## Install packages
rm(list = ls()) #clear environment
CRAN_packages <- c("tibble", "dplyr", "vcfR") #CRAN packages
Bioconductor_packages <- "SeqArray" #Bioconductor package
for(p in CRAN_packages) if(!requireNamespace(p, quietly = TRUE)) install.packages(p) #install missing CRAN packages
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") #install BiocManager
if(!requireNamespace(Bioconductor_packages, quietly = TRUE)) BiocManager::install(Bioconductor_packages, ask = FALSE, update = FALSE) #install missing Bioconductor package
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv") #import (NicheDiv::transform.skewed.variables function from my other package (under development)
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")


## Set directories
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
intermediate_files_folder <- "Empirical_examples/Intermediate_files"
figure_files_folder <- "Figure_files"


## Import and prepare data
Polygonia_SNP <- process.SNP.data.SOM(vcf.path = "./Empirical_examples/Dupuis_et_al_2018/Polygonia_961SNPs.vcf", #filter loci and individuals and create SNP matrix dataframe
                                      missing.loci.cutoff.lenient = 0.7,
                                      missing.loci.cutoff.final = 0.5,
                                      missing.individuals.cutoff = 0.5)
rownames(Polygonia_SNP) <- sub(".*?(\\d+)$", "\\1", rownames(Polygonia_SNP)) #only keep numeric identifier as rownames
Polygonia_COI <- process.SNP.data.SOM(nexus.path = "./Empirical_examples/Dupuis_et_al_2018/Polygonia_COI.nex",
                                      missing.loci.cutoff.lenient = 0.7, 
                                      missing.loci.cutoff.final = 0.5,
                                      missing.individuals.cutoff = 0.5)
Polygonia_COI_numeric_rownames <- sub(".*?(\\d+)$", "\\1", rownames(Polygonia_COI)) #extract numeric code from each rowname (e.g., "pf_8301" -> "8301")
Polygonia_COI <- Polygonia_COI[!duplicated(Polygonia_COI_numeric_rownames), , drop = FALSE] #keep only first occurrence for each numeric code (remove duplicates)
rownames(Polygonia_COI) <- Polygonia_COI_numeric_rownames[!duplicated(Polygonia_COI_numeric_rownames)] #set rownames to unique numeric codes
Polygonia_RGB <- read.delim("./Empirical_examples/Dupuis_et_al_2018/Polygonia_RGB_characters.txt", stringsAsFactors = FALSE)
rownames(Polygonia_RGB) <- Polygonia_RGB$Species
Polygonia_RGB <- magrittr::`%>%`(Polygonia_RGB, dplyr::select(-Name, -Species)) #remove columns
Polygonia_wing_scores <- read.delim("./Empirical_examples/Dupuis_et_al_2018/Polygonia_visually_scored.txt", stringsAsFactors = FALSE)
rownames(Polygonia_wing_scores) <- Polygonia_wing_scores$Name
Polygonia_wing_scores <- Polygonia_wing_scores |>
  dplyr::select(-Name, -Species) |> #remove columns
  dplyr::rename(Wing_character_1 = Ch1,
                Wing_character_2 = Ch2,
                Wing_character_3 = Ch3,
                Wing_character_4 = Ch4,
                Wing_character_5 = Ch5,
                Wing_character_6 = Ch6,
                Wing_character_7 = Ch7,
                Wing_character_8 = Ch8,
                Wing_character_9 = Ch9,
                Wing_character_10 = Ch10)
Polygonia_wing_scores$Wing_character_8 <- factor(Polygonia_wing_scores$Wing_character_8, levels = 1:4)
Wing_character_8_states <- stats::model.matrix(~ Wing_character_8 - 1, data = Polygonia_wing_scores) #make wing character 8 binary (since it is not ordinal)
colnames(Wing_character_8_states) <- paste0("Wing_character_8_state_", 1:4)
Polygonia_wing_scores <- cbind(Polygonia_wing_scores, Wing_character_8_states)
Polygonia_wing_scores$Wing_character_8 <- NULL
Polygonia_metadata <- read.csv("./Empirical_examples/Dupuis_et_al_2018/Polygonia_metadata.csv", header = TRUE, sep = ";")
rownames(Polygonia_metadata) <- Polygonia_metadata$ID

Polygonia_spatial <- dplyr::select(Polygonia_metadata, Latitude, Longitude) #create dataframe with Lat and Long
Polygonia_spatial$Elevation <- NA #initialize elevation column with NA
Polygonia_spatial_sf <- sf::st_as_sf(Polygonia_spatial[!is.na(Polygonia_spatial$Latitude) & !is.na(Polygonia_spatial$Longitude), ], coords = c("Longitude", "Latitude"), crs = 4326) #extract elevation
Polygonia_spatial$Elevation[!is.na(Polygonia_spatial$Latitude) & !is.na(Polygonia_spatial$Longitude)] <-
  elevatr::get_elev_point(locations = Polygonia_spatial_sf, prj = sf::st_crs(Polygonia_spatial_sf)$proj4string, src = "aws")$elevation
Polygonia_morphotype <- dplyr::select(Polygonia_metadata, Morphotype)
Polygonia_metadata <- dplyr::select(Polygonia_metadata, Species, ID)
rownames(Polygonia_RGB) <- as.character(rownames(Polygonia_RGB))
rownames(Polygonia_wing_scores) <- as.character(rownames(Polygonia_wing_scores))
rownames(Polygonia_morphotype) <- as.character(rownames(Polygonia_morphotype))
Polygonia_RGB_wing_scores <- merge(Polygonia_RGB,
                                   Polygonia_wing_scores,
                                   by = "row.names",
                                   all = FALSE)
rownames(Polygonia_RGB_wing_scores) <- Polygonia_RGB_wing_scores$Row.names
Polygonia_RGB_wing_scores$Row.names <- NULL
Polygonia_morphology <- merge(Polygonia_RGB_wing_scores,
                              Polygonia_morphotype,
                              by = "row.names",
                              all = FALSE)
rownames(Polygonia_morphology) <- Polygonia_morphology$Row.names
Polygonia_morphology$Row.names <- NULL
Polygonia_morphology <- make.cols.binary.SOM(Polygonia_morphology, #convert Morphotype to binary columns and remove original
                                             make.binary.cols = "Morphotype",
                                             append.to.original = TRUE)
Polygonia_morphology$Morphotype <- NULL
non.continuous.cols <- grepl("^Wing_character_|^Morphotype_", colnames(Polygonia_morphology)) #identify non-continuous traits
Polygonia_morphology_categorical <- Polygonia_morphology[, non.continuous.cols, drop = FALSE] #non-continuous traits
Polygonia_morphology <- Polygonia_morphology[, !non.continuous.cols, drop = FALSE] #continuous traits
Polygonia_morphology <- remove.lowCV.multicollinearity.SOM(Polygonia_morphology, #remove highly correlated and low-variance variables
                                                           CV.threshold = 0.05,
                                                           cor.threshold = 0.9)
Polygonia_environmental <- read.csv("./Empirical_examples/Dupuis_et_al_2018/Polygonia_environmental.csv",
                                    row.names = 1, header = TRUE)
Polygonia_environmental <- dplyr::select(Polygonia_environmental, -Latitude, -Longitude, -Elevation)
Polygonia_environmental_rownames <- rownames(Polygonia_environmental) #save rownames
Polygonia_environmental <- as.data.frame(lapply(Polygonia_environmental, as.numeric)) #ensure all columns are numeric
rownames(Polygonia_environmental) <- Polygonia_environmental_rownames #reassign saved row names
Polygonia_environmental <- (NicheDiv::transform.skewed.variables(Polygonia_environmental))$transformed #transform skewed variables
Polygonia_environmental <- remove.lowCV.multicollinearity.SOM(Polygonia_environmental, #remove highly correlated and low-variance variables
                                                              CV.threshold = 0.05,
                                                              cor.threshold = 0.9)
for (Polygonia_shared_data in c("Polygonia_morphology", "Polygonia_morphology_categorical", "Polygonia_SNP", "Polygonia_COI", "Polygonia_spatial", "Polygonia_environmental", "Polygonia_metadata")) {
  Polygonia_shared_data_mat <- get(Polygonia_shared_data)
  rownames(Polygonia_shared_data_mat) <- make.unique(as.character(rownames(Polygonia_shared_data_mat)))
  assign(Polygonia_shared_data, Polygonia_shared_data_mat, envir = .GlobalEnv)}
Polygonia_common_IDs <- Reduce(intersect, list(
  rownames(Polygonia_morphology),
  rownames(Polygonia_morphology_categorical),
  rownames(Polygonia_SNP),
  rownames(Polygonia_COI),
  rownames(Polygonia_spatial),
  rownames(Polygonia_environmental),
  rownames(Polygonia_metadata)))
Polygonia_morphology2 <- Polygonia_morphology[Polygonia_common_IDs, , drop = FALSE]
Polygonia_morphology_categorical2 <- Polygonia_morphology_categorical[Polygonia_common_IDs, , drop = FALSE]
Polygonia_SNP2 <- Polygonia_SNP[Polygonia_common_IDs, , drop = FALSE]
Polygonia_COI2 <- Polygonia_COI[Polygonia_common_IDs, , drop = FALSE]
Polygonia_spatial2 <- Polygonia_spatial[Polygonia_common_IDs, , drop = FALSE]
Polygonia_environmental2 <- Polygonia_environmental[Polygonia_common_IDs, , drop = FALSE]
Polygonia_metadata2 <- Polygonia_metadata[Polygonia_common_IDs, , drop = FALSE]
Polygonia_morphology <- Polygonia_morphology2[rowSums(!is.na(Polygonia_morphology2)) > 0, , drop = FALSE]
Polygonia_morphology_categorical <- Polygonia_morphology_categorical2[rowSums(!is.na(Polygonia_morphology_categorical2)) > 0, , drop = FALSE]
Polygonia_SNP <- Polygonia_SNP2[rowSums(!is.na(Polygonia_SNP2)) > 0, , drop = FALSE]
Polygonia_COI <- Polygonia_COI2[rowSums(!is.na(Polygonia_COI2)) > 0, , drop = FALSE]
Polygonia_spatial <- Polygonia_spatial2[rowSums(!is.na(Polygonia_spatial2)) > 0, , drop = FALSE]
Polygonia_environmental <- Polygonia_environmental2[rowSums(!is.na(Polygonia_environmental2)) > 0, , drop = FALSE]
Polygonia_metadata <- Polygonia_metadata2[rowSums(!is.na(Polygonia_metadata2)) > 0, , drop = FALSE]
Polygonia_species_vec <- Polygonia_metadata[rownames(Polygonia_morphology), "Species"]
Polygonia_new_rownames <- paste(rownames(Polygonia_morphology), Polygonia_species_vec, sep = "_")
rownames(Polygonia_morphology) <- Polygonia_new_rownames
rownames(Polygonia_morphology_categorical) <- Polygonia_new_rownames
rownames(Polygonia_SNP) <- Polygonia_new_rownames
rownames(Polygonia_COI) <- Polygonia_new_rownames
rownames(Polygonia_spatial) <- Polygonia_new_rownames
rownames(Polygonia_metadata) <- Polygonia_new_rownames
rownames(Polygonia_environmental) <- Polygonia_new_rownames
Polygonia_all_data <- list(Morphology = Polygonia_morphology,
                           Morphology_2 = Polygonia_morphology_categorical,
                           SNP = Polygonia_SNP,
                           COI = Polygonia_COI,
                           Environmental = Polygonia_environmental,
                           Spatial = Polygonia_spatial)
Polygonia_SOM_tr <- train.SOM(input_data = Polygonia_all_data, #200 samples, 4.1min
                              save.SOM.results = TRUE,
                              save.SOM.results.name = file.path(intermediate_files_folder, "Polygonia_SOM_tr.Rdata"),
                              max.NA.row = 0.5,
                              max.NA.col = 0.5)
Polygonia_SOM <- clustering.SOM(Polygonia_SOM_tr, #3.0min
                                clustering.method = "kmeans+BICelbow",
                                save.SOM.results.name = file.path(intermediate_files_folder, "Polygonia_SOM_kmeansBICelbow.Rdata"))


## Figure 3A
plot_width_cm <- 8.83
plot_height_cm <- 5.93
bottom_tick_label_gap <- 0.5
side_tick_label_gap <- 0.65
y_axis_breaks <- c(0.00, 0.01, 0.02, 0.03)
y_axis_labels <- sprintf("%.2f", y_axis_breaks)
bottom_margin <- 2
left_margin <- 2
top_margin <- 1
right_margin <- 0.5
bottom_y_axis_padding_fraction <- 0.14
top_y_axis_padding_fraction <- 0
lines_alpha <- 0.25
lines_thickness <- 0.95
figure_name <- "Figure_3A.svg"
if ("learning_values" %in% names(Polygonia_SOM)) {
  learning_values_list <- list(Polygonia_SOM$learning_values)
} else {
  learning_values_list <- Polygonia_SOM$learning_values_list
}
learning_values_list <- lapply(learning_values_list, function(learning_values_matrix) {
  if (is.data.frame(learning_values_matrix)) return(as.matrix(learning_values_matrix))
  learning_values_matrix
})
layer_names <- Polygonia_SOM$input_data_names
finite_learning_values <- unlist(lapply(learning_values_list, function(learning_values_matrix) learning_values_matrix[is.finite(learning_values_matrix)]))
global_y_limits <- range(finite_learning_values, na.rm = TRUE)
if (diff(global_y_limits) == 0) {
  global_y_limits <- global_y_limits + c(-0.5, 0.5)
} else {
  y_axis_range <- diff(global_y_limits)
  bottom_y_axis_padding <- bottom_y_axis_padding_fraction * y_axis_range
  top_y_axis_padding <- top_y_axis_padding_fraction * y_axis_range
  global_y_limits <- c(global_y_limits[1] - bottom_y_axis_padding, global_y_limits[2] + top_y_axis_padding)
}
svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72), family = "Arial")
base_font_size <- par("ps")
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
first_learning_values_matrix <- learning_values_list[[1]]
number_of_training_steps <- nrow(first_learning_values_matrix)
par(mfrow = c(1, 1), mar = c(bottom_margin, left_margin, top_margin, right_margin), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1, family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black")
plot(NULL, xlim = c(1, number_of_training_steps), ylim = global_y_limits, xlab = "", ylab = "", main = "", axes = FALSE)
axis(1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size, col = "black", col.axis = "black", font.axis = 1, family = "Arial")
axis(2, at = y_axis_breaks, labels = y_axis_labels, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size, col = "black", col.axis = "black", font.axis = 1, family = "Arial")
layer_colors <- setNames(viridis::turbo(length(layer_names)), layer_names)
for (layer_index in seq_along(learning_values_list)) {
  current_learning_values_matrix <- learning_values_list[[layer_index]]
  current_layer_name <- layer_names[layer_index]
  for (replicate_index in seq_len(ncol(current_learning_values_matrix))) {
    lines(seq_len(nrow(current_learning_values_matrix)), current_learning_values_matrix[, replicate_index], col = adjustcolor(layer_colors[current_layer_name], alpha.f = lines_alpha), lwd = lines_thickness)
  }
}
dev.off()

plot_width_cm <- 6
plot_height_cm <- 1.6
legend_line_thickness <- 2.5
legend_text_font_size <- 7.1
figure_name <- "Figure_3A_legend.svg"

layer_names <- Polygonia_SOM$input_data_names
legend_labels <- sub("^Morphology_", "Morphology ", layer_names)
layer_colors <- setNames(viridis::turbo(length(layer_names)), layer_names)
svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72), family = "Arial", bg = "transparent")
base_font_size <- par("ps")
legend_text_relative_font_size <- (legend_text_font_size * (96 / 72)) / base_font_size
par(mar = c(0, 0, 0, 0), bg = NA, family = "Arial", fg = "black", col = "black")
plot.new()
legend("center",
       legend = legend_labels,
       col = layer_colors[layer_names],
       lty = 1,
       lwd = legend_line_thickness,
       ncol = 2,
       bty = "n",
       cex = legend_text_relative_font_size,
       text.col = "black",
       text.font = 1,
       x.intersp = 0.7,
       y.intersp = 0.9)
dev.off()


## Figure 3B
plot_width_cm <- 6.11
plot_height_cm <- 9.79
panel_heights <- c(1, 1, 1.1)
panel_gap <- 0.85
side_tick_label_gap <- 0.65
bottom_tick_label_gap <- 0.5
figure_name <- "Figure_3B.svg"

svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72))
max_k <- Polygonia_SOM$max_k
optim_k_vals <- as.numeric(Polygonia_SOM$optim_k_vals)
BIC_values <- Polygonia_SOM$BIC_values[seq_len(max_k), , drop = FALSE]
support_values <- if (!is.null(Polygonia_SOM$support_values)) Polygonia_SOM$support_values[seq_len(max_k), , drop = FALSE] else BIC_values
k_colors <- viridis::magma(max_k)
base_font_size <- par("ps")
bottom_numbers_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
layout(matrix(1:3, ncol = 1), heights = panel_heights)
par(bty = "n", oma = c(0, 0, 0, 0), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mar = c(panel_gap, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
finite_k_rows <- apply(support_values, 1, function(x) any(is.finite(x) & !is.na(x)))
plotted_k_values <- seq_len(max_k)[finite_k_rows]
values_for_plot <- t(support_values[finite_k_rows, , drop = FALSE])
boxplot(values_for_plot, at = plotted_k_values, xlim = c(0.5, max_k + 0.5), outline = FALSE, notch = FALSE, axes = FALSE, ylab = "", main = "", whisklty = 1, staplelty = 1, col = k_colors[plotted_k_values])
axis(1, at = seq_len(max_k), labels = FALSE)
y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = y_axis_breaks, labels = round(y_axis_breaks, 1), las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
par(mar = c(panel_gap, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
delta_BIC_matrix <- apply(BIC_values, 2, function(x) { previous_BIC <- x[-length(x)]; current_BIC <- x[-1]; delta_BIC <- previous_BIC - current_BIC; delta_BIC[!is.finite(previous_BIC) | !is.finite(current_BIC)] <- NA_real_; delta_BIC })
if (is.null(dim(delta_BIC_matrix))) delta_BIC_matrix <- matrix(delta_BIC_matrix, ncol = 1)
finite_delta_rows <- apply(delta_BIC_matrix, 1, function(x) any(is.finite(x) & !is.na(x)))
plotted_delta_k_values <- seq.int(2, max_k)[finite_delta_rows]
delta_BIC_for_plot <- t(delta_BIC_matrix[finite_delta_rows, , drop = FALSE])
boxplot(delta_BIC_for_plot, at = plotted_delta_k_values, xlim = c(0.5, max_k + 0.5), outline = FALSE, notch = FALSE, axes = FALSE, ylab = "", main = "", whisklty = 1, staplelty = 1, col = k_colors[plotted_delta_k_values])
axis(1, at = seq_len(max_k), labels = FALSE)
y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = y_axis_breaks, labels = round(y_axis_breaks, 1), las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
par(mar = c(2, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
k_frequency_values <- table(factor(optim_k_vals, levels = seq_len(max_k))) / length(optim_k_vals)
bar_midpoints <- barplot(k_frequency_values, ylim = c(0, 1), col = k_colors, axes = FALSE, axisnames = FALSE, ylab = "", main = "")
axis(1, at = bar_midpoints, labels = seq_len(max_k), mgp = c(3, bottom_tick_label_gap, 0), font.axis = 1, cex.axis = bottom_numbers_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
dev.off()


## Figure 3C & Figure 3D
plot_width_cm <- 12
plot_height_cm <- 4
panel_gap_cm <- 0.52
boundary_color_clusters <- "red"
boundary_line_width_clusters <- 3
point_color_clusters <- "white"
point_shape_clusters <- 19
point_size_clusters <- 0.9
cluster_shape_clusters <- "straight"
cluster_shape_neighbor_dist <- "straight"
legend_font_size <- 7
figure_name <- "Figure_3CD.svg"

calc_unit_neighbor_dist <- function(som_model) {
  number_of_units <- nrow(som_model$grid$pts)
  codebook_distance_matrix <- as.matrix(kohonen::object.distances(som_model, type = "codes"))
  grid_distance_matrix <- as.matrix(kohonen::unit.distances(som_model$grid))
  neighbor_matrix <- abs(grid_distance_matrix - 1) <= 0.001
  codebook_distance_matrix[!neighbor_matrix] <- NA_real_
  unit_mean_neighbor_distances <- colMeans(codebook_distance_matrix, na.rm = TRUE)
  unit_mean_neighbor_distances[!is.finite(unit_mean_neighbor_distances)] <- NA_real_
  unit_mean_neighbor_distances
}
count_SOM_clusters <- function(cluster_vector) {
  cluster_vector <- as.integer(cluster_vector)
  cluster_vector <- cluster_vector[is.finite(cluster_vector) & !is.na(cluster_vector) & cluster_vector >= 1]
  if (length(cluster_vector) == 0) return(NA_integer_)
  length(unique(cluster_vector))
}
choose_representative_replicate <- function(som_models, som_clusters) {
  number_of_replicates <- length(som_clusters)
  if (number_of_replicates == 1) return(1L)
  sample_cluster_assignments <- vector("list", number_of_replicates)
  for (replicate_index in seq_len(number_of_replicates)) {
    unit_classif <- as.integer(som_models[[replicate_index]]$unit.classif)
    unit_cluster_labels <- as.integer(som_clusters[[replicate_index]])
    sample_cluster_assignments[[replicate_index]] <- unit_cluster_labels[unit_classif]
  }
  k_values <- vapply(som_clusters, count_SOM_clusters, integer(1))
  k_frequency <- table(k_values[!is.na(k_values)])
  modal_k_values <- as.integer(names(k_frequency)[k_frequency == max(k_frequency)])
  selected_k <- min(modal_k_values)
  candidate_replicates <- which(k_values == selected_k)
  if (length(candidate_replicates) == 1) return(candidate_replicates)
  pairwise_adjusted_rand_index <- matrix(NA_real_, nrow = length(candidate_replicates), ncol = length(candidate_replicates))
  diag(pairwise_adjusted_rand_index) <- NA_real_
  for (candidate_index_1 in seq_len(length(candidate_replicates) - 1)) {
    for (candidate_index_2 in seq.int(candidate_index_1 + 1, length(candidate_replicates))) {
      replicate_index_1 <- candidate_replicates[candidate_index_1]
      replicate_index_2 <- candidate_replicates[candidate_index_2]
      current_adjusted_rand_index <- mclust::adjustedRandIndex(sample_cluster_assignments[[replicate_index_1]], sample_cluster_assignments[[replicate_index_2]])
      pairwise_adjusted_rand_index[candidate_index_1, candidate_index_2] <- current_adjusted_rand_index
      pairwise_adjusted_rand_index[candidate_index_2, candidate_index_1] <- current_adjusted_rand_index
    }
  }
  mean_adjusted_rand_index <- rowMeans(pairwise_adjusted_rand_index, na.rm = TRUE)
  if (all(!is.finite(mean_adjusted_rand_index) | is.na(mean_adjusted_rand_index))) return(candidate_replicates[1])
  representative_candidate_index <- which.max(replace(mean_adjusted_rand_index, !is.finite(mean_adjusted_rand_index) | is.na(mean_adjusted_rand_index), -Inf))
  candidate_replicates[representative_candidate_index]
}
som_models_use <- Polygonia_SOM$som_models
som_clusters_use <- Polygonia_SOM$som_clusters
replicate_index <- choose_representative_replicate(som_models = som_models_use, som_clusters = som_clusters_use)
som_model <- som_models_use[[replicate_index]]
som_cluster <- as.integer(som_clusters_use[[replicate_index]])
neighbor_distances <- calc_unit_neighbor_dist(som_model)
som_cluster <- as.integer(factor(som_cluster, levels = sort(unique(som_cluster))))
number_of_clusters <- length(unique(som_cluster))
cluster_colors <- viridis::viridis(number_of_clusters)
SOM_cluster_plot_col <- cluster_colors[som_cluster]
device_width_in <- (plot_width_cm / 2.54) * (96 / 72)
device_height_in <- (plot_height_cm / 2.54) * (96 / 72)
measurement_file <- tempfile(fileext = ".svg")
svg(measurement_file, width = device_width_in, height = device_height_in, family = "Arial")
base_font_size <- par("ps")
legend_relative_font_size <- (legend_font_size * (96 / 72)) / base_font_size
line_height_in <- par("csi") * par("mex")
panel_gap_in <- (panel_gap_cm / 2.54) * (96 / 72)
neighbor_panel_margin_width_in <- (4.6 + 0.6) * line_height_in
cluster_panel_margin_width_in <- (0.6 + 0.6) * line_height_in
map_plot_width_in <- (device_width_in - panel_gap_in - neighbor_panel_margin_width_in - cluster_panel_margin_width_in) / 2
neighbor_panel_width_in <- map_plot_width_in + neighbor_panel_margin_width_in
neighbor_panel_end_initial <- neighbor_panel_width_in / device_width_in
cluster_panel_start_initial <- (neighbor_panel_width_in + panel_gap_in) / device_width_in
par(family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(fig = c(0, neighbor_panel_end_initial, 0, 1), new = FALSE)
plot(x = som_model, type = "property", property = neighbor_distances, main = "", shape = cluster_shape_neighbor_dist, cex = legend_relative_font_size, palette.name = function(n) rev(viridis::cividis(n)))
left_map_right <- grconvertX(max(som_model$grid$pts[, 1]) + 0.5, from = "user", to = "ndc")
par(fig = c(cluster_panel_start_initial, 1, 0, 1), new = TRUE)
plot(x = som_model, shape = cluster_shape_clusters, type = "mapping", bgcol = SOM_cluster_plot_col, main = "", pchs = point_shape_clusters, cex = point_size_clusters, col = point_color_clusters)
right_map_left <- grconvertX(min(som_model$grid$pts[, 1]) - 0.5, from = "user", to = "ndc")
dev.off()
unlink(measurement_file)
actual_gap_fraction <- right_map_left - left_map_right
desired_gap_fraction <- panel_gap_cm / plot_width_cm
panel_shift_fraction <- (actual_gap_fraction - desired_gap_fraction) / 2
neighbor_panel_start <- panel_shift_fraction
neighbor_panel_end <- neighbor_panel_end_initial + panel_shift_fraction
cluster_panel_start <- cluster_panel_start_initial - panel_shift_fraction
cluster_panel_end <- 1 - panel_shift_fraction
svg(file.path(figure_files_folder, figure_name), width = device_width_in, height = device_height_in, family = "Arial")
base_font_size <- par("ps")
legend_relative_font_size <- (legend_font_size * (96 / 72)) / base_font_size
par(family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(fig = c(neighbor_panel_start, neighbor_panel_end, 0, 1), new = FALSE)
plot(x = som_model, type = "property", property = neighbor_distances, main = "", shape = cluster_shape_neighbor_dist, cex = legend_relative_font_size, palette.name = function(n) rev(viridis::cividis(n)))
par(fig = c(cluster_panel_start, cluster_panel_end, 0, 1), new = TRUE)
plot(x = som_model, shape = cluster_shape_clusters, type = "mapping", bgcol = SOM_cluster_plot_col, main = "", pchs = point_shape_clusters, cex = point_size_clusters, col = point_color_clusters)
if (number_of_clusters > 1) kohonen::add.cluster.boundaries(som_model, som_cluster, lwd = boundary_line_width_clusters, col = boundary_color_clusters)
dev.off()


## Figure 3E
plot_width_cm <- 15.36
plot_height_cm <- 4.8
axis.ticks.font.size <- 7.1
side_tick_label_gap <- 0.6
axis_bar_gap <- 4.3
font_family <- "Arial"
font_color <- "black"
plot_file_name <- "Figure_3E.svg"

Polygonia_species <- Polygonia_metadata[rownames(Polygonia_SOM$ancestry_matrix), "Species"]
Polygonia_order <- order(as.character(Polygonia_species))
Polygonia_ancestry_plot <- Polygonia_SOM$ancestry_matrix[Polygonia_order, , drop = FALSE]
Polygonia_species_plot <- Polygonia_species[Polygonia_order]
svg_scaling_factor <- 96 / 72
cluster_colors <- viridis::viridis(ncol(Polygonia_ancestry_plot))
plotting_assignment_coefficients <- apply(cbind(0, Polygonia_ancestry_plot), 1, cumsum)
svg(file.path(figure_files_folder, plot_file_name),
    width = (plot_width_cm / 2.54) * svg_scaling_factor,
    height = (plot_height_cm / 2.54) * svg_scaling_factor,
    family = font_family)
base_font_size <- par("ps")
axis_ticks_relative_font_size <- (axis.ticks.font.size * svg_scaling_factor) / base_font_size
par(mar = c(2, 2, 1.5, 1.5),
    family = font_family,
    fg = font_color,
    col.axis = font_color,
    col.lab = font_color,
    col.main = font_color,
    bty = "n",
    cex = 1,
    cex.axis = 1,
    cex.lab = 1,
    cex.main = 1)
plot(0,
     xlim = c(-axis_bar_gap, nrow(Polygonia_ancestry_plot)),
     ylim = c(0, 1),
     type = "n",
     ylab = "",
     xlab = "",
     xaxt = "n",
     yaxt = "n",
     xaxs = "i",
     yaxs = "i",
     frame.plot = FALSE)
axis(side = 2,
     las = 3,
     mgp = c(3, side_tick_label_gap, 0),
     col = font_color,
     col.axis = font_color,
     cex.axis = axis_ticks_relative_font_size)
for (cluster_index in seq_len(ncol(Polygonia_ancestry_plot))) {
  for (individual_index in seq_len(nrow(Polygonia_ancestry_plot))) {
    polygon(x = c(individual_index - 1, individual_index, individual_index, individual_index - 1),
            y = c(plotting_assignment_coefficients[cluster_index, individual_index],
                  plotting_assignment_coefficients[cluster_index, individual_index],
                  plotting_assignment_coefficients[cluster_index + 1, individual_index],
                  plotting_assignment_coefficients[cluster_index + 1, individual_index]),
            col = cluster_colors[cluster_index],
            border = cluster_colors[cluster_index],
            lwd = 0.5)
  }
}
dev.off()



#### Figure 4 ##################################################################

## Install packages
rm(list = ls()) #clear environment
CRAN_packages <- c("tibble", "dplyr", "vcfR") #CRAN packages
Bioconductor_packages <- "SeqArray" #Bioconductor package
for(p in CRAN_packages) if(!requireNamespace(p, quietly = TRUE)) install.packages(p) #install missing CRAN packages
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") #install BiocManager
if(!requireNamespace(Bioconductor_packages, quietly = TRUE)) BiocManager::install(Bioconductor_packages, ask = FALSE, update = FALSE) #install missing Bioconductor package
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv") #import (NicheDiv::transform.skewed.variables function from my other package (under development)
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")


## Set directories
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
intermediate_files_folder <- "Empirical_examples/Intermediate_files"
figure_files_folder <- "Figure_files"


## Import and prepare data
Polygonia_SNP <- process.SNP.data.SOM(vcf.path = "./Empirical_examples/Dupuis_et_al_2018/Polygonia_961SNPs.vcf", #filter loci and individuals and create SNP matrix dataframe
                                      missing.loci.cutoff.lenient = 0.7,
                                      missing.loci.cutoff.final = 0.5,
                                      missing.individuals.cutoff = 0.5)
rownames(Polygonia_SNP) <- sub(".*?(\\d+)$", "\\1", rownames(Polygonia_SNP)) #only keep numeric identifier as rownames
Polygonia_COI <- process.SNP.data.SOM(nexus.path = "./Empirical_examples/Dupuis_et_al_2018/Polygonia_COI.nex",
                                      missing.loci.cutoff.lenient = 0.7, 
                                      missing.loci.cutoff.final = 0.5,
                                      missing.individuals.cutoff = 0.5)
Polygonia_COI_numeric_rownames <- sub(".*?(\\d+)$", "\\1", rownames(Polygonia_COI)) #extract numeric code from each rowname (e.g., "pf_8301" -> "8301")
Polygonia_COI <- Polygonia_COI[!duplicated(Polygonia_COI_numeric_rownames), , drop = FALSE] #keep only first occurrence for each numeric code (remove duplicates)
rownames(Polygonia_COI) <- Polygonia_COI_numeric_rownames[!duplicated(Polygonia_COI_numeric_rownames)] #set rownames to unique numeric codes
Polygonia_RGB <- read.delim("./Empirical_examples/Dupuis_et_al_2018/Polygonia_RGB_characters.txt", stringsAsFactors = FALSE)
rownames(Polygonia_RGB) <- Polygonia_RGB$Species
Polygonia_RGB <- magrittr::`%>%`(Polygonia_RGB, dplyr::select(-Name, -Species)) #remove columns
Polygonia_wing_scores <- read.delim("./Empirical_examples/Dupuis_et_al_2018/Polygonia_visually_scored.txt", stringsAsFactors = FALSE)
rownames(Polygonia_wing_scores) <- Polygonia_wing_scores$Name
Polygonia_wing_scores <- Polygonia_wing_scores |>
  dplyr::select(-Name, -Species) |> #remove columns
  dplyr::rename(Wing_character_1 = Ch1,
                Wing_character_2 = Ch2,
                Wing_character_3 = Ch3,
                Wing_character_4 = Ch4,
                Wing_character_5 = Ch5,
                Wing_character_6 = Ch6,
                Wing_character_7 = Ch7,
                Wing_character_8 = Ch8,
                Wing_character_9 = Ch9,
                Wing_character_10 = Ch10)
Polygonia_wing_scores$Wing_character_8 <- factor(Polygonia_wing_scores$Wing_character_8, levels = 1:4)
Wing_character_8_states <- stats::model.matrix(~ Wing_character_8 - 1, data = Polygonia_wing_scores) #make wing character 8 binary (since it is not ordinal)
colnames(Wing_character_8_states) <- paste0("Wing_character_8_state_", 1:4)
Polygonia_wing_scores <- cbind(Polygonia_wing_scores, Wing_character_8_states)
Polygonia_wing_scores$Wing_character_8 <- NULL
Polygonia_metadata <- read.csv("./Empirical_examples/Dupuis_et_al_2018/Polygonia_metadata.csv", header = TRUE, sep = ";")
rownames(Polygonia_metadata) <- Polygonia_metadata$ID

Polygonia_spatial <- dplyr::select(Polygonia_metadata, Latitude, Longitude) #create dataframe with Lat and Long
Polygonia_spatial$Elevation <- NA #initialize elevation column with NA
Polygonia_spatial_sf <- sf::st_as_sf(Polygonia_spatial[!is.na(Polygonia_spatial$Latitude) & !is.na(Polygonia_spatial$Longitude), ], coords = c("Longitude", "Latitude"), crs = 4326) #extract elevation
Polygonia_spatial$Elevation[!is.na(Polygonia_spatial$Latitude) & !is.na(Polygonia_spatial$Longitude)] <-
  elevatr::get_elev_point(locations = Polygonia_spatial_sf, prj = sf::st_crs(Polygonia_spatial_sf)$proj4string, src = "aws")$elevation
Polygonia_morphotype <- dplyr::select(Polygonia_metadata, Morphotype)
Polygonia_metadata <- dplyr::select(Polygonia_metadata, Species, ID)
rownames(Polygonia_RGB) <- as.character(rownames(Polygonia_RGB))
rownames(Polygonia_wing_scores) <- as.character(rownames(Polygonia_wing_scores))
rownames(Polygonia_morphotype) <- as.character(rownames(Polygonia_morphotype))
Polygonia_RGB_wing_scores <- merge(Polygonia_RGB,
                                   Polygonia_wing_scores,
                                   by = "row.names",
                                   all = FALSE)
rownames(Polygonia_RGB_wing_scores) <- Polygonia_RGB_wing_scores$Row.names
Polygonia_RGB_wing_scores$Row.names <- NULL
Polygonia_morphology <- merge(Polygonia_RGB_wing_scores,
                              Polygonia_morphotype,
                              by = "row.names",
                              all = FALSE)
rownames(Polygonia_morphology) <- Polygonia_morphology$Row.names
Polygonia_morphology$Row.names <- NULL
Polygonia_morphology <- make.cols.binary.SOM(Polygonia_morphology, #convert Morphotype to binary columns and remove original
                                             make.binary.cols = "Morphotype",
                                             append.to.original = TRUE)
Polygonia_morphology$Morphotype <- NULL
non.continuous.cols <- grepl("^Wing_character_|^Morphotype_", colnames(Polygonia_morphology)) #identify non-continuous traits
Polygonia_morphology_categorical <- Polygonia_morphology[, non.continuous.cols, drop = FALSE] #non-continuous traits
Polygonia_morphology <- Polygonia_morphology[, !non.continuous.cols, drop = FALSE] #continuous traits
Polygonia_morphology <- remove.lowCV.multicollinearity.SOM(Polygonia_morphology, #remove highly correlated and low-variance variables
                                                           CV.threshold = 0.05,
                                                           cor.threshold = 0.9)
Polygonia_environmental <- read.csv("./Empirical_examples/Dupuis_et_al_2018/Polygonia_environmental.csv",
                                    row.names = 1, header = TRUE)
Polygonia_environmental <- dplyr::select(Polygonia_environmental, -Latitude, -Longitude, -Elevation)
Polygonia_environmental_rownames <- rownames(Polygonia_environmental) #save rownames
Polygonia_environmental <- as.data.frame(lapply(Polygonia_environmental, as.numeric)) #ensure all columns are numeric
rownames(Polygonia_environmental) <- Polygonia_environmental_rownames #reassign saved row names
Polygonia_environmental <- (NicheDiv::transform.skewed.variables(Polygonia_environmental))$transformed #transform skewed variables
Polygonia_environmental <- remove.lowCV.multicollinearity.SOM(Polygonia_environmental, #remove highly correlated and low-variance variables
                                                              CV.threshold = 0.05,
                                                              cor.threshold = 0.9)
for (Polygonia_shared_data in c("Polygonia_morphology", "Polygonia_morphology_categorical", "Polygonia_SNP", "Polygonia_COI", "Polygonia_spatial", "Polygonia_environmental", "Polygonia_metadata")) {
  Polygonia_shared_data_mat <- get(Polygonia_shared_data)
  rownames(Polygonia_shared_data_mat) <- make.unique(as.character(rownames(Polygonia_shared_data_mat)))
  assign(Polygonia_shared_data, Polygonia_shared_data_mat, envir = .GlobalEnv)}
Polygonia_common_IDs <- Reduce(intersect, list(
  rownames(Polygonia_morphology),
  rownames(Polygonia_morphology_categorical),
  rownames(Polygonia_SNP),
  rownames(Polygonia_COI),
  rownames(Polygonia_spatial),
  rownames(Polygonia_environmental),
  rownames(Polygonia_metadata)))
Polygonia_morphology2 <- Polygonia_morphology[Polygonia_common_IDs, , drop = FALSE]
Polygonia_morphology_categorical2 <- Polygonia_morphology_categorical[Polygonia_common_IDs, , drop = FALSE]
Polygonia_SNP2 <- Polygonia_SNP[Polygonia_common_IDs, , drop = FALSE]
Polygonia_COI2 <- Polygonia_COI[Polygonia_common_IDs, , drop = FALSE]
Polygonia_spatial2 <- Polygonia_spatial[Polygonia_common_IDs, , drop = FALSE]
Polygonia_environmental2 <- Polygonia_environmental[Polygonia_common_IDs, , drop = FALSE]
Polygonia_metadata2 <- Polygonia_metadata[Polygonia_common_IDs, , drop = FALSE]
Polygonia_morphology <- Polygonia_morphology2[rowSums(!is.na(Polygonia_morphology2)) > 0, , drop = FALSE]
Polygonia_morphology_categorical <- Polygonia_morphology_categorical2[rowSums(!is.na(Polygonia_morphology_categorical2)) > 0, , drop = FALSE]
Polygonia_SNP <- Polygonia_SNP2[rowSums(!is.na(Polygonia_SNP2)) > 0, , drop = FALSE]
Polygonia_COI <- Polygonia_COI2[rowSums(!is.na(Polygonia_COI2)) > 0, , drop = FALSE]
Polygonia_spatial <- Polygonia_spatial2[rowSums(!is.na(Polygonia_spatial2)) > 0, , drop = FALSE]
Polygonia_environmental <- Polygonia_environmental2[rowSums(!is.na(Polygonia_environmental2)) > 0, , drop = FALSE]
Polygonia_metadata <- Polygonia_metadata2[rowSums(!is.na(Polygonia_metadata2)) > 0, , drop = FALSE]
Polygonia_species_vec <- Polygonia_metadata[rownames(Polygonia_morphology), "Species"]
Polygonia_new_rownames <- paste(rownames(Polygonia_morphology), Polygonia_species_vec, sep = "_")
rownames(Polygonia_morphology) <- Polygonia_new_rownames
rownames(Polygonia_morphology_categorical) <- Polygonia_new_rownames
rownames(Polygonia_SNP) <- Polygonia_new_rownames
rownames(Polygonia_COI) <- Polygonia_new_rownames
rownames(Polygonia_spatial) <- Polygonia_new_rownames
rownames(Polygonia_metadata) <- Polygonia_new_rownames
rownames(Polygonia_environmental) <- Polygonia_new_rownames
Polygonia_all_data <- list(Morphology = Polygonia_morphology,
                           Morphology_2 = Polygonia_morphology_categorical,
                           SNP = Polygonia_SNP,
                           COI = Polygonia_COI,
                           Environmental = Polygonia_environmental,
                           Spatial = Polygonia_spatial)
Polygonia_SOM_tr <- train.SOM(input_data = Polygonia_all_data, #200 samples, 4.1min
                              save.SOM.results = TRUE,
                              save.SOM.results.name = file.path(intermediate_files_folder, "Polygonia_SOM_tr.Rdata"),
                              max.NA.row = 0.5,
                              max.NA.col = 0.5)
Polygonia_SOM <- clustering.SOM(Polygonia_SOM_tr, #3.0min
                                clustering.method = "kmeans+BICelbow",
                                save.SOM.results.name = file.path(intermediate_files_folder, "Polygonia_SOM_kmeansBICelbow.Rdata"))


## Figure 4A
plot_width_cm <- 15.88
plot_height_cm <- 11.51
row_gap <- 1.45
column_gap <- 3.5
bottom_tick_label_gap <- 0.6
top_margin_mm <- 0
left_margin_mm <- 3.4
right_margin_mm <- 0
layer_label_font_size <- 9
bar_label_font_size <- 7.1
axis_ticks_font_size <- 7.1
bars_threshold_N <- 20
importance_threshold <- 0.001
figure_name <- "Figure_4A.svg"
variable_label_abbreviations <- c("Wing_character_1" = "WC1",
                                  "Wing_character_2" = "WC2",
                                  "Wing_character_3" = "WC3",
                                  "Wing_character_4" = "WC4",
                                  "Wing_character_5" = "WC5",
                                  "Wing_character_6" = "WC6",
                                  "Wing_character_7" = "WC7",
                                  "Wing_character_9" = "WC9",
                                  "Wing_character_10" = "WC10",
                                  "Wing_character_8_state_1" = "WC8-1",
                                  "Wing_character_8_state_2" = "WC8-2",
                                  "Wing_character_8_state_3" = "WC8-3",
                                  "Wing_character_8_state_4" = "WC8-4",
                                  "Morphotype_Contrasted" = "MT cont",
                                  "Morphotype_Smeared" = "MT smear",
                                  "Latitude" = "Lat",
                                  "Longitude" = "Long",
                                  "Elevation" = "Elev")
calculate.etasquared.per.variable <- function(codebook_matrix, neuron_cluster_vector, som_model, baseline_weight = 1, calculation_block_size = 10000) {
  codebook_matrix <- as.matrix(codebook_matrix)
  storage.mode(codebook_matrix) <- "numeric"
  neuron_cluster_vector <- as.integer(neuron_cluster_vector)
  valid_neuron_rows <- is.finite(neuron_cluster_vector) & !is.na(neuron_cluster_vector)
  number_of_SOM_units <- length(neuron_cluster_vector)
  neuron_sample_counts <- tabulate(som_model$unit.classif, nbins = number_of_SOM_units)
  codebook_matrix <- codebook_matrix[valid_neuron_rows, , drop = FALSE]
  neuron_cluster_vector <- neuron_cluster_vector[valid_neuron_rows]
  neuron_sample_counts <- neuron_sample_counts[valid_neuron_rows]
  neuron_weights <- neuron_sample_counts + baseline_weight
  number_of_variables <- ncol(codebook_matrix)
  etasquared_values <- rep(NA_real_, number_of_variables)
  names(etasquared_values) <- colnames(codebook_matrix)
  if (nrow(codebook_matrix) < 2 || length(unique(neuron_cluster_vector)) < 2 || !is.finite(sum(neuron_weights)) || sum(neuron_weights) <= 0) return(etasquared_values)
  total_neuron_weight <- sum(neuron_weights)
  cluster_weight_sums <- rowsum(matrix(neuron_weights, ncol = 1), group = neuron_cluster_vector, reorder = FALSE)[, 1]
  block_start_indices <- seq.int(from = 1, to = number_of_variables, by = calculation_block_size)
  for (block_start_index in block_start_indices) {
    block_end_index <- min(block_start_index + calculation_block_size - 1, number_of_variables)
    block_variable_indices <- block_start_index:block_end_index
    codebook_block <- codebook_matrix[, block_variable_indices, drop = FALSE]
    finite_block_variable_positions <- which(colSums(!is.finite(codebook_block)) == 0)
    if (length(finite_block_variable_positions) > 0) {
      finite_codebook_block <- codebook_block[, finite_block_variable_positions, drop = FALSE]
      weighted_codebook_block <- finite_codebook_block * neuron_weights
      weighted_variable_sums <- colSums(weighted_codebook_block)
      weighted_squared_variable_sums <- colSums(finite_codebook_block * weighted_codebook_block)
      total_sum_of_squares <- weighted_squared_variable_sums - weighted_variable_sums^2 / total_neuron_weight
      cluster_weighted_variable_sums <- rowsum(weighted_codebook_block, group = neuron_cluster_vector, reorder = FALSE)
      cluster_mean_components <- sweep(cluster_weighted_variable_sums^2, MARGIN = 1, STATS = cluster_weight_sums, FUN = "/")
      between_cluster_sum_of_squares <- colSums(cluster_mean_components) - weighted_variable_sums^2 / total_neuron_weight
      finite_etasquared_values <- rep(NA_real_, length(finite_block_variable_positions))
      degenerate_variables <- is.finite(total_sum_of_squares) & total_sum_of_squares <= 0
      valid_variables <- is.finite(total_sum_of_squares) & total_sum_of_squares > 0 & is.finite(between_cluster_sum_of_squares)
      finite_etasquared_values[degenerate_variables] <- 0
      finite_etasquared_values[valid_variables] <- between_cluster_sum_of_squares[valid_variables] / total_sum_of_squares[valid_variables]
      numerical_tolerance <- sqrt(.Machine$double.eps)
      finite_etasquared_values[finite_etasquared_values < 0 & finite_etasquared_values > -numerical_tolerance] <- 0
      finite_etasquared_values[finite_etasquared_values > 1 & finite_etasquared_values < 1 + numerical_tolerance] <- 1
      complete_variable_indices <- block_variable_indices[finite_block_variable_positions]
      etasquared_values[complete_variable_indices] <- finite_etasquared_values
    }
    incomplete_block_variable_positions <- which(colSums(!is.finite(codebook_block)) > 0)
    if (length(incomplete_block_variable_positions) > 0) {
      incomplete_variable_indices <- block_variable_indices[incomplete_block_variable_positions]
      etasquared_values[incomplete_variable_indices] <- vapply(incomplete_variable_indices, function(variable_index) {
        variable_values <- codebook_matrix[, variable_index]
        valid_variable_rows <- is.finite(variable_values) & !is.na(variable_values) & is.finite(neuron_weights) & !is.na(neuron_weights)
        variable_values <- variable_values[valid_variable_rows]
        cluster_labels <- neuron_cluster_vector[valid_variable_rows]
        variable_weights <- neuron_weights[valid_variable_rows]
        if (length(variable_values) < 2) return(NA_real_)
        if (length(unique(cluster_labels)) < 2) return(NA_real_)
        if (sum(variable_weights) <= 0) return(NA_real_)
        weighted_grand_mean <- sum(variable_weights * variable_values) / sum(variable_weights)
        total_sum_of_squares <- sum(variable_weights * (variable_values - weighted_grand_mean)^2)
        if (!is.finite(total_sum_of_squares) || total_sum_of_squares <= 0) return(0)
        cluster_weighted_sums <- tapply(variable_weights * variable_values, cluster_labels, sum)
        cluster_weight_sums_current <- tapply(variable_weights, cluster_labels, sum)
        cluster_weighted_means <- cluster_weighted_sums / cluster_weight_sums_current
        between_cluster_sum_of_squares <- sum(cluster_weight_sums_current * (cluster_weighted_means - weighted_grand_mean)^2)
        as.numeric(between_cluster_sum_of_squares / total_sum_of_squares)
      }, numeric(1))
    }
  }
  etasquared_values
}
matrix_names <- Polygonia_SOM$input_data_names
first_codebook_list <- kohonen::getCodes(Polygonia_SOM$som_models[[1]])
if (!is.list(first_codebook_list)) first_codebook_list <- list(first_codebook_list)
number_of_layers <- length(first_codebook_list)
if (length(matrix_names) != number_of_layers) matrix_names <- paste0("layer", seq_len(number_of_layers))
matrix_names_plot <- matrix_names
matrix_names_plot[matrix_names_plot == "Morphology_2"] <- "Morphology 2"
replicate_k_values <- vapply(Polygonia_SOM$som_clusters, function(cluster_vector) length(unique(cluster_vector[is.finite(cluster_vector) & !is.na(cluster_vector)])), integer(1))
retained_replicate_indices <- which(is.finite(replicate_k_values) & !is.na(replicate_k_values) & replicate_k_values >= 2L)
first_retained_codebook_list <- kohonen::getCodes(Polygonia_SOM$som_models[[retained_replicate_indices[1]]])
if (!is.list(first_retained_codebook_list)) first_retained_codebook_list <- list(first_retained_codebook_list)
for (layer_index in seq_len(number_of_layers)) {
  if (is.null(colnames(first_retained_codebook_list[[layer_index]]))) colnames(first_retained_codebook_list[[layer_index]]) <- paste0("V", seq_len(ncol(first_retained_codebook_list[[layer_index]])))
}
Polygonia_SOM_variable_importance <- vector("list", number_of_layers)
names(Polygonia_SOM_variable_importance) <- matrix_names
for (layer_index in seq_len(number_of_layers)) {
  Polygonia_SOM_variable_importance[[layer_index]] <- matrix(NA_real_, nrow = length(retained_replicate_indices), ncol = ncol(first_retained_codebook_list[[layer_index]]), dimnames = list(paste0("R", retained_replicate_indices), colnames(first_retained_codebook_list[[layer_index]])))
}
for (retained_replicate_position in seq_along(retained_replicate_indices)) {
  retained_replicate_index <- retained_replicate_indices[retained_replicate_position]
  som_model <- Polygonia_SOM$som_models[[retained_replicate_index]]
  codebook_list <- kohonen::getCodes(som_model)
  if (!is.list(codebook_list)) codebook_list <- list(codebook_list)
  som_cluster_vector <- Polygonia_SOM$som_clusters[[retained_replicate_index]]
  for (layer_index in seq_len(number_of_layers)) {
    if (is.null(colnames(codebook_list[[layer_index]]))) colnames(codebook_list[[layer_index]]) <- paste0("V", seq_len(ncol(codebook_list[[layer_index]])))
    Polygonia_SOM_variable_importance[[layer_index]][retained_replicate_position, ] <- calculate.etasquared.per.variable(codebook_matrix = codebook_list[[layer_index]], neuron_cluster_vector = som_cluster_vector, som_model = som_model, calculation_block_size = 10000)
  }
}
svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72), family = "Arial")
base_font_size <- par("ps")
layer_label_relative_font_size <- (layer_label_font_size * (96 / 72)) / base_font_size
bar_label_relative_font_size <- (bar_label_font_size * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (axis_ticks_font_size * (96 / 72)) / base_font_size
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
left_margin_lines <- (left_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
par(mfrow = if (number_of_layers <= 3) c(1, number_of_layers) else if (number_of_layers == 4) c(2, 2) else if (number_of_layers <= 6) c(2, 3) else if (number_of_layers <= 8) c(2, 4) else if (number_of_layers == 9) c(3, 3) else c(ceiling(number_of_layers / 3), 3), oma = c(0, left_margin_lines, top_margin_lines, right_margin_lines), mar = c(2, column_gap / 2, row_gap, column_gap / 2), mgp = c(3, bottom_tick_label_gap, 0), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1, bty = "o", family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black")
layer_colors <- setNames(viridis::turbo(length(matrix_names)), matrix_names)
for (layer_index in seq_along(Polygonia_SOM_variable_importance)) {
  variable_importance_matrix <- Polygonia_SOM_variable_importance[[layer_index]]
  if (is.null(variable_importance_matrix) || nrow(variable_importance_matrix) == 0 || ncol(variable_importance_matrix) == 0) {
    plot.new()
    next
  }
  median_metric_per_variable <- matrixStats::colMedians(variable_importance_matrix, na.rm = TRUE, useNames = TRUE)
  median_metric_per_variable[!is.finite(median_metric_per_variable)] <- NA_real_
  retained_variable_indices <- which(is.finite(median_metric_per_variable) & !is.na(median_metric_per_variable) & median_metric_per_variable > importance_threshold)
  if (length(retained_variable_indices) == 0) {
    plot.new()
    next
  }
  variable_importance_matrix <- variable_importance_matrix[, retained_variable_indices, drop = FALSE]
  retained_variable_medians <- median_metric_per_variable[colnames(variable_importance_matrix)]
  variable_sort_indices <- order(retained_variable_medians, decreasing = FALSE)
  variable_importance_matrix <- variable_importance_matrix[, variable_sort_indices, drop = FALSE]
  variable_labels <- colnames(variable_importance_matrix)
  abbreviation_matches <- match(variable_labels, names(variable_label_abbreviations))
  variable_labels[!is.na(abbreviation_matches)] <- variable_label_abbreviations[abbreviation_matches[!is.na(abbreviation_matches)]]
  number_of_bars <- ncol(variable_importance_matrix)
  boxplot(variable_importance_matrix, horizontal = TRUE, las = 1, notch = FALSE, outline = FALSE, col = rep(layer_colors[matrix_names[layer_index]], number_of_bars), axes = FALSE, whisklty = if (number_of_bars > bars_threshold_N) 0 else 1, staplelty = if (number_of_bars > bars_threshold_N) 0 else 1, names = rep("", number_of_bars))
  axis(1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size, col = "black", col.axis = "black", font.axis = 1)
  axis(2, at = seq_len(number_of_bars), labels = if (number_of_bars <= bars_threshold_N) variable_labels else FALSE, las = 1, tick = FALSE, cex.axis = bar_label_relative_font_size, col.axis = "black", font.axis = 1)
  mtext(matrix_names_plot[layer_index], side = 3, line = 0.3, cex = layer_label_relative_font_size, font = 1, family = "Arial", col = "black")
}
dev.off()


## Figure 4B
plot_width_cm <- 15.13
plot_height_cm <- 6.38
panel_gap <- 1.55
side_tick_label_gap <- 0.65
bottom_tick_label_gap <- 0.74
bottom_margin <- 5.5
top_margin_mm <- 1.5
right_margin_mm <- 1.5
layer_label_font_size <- 7.1
axis_ticks_font_size <- 7.1
point_cex <- 0.8
point_alpha <- 0.65
figure_name <- "Figure_4B.svg"
load(file.path(intermediate_files_folder, "Polygonia_SOM_lolo.Rdata"))
layer_summary_table <- leave.one.layer.out.results$layer.summary
replicate_matched_results_table <- leave.one.layer.out.results$replicate.matched.results
successful_replicate_matched_results_table <- replicate_matched_results_table[is.na(replicate_matched_results_table$error) | replicate_matched_results_table$error == "", , drop = FALSE]
SOM_layer_names <- as.character(layer_summary_table$layer)
SOM_layer_names_plot <- SOM_layer_names
SOM_layer_names_plot[SOM_layer_names_plot == "Morphology_2"] <- "Morphology 2"
successful_replicate_matched_results_table$layer <- factor(successful_replicate_matched_results_table$layer, levels = SOM_layer_names)
layer_colors <- setNames(viridis::turbo(length(SOM_layer_names)), SOM_layer_names)
show_assignment_margin_plot <- any(is.finite(successful_replicate_matched_results_table$delta.mean.assignment.margin) & !is.na(successful_replicate_matched_results_table$delta.mean.assignment.margin))
svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72), family = "Arial")
base_font_size <- par("ps")
layer_label_relative_font_size <- (layer_label_font_size * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (axis_ticks_font_size * (96 / 72)) / base_font_size
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
layout(matrix(seq_len(if (show_assignment_margin_plot) 3 else 2), nrow = 1))
par(oma = c(0, 0, top_margin_lines, right_margin_lines), bty = "o", family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mar = c(bottom_margin, 2, 0.5, panel_gap / 2))
boxplot(absolute.k.deviation ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names_plot, las = 2, mgp = c(3, bottom_tick_label_gap, 0), font.axis = 1, cex.axis = layer_label_relative_font_size, col = "black", col.axis = "black")
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), font.axis = 1, cex.axis = axis_ticks_relative_font_size, col = "black", col.axis = "black")
for (layer_index in seq_along(SOM_layer_names)) {
  current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == SOM_layer_names[layer_index], "absolute.k.deviation"]
  current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
  if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = point_cex, col = adjustcolor(layer_colors[SOM_layer_names[layer_index]], alpha.f = point_alpha))
}
box()
par(mar = c(bottom_margin, 2 + panel_gap / 2, 0.5, panel_gap / 2))
boxplot(pairwise.coassignment.change ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names_plot, las = 2, mgp = c(3, bottom_tick_label_gap, 0), font.axis = 1, cex.axis = layer_label_relative_font_size, col = "black", col.axis = "black")
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), font.axis = 1, cex.axis = axis_ticks_relative_font_size, col = "black", col.axis = "black")
for (layer_index in seq_along(SOM_layer_names)) {
  current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == SOM_layer_names[layer_index], "pairwise.coassignment.change"]
  current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
  if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = point_cex, col = adjustcolor(layer_colors[SOM_layer_names[layer_index]], alpha.f = point_alpha))
}
box()
if (show_assignment_margin_plot) {
  par(mar = c(bottom_margin, 2 + panel_gap / 2, 0.5, 0.5))
  boxplot(delta.mean.assignment.margin ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
  axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names_plot, las = 2, mgp = c(3, bottom_tick_label_gap, 0), font.axis = 1, cex.axis = layer_label_relative_font_size, col = "black", col.axis = "black")
  axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), font.axis = 1, cex.axis = axis_ticks_relative_font_size, col = "black", col.axis = "black")
  for (layer_index in seq_along(SOM_layer_names)) {
    current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == SOM_layer_names[layer_index], "delta.mean.assignment.margin"]
    current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
    if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = point_cex, col = adjustcolor(layer_colors[SOM_layer_names[layer_index]], alpha.f = point_alpha))
  }
  box()
}
dev.off()



#### Figure 5 ##################################################################

## Install packages
rm(list = ls()) #clear environment
CRAN_packages <- c("tibble", "dplyr", "vcfR") #CRAN packages
Bioconductor_packages <- "SeqArray" #Bioconductor package
for(p in CRAN_packages) if(!requireNamespace(p, quietly = TRUE)) install.packages(p) #install missing CRAN packages
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") #install BiocManager
if(!requireNamespace(Bioconductor_packages, quietly = TRUE)) BiocManager::install(Bioconductor_packages, ask = FALSE, update = FALSE) #install missing Bioconductor package
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv") #import (NicheDiv::transform.skewed.variables function from my other package (under development)
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")


## Set directories
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
intermediate_files_folder <- "Empirical_examples/Intermediate_files"
figure_files_folder <- "Figure_files"


## Import and prepare data
Aeneus_data <- read.csv(file = "./Empirical_examples/Pyron_et_al_2024/aeneus56.csv",
                        row.names = 1,
                        header = TRUE, 
                        colClasses = c(huc2 = "character", huc4 = "character", huc6 = "character", huc8 = "character", huc10 = "character", huc12 = "character"))
Aeneus_SNP <- process.SNP.data.SOM(vcf.path = "./Empirical_examples/Pyron_et_al_2024/aeneus56.vcf.gz", #filter loci and individuals and create SNP matrix dataframe
                                   missing.loci.cutoff.lenient = 0.7,
                                   missing.loci.cutoff.final = 0.5,
                                   missing.individuals.cutoff = 0.5)
rownames(Aeneus_SNP) <- Aeneus_data$Sample[match(rownames(Aeneus_SNP), rownames(Aeneus_data))] #rename alleles
Aeneus_spatial <- Aeneus_data[, c("Lat", "Long", "Elev")] #extract coordinates and elevation
Aeneus_spatial <- dplyr::rename(Aeneus_spatial, Latitude = Lat, Longitude = Long, Elevation = Elev) #rename variables
rownames(Aeneus_spatial) <- Aeneus_data$Sample #assign rownames
Aeneus_environmental <- read.csv("./Empirical_examples/Pyron_et_al_2024/Aeneus56_environmental.csv", row.names = 1, header = TRUE) #read CSV with rownames
Aeneus_environmental <- Aeneus_environmental[, !names(Aeneus_environmental) %in% c("Latitude", "Longitude", "Elevation")] #remove spatial variables
Aeneus_environmental <- as.data.frame(lapply(Aeneus_environmental, as.numeric)) #ensure numeric
rownames(Aeneus_environmental) <- Aeneus_data$Sample #assign rownames
rownames(Aeneus_data) <- Aeneus_data$Sample #assign rownames
Aeneus_watershed <- make.cols.binary.SOM(dataframe = Aeneus_data, make.binary.cols = c("lvl4", "huc2", "huc4", "huc6", "huc8", "huc10", "huc12"))
Aeneus_watershed <- Aeneus_watershed[rownames(Aeneus_data), , drop = FALSE]
Aeneus_environmental <- (NicheDiv::transform.skewed.variables(Aeneus_environmental))$transformed #transform skewed variables
Aeneus_environmental <- remove.lowCV.multicollinearity.SOM(Aeneus_environmental,  CV.threshold = 0.05, cor.threshold = 0.9)
Aeneus_trait_names <- c("SVL", "TL", "AG", "CW", "FL", "HL", "SG", "TW", "TO", "FI", "HW", "ED", "IN", "ES", "ON", "IO", "IC") #define trait columns
Aeneus_trait_data <- Aeneus_data[, Aeneus_trait_names] #extract variables
rownames(Aeneus_trait_data) <- Aeneus_data$Sample #assign sample IDs as rownames
Aeneus_trait_data <- Aeneus_trait_data[rowSums(!is.na(Aeneus_trait_data)) > 0, ] #remove samples with only NA values
Aeneus_log_traits <- log(Aeneus_trait_data) #log-transform all traits
Aeneus_filtered_log_traits <- remove.lowCV.multicollinearity.SOM(Aeneus_log_traits, CV.threshold = 0.05, cor.threshold = 0.9, exclude.cols = "SVL")
rownames(Aeneus_filtered_log_traits) <- rownames(Aeneus_trait_data) #set rownames for filtered log traits
Aeneus_SVL <- Aeneus_filtered_log_traits[, "SVL"] #extract SVL and residualize all others
Aeneus_residuals_mat <- sapply(colnames(Aeneus_filtered_log_traits)[colnames(Aeneus_filtered_log_traits) != "SVL"], function(trait) stats::resid(stats::lm(Aeneus_filtered_log_traits[, trait] ~ Aeneus_SVL, na.action = stats::na.exclude))) #regress each trait on SVL and retain NA alignment
rownames(Aeneus_residuals_mat) <- rownames(Aeneus_filtered_log_traits) #assign rownames to residual matrix
Aeneus_morphology <- as.data.frame(cbind(SVL = Aeneus_SVL, Aeneus_residuals_mat)) #combine log(SVL) and residuals
Aeneus_morphology <- Aeneus_morphology[rownames(Aeneus_trait_data), ] #keep only samples with trait data
Aeneus_SOM_data <- list(Alleles = Aeneus_SNP,
                        Spatial = Aeneus_spatial,
                        Environmental = Aeneus_environmental,
                        Watershed = Aeneus_watershed,
                        Morphology = Aeneus_morphology)
Aeneus_SOM_tr <- train.SOM(input_data = Aeneus_SOM_data, #40 samples, 2.8min
                           save.SOM.results = TRUE,
                           save.SOM.results.name = file.path(intermediate_files_folder, "Aeneus_SOM_tr.Rdata"),
                           max.NA.row = 0.5,
                           max.NA.col = 0.5)
Aeneus_SOM <- clustering.SOM(Aeneus_SOM_tr, #3.9min
                             clustering.method = "kmeans+BICelbow",
                             save.SOM.results.name = file.path(intermediate_files_folder, "Aeneus_SOM_kmeansBICelbow.Rdata"))
Aeneus_SOM$optim_k_summary


## Figure 5A & Figure 5B
plot_width_cm <- 12.16
plot_height_cm <- 4
panel_gap_cm <- 0.52
boundary_color_clusters <- "red"
boundary_line_width_clusters <- 3
point_color_clusters <- "white"
point_shape_clusters <- 19
point_size_clusters <- 0.9
cluster_shape_clusters <- "straight"
cluster_shape_neighbor_dist <- "straight"
legend_font_size <- 7
figure_name <- "Figure_5AB.svg"

calc_unit_neighbor_dist <- function(som_model) {
  number_of_units <- nrow(som_model$grid$pts)
  codebook_distance_matrix <- as.matrix(kohonen::object.distances(som_model, type = "codes"))
  grid_distance_matrix <- as.matrix(kohonen::unit.distances(som_model$grid))
  neighbor_matrix <- abs(grid_distance_matrix - 1) <= 0.001
  codebook_distance_matrix[!neighbor_matrix] <- NA_real_
  unit_mean_neighbor_distances <- colMeans(codebook_distance_matrix, na.rm = TRUE)
  unit_mean_neighbor_distances[!is.finite(unit_mean_neighbor_distances)] <- NA_real_
  unit_mean_neighbor_distances
}
count_SOM_clusters <- function(cluster_vector) {
  cluster_vector <- as.integer(cluster_vector)
  cluster_vector <- cluster_vector[is.finite(cluster_vector) & !is.na(cluster_vector) & cluster_vector >= 1]
  if (length(cluster_vector) == 0) return(NA_integer_)
  length(unique(cluster_vector))
}
choose_representative_replicate <- function(som_models, som_clusters) {
  number_of_replicates <- length(som_clusters)
  if (number_of_replicates == 1) return(1L)
  sample_cluster_assignments <- vector("list", number_of_replicates)
  for (replicate_index in seq_len(number_of_replicates)) {
    unit_classif <- as.integer(som_models[[replicate_index]]$unit.classif)
    unit_cluster_labels <- as.integer(som_clusters[[replicate_index]])
    sample_cluster_assignments[[replicate_index]] <- unit_cluster_labels[unit_classif]
  }
  k_values <- vapply(som_clusters, count_SOM_clusters, integer(1))
  k_frequency <- table(k_values[!is.na(k_values)])
  modal_k_values <- as.integer(names(k_frequency)[k_frequency == max(k_frequency)])
  selected_k <- min(modal_k_values)
  candidate_replicates <- which(k_values == selected_k)
  if (length(candidate_replicates) == 1) return(candidate_replicates)
  pairwise_adjusted_rand_index <- matrix(NA_real_, nrow = length(candidate_replicates), ncol = length(candidate_replicates))
  diag(pairwise_adjusted_rand_index) <- NA_real_
  for (candidate_index_1 in seq_len(length(candidate_replicates) - 1)) {
    for (candidate_index_2 in seq.int(candidate_index_1 + 1, length(candidate_replicates))) {
      replicate_index_1 <- candidate_replicates[candidate_index_1]
      replicate_index_2 <- candidate_replicates[candidate_index_2]
      current_adjusted_rand_index <- mclust::adjustedRandIndex(sample_cluster_assignments[[replicate_index_1]], sample_cluster_assignments[[replicate_index_2]])
      pairwise_adjusted_rand_index[candidate_index_1, candidate_index_2] <- current_adjusted_rand_index
      pairwise_adjusted_rand_index[candidate_index_2, candidate_index_1] <- current_adjusted_rand_index
    }
  }
  mean_adjusted_rand_index <- rowMeans(pairwise_adjusted_rand_index, na.rm = TRUE)
  if (all(!is.finite(mean_adjusted_rand_index) | is.na(mean_adjusted_rand_index))) return(candidate_replicates[1])
  representative_candidate_index <- which.max(replace(mean_adjusted_rand_index, !is.finite(mean_adjusted_rand_index) | is.na(mean_adjusted_rand_index), -Inf))
  candidate_replicates[representative_candidate_index]
}
som_models_use <- Aeneus_SOM$som_models
som_clusters_use <- Aeneus_SOM$som_clusters
replicate_index <- choose_representative_replicate(som_models = som_models_use, som_clusters = som_clusters_use)
som_model <- som_models_use[[replicate_index]]
som_cluster <- as.integer(som_clusters_use[[replicate_index]])
neighbor_distances <- calc_unit_neighbor_dist(som_model)
som_cluster <- as.integer(factor(som_cluster, levels = sort(unique(som_cluster))))
number_of_clusters <- length(unique(som_cluster))
cluster_colors <- viridis::viridis(number_of_clusters)
SOM_cluster_plot_col <- cluster_colors[som_cluster]
device_width_in <- (plot_width_cm / 2.54) * (96 / 72)
device_height_in <- (plot_height_cm / 2.54) * (96 / 72)
measurement_file <- tempfile(fileext = ".svg")
svg(measurement_file, width = device_width_in, height = device_height_in, family = "Arial")
base_font_size <- par("ps")
legend_relative_font_size <- (legend_font_size * (96 / 72)) / base_font_size
line_height_in <- par("csi") * par("mex")
panel_gap_in <- (panel_gap_cm / 2.54) * (96 / 72)
neighbor_panel_margin_width_in <- (4.6 + 0.6) * line_height_in
cluster_panel_margin_width_in <- (0.6 + 0.6) * line_height_in
map_plot_width_in <- (device_width_in - panel_gap_in - neighbor_panel_margin_width_in - cluster_panel_margin_width_in) / 2
neighbor_panel_width_in <- map_plot_width_in + neighbor_panel_margin_width_in
neighbor_panel_end_initial <- neighbor_panel_width_in / device_width_in
cluster_panel_start_initial <- (neighbor_panel_width_in + panel_gap_in) / device_width_in
par(family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(fig = c(0, neighbor_panel_end_initial, 0, 1), new = FALSE)
plot(x = som_model, type = "property", property = neighbor_distances, main = "", shape = cluster_shape_neighbor_dist, cex = legend_relative_font_size, palette.name = function(n) rev(viridis::cividis(n)))
left_map_right <- grconvertX(max(som_model$grid$pts[, 1]) + 0.5, from = "user", to = "ndc")
par(fig = c(cluster_panel_start_initial, 1, 0, 1), new = TRUE)
plot(x = som_model, shape = cluster_shape_clusters, type = "mapping", bgcol = SOM_cluster_plot_col, main = "", pchs = point_shape_clusters, cex = point_size_clusters, col = point_color_clusters)
right_map_left <- grconvertX(min(som_model$grid$pts[, 1]) - 0.5, from = "user", to = "ndc")
dev.off()
unlink(measurement_file)
actual_gap_fraction <- right_map_left - left_map_right
desired_gap_fraction <- panel_gap_cm / plot_width_cm
panel_shift_fraction <- (actual_gap_fraction - desired_gap_fraction) / 2
neighbor_panel_start <- panel_shift_fraction
neighbor_panel_end <- neighbor_panel_end_initial + panel_shift_fraction
cluster_panel_start <- cluster_panel_start_initial - panel_shift_fraction
cluster_panel_end <- 1 - panel_shift_fraction
svg(file.path(figure_files_folder, figure_name), width = device_width_in, height = device_height_in, family = "Arial")
base_font_size <- par("ps")
legend_relative_font_size <- (legend_font_size * (96 / 72)) / base_font_size
par(family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(fig = c(neighbor_panel_start, neighbor_panel_end, 0, 1), new = FALSE)
plot(x = som_model, type = "property", property = neighbor_distances, main = "", shape = cluster_shape_neighbor_dist, cex = legend_relative_font_size, palette.name = function(n) rev(viridis::cividis(n)))
par(fig = c(cluster_panel_start, cluster_panel_end, 0, 1), new = TRUE)
plot(x = som_model, shape = cluster_shape_clusters, type = "mapping", bgcol = SOM_cluster_plot_col, main = "", pchs = point_shape_clusters, cex = point_size_clusters, col = point_color_clusters)
if (number_of_clusters > 1) kohonen::add.cluster.boundaries(som_model, som_cluster, lwd = boundary_line_width_clusters, col = boundary_color_clusters)
dev.off()


## Figure 5C
plot_width_cm <- 6.08
plot_height_cm <- 10.09
panel_heights <- c(1, 1, 1.1)
panel_gap <- 0.85
side_tick_label_gap <- 0.65
bottom_tick_label_gap <- 0.5
figure_name <- "Figure_5C.svg"

svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72))
max_k <- Aeneus_SOM$max_k
optim_k_vals <- as.numeric(Aeneus_SOM$optim_k_vals)
BIC_values <- Aeneus_SOM$BIC_values[seq_len(max_k), , drop = FALSE]
support_values <- if (!is.null(Aeneus_SOM$support_values)) Aeneus_SOM$support_values[seq_len(max_k), , drop = FALSE] else BIC_values
k_colors <- viridis::magma(max_k)
base_font_size <- par("ps")
bottom_numbers_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
layout(matrix(1:3, ncol = 1), heights = panel_heights)
par(bty = "n", oma = c(0, 0, 0, 0), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mar = c(panel_gap, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
finite_k_rows <- apply(support_values, 1, function(x) any(is.finite(x) & !is.na(x)))
plotted_k_values <- seq_len(max_k)[finite_k_rows]
values_for_plot <- t(support_values[finite_k_rows, , drop = FALSE])
boxplot(values_for_plot, at = plotted_k_values, xlim = c(0.5, max_k + 0.5), outline = FALSE, notch = FALSE, axes = FALSE, ylab = "", main = "", whisklty = 1, staplelty = 1, col = k_colors[plotted_k_values])
axis(1, at = seq_len(max_k), labels = FALSE)
y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = y_axis_breaks, labels = round(y_axis_breaks, 1), las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
par(mar = c(panel_gap, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
delta_BIC_matrix <- apply(BIC_values, 2, function(x) { previous_BIC <- x[-length(x)]; current_BIC <- x[-1]; delta_BIC <- previous_BIC - current_BIC; delta_BIC[!is.finite(previous_BIC) | !is.finite(current_BIC)] <- NA_real_; delta_BIC })
if (is.null(dim(delta_BIC_matrix))) delta_BIC_matrix <- matrix(delta_BIC_matrix, ncol = 1)
finite_delta_rows <- apply(delta_BIC_matrix, 1, function(x) any(is.finite(x) & !is.na(x)))
plotted_delta_k_values <- seq.int(2, max_k)[finite_delta_rows]
delta_BIC_for_plot <- t(delta_BIC_matrix[finite_delta_rows, , drop = FALSE])
boxplot(delta_BIC_for_plot, at = plotted_delta_k_values, xlim = c(0.5, max_k + 0.5), outline = FALSE, notch = FALSE, axes = FALSE, ylab = "", main = "", whisklty = 1, staplelty = 1, col = k_colors[plotted_delta_k_values])
axis(1, at = seq_len(max_k), labels = FALSE)
y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = y_axis_breaks, labels = round(y_axis_breaks, 1), las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
par(mar = c(2, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
k_frequency_values <- table(factor(optim_k_vals, levels = seq_len(max_k))) / length(optim_k_vals)
bar_midpoints <- barplot(k_frequency_values, ylim = c(0, 1), col = k_colors, axes = FALSE, axisnames = FALSE, ylab = "", main = "")
axis(1, at = bar_midpoints, labels = seq_len(max_k), mgp = c(3, bottom_tick_label_gap, 0), font.axis = 1, cex.axis = bottom_numbers_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
dev.off()


## Figure 5D
plot_width_cm <- 8.72
plot_height_cm <- 7.55
panel_gap <- 1.55
side_tick_label_gap <- 0.65
bottom_tick_label_gap <- 0.75
top_margin_mm <- 2
right_margin_mm <- 2
figure_name <- "Figure_5D.svg"

load(file.path(intermediate_files_folder, "Aeneus_SOM_lolo.Rdata"))
layer_summary_table <- leave.one.layer.out.results$layer.summary
replicate_matched_results_table <- leave.one.layer.out.results$replicate.matched.results
successful_replicate_matched_results_table <- replicate_matched_results_table[is.na(replicate_matched_results_table$error) | replicate_matched_results_table$error == "", , drop = FALSE]
SOM_layer_names <- as.character(layer_summary_table$layer)
SOM_layer_names_plot <- SOM_layer_names
SOM_layer_names_plot[SOM_layer_names_plot == "Alleles"] <- "SNP"
successful_replicate_matched_results_table$layer <- factor(successful_replicate_matched_results_table$layer, levels = SOM_layer_names)
layer_colors <- setNames(viridis::turbo(length(SOM_layer_names)), SOM_layer_names)
show_assignment_margin_plot <- any(is.finite(successful_replicate_matched_results_table$delta.mean.assignment.margin) & !is.na(successful_replicate_matched_results_table$delta.mean.assignment.margin))
svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72))
base_font_size <- par("ps")
axis_labels_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
if (show_assignment_margin_plot) {
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1))
} else {
  layout(matrix(1:2, nrow = 1), widths = c(1, 1))
}
par(bty = "o", oma = c(7, 0, top_margin_lines, right_margin_lines), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mar = c(0, 2, 0, panel_gap), mgp = c(3, side_tick_label_gap, 0))
boxplot(absolute.k.deviation ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names_plot, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
for (layer_index in seq_along(SOM_layer_names)) {
  current_layer_name <- SOM_layer_names[layer_index]
  current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "absolute.k.deviation"]
  current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
  if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
}
box()
par(mar = c(0, 2, 0, if (show_assignment_margin_plot) panel_gap else 0), mgp = c(3, side_tick_label_gap, 0))
boxplot(pairwise.coassignment.change ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names_plot, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
for (layer_index in seq_along(SOM_layer_names)) {
  current_layer_name <- SOM_layer_names[layer_index]
  current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "pairwise.coassignment.change"]
  current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
  if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
}
box()
if (show_assignment_margin_plot) {
  par(mar = c(0, 2, 0, 0), mgp = c(3, side_tick_label_gap, 0))
  boxplot(delta.mean.assignment.margin ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
  axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names_plot, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
  axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
  for (layer_index in seq_along(SOM_layer_names)) {
    current_layer_name <- SOM_layer_names[layer_index]
    current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "delta.mean.assignment.margin"]
    current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
    if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
  }
  box()
}
dev.off()


## Figure 5E
plot_width_cm <- 9.48
plot_height_cm <- 7.38
row_gap <- 1.45
column_gap <- 1.45
bottom_tick_label_gap <- 0.6
top_margin_mm <- 2
left_margin_mm <- 2
right_margin_mm <- 2
figure_name <- "Figure_5E.svg"
bars_threshold_N <- 20
importance_threshold <- 0.0001
matrix_names <- Aeneus_SOM$input_data_names
matrix_names_plot <- matrix_names
matrix_names_plot[matrix_names_plot == "Alleles"] <- "SNP"
first_codebook_list <- kohonen::getCodes(Aeneus_SOM$som_models[[1]])
if (!is.list(first_codebook_list)) first_codebook_list <- list(first_codebook_list)
number_of_layers <- length(first_codebook_list)
if (length(matrix_names) != number_of_layers) matrix_names <- paste0("layer", seq_len(number_of_layers))
if (length(matrix_names_plot) != number_of_layers) matrix_names_plot <- matrix_names
for (layer_index in seq_len(number_of_layers)) {
  if (is.null(colnames(first_codebook_list[[layer_index]]))) colnames(first_codebook_list[[layer_index]]) <- paste0("V", seq_len(ncol(first_codebook_list[[layer_index]])))
}
retained_replicate_indices <- seq_along(Aeneus_SOM$som_models)
Aeneus_SOM_variable_importance <- vector("list", number_of_layers)
names(Aeneus_SOM_variable_importance) <- matrix_names
for (layer_index in seq_len(number_of_layers)) {
  Aeneus_SOM_variable_importance[[layer_index]] <- matrix(NA_real_, nrow = length(retained_replicate_indices), ncol = ncol(first_codebook_list[[layer_index]]), dimnames = list(paste0("R", retained_replicate_indices), colnames(first_codebook_list[[layer_index]])))
}
for (retained_replicate_position in seq_along(retained_replicate_indices)) {
  retained_replicate_index <- retained_replicate_indices[retained_replicate_position]
  som_model <- Aeneus_SOM$som_models[[retained_replicate_index]]
  codebook_list <- kohonen::getCodes(som_model)
  if (!is.list(codebook_list)) codebook_list <- list(codebook_list)
  for (layer_index in seq_len(number_of_layers)) {
    if (is.null(colnames(codebook_list[[layer_index]]))) colnames(codebook_list[[layer_index]]) <- paste0("V", seq_len(ncol(codebook_list[[layer_index]])))
    codebook_matrix <- as.matrix(codebook_list[[layer_index]])
    neuron_sample_counts <- tabulate(som_model$unit.classif, nbins = nrow(codebook_matrix))
    neuron_weights <- neuron_sample_counts + 1
    Aeneus_SOM_variable_importance[[layer_index]][retained_replicate_position, ] <- apply(codebook_matrix, 2, function(variable_values) {
      valid_rows <- is.finite(variable_values) & !is.na(variable_values) & is.finite(neuron_weights) & !is.na(neuron_weights)
      variable_values <- variable_values[valid_rows]
      variable_weights <- neuron_weights[valid_rows]
      if (length(variable_values) < 2 || sum(variable_weights) <= 0) return(NA_real_)
      weighted_variable_mean <- sum(variable_weights * variable_values) / sum(variable_weights)
      sum(variable_weights * (variable_values - weighted_variable_mean)^2) / sum(variable_weights)
    })
  }
}
svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72), family = "Arial")
base_font_size <- par("ps")
layer_label_relative_font_size <- (9 * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
left_margin_lines <- (left_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
par(mfrow = if (number_of_layers <= 3) c(1, number_of_layers) else if (number_of_layers == 4) c(2, 2) else if (number_of_layers <= 6) c(2, 3) else if (number_of_layers <= 8) c(2, 4) else if (number_of_layers == 9) c(3, 3) else c(ceiling(number_of_layers / 3), 3), oma = c(0, left_margin_lines, top_margin_lines, right_margin_lines), mar = c(2.2, column_gap / 2, row_gap, column_gap / 2), mgp = c(3, bottom_tick_label_gap, 0), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1, bty = "o", family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black")
layer_colors <- setNames(viridis::turbo(length(matrix_names)), matrix_names)
for (layer_index in seq_along(Aeneus_SOM_variable_importance)) {
  variable_importance_matrix <- Aeneus_SOM_variable_importance[[layer_index]]
  if (is.null(variable_importance_matrix) || nrow(variable_importance_matrix) == 0 || ncol(variable_importance_matrix) == 0) {
    plot.new()
    next
  }
  median_metric_per_variable <- matrixStats::colMedians(variable_importance_matrix, na.rm = TRUE, useNames = TRUE)
  median_metric_per_variable[!is.finite(median_metric_per_variable)] <- NA_real_
  retained_variable_indices <- which(is.finite(median_metric_per_variable) & !is.na(median_metric_per_variable) & median_metric_per_variable > importance_threshold)
  if (length(retained_variable_indices) == 0) {
    plot.new()
    next
  }
  variable_importance_matrix <- variable_importance_matrix[, retained_variable_indices, drop = FALSE]
  retained_variable_medians <- median_metric_per_variable[colnames(variable_importance_matrix)]
  variable_sort_indices <- order(retained_variable_medians, decreasing = FALSE)
  variable_importance_matrix <- variable_importance_matrix[, variable_sort_indices, drop = FALSE]
  number_of_bars <- ncol(variable_importance_matrix)
  boxplot(variable_importance_matrix, horizontal = TRUE, las = 1, notch = FALSE, outline = FALSE, col = rep(layer_colors[matrix_names[layer_index]], number_of_bars), border = "black", axes = FALSE, whisklty = if (number_of_bars > bars_threshold_N) 0 else 1, staplelty = if (number_of_bars > bars_threshold_N) 0 else 1, names = rep("", number_of_bars))
  axis(1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size, col = "black", col.axis = "black", font.axis = 1)
  mtext(matrix_names_plot[layer_index], side = 3, line = 0.3, cex = layer_label_relative_font_size, font = 1, family = "Arial", col = "black")
}
dev.off()
top5_variables_with_values <- lapply(Aeneus_SOM_variable_importance, function(variable_importance_matrix) {
  median_importance <- matrixStats::colMedians(variable_importance_matrix, na.rm = TRUE, useNames = TRUE)
  median_importance <- median_importance[is.finite(median_importance) & !is.na(median_importance) & median_importance > importance_threshold]
  head(sort(median_importance, decreasing = TRUE), 5)
})
top5_variables_with_values


## Figure 5F
plot_width_cm <- 9.45
plot_height_cm <- 12.5
pie_size <- 2.2
lon_buffer_range <- 4.5
lat_buffer_range <- 5.7
scale_size <- 0.3
bottom_tick_label_gap <- 0.55
side_tick_label_gap <- 0.8
north_arrow_position <- c(0.93, 0.14)
north_arrow_length <- 1.05
north_arrow_N_position <- 0.39
scale_position <- c(0.045, 0.095)
figure_name <- "Figure_5F.svg"
Coordinates <- Aeneus_spatial[, c("Latitude", "Longitude")]

ancestry_matrix <- as.matrix(Aeneus_SOM$ancestry_matrix)
Coordinates <- as.data.frame(Coordinates, stringsAsFactors = FALSE)
Coordinates$Latitude <- as.numeric(Coordinates$Latitude)
Coordinates$Longitude <- as.numeric(Coordinates$Longitude)
coordinate_sample_names <- rownames(Coordinates)
ancestry_sample_names <- rownames(ancestry_matrix)
matched_sample_names <- intersect(ancestry_sample_names, coordinate_sample_names)
Coordinates <- Coordinates[matched_sample_names, , drop = FALSE]
ancestry_matrix <- ancestry_matrix[matched_sample_names, , drop = FALSE]
rows_with_missing_coordinates <- which(!is.finite(Coordinates$Latitude) | !is.finite(Coordinates$Longitude))
if (length(rows_with_missing_coordinates) > 0) {
  Coordinates <- Coordinates[-rows_with_missing_coordinates, , drop = FALSE]
  ancestry_matrix <- ancestry_matrix[-rows_with_missing_coordinates, , drop = FALSE]
}
ancestry_row_sums <- rowSums(ancestry_matrix, na.rm = TRUE)
ancestry_matrix <- ancestry_matrix / ancestry_row_sums
ancestry_proportions <- as.data.frame(ancestry_matrix)
cluster_colors <- viridis::viridis(ncol(ancestry_matrix))
longitude_minimum <- min(Coordinates$Longitude) - lon_buffer_range
longitude_maximum <- max(Coordinates$Longitude) + lon_buffer_range
latitude_minimum <- min(Coordinates$Latitude) - lat_buffer_range
latitude_maximum <- max(Coordinates$Latitude) + lat_buffer_range
pie_radius <- pie_size * 0.01 * max(longitude_maximum - longitude_minimum, latitude_maximum - latitude_minimum)
add_admixture_pie <- function(longitude, latitude, ancestry_proportions, cluster_colors, x.radius, y.radius, border.color = "black", line.width = 0.8, number.of.points = 80) {
  ancestry_proportions <- as.numeric(ancestry_proportions)
  ancestry_proportions[is.na(ancestry_proportions) | !is.finite(ancestry_proportions) | ancestry_proportions < 0] <- 0
  if (sum(ancestry_proportions) <= 0) return(invisible(NULL))
  ancestry_proportions <- ancestry_proportions / sum(ancestry_proportions)
  slice_start_angles <- c(0, cumsum(ancestry_proportions)[-length(ancestry_proportions)]) * 2 * pi
  slice_end_angles <- cumsum(ancestry_proportions) * 2 * pi
  for (slice_index in seq_along(ancestry_proportions)) {
    if (ancestry_proportions[slice_index] <= 0) next
    slice_angles <- seq(slice_start_angles[slice_index], slice_end_angles[slice_index], length.out = max(3, ceiling(number.of.points * ancestry_proportions[slice_index])))
    polygon(x = c(longitude, longitude + x.radius * cos(slice_angles), longitude), y = c(latitude, latitude + y.radius * sin(slice_angles), latitude), col = cluster_colors[slice_index], border = border.color, lwd = line.width)
  }
  circle_angles <- seq(0, 2 * pi, length.out = number.of.points)
  lines(longitude + x.radius * cos(circle_angles), latitude + y.radius * sin(circle_angles), col = border.color, lwd = line.width)
  invisible(NULL)
}
add_map_legend <- function(legend.position, legend.title, legend.labels, legend.colors, legend.text.relative.font.size, legend.title.relative.font.size, legend.text.font, legend.symbol.size, legend.box) {
  plot_coordinate_limits <- par("usr")
  plot_longitude_range <- plot_coordinate_limits[2] - plot_coordinate_limits[1]
  plot_latitude_range <- plot_coordinate_limits[4] - plot_coordinate_limits[3]
  legend_text_width <- max(strwidth(legend.labels, units = "user", cex = legend.text.relative.font.size, font = legend.text.font, family = "Arial"))
  legend_text_height <- strheight("M", units = "user", cex = legend.text.relative.font.size, font = legend.text.font, family = "Arial")
  legend_title_width <- 0
  legend_title_height <- 0
  legend_title_gap <- 0
  if (!is.null(legend.title)) {
    legend_title_width <- strwidth(legend.title, units = "user", cex = legend.title.relative.font.size, font = 2, family = "Arial")
    legend_title_height <- strheight("M", units = "user", cex = legend.title.relative.font.size, font = 2, family = "Arial")
    legend_title_gap <- 0.5 * legend_text_height
  }
  legend_symbol_width <- strwidth("M", units = "user", cex = legend.text.relative.font.size, family = "Arial") * legend.symbol.size
  legend_symbol_gap <- 0.7 * strwidth("M", units = "user", cex = legend.text.relative.font.size, family = "Arial")
  legend_padding_x <- 0.9 * strwidth("M", units = "user", cex = legend.text.relative.font.size, family = "Arial")
  legend_padding_y <- 0.7 * legend_text_height
  legend_line_gap <- 0.35 * legend_text_height
  legend_width <- max(legend_title_width, legend_symbol_width + legend_symbol_gap + legend_text_width) + 2 * legend_padding_x
  legend_height <- 2 * legend_padding_y + legend_title_height + legend_title_gap + length(legend.labels) * legend_text_height + (length(legend.labels) - 1) * legend_line_gap
  if (legend.position == "topright") {
    legend_left <- plot_coordinate_limits[2] - legend_width
    legend_bottom <- plot_coordinate_limits[4] - legend_height
  } else if (legend.position == "topleft") {
    legend_left <- plot_coordinate_limits[1]
    legend_bottom <- plot_coordinate_limits[4] - legend_height
  } else if (legend.position == "bottomright") {
    legend_left <- plot_coordinate_limits[2] - legend_width
    legend_bottom <- plot_coordinate_limits[3]
  } else if (legend.position == "bottomleft") {
    legend_left <- plot_coordinate_limits[1]
    legend_bottom <- plot_coordinate_limits[3]
  } else if (legend.position == "right") {
    legend_left <- plot_coordinate_limits[2] - legend_width
    legend_bottom <- plot_coordinate_limits[3] + 0.5 * plot_latitude_range - 0.5 * legend_height
  } else if (legend.position == "left") {
    legend_left <- plot_coordinate_limits[1]
    legend_bottom <- plot_coordinate_limits[3] + 0.5 * plot_latitude_range - 0.5 * legend_height
  } else if (legend.position == "top") {
    legend_left <- plot_coordinate_limits[1] + 0.5 * plot_longitude_range - 0.5 * legend_width
    legend_bottom <- plot_coordinate_limits[4] - legend_height
  } else if (legend.position == "bottom") {
    legend_left <- plot_coordinate_limits[1] + 0.5 * plot_longitude_range - 0.5 * legend_width
    legend_bottom <- plot_coordinate_limits[3]
  } else {
    legend_left <- plot_coordinate_limits[1] + 0.5 * plot_longitude_range - 0.5 * legend_width
    legend_bottom <- plot_coordinate_limits[3] + 0.5 * plot_latitude_range - 0.5 * legend_height
  }
  legend_right <- legend_left + legend_width
  legend_top <- legend_bottom + legend_height
  if (legend.box) rect(legend_left, legend_bottom, legend_right, legend_top, col = "white", border = "black")
  current_legend_y_position <- legend_top - legend_padding_y
  if (!is.null(legend.title)) {
    text(x = legend_left + legend_padding_x, y = current_legend_y_position - 0.5 * legend_title_height, labels = legend.title, adj = c(0, 0.5), font = 2, cex = legend.title.relative.font.size, family = "Arial", col = "black")
    current_legend_y_position <- current_legend_y_position - legend_title_height - legend_title_gap
  }
  legend_symbol_x_position <- legend_left + legend_padding_x + 0.5 * legend_symbol_width
  legend_text_x_position <- legend_left + legend_padding_x + legend_symbol_width + legend_symbol_gap
  for (legend_index in seq_along(legend.labels)) {
    legend_entry_y_position <- current_legend_y_position - 0.5 * legend_text_height - (legend_index - 1) * (legend_text_height + legend_line_gap)
    points(x = legend_symbol_x_position, y = legend_entry_y_position, pch = 21, cex = legend.symbol.size, bg = legend.colors[legend_index], col = "black")
    text(x = legend_text_x_position, y = legend_entry_y_position, labels = legend.labels[legend_index], adj = c(0, 0.5), cex = legend.text.relative.font.size, font = legend.text.font, family = "Arial", col = "black")
  }
  invisible(NULL)
}
svg(file.path(figure_files_folder, figure_name), width = plot_width_cm / 2.54, height = plot_height_cm / 2.54, family = "Arial")
base_font_size <- par("ps")
scale_text_relative_font_size <- 7.1 / base_font_size
axis_labels_relative_font_size <- 9.1 / base_font_size
legend_text_relative_font_size <- 9.1 / base_font_size
legend_title_relative_font_size <- 9.1 / base_font_size
par(mfrow = c(1, 1), oma = c(2, 2, 1, 1), mar = c(1, 1, 1, 1), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1, family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", bty = "o")
current_plot_region_size_inches <- par("pin")
longitude_range <- longitude_maximum - longitude_minimum
latitude_range <- latitude_maximum - latitude_minimum
mean_map_latitude <- mean(c(latitude_minimum, latitude_maximum))
longitude_latitude_correction <- cos(mean_map_latitude * pi / 180)
if (!is.finite(longitude_latitude_correction) || longitude_latitude_correction <= 0) longitude_latitude_correction <- 1
target_height_width_ratio <- latitude_range / (longitude_range * longitude_latitude_correction)
adjusted_plot_width_inches <- current_plot_region_size_inches[1]
adjusted_plot_height_inches <- adjusted_plot_width_inches * target_height_width_ratio
if (adjusted_plot_height_inches > current_plot_region_size_inches[2]) {
  adjusted_plot_height_inches <- current_plot_region_size_inches[2]
  adjusted_plot_width_inches <- adjusted_plot_height_inches / target_height_width_ratio
}
par(pin = c(adjusted_plot_width_inches, adjusted_plot_height_inches))
plot.new()
plot.window(xlim = c(longitude_minimum, longitude_maximum), ylim = c(latitude_minimum, latitude_maximum), xaxs = "i", yaxs = "i")
plot_coordinate_limits <- par("usr")
plot_region_size_inches <- par("pin")
x_units_per_inch <- plot_region_size_inches[1] / diff(plot_coordinate_limits[1:2])
y_units_per_inch <- plot_region_size_inches[2] / diff(plot_coordinate_limits[3:4])
if (!is.finite(x_units_per_inch) || !is.finite(y_units_per_inch) || x_units_per_inch <= 0 || y_units_per_inch <= 0) {
  pie_radius_x <- pie_radius
  pie_radius_y <- pie_radius
} else {
  pie_radius_x <- pie_radius * y_units_per_inch / x_units_per_inch
  pie_radius_y <- pie_radius
}
maps::map("world", fill = TRUE, col = "lightgrey", border = "black", xlim = c(longitude_minimum, longitude_maximum), ylim = c(latitude_minimum, latitude_maximum), add = TRUE)
try(maps::map("state", add = TRUE, xlim = c(longitude_minimum, longitude_maximum), ylim = c(latitude_minimum, latitude_maximum), col = "black", lwd = 0.5), silent = TRUE)
longitude_clip_buffer <- 0.2 * (longitude_maximum - longitude_minimum)
latitude_clip_buffer <- 0.2 * (latitude_maximum - latitude_minimum)
clip(longitude_minimum - longitude_clip_buffer, longitude_maximum + longitude_clip_buffer, latitude_minimum - latitude_clip_buffer, latitude_maximum + latitude_clip_buffer)
axis(1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size, col = "black", col.axis = "black", font.axis = 1)
axis(2, las = 2, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size, col = "black", col.axis = "black", font.axis = 1)
box(col = "black")
for (sample_index in seq_len(nrow(ancestry_proportions))) {
  add_admixture_pie(longitude = Coordinates$Longitude[sample_index], latitude = Coordinates$Latitude[sample_index], ancestry_proportions = as.numeric(ancestry_proportions[sample_index, ]), cluster_colors = cluster_colors, x.radius = pie_radius_x, y.radius = pie_radius_y)
}
legend_labels <- paste("Cluster", seq_along(cluster_colors))
add_map_legend(legend.position = "topright", legend.title = "Cluster", legend.labels = legend_labels, legend.colors = cluster_colors, legend.text.relative.font.size = legend_text_relative_font_size, legend.title.relative.font.size = legend_title_relative_font_size, legend.text.font = 1, legend.symbol.size = 1.6, legend.box = TRUE)
scale_position_longitude <- scale_position[1] * (longitude_maximum - longitude_minimum) + longitude_minimum
scale_position_latitude <- scale_position[2] * (latitude_maximum - latitude_minimum) + latitude_minimum
maps::map.scale(x = scale_position_longitude, y = scale_position_latitude, cex = scale_text_relative_font_size, relwidth = scale_size, ratio = FALSE)
north_arrow_longitude <- north_arrow_position[1] * (longitude_maximum - longitude_minimum) + longitude_minimum
north_arrow_latitude <- north_arrow_position[2] * (latitude_maximum - latitude_minimum) + latitude_minimum
arrows(x0 = north_arrow_longitude, y0 = north_arrow_latitude, x1 = north_arrow_longitude, y1 = north_arrow_latitude + north_arrow_length, length = 0.13, col = "black", lwd = 1.7)
text(x = north_arrow_longitude, y = north_arrow_latitude + north_arrow_length + north_arrow_N_position, labels = "N", cex = 0.8, col = "black", family = "Arial", font = 1)
dev.off()




#### Supplementary Figure S1 ###################################################

## Set environment
rm(list = ls()) #clear environment
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")


## Install and load required packages
required_packages <- c("emmeans",
                       "glmmTMB",
                       "nlme",
                       "patchwork",
                       "svglite",
                       "tidyverse")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
  library(package_name, character.only = TRUE)
}


## Set main plotting parameters
plot_width_cm <- 16
plot_height_cm <- 20

plot_font_family <- "Arial"
plot_text_color <- "black"

plot_base_size <- 9
plot_axis_title_size <- 9
plot_axis_text_size <- 7.1
plot_column_title_size <- 9

plot_raw_point_alpha <- 0.35
plot_raw_point_size <- 1
plot_estimate_point_size <- 2.6
plot_estimate_linewidth <- 0.7
plot_errorbar_linewidth <- 1.3
plot_errorbar_width <- 0
plot_categorical_jitter_width <- 0.06
plot_SE_multiplier <- 1
X_text_angle_plot <- 90
minimum_column_width_units <- 3.5
panel_margins_mm <- 2.5

plot_color_clustering_methods <- "#56B4E9"
plot_color_neighborhood <- "#009E73"
plot_color_N_steps <- "#E69F00"
plot_color_NA <- "#D55E00"
plot_color_learning_rate_tuning <- "#CC80A7"


## Set input and output
simulation_k3_dir <- file.path("Simulations", "Simulation_set_1", "k3")
sim_results_clustering_methods_csv <- "Sim_results_clustering_methods.csv"
sim_results_neighborhoods_csv <- "Sim_results_neighborhoods.csv"
sim_results_N_steps_csv <- "Sim_results_N_steps.csv"
sim_results_NA_csv <- "Sim_results_NA.csv"
sim_results_learning_rate_tuning_csv <- "Sim_results_learning_rate_tuning.csv"

plot_dir <- "Figure_files"
plot_file_name <- "Supplementary_Figure_S1.svg"


## Combine directories and file names
sim_results_clustering_methods_csv <- file.path(simulation_k3_dir, sim_results_clustering_methods_csv)
sim_results_neighborhoods_csv <- file.path(simulation_k3_dir, sim_results_neighborhoods_csv)
sim_results_N_steps_csv <- file.path(simulation_k3_dir, sim_results_N_steps_csv)
sim_results_NA_csv <- file.path(simulation_k3_dir, sim_results_NA_csv)
sim_results_learning_rate_tuning_csv <- file.path(simulation_k3_dir, sim_results_learning_rate_tuning_csv)
combined_svg <- file.path(plot_dir, plot_file_name)


## Create output directory if needed
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


## Check that all completed result files exist
required_result_files <- c(sim_results_clustering_methods_csv,
                           sim_results_neighborhoods_csv,
                           sim_results_N_steps_csv,
                           sim_results_NA_csv,
                           sim_results_learning_rate_tuning_csv)
missing_result_files <- required_result_files[!file.exists(required_result_files)]
if (length(missing_result_files) > 0) stop("The following completed result files are missing: ", paste(missing_result_files, collapse = ", "))


## Define tests
test_specs <- list(
  clustering_methods = list(
    title = "Clustering\nmethod",
    condition_variable = "clustering_method",
    condition_levels = c("kmeans+BICelbow",
                         "kmeans+BICthreshold",
                         "GMM+BICthreshold",
                         "hierarchical+DB",
                         "HDBSCAN",
                         "OPTICS+Silhouette"),
    estimate_color = plot_color_clustering_methods,
    connect_estimates = FALSE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  neighborhoods = list(
    title = "Neighborhood\nfunction",
    condition_variable = "neighborhood_function",
    condition_levels = c("gaussian", "bubble"),
    estimate_color = plot_color_neighborhood,
    connect_estimates = FALSE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  N_steps = list(
    title = "Number of training\niterations",
    condition_variable = "N_steps",
    condition_levels = c(20,
                         50,
                         100,
                         200,
                         500,
                         1000,
                         2000),
    estimate_color = plot_color_N_steps,
    connect_estimates = TRUE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  learning_rate_tuning = list(
    title = "Learning-rate\ntuning",
    condition_variable = "learning_rate_tuning",
    condition_levels = c(FALSE, TRUE),
    estimate_color = plot_color_learning_rate_tuning,
    connect_estimates = FALSE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  missing_data = list(
    title = "Missing data\nproportion",
    condition_variable = "missing_data_prop",
    condition_levels = c(0.0,
                         0.1,
                         0.2,
                         0.3,
                         0.4,
                         0.5,
                         0.6,
                         0.7,
                         0.8),
    estimate_color = plot_color_NA,
    connect_estimates = TRUE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  )
)


## Define metrics
metric_specs <- list(
  ARI = list(
    response_variable = "ARI",
    metric_label = "Adjusted Rand Index"
  ),
  QE = list(
    response_variable = "QE",
    metric_label = "Quantization error"
  ),
  TE = list(
    response_variable = "TE",
    metric_label = "Topographic error"
  ),
  Time = list(
    response_variable = "Time",
    metric_label = "Time (seconds)"
  )
)


## Calculate category-based column widths
test_category_counts <- vapply(test_specs, function(test_spec) length(test_spec$condition_levels), numeric(1))
column_width_units <- pmax(test_category_counts, minimum_column_width_units)


## Create additional functions
safe.mean <- function(input_vector) {
  if (all(is.na(input_vector))) return(NA_real_)
  mean(input_vector, na.rm = TRUE)
}
prepare.continuous.benchmark.data <- function(input_data,
                                              condition_variable,
                                              response_variable) {
  model_data <- input_data[ , c("sim", condition_variable, response_variable), drop = FALSE]
  colnames(model_data)[colnames(model_data) == condition_variable] <- "condition"
  colnames(model_data)[colnames(model_data) == response_variable] <- "response"
  model_data$response <- suppressWarnings(as.numeric(model_data$response))
  model_data <- model_data[!is.na(model_data$condition) & !is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data <- model_data %>%
    dplyr::group_by(sim, condition) %>%
    dplyr::summarise(response = safe.mean(response), .groups = "drop")
  model_data <- model_data[!is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data
}


## Create plotting functions
create.benchmark.plot.theme <- function() {
  theme_classic(base_size = plot_base_size, base_family = plot_font_family) +
    theme(text = element_text(family = plot_font_family, colour = plot_text_color),
          axis.title = element_text(family = plot_font_family, size = plot_axis_title_size),
          axis.text = element_text(family = plot_font_family, size = plot_axis_text_size),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          axis.ticks.length = grid::unit(0.8, "mm"),
          axis.text.x.top = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.line.x.top = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          plot.title = element_text(family = plot_font_family,
                                    size = plot_column_title_size,
                                    face = "bold",
                                    hjust = 0.5,
                                    lineheight = 0.9,
                                    margin = margin(b = 1.2, unit = "mm")),
          plot.margin = margin(1, 1, 1, 1, unit = "mm"))
}
plot.categorical.gaussian.benchmark.model <- function(input_data,
                                                      condition_variable,
                                                      response_variable,
                                                      condition_levels,
                                                      estimate_color,
                                                      connect_estimates = FALSE) {
  model_data <- prepare.continuous.benchmark.data(input_data = input_data, condition_variable = condition_variable, response_variable = response_variable)
  model_data$condition_factor <- factor(as.character(model_data$condition), levels = as.character(condition_levels))
  model_data <- model_data[!is.na(model_data$condition_factor), , drop = FALSE]
  if (length(unique(model_data$sim)) < 2) stop("Gaussian mixed model requires at least two simulation replicates")
  model_output <- nlme::lme(response ~ condition_factor, random = ~ 1 | sim, data = model_data, na.action = na.omit, method = "REML")
  emmeans_output <- emmeans::emmeans(model_output, ~ condition_factor)
  emmeans_output_df <- as.data.frame(emmeans_output)
  emmeans_output_df$lower_SE <- emmeans_output_df$emmean - plot_SE_multiplier * emmeans_output_df$SE
  emmeans_output_df$upper_SE <- emmeans_output_df$emmean + plot_SE_multiplier * emmeans_output_df$SE
  model_plot <- ggplot() +
    geom_point(data = model_data,
               aes(x = condition_factor, y = response),
               alpha = plot_raw_point_alpha,
               size = plot_raw_point_size,
               position = position_jitter(width = plot_categorical_jitter_width, height = 0)) +
    geom_errorbar(data = emmeans_output_df,
                  aes(x = condition_factor, ymin = lower_SE, ymax = upper_SE),
                  colour = estimate_color,
                  linewidth = plot_errorbar_linewidth,
                  width = plot_errorbar_width)
  if (connect_estimates) {
    model_plot <- model_plot +
      geom_line(data = emmeans_output_df,
                aes(x = condition_factor, y = emmean, group = 1),
                colour = estimate_color,
                linewidth = plot_estimate_linewidth)
  }
  model_plot <- model_plot +
    geom_point(data = emmeans_output_df,
               aes(x = condition_factor, y = emmean),
               colour = estimate_color,
               size = plot_estimate_point_size) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
    create.benchmark.plot.theme()
  list(model = model_output,
       emmeans = emmeans_output,
       emmeans_df = emmeans_output_df,
       plot = model_plot)
}
format.benchmark.panel <- function(input_plot,
                                   test_title = NULL,
                                   metric_label = NULL,
                                   show_x_axis = FALSE,
                                   x_text_angle = 0,
                                   x_text_hjust = 1,
                                   x_text_vjust = 0.5,
                                   panel_tag = NULL) {
  output_plot <- input_plot +
    labs(title = test_title,
         x = NULL,
         y = metric_label,
         tag = panel_tag) +
    theme(axis.title.x = element_blank(),
          axis.text.y = element_text(family = plot_font_family, size = plot_axis_text_size),
          axis.ticks.y = element_line(linewidth = 0.3),
          axis.line.y = element_line(linewidth = 0.3))
  if (is.null(test_title)) output_plot <- output_plot + theme(plot.title = element_blank())
  if (is.null(metric_label)) {
    output_plot <- output_plot + theme(axis.title.y = element_blank())
  } else {
    output_plot <- output_plot + theme(axis.title.y = element_text(family = plot_font_family, size = plot_axis_title_size, face = "bold", margin = margin(r = 2, unit = "mm")))
  }
  if (show_x_axis) {
    output_plot <- output_plot +
      theme(axis.text.x = element_text(family = plot_font_family, size = plot_axis_text_size, angle = x_text_angle, hjust = x_text_hjust, vjust = x_text_vjust),
            axis.ticks.x = element_line(linewidth = 0.3),
            axis.line.x = element_line(linewidth = 0.3))
  } else {
    output_plot <- output_plot +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_line(linewidth = 0.3))
  }
  output_plot
}


## Load and check completed analyses
analysis_data <- list(
  clustering_methods = read.csv(sim_results_clustering_methods_csv, stringsAsFactors = FALSE),
  neighborhoods = read.csv(sim_results_neighborhoods_csv, stringsAsFactors = FALSE),
  N_steps = read.csv(sim_results_N_steps_csv, stringsAsFactors = FALSE),
  missing_data = read.csv(sim_results_NA_csv, stringsAsFactors = FALSE),
  learning_rate_tuning = read.csv(sim_results_learning_rate_tuning_csv, stringsAsFactors = FALSE)
)
required_metric_columns <- vapply(metric_specs, function(metric_spec) metric_spec$response_variable, character(1))
for (test_name in names(test_specs)) {
  condition_variable <- test_specs[[test_name]]$condition_variable
  required_columns <- c("sim", condition_variable, required_metric_columns)
  missing_columns <- setdiff(required_columns, colnames(analysis_data[[test_name]]))
  if (length(missing_columns) > 0) stop("The result table for ", test_specs[[test_name]]$title, " is missing the following columns: ", paste(missing_columns, collapse = ", "))
}


## Create plots for all remaining metrics
metric_plots <- lapply(metric_specs, function(metric_spec) vector(mode = "list", length = length(test_specs)))
metric_models <- lapply(metric_specs, function(metric_spec) list())
metric_emmeans <- lapply(metric_specs, function(metric_spec) list())
for (metric_index in seq_along(metric_specs)) {
  metric_name <- names(metric_specs)[metric_index]
  metric_spec <- metric_specs[[metric_name]]
  for (test_index in seq_along(test_specs)) {
    test_name <- names(test_specs)[test_index]
    test_spec <- test_specs[[test_name]]
    model_plot_output <- plot.categorical.gaussian.benchmark.model(
      input_data = analysis_data[[test_name]],
      condition_variable = test_spec$condition_variable,
      response_variable = metric_spec$response_variable,
      condition_levels = test_spec$condition_levels,
      estimate_color = test_spec$estimate_color,
      connect_estimates = test_spec$connect_estimates
    )
    test_title <- if (metric_index == 1) test_spec$title else NULL
    metric_label <- if (test_index == 1) metric_spec$metric_label else NULL
    show_x_axis <- metric_index == length(metric_specs)
    metric_plots[[metric_name]][[test_index]] <- format.benchmark.panel(
      input_plot = model_plot_output$plot,
      test_title = test_title,
      metric_label = metric_label,
      show_x_axis = show_x_axis,
      x_text_angle = test_spec$x_text_angle,
      x_text_hjust = test_spec$x_text_hjust,
      x_text_vjust = test_spec$x_text_vjust,
      panel_tag = NULL
    )
    metric_models[[metric_name]][[test_name]] <- model_plot_output$model
    metric_emmeans[[metric_name]][[test_name]] <- model_plot_output$emmeans
  }
}


## Combine plots using category-based column widths
combined_plot_list <- do.call(c, metric_plots)
combined_plot <- patchwork::wrap_plots(plotlist = combined_plot_list,
                                       ncol = length(test_specs),
                                       nrow = length(metric_specs),
                                       byrow = TRUE) +
  patchwork::plot_layout(widths = unname(column_width_units), heights = rep(1, length(metric_specs))) &
  theme(text = element_text(family = plot_font_family),
        plot.margin = margin(panel_margins_mm, panel_margins_mm, panel_margins_mm, panel_margins_mm, unit = "mm"))


## Show and save combined plot
combined_plot
ggsave(filename = combined_svg,
       plot = combined_plot,
       device = svglite::svglite,
       width = plot_width_cm,
       height = plot_height_cm,
       units = "cm",
       bg = "white",
       limitsize = FALSE)
message("Combined SVG saved to: ", combined_svg)




#### Supplementary Figure S2 ###################################################

## Set environment
rm(list = ls()) #clear environment
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")


## Install and load required packages
required_packages <- c("emmeans",
                       "glmmTMB",
                       "nlme",
                       "patchwork",
                       "svglite",
                       "tidyverse")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
  library(package_name, character.only = TRUE)
}


## Set main plotting parameters
plot_width_cm <- 16
plot_height_cm <- 15

plot_font_family <- "Arial"
plot_text_color <- "black"

plot_base_size <- 9
plot_axis_title_size <- 9
plot_axis_text_size <- 7.1
plot_column_title_size <- 9

plot_raw_point_alpha <- 0.35
plot_raw_point_size <- 1
plot_estimate_point_size <- 2.6
plot_estimate_linewidth <- 0.7
plot_errorbar_linewidth <- 1.3
plot_errorbar_width <- 0
plot_categorical_jitter_width <- 0.06
plot_SE_multiplier <- 1
X_text_angle_plot <- 90
minimum_column_width_units <- 3.5
panel_margins_mm <- 2.5

plot_color_clustering_methods <- "#56B4E9"
plot_color_neighborhood <- "#009E73"
plot_color_N_steps <- "#E69F00"
plot_color_NA <- "#D55E00"
plot_color_learning_rate_tuning <- "#CC80A7"


## Set input and output
simulation_k7_dir <- file.path("Simulations", "Simulation_set_1", "k7")
sim_results_clustering_methods_csv <- "Sim_results_clustering_methods.csv"
sim_results_neighborhoods_csv <- "Sim_results_neighborhoods.csv"
sim_results_N_steps_csv <- "Sim_results_N_steps.csv"
sim_results_NA_csv <- "Sim_results_NA.csv"
sim_results_learning_rate_tuning_csv <- "Sim_results_learning_rate_tuning.csv"

plot_dir <- "Figure_files"
plot_file_name <- "Supplementary_Figure_S2.svg"


## Combine directories and file names
sim_results_clustering_methods_csv <- file.path(simulation_k7_dir, sim_results_clustering_methods_csv)
sim_results_neighborhoods_csv <- file.path(simulation_k7_dir, sim_results_neighborhoods_csv)
sim_results_N_steps_csv <- file.path(simulation_k7_dir, sim_results_N_steps_csv)
sim_results_NA_csv <- file.path(simulation_k7_dir, sim_results_NA_csv)
sim_results_learning_rate_tuning_csv <- file.path(simulation_k7_dir, sim_results_learning_rate_tuning_csv)
combined_svg <- file.path(plot_dir, plot_file_name)


## Create output directory if needed
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


## Check that all completed result files exist
required_result_files <- c(sim_results_clustering_methods_csv,
                           sim_results_neighborhoods_csv,
                           sim_results_N_steps_csv,
                           sim_results_NA_csv,
                           sim_results_learning_rate_tuning_csv)
missing_result_files <- required_result_files[!file.exists(required_result_files)]
if (length(missing_result_files) > 0) stop("The following completed result files are missing: ", paste(missing_result_files, collapse = ", "))


## Define tests
test_specs <- list(
  clustering_methods = list(
    title = "Clustering\nmethod",
    condition_variable = "clustering_method",
    condition_levels = c("kmeans+BICelbow",
                         "kmeans+BICthreshold",
                         "GMM+BICthreshold",
                         "hierarchical+DB",
                         "HDBSCAN",
                         "OPTICS+Silhouette"),
    estimate_color = plot_color_clustering_methods,
    connect_estimates = FALSE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  neighborhoods = list(
    title = "Neighborhood\nfunction",
    condition_variable = "neighborhood_function",
    condition_levels = c("gaussian", "bubble"),
    estimate_color = plot_color_neighborhood,
    connect_estimates = FALSE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  N_steps = list(
    title = "N training\niterations",
    condition_variable = "N_steps",
    condition_levels = c(20,
                         50,
                         100,
                         200,
                         500,
                         1000,
                         2000),
    estimate_color = plot_color_N_steps,
    connect_estimates = TRUE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  learning_rate_tuning = list(
    title = "Learning-rate\ntuning",
    condition_variable = "learning_rate_tuning",
    condition_levels = c(FALSE, TRUE),
    estimate_color = plot_color_learning_rate_tuning,
    connect_estimates = FALSE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  missing_data = list(
    title = "Missing data\nproportion",
    condition_variable = "missing_data_prop",
    condition_levels = c(0.0,
                         0.1,
                         0.2,
                         0.3,
                         0.4,
                         0.5,
                         0.6,
                         0.7,
                         0.8),
    estimate_color = plot_color_NA,
    connect_estimates = TRUE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  )
)


## Calculate category-based column widths
test_category_counts <- vapply(test_specs, function(test_spec) length(test_spec$condition_levels), numeric(1))
column_width_units <- pmax(test_category_counts, minimum_column_width_units)


## Create additional functions
safe.mean <- function(input_vector) {
  if (all(is.na(input_vector))) return(NA_real_)
  mean(input_vector, na.rm = TRUE)
}
convert.to.logical <- function(input_vector) {
  if (is.logical(input_vector)) return(input_vector)
  if (is.numeric(input_vector)) return(ifelse(is.na(input_vector), NA, input_vector != 0))
  input_vector <- tolower(trimws(as.character(input_vector)))
  output_vector <- rep(NA, length(input_vector))
  output_vector[input_vector %in% c("true", "t", "1", "yes", "y")] <- TRUE
  output_vector[input_vector %in% c("false", "f", "0", "no", "n")] <- FALSE
  output_vector
}
prepare.continuous.benchmark.data <- function(input_data,
                                              condition_variable,
                                              response_variable) {
  model_data <- input_data[ , c("sim", condition_variable, response_variable), drop = FALSE]
  colnames(model_data)[colnames(model_data) == condition_variable] <- "condition"
  colnames(model_data)[colnames(model_data) == response_variable] <- "response"
  model_data$response <- suppressWarnings(as.numeric(model_data$response))
  model_data <- model_data[!is.na(model_data$condition) & !is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data <- model_data %>%
    dplyr::group_by(sim, condition) %>%
    dplyr::summarise(response = safe.mean(response), .groups = "drop")
  model_data <- model_data[!is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data
}
prepare.binary.benchmark.data <- function(input_data,
                                          condition_variable,
                                          response_variable) {
  model_data <- input_data[ , c("sim", condition_variable, response_variable), drop = FALSE]
  colnames(model_data)[colnames(model_data) == condition_variable] <- "condition"
  colnames(model_data)[colnames(model_data) == response_variable] <- "response"
  model_data$response <- as.numeric(convert.to.logical(model_data$response))
  model_data <- model_data[!is.na(model_data$condition) & !is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data
}
extract.response.emmeans.columns <- function(emmeans_output_df) {
  response_column <- if ("prob" %in% colnames(emmeans_output_df)) {
    "prob"
  } else if ("response" %in% colnames(emmeans_output_df)) {
    "response"
  } else {
    stop("Could not find the response column in the emmeans output")
  }
  emmeans_output_df$estimated_response <- emmeans_output_df[[response_column]]
  emmeans_output_df$lower_SE <- pmax(0, emmeans_output_df$estimated_response - plot_SE_multiplier * emmeans_output_df$SE)
  emmeans_output_df$upper_SE <- pmin(1, emmeans_output_df$estimated_response + plot_SE_multiplier * emmeans_output_df$SE)
  emmeans_output_df
}


## Create plotting functions
create.benchmark.plot.theme <- function() {
  theme_classic(base_size = plot_base_size, base_family = plot_font_family) +
    theme(text = element_text(family = plot_font_family, colour = plot_text_color),
          axis.title = element_text(family = plot_font_family, size = plot_axis_title_size),
          axis.text = element_text(family = plot_font_family, size = plot_axis_text_size),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          axis.ticks.length = grid::unit(0.8, "mm"),
          axis.text.x.top = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.line.x.top = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          plot.title = element_text(family = plot_font_family,
                                    size = plot_column_title_size,
                                    face = "bold",
                                    hjust = 0.5,
                                    lineheight = 0.9,
                                    margin = margin(b = 1.2, unit = "mm")),
          plot.margin = margin(1, 1, 1, 1, unit = "mm"))
}
plot.categorical.binary.benchmark.model <- function(input_data,
                                                    condition_variable,
                                                    response_variable,
                                                    condition_levels,
                                                    estimate_color,
                                                    connect_estimates = FALSE) {
  model_data <- prepare.binary.benchmark.data(input_data = input_data, condition_variable = condition_variable, response_variable = response_variable)
  model_data$condition_factor <- factor(as.character(model_data$condition), levels = as.character(condition_levels))
  model_data <- model_data[!is.na(model_data$condition_factor), , drop = FALSE]
  model_data$sim_condition <- interaction(model_data$sim, model_data$condition_factor, drop = TRUE)
  model_data$sim_condition <- factor(model_data$sim_condition)
  if (length(unique(model_data$sim)) < 2) stop("Binary mixed model requires at least two simulation replicates")
  model_output <- glmmTMB::glmmTMB(response ~ condition_factor + (1 | sim) + (1 | sim_condition),family = binomial, data = model_data)
  emmeans_output <- emmeans::emmeans(model_output, ~ condition_factor, type = "response")
  emmeans_output_df <- as.data.frame(emmeans_output)
  emmeans_output_df <- extract.response.emmeans.columns(emmeans_output_df)
  raw_plot_data <- model_data %>%
    dplyr::group_by(sim, condition_factor) %>%
    dplyr::summarise(response = safe.mean(response), .groups = "drop")
  model_plot <- ggplot() +
    geom_point(data = raw_plot_data,
               aes(x = condition_factor, y = response),
               alpha = plot_raw_point_alpha,
               size = plot_raw_point_size,
               position = position_jitter(width = plot_categorical_jitter_width, height = 0)) +
    geom_errorbar(data = emmeans_output_df,
                  aes(x = condition_factor, ymin = lower_SE, ymax = upper_SE),
                  colour = estimate_color,
                  linewidth = plot_errorbar_linewidth,
                  width = plot_errorbar_width)
  if (connect_estimates) {
    model_plot <- model_plot +
      geom_line(data = emmeans_output_df,
                aes(x = condition_factor, y = estimated_response, group = 1),
                colour = estimate_color,
                linewidth = plot_estimate_linewidth)
  }
  model_plot <- model_plot +
    geom_point(data = emmeans_output_df,
               aes(x = condition_factor, y = estimated_response),
               colour = estimate_color,
               size = plot_estimate_point_size) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
    create.benchmark.plot.theme()
  list(model = model_output,
       emmeans = emmeans_output,
       emmeans_df = emmeans_output_df,
       plot = model_plot)
}
plot.categorical.gaussian.benchmark.model <- function(input_data,
                                                      condition_variable,
                                                      response_variable,
                                                      condition_levels,
                                                      estimate_color,
                                                      connect_estimates = FALSE) {
  model_data <- prepare.continuous.benchmark.data(input_data = input_data, condition_variable = condition_variable, response_variable = response_variable)
  model_data$condition_factor <- factor(as.character(model_data$condition), levels = as.character(condition_levels))
  model_data <- model_data[!is.na(model_data$condition_factor), , drop = FALSE]
  if (length(unique(model_data$sim)) < 2) stop("Gaussian mixed model requires at least two simulation replicates")
  model_output <- nlme::lme(response ~ condition_factor, random = ~ 1 | sim, data = model_data, na.action = na.omit, method = "REML")
  emmeans_output <- emmeans::emmeans(model_output, ~ condition_factor)
  emmeans_output_df <- as.data.frame(emmeans_output)
  emmeans_output_df$lower_SE <- emmeans_output_df$emmean - plot_SE_multiplier * emmeans_output_df$SE
  emmeans_output_df$upper_SE <- emmeans_output_df$emmean + plot_SE_multiplier * emmeans_output_df$SE
  model_plot <- ggplot() +
    geom_point(data = model_data,
               aes(x = condition_factor, y = response),
               alpha = plot_raw_point_alpha,
               size = plot_raw_point_size,
               position = position_jitter(width = plot_categorical_jitter_width, height = 0)) +
    geom_errorbar(data = emmeans_output_df,
                  aes(x = condition_factor, ymin = lower_SE, ymax = upper_SE),
                  colour = estimate_color,
                  linewidth = plot_errorbar_linewidth,
                  width = plot_errorbar_width)
  if (connect_estimates) {
    model_plot <- model_plot +
      geom_line(data = emmeans_output_df,
                aes(x = condition_factor, y = emmean, group = 1),
                colour = estimate_color,
                linewidth = plot_estimate_linewidth)
  }
  model_plot <- model_plot +
    geom_point(data = emmeans_output_df,
               aes(x = condition_factor, y = emmean),
               colour = estimate_color,
               size = plot_estimate_point_size) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
    create.benchmark.plot.theme()
  list(model = model_output,
       emmeans = emmeans_output,
       emmeans_df = emmeans_output_df,
       plot = model_plot)
}
format.benchmark.panel <- function(input_plot,
                                   test_title = NULL,
                                   metric_label = NULL,
                                   show_x_axis = FALSE,
                                   x_text_angle = 0,
                                   x_text_hjust = 1,
                                   x_text_vjust = 0.5,
                                   panel_tag = NULL) {
  output_plot <- input_plot +
    labs(title = test_title,
         x = NULL,
         y = metric_label,
         tag = panel_tag) +
    theme(axis.title.x = element_blank(),
          axis.text.y = element_text(family = plot_font_family, size = plot_axis_text_size),
          axis.ticks.y = element_line(linewidth = 0.3),
          axis.line.y = element_line(linewidth = 0.3))
  if (is.null(test_title)) output_plot <- output_plot + theme(plot.title = element_blank())
  if (is.null(metric_label)) {
    output_plot <- output_plot + theme(axis.title.y = element_blank())
  } else {
    output_plot <- output_plot + theme(axis.title.y = element_text(family = plot_font_family, size = plot_axis_title_size, face = "bold", margin = margin(r = 2, unit = "mm")))
  }
  if (show_x_axis) {
    output_plot <- output_plot +
      theme(axis.text.x = element_text(family = plot_font_family, size = plot_axis_text_size, angle = x_text_angle, hjust = x_text_hjust, vjust = x_text_vjust),
            axis.ticks.x = element_line(linewidth = 0.3),
            axis.line.x = element_line(linewidth = 0.3))
  } else {
    output_plot <- output_plot +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_line(linewidth = 0.3))
  }
  output_plot
}


## Load and check completed analyses
analysis_data <- list(
  clustering_methods = read.csv(sim_results_clustering_methods_csv, stringsAsFactors = FALSE),
  neighborhoods = read.csv(sim_results_neighborhoods_csv, stringsAsFactors = FALSE),
  N_steps = read.csv(sim_results_N_steps_csv, stringsAsFactors = FALSE),
  missing_data = read.csv(sim_results_NA_csv, stringsAsFactors = FALSE),
  learning_rate_tuning = read.csv(sim_results_learning_rate_tuning_csv, stringsAsFactors = FALSE)
)
for (test_name in names(test_specs)) {
  condition_variable <- test_specs[[test_name]]$condition_variable
  required_columns <- c("sim", condition_variable, "K_correct", "Acc")
  missing_columns <- setdiff(required_columns, colnames(analysis_data[[test_name]]))
  if (length(missing_columns) > 0) stop("The result table for ", test_specs[[test_name]]$title, " is missing the following columns: ", paste(missing_columns, collapse = ", "))
}


## Create top-row K-correct plots
K_correct_plots <- vector(mode = "list", length = length(test_specs))
K_correct_models <- vector(mode = "list", length = length(test_specs))
K_correct_emmeans <- vector(mode = "list", length = length(test_specs))
for (test_index in seq_along(test_specs)) {
  test_name <- names(test_specs)[test_index]
  test_spec <- test_specs[[test_name]]
  model_plot_output <- plot.categorical.binary.benchmark.model(
    input_data = analysis_data[[test_name]],
    condition_variable = test_spec$condition_variable,
    response_variable = "K_correct",
    condition_levels = test_spec$condition_levels,
    estimate_color = test_spec$estimate_color,
    connect_estimates = test_spec$connect_estimates
  )
  metric_label <- if (test_index == 1) "Proportion of correct K" else NULL
  K_correct_plots[[test_index]] <- format.benchmark.panel(input_plot = model_plot_output$plot,
                                                          test_title = test_spec$title,
                                                          metric_label = metric_label,
                                                          show_x_axis = FALSE,
                                                          panel_tag = NULL)
  K_correct_models[[test_index]] <- model_plot_output$model
  K_correct_emmeans[[test_index]] <- model_plot_output$emmeans
  names(K_correct_models)[test_index] <- test_name
  names(K_correct_emmeans)[test_index] <- test_name
}


## Create bottom-row assignment-accuracy plots
Acc_plots <- vector(mode = "list", length = length(test_specs))
Acc_models <- vector(mode = "list", length = length(test_specs))
Acc_emmeans <- vector(mode = "list", length = length(test_specs))
for (test_index in seq_along(test_specs)) {
  test_name <- names(test_specs)[test_index]
  test_spec <- test_specs[[test_name]]
  model_plot_output <- plot.categorical.gaussian.benchmark.model(input_data = analysis_data[[test_name]],
                                                                 condition_variable = test_spec$condition_variable,
                                                                 response_variable = "Acc",
                                                                 condition_levels = test_spec$condition_levels,
                                                                 estimate_color = test_spec$estimate_color,
                                                                 connect_estimates = test_spec$connect_estimates)
  metric_label <- if (test_index == 1) "Assignment accuracy" else NULL
  Acc_plots[[test_index]] <- format.benchmark.panel(input_plot = model_plot_output$plot,
                                                    metric_label = metric_label,
                                                    show_x_axis = TRUE,
                                                    x_text_angle = test_spec$x_text_angle,
                                                    x_text_hjust = test_spec$x_text_hjust,
                                                    x_text_vjust = test_spec$x_text_vjust,
                                                    panel_tag = NULL)
  Acc_models[[test_index]] <- model_plot_output$model
  Acc_emmeans[[test_index]] <- model_plot_output$emmeans
  names(Acc_models)[test_index] <- test_name
  names(Acc_emmeans)[test_index] <- test_name
}


## Combine plots using category-based column widths
combined_plot <- patchwork::wrap_plots(plotlist = c(K_correct_plots, Acc_plots),
                                       ncol = length(test_specs),
                                       nrow = 2,
                                       byrow = TRUE) +
  patchwork::plot_layout(widths = unname(column_width_units), heights = c(1, 1)) &
  theme(text = element_text(family = plot_font_family),
        plot.margin = margin(panel_margins_mm, panel_margins_mm, panel_margins_mm, panel_margins_mm, unit = "mm"))


## Show and save combined plot
combined_plot
ggsave(filename = combined_svg,
       plot = combined_plot,
       device = svglite::svglite,
       width = plot_width_cm,
       height = plot_height_cm,
       units = "cm",
       bg = "white",
       limitsize = FALSE)
message("Combined SVG saved to: ", combined_svg)




#### Supplementary Figure S3 ###################################################

## Set environment
rm(list = ls()) #clear environment
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")


## Install and load required packages
required_packages <- c("emmeans",
                       "glmmTMB",
                       "nlme",
                       "patchwork",
                       "svglite",
                       "tidyverse")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
  library(package_name, character.only = TRUE)
}


## Set main plotting parameters
plot_width_cm <- 16
plot_height_cm <- 20

plot_font_family <- "Arial"
plot_text_color <- "black"

plot_base_size <- 9
plot_axis_title_size <- 9
plot_axis_text_size <- 7.1
plot_column_title_size <- 9

plot_raw_point_alpha <- 0.35
plot_raw_point_size <- 1
plot_estimate_point_size <- 2.6
plot_estimate_linewidth <- 0.7
plot_errorbar_linewidth <- 1.3
plot_errorbar_width <- 0
plot_categorical_jitter_width <- 0.06
plot_SE_multiplier <- 1
X_text_angle_plot <- 90
minimum_column_width_units <- 3.5
panel_margins_mm <- 2.5

plot_color_clustering_methods <- "#56B4E9"
plot_color_neighborhood <- "#009E73"
plot_color_N_steps <- "#E69F00"
plot_color_NA <- "#D55E00"
plot_color_learning_rate_tuning <- "#CC80A7"


## Set input and output
simulation_k7_dir <- file.path("Simulations", "Simulation_set_1", "k7")
sim_results_clustering_methods_csv <- "Sim_results_clustering_methods.csv"
sim_results_neighborhoods_csv <- "Sim_results_neighborhoods.csv"
sim_results_N_steps_csv <- "Sim_results_N_steps.csv"
sim_results_NA_csv <- "Sim_results_NA.csv"
sim_results_learning_rate_tuning_csv <- "Sim_results_learning_rate_tuning.csv"

plot_dir <- "Figure_files"
plot_file_name <- "Supplementary_Figure_S3.svg"


## Combine directories and file names
sim_results_clustering_methods_csv <- file.path(simulation_k7_dir, sim_results_clustering_methods_csv)
sim_results_neighborhoods_csv <- file.path(simulation_k7_dir, sim_results_neighborhoods_csv)
sim_results_N_steps_csv <- file.path(simulation_k7_dir, sim_results_N_steps_csv)
sim_results_NA_csv <- file.path(simulation_k7_dir, sim_results_NA_csv)
sim_results_learning_rate_tuning_csv <- file.path(simulation_k7_dir, sim_results_learning_rate_tuning_csv)
combined_svg <- file.path(plot_dir, plot_file_name)


## Create output directory if needed
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


## Check that all completed result files exist
required_result_files <- c(sim_results_clustering_methods_csv,
                           sim_results_neighborhoods_csv,
                           sim_results_N_steps_csv,
                           sim_results_NA_csv,
                           sim_results_learning_rate_tuning_csv)
missing_result_files <- required_result_files[!file.exists(required_result_files)]
if (length(missing_result_files) > 0) stop("The following completed result files are missing: ", paste(missing_result_files, collapse = ", "))


## Define tests
test_specs <- list(
  clustering_methods = list(
    title = "Clustering\nmethod",
    condition_variable = "clustering_method",
    condition_levels = c("kmeans+BICelbow",
                         "kmeans+BICthreshold",
                         "GMM+BICthreshold",
                         "hierarchical+DB",
                         "HDBSCAN",
                         "OPTICS+Silhouette"),
    estimate_color = plot_color_clustering_methods,
    connect_estimates = FALSE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  neighborhoods = list(
    title = "Neighborhood\nfunction",
    condition_variable = "neighborhood_function",
    condition_levels = c("gaussian", "bubble"),
    estimate_color = plot_color_neighborhood,
    connect_estimates = FALSE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  N_steps = list(
    title = "Number of training\niterations",
    condition_variable = "N_steps",
    condition_levels = c(20,
                         50,
                         100,
                         200,
                         500,
                         1000,
                         2000),
    estimate_color = plot_color_N_steps,
    connect_estimates = TRUE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  learning_rate_tuning = list(
    title = "Learning-rate\ntuning",
    condition_variable = "learning_rate_tuning",
    condition_levels = c(FALSE, TRUE),
    estimate_color = plot_color_learning_rate_tuning,
    connect_estimates = FALSE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  ),
  missing_data = list(
    title = "Missing data\nproportion",
    condition_variable = "missing_data_prop",
    condition_levels = c(0.0,
                         0.1,
                         0.2,
                         0.3,
                         0.4,
                         0.5,
                         0.6,
                         0.7,
                         0.8),
    estimate_color = plot_color_NA,
    connect_estimates = TRUE,
    x_text_angle = X_text_angle_plot,
    x_text_hjust = 1,
    x_text_vjust = 0.5
  )
)


## Define metrics
metric_specs <- list(
  ARI = list(
    response_variable = "ARI",
    metric_label = "Adjusted Rand Index"
  ),
  QE = list(
    response_variable = "QE",
    metric_label = "Quantization error"
  ),
  TE = list(
    response_variable = "TE",
    metric_label = "Topographic error"
  ),
  Time = list(
    response_variable = "Time",
    metric_label = "Time (seconds)"
  )
)


## Calculate category-based column widths
test_category_counts <- vapply(test_specs, function(test_spec) length(test_spec$condition_levels), numeric(1))
column_width_units <- pmax(test_category_counts, minimum_column_width_units)


## Create additional functions
safe.mean <- function(input_vector) {
  if (all(is.na(input_vector))) return(NA_real_)
  mean(input_vector, na.rm = TRUE)
}
prepare.continuous.benchmark.data <- function(input_data,
                                              condition_variable,
                                              response_variable) {
  model_data <- input_data[ , c("sim", condition_variable, response_variable), drop = FALSE]
  colnames(model_data)[colnames(model_data) == condition_variable] <- "condition"
  colnames(model_data)[colnames(model_data) == response_variable] <- "response"
  model_data$response <- suppressWarnings(as.numeric(model_data$response))
  model_data <- model_data[!is.na(model_data$condition) & !is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data <- model_data %>%
    dplyr::group_by(sim, condition) %>%
    dplyr::summarise(response = safe.mean(response), .groups = "drop")
  model_data <- model_data[!is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data
}


## Create plotting functions
create.benchmark.plot.theme <- function() {
  theme_classic(base_size = plot_base_size, base_family = plot_font_family) +
    theme(text = element_text(family = plot_font_family, colour = plot_text_color),
          axis.title = element_text(family = plot_font_family, size = plot_axis_title_size),
          axis.text = element_text(family = plot_font_family, size = plot_axis_text_size),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          axis.ticks.length = grid::unit(0.8, "mm"),
          axis.text.x.top = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.line.x.top = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          plot.title = element_text(family = plot_font_family,
                                    size = plot_column_title_size,
                                    face = "bold",
                                    hjust = 0.5,
                                    lineheight = 0.9,
                                    margin = margin(b = 1.2, unit = "mm")),
          plot.margin = margin(1, 1, 1, 1, unit = "mm"))
}
plot.categorical.gaussian.benchmark.model <- function(input_data,
                                                      condition_variable,
                                                      response_variable,
                                                      condition_levels,
                                                      estimate_color,
                                                      connect_estimates = FALSE) {
  model_data <- prepare.continuous.benchmark.data(input_data = input_data, condition_variable = condition_variable, response_variable = response_variable)
  model_data$condition_factor <- factor(as.character(model_data$condition), levels = as.character(condition_levels))
  model_data <- model_data[!is.na(model_data$condition_factor), , drop = FALSE]
  if (length(unique(model_data$sim)) < 2) stop("Gaussian mixed model requires at least two simulation replicates")
  model_output <- nlme::lme(response ~ condition_factor, random = ~ 1 | sim, data = model_data, na.action = na.omit, method = "REML")
  emmeans_output <- emmeans::emmeans(model_output, ~ condition_factor)
  emmeans_output_df <- as.data.frame(emmeans_output)
  emmeans_output_df$lower_SE <- emmeans_output_df$emmean - plot_SE_multiplier * emmeans_output_df$SE
  emmeans_output_df$upper_SE <- emmeans_output_df$emmean + plot_SE_multiplier * emmeans_output_df$SE
  model_plot <- ggplot() +
    geom_point(data = model_data,
               aes(x = condition_factor, y = response),
               alpha = plot_raw_point_alpha,
               size = plot_raw_point_size,
               position = position_jitter(width = plot_categorical_jitter_width, height = 0)) +
    geom_errorbar(data = emmeans_output_df,
                  aes(x = condition_factor, ymin = lower_SE, ymax = upper_SE),
                  colour = estimate_color,
                  linewidth = plot_errorbar_linewidth,
                  width = plot_errorbar_width)
  if (connect_estimates) {
    model_plot <- model_plot +
      geom_line(data = emmeans_output_df,
                aes(x = condition_factor, y = emmean, group = 1),
                colour = estimate_color,
                linewidth = plot_estimate_linewidth)
  }
  model_plot <- model_plot +
    geom_point(data = emmeans_output_df,
               aes(x = condition_factor, y = emmean),
               colour = estimate_color,
               size = plot_estimate_point_size) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
    create.benchmark.plot.theme()
  list(model = model_output,
       emmeans = emmeans_output,
       emmeans_df = emmeans_output_df,
       plot = model_plot)
}
format.benchmark.panel <- function(input_plot,
                                   test_title = NULL,
                                   metric_label = NULL,
                                   show_x_axis = FALSE,
                                   x_text_angle = 0,
                                   x_text_hjust = 1,
                                   x_text_vjust = 0.5,
                                   panel_tag = NULL) {
  output_plot <- input_plot +
    labs(title = test_title,
         x = NULL,
         y = metric_label,
         tag = panel_tag) +
    theme(axis.title.x = element_blank(),
          axis.text.y = element_text(family = plot_font_family, size = plot_axis_text_size),
          axis.ticks.y = element_line(linewidth = 0.3),
          axis.line.y = element_line(linewidth = 0.3))
  if (is.null(test_title)) output_plot <- output_plot + theme(plot.title = element_blank())
  if (is.null(metric_label)) {
    output_plot <- output_plot + theme(axis.title.y = element_blank())
  } else {
    output_plot <- output_plot + theme(axis.title.y = element_text(family = plot_font_family, size = plot_axis_title_size, face = "bold", margin = margin(r = 2, unit = "mm")))
  }
  if (show_x_axis) {
    output_plot <- output_plot +
      theme(axis.text.x = element_text(family = plot_font_family, size = plot_axis_text_size, angle = x_text_angle, hjust = x_text_hjust, vjust = x_text_vjust),
            axis.ticks.x = element_line(linewidth = 0.3),
            axis.line.x = element_line(linewidth = 0.3))
  } else {
    output_plot <- output_plot +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_line(linewidth = 0.3))
  }
  output_plot
}


## Load and check completed analyses
analysis_data <- list(
  clustering_methods = read.csv(sim_results_clustering_methods_csv, stringsAsFactors = FALSE),
  neighborhoods = read.csv(sim_results_neighborhoods_csv, stringsAsFactors = FALSE),
  N_steps = read.csv(sim_results_N_steps_csv, stringsAsFactors = FALSE),
  missing_data = read.csv(sim_results_NA_csv, stringsAsFactors = FALSE),
  learning_rate_tuning = read.csv(sim_results_learning_rate_tuning_csv, stringsAsFactors = FALSE)
)
required_metric_columns <- vapply(metric_specs, function(metric_spec) metric_spec$response_variable, character(1))
for (test_name in names(test_specs)) {
  condition_variable <- test_specs[[test_name]]$condition_variable
  required_columns <- c("sim", condition_variable, required_metric_columns)
  missing_columns <- setdiff(required_columns, colnames(analysis_data[[test_name]]))
  if (length(missing_columns) > 0) stop("The result table for ", test_specs[[test_name]]$title, " is missing the following columns: ", paste(missing_columns, collapse = ", "))
}


## Create plots for all remaining metrics
metric_plots <- lapply(metric_specs, function(metric_spec) vector(mode = "list", length = length(test_specs)))
metric_models <- lapply(metric_specs, function(metric_spec) list())
metric_emmeans <- lapply(metric_specs, function(metric_spec) list())
for (metric_index in seq_along(metric_specs)) {
  metric_name <- names(metric_specs)[metric_index]
  metric_spec <- metric_specs[[metric_name]]
  for (test_index in seq_along(test_specs)) {
    test_name <- names(test_specs)[test_index]
    test_spec <- test_specs[[test_name]]
    model_plot_output <- plot.categorical.gaussian.benchmark.model(
      input_data = analysis_data[[test_name]],
      condition_variable = test_spec$condition_variable,
      response_variable = metric_spec$response_variable,
      condition_levels = test_spec$condition_levels,
      estimate_color = test_spec$estimate_color,
      connect_estimates = test_spec$connect_estimates
    )
    test_title <- if (metric_index == 1) test_spec$title else NULL
    metric_label <- if (test_index == 1) metric_spec$metric_label else NULL
    show_x_axis <- metric_index == length(metric_specs)
    metric_plots[[metric_name]][[test_index]] <- format.benchmark.panel(
      input_plot = model_plot_output$plot,
      test_title = test_title,
      metric_label = metric_label,
      show_x_axis = show_x_axis,
      x_text_angle = test_spec$x_text_angle,
      x_text_hjust = test_spec$x_text_hjust,
      x_text_vjust = test_spec$x_text_vjust,
      panel_tag = NULL
    )
    metric_models[[metric_name]][[test_name]] <- model_plot_output$model
    metric_emmeans[[metric_name]][[test_name]] <- model_plot_output$emmeans
  }
}


## Combine plots using category-based column widths
combined_plot_list <- do.call(c, metric_plots)
combined_plot <- patchwork::wrap_plots(plotlist = combined_plot_list,
                                       ncol = length(test_specs),
                                       nrow = length(metric_specs),
                                       byrow = TRUE) +
  patchwork::plot_layout(widths = unname(column_width_units), heights = rep(1, length(metric_specs))) &
  theme(text = element_text(family = plot_font_family),
        plot.margin = margin(panel_margins_mm, panel_margins_mm, panel_margins_mm, panel_margins_mm, unit = "mm"))


## Show and save combined plot
combined_plot
ggsave(filename = combined_svg,
       plot = combined_plot,
       device = svglite::svglite,
       width = plot_width_cm,
       height = plot_height_cm,
       units = "cm",
       bg = "white",
       limitsize = FALSE)
message("Combined SVG saved to: ", combined_svg)




#### Supplementary Figure S4 ###################################################


## Set environment
rm(list = ls())
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")


## Install and load required packages
required_packages <- c("clue",
                       "emmeans",
                       "glmmTMB",
                       "kohonen",
                       "MASS",
                       "mclust",
                       "nlme",
                       "sn",
                       "tidyverse",
                       "viridis")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
  library(package_name, character.only = TRUE)
}


## Set parameters
simulation_k7_dir <- file.path("Simulations", "Simulation_set_1", "k7")
sim_data_base_file <- file.path(simulation_k7_dir, "Sim_data_base_complete.rds")
intermediate_results_dir <- file.path(simulation_k7_dir, "Intermediate_results")

overwrite_results <- FALSE

plot_width_cm <- 8.27
plot_height_cm <- 8.14
plot_font_family <- "Arial"
plot_text_color <- "black"
plot_base_size <- 9
plot_axis_title_size <- 9
plot_axis_text_size <- 7.1
plot_tick_number_axis_gap <- 0.65
plot_axis_title_number_gap <- 1.33

plot_dir <- "Figure_files"
plot_BIC_file_name <- "Supplementary_Figure_S4A.svg"
plot_delta_BIC_file_name <- "Supplementary_Figure_S4B.svg"

N_individuals <- 140
N_steps_SOM <- 100
N_replicates_SOM <- 110
SOM_grid_size <- rep(floor(sqrt(4 * sqrt(N_individuals))), 2)
max_NA_row_SOM <- 1
max_NA_col_SOM <- 1
BIC_threshold_SOM <- 6
learning_rate_tuning <- FALSE
max_k_SOM <- 12


## Load K = 7 simulations
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!file.exists(sim_data_base_file)) stop("K = 7 base simulation file not found: ", sim_data_base_file)
simulation_results_base <- readRDS(sim_data_base_file)


## Select random simulation replicate
set.seed(1)
simulation_index <- sample(seq_along(simulation_results_base), 1)
simulation_index
simulation_data <- simulation_results_base[[simulation_index]]


## Train and cluster SOM
som_input <- list(SNP = simulation_data$SNP,
                  Morphology = simulation_data$Morphology,
                  Climate = simulation_data$Climate,
                  Host = simulation_data$Host)
som_output <- train.SOM(input_data = som_input,
                        overwrite.SOM.results = overwrite_results,
                        learning.rate.tuning = learning_rate_tuning,
                        training.neighborhoods = "gaussian",
                        N.steps = N_steps_SOM,
                        grid.size = SOM_grid_size,
                        N.replicates = N_replicates_SOM,
                        max.NA.row = max_NA_row_SOM,
                        max.NA.col = max_NA_col_SOM,
                        save.SOM.results.name = file.path(intermediate_results_dir, paste0("kmeans+BICelbow_k7_BIC_test_sim", simulation_index, "_training.Rdata")),
                        save.SOM.results = TRUE)
clustering_output <- clustering.SOM(SOM.output = som_output,
                                    overwrite.SOM.results = overwrite_results,
                                    max.k = max_k_SOM,
                                    BIC.thresh = BIC_threshold_SOM,
                                    save.SOM.results = TRUE,
                                    save.SOM.results.name = file.path(intermediate_results_dir, paste0("kmeans+BICelbow_k7_BIC_test_sim", simulation_index, "_clustering.Rdata")),
                                    clustering.method = "kmeans+BICelbow")
clustering_output$optim_k_summary


## Prepare BIC values
if (is.null(clustering_output$BIC_values)) stop("clustering.SOM did not return BIC_values")
BIC_values <- clustering_output$BIC_values[seq_len(max_k_SOM), , drop = FALSE]
k_colors <- viridis::magma(max_k_SOM)
finite_BIC_rows <- apply(BIC_values, 1, function(BIC_values_for_K) any(is.finite(BIC_values_for_K) & !is.na(BIC_values_for_K)))
plotted_BIC_k_values <- seq_len(max_k_SOM)[finite_BIC_rows]
BIC_values_for_plot <- t(BIC_values[finite_BIC_rows, , drop = FALSE])


## Calculate delta BIC
delta_BIC_matrix <- apply(BIC_values, 2, function(BIC_replicate_values) {
  previous_BIC_values <- BIC_replicate_values[-length(BIC_replicate_values)]
  current_BIC_values <- BIC_replicate_values[-1]
  delta_BIC_values <- previous_BIC_values - current_BIC_values
  delta_BIC_values[!is.finite(previous_BIC_values) | !is.finite(current_BIC_values)] <- NA_real_
  return(delta_BIC_values)
})
if (is.null(dim(delta_BIC_matrix))) delta_BIC_matrix <- matrix(delta_BIC_matrix, ncol = 1)
rownames(delta_BIC_matrix) <- paste0("k", 2:max_k_SOM, "-k", 1:(max_k_SOM - 1))
finite_delta_BIC_rows <- apply(delta_BIC_matrix, 1, function(delta_BIC_values_for_K) any(is.finite(delta_BIC_values_for_K) & !is.na(delta_BIC_values_for_K)))
plotted_delta_BIC_k_values <- seq.int(2, max_k_SOM)[finite_delta_BIC_rows]
delta_BIC_values_for_plot <- t(delta_BIC_matrix[finite_delta_BIC_rows, , drop = FALSE])


## Plot BIC
svg_scaling_factor <- 96 / 72
svg(file.path(plot_dir, plot_BIC_file_name), width = (plot_width_cm / 2.54) * svg_scaling_factor, height = (plot_height_cm / 2.54) * svg_scaling_factor, pointsize = plot_base_size, family = plot_font_family)
par(family = plot_font_family, fg = plot_text_color, col.axis = plot_text_color, bty = "n", mar = c(4, 4.5, 1, 1))
boxplot(BIC_values_for_plot,
        at = plotted_BIC_k_values,
        xlim = c(0.5, max_k_SOM + 0.5),
        outline = FALSE,
        notch = FALSE,
        axes = FALSE,
        ylab = "",
        main = "",
        whisklty = 1,
        staplelty = 1,
        col = k_colors[plotted_BIC_k_values])
axis(1, at = seq_len(max_k_SOM), labels = seq_len(max_k_SOM), mgp = c(3, plot_tick_number_axis_gap, 0), cex.axis = (plot_axis_text_size * svg_scaling_factor) / plot_base_size, col = plot_text_color, col.axis = plot_text_color)
BIC_y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = BIC_y_axis_breaks, labels = round(BIC_y_axis_breaks, 1), las = 3, mgp = c(3, plot_tick_number_axis_gap, 0), cex.axis = (plot_axis_text_size * svg_scaling_factor) / plot_base_size, col = plot_text_color, col.axis = plot_text_color)
mtext("BIC", side = 2, line = plot_tick_number_axis_gap + plot_axis_title_number_gap, font = 2, cex = (plot_axis_title_size * svg_scaling_factor) / plot_base_size, col = plot_text_color)
dev.off()


## Plot delta BIC
svg(file.path(plot_dir, plot_delta_BIC_file_name), width = (plot_width_cm / 2.54) * svg_scaling_factor, height = (plot_height_cm / 2.54) * svg_scaling_factor, pointsize = plot_base_size, family = plot_font_family)
par(family = plot_font_family, fg = plot_text_color, col.axis = plot_text_color, bty = "n", mar = c(4, 4.5, 1, 1))
boxplot(delta_BIC_values_for_plot,
        at = plotted_delta_BIC_k_values,
        xlim = c(0.5, max_k_SOM + 0.5),
        outline = FALSE,
        notch = FALSE,
        axes = FALSE,
        ylab = "",
        main = "",
        whisklty = 1,
        staplelty = 1,
        col = k_colors[plotted_delta_BIC_k_values])
axis(1, at = seq_len(max_k_SOM), labels = seq_len(max_k_SOM), mgp = c(3, plot_tick_number_axis_gap, 0), cex.axis = (plot_axis_text_size * svg_scaling_factor) / plot_base_size, col = plot_text_color, col.axis = plot_text_color)
delta_BIC_y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = delta_BIC_y_axis_breaks, labels = round(delta_BIC_y_axis_breaks, 1), las = 3, mgp = c(3, plot_tick_number_axis_gap, 0), cex.axis = (plot_axis_text_size * svg_scaling_factor) / plot_base_size, col = plot_text_color, col.axis = plot_text_color)
mtext("delta BIC", side = 2, line = plot_tick_number_axis_gap + plot_axis_title_number_gap, font = 2, cex = (plot_axis_title_size * svg_scaling_factor) / plot_base_size, col = plot_text_color)
dev.off()




#### Supplementary Figure S5 ##################################################################

## Set environment
rm(list = ls()) #clear environment
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
invisible(gc())


## Install and load required packages
required_packages <- c("svglite", "tidyverse", "viridisLite")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
  library(package_name, character.only = TRUE)
}


## Set main plotting parameters
plot_width_cm <- 15.36
plot_height_cm <- 19.24

plot_font_family <- "Arial"
plot_text_color <- "black"

plot_base_size <- 9
plot_axis_title_size <- 9
plot_axis_text_size <- 7.1
plot_column_title_size <- 9
plot_row_title_size <- 9
plot_legend_title_size <- 9
plot_legend_text_size <- 9

plot_axis_title_margin_mm <- 3
legend_box_margin_mm <- 1
panel_spacing_mm <- 2.5
plot_column_title_margin_mm <- 2.5
plot_x_axis_left_expansion <- 10000
plot_x_axis_right_expansion <- 10000

plot_point_size <- 1
plot_point_alpha <- 0.7
plot_fitted_line_width <- 1
plot_threshold_line_width <- 0.6
plot_threshold_line_type <- "dashed"
plot_K2_threshold_proportion <- 0.5
plot_reference_line_positions <- c(50000, 100000, 150000, 200000)
plot_reference_line_width <- 0.1
legend_box_line_width <- 0.3
migration_values <- c(0, 1e-6, 4e-6, 7e-6)
migration_labels <- c("0", "1e-6", "4e-6", "7e-6")
migration_colors <- stats::setNames(c("#2B0A3D", "#31668EFF", "#30B57BFF", "#F0DD1CFF"), migration_labels)

## Set input and output
simulation_dir <- file.path("Simulations", "Simulation_set_2")
results_root_dir <- file.path(simulation_dir, "fastsimcoal2_results")
balanced_results_dir <- file.path(results_root_dir, "symmetric_16")
unbalanced_results_dir <- file.path(results_root_dir, "asymmetric_16")
balanced_SOM_results_csv <- file.path(balanced_results_dir, "fastsimcoal2_SOM_results.csv")
balanced_SOM_optim_k_csv <- file.path(balanced_results_dir, "fastsimcoal2_SOM_optim_k_summary.csv")
balanced_STRUCTURE_csv <- file.path(balanced_results_dir, "STRUCTURE_results", "STRUCTURE_binomial_table_by_vcf.csv")
unbalanced_SOM_results_csv <- file.path(unbalanced_results_dir, "fastsimcoal2_SOM_results.csv")
unbalanced_SOM_optim_k_csv <- file.path(unbalanced_results_dir, "fastsimcoal2_SOM_optim_k_summary.csv")
unbalanced_STRUCTURE_csv <- file.path(unbalanced_results_dir, "STRUCTURE_results", "STRUCTURE_binomial_table_by_vcf.csv")
plot_dir <- "Figure_files"
plot_file_name <- "Supplementary_figure_S5.svg"
combined_svg <- file.path(plot_dir, plot_file_name)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


## Check that all completed result files exist
required_result_files <- c(balanced_SOM_results_csv,
                           balanced_SOM_optim_k_csv,
                           balanced_STRUCTURE_csv,
                           unbalanced_SOM_results_csv,
                           unbalanced_SOM_optim_k_csv,
                           unbalanced_STRUCTURE_csv)
missing_result_files <- required_result_files[!file.exists(required_result_files)]
if (length(missing_result_files) > 0) stop("The following completed result files are missing: ", paste(missing_result_files, collapse = ", "))


## Create additional functions
check.required.columns <- function(input_table,
                                   required_columns,
                                   table_name) {
  missing_columns <- setdiff(required_columns, colnames(input_table))
  if (length(missing_columns) > 0) stop(table_name, " is missing the following columns: ", paste(missing_columns, collapse = ", "))
}
standardize.migration.labels <- function(input_data) {
  migration_indices <- vapply(input_data$mig, function(current_migration) {
    migration_difference <- abs(migration_values - current_migration)
    matching_index <- which.min(migration_difference)
    if (length(matching_index) != 1 || !is.finite(migration_difference[matching_index]) || migration_difference[matching_index] > sqrt(.Machine$double.eps)) return(NA_integer_)
    matching_index
  }, integer(1))
  if (any(is.na(migration_indices))) stop("At least one migration rate could not be matched to the expected migration rates")
  input_data$mig.tag <- factor(migration_labels[migration_indices], levels = migration_labels)
  input_data
}
create.analysis.binomial.table <- function(SOM_results_csv,
                                           SOM_optim_k_csv,
                                           STRUCTURE_csv,
                                           design_label) {
  result_table <- read.csv(SOM_results_csv, stringsAsFactors = FALSE)
  optim_k_result_table <- read.csv(SOM_optim_k_csv, stringsAsFactors = FALSE)
  STRUCTURE_binomial_table <- read.csv(STRUCTURE_csv, stringsAsFactors = FALSE)
  check.required.columns(result_table, c("file", "status", "mig", "tdiv", "deNovo.kmeans.best.k", "sNMF.best.k"), paste0(design_label, " SOM result table"))
  check.required.columns(optim_k_result_table, c("file", "mig", "tdiv", "Count", "k.label"), paste0(design_label, " SOM optim-k table"))
  check.required.columns(STRUCTURE_binomial_table, c("file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2"), paste0(design_label, " STRUCTURE binomial table"))
  plot_result_table <- result_table[result_table$status == "ok", ]
  plot_result_table <- plot_result_table[order(plot_result_table$mig, plot_result_table$tdiv), ]
  if (nrow(plot_result_table) == 0) stop("No successful SOM result rows are available for ", design_label)
  SOM_optim_k_count_table <- optim_k_result_table[!is.na(optim_k_result_table$Count) & optim_k_result_table$file %in% plot_result_table$file, ]
  SOM_optim_k_count_table$Count <- as.numeric(as.character(SOM_optim_k_count_table$Count))
  SOM_total_count_table <- stats::aggregate(Count ~ file + mig + tdiv, data = SOM_optim_k_count_table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE))
  colnames(SOM_total_count_table)[colnames(SOM_total_count_table) == "Count"] <- "n.total"
  SOM_k2_count_table <- SOM_optim_k_count_table[SOM_optim_k_count_table$k.label == "k2", ]
  if (nrow(SOM_k2_count_table) > 0) {
    SOM_k2_count_table <- stats::aggregate(Count ~ file + mig + tdiv, data = SOM_k2_count_table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE))
    colnames(SOM_k2_count_table)[colnames(SOM_k2_count_table) == "Count"] <- "n.k2"
  } else {
    SOM_k2_count_table <- SOM_total_count_table[, c("file", "mig", "tdiv")]
    SOM_k2_count_table$n.k2 <- 0
  }
  SOM_binomial_table <- merge(SOM_total_count_table, SOM_k2_count_table, by = c("file", "mig", "tdiv"), all.x = TRUE)
  SOM_binomial_table$n.k2[is.na(SOM_binomial_table$n.k2)] <- 0
  SOM_binomial_table$n.not.k2 <- SOM_binomial_table$n.total - SOM_binomial_table$n.k2
  SOM_binomial_table$proportion.k2 <- SOM_binomial_table$n.k2 / SOM_binomial_table$n.total
  SOM_binomial_table$method <- "SOM"
  SOM_binomial_table <- SOM_binomial_table[, c("method", "file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")]
  DAPC_binomial_source_table <- result_table[result_table$status == "ok" & !is.na(result_table$deNovo.kmeans.best.k), ]
  if ("deNovo.kmeans.status" %in% colnames(DAPC_binomial_source_table)) DAPC_binomial_source_table <- DAPC_binomial_source_table[DAPC_binomial_source_table$deNovo.kmeans.status == "ok", ]
  DAPC_best_k <- as.integer(as.character(DAPC_binomial_source_table$deNovo.kmeans.best.k))
  DAPC_binomial_table <- data.frame(method = "DAPC",
                                    file = DAPC_binomial_source_table$file,
                                    mig = DAPC_binomial_source_table$mig,
                                    tdiv = DAPC_binomial_source_table$tdiv,
                                    n.k2 = as.integer(DAPC_best_k == 2L),
                                    n.not.k2 = as.integer(DAPC_best_k != 2L),
                                    proportion.k2 = as.numeric(DAPC_best_k == 2L),
                                    stringsAsFactors = FALSE)
  sNMF_binomial_source_table <- result_table[result_table$status == "ok" & !is.na(result_table$sNMF.best.k), ]
  if ("sNMF.status" %in% colnames(sNMF_binomial_source_table)) sNMF_binomial_source_table <- sNMF_binomial_source_table[sNMF_binomial_source_table$sNMF.status == "ok", ]
  sNMF_best_k <- as.integer(as.character(sNMF_binomial_source_table$sNMF.best.k))
  sNMF_binomial_table <- data.frame(method = "sNMF",
                                    file = sNMF_binomial_source_table$file,
                                    mig = sNMF_binomial_source_table$mig,
                                    tdiv = sNMF_binomial_source_table$tdiv,
                                    n.k2 = as.integer(sNMF_best_k == 2L),
                                    n.not.k2 = as.integer(sNMF_best_k != 2L),
                                    proportion.k2 = as.numeric(sNMF_best_k == 2L),
                                    stringsAsFactors = FALSE)
  STRUCTURE_binomial_table <- STRUCTURE_binomial_table[, c("file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")]
  STRUCTURE_binomial_table$method <- "STRUCTURE"
  STRUCTURE_binomial_table <- STRUCTURE_binomial_table[, c("method", "file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")]
  combined_binomial_table <- rbind(SOM_binomial_table, DAPC_binomial_table, sNMF_binomial_table, STRUCTURE_binomial_table)
  combined_binomial_table$design <- design_label
  combined_binomial_table <- standardize.migration.labels(combined_binomial_table)
  combined_binomial_table
}
fit.K2.binomial.models <- function(input_data) {
  prediction_list <- list()
  threshold_list <- list()
  prediction_index <- 1
  threshold_index <- 1
  for (current_design in levels(input_data$design)) {
    for (current_method in levels(input_data$method)) {
      for (current_migration in levels(input_data$mig.tag)) {
        current_binomial_table <- input_data[input_data$design == current_design & input_data$method == current_method & input_data$mig.tag == current_migration, ]
        current_binomial_table <- current_binomial_table[is.finite(current_binomial_table$tdiv), ]
        if (nrow(current_binomial_table) == 0) next
        current_binomial_table <- current_binomial_table[order(current_binomial_table$tdiv), ]
        prediction_tdiv <- seq(min(current_binomial_table$tdiv), max(current_binomial_table$tdiv), length.out = 1000)
        if (all(current_binomial_table$n.k2 == 0)) {
          current_prediction_table <- data.frame(design = current_design,
                                                 method = current_method,
                                                 mig.tag = current_migration,
                                                 tdiv = prediction_tdiv,
                                                 fitted.proportion.k2 = 0,
                                                 stringsAsFactors = FALSE)
        } else if (all(current_binomial_table$n.not.k2 == 0)) {
          current_prediction_table <- data.frame(design = current_design,
                                                 method = current_method,
                                                 mig.tag = current_migration,
                                                 tdiv = prediction_tdiv,
                                                 fitted.proportion.k2 = 1,
                                                 stringsAsFactors = FALSE)
        } else {
          current_fit <- stats::glm(cbind(n.k2, n.not.k2) ~ tdiv, data = current_binomial_table, family = stats::binomial(link = "logit"), control = stats::glm.control(maxit = 100))
          current_prediction_table <- data.frame(design = current_design,
                                                 method = current_method,
                                                 mig.tag = current_migration,
                                                 tdiv = prediction_tdiv,
                                                 fitted.proportion.k2 = stats::predict(current_fit, newdata = data.frame(tdiv = prediction_tdiv), type = "response"),
                                                 stringsAsFactors = FALSE)
          current_coefficients <- stats::coef(current_fit)
          current_intercept <- unname(current_coefficients["(Intercept)"])
          current_tdiv_slope <- unname(current_coefficients["tdiv"])
          if (is.finite(current_intercept) && is.finite(current_tdiv_slope) && current_tdiv_slope != 0) {
            current_threshold_tdiv <- (stats::qlogis(plot_K2_threshold_proportion) - current_intercept) / current_tdiv_slope
            threshold_within_range <- current_threshold_tdiv >= min(current_binomial_table$tdiv) && current_threshold_tdiv <= max(current_binomial_table$tdiv)
            if (threshold_within_range) {
              threshold_list[[threshold_index]] <- data.frame(design = current_design,
                                                              method = current_method,
                                                              mig.tag = current_migration,
                                                              tdiv.at.threshold.k2 = current_threshold_tdiv,
                                                              stringsAsFactors = FALSE)
              threshold_index <- threshold_index + 1
            }
          }
        }
        prediction_list[[prediction_index]] <- current_prediction_table
        prediction_index <- prediction_index + 1
      }
    }
  }
  prediction_table <- do.call(rbind, prediction_list)
  if (length(threshold_list) > 0) {
    threshold_table <- do.call(rbind, threshold_list)
  } else {
    threshold_table <- data.frame(design = character(0), method = character(0), mig.tag = character(0), tdiv.at.threshold.k2 = numeric(0))
  }
  prediction_table$design <- factor(prediction_table$design, levels = levels(input_data$design))
  prediction_table$method <- factor(prediction_table$method, levels = levels(input_data$method))
  prediction_table$mig.tag <- factor(prediction_table$mig.tag, levels = migration_labels)
  threshold_table$design <- factor(threshold_table$design, levels = levels(input_data$design))
  threshold_table$method <- factor(threshold_table$method, levels = levels(input_data$method))
  threshold_table$mig.tag <- factor(threshold_table$mig.tag, levels = migration_labels)
  list(prediction_table = prediction_table, threshold_table = threshold_table)
}


## Load balanced analysis
balanced_binomial_table <- create.analysis.binomial.table(SOM_results_csv = balanced_SOM_results_csv,
                                                          SOM_optim_k_csv = balanced_SOM_optim_k_csv,
                                                          STRUCTURE_csv = balanced_STRUCTURE_csv,
                                                          design_label = "Balanced")


## Load unbalanced analysis
unbalanced_binomial_table <- create.analysis.binomial.table(SOM_results_csv = unbalanced_SOM_results_csv,
                                                            SOM_optim_k_csv = unbalanced_SOM_optim_k_csv,
                                                            STRUCTURE_csv = unbalanced_STRUCTURE_csv,
                                                            design_label = "Unbalanced")


## Combine balanced and unbalanced analyses
combined_binomial_table <- rbind(balanced_binomial_table, unbalanced_binomial_table)
combined_binomial_table$design <- factor(combined_binomial_table$design, levels = c("Balanced", "Unbalanced"))
combined_binomial_table$method <- factor(combined_binomial_table$method, levels = c("SOM", "DAPC", "sNMF", "STRUCTURE"))


## Fit binomial models
fitted_K2_output <- fit.K2.binomial.models(combined_binomial_table)
fitted_prediction_table <- fitted_K2_output$prediction_table
threshold_table <- fitted_K2_output$threshold_table


## Create figure
combined_plot <- ggplot() +
  geom_vline(xintercept = plot_reference_line_positions,
             linewidth = plot_reference_line_width,
             colour = plot_text_color) +
  geom_point(data = combined_binomial_table,
             aes(x = tdiv, y = proportion.k2, colour = mig.tag, group = mig.tag),
             size = plot_point_size,
             alpha = plot_point_alpha) +
  geom_line(data = fitted_prediction_table,
            aes(x = tdiv, y = fitted.proportion.k2, colour = mig.tag, group = mig.tag),
            linewidth = plot_fitted_line_width) +
  geom_vline(data = threshold_table,
             aes(xintercept = tdiv.at.threshold.k2, colour = mig.tag),
             linewidth = plot_threshold_line_width,
             linetype = plot_threshold_line_type,
             show.legend = FALSE) +
  facet_grid(rows = vars(method),
             cols = vars(design),
             labeller = labeller(design = c("Balanced" = "Even sampling", "Unbalanced" = "Uneven sampling")),
             axes = "all",
             axis.labels = "margins") +
  scale_color_manual(values = migration_colors,
                     breaks = migration_labels,
                     labels = migration_labels,
                     drop = FALSE) +
  scale_x_continuous(expand = expansion(add = c(plot_x_axis_left_expansion, plot_x_axis_right_expansion))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.25),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "Divergence time (generations)",
       y = "Proportion choosing k = 2",
       colour = "Migration rate") +
  theme_classic(base_size = plot_base_size,
                base_family = plot_font_family) +
  theme(text = element_text(family = plot_font_family,
                            colour = plot_text_color),
        axis.title.x = element_text(family = plot_font_family,
                                    size = plot_axis_title_size,
                                    colour = plot_text_color,
                                    margin = margin(t = plot_axis_title_margin_mm, unit = "mm")),
        axis.title.y = element_text(family = plot_font_family,
                                    size = plot_axis_title_size,
                                    colour = plot_text_color,
                                    margin = margin(r = plot_axis_title_margin_mm, unit = "mm")),
        axis.text = element_text(family = plot_font_family,
                                 size = plot_axis_text_size,
                                 colour = plot_text_color),
        strip.background = element_blank(),
        strip.text.x = element_text(family = plot_font_family,
                                    size = plot_column_title_size,
                                    face = "bold",
                                    colour = plot_text_color,
                                    margin = margin(b = plot_column_title_margin_mm, unit = "mm")),
        strip.text.y = element_text(family = plot_font_family,
                                    size = plot_row_title_size,
                                    face = "bold",
                                    colour = plot_text_color),
        legend.position = "bottom",
        legend.title = element_text(family = plot_font_family,
                                    size = plot_legend_title_size,
                                    face = "bold",
                                    colour = plot_text_color),
        legend.text = element_text(family = plot_font_family,
                                   size = plot_legend_text_size,
                                   colour = plot_text_color),
        legend.box.background = element_rect(colour = plot_text_color,
                                             fill = NA,
                                             linewidth = legend_box_line_width),
        legend.box.margin = margin(legend_box_margin_mm,
                                   legend_box_margin_mm,
                                   legend_box_margin_mm,
                                   legend_box_margin_mm,
                                   unit = "mm"),
        panel.spacing = grid::unit(panel_spacing_mm, "mm"),
        plot.margin = margin(1, 1, 1, 1, unit = "mm")) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))


## Show and save figure
combined_plot
ggsave(filename = combined_svg,
       plot = combined_plot,
       device = svglite::svglite,
       width = plot_width_cm,
       height = plot_height_cm,
       units = "cm",
       bg = "white",
       limitsize = FALSE)
message("Combined SVG saved to: ", combined_svg)




#### Supplementary Figure S6 ##################################################################


## Set environment
rm(list = ls()) #clear environment
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
invisible(gc())


## Install and load required packages
required_packages <- c("svglite", "tidyverse", "viridisLite")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
  library(package_name, character.only = TRUE)
}


## Set main plotting parameters
plot_width_cm <- 15.36
plot_height_cm <- 19.24

plot_font_family <- "Arial"
plot_text_color <- "black"

plot_base_size <- 9
plot_axis_title_size <- 9
plot_axis_text_size <- 7.1
plot_column_title_size <- 9
plot_row_title_size <- 9
plot_legend_title_size <- 9
plot_legend_text_size <- 9

plot_axis_title_margin_mm <- 3
legend_box_margin_mm <- 1
panel_spacing_mm <- 2.5
plot_column_title_margin_mm <- 2.5
plot_x_axis_left_expansion <- 10000
plot_x_axis_right_expansion <- 10000

plot_point_size <- 1
plot_point_alpha <- 0.7
plot_fitted_line_width <- 1
plot_threshold_line_width <- 0.6
plot_threshold_line_type <- "dashed"
plot_K2_threshold_proportion <- 0.5
plot_reference_line_positions <- c(50000, 100000, 150000, 200000)
plot_reference_line_width <- 0.1
legend_box_line_width <- 0.3
migration_values <- c(0, 1e-6, 4e-6, 7e-6)
migration_labels <- c("0", "1e-6", "4e-6", "7e-6")
migration_colors <- stats::setNames(c("#2B0A3D", "#31668EFF", "#30B57BFF", "#F0DD1CFF"), migration_labels)


## Set input and output
simulation_dir <- file.path("Simulations", "Simulation_set_2")
results_root_dir <- file.path(simulation_dir, "fastsimcoal2_results")
balanced_results_dir <- file.path(results_root_dir, "symmetric_8")
unbalanced_results_dir <- file.path(results_root_dir, "asymmetric_8")
balanced_SOM_results_csv <- file.path(balanced_results_dir, "fastsimcoal2_SOM_results.csv")
balanced_SOM_optim_k_csv <- file.path(balanced_results_dir, "fastsimcoal2_SOM_optim_k_summary.csv")
balanced_STRUCTURE_csv <- file.path(balanced_results_dir, "STRUCTURE_results", "STRUCTURE_binomial_table_by_vcf.csv")
unbalanced_SOM_results_csv <- file.path(unbalanced_results_dir, "fastsimcoal2_SOM_results.csv")
unbalanced_SOM_optim_k_csv <- file.path(unbalanced_results_dir, "fastsimcoal2_SOM_optim_k_summary.csv")
unbalanced_STRUCTURE_csv <- file.path(unbalanced_results_dir, "STRUCTURE_results", "STRUCTURE_binomial_table_by_vcf.csv")
plot_dir <- "Figure_files"
plot_file_name <- "Supplementary_figure_S6.svg"
combined_svg <- file.path(plot_dir, plot_file_name)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


## Check that all completed result files exist
required_result_files <- c(balanced_SOM_results_csv,
                           balanced_SOM_optim_k_csv,
                           balanced_STRUCTURE_csv,
                           unbalanced_SOM_results_csv,
                           unbalanced_SOM_optim_k_csv,
                           unbalanced_STRUCTURE_csv)
missing_result_files <- required_result_files[!file.exists(required_result_files)]
if (length(missing_result_files) > 0) stop("The following completed result files are missing: ", paste(missing_result_files, collapse = ", "))


## Create additional functions
check.required.columns <- function(input_table,
                                   required_columns,
                                   table_name) {
  missing_columns <- setdiff(required_columns, colnames(input_table))
  if (length(missing_columns) > 0) stop(table_name, " is missing the following columns: ", paste(missing_columns, collapse = ", "))
}
standardize.migration.labels <- function(input_data) {
  migration_indices <- vapply(input_data$mig, function(current_migration) {
    migration_difference <- abs(migration_values - current_migration)
    matching_index <- which.min(migration_difference)
    if (length(matching_index) != 1 || !is.finite(migration_difference[matching_index]) || migration_difference[matching_index] > sqrt(.Machine$double.eps)) return(NA_integer_)
    matching_index
  }, integer(1))
  if (any(is.na(migration_indices))) stop("At least one migration rate could not be matched to the expected migration rates")
  input_data$mig.tag <- factor(migration_labels[migration_indices], levels = migration_labels)
  input_data
}
create.analysis.binomial.table <- function(SOM_results_csv,
                                           SOM_optim_k_csv,
                                           STRUCTURE_csv,
                                           design_label) {
  result_table <- read.csv(SOM_results_csv, stringsAsFactors = FALSE)
  optim_k_result_table <- read.csv(SOM_optim_k_csv, stringsAsFactors = FALSE)
  STRUCTURE_binomial_table <- read.csv(STRUCTURE_csv, stringsAsFactors = FALSE)
  check.required.columns(result_table, c("file", "status", "mig", "tdiv", "deNovo.kmeans.best.k", "sNMF.best.k"), paste0(design_label, " SOM result table"))
  check.required.columns(optim_k_result_table, c("file", "mig", "tdiv", "Count", "k.label"), paste0(design_label, " SOM optim-k table"))
  check.required.columns(STRUCTURE_binomial_table, c("file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2"), paste0(design_label, " STRUCTURE binomial table"))
  plot_result_table <- result_table[result_table$status == "ok", ]
  plot_result_table <- plot_result_table[order(plot_result_table$mig, plot_result_table$tdiv), ]
  if (nrow(plot_result_table) == 0) stop("No successful SOM result rows are available for ", design_label)
  SOM_optim_k_count_table <- optim_k_result_table[!is.na(optim_k_result_table$Count) & optim_k_result_table$file %in% plot_result_table$file, ]
  SOM_optim_k_count_table$Count <- as.numeric(as.character(SOM_optim_k_count_table$Count))
  SOM_total_count_table <- stats::aggregate(Count ~ file + mig + tdiv, data = SOM_optim_k_count_table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE))
  colnames(SOM_total_count_table)[colnames(SOM_total_count_table) == "Count"] <- "n.total"
  SOM_k2_count_table <- SOM_optim_k_count_table[SOM_optim_k_count_table$k.label == "k2", ]
  if (nrow(SOM_k2_count_table) > 0) {
    SOM_k2_count_table <- stats::aggregate(Count ~ file + mig + tdiv, data = SOM_k2_count_table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE))
    colnames(SOM_k2_count_table)[colnames(SOM_k2_count_table) == "Count"] <- "n.k2"
  } else {
    SOM_k2_count_table <- SOM_total_count_table[, c("file", "mig", "tdiv")]
    SOM_k2_count_table$n.k2 <- 0
  }
  SOM_binomial_table <- merge(SOM_total_count_table, SOM_k2_count_table, by = c("file", "mig", "tdiv"), all.x = TRUE)
  SOM_binomial_table$n.k2[is.na(SOM_binomial_table$n.k2)] <- 0
  SOM_binomial_table$n.not.k2 <- SOM_binomial_table$n.total - SOM_binomial_table$n.k2
  SOM_binomial_table$proportion.k2 <- SOM_binomial_table$n.k2 / SOM_binomial_table$n.total
  SOM_binomial_table$method <- "SOM"
  SOM_binomial_table <- SOM_binomial_table[, c("method", "file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")]
  DAPC_binomial_source_table <- result_table[result_table$status == "ok" & !is.na(result_table$deNovo.kmeans.best.k), ]
  if ("deNovo.kmeans.status" %in% colnames(DAPC_binomial_source_table)) DAPC_binomial_source_table <- DAPC_binomial_source_table[DAPC_binomial_source_table$deNovo.kmeans.status == "ok", ]
  DAPC_best_k <- as.integer(as.character(DAPC_binomial_source_table$deNovo.kmeans.best.k))
  DAPC_binomial_table <- data.frame(method = "DAPC",
                                    file = DAPC_binomial_source_table$file,
                                    mig = DAPC_binomial_source_table$mig,
                                    tdiv = DAPC_binomial_source_table$tdiv,
                                    n.k2 = as.integer(DAPC_best_k == 2L),
                                    n.not.k2 = as.integer(DAPC_best_k != 2L),
                                    proportion.k2 = as.numeric(DAPC_best_k == 2L),
                                    stringsAsFactors = FALSE)
  sNMF_binomial_source_table <- result_table[result_table$status == "ok" & !is.na(result_table$sNMF.best.k), ]
  if ("sNMF.status" %in% colnames(sNMF_binomial_source_table)) sNMF_binomial_source_table <- sNMF_binomial_source_table[sNMF_binomial_source_table$sNMF.status == "ok", ]
  sNMF_best_k <- as.integer(as.character(sNMF_binomial_source_table$sNMF.best.k))
  sNMF_binomial_table <- data.frame(method = "sNMF",
                                    file = sNMF_binomial_source_table$file,
                                    mig = sNMF_binomial_source_table$mig,
                                    tdiv = sNMF_binomial_source_table$tdiv,
                                    n.k2 = as.integer(sNMF_best_k == 2L),
                                    n.not.k2 = as.integer(sNMF_best_k != 2L),
                                    proportion.k2 = as.numeric(sNMF_best_k == 2L),
                                    stringsAsFactors = FALSE)
  STRUCTURE_binomial_table <- STRUCTURE_binomial_table[, c("file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")]
  STRUCTURE_binomial_table$method <- "STRUCTURE"
  STRUCTURE_binomial_table <- STRUCTURE_binomial_table[, c("method", "file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")]
  combined_binomial_table <- rbind(SOM_binomial_table, DAPC_binomial_table, sNMF_binomial_table, STRUCTURE_binomial_table)
  combined_binomial_table$design <- design_label
  combined_binomial_table <- standardize.migration.labels(combined_binomial_table)
  combined_binomial_table
}
fit.K2.binomial.models <- function(input_data) {
  prediction_list <- list()
  threshold_list <- list()
  prediction_index <- 1
  threshold_index <- 1
  for (current_design in levels(input_data$design)) {
    for (current_method in levels(input_data$method)) {
      for (current_migration in levels(input_data$mig.tag)) {
        current_binomial_table <- input_data[input_data$design == current_design & input_data$method == current_method & input_data$mig.tag == current_migration, ]
        current_binomial_table <- current_binomial_table[is.finite(current_binomial_table$tdiv), ]
        if (nrow(current_binomial_table) == 0) next
        current_binomial_table <- current_binomial_table[order(current_binomial_table$tdiv), ]
        prediction_tdiv <- seq(min(current_binomial_table$tdiv), max(current_binomial_table$tdiv), length.out = 1000)
        if (all(current_binomial_table$n.k2 == 0)) {
          current_prediction_table <- data.frame(design = current_design,
                                                 method = current_method,
                                                 mig.tag = current_migration,
                                                 tdiv = prediction_tdiv,
                                                 fitted.proportion.k2 = 0,
                                                 stringsAsFactors = FALSE)
        } else if (all(current_binomial_table$n.not.k2 == 0)) {
          current_prediction_table <- data.frame(design = current_design,
                                                 method = current_method,
                                                 mig.tag = current_migration,
                                                 tdiv = prediction_tdiv,
                                                 fitted.proportion.k2 = 1,
                                                 stringsAsFactors = FALSE)
        } else {
          current_fit <- stats::glm(cbind(n.k2, n.not.k2) ~ tdiv, data = current_binomial_table, family = stats::binomial(link = "logit"), control = stats::glm.control(maxit = 100))
          current_prediction_table <- data.frame(design = current_design,
                                                 method = current_method,
                                                 mig.tag = current_migration,
                                                 tdiv = prediction_tdiv,
                                                 fitted.proportion.k2 = stats::predict(current_fit, newdata = data.frame(tdiv = prediction_tdiv), type = "response"),
                                                 stringsAsFactors = FALSE)
          current_coefficients <- stats::coef(current_fit)
          current_intercept <- unname(current_coefficients["(Intercept)"])
          current_tdiv_slope <- unname(current_coefficients["tdiv"])
          if (is.finite(current_intercept) && is.finite(current_tdiv_slope) && current_tdiv_slope != 0) {
            current_threshold_tdiv <- (stats::qlogis(plot_K2_threshold_proportion) - current_intercept) / current_tdiv_slope
            threshold_within_range <- current_threshold_tdiv >= min(current_binomial_table$tdiv) && current_threshold_tdiv <= max(current_binomial_table$tdiv)
            if (threshold_within_range) {
              threshold_list[[threshold_index]] <- data.frame(design = current_design,
                                                              method = current_method,
                                                              mig.tag = current_migration,
                                                              tdiv.at.threshold.k2 = current_threshold_tdiv,
                                                              stringsAsFactors = FALSE)
              threshold_index <- threshold_index + 1
            }
          }
        }
        prediction_list[[prediction_index]] <- current_prediction_table
        prediction_index <- prediction_index + 1
      }
    }
  }
  prediction_table <- do.call(rbind, prediction_list)
  if (length(threshold_list) > 0) {
    threshold_table <- do.call(rbind, threshold_list)
  } else {
    threshold_table <- data.frame(design = character(0), method = character(0), mig.tag = character(0), tdiv.at.threshold.k2 = numeric(0))
  }
  prediction_table$design <- factor(prediction_table$design, levels = levels(input_data$design))
  prediction_table$method <- factor(prediction_table$method, levels = levels(input_data$method))
  prediction_table$mig.tag <- factor(prediction_table$mig.tag, levels = migration_labels)
  threshold_table$design <- factor(threshold_table$design, levels = levels(input_data$design))
  threshold_table$method <- factor(threshold_table$method, levels = levels(input_data$method))
  threshold_table$mig.tag <- factor(threshold_table$mig.tag, levels = migration_labels)
  list(prediction_table = prediction_table, threshold_table = threshold_table)
}


## Load balanced analysis
balanced_binomial_table <- create.analysis.binomial.table(SOM_results_csv = balanced_SOM_results_csv,
                                                          SOM_optim_k_csv = balanced_SOM_optim_k_csv,
                                                          STRUCTURE_csv = balanced_STRUCTURE_csv,
                                                          design_label = "Balanced")


## Load unbalanced analysis
unbalanced_binomial_table <- create.analysis.binomial.table(SOM_results_csv = unbalanced_SOM_results_csv,
                                                            SOM_optim_k_csv = unbalanced_SOM_optim_k_csv,
                                                            STRUCTURE_csv = unbalanced_STRUCTURE_csv,
                                                            design_label = "Unbalanced")


## Combine balanced and unbalanced analyses
combined_binomial_table <- rbind(balanced_binomial_table, unbalanced_binomial_table)
combined_binomial_table$design <- factor(combined_binomial_table$design, levels = c("Balanced", "Unbalanced"))
combined_binomial_table$method <- factor(combined_binomial_table$method, levels = c("SOM", "DAPC", "sNMF", "STRUCTURE"))


## Fit binomial models
fitted_K2_output <- fit.K2.binomial.models(combined_binomial_table)
fitted_prediction_table <- fitted_K2_output$prediction_table
threshold_table <- fitted_K2_output$threshold_table


## Create figure
combined_plot <- ggplot() +
  geom_vline(xintercept = plot_reference_line_positions,
             linewidth = plot_reference_line_width,
             colour = plot_text_color) +
  geom_point(data = combined_binomial_table,
             aes(x = tdiv, y = proportion.k2, colour = mig.tag, group = mig.tag),
             size = plot_point_size,
             alpha = plot_point_alpha) +
  geom_line(data = fitted_prediction_table,
            aes(x = tdiv, y = fitted.proportion.k2, colour = mig.tag, group = mig.tag),
            linewidth = plot_fitted_line_width) +
  geom_vline(data = threshold_table,
             aes(xintercept = tdiv.at.threshold.k2, colour = mig.tag),
             linewidth = plot_threshold_line_width,
             linetype = plot_threshold_line_type,
             show.legend = FALSE) +
  facet_grid(rows = vars(method),
             cols = vars(design),
             labeller = labeller(design = c("Balanced" = "Even sampling", "Unbalanced" = "Uneven sampling")),
             axes = "all",
             axis.labels = "margins") +
  scale_color_manual(values = migration_colors,
                     breaks = migration_labels,
                     labels = migration_labels,
                     drop = FALSE) +
  scale_x_continuous(expand = expansion(add = c(plot_x_axis_left_expansion, plot_x_axis_right_expansion))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.25),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "Divergence time (generations)",
       y = "Proportion choosing k = 2",
       colour = "Migration rate") +
  theme_classic(base_size = plot_base_size,
                base_family = plot_font_family) +
  theme(text = element_text(family = plot_font_family,
                            colour = plot_text_color),
        axis.title.x = element_text(family = plot_font_family,
                                    size = plot_axis_title_size,
                                    colour = plot_text_color,
                                    margin = margin(t = plot_axis_title_margin_mm, unit = "mm")),
        axis.title.y = element_text(family = plot_font_family,
                                    size = plot_axis_title_size,
                                    colour = plot_text_color,
                                    margin = margin(r = plot_axis_title_margin_mm, unit = "mm")),
        axis.text = element_text(family = plot_font_family,
                                 size = plot_axis_text_size,
                                 colour = plot_text_color),
        strip.background = element_blank(),
        strip.text.x = element_text(family = plot_font_family,
                                    size = plot_column_title_size,
                                    face = "bold",
                                    colour = plot_text_color,
                                    margin = margin(b = plot_column_title_margin_mm, unit = "mm")),
        strip.text.y = element_text(family = plot_font_family,
                                    size = plot_row_title_size,
                                    face = "bold",
                                    colour = plot_text_color),
        legend.position = "bottom",
        legend.title = element_text(family = plot_font_family,
                                    size = plot_legend_title_size,
                                    face = "bold",
                                    colour = plot_text_color),
        legend.text = element_text(family = plot_font_family,
                                   size = plot_legend_text_size,
                                   colour = plot_text_color),
        legend.box.background = element_rect(colour = plot_text_color,
                                             fill = NA,
                                             linewidth = legend_box_line_width),
        legend.box.margin = margin(legend_box_margin_mm,
                                   legend_box_margin_mm,
                                   legend_box_margin_mm,
                                   legend_box_margin_mm,
                                   unit = "mm"),
        panel.spacing = grid::unit(panel_spacing_mm, "mm"),
        plot.margin = margin(1, 1, 1, 1, unit = "mm")) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))


## Show and save figure
combined_plot
ggsave(filename = combined_svg,
       plot = combined_plot,
       device = svglite::svglite,
       width = plot_width_cm,
       height = plot_height_cm,
       units = "cm",
       bg = "white",
       limitsize = FALSE)
message("Combined SVG saved to: ", combined_svg)




#### Supplementary Figure S7 ##################################################################

## Set environment
rm(list = ls()) #clear environment
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
invisible(gc())


## Install and load required packages
required_packages <- c("svglite", "tidyverse", "viridisLite")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
  library(package_name, character.only = TRUE)
}


## Set main plotting parameters
plot_width_cm <- 15.36
plot_height_cm <- 19.24

plot_font_family <- "Arial"
plot_text_color <- "black"

plot_base_size <- 9
plot_axis_title_size <- 9
plot_axis_text_size <- 7.1
plot_column_title_size <- 9
plot_row_title_size <- 9
plot_legend_title_size <- 9
plot_legend_text_size <- 9

plot_axis_title_margin_mm <- 3
legend_box_margin_mm <- 1
panel_spacing_mm <- 2.5
plot_column_title_margin_mm <- 2.5
plot_x_axis_left_expansion <- 10000
plot_x_axis_right_expansion <- 10000

plot_point_size <- 1
plot_point_alpha <- 0.7
plot_fitted_line_width <- 1
plot_threshold_line_width <- 0.6
plot_threshold_line_type <- "dashed"
plot_K2_threshold_proportion <- 0.5
plot_reference_line_positions <- c(50000, 100000, 150000, 200000)
plot_reference_line_width <- 0.1
legend_box_line_width <- 0.3
migration_values <- c(0, 1e-6, 4e-6, 7e-6)
migration_labels <- c("0", "1e-6", "4e-6", "7e-6")
migration_colors <- stats::setNames(c("#2B0A3D", "#31668EFF", "#30B57BFF", "#F0DD1CFF"), migration_labels)


## Set input and output
simulation_dir <- file.path("Simulations", "Simulation_set_2")
results_root_dir <- file.path(simulation_dir, "fastsimcoal2_results")
balanced_results_dir <- file.path(results_root_dir, "symmetric_24")
unbalanced_results_dir <- file.path(results_root_dir, "asymmetric_24")
balanced_SOM_results_csv <- file.path(balanced_results_dir, "fastsimcoal2_SOM_results.csv")
balanced_SOM_optim_k_csv <- file.path(balanced_results_dir, "fastsimcoal2_SOM_optim_k_summary.csv")
balanced_STRUCTURE_csv <- file.path(balanced_results_dir, "STRUCTURE_results", "STRUCTURE_binomial_table_by_vcf.csv")
unbalanced_SOM_results_csv <- file.path(unbalanced_results_dir, "fastsimcoal2_SOM_results.csv")
unbalanced_SOM_optim_k_csv <- file.path(unbalanced_results_dir, "fastsimcoal2_SOM_optim_k_summary.csv")
unbalanced_STRUCTURE_csv <- file.path(unbalanced_results_dir, "STRUCTURE_results", "STRUCTURE_binomial_table_by_vcf.csv")
plot_dir <- "Figure_files"
plot_file_name <- "Supplementary_figure_S7.svg"
combined_svg <- file.path(plot_dir, plot_file_name)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


## Check that all completed result files exist
required_result_files <- c(balanced_SOM_results_csv,
                           balanced_SOM_optim_k_csv,
                           balanced_STRUCTURE_csv,
                           unbalanced_SOM_results_csv,
                           unbalanced_SOM_optim_k_csv,
                           unbalanced_STRUCTURE_csv)
missing_result_files <- required_result_files[!file.exists(required_result_files)]
if (length(missing_result_files) > 0) stop("The following completed result files are missing: ", paste(missing_result_files, collapse = ", "))


## Create additional functions
check.required.columns <- function(input_table,
                                   required_columns,
                                   table_name) {
  missing_columns <- setdiff(required_columns, colnames(input_table))
  if (length(missing_columns) > 0) stop(table_name, " is missing the following columns: ", paste(missing_columns, collapse = ", "))
}
standardize.migration.labels <- function(input_data) {
  migration_indices <- vapply(input_data$mig, function(current_migration) {
    migration_difference <- abs(migration_values - current_migration)
    matching_index <- which.min(migration_difference)
    if (length(matching_index) != 1 ||
        !is.finite(migration_difference[matching_index]) ||
        migration_difference[matching_index] > sqrt(.Machine$double.eps)) return(NA_integer_)
    matching_index
  }, integer(1))
  if (any(is.na(migration_indices))) stop("At least one migration rate could not be matched to the expected migration rates")
  input_data$mig.tag <- factor(migration_labels[migration_indices], levels = migration_labels)
  input_data
}
create.analysis.binomial.table <- function(SOM_results_csv,
                                           SOM_optim_k_csv,
                                           STRUCTURE_csv,
                                           design_label) {
  result_table <- read.csv(SOM_results_csv, stringsAsFactors = FALSE)
  optim_k_result_table <- read.csv(SOM_optim_k_csv, stringsAsFactors = FALSE)
  STRUCTURE_binomial_table <- read.csv(STRUCTURE_csv, stringsAsFactors = FALSE)
  check.required.columns(result_table,
                         c("file", "status", "mig", "tdiv", "fst", "deNovo.kmeans.best.k", "sNMF.best.k"),
                         paste0(design_label, " SOM result table"))
  check.required.columns(optim_k_result_table,
                         c("file", "mig", "tdiv", "Count", "k.label"),
                         paste0(design_label, " SOM optim-k table"))
  check.required.columns(STRUCTURE_binomial_table,
                         c("file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2"),
                         paste0(design_label, " STRUCTURE binomial table"))
  plot_result_table <- result_table[result_table$status == "ok", ]
  plot_result_table <- plot_result_table[order(plot_result_table$mig, plot_result_table$tdiv), ]
  if (nrow(plot_result_table) == 0) stop("No successful SOM result rows are available for ", design_label)
  SOM_optim_k_count_table <- optim_k_result_table[!is.na(optim_k_result_table$Count) & optim_k_result_table$file %in% plot_result_table$file, ]
  SOM_optim_k_count_table$Count <- as.numeric(as.character(SOM_optim_k_count_table$Count))
  SOM_total_count_table <- stats::aggregate(Count ~ file + mig + tdiv,
                                            data = SOM_optim_k_count_table,
                                            FUN = function(x) sum(as.numeric(x), na.rm = TRUE))
  colnames(SOM_total_count_table)[colnames(SOM_total_count_table) == "Count"] <- "n.total"
  SOM_k2_count_table <- SOM_optim_k_count_table[SOM_optim_k_count_table$k.label == "k2", ]
  if (nrow(SOM_k2_count_table) > 0) {
    SOM_k2_count_table <- stats::aggregate(Count ~ file + mig + tdiv,
                                           data = SOM_k2_count_table,
                                           FUN = function(x) sum(as.numeric(x), na.rm = TRUE))
    colnames(SOM_k2_count_table)[colnames(SOM_k2_count_table) == "Count"] <- "n.k2"
  } else {
    SOM_k2_count_table <- SOM_total_count_table[, c("file", "mig", "tdiv")]
    SOM_k2_count_table$n.k2 <- 0
  }
  SOM_binomial_table <- merge(SOM_total_count_table,
                              SOM_k2_count_table,
                              by = c("file", "mig", "tdiv"),
                              all.x = TRUE)
  SOM_binomial_table$n.k2[is.na(SOM_binomial_table$n.k2)] <- 0
  SOM_binomial_table$n.not.k2 <- SOM_binomial_table$n.total - SOM_binomial_table$n.k2
  SOM_binomial_table$proportion.k2 <- SOM_binomial_table$n.k2 / SOM_binomial_table$n.total
  SOM_binomial_table$method <- "SOM"
  SOM_binomial_table <- SOM_binomial_table[, c("method", "file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")]
  DAPC_binomial_source_table <- result_table[result_table$status == "ok" & !is.na(result_table$deNovo.kmeans.best.k), ]
  if ("deNovo.kmeans.status" %in% colnames(DAPC_binomial_source_table)) DAPC_binomial_source_table <- DAPC_binomial_source_table[DAPC_binomial_source_table$deNovo.kmeans.status == "ok", ]
  DAPC_best_k <- as.integer(as.character(DAPC_binomial_source_table$deNovo.kmeans.best.k))
  DAPC_binomial_table <- data.frame(method = "DAPC",
                                    file = DAPC_binomial_source_table$file,
                                    mig = DAPC_binomial_source_table$mig,
                                    tdiv = DAPC_binomial_source_table$tdiv,
                                    n.k2 = as.integer(DAPC_best_k == 2L),
                                    n.not.k2 = as.integer(DAPC_best_k != 2L),
                                    proportion.k2 = as.numeric(DAPC_best_k == 2L),
                                    stringsAsFactors = FALSE)
  sNMF_binomial_source_table <- result_table[result_table$status == "ok" & !is.na(result_table$sNMF.best.k), ]
  if ("sNMF.status" %in% colnames(sNMF_binomial_source_table)) sNMF_binomial_source_table <- sNMF_binomial_source_table[sNMF_binomial_source_table$sNMF.status == "ok", ]
  sNMF_best_k <- as.integer(as.character(sNMF_binomial_source_table$sNMF.best.k))
  sNMF_binomial_table <- data.frame(method = "sNMF",
                                    file = sNMF_binomial_source_table$file,
                                    mig = sNMF_binomial_source_table$mig,
                                    tdiv = sNMF_binomial_source_table$tdiv,
                                    n.k2 = as.integer(sNMF_best_k == 2L),
                                    n.not.k2 = as.integer(sNMF_best_k != 2L),
                                    proportion.k2 = as.numeric(sNMF_best_k == 2L),
                                    stringsAsFactors = FALSE)
  STRUCTURE_binomial_table <- STRUCTURE_binomial_table[, c("file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")]
  STRUCTURE_binomial_table$method <- "STRUCTURE"
  STRUCTURE_binomial_table <- STRUCTURE_binomial_table[, c("method", "file", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")]
  combined_binomial_table <- rbind(SOM_binomial_table,
                                   DAPC_binomial_table,
                                   sNMF_binomial_table,
                                   STRUCTURE_binomial_table)
  combined_binomial_table$design <- design_label
  combined_binomial_table <- standardize.migration.labels(combined_binomial_table)
  list(result_table = result_table,
       combined_binomial_table = combined_binomial_table)
}
calculate.K2.thresholds <- function(input_data) {
  threshold_list <- list()
  threshold_index <- 1
  for (current_design in levels(input_data$design)) {
    for (current_method in levels(input_data$method)) {
      for (current_migration in levels(input_data$mig.tag)) {
        current_binomial_table <- input_data[input_data$design == current_design &
                                               input_data$method == current_method &
                                               input_data$mig.tag == current_migration, ]
        current_binomial_table <- current_binomial_table[is.finite(current_binomial_table$tdiv), ]
        if (nrow(current_binomial_table) == 0 ||
            all(current_binomial_table$n.k2 == 0) ||
            all(current_binomial_table$n.not.k2 == 0)) next
        current_fit <- stats::glm(cbind(n.k2, n.not.k2) ~ tdiv,
                                  data = current_binomial_table,
                                  family = stats::binomial(link = "logit"),
                                  control = stats::glm.control(maxit = 100))
        current_coefficients <- stats::coef(current_fit)
        current_intercept <- unname(current_coefficients["(Intercept)"])
        current_tdiv_slope <- unname(current_coefficients["tdiv"])
        if (!is.finite(current_intercept) || !is.finite(current_tdiv_slope) || current_tdiv_slope == 0) next
        current_threshold_tdiv <- (stats::qlogis(plot_K2_threshold_proportion) - current_intercept) / current_tdiv_slope
        threshold_within_range <- current_threshold_tdiv >= min(current_binomial_table$tdiv) &&
          current_threshold_tdiv <= max(current_binomial_table$tdiv)
        if (!threshold_within_range) next
        threshold_list[[threshold_index]] <- data.frame(design = current_design,
                                                        method = current_method,
                                                        mig.tag = current_migration,
                                                        tdiv.at.threshold.k2 = current_threshold_tdiv,
                                                        stringsAsFactors = FALSE)
        threshold_index <- threshold_index + 1
      }
    }
  }
  if (length(threshold_list) > 0) {
    threshold_table <- do.call(rbind, threshold_list)
  } else {
    threshold_table <- data.frame(design = character(0),
                                  method = character(0),
                                  mig.tag = character(0),
                                  tdiv.at.threshold.k2 = numeric(0))
  }
  threshold_table$design <- factor(threshold_table$design, levels = levels(input_data$design))
  threshold_table$method <- factor(threshold_table$method, levels = levels(input_data$method))
  threshold_table$mig.tag <- factor(threshold_table$mig.tag, levels = migration_labels)
  threshold_table
}
create.Fst.plot.table <- function(result_table,
                                  design_label) {
  current_Fst_table <- result_table[result_table$status == "ok" & is.finite(result_table$fst), ]
  current_Fst_table <- current_Fst_table[order(current_Fst_table$mig, current_Fst_table$tdiv), ]
  if (nrow(current_Fst_table) == 0) stop("No rows with finite Fst values are available for ", design_label)
  current_Fst_table <- standardize.migration.labels(current_Fst_table)
  current_Fst_table$design <- design_label
  method_levels <- c("SOM", "DAPC", "sNMF", "STRUCTURE")
  method_Fst_tables <- lapply(method_levels, function(current_method) {
    current_method_table <- current_Fst_table
    current_method_table$method <- current_method
    current_method_table
  })
  do.call(rbind, method_Fst_tables)
}


## Load balanced analysis
balanced_analysis_output <- create.analysis.binomial.table(SOM_results_csv = balanced_SOM_results_csv,
                                                           SOM_optim_k_csv = balanced_SOM_optim_k_csv,
                                                           STRUCTURE_csv = balanced_STRUCTURE_csv,
                                                           design_label = "Balanced")
balanced_result_table <- balanced_analysis_output$result_table
balanced_binomial_table <- balanced_analysis_output$combined_binomial_table


## Load unbalanced analysis
unbalanced_analysis_output <- create.analysis.binomial.table(SOM_results_csv = unbalanced_SOM_results_csv,
                                                             SOM_optim_k_csv = unbalanced_SOM_optim_k_csv,
                                                             STRUCTURE_csv = unbalanced_STRUCTURE_csv,
                                                             design_label = "Unbalanced")
unbalanced_result_table <- unbalanced_analysis_output$result_table
unbalanced_binomial_table <- unbalanced_analysis_output$combined_binomial_table


## Combine balanced and unbalanced analyses
combined_binomial_table <- rbind(balanced_binomial_table,
                                 unbalanced_binomial_table)
combined_binomial_table$design <- factor(combined_binomial_table$design,
                                         levels = c("Balanced", "Unbalanced"))
combined_binomial_table$method <- factor(combined_binomial_table$method,
                                         levels = c("SOM", "DAPC", "sNMF", "STRUCTURE"))


## Prepare Fst plot table
balanced_Fst_table <- create.Fst.plot.table(result_table = balanced_result_table,
                                            design_label = "Balanced")
unbalanced_Fst_table <- create.Fst.plot.table(result_table = unbalanced_result_table,
                                              design_label = "Unbalanced")
Fst_plot_table <- rbind(balanced_Fst_table,
                        unbalanced_Fst_table)
Fst_plot_table$design <- factor(Fst_plot_table$design,
                                levels = c("Balanced", "Unbalanced"))
Fst_plot_table$method <- factor(Fst_plot_table$method,
                                levels = c("SOM", "DAPC", "sNMF", "STRUCTURE"))


## Calculate K2 thresholds
threshold_table <- calculate.K2.thresholds(combined_binomial_table)


## Create figure
combined_plot <- ggplot() +
  geom_vline(xintercept = plot_reference_line_positions,
             linewidth = plot_reference_line_width,
             colour = plot_text_color) +
  geom_point(data = Fst_plot_table,
             aes(x = tdiv,
                 y = fst,
                 colour = mig.tag,
                 group = mig.tag),
             size = plot_point_size,
             alpha = plot_point_alpha) +
  geom_smooth(data = Fst_plot_table,
              aes(x = tdiv,
                  y = fst,
                  colour = mig.tag,
                  group = mig.tag),
              method = "loess",
              formula = y ~ x,
              linewidth = plot_fitted_line_width,
              se = FALSE) +
  geom_vline(data = threshold_table,
             aes(xintercept = tdiv.at.threshold.k2,
                 colour = mig.tag),
             linewidth = plot_threshold_line_width,
             linetype = plot_threshold_line_type,
             show.legend = FALSE) +
  facet_grid(rows = vars(method),
             cols = vars(design),
             labeller = labeller(design = c("Balanced" = "Balanced sampling",
                                            "Unbalanced" = "Unbalanced sampling")),
             axes = "all",
             axis.labels = "margins") +
  scale_color_manual(values = migration_colors,
                     breaks = migration_labels,
                     labels = migration_labels,
                     drop = FALSE) +
  scale_x_continuous(expand = expansion(add = c(plot_x_axis_left_expansion,
                                                plot_x_axis_right_expansion))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "Divergence time (generations)",
       y = "Weir and Cockerham Fst",
       colour = "Migration rate") +
  theme_classic(base_size = plot_base_size,
                base_family = plot_font_family) +
  theme(text = element_text(family = plot_font_family,
                            colour = plot_text_color),
        axis.title.x = element_text(family = plot_font_family,
                                    size = plot_axis_title_size,
                                    colour = plot_text_color,
                                    margin = margin(t = plot_axis_title_margin_mm, unit = "mm")),
        axis.title.y = element_text(family = plot_font_family,
                                    size = plot_axis_title_size,
                                    colour = plot_text_color,
                                    margin = margin(r = plot_axis_title_margin_mm, unit = "mm")),
        axis.text = element_text(family = plot_font_family,
                                 size = plot_axis_text_size,
                                 colour = plot_text_color),
        strip.background = element_blank(),
        strip.text.x = element_text(family = plot_font_family,
                                    size = plot_column_title_size,
                                    face = "bold",
                                    colour = plot_text_color,
                                    margin = margin(b = plot_column_title_margin_mm, unit = "mm")),
        strip.text.y = element_text(family = plot_font_family,
                                    size = plot_row_title_size,
                                    face = "bold",
                                    colour = plot_text_color),
        legend.position = "bottom",
        legend.title = element_text(family = plot_font_family,
                                    size = plot_legend_title_size,
                                    face = "bold",
                                    colour = plot_text_color),
        legend.text = element_text(family = plot_font_family,
                                   size = plot_legend_text_size,
                                   colour = plot_text_color),
        legend.box.background = element_rect(colour = plot_text_color,
                                             fill = NA,
                                             linewidth = legend_box_line_width),
        legend.box.margin = margin(legend_box_margin_mm,
                                   legend_box_margin_mm,
                                   legend_box_margin_mm,
                                   legend_box_margin_mm,
                                   unit = "mm"),
        panel.spacing = grid::unit(panel_spacing_mm, "mm"),
        plot.margin = margin(1, 1, 1, 1, unit = "mm")) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))


## Show and save figure
combined_plot
ggsave(filename = combined_svg,
       plot = combined_plot,
       device = svglite::svglite,
       width = plot_width_cm,
       height = plot_height_cm,
       units = "cm",
       bg = "white",
       limitsize = FALSE)
message("Combined SVG saved to: ", combined_svg)




#### Supplementary Figure S8 ###################################################

## Install packages
rm(list = ls()) #clear environment
CRAN_packages <- c("tibble", "dplyr", "vcfR") #CRAN packages
Bioconductor_packages <- "SeqArray" #Bioconductor package
for(p in CRAN_packages) if(!requireNamespace(p, quietly = TRUE)) install.packages(p) #install missing CRAN packages
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") #install BiocManager
if(!requireNamespace(Bioconductor_packages, quietly = TRUE)) BiocManager::install(Bioconductor_packages, ask = FALSE, update = FALSE) #install missing Bioconductor package
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv") #import (NicheDiv::transform.skewed.variables function from my other package (under development)
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")


## Set directories
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
intermediate_files_folder <- "Empirical_examples/Intermediate_files"
figure_files_folder <- "Figure_files"


## Import and prepare data
Monticola71_data <- read.csv(file = "./Empirical_examples/Pyron_2023/monticola71.csv",
                             row.names = 1,
                             header = TRUE,
                             colClasses = c(huc2 = "character",
                                            huc4 = "character",
                                            huc6 = "character",
                                            huc8 = "character",
                                            huc10 = "character",
                                            huc12 = "character"))
Monticola71_SNP <- process.SNP.data.SOM(vcf.path = "./Empirical_examples/Pyron_2023/Monticola71.vcf.gz", #filter loci and individuals and create SNP matrix dataframe
                                        missing.loci.cutoff.lenient = 0.7,
                                        missing.loci.cutoff.final = 0.5,
                                        missing.individuals.cutoff = 0.6)
Monticola71_snp_to_sample <- Monticola71_data$Sample[match(rownames(Monticola71_SNP), rownames(Monticola71_data))] #returns RAP**** names
rownames(Monticola71_SNP) <- Monticola71_snp_to_sample #rename SNP matrix to RAP codes
Monticola71_spatial <- Monticola71_data[, c("Lat", "Long", "Elev"), drop = FALSE] #extract numeric spatial data
Monticola71_spatial <- dplyr::rename(Monticola71_spatial, Latitude = Lat, Longitude = Long, Elevation = Elev)
rownames(Monticola71_spatial) <- Monticola71_data$Sample #assign rownames
nrow(Monticola71_spatial) #number of samples: 71
Monticola71_environmental <- read.csv("./Empirical_examples/Pyron_2023/Monticola71_environmental.csv", header = TRUE) #read CSV
rownames(Monticola71_environmental) <- Monticola71_environmental$Sample
Monticola71_environmental <- Monticola71_environmental[, !names(Monticola71_environmental) %in% c("Sample", "ID")] #remove ID columns
Monticola71_environmental <- Monticola71_environmental[, !names(Monticola71_environmental) %in% c("Latitude", "Longitude", "Elevation")] #remove spatial variables
Monticola71_environmental[] <- lapply(Monticola71_environmental, as.numeric) #ensure numeric
rownames(Monticola71_data) <- Monticola71_data$Sample #assign rownames
Monticola71_watershed <- make.cols.binary.SOM(dataframe = Monticola71_data, make.binary.cols = c("lvl4", "huc2", "huc4", "huc6", "huc8", "huc10", "huc12"))
Monticola71_watershed <- Monticola71_watershed[rownames(Monticola71_data), , drop = FALSE]
Monticola71_environmental <- (NicheDiv::transform.skewed.variables(Monticola71_environmental))$transformed #transform skewed variables
Monticola71_environmental <- remove.lowCV.multicollinearity.SOM(Monticola71_environmental, #remove highly correlated and low-variance variables
                                                                CV.threshold = 0.05,
                                                                cor.threshold = 0.9)
Monticola71_trait_names <- c("SVL", "TL", "AG", "CW", "FL", "HL", "SG", "TW", "TO", "FI", "HW", "ED", "IN", "ES", "ON", "IO", "IC")
Monticola71_trait_data <- Monticola71_data[, Monticola71_trait_names] #extract variables
Monticola71_log_traits <- log(Monticola71_trait_data) #no log-transformation
Monticola71_filtered_log_traits <- remove.lowCV.multicollinearity.SOM(Monticola71_log_traits, #filter log-transformed traits by CV and correlation, excluding SVL from removal
                                                                      CV.threshold = 0.05,
                                                                      cor.threshold = 0.9,
                                                                      exclude.cols = "SVL")
Monticola71_SVL <- Monticola71_filtered_log_traits[, "SVL"] #extract SVL and residualize all others
Monticola71_residuals_mat <- sapply(colnames(Monticola71_filtered_log_traits)[colnames(Monticola71_filtered_log_traits) != "SVL"], function(trait) {stats::resid(stats::lm(Monticola71_filtered_log_traits[, trait] ~ Monticola71_SVL))}) #regress each trait on SVL and store residuals
rownames(Monticola71_filtered_log_traits) <- Monticola71_data$Sample #set rownames for log-transformed traits
rownames(Monticola71_residuals_mat) <- Monticola71_data$Sample #set rownames for residualized traits
Monticola71_morphology <- as.data.frame(cbind(SVL = Monticola71_SVL, Monticola71_residuals_mat)) #combine log(SVL) and residuals
Monticola71_SOM_data <- list(SNP = Monticola71_SNP,
                             Spatial = Monticola71_spatial,
                             Environmental = Monticola71_environmental,
                             Watershed = Monticola71_watershed,
                             Morphology = Monticola71_morphology)
Monticola71_SOM_tr <- train.SOM(input_data = Monticola71_SOM_data, #71 samples, 14.7min
                                save.SOM.results = TRUE,
                                save.SOM.results.name = file.path(intermediate_files_folder, "Monticola71_SOM_tr.Rdata"),
                                max.NA.row = 0.6,
                                max.NA.col = 0.5)
Monticola71_SOM <- clustering.SOM(Monticola71_SOM_tr, #15.8min
                                  clustering.method = "kmeans+BICelbow",
                                  save.SOM.results.name = file.path(intermediate_files_folder, "Monticola71_SOM_kmeansBICelbow.Rdata"))
Monticola71_SOM$optim_k_summary #k2 100%


## Supplementary Figure S8A & Supplementary Figure S8B
plot_width_cm <- 12.38
plot_height_cm <- 4
panel_gap_cm <- 0.52
boundary_color_clusters <- "red"
boundary_line_width_clusters <- 3
point_color_clusters <- "white"
point_shape_clusters <- 19
point_size_clusters <- 0.9
cluster_shape_clusters <- "straight"
cluster_shape_neighbor_dist <- "straight"
legend_font_size <- 7
figure_name <- "Supplementary_Figure_S8AB.svg"
calc_unit_neighbor_dist <- function(som_model) {
  number_of_units <- nrow(som_model$grid$pts)
  codebook_distance_matrix <- as.matrix(kohonen::object.distances(som_model, type = "codes"))
  grid_distance_matrix <- as.matrix(kohonen::unit.distances(som_model$grid))
  neighbor_matrix <- abs(grid_distance_matrix - 1) <= 0.001
  codebook_distance_matrix[!neighbor_matrix] <- NA_real_
  unit_mean_neighbor_distances <- colMeans(codebook_distance_matrix, na.rm = TRUE)
  unit_mean_neighbor_distances[!is.finite(unit_mean_neighbor_distances)] <- NA_real_
  unit_mean_neighbor_distances
}
count_SOM_clusters <- function(cluster_vector) {
  cluster_vector <- as.integer(cluster_vector)
  cluster_vector <- cluster_vector[is.finite(cluster_vector) & !is.na(cluster_vector) & cluster_vector >= 1]
  if (length(cluster_vector) == 0) return(NA_integer_)
  length(unique(cluster_vector))
}
choose_representative_replicate <- function(som_models, som_clusters) {
  number_of_replicates <- length(som_clusters)
  if (number_of_replicates == 1) return(1L)
  sample_cluster_assignments <- vector("list", number_of_replicates)
  for (replicate_index in seq_len(number_of_replicates)) {
    unit_classif <- as.integer(som_models[[replicate_index]]$unit.classif)
    unit_cluster_labels <- as.integer(som_clusters[[replicate_index]])
    sample_cluster_assignments[[replicate_index]] <- unit_cluster_labels[unit_classif]
  }
  k_values <- vapply(som_clusters, count_SOM_clusters, integer(1))
  k_frequency <- table(k_values[!is.na(k_values)])
  modal_k_values <- as.integer(names(k_frequency)[k_frequency == max(k_frequency)])
  selected_k <- min(modal_k_values)
  candidate_replicates <- which(k_values == selected_k)
  if (length(candidate_replicates) == 1) return(candidate_replicates)
  pairwise_adjusted_rand_index <- matrix(NA_real_, nrow = length(candidate_replicates), ncol = length(candidate_replicates))
  diag(pairwise_adjusted_rand_index) <- NA_real_
  for (candidate_index_1 in seq_len(length(candidate_replicates) - 1)) {
    for (candidate_index_2 in seq.int(candidate_index_1 + 1, length(candidate_replicates))) {
      replicate_index_1 <- candidate_replicates[candidate_index_1]
      replicate_index_2 <- candidate_replicates[candidate_index_2]
      current_adjusted_rand_index <- mclust::adjustedRandIndex(sample_cluster_assignments[[replicate_index_1]], sample_cluster_assignments[[replicate_index_2]])
      pairwise_adjusted_rand_index[candidate_index_1, candidate_index_2] <- current_adjusted_rand_index
      pairwise_adjusted_rand_index[candidate_index_2, candidate_index_1] <- current_adjusted_rand_index
    }
  }
  mean_adjusted_rand_index <- rowMeans(pairwise_adjusted_rand_index, na.rm = TRUE)
  if (all(!is.finite(mean_adjusted_rand_index) | is.na(mean_adjusted_rand_index))) return(candidate_replicates[1])
  representative_candidate_index <- which.max(replace(mean_adjusted_rand_index, !is.finite(mean_adjusted_rand_index) | is.na(mean_adjusted_rand_index), -Inf))
  candidate_replicates[representative_candidate_index]
}
som_models_use <- Monticola71_SOM$som_models
som_clusters_use <- Monticola71_SOM$som_clusters
replicate_index <- choose_representative_replicate(som_models = som_models_use, som_clusters = som_clusters_use)
som_model <- som_models_use[[replicate_index]]
som_cluster <- as.integer(som_clusters_use[[replicate_index]])
neighbor_distances <- calc_unit_neighbor_dist(som_model)
som_cluster <- as.integer(factor(som_cluster, levels = sort(unique(som_cluster))))
number_of_clusters <- length(unique(som_cluster))
cluster_colors <- viridis::viridis(number_of_clusters)
SOM_cluster_plot_col <- cluster_colors[som_cluster]
device_width_in <- (plot_width_cm / 2.54) * (96 / 72)
device_height_in <- (plot_height_cm / 2.54) * (96 / 72)
measurement_file <- tempfile(fileext = ".svg")
svg(measurement_file, width = device_width_in, height = device_height_in, family = "Arial")
base_font_size <- par("ps")
legend_relative_font_size <- (legend_font_size * (96 / 72)) / base_font_size
line_height_in <- par("csi") * par("mex")
panel_gap_in <- (panel_gap_cm / 2.54) * (96 / 72)
neighbor_panel_margin_width_in <- (4.6 + 0.6) * line_height_in
cluster_panel_margin_width_in <- (0.6 + 0.6) * line_height_in
map_plot_width_in <- (device_width_in - panel_gap_in - neighbor_panel_margin_width_in - cluster_panel_margin_width_in) / 2
neighbor_panel_width_in <- map_plot_width_in + neighbor_panel_margin_width_in
neighbor_panel_end_initial <- neighbor_panel_width_in / device_width_in
cluster_panel_start_initial <- (neighbor_panel_width_in + panel_gap_in) / device_width_in
par(family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(fig = c(0, neighbor_panel_end_initial, 0, 1), new = FALSE)
plot(x = som_model, type = "property", property = neighbor_distances, main = "", shape = cluster_shape_neighbor_dist, cex = legend_relative_font_size, palette.name = function(n) rev(viridis::cividis(n)))
left_map_right <- grconvertX(max(som_model$grid$pts[, 1]) + 0.5, from = "user", to = "ndc")
par(fig = c(cluster_panel_start_initial, 1, 0, 1), new = TRUE)
plot(x = som_model, shape = cluster_shape_clusters, type = "mapping", bgcol = SOM_cluster_plot_col, main = "", pchs = point_shape_clusters, cex = point_size_clusters, col = point_color_clusters)
right_map_left <- grconvertX(min(som_model$grid$pts[, 1]) - 0.5, from = "user", to = "ndc")
dev.off()
unlink(measurement_file)
actual_gap_fraction <- right_map_left - left_map_right
desired_gap_fraction <- panel_gap_cm / plot_width_cm
panel_shift_fraction <- (actual_gap_fraction - desired_gap_fraction) / 2
neighbor_panel_start <- panel_shift_fraction
neighbor_panel_end <- neighbor_panel_end_initial + panel_shift_fraction
cluster_panel_start <- cluster_panel_start_initial - panel_shift_fraction
cluster_panel_end <- 1 - panel_shift_fraction
svg(file.path(figure_files_folder, figure_name), width = device_width_in, height = device_height_in, family = "Arial")
base_font_size <- par("ps")
legend_relative_font_size <- (legend_font_size * (96 / 72)) / base_font_size
par(family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(fig = c(neighbor_panel_start, neighbor_panel_end, 0, 1), new = FALSE)
plot(x = som_model, type = "property", property = neighbor_distances, main = "", shape = cluster_shape_neighbor_dist, cex = legend_relative_font_size, palette.name = function(n) rev(viridis::cividis(n)))
par(fig = c(cluster_panel_start, cluster_panel_end, 0, 1), new = TRUE)
plot(x = som_model, shape = cluster_shape_clusters, type = "mapping", bgcol = SOM_cluster_plot_col, main = "", pchs = point_shape_clusters, cex = point_size_clusters, col = point_color_clusters)
if (number_of_clusters > 1) kohonen::add.cluster.boundaries(som_model, som_cluster, lwd = boundary_line_width_clusters, col = boundary_color_clusters)
dev.off()


## Supplementary Figure S8C
plot_width_cm <- 5.68
plot_height_cm <- 10.09
panel_heights <- c(1, 1, 1.1)
panel_gap <- 0.85
side_tick_label_gap <- 0.65
bottom_tick_label_gap <- 0.5
figure_name <- "Supplementary_Figure_S8C.svg"

svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72))
max_k <- Monticola71_SOM$max_k
optim_k_vals <- as.numeric(Monticola71_SOM$optim_k_vals)
BIC_values <- Monticola71_SOM$BIC_values[seq_len(max_k), , drop = FALSE]
support_values <- if (!is.null(Monticola71_SOM$support_values)) Monticola71_SOM$support_values[seq_len(max_k), , drop = FALSE] else BIC_values
k_colors <- viridis::magma(max_k)
base_font_size <- par("ps")
bottom_numbers_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
layout(matrix(1:3, ncol = 1), heights = panel_heights)
par(bty = "n", oma = c(0, 0, 0, 0), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mar = c(panel_gap, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
finite_k_rows <- apply(support_values, 1, function(x) any(is.finite(x) & !is.na(x)))
plotted_k_values <- seq_len(max_k)[finite_k_rows]
values_for_plot <- t(support_values[finite_k_rows, , drop = FALSE])
boxplot(values_for_plot, at = plotted_k_values, xlim = c(0.5, max_k + 0.5), outline = FALSE, notch = FALSE, axes = FALSE, ylab = "", main = "", whisklty = 1, staplelty = 1, col = k_colors[plotted_k_values])
axis(1, at = seq_len(max_k), labels = FALSE)
y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = y_axis_breaks, labels = round(y_axis_breaks, 1), las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
par(mar = c(panel_gap, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
delta_BIC_matrix <- apply(BIC_values, 2, function(x) { previous_BIC <- x[-length(x)]; current_BIC <- x[-1]; delta_BIC <- previous_BIC - current_BIC; delta_BIC[!is.finite(previous_BIC) | !is.finite(current_BIC)] <- NA_real_; delta_BIC })
if (is.null(dim(delta_BIC_matrix))) delta_BIC_matrix <- matrix(delta_BIC_matrix, ncol = 1)
finite_delta_rows <- apply(delta_BIC_matrix, 1, function(x) any(is.finite(x) & !is.na(x)))
plotted_delta_k_values <- seq.int(2, max_k)[finite_delta_rows]
delta_BIC_for_plot <- t(delta_BIC_matrix[finite_delta_rows, , drop = FALSE])
boxplot(delta_BIC_for_plot, at = plotted_delta_k_values, xlim = c(0.5, max_k + 0.5), outline = FALSE, notch = FALSE, axes = FALSE, ylab = "", main = "", whisklty = 1, staplelty = 1, col = k_colors[plotted_delta_k_values])
axis(1, at = seq_len(max_k), labels = FALSE)
y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = y_axis_breaks, labels = round(y_axis_breaks, 1), las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
par(mar = c(2, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
k_frequency_values <- table(factor(optim_k_vals, levels = seq_len(max_k))) / length(optim_k_vals)
bar_midpoints <- barplot(k_frequency_values, ylim = c(0, 1), col = k_colors, axes = FALSE, axisnames = FALSE, ylab = "", main = "")
axis(1, at = bar_midpoints, labels = seq_len(max_k), mgp = c(3, bottom_tick_label_gap, 0), font.axis = 1, cex.axis = bottom_numbers_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
dev.off()



## Supplementary Figure S8D
plot_width_cm <- 9.1
plot_height_cm <- 7.88
panel_gap <- 1.55
side_tick_label_gap <- 0.65
bottom_tick_label_gap <- 0.75
top_margin_mm <- 2
right_margin_mm <- 2
figure_name <- "Supplementary_Figure_S8D.svg"
load(file.path(intermediate_files_folder, "Monticola71_SOM_lolo.Rdata"))
layer_summary_table <- leave.one.layer.out.results$layer.summary
replicate_matched_results_table <- leave.one.layer.out.results$replicate.matched.results
successful_replicate_matched_results_table <- replicate_matched_results_table[is.na(replicate_matched_results_table$error) | replicate_matched_results_table$error == "", , drop = FALSE]
SOM_layer_names <- as.character(layer_summary_table$layer)
successful_replicate_matched_results_table$layer <- factor(successful_replicate_matched_results_table$layer, levels = SOM_layer_names)
layer_colors <- setNames(viridis::turbo(length(SOM_layer_names)), SOM_layer_names)
show_assignment_margin_plot <- any(is.finite(successful_replicate_matched_results_table$delta.mean.assignment.margin) & !is.na(successful_replicate_matched_results_table$delta.mean.assignment.margin))
svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72))
base_font_size <- par("ps")
axis_labels_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
if (show_assignment_margin_plot) {
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1))
} else {
  layout(matrix(1:2, nrow = 1), widths = c(1, 1))
}
par(bty = "o", oma = c(7, 0, top_margin_lines, right_margin_lines), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mar = c(0, 2, 0, panel_gap), mgp = c(3, side_tick_label_gap, 0))
boxplot(absolute.k.deviation ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
for (layer_index in seq_along(SOM_layer_names)) {
  current_layer_name <- SOM_layer_names[layer_index]
  current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "absolute.k.deviation"]
  current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
  if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
}
box()
par(mar = c(0, 2, 0, if (show_assignment_margin_plot) panel_gap else 0), mgp = c(3, side_tick_label_gap, 0))
boxplot(pairwise.coassignment.change ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
for (layer_index in seq_along(SOM_layer_names)) {
  current_layer_name <- SOM_layer_names[layer_index]
  current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "pairwise.coassignment.change"]
  current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
  if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
}
box()
if (show_assignment_margin_plot) {
  par(mar = c(0, 2, 0, 0), mgp = c(3, side_tick_label_gap, 0))
  boxplot(delta.mean.assignment.margin ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
  axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
  axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
  for (layer_index in seq_along(SOM_layer_names)) {
    current_layer_name <- SOM_layer_names[layer_index]
    current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "delta.mean.assignment.margin"]
    current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
    if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
  }
  box()
}
dev.off()


## Supplementary Figure S8E
plot_width_cm <- 10.42
plot_height_cm <- 7.38
row_gap <- 1.45
column_gap <- 2.5
bottom_tick_label_gap <- 0.6
top_margin_mm <- 2
left_margin_mm <- 2
right_margin_mm <- 2
figure_name <- "Supplementary_Figure_S8E.svg"
bars_threshold_N <- 20
importance_threshold <- 0.0001
layer_label_font_size <- 9
bar_label_font_size <- 7.1
axis_ticks_font_size <- 7.1
variable_label_abbreviations <- c("Latitude" = "Lat",
                                  "Longitude" = "Long",
                                  "Elevation" = "Elev")

matrix_names <- Monticola71_SOM$input_data_names
first_codebook_list <- kohonen::getCodes(Monticola71_SOM$som_models[[1]])
if (!is.list(first_codebook_list)) first_codebook_list <- list(first_codebook_list)
number_of_layers <- length(first_codebook_list)
if (length(matrix_names) != number_of_layers) matrix_names <- paste0("layer", seq_len(number_of_layers))
for (layer_index in seq_len(number_of_layers)) {
  if (is.null(colnames(first_codebook_list[[layer_index]]))) colnames(first_codebook_list[[layer_index]]) <- paste0("V", seq_len(ncol(first_codebook_list[[layer_index]])))
}
retained_replicate_indices <- seq_along(Monticola71_SOM$som_models)
Monticola71_SOM_variable_importance <- vector("list", number_of_layers)
names(Monticola71_SOM_variable_importance) <- matrix_names
for (layer_index in seq_len(number_of_layers)) {
  Monticola71_SOM_variable_importance[[layer_index]] <- matrix(NA_real_, nrow = length(retained_replicate_indices), ncol = ncol(first_codebook_list[[layer_index]]), dimnames = list(paste0("R", retained_replicate_indices), colnames(first_codebook_list[[layer_index]])))
}
for (retained_replicate_position in seq_along(retained_replicate_indices)) {
  retained_replicate_index <- retained_replicate_indices[retained_replicate_position]
  som_model <- Monticola71_SOM$som_models[[retained_replicate_index]]
  neuron_cluster_vector <- Monticola71_SOM$som_clusters[[retained_replicate_index]]
  codebook_list <- kohonen::getCodes(som_model)
  if (!is.list(codebook_list)) codebook_list <- list(codebook_list)
  for (layer_index in seq_len(number_of_layers)) {
    if (is.null(colnames(codebook_list[[layer_index]]))) colnames(codebook_list[[layer_index]]) <- paste0("V", seq_len(ncol(codebook_list[[layer_index]])))
    codebook_matrix <- as.matrix(codebook_list[[layer_index]])
    valid_cluster_rows <- is.finite(neuron_cluster_vector) & !is.na(neuron_cluster_vector)
    neuron_sample_counts <- tabulate(som_model$unit.classif, nbins = length(neuron_cluster_vector))
    codebook_matrix <- codebook_matrix[valid_cluster_rows, , drop = FALSE]
    cluster_labels <- neuron_cluster_vector[valid_cluster_rows]
    neuron_sample_counts <- neuron_sample_counts[valid_cluster_rows]
    neuron_weights <- neuron_sample_counts + 1
    Monticola71_SOM_variable_importance[[layer_index]][retained_replicate_position, ] <- apply(codebook_matrix, 2, function(variable_values) {
      valid_rows <- is.finite(variable_values) & !is.na(variable_values) & is.finite(neuron_weights) & !is.na(neuron_weights)
      variable_values <- variable_values[valid_rows]
      variable_cluster_labels <- cluster_labels[valid_rows]
      variable_weights <- neuron_weights[valid_rows]
      if (length(variable_values) < 2 || length(unique(variable_cluster_labels)) < 2 || sum(variable_weights) <= 0) return(NA_real_)
      weighted_variable_mean <- sum(variable_weights * variable_values) / sum(variable_weights)
      total_sum_of_squares <- sum(variable_weights * (variable_values - weighted_variable_mean)^2)
      if (!is.finite(total_sum_of_squares) || total_sum_of_squares <= 0) return(0)
      cluster_means <- tapply(variable_weights * variable_values, variable_cluster_labels, sum) / tapply(variable_weights, variable_cluster_labels, sum)
      cluster_sizes <- tapply(variable_weights, variable_cluster_labels, sum)
      between_cluster_sum_of_squares <- sum(cluster_sizes * (cluster_means - weighted_variable_mean)^2)
      as.numeric(between_cluster_sum_of_squares / total_sum_of_squares)
    })
  }
}
svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72), family = "Arial")
base_font_size <- par("ps")
layer_label_relative_font_size <- (layer_label_font_size * (96 / 72)) / base_font_size
bar_label_relative_font_size <- (bar_label_font_size * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (axis_ticks_font_size * (96 / 72)) / base_font_size
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
left_margin_lines <- (left_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
par(mfrow = if (number_of_layers <= 3) c(1, number_of_layers) else if (number_of_layers == 4) c(2, 2) else if (number_of_layers <= 6) c(2, 3) else if (number_of_layers <= 8) c(2, 4) else if (number_of_layers == 9) c(3, 3) else c(ceiling(number_of_layers / 3), 3), oma = c(0, left_margin_lines, top_margin_lines, right_margin_lines), mar = c(2.2, column_gap / 2, row_gap, column_gap / 2), mgp = c(3, bottom_tick_label_gap, 0), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1, bty = "o", family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black")
layer_colors <- setNames(viridis::turbo(length(matrix_names)), matrix_names)
for (layer_index in seq_along(Monticola71_SOM_variable_importance)) {
  variable_importance_matrix <- Monticola71_SOM_variable_importance[[layer_index]]
  if (is.null(variable_importance_matrix) || nrow(variable_importance_matrix) == 0 || ncol(variable_importance_matrix) == 0) {
    plot.new()
    next
  }
  median_metric_per_variable <- matrixStats::colMedians(variable_importance_matrix, na.rm = TRUE, useNames = TRUE)
  median_metric_per_variable[!is.finite(median_metric_per_variable)] <- NA_real_
  retained_variable_indices <- which(is.finite(median_metric_per_variable) & !is.na(median_metric_per_variable) & median_metric_per_variable > importance_threshold)
  if (length(retained_variable_indices) == 0) {
    plot.new()
    next
  }
  variable_importance_matrix <- variable_importance_matrix[, retained_variable_indices, drop = FALSE]
  retained_variable_medians <- median_metric_per_variable[colnames(variable_importance_matrix)]
  variable_sort_indices <- order(retained_variable_medians, decreasing = FALSE)
  variable_importance_matrix <- variable_importance_matrix[, variable_sort_indices, drop = FALSE]
  variable_labels <- colnames(variable_importance_matrix)
  abbreviation_matches <- match(variable_labels, names(variable_label_abbreviations))
  variable_labels[!is.na(abbreviation_matches)] <- variable_label_abbreviations[abbreviation_matches[!is.na(abbreviation_matches)]]
  number_of_bars <- ncol(variable_importance_matrix)
  boxplot(variable_importance_matrix, horizontal = TRUE, las = 1, notch = FALSE, outline = FALSE, col = rep(layer_colors[matrix_names[layer_index]], number_of_bars), border = "black", axes = FALSE, whisklty = if (number_of_bars > bars_threshold_N) 0 else 1, staplelty = if (number_of_bars > bars_threshold_N) 0 else 1, names = rep("", number_of_bars))
  axis(1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size, col = "black", col.axis = "black", font.axis = 1)
  axis(2, at = seq_len(number_of_bars), labels = if (number_of_bars <= bars_threshold_N) variable_labels else FALSE, las = 1, tick = FALSE, cex.axis = bar_label_relative_font_size, col.axis = "black", font.axis = 1)
  mtext(matrix_names[layer_index], side = 3, line = 0.3, cex = layer_label_relative_font_size, font = 1, family = "Arial", col = "black")
}
dev.off()


## Supplementary Figure S8F
plot_width_cm <- 8.93
plot_height_cm <- 12.7
pie_size <- 2
lon_buffer_range <- 1.5
lat_buffer_range <- 3.1
scale_size <- 0.3
bottom_tick_label_gap <- 0.55
side_tick_label_gap <- 0.8
north_arrow_position <- c(0.93, 0.13)
north_arrow_length <- 1.05
north_arrow_N_position <- 0.39
scale_position <- c(0.08, 0.085)
figure_name <- "Supplementary_Figure_S8F.svg"

Coordinates <- Monticola71_spatial[, c("Latitude", "Longitude")]
ancestry_matrix <- as.matrix(Monticola71_SOM$ancestry_matrix)
Coordinates <- as.data.frame(Coordinates, stringsAsFactors = FALSE)
Coordinates$Latitude <- as.numeric(Coordinates$Latitude)
Coordinates$Longitude <- as.numeric(Coordinates$Longitude)
coordinate_sample_names <- rownames(Coordinates)
ancestry_sample_names <- rownames(ancestry_matrix)
matched_sample_names <- intersect(ancestry_sample_names, coordinate_sample_names)
Coordinates <- Coordinates[matched_sample_names, , drop = FALSE]
ancestry_matrix <- ancestry_matrix[matched_sample_names, , drop = FALSE]
rows_with_missing_coordinates <- which(!is.finite(Coordinates$Latitude) | !is.finite(Coordinates$Longitude))
if (length(rows_with_missing_coordinates) > 0) {
  Coordinates <- Coordinates[-rows_with_missing_coordinates, , drop = FALSE]
  ancestry_matrix <- ancestry_matrix[-rows_with_missing_coordinates, , drop = FALSE]
}
ancestry_row_sums <- rowSums(ancestry_matrix, na.rm = TRUE)
ancestry_matrix <- ancestry_matrix / ancestry_row_sums
ancestry_proportions <- as.data.frame(ancestry_matrix)
cluster_colors <- viridis::viridis(ncol(ancestry_matrix))
longitude_minimum <- min(Coordinates$Longitude) - lon_buffer_range
longitude_maximum <- max(Coordinates$Longitude) + lon_buffer_range
latitude_minimum <- min(Coordinates$Latitude) - lat_buffer_range
latitude_maximum <- max(Coordinates$Latitude) + lat_buffer_range
pie_radius <- pie_size * 0.01 * max(longitude_maximum - longitude_minimum, latitude_maximum - latitude_minimum)
add_admixture_pie <- function(longitude, latitude, ancestry_proportions, cluster_colors, x.radius, y.radius, border.color = "black", line.width = 0.8, number.of.points = 80) {
  ancestry_proportions <- as.numeric(ancestry_proportions)
  ancestry_proportions[is.na(ancestry_proportions) | !is.finite(ancestry_proportions) | ancestry_proportions < 0] <- 0
  if (sum(ancestry_proportions) <= 0) return(invisible(NULL))
  ancestry_proportions <- ancestry_proportions / sum(ancestry_proportions)
  slice_start_angles <- c(0, cumsum(ancestry_proportions)[-length(ancestry_proportions)]) * 2 * pi
  slice_end_angles <- cumsum(ancestry_proportions) * 2 * pi
  for (slice_index in seq_along(ancestry_proportions)) {
    if (ancestry_proportions[slice_index] <= 0) next
    slice_angles <- seq(slice_start_angles[slice_index], slice_end_angles[slice_index], length.out = max(3, ceiling(number.of.points * ancestry_proportions[slice_index])))
    polygon(x = c(longitude, longitude + x.radius * cos(slice_angles), longitude), y = c(latitude, latitude + y.radius * sin(slice_angles), latitude), col = cluster_colors[slice_index], border = border.color, lwd = line.width)
  }
  circle_angles <- seq(0, 2 * pi, length.out = number.of.points)
  lines(longitude + x.radius * cos(circle_angles), latitude + y.radius * sin(circle_angles), col = border.color, lwd = line.width)
  invisible(NULL)
}
add_map_legend <- function(legend.position, legend.title, legend.labels, legend.colors, legend.text.relative.font.size, legend.title.relative.font.size, legend.text.font, legend.symbol.size, legend.box) {
  plot_coordinate_limits <- par("usr")
  plot_longitude_range <- plot_coordinate_limits[2] - plot_coordinate_limits[1]
  plot_latitude_range <- plot_coordinate_limits[4] - plot_coordinate_limits[3]
  legend_text_width <- max(strwidth(legend.labels, units = "user", cex = legend.text.relative.font.size, font = legend.text.font, family = "Arial"))
  legend_text_height <- strheight("M", units = "user", cex = legend.text.relative.font.size, font = legend.text.font, family = "Arial")
  legend_title_width <- 0
  legend_title_height <- 0
  legend_title_gap <- 0
  if (!is.null(legend.title)) {
    legend_title_width <- strwidth(legend.title, units = "user", cex = legend.title.relative.font.size, font = 2, family = "Arial")
    legend_title_height <- strheight("M", units = "user", cex = legend.title.relative.font.size, font = 2, family = "Arial")
    legend_title_gap <- 0.5 * legend_text_height
  }
  legend_symbol_width <- strwidth("M", units = "user", cex = legend.text.relative.font.size, family = "Arial") * legend.symbol.size
  legend_symbol_gap <- 0.7 * strwidth("M", units = "user", cex = legend.text.relative.font.size, family = "Arial")
  legend_padding_x <- 0.9 * strwidth("M", units = "user", cex = legend.text.relative.font.size, family = "Arial")
  legend_padding_y <- 0.7 * legend_text_height
  legend_line_gap <- 0.35 * legend_text_height
  legend_width <- max(legend_title_width, legend_symbol_width + legend_symbol_gap + legend_text_width) + 2 * legend_padding_x
  legend_height <- 2 * legend_padding_y + legend_title_height + legend_title_gap + length(legend.labels) * legend_text_height + (length(legend.labels) - 1) * legend_line_gap
  if (legend.position == "topright") {
    legend_left <- plot_coordinate_limits[2] - legend_width
    legend_bottom <- plot_coordinate_limits[4] - legend_height
  } else if (legend.position == "topleft") {
    legend_left <- plot_coordinate_limits[1]
    legend_bottom <- plot_coordinate_limits[4] - legend_height
  } else if (legend.position == "bottomright") {
    legend_left <- plot_coordinate_limits[2] - legend_width
    legend_bottom <- plot_coordinate_limits[3]
  } else if (legend.position == "bottomleft") {
    legend_left <- plot_coordinate_limits[1]
    legend_bottom <- plot_coordinate_limits[3]
  } else if (legend.position == "right") {
    legend_left <- plot_coordinate_limits[2] - legend_width
    legend_bottom <- plot_coordinate_limits[3] + 0.5 * plot_latitude_range - 0.5 * legend_height
  } else if (legend.position == "left") {
    legend_left <- plot_coordinate_limits[1]
    legend_bottom <- plot_coordinate_limits[3] + 0.5 * plot_latitude_range - 0.5 * legend_height
  } else if (legend.position == "top") {
    legend_left <- plot_coordinate_limits[1] + 0.5 * plot_longitude_range - 0.5 * legend_width
    legend_bottom <- plot_coordinate_limits[4] - legend_height
  } else if (legend.position == "bottom") {
    legend_left <- plot_coordinate_limits[1] + 0.5 * plot_longitude_range - 0.5 * legend_width
    legend_bottom <- plot_coordinate_limits[3]
  } else {
    legend_left <- plot_coordinate_limits[1] + 0.5 * plot_longitude_range - 0.5 * legend_width
    legend_bottom <- plot_coordinate_limits[3] + 0.5 * plot_latitude_range - 0.5 * legend_height
  }
  legend_right <- legend_left + legend_width
  legend_top <- legend_bottom + legend_height
  if (legend.box) rect(legend_left, legend_bottom, legend_right, legend_top, col = "white", border = "black")
  current_legend_y_position <- legend_top - legend_padding_y
  if (!is.null(legend.title)) {
    text(x = legend_left + legend_padding_x, y = current_legend_y_position - 0.5 * legend_title_height, labels = legend.title, adj = c(0, 0.5), font = 2, cex = legend.title.relative.font.size, family = "Arial", col = "black")
    current_legend_y_position <- current_legend_y_position - legend_title_height - legend_title_gap
  }
  legend_symbol_x_position <- legend_left + legend_padding_x + 0.5 * legend_symbol_width
  legend_text_x_position <- legend_left + legend_padding_x + legend_symbol_width + legend_symbol_gap
  for (legend_index in seq_along(legend.labels)) {
    legend_entry_y_position <- current_legend_y_position - 0.5 * legend_text_height - (legend_index - 1) * (legend_text_height + legend_line_gap)
    points(x = legend_symbol_x_position, y = legend_entry_y_position, pch = 21, cex = legend.symbol.size, bg = legend.colors[legend_index], col = "black")
    text(x = legend_text_x_position, y = legend_entry_y_position, labels = legend.labels[legend_index], adj = c(0, 0.5), cex = legend.text.relative.font.size, font = legend.text.font, family = "Arial", col = "black")
  }
  invisible(NULL)
}
svg(file.path(figure_files_folder, figure_name), width = plot_width_cm / 2.54, height = plot_height_cm / 2.54, family = "Arial")
base_font_size <- par("ps")
scale_text_relative_font_size <- 7.1 / base_font_size
axis_labels_relative_font_size <- 9.1 / base_font_size
legend_text_relative_font_size <- 9.1 / base_font_size
legend_title_relative_font_size <- 9.1 / base_font_size
par(mfrow = c(1, 1), oma = c(2, 2, 1, 1), mar = c(1, 1, 1, 1), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1, family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", bty = "o")
current_plot_region_size_inches <- par("pin")
longitude_range <- longitude_maximum - longitude_minimum
latitude_range <- latitude_maximum - latitude_minimum
mean_map_latitude <- mean(c(latitude_minimum, latitude_maximum))
longitude_latitude_correction <- cos(mean_map_latitude * pi / 180)
if (!is.finite(longitude_latitude_correction) || longitude_latitude_correction <= 0) longitude_latitude_correction <- 1
target_height_width_ratio <- latitude_range / (longitude_range * longitude_latitude_correction)
adjusted_plot_width_inches <- current_plot_region_size_inches[1]
adjusted_plot_height_inches <- adjusted_plot_width_inches * target_height_width_ratio
if (adjusted_plot_height_inches > current_plot_region_size_inches[2]) {
  adjusted_plot_height_inches <- current_plot_region_size_inches[2]
  adjusted_plot_width_inches <- adjusted_plot_height_inches / target_height_width_ratio
}
par(pin = c(adjusted_plot_width_inches, adjusted_plot_height_inches))
plot.new()
plot.window(xlim = c(longitude_minimum, longitude_maximum), ylim = c(latitude_minimum, latitude_maximum), xaxs = "i", yaxs = "i")
plot_coordinate_limits <- par("usr")
plot_region_size_inches <- par("pin")
x_units_per_inch <- plot_region_size_inches[1] / diff(plot_coordinate_limits[1:2])
y_units_per_inch <- plot_region_size_inches[2] / diff(plot_coordinate_limits[3:4])
if (!is.finite(x_units_per_inch) || !is.finite(y_units_per_inch) || x_units_per_inch <= 0 || y_units_per_inch <= 0) {
  pie_radius_x <- pie_radius
  pie_radius_y <- pie_radius
} else {
  pie_radius_x <- pie_radius * y_units_per_inch / x_units_per_inch
  pie_radius_y <- pie_radius
}
maps::map("world", fill = TRUE, col = "lightgrey", border = "black", xlim = c(longitude_minimum, longitude_maximum), ylim = c(latitude_minimum, latitude_maximum), add = TRUE)
try(maps::map("state", add = TRUE, xlim = c(longitude_minimum, longitude_maximum), ylim = c(latitude_minimum, latitude_maximum), col = "black", lwd = 0.5), silent = TRUE)
longitude_clip_buffer <- 0.2 * (longitude_maximum - longitude_minimum)
latitude_clip_buffer <- 0.2 * (latitude_maximum - latitude_minimum)
clip(longitude_minimum - longitude_clip_buffer, longitude_maximum + longitude_clip_buffer, latitude_minimum - latitude_clip_buffer, latitude_maximum + latitude_clip_buffer)
axis(1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size, col = "black", col.axis = "black", font.axis = 1)
axis(2, las = 2, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size, col = "black", col.axis = "black", font.axis = 1)
box(col = "black")
for (sample_index in seq_len(nrow(ancestry_proportions))) {
  add_admixture_pie(longitude = Coordinates$Longitude[sample_index], latitude = Coordinates$Latitude[sample_index], ancestry_proportions = as.numeric(ancestry_proportions[sample_index, ]), cluster_colors = cluster_colors, x.radius = pie_radius_x, y.radius = pie_radius_y)
}
species_by_sample <- as.character(Monticola71_data[matched_sample_names, "Species"])
species_names <- unique(species_by_sample)
species_cluster_means <- t(vapply(species_names, function(species_name) colMeans(ancestry_matrix[species_by_sample == species_name, , drop = FALSE], na.rm = TRUE), numeric(ncol(ancestry_matrix))))
cluster_species_names <- species_names[apply(species_cluster_means, 2, which.max)]
cluster_species_names <- gsub("_", " ", cluster_species_names)
cluster_species_names <- sub("^Desmognathus ", "", cluster_species_names)
cluster_species_names <- paste("D.", cluster_species_names)
legend_labels <- parse(text = paste0("italic('", cluster_species_names, "')"))
add_map_legend(legend.position = "topright", legend.title = "Species", legend.labels = legend_labels, legend.colors = cluster_colors, legend.text.relative.font.size = legend_text_relative_font_size, legend.title.relative.font.size = legend_title_relative_font_size, legend.text.font = 1, legend.symbol.size = 1.6, legend.box = TRUE)
scale_position_longitude <- scale_position[1] * (longitude_maximum - longitude_minimum) + longitude_minimum
scale_position_latitude <- scale_position[2] * (latitude_maximum - latitude_minimum) + latitude_minimum
maps::map.scale(x = scale_position_longitude, y = scale_position_latitude, cex = scale_text_relative_font_size, relwidth = scale_size, ratio = FALSE)
north_arrow_longitude <- north_arrow_position[1] * (longitude_maximum - longitude_minimum) + longitude_minimum
north_arrow_latitude <- north_arrow_position[2] * (latitude_maximum - latitude_minimum) + latitude_minimum
arrows(x0 = north_arrow_longitude, y0 = north_arrow_latitude, x1 = north_arrow_longitude, y1 = north_arrow_latitude + north_arrow_length, length = 0.13, col = "black", lwd = 1.7)
text(x = north_arrow_longitude, y = north_arrow_latitude + north_arrow_length + north_arrow_N_position, labels = "N", cex = 0.8, col = "black", family = "Arial", font = 1)
dev.off()




#### Supplementary Figure S9 ###################################################

## Install packages
rm(list = ls()) #clear environment
CRAN_packages <- c("tibble", "dplyr", "vcfR") #CRAN packages
Bioconductor_packages <- "SeqArray" #Bioconductor package
for(p in CRAN_packages) if(!requireNamespace(p, quietly = TRUE)) install.packages(p) #install missing CRAN packages
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") #install BiocManager
if(!requireNamespace(Bioconductor_packages, quietly = TRUE)) BiocManager::install(Bioconductor_packages, ask = FALSE, update = FALSE) #install missing Bioconductor package
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv") #import (NicheDiv::transform.skewed.variables function from my other package (under development)
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")


## Set directories
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
intermediate_files_folder <- "Empirical_examples/Intermediate_files"
figure_files_folder <- "Figure_files"


## Import and prepare data
Pascagoula_data <- read.csv(file = "./Empirical_examples/Pyron_et_al_2022/pascagoula22.csv",
                            row.names = 1,
                            header = TRUE, 
                            colClasses = c(huc2 = "character",
                                           huc4 = "character",
                                           huc6 = "character",
                                           huc8 = "character",
                                           huc10 = "character",
                                           huc12 = "character"))
Pascagoula_SNP <- process.SNP.data.SOM(vcf.path = "./Empirical_examples/Pyron_et_al_2022/pascagoula22.vcf.gz", #filter loci and individuals and create SNP matrix dataframe
                                       missing.loci.cutoff.lenient = 0.7,
                                       missing.loci.cutoff.final = 0.5,
                                       missing.individuals.cutoff = 0.5)
rownames(Pascagoula_SNP) <- Pascagoula_data$Sample[match(rownames(Pascagoula_SNP), rownames(Pascagoula_data))] #rename alleles
Pascagoula_spatial <- Pascagoula_data[, c("Lat", "Long", "Elev")] #extract coordinates and elevation
Pascagoula_spatial <- dplyr::rename(Pascagoula_spatial, Latitude = Lat, Longitude = Long, Elevation = Elev)
rownames(Pascagoula_spatial) <- Pascagoula_data$Sample #assign rownames
Pascagoula_environmental <- read.csv("./Empirical_examples/Pyron_et_al_2022/Pascagoula22_environmental.csv", header = TRUE) #read CSV
rownames(Pascagoula_environmental) <- Pascagoula_environmental$Sample
Pascagoula_environmental <- Pascagoula_environmental[, !names(Pascagoula_environmental) %in% c("Sample", "ID")] #remove ID columns
Pascagoula_environmental <- Pascagoula_environmental[, !names(Pascagoula_environmental) %in% c("Latitude", "Longitude", "Elevation")] #remove spatial variables
Pascagoula_environmental[] <- lapply(Pascagoula_environmental, as.numeric) #ensure numeric
rownames(Pascagoula_data) <- Pascagoula_data$Sample #assign rownames
Pascagoula_watershed <- make.cols.binary.SOM(dataframe = Pascagoula_data, make.binary.cols = c("lvl4", "huc2", "huc4", "huc6", "huc8", "huc10", "huc12"))
Pascagoula_watershed <- Pascagoula_watershed[rownames(Pascagoula_data), , drop = FALSE]
Pascagoula_environmental <- (NicheDiv::transform.skewed.variables(Pascagoula_environmental))$transformed #transform skewed variables
Pascagoula_environmental <- remove.lowCV.multicollinearity.SOM(Pascagoula_environmental, #remove highly correlated and low-variance variables
                                                               CV.threshold = 0.05,
                                                               cor.threshold = 0.9)
Pascagoula_trait_names <- c("SVL", "TL", "AG", "CW", "FL", "HL", "SG", "TW", "TO", "FI", "HW", "ED", "IN", "ES", "ON", "IO", "IC")
Pascagoula_trait_data <- Pascagoula_data[, Pascagoula_trait_names] #extract variables
rownames(Pascagoula_trait_data) <- Pascagoula_data$Sample #assign rownames
Pascagoula_log_traits <- log(Pascagoula_trait_data) #log-transform all traits
Pascagoula_filtered_log_traits <- remove.lowCV.multicollinearity.SOM(Pascagoula_log_traits, #filter log-transformed traits by CV and correlation, excluding SVL from removal
                                                                     CV.threshold = 0.05,
                                                                     cor.threshold = 0.9,
                                                                     exclude.cols = "SVL")
Pascagoula_SVL <- Pascagoula_filtered_log_traits[, "SVL"] #extract SVL and residualize all others
Pascagoula_residuals_mat <- sapply(colnames(Pascagoula_filtered_log_traits)[colnames(Pascagoula_filtered_log_traits) != "SVL"], function(trait) {stats::resid(stats::lm(Pascagoula_filtered_log_traits[, trait] ~ Pascagoula_SVL))}) #regress each trait on SVL and store residuals
rownames(Pascagoula_filtered_log_traits) <- Pascagoula_data$Sample #set rownames for log-transformed traits
rownames(Pascagoula_residuals_mat) <- Pascagoula_data$Sample #set rownames for residualized traits
Pascagoula_morphology <- as.data.frame(cbind(SVL = Pascagoula_SVL, Pascagoula_residuals_mat)) #combine log(SVL) and residuals
Pascagoula_SOM_data <- list(Alleles = Pascagoula_SNP,
                            Spatial = Pascagoula_spatial,
                            Environmental = Pascagoula_environmental,
                            Watershed = Pascagoula_watershed,
                            Morphology = Pascagoula_morphology)
Pascagoula_SOM_tr <- train.SOM(input_data = Pascagoula_SOM_data, #22 samples, 0.7min
                               max.NA.row = 0.5,
                               max.NA.col = 0.5,
                               save.SOM.results = TRUE,
                               save.SOM.results.name = file.path(intermediate_files_folder, "Pascagoula_SOM_tr.Rdata"),
                               grid.multiplier = 4)
Pascagoula_SOM <- clustering.SOM(Pascagoula_SOM_tr, #2.6min
                                 clustering.method = "kmeans+BICelbow",
                                 save.SOM.results.name = file.path(intermediate_files_folder, "Pascagoula_SOM_kmeansBICelbow.Rdata"))
Pascagoula_SOM$optim_k_summary #k2 99%


## Supplementary Figure S9A & Supplementary Figure S9B
plot_width_cm <- 12.84
plot_height_cm <- 4
panel_gap_cm <- 0.52
boundary_color_clusters <- "red"
boundary_line_width_clusters <- 3
point_color_clusters <- "white"
point_shape_clusters <- 19
point_size_clusters <- 0.9
cluster_shape_clusters <- "straight"
cluster_shape_neighbor_dist <- "straight"
legend_font_size <- 7
figure_name <- "Supplementary_Figure_S9AB.svg"

calc_unit_neighbor_dist <- function(som_model) {
  number_of_units <- nrow(som_model$grid$pts)
  codebook_distance_matrix <- as.matrix(kohonen::object.distances(som_model, type = "codes"))
  grid_distance_matrix <- as.matrix(kohonen::unit.distances(som_model$grid))
  neighbor_matrix <- abs(grid_distance_matrix - 1) <= 0.001
  codebook_distance_matrix[!neighbor_matrix] <- NA_real_
  unit_mean_neighbor_distances <- colMeans(codebook_distance_matrix, na.rm = TRUE)
  unit_mean_neighbor_distances[!is.finite(unit_mean_neighbor_distances)] <- NA_real_
  unit_mean_neighbor_distances
}
count_SOM_clusters <- function(cluster_vector) {
  cluster_vector <- as.integer(cluster_vector)
  cluster_vector <- cluster_vector[is.finite(cluster_vector) & !is.na(cluster_vector) & cluster_vector >= 1]
  if (length(cluster_vector) == 0) return(NA_integer_)
  length(unique(cluster_vector))
}
choose_representative_replicate <- function(som_models, som_clusters) {
  number_of_replicates <- length(som_clusters)
  if (number_of_replicates == 1) return(1L)
  sample_cluster_assignments <- vector("list", number_of_replicates)
  for (replicate_index in seq_len(number_of_replicates)) {
    unit_classif <- as.integer(som_models[[replicate_index]]$unit.classif)
    unit_cluster_labels <- as.integer(som_clusters[[replicate_index]])
    sample_cluster_assignments[[replicate_index]] <- unit_cluster_labels[unit_classif]
  }
  k_values <- vapply(som_clusters, count_SOM_clusters, integer(1))
  k_frequency <- table(k_values[!is.na(k_values)])
  modal_k_values <- as.integer(names(k_frequency)[k_frequency == max(k_frequency)])
  selected_k <- min(modal_k_values)
  candidate_replicates <- which(k_values == selected_k)
  if (length(candidate_replicates) == 1) return(candidate_replicates)
  pairwise_adjusted_rand_index <- matrix(NA_real_, nrow = length(candidate_replicates), ncol = length(candidate_replicates))
  diag(pairwise_adjusted_rand_index) <- NA_real_
  for (candidate_index_1 in seq_len(length(candidate_replicates) - 1)) {
    for (candidate_index_2 in seq.int(candidate_index_1 + 1, length(candidate_replicates))) {
      replicate_index_1 <- candidate_replicates[candidate_index_1]
      replicate_index_2 <- candidate_replicates[candidate_index_2]
      current_adjusted_rand_index <- mclust::adjustedRandIndex(sample_cluster_assignments[[replicate_index_1]], sample_cluster_assignments[[replicate_index_2]])
      pairwise_adjusted_rand_index[candidate_index_1, candidate_index_2] <- current_adjusted_rand_index
      pairwise_adjusted_rand_index[candidate_index_2, candidate_index_1] <- current_adjusted_rand_index
    }
  }
  mean_adjusted_rand_index <- rowMeans(pairwise_adjusted_rand_index, na.rm = TRUE)
  if (all(!is.finite(mean_adjusted_rand_index) | is.na(mean_adjusted_rand_index))) return(candidate_replicates[1])
  representative_candidate_index <- which.max(replace(mean_adjusted_rand_index, !is.finite(mean_adjusted_rand_index) | is.na(mean_adjusted_rand_index), -Inf))
  candidate_replicates[representative_candidate_index]
}
som_models_use <- Pascagoula_SOM$som_models
som_clusters_use <- Pascagoula_SOM$som_clusters
replicate_index <- choose_representative_replicate(som_models = som_models_use, som_clusters = som_clusters_use)
som_model <- som_models_use[[replicate_index]]
som_cluster <- as.integer(som_clusters_use[[replicate_index]])
neighbor_distances <- calc_unit_neighbor_dist(som_model)
som_cluster <- as.integer(factor(som_cluster, levels = sort(unique(som_cluster))))
number_of_clusters <- length(unique(som_cluster))
cluster_colors <- viridis::viridis(number_of_clusters)
SOM_cluster_plot_col <- cluster_colors[som_cluster]
device_width_in <- (plot_width_cm / 2.54) * (96 / 72)
device_height_in <- (plot_height_cm / 2.54) * (96 / 72)
measurement_file <- tempfile(fileext = ".svg")
svg(measurement_file, width = device_width_in, height = device_height_in, family = "Arial")
base_font_size <- par("ps")
legend_relative_font_size <- (legend_font_size * (96 / 72)) / base_font_size
line_height_in <- par("csi") * par("mex")
panel_gap_in <- (panel_gap_cm / 2.54) * (96 / 72)
neighbor_panel_margin_width_in <- (4.6 + 0.6) * line_height_in
cluster_panel_margin_width_in <- (0.6 + 0.6) * line_height_in
map_plot_width_in <- (device_width_in - panel_gap_in - neighbor_panel_margin_width_in - cluster_panel_margin_width_in) / 2
neighbor_panel_width_in <- map_plot_width_in + neighbor_panel_margin_width_in
neighbor_panel_end_initial <- neighbor_panel_width_in / device_width_in
cluster_panel_start_initial <- (neighbor_panel_width_in + panel_gap_in) / device_width_in
par(family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(fig = c(0, neighbor_panel_end_initial, 0, 1), new = FALSE)
plot(x = som_model, type = "property", property = neighbor_distances, main = "", shape = cluster_shape_neighbor_dist, cex = legend_relative_font_size, palette.name = function(n) rev(viridis::cividis(n)))
left_map_right <- grconvertX(max(som_model$grid$pts[, 1]) + 0.5, from = "user", to = "ndc")
par(fig = c(cluster_panel_start_initial, 1, 0, 1), new = TRUE)
plot(x = som_model, shape = cluster_shape_clusters, type = "mapping", bgcol = SOM_cluster_plot_col, main = "", pchs = point_shape_clusters, cex = point_size_clusters, col = point_color_clusters)
right_map_left <- grconvertX(min(som_model$grid$pts[, 1]) - 0.5, from = "user", to = "ndc")
dev.off()
unlink(measurement_file)
actual_gap_fraction <- right_map_left - left_map_right
desired_gap_fraction <- panel_gap_cm / plot_width_cm
panel_shift_fraction <- (actual_gap_fraction - desired_gap_fraction) / 2
neighbor_panel_start <- panel_shift_fraction
neighbor_panel_end <- neighbor_panel_end_initial + panel_shift_fraction
cluster_panel_start <- cluster_panel_start_initial - panel_shift_fraction
cluster_panel_end <- 1 - panel_shift_fraction
svg(file.path(figure_files_folder, figure_name), width = device_width_in, height = device_height_in, family = "Arial")
base_font_size <- par("ps")
legend_relative_font_size <- (legend_font_size * (96 / 72)) / base_font_size
par(family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(fig = c(neighbor_panel_start, neighbor_panel_end, 0, 1), new = FALSE)
plot(x = som_model, type = "property", property = neighbor_distances, main = "", shape = cluster_shape_neighbor_dist, cex = legend_relative_font_size, palette.name = function(n) rev(viridis::cividis(n)))
par(fig = c(cluster_panel_start, cluster_panel_end, 0, 1), new = TRUE)
plot(x = som_model, shape = cluster_shape_clusters, type = "mapping", bgcol = SOM_cluster_plot_col, main = "", pchs = point_shape_clusters, cex = point_size_clusters, col = point_color_clusters)
if (number_of_clusters > 1) kohonen::add.cluster.boundaries(som_model, som_cluster, lwd = boundary_line_width_clusters, col = boundary_color_clusters)
dev.off()


## Supplementary Figure S9C
plot_width_cm <- 5.68
plot_height_cm <- 10.09
panel_heights <- c(1, 1, 1.1)
panel_gap <- 0.85
side_tick_label_gap <- 0.65
bottom_tick_label_gap <- 0.5
figure_name <- "Supplementary_Figure_S9C.svg"

svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72))
max_k <- Pascagoula_SOM$max_k
optim_k_vals <- as.numeric(Pascagoula_SOM$optim_k_vals)
BIC_values <- Pascagoula_SOM$BIC_values[seq_len(max_k), , drop = FALSE]
support_values <- if (!is.null(Pascagoula_SOM$support_values)) Pascagoula_SOM$support_values[seq_len(max_k), , drop = FALSE] else BIC_values
k_colors <- viridis::magma(max_k)
base_font_size <- par("ps")
bottom_numbers_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
layout(matrix(1:3, ncol = 1), heights = panel_heights)
par(bty = "n", oma = c(0, 0, 0, 0), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mar = c(panel_gap, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
finite_k_rows <- apply(support_values, 1, function(x) any(is.finite(x) & !is.na(x)))
plotted_k_values <- seq_len(max_k)[finite_k_rows]
values_for_plot <- t(support_values[finite_k_rows, , drop = FALSE])
boxplot(values_for_plot, at = plotted_k_values, xlim = c(0.5, max_k + 0.5), outline = FALSE, notch = FALSE, axes = FALSE, ylab = "", main = "", whisklty = 1, staplelty = 1, col = k_colors[plotted_k_values])
axis(1, at = seq_len(max_k), labels = FALSE)
y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = y_axis_breaks, labels = round(y_axis_breaks, 1), las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
par(mar = c(panel_gap, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
delta_BIC_matrix <- apply(BIC_values, 2, function(x) { previous_BIC <- x[-length(x)]; current_BIC <- x[-1]; delta_BIC <- previous_BIC - current_BIC; delta_BIC[!is.finite(previous_BIC) | !is.finite(current_BIC)] <- NA_real_; delta_BIC })
if (is.null(dim(delta_BIC_matrix))) delta_BIC_matrix <- matrix(delta_BIC_matrix, ncol = 1)
finite_delta_rows <- apply(delta_BIC_matrix, 1, function(x) any(is.finite(x) & !is.na(x)))
plotted_delta_k_values <- seq.int(2, max_k)[finite_delta_rows]
delta_BIC_for_plot <- t(delta_BIC_matrix[finite_delta_rows, , drop = FALSE])
boxplot(delta_BIC_for_plot, at = plotted_delta_k_values, xlim = c(0.5, max_k + 0.5), outline = FALSE, notch = FALSE, axes = FALSE, ylab = "", main = "", whisklty = 1, staplelty = 1, col = k_colors[plotted_delta_k_values])
axis(1, at = seq_len(max_k), labels = FALSE)
y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = y_axis_breaks, labels = round(y_axis_breaks, 1), las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
par(mar = c(2, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
k_frequency_values <- table(factor(optim_k_vals, levels = seq_len(max_k))) / length(optim_k_vals)
bar_midpoints <- barplot(k_frequency_values, ylim = c(0, 1), col = k_colors, axes = FALSE, axisnames = FALSE, ylab = "", main = "")
axis(1, at = bar_midpoints, labels = seq_len(max_k), mgp = c(3, bottom_tick_label_gap, 0), font.axis = 1, cex.axis = bottom_numbers_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
dev.off()


## Supplementary Figure S9D
plot_width_cm <- 9.1
plot_height_cm <- 7.09
panel_gap <- 1.55
side_tick_label_gap <- 0.65
bottom_tick_label_gap <- 0.75
bottom_margin_mm <- 30
top_margin_mm <- 2
right_margin_mm <- 2
figure_name <- "Supplementary_Figure_S9D.svg"

load(file.path(intermediate_files_folder, "Pascagoula_SOM_lolo.Rdata"))
layer_summary_table <- leave.one.layer.out.results$layer.summary
replicate_matched_results_table <- leave.one.layer.out.results$replicate.matched.results
successful_replicate_matched_results_table <- replicate_matched_results_table[is.na(replicate_matched_results_table$error) | replicate_matched_results_table$error == "", , drop = FALSE]
SOM_layer_names <- as.character(layer_summary_table$layer)
SOM_layer_names_plot <- SOM_layer_names
SOM_layer_names_plot[SOM_layer_names_plot == "Alleles"] <- "SNP"
successful_replicate_matched_results_table$layer <- factor(successful_replicate_matched_results_table$layer, levels = SOM_layer_names)
layer_colors <- setNames(viridis::turbo(length(SOM_layer_names)), SOM_layer_names)
show_assignment_margin_plot <- any(is.finite(successful_replicate_matched_results_table$delta.mean.assignment.margin) & !is.na(successful_replicate_matched_results_table$delta.mean.assignment.margin))
svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72))
base_font_size <- par("ps")
axis_labels_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
bottom_margin_lines <- (bottom_margin_mm / 25.4) / (par("csi") * par("mex"))
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
if (show_assignment_margin_plot) {
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1))
} else {
  layout(matrix(1:2, nrow = 1), widths = c(1, 1))
}
par(bty = "o", oma = c(bottom_margin_lines, 0, top_margin_lines, right_margin_lines), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mar = c(0, 2, 0, panel_gap), mgp = c(3, side_tick_label_gap, 0))
boxplot(absolute.k.deviation ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names_plot, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
for (layer_index in seq_along(SOM_layer_names)) {
  current_layer_name <- SOM_layer_names[layer_index]
  current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "absolute.k.deviation"]
  current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
  if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
}
box()
par(mar = c(0, 2, 0, if (show_assignment_margin_plot) panel_gap else 0), mgp = c(3, side_tick_label_gap, 0))
boxplot(pairwise.coassignment.change ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names_plot, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
for (layer_index in seq_along(SOM_layer_names)) {
  current_layer_name <- SOM_layer_names[layer_index]
  current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "pairwise.coassignment.change"]
  current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
  if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
}
box()
if (show_assignment_margin_plot) {
  par(mar = c(0, 2, 0, 0), mgp = c(3, side_tick_label_gap, 0))
  boxplot(delta.mean.assignment.margin ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
  axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names_plot, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
  axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
  for (layer_index in seq_along(SOM_layer_names)) {
    current_layer_name <- SOM_layer_names[layer_index]
    current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "delta.mean.assignment.margin"]
    current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
    if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
  }
  box()
}
dev.off()


## Supplementary Figure S8E
plot_width_cm <- 10.09
plot_height_cm <- 7.38
row_gap <- 1.45
column_gap <- 1.2
bottom_tick_label_gap <- 0.6
top_margin_mm <- 2
left_margin_mm <- 2
right_margin_mm <- 2
figure_name <- "Supplementary_Figure_S8E.svg"
bars_threshold_N <- 20
importance_threshold <- 0.0001
layer_label_font_size <- 9
axis_ticks_font_size <- 7.1
overwrite_variable_importance <- FALSE
variable_importance_file <- file.path(intermediate_files_folder, "Monticola71_SOM_variable_importance.rds")

matrix_names <- Monticola71_SOM$input_data_names
first_codebook_list <- kohonen::getCodes(Monticola71_SOM$som_models[[1]])
if (!is.list(first_codebook_list)) first_codebook_list <- list(first_codebook_list)
number_of_layers <- length(first_codebook_list)
if (length(matrix_names) != number_of_layers) matrix_names <- paste0("layer", seq_len(number_of_layers))

if (file.exists(variable_importance_file) && !overwrite_variable_importance) {
  Monticola71_SOM_variable_importance <- readRDS(variable_importance_file)
} else {
  for (layer_index in seq_len(number_of_layers)) {
    if (is.null(colnames(first_codebook_list[[layer_index]]))) colnames(first_codebook_list[[layer_index]]) <- paste0("V", seq_len(ncol(first_codebook_list[[layer_index]])))
  }
  retained_replicate_indices <- seq_along(Monticola71_SOM$som_models)
  Monticola71_SOM_variable_importance <- vector("list", number_of_layers)
  names(Monticola71_SOM_variable_importance) <- matrix_names
  for (layer_index in seq_len(number_of_layers)) {
    Monticola71_SOM_variable_importance[[layer_index]] <- matrix(NA_real_, nrow = length(retained_replicate_indices), ncol = ncol(first_codebook_list[[layer_index]]), dimnames = list(paste0("R", retained_replicate_indices), colnames(first_codebook_list[[layer_index]])))
  }
  for (retained_replicate_position in seq_along(retained_replicate_indices)) {
    retained_replicate_index <- retained_replicate_indices[retained_replicate_position]
    som_model <- Monticola71_SOM$som_models[[retained_replicate_index]]
    neuron_cluster_vector <- Monticola71_SOM$som_clusters[[retained_replicate_index]]
    codebook_list <- kohonen::getCodes(som_model)
    if (!is.list(codebook_list)) codebook_list <- list(codebook_list)
    for (layer_index in seq_len(number_of_layers)) {
      if (is.null(colnames(codebook_list[[layer_index]]))) colnames(codebook_list[[layer_index]]) <- paste0("V", seq_len(ncol(codebook_list[[layer_index]])))
      codebook_matrix <- as.matrix(codebook_list[[layer_index]])
      valid_cluster_rows <- is.finite(neuron_cluster_vector) & !is.na(neuron_cluster_vector)
      neuron_sample_counts <- tabulate(som_model$unit.classif, nbins = length(neuron_cluster_vector))
      codebook_matrix <- codebook_matrix[valid_cluster_rows, , drop = FALSE]
      cluster_labels <- neuron_cluster_vector[valid_cluster_rows]
      neuron_sample_counts <- neuron_sample_counts[valid_cluster_rows]
      neuron_weights <- neuron_sample_counts + 1
      Monticola71_SOM_variable_importance[[layer_index]][retained_replicate_position, ] <- apply(codebook_matrix, 2, function(variable_values) {
        valid_rows <- is.finite(variable_values) & !is.na(variable_values) & is.finite(neuron_weights) & !is.na(neuron_weights)
        variable_values <- variable_values[valid_rows]
        variable_cluster_labels <- cluster_labels[valid_rows]
        variable_weights <- neuron_weights[valid_rows]
        if (length(variable_values) < 2 || length(unique(variable_cluster_labels)) < 2 || sum(variable_weights) <= 0) return(NA_real_)
        weighted_variable_mean <- sum(variable_weights * variable_values) / sum(variable_weights)
        total_sum_of_squares <- sum(variable_weights * (variable_values - weighted_variable_mean)^2)
        if (!is.finite(total_sum_of_squares) || total_sum_of_squares <= 0) return(0)
        cluster_means <- tapply(variable_weights * variable_values, variable_cluster_labels, sum) / tapply(variable_weights, variable_cluster_labels, sum)
        cluster_sizes <- tapply(variable_weights, variable_cluster_labels, sum)
        between_cluster_sum_of_squares <- sum(cluster_sizes * (cluster_means - weighted_variable_mean)^2)
        as.numeric(between_cluster_sum_of_squares / total_sum_of_squares)
      })
    }
  }
  saveRDS(Monticola71_SOM_variable_importance, variable_importance_file)
}

svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72), family = "Arial")
base_font_size <- par("ps")
layer_label_relative_font_size <- (layer_label_font_size * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (axis_ticks_font_size * (96 / 72)) / base_font_size
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
left_margin_lines <- (left_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
par(mfrow = if (number_of_layers <= 3) c(1, number_of_layers) else if (number_of_layers == 4) c(2, 2) else if (number_of_layers <= 6) c(2, 3) else if (number_of_layers <= 8) c(2, 4) else if (number_of_layers == 9) c(3, 3) else c(ceiling(number_of_layers / 3), 3), oma = c(0, left_margin_lines, top_margin_lines, right_margin_lines), mar = c(2.2, column_gap / 2, row_gap, column_gap / 2), mgp = c(3, bottom_tick_label_gap, 0), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1, bty = "o", family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black")
layer_colors <- setNames(viridis::turbo(length(matrix_names)), matrix_names)
for (layer_index in seq_along(Monticola71_SOM_variable_importance)) {
  variable_importance_matrix <- Monticola71_SOM_variable_importance[[layer_index]]
  if (is.null(variable_importance_matrix) || nrow(variable_importance_matrix) == 0 || ncol(variable_importance_matrix) == 0) {
    plot.new()
    next
  }
  median_metric_per_variable <- matrixStats::colMedians(variable_importance_matrix, na.rm = TRUE, useNames = TRUE)
  median_metric_per_variable[!is.finite(median_metric_per_variable)] <- NA_real_
  retained_variable_indices <- which(is.finite(median_metric_per_variable) & !is.na(median_metric_per_variable) & median_metric_per_variable > importance_threshold)
  if (length(retained_variable_indices) == 0) {
    plot.new()
    next
  }
  variable_importance_matrix <- variable_importance_matrix[, retained_variable_indices, drop = FALSE]
  retained_variable_medians <- median_metric_per_variable[colnames(variable_importance_matrix)]
  variable_sort_indices <- order(retained_variable_medians, decreasing = FALSE)
  variable_importance_matrix <- variable_importance_matrix[, variable_sort_indices, drop = FALSE]
  number_of_bars <- ncol(variable_importance_matrix)
  boxplot(variable_importance_matrix, horizontal = TRUE, las = 1, notch = FALSE, outline = FALSE, col = rep(layer_colors[matrix_names[layer_index]], number_of_bars), border = "black", axes = FALSE, whisklty = if (number_of_bars > bars_threshold_N) 0 else 1, staplelty = if (number_of_bars > bars_threshold_N) 0 else 1, names = rep("", number_of_bars))
  axis(1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size, col = "black", col.axis = "black", font.axis = 1)
  mtext(matrix_names[layer_index], side = 3, line = 0.3, cex = layer_label_relative_font_size, font = 1, family = "Arial", col = "black")
}
dev.off()

top5_variables_with_values <- lapply(Monticola71_SOM_variable_importance, function(variable_importance_matrix) {
  median_importance <- matrixStats::colMedians(variable_importance_matrix, na.rm = TRUE, useNames = TRUE)
  median_importance <- median_importance[is.finite(median_importance) & !is.na(median_importance) & median_importance > importance_threshold]
  head(sort(median_importance, decreasing = TRUE), 5)
})
top5_variables_with_values


## Supplementary Figure S9F
plot_width_cm <- 8.93
plot_height_cm <- 12.7
pie_size <- 2.2
lon_buffer_range <- 2.1
lat_buffer_range <- 2.48
scale_size <- 0.3
bottom_tick_label_gap <- 0.55
side_tick_label_gap <- 0.8
north_arrow_position <- c(0.93, 0.1)
north_arrow_length <- 0.8
north_arrow_N_position <- 0.35
scale_position <- c(0.08, 0.095)
figure_name <- "Supplementary_Figure_S9F.svg"

Coordinates <- Pascagoula_spatial[, c("Latitude", "Longitude")]
ancestry_matrix <- as.matrix(Pascagoula_SOM$ancestry_matrix)
Coordinates <- as.data.frame(Coordinates, stringsAsFactors = FALSE)
Coordinates$Latitude <- as.numeric(Coordinates$Latitude)
Coordinates$Longitude <- as.numeric(Coordinates$Longitude)
coordinate_sample_names <- rownames(Coordinates)
ancestry_sample_names <- rownames(ancestry_matrix)
matched_sample_names <- intersect(ancestry_sample_names, coordinate_sample_names)
Coordinates <- Coordinates[matched_sample_names, , drop = FALSE]
ancestry_matrix <- ancestry_matrix[matched_sample_names, , drop = FALSE]
rows_with_missing_coordinates <- which(!is.finite(Coordinates$Latitude) | !is.finite(Coordinates$Longitude))
if (length(rows_with_missing_coordinates) > 0) {
  Coordinates <- Coordinates[-rows_with_missing_coordinates, , drop = FALSE]
  ancestry_matrix <- ancestry_matrix[-rows_with_missing_coordinates, , drop = FALSE]
}
ancestry_row_sums <- rowSums(ancestry_matrix, na.rm = TRUE)
ancestry_matrix <- ancestry_matrix / ancestry_row_sums
ancestry_proportions <- as.data.frame(ancestry_matrix)
cluster_colors <- viridis::viridis(ncol(ancestry_matrix))
longitude_minimum <- min(Coordinates$Longitude) - lon_buffer_range
longitude_maximum <- max(Coordinates$Longitude) + lon_buffer_range
latitude_minimum <- min(Coordinates$Latitude) - lat_buffer_range
latitude_maximum <- max(Coordinates$Latitude) + lat_buffer_range
pie_radius <- pie_size * 0.01 * max(longitude_maximum - longitude_minimum, latitude_maximum - latitude_minimum)
add_admixture_pie <- function(longitude, latitude, ancestry_proportions, cluster_colors, x.radius, y.radius, border.color = "black", line.width = 0.8, number.of.points = 80) {
  ancestry_proportions <- as.numeric(ancestry_proportions)
  ancestry_proportions[is.na(ancestry_proportions) | !is.finite(ancestry_proportions) | ancestry_proportions < 0] <- 0
  if (sum(ancestry_proportions) <= 0) return(invisible(NULL))
  ancestry_proportions <- ancestry_proportions / sum(ancestry_proportions)
  slice_start_angles <- c(0, cumsum(ancestry_proportions)[-length(ancestry_proportions)]) * 2 * pi
  slice_end_angles <- cumsum(ancestry_proportions) * 2 * pi
  for (slice_index in seq_along(ancestry_proportions)) {
    if (ancestry_proportions[slice_index] <= 0) next
    slice_angles <- seq(slice_start_angles[slice_index], slice_end_angles[slice_index], length.out = max(3, ceiling(number.of.points * ancestry_proportions[slice_index])))
    polygon(x = c(longitude, longitude + x.radius * cos(slice_angles), longitude), y = c(latitude, latitude + y.radius * sin(slice_angles), latitude), col = cluster_colors[slice_index], border = border.color, lwd = line.width)
  }
  circle_angles <- seq(0, 2 * pi, length.out = number.of.points)
  lines(longitude + x.radius * cos(circle_angles), latitude + y.radius * sin(circle_angles), col = border.color, lwd = line.width)
  invisible(NULL)
}
add_map_legend <- function(legend.position, legend.title, legend.labels, legend.colors, legend.text.relative.font.size, legend.title.relative.font.size, legend.text.font, legend.symbol.size, legend.box) {
  plot_coordinate_limits <- par("usr")
  plot_longitude_range <- plot_coordinate_limits[2] - plot_coordinate_limits[1]
  plot_latitude_range <- plot_coordinate_limits[4] - plot_coordinate_limits[3]
  legend_text_width <- max(strwidth(legend.labels, units = "user", cex = legend.text.relative.font.size, font = legend.text.font, family = "Arial"))
  legend_text_height <- strheight("M", units = "user", cex = legend.text.relative.font.size, font = legend.text.font, family = "Arial")
  legend_title_width <- 0
  legend_title_height <- 0
  legend_title_gap <- 0
  if (!is.null(legend.title)) {
    legend_title_width <- strwidth(legend.title, units = "user", cex = legend.title.relative.font.size, font = 2, family = "Arial")
    legend_title_height <- strheight("M", units = "user", cex = legend.title.relative.font.size, font = 2, family = "Arial")
    legend_title_gap <- 0.5 * legend_text_height
  }
  legend_symbol_width <- strwidth("M", units = "user", cex = legend.text.relative.font.size, family = "Arial") * legend.symbol.size
  legend_symbol_gap <- 0.7 * strwidth("M", units = "user", cex = legend.text.relative.font.size, family = "Arial")
  legend_padding_x <- 0.9 * strwidth("M", units = "user", cex = legend.text.relative.font.size, family = "Arial")
  legend_padding_y <- 0.7 * legend_text_height
  legend_line_gap <- 0.35 * legend_text_height
  legend_width <- max(legend_title_width, legend_symbol_width + legend_symbol_gap + legend_text_width) + 2 * legend_padding_x
  legend_height <- 2 * legend_padding_y + legend_title_height + legend_title_gap + length(legend.labels) * legend_text_height + (length(legend.labels) - 1) * legend_line_gap
  if (legend.position == "topright") {
    legend_left <- plot_coordinate_limits[2] - legend_width
    legend_bottom <- plot_coordinate_limits[4] - legend_height
  } else if (legend.position == "topleft") {
    legend_left <- plot_coordinate_limits[1]
    legend_bottom <- plot_coordinate_limits[4] - legend_height
  } else if (legend.position == "bottomright") {
    legend_left <- plot_coordinate_limits[2] - legend_width
    legend_bottom <- plot_coordinate_limits[3]
  } else if (legend.position == "bottomleft") {
    legend_left <- plot_coordinate_limits[1]
    legend_bottom <- plot_coordinate_limits[3]
  } else if (legend.position == "right") {
    legend_left <- plot_coordinate_limits[2] - legend_width
    legend_bottom <- plot_coordinate_limits[3] + 0.5 * plot_latitude_range - 0.5 * legend_height
  } else if (legend.position == "left") {
    legend_left <- plot_coordinate_limits[1]
    legend_bottom <- plot_coordinate_limits[3] + 0.5 * plot_latitude_range - 0.5 * legend_height
  } else if (legend.position == "top") {
    legend_left <- plot_coordinate_limits[1] + 0.5 * plot_longitude_range - 0.5 * legend_width
    legend_bottom <- plot_coordinate_limits[4] - legend_height
  } else if (legend.position == "bottom") {
    legend_left <- plot_coordinate_limits[1] + 0.5 * plot_longitude_range - 0.5 * legend_width
    legend_bottom <- plot_coordinate_limits[3]
  } else {
    legend_left <- plot_coordinate_limits[1] + 0.5 * plot_longitude_range - 0.5 * legend_width
    legend_bottom <- plot_coordinate_limits[3] + 0.5 * plot_latitude_range - 0.5 * legend_height
  }
  legend_right <- legend_left + legend_width
  legend_top <- legend_bottom + legend_height
  if (legend.box) rect(legend_left, legend_bottom, legend_right, legend_top, col = "white", border = "black")
  current_legend_y_position <- legend_top - legend_padding_y
  if (!is.null(legend.title)) {
    text(x = legend_left + legend_padding_x, y = current_legend_y_position - 0.5 * legend_title_height, labels = legend.title, adj = c(0, 0.5), font = 2, cex = legend.title.relative.font.size, family = "Arial", col = "black")
    current_legend_y_position <- current_legend_y_position - legend_title_height - legend_title_gap
  }
  legend_symbol_x_position <- legend_left + legend_padding_x + 0.5 * legend_symbol_width
  legend_text_x_position <- legend_left + legend_padding_x + legend_symbol_width + legend_symbol_gap
  for (legend_index in seq_along(legend.labels)) {
    legend_entry_y_position <- current_legend_y_position - 0.5 * legend_text_height - (legend_index - 1) * (legend_text_height + legend_line_gap)
    points(x = legend_symbol_x_position, y = legend_entry_y_position, pch = 21, cex = legend.symbol.size, bg = legend.colors[legend_index], col = "black")
    text(x = legend_text_x_position, y = legend_entry_y_position, labels = legend.labels[legend_index], adj = c(0, 0.5), cex = legend.text.relative.font.size, font = legend.text.font, family = "Arial", col = "black")
  }
  invisible(NULL)
}
svg(file.path(figure_files_folder, figure_name), width = plot_width_cm / 2.54, height = plot_height_cm / 2.54, family = "Arial")
base_font_size <- par("ps")
scale_text_relative_font_size <- 7.1 / base_font_size
axis_labels_relative_font_size <- 9.1 / base_font_size
legend_text_relative_font_size <- 9.1 / base_font_size
legend_title_relative_font_size <- 9.1 / base_font_size
par(mfrow = c(1, 1), oma = c(2, 2, 1, 1), mar = c(1, 1, 1, 1), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1, family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", bty = "o")
current_plot_region_size_inches <- par("pin")
longitude_range <- longitude_maximum - longitude_minimum
latitude_range <- latitude_maximum - latitude_minimum
mean_map_latitude <- mean(c(latitude_minimum, latitude_maximum))
longitude_latitude_correction <- cos(mean_map_latitude * pi / 180)
if (!is.finite(longitude_latitude_correction) || longitude_latitude_correction <= 0) longitude_latitude_correction <- 1
target_height_width_ratio <- latitude_range / (longitude_range * longitude_latitude_correction)
adjusted_plot_width_inches <- current_plot_region_size_inches[1]
adjusted_plot_height_inches <- adjusted_plot_width_inches * target_height_width_ratio
if (adjusted_plot_height_inches > current_plot_region_size_inches[2]) {
  adjusted_plot_height_inches <- current_plot_region_size_inches[2]
  adjusted_plot_width_inches <- adjusted_plot_height_inches / target_height_width_ratio
}
par(pin = c(adjusted_plot_width_inches, adjusted_plot_height_inches))
plot.new()
plot.window(xlim = c(longitude_minimum, longitude_maximum), ylim = c(latitude_minimum, latitude_maximum), xaxs = "i", yaxs = "i")
plot_coordinate_limits <- par("usr")
plot_region_size_inches <- par("pin")
x_units_per_inch <- plot_region_size_inches[1] / diff(plot_coordinate_limits[1:2])
y_units_per_inch <- plot_region_size_inches[2] / diff(plot_coordinate_limits[3:4])
if (!is.finite(x_units_per_inch) || !is.finite(y_units_per_inch) || x_units_per_inch <= 0 || y_units_per_inch <= 0) {
  pie_radius_x <- pie_radius
  pie_radius_y <- pie_radius
} else {
  pie_radius_x <- pie_radius * y_units_per_inch / x_units_per_inch
  pie_radius_y <- pie_radius
}
maps::map("world", fill = TRUE, col = "lightgrey", border = "black", xlim = c(longitude_minimum, longitude_maximum), ylim = c(latitude_minimum, latitude_maximum), add = TRUE)
try(maps::map("state", add = TRUE, xlim = c(longitude_minimum, longitude_maximum), ylim = c(latitude_minimum, latitude_maximum), col = "black", lwd = 0.5), silent = TRUE)
longitude_clip_buffer <- 0.2 * (longitude_maximum - longitude_minimum)
latitude_clip_buffer <- 0.2 * (latitude_maximum - latitude_minimum)
clip(longitude_minimum - longitude_clip_buffer, longitude_maximum + longitude_clip_buffer, latitude_minimum - latitude_clip_buffer, latitude_maximum + latitude_clip_buffer)
axis(1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size, col = "black", col.axis = "black", font.axis = 1)
axis(2, las = 2, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size, col = "black", col.axis = "black", font.axis = 1)
box(col = "black")
for (sample_index in seq_len(nrow(ancestry_proportions))) {
  add_admixture_pie(longitude = Coordinates$Longitude[sample_index], latitude = Coordinates$Latitude[sample_index], ancestry_proportions = as.numeric(ancestry_proportions[sample_index, ]), cluster_colors = cluster_colors, x.radius = pie_radius_x, y.radius = pie_radius_y)
}
species_by_sample <- as.character(Pascagoula_data[matched_sample_names, "Species"])
species_names <- unique(species_by_sample)
species_cluster_means <- t(vapply(species_names, function(species_name) colMeans(ancestry_matrix[species_by_sample == species_name, , drop = FALSE], na.rm = TRUE), numeric(ncol(ancestry_matrix))))
cluster_species_names <- species_names[apply(species_cluster_means, 2, which.max)]
cluster_species_names <- paste("D.", gsub("_", " ", cluster_species_names))
legend_labels <- parse(text = paste0("italic('", cluster_species_names, "')"))
add_map_legend(legend.position = "topright", legend.title = "Species", legend.labels = legend_labels, legend.colors = cluster_colors, legend.text.relative.font.size = legend_text_relative_font_size, legend.title.relative.font.size = legend_title_relative_font_size, legend.text.font = 1, legend.symbol.size = 1.6, legend.box = TRUE)
scale_position_longitude <- scale_position[1] * (longitude_maximum - longitude_minimum) + longitude_minimum
scale_position_latitude <- scale_position[2] * (latitude_maximum - latitude_minimum) + latitude_minimum
maps::map.scale(x = scale_position_longitude, y = scale_position_latitude, cex = scale_text_relative_font_size, relwidth = scale_size, ratio = FALSE)
north_arrow_longitude <- north_arrow_position[1] * (longitude_maximum - longitude_minimum) + longitude_minimum
north_arrow_latitude <- north_arrow_position[2] * (latitude_maximum - latitude_minimum) + latitude_minimum
arrows(x0 = north_arrow_longitude, y0 = north_arrow_latitude, x1 = north_arrow_longitude, y1 = north_arrow_latitude + north_arrow_length, length = 0.13, col = "black", lwd = 1.7)
text(x = north_arrow_longitude, y = north_arrow_latitude + north_arrow_length + north_arrow_N_position, labels = "N", cex = 0.8, col = "black", family = "Arial", font = 1)
dev.off()


#### Supplementary Figure S10 ##################################################

## Install packages
rm(list = ls()) #clear environment
CRAN_packages <- c("tibble", "dplyr", "vcfR", "kohonen", "mclust", "viridis", "matrixStats") #CRAN packages
Bioconductor_packages <- "SeqArray" #Bioconductor package
for(p in CRAN_packages) if(!requireNamespace(p, quietly = TRUE)) install.packages(p) #install missing CRAN packages
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") #install BiocManager
if(!requireNamespace(Bioconductor_packages, quietly = TRUE)) BiocManager::install(Bioconductor_packages, ask = FALSE, update = FALSE) #install missing Bioconductor package
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv") #import (NicheDiv::transform.skewed.variables function from my other package (under development)
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")


## Set directories
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
intermediate_files_folder <- "Empirical_examples/Intermediate_files"
figure_files_folder <- "Figure_files"


## Import and prepare data
Viburnum_SNP <- process.SNP.data.SOM(vcf.path = "./Empirical_examples/Spriggs_et_al_2018/nudum-c88-d6-min50.vcf.gz",
                                     missing.loci.cutoff.lenient = 0.7,
                                     missing.loci.cutoff.final = 0.5,
                                     missing.individuals.cutoff = 0.5)
Viburnum_morphology <- read.delim("./Empirical_examples/Spriggs_et_al_2018/morphological_trait_data2.txt", stringsAsFactors = FALSE)
Viburnum_morphology <- Viburnum_morphology[!duplicated(Viburnum_morphology$Individual), ] #remove duplicate IDs
rownames(Viburnum_morphology) <- Viburnum_morphology$Individual #add rownames
Viburnum_morphology$Individual <- NULL
Viburnum_morphology$State <- NULL
Viburnum_morphology$County <- NULL
Viburnum_morphology <- dplyr::rename(Viburnum_morphology, #rename variables
                                     Petiole_length = petiole_length,
                                     Leaf_Area = Area,
                                     Leaf_perimeter = Perimeter,
                                     Leaf_major_ellipse_axis = Major,
                                     Leaf_minor_ellipse_axis = Minor,
                                     Leaf_circularity = Circ.,
                                     Leaf_aspect_ratio = AR,
                                     Leaf_glossy_surface = leaves_lustrous,
                                     Leaf_dark_coloration = drys_dark,
                                     Leaf_margin_type = Leaf_margin,
                                     Mean_peduncle_length = mean_peduncle,
                                     Mean_cyme_length = mean_cyme,
                                     Blade_petiole_length_difference = length_difference)
Viburnum_log_transform_traits <- c(
  "Petiole_length",
  "Leaf_Area",
  "Leaf_perimeter",
  "Leaf_major_ellipse_axis",
  "Leaf_minor_ellipse_axis",
  "Mean_peduncle_length",
  "Mean_cyme_length")
Viburnum_morphology[, Viburnum_log_transform_traits] <- lapply(Viburnum_morphology[, Viburnum_log_transform_traits, drop = FALSE], function(trait_values) { #log transform positive size traits
  trait_values <- as.numeric(trait_values)
  if(any(trait_values <= 0, na.rm = TRUE)) trait_values
  else log(trait_values)
}) #log transform positive size traits
Viburnum_morphology <- remove.lowCV.multicollinearity.SOM(Viburnum_morphology, #remove highly correlated and low-variance variables
                                                          CV.threshold = 0.05,
                                                          cor.threshold = 0.9,
                                                          exclude.cols = c("Leaf_glossy_surface",
                                                                           "Leaf_dark_coloration",
                                                                           "Leaf_margin_type"))
Viburnum_metadata <- read.delim("./Empirical_examples/Spriggs_et_al_2018/morphological_trait_data2.txt", stringsAsFactors = FALSE)
Viburnum_metadata <- Viburnum_metadata[!duplicated(Viburnum_metadata$Individual), ] #remove duplicate IDs
rownames(Viburnum_metadata) <- Viburnum_metadata$Individual #add rownames
Viburnum_metadata <- Viburnum_metadata[, c("State", "County"), drop = FALSE] #only keep State and County columns
nrow(Viburnum_metadata) #number of samples: 145
Viburnum_sample_names <- rownames(Viburnum_SNP)
extract.sample.id.and.species <- function(sample_name) {
  if (grepl("^carlesii_", sample_name)) {
    species_name <- "carlesii"
    sample_id <- sub("^carlesii_", "", sample_name)
  } else if (grepl("^rhytidophyllum_", sample_name)) {
    species_name <- "rhytidophyllum"
    sample_id <- sub("^rhytidophyllum_", "", sample_name)
  } else {
    parts <- strsplit(sample_name, "_")[[1]]
    if (length(parts) > 2 && grepl("^[a-z]+$", parts[length(parts) - 1]) && grepl("^[A-Z]+$", parts[length(parts)])) {
      species_name <- paste(parts[(length(parts) - 1):length(parts)], collapse = "_")
      sample_id <- paste(parts[1:(length(parts) - 2)], collapse = "_")
    } else {
      species_name <- parts[length(parts)]
      sample_id <- paste(parts[1:(length(parts) - 1)], collapse = "_")}}
  data.frame(Sample_ID = sample_id, Species_Name = species_name, stringsAsFactors = FALSE)}
Viburnum_id_species_df <- do.call(rbind, lapply(Viburnum_sample_names, extract.sample.id.and.species))
rownames(Viburnum_SNP) <- Viburnum_id_species_df$Sample_ID
Viburnum_SNP <- Viburnum_SNP[Viburnum_id_species_df$Species_Name == "nudum", , drop = FALSE] #remove outgroup samples from SNP dataset
Viburnum_shared_ids <- Reduce(intersect, list(rownames(Viburnum_SNP), rownames(Viburnum_morphology), rownames(Viburnum_metadata)))
Viburnum_SNP <- Viburnum_SNP[Viburnum_shared_ids, , drop = FALSE]
Viburnum_morphology <- Viburnum_morphology[Viburnum_shared_ids, , drop = FALSE]
Viburnum_metadata <- Viburnum_metadata[Viburnum_shared_ids, , drop = FALSE]
Viburnum_SOM_data <- list(Morphology = Viburnum_morphology,
                          SNP = Viburnum_SNP)
Viburnum_SOM_tr <- train.SOM(Viburnum_SOM_data, #46 samples, 8.4min
                             max.NA.row = 0.5,
                             max.NA.col = 0.5,
                             save.SOM.results.name = file.path(intermediate_files_folder, "Viburnum_SOM_tr.Rdata"),
                             save.SOM.results = TRUE)
Viburnum_SOM <- clustering.SOM(Viburnum_SOM_tr, #21.7min
                               clustering.method = "kmeans+BICelbow",
                               save.SOM.results.name = file.path(intermediate_files_folder, "Viburnum_SOM_kmeansBICelbow.Rdata"))
Viburnum_SOM$optim_k_summary #k2 100%


## Supplementary Figure S10A & Supplementary Figure S10B
plot_width_cm <- 12.49
plot_height_cm <- 4
panel_gap_cm <- 0.52
boundary_color_clusters <- "red"
boundary_line_width_clusters <- 3
point_color_clusters <- "white"
point_shape_clusters <- 19
point_size_clusters <- 0.9
cluster_shape_clusters <- "straight"
cluster_shape_neighbor_dist <- "straight"
legend_font_size <- 7
figure_name <- "Supplementary_Figure_S10AB.svg"

calc_unit_neighbor_dist <- function(som_model) {
  number_of_units <- nrow(som_model$grid$pts)
  codebook_distance_matrix <- as.matrix(kohonen::object.distances(som_model, type = "codes"))
  grid_distance_matrix <- as.matrix(kohonen::unit.distances(som_model$grid))
  neighbor_matrix <- abs(grid_distance_matrix - 1) <= 0.001
  codebook_distance_matrix[!neighbor_matrix] <- NA_real_
  unit_mean_neighbor_distances <- colMeans(codebook_distance_matrix, na.rm = TRUE)
  unit_mean_neighbor_distances[!is.finite(unit_mean_neighbor_distances)] <- NA_real_
  unit_mean_neighbor_distances
}
count_SOM_clusters <- function(cluster_vector) {
  cluster_vector <- as.integer(cluster_vector)
  cluster_vector <- cluster_vector[is.finite(cluster_vector) & !is.na(cluster_vector) & cluster_vector >= 1]
  if (length(cluster_vector) == 0) return(NA_integer_)
  length(unique(cluster_vector))
}
choose_representative_replicate <- function(som_models, som_clusters) {
  number_of_replicates <- length(som_clusters)
  if (number_of_replicates == 1) return(1L)
  sample_cluster_assignments <- vector("list", number_of_replicates)
  for (replicate_index in seq_len(number_of_replicates)) {
    unit_classif <- as.integer(som_models[[replicate_index]]$unit.classif)
    unit_cluster_labels <- as.integer(som_clusters[[replicate_index]])
    sample_cluster_assignments[[replicate_index]] <- unit_cluster_labels[unit_classif]
  }
  k_values <- vapply(som_clusters, count_SOM_clusters, integer(1))
  k_frequency <- table(k_values[!is.na(k_values)])
  modal_k_values <- as.integer(names(k_frequency)[k_frequency == max(k_frequency)])
  selected_k <- min(modal_k_values)
  candidate_replicates <- which(k_values == selected_k)
  if (length(candidate_replicates) == 1) return(candidate_replicates)
  pairwise_adjusted_rand_index <- matrix(NA_real_, nrow = length(candidate_replicates), ncol = length(candidate_replicates))
  diag(pairwise_adjusted_rand_index) <- NA_real_
  for (candidate_index_1 in seq_len(length(candidate_replicates) - 1)) {
    for (candidate_index_2 in seq.int(candidate_index_1 + 1, length(candidate_replicates))) {
      replicate_index_1 <- candidate_replicates[candidate_index_1]
      replicate_index_2 <- candidate_replicates[candidate_index_2]
      current_adjusted_rand_index <- mclust::adjustedRandIndex(sample_cluster_assignments[[replicate_index_1]], sample_cluster_assignments[[replicate_index_2]])
      pairwise_adjusted_rand_index[candidate_index_1, candidate_index_2] <- current_adjusted_rand_index
      pairwise_adjusted_rand_index[candidate_index_2, candidate_index_1] <- current_adjusted_rand_index
    }
  }
  mean_adjusted_rand_index <- rowMeans(pairwise_adjusted_rand_index, na.rm = TRUE)
  if (all(!is.finite(mean_adjusted_rand_index) | is.na(mean_adjusted_rand_index))) return(candidate_replicates[1])
  representative_candidate_index <- which.max(replace(mean_adjusted_rand_index, !is.finite(mean_adjusted_rand_index) | is.na(mean_adjusted_rand_index), -Inf))
  candidate_replicates[representative_candidate_index]
}
som_models_use <- Viburnum_SOM$som_models
som_clusters_use <- Viburnum_SOM$som_clusters
replicate_index <- choose_representative_replicate(som_models = som_models_use, som_clusters = som_clusters_use)
som_model <- som_models_use[[replicate_index]]
som_cluster <- as.integer(som_clusters_use[[replicate_index]])
neighbor_distances <- calc_unit_neighbor_dist(som_model)
som_cluster <- as.integer(factor(som_cluster, levels = sort(unique(som_cluster))))
number_of_clusters <- length(unique(som_cluster))
cluster_colors <- viridis::viridis(number_of_clusters)
SOM_cluster_plot_col <- cluster_colors[som_cluster]
device_width_in <- (plot_width_cm / 2.54) * (96 / 72)
device_height_in <- (plot_height_cm / 2.54) * (96 / 72)
measurement_file <- tempfile(fileext = ".svg")
svg(measurement_file, width = device_width_in, height = device_height_in, family = "Arial")
base_font_size <- par("ps")
legend_relative_font_size <- (legend_font_size * (96 / 72)) / base_font_size
line_height_in <- par("csi") * par("mex")
panel_gap_in <- (panel_gap_cm / 2.54) * (96 / 72)
neighbor_panel_margin_width_in <- (4.6 + 0.6) * line_height_in
cluster_panel_margin_width_in <- (0.6 + 0.6) * line_height_in
map_plot_width_in <- (device_width_in - panel_gap_in - neighbor_panel_margin_width_in - cluster_panel_margin_width_in) / 2
neighbor_panel_width_in <- map_plot_width_in + neighbor_panel_margin_width_in
neighbor_panel_end_initial <- neighbor_panel_width_in / device_width_in
cluster_panel_start_initial <- (neighbor_panel_width_in + panel_gap_in) / device_width_in
par(family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(fig = c(0, neighbor_panel_end_initial, 0, 1), new = FALSE)
plot(x = som_model, type = "property", property = neighbor_distances, main = "", shape = cluster_shape_neighbor_dist, cex = legend_relative_font_size, palette.name = function(n) rev(viridis::cividis(n)))
left_map_right <- grconvertX(max(som_model$grid$pts[, 1]) + 0.5, from = "user", to = "ndc")
par(fig = c(cluster_panel_start_initial, 1, 0, 1), new = TRUE)
plot(x = som_model, shape = cluster_shape_clusters, type = "mapping", bgcol = SOM_cluster_plot_col, main = "", pchs = point_shape_clusters, cex = point_size_clusters, col = point_color_clusters)
right_map_left <- grconvertX(min(som_model$grid$pts[, 1]) - 0.5, from = "user", to = "ndc")
dev.off()
unlink(measurement_file)
actual_gap_fraction <- right_map_left - left_map_right
desired_gap_fraction <- panel_gap_cm / plot_width_cm
panel_shift_fraction <- (actual_gap_fraction - desired_gap_fraction) / 2
neighbor_panel_start <- panel_shift_fraction
neighbor_panel_end <- neighbor_panel_end_initial + panel_shift_fraction
cluster_panel_start <- cluster_panel_start_initial - panel_shift_fraction
cluster_panel_end <- 1 - panel_shift_fraction
svg(file.path(figure_files_folder, figure_name), width = device_width_in, height = device_height_in, family = "Arial")
base_font_size <- par("ps")
legend_relative_font_size <- (legend_font_size * (96 / 72)) / base_font_size
par(family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(fig = c(neighbor_panel_start, neighbor_panel_end, 0, 1), new = FALSE)
plot(x = som_model, type = "property", property = neighbor_distances, main = "", shape = cluster_shape_neighbor_dist, cex = legend_relative_font_size, palette.name = function(n) rev(viridis::cividis(n)))
par(fig = c(cluster_panel_start, cluster_panel_end, 0, 1), new = TRUE)
plot(x = som_model, shape = cluster_shape_clusters, type = "mapping", bgcol = SOM_cluster_plot_col, main = "", pchs = point_shape_clusters, cex = point_size_clusters, col = point_color_clusters)
if (number_of_clusters > 1) kohonen::add.cluster.boundaries(som_model, som_cluster, lwd = boundary_line_width_clusters, col = boundary_color_clusters)
dev.off()


## Supplementary Figure S10C
plot_width_cm <- 5.68
plot_height_cm <- 8.98
panel_heights <- c(1, 1, 1.1)
panel_gap <- 0.85
side_tick_label_gap <- 0.65
bottom_tick_label_gap <- 0.5
figure_name <- "Supplementary_Figure_S10C.svg"

svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72))
max_k <- Viburnum_SOM$max_k
optim_k_vals <- as.numeric(Viburnum_SOM$optim_k_vals)
BIC_values <- Viburnum_SOM$BIC_values[seq_len(max_k), , drop = FALSE]
support_values <- if (!is.null(Viburnum_SOM$support_values)) Viburnum_SOM$support_values[seq_len(max_k), , drop = FALSE] else BIC_values
k_colors <- viridis::magma(max_k)
base_font_size <- par("ps")
bottom_numbers_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
layout(matrix(1:3, ncol = 1), heights = panel_heights)
par(bty = "n", oma = c(0, 0, 0, 0), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mar = c(panel_gap, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
finite_k_rows <- apply(support_values, 1, function(x) any(is.finite(x) & !is.na(x)))
plotted_k_values <- seq_len(max_k)[finite_k_rows]
values_for_plot <- t(support_values[finite_k_rows, , drop = FALSE])
boxplot(values_for_plot, at = plotted_k_values, xlim = c(0.5, max_k + 0.5), outline = FALSE, notch = FALSE, axes = FALSE, ylab = "", main = "", whisklty = 1, staplelty = 1, col = k_colors[plotted_k_values])
axis(1, at = seq_len(max_k), labels = FALSE)
y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = y_axis_breaks, labels = round(y_axis_breaks, 1), las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
par(mar = c(panel_gap, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
delta_BIC_matrix <- apply(BIC_values, 2, function(x) { previous_BIC <- x[-length(x)]; current_BIC <- x[-1]; delta_BIC <- previous_BIC - current_BIC; delta_BIC[!is.finite(previous_BIC) | !is.finite(current_BIC)] <- NA_real_; delta_BIC })
if (is.null(dim(delta_BIC_matrix))) delta_BIC_matrix <- matrix(delta_BIC_matrix, ncol = 1)
finite_delta_rows <- apply(delta_BIC_matrix, 1, function(x) any(is.finite(x) & !is.na(x)))
plotted_delta_k_values <- seq.int(2, max_k)[finite_delta_rows]
delta_BIC_for_plot <- t(delta_BIC_matrix[finite_delta_rows, , drop = FALSE])
boxplot(delta_BIC_for_plot, at = plotted_delta_k_values, xlim = c(0.5, max_k + 0.5), outline = FALSE, notch = FALSE, axes = FALSE, ylab = "", main = "", whisklty = 1, staplelty = 1, col = k_colors[plotted_delta_k_values])
axis(1, at = seq_len(max_k), labels = FALSE)
y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = y_axis_breaks, labels = round(y_axis_breaks, 1), las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
par(mar = c(2, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
k_frequency_values <- table(factor(optim_k_vals, levels = seq_len(max_k))) / length(optim_k_vals)
bar_midpoints <- barplot(k_frequency_values, ylim = c(0, 1), col = k_colors, axes = FALSE, axisnames = FALSE, ylab = "", main = "")
axis(1, at = bar_midpoints, labels = seq_len(max_k), mgp = c(3, bottom_tick_label_gap, 0), font.axis = 1, cex.axis = bottom_numbers_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
dev.off()


## Supplementary Figure S10D
plot_width_cm <- 9.1
plot_height_cm <- 5.7
panel_gap <- 1.55
side_tick_label_gap <- 0.65
bottom_tick_label_gap <- 0.75
bottom_margin_mm <- 23
top_margin_mm <- 2
right_margin_mm <- 2
figure_name <- "Supplementary_Figure_S10D.svg"
leave_one_layer_out_file <- file.path(intermediate_files_folder, "Viburnum_SOM_lolo.Rdata")

if (!file.exists(leave_one_layer_out_file)) {
  plot.layer.importance.leaveoneout.SOM(Viburnum_SOM,
                                        save.leave.one.layer.out.results = TRUE,
                                        save.leave.one.layer.out.results.name = leave_one_layer_out_file)
}
load(leave_one_layer_out_file)
layer_summary_table <- leave.one.layer.out.results$layer.summary
replicate_matched_results_table <- leave.one.layer.out.results$replicate.matched.results
successful_replicate_matched_results_table <- replicate_matched_results_table[is.na(replicate_matched_results_table$error) | replicate_matched_results_table$error == "", , drop = FALSE]
SOM_layer_names <- as.character(layer_summary_table$layer)
successful_replicate_matched_results_table$layer <- factor(successful_replicate_matched_results_table$layer, levels = SOM_layer_names)
layer_colors <- setNames(viridis::turbo(length(SOM_layer_names)), SOM_layer_names)
show_assignment_margin_plot <- any(is.finite(successful_replicate_matched_results_table$delta.mean.assignment.margin) & !is.na(successful_replicate_matched_results_table$delta.mean.assignment.margin))
svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72))
base_font_size <- par("ps")
axis_labels_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
bottom_margin_lines <- (bottom_margin_mm / 25.4) / (par("csi") * par("mex"))
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
if (show_assignment_margin_plot) {
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1))
} else {
  layout(matrix(1:2, nrow = 1), widths = c(1, 1))
}
par(bty = "o", oma = c(bottom_margin_lines, 0, top_margin_lines, right_margin_lines), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mar = c(0, 2, 0, panel_gap), mgp = c(3, side_tick_label_gap, 0))
boxplot(absolute.k.deviation ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
for (layer_index in seq_along(SOM_layer_names)) {
  current_layer_name <- SOM_layer_names[layer_index]
  current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "absolute.k.deviation"]
  current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
  if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
}
box()
par(mar = c(0, 2, 0, if (show_assignment_margin_plot) panel_gap else 0), mgp = c(3, side_tick_label_gap, 0))
boxplot(pairwise.coassignment.change ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
for (layer_index in seq_along(SOM_layer_names)) {
  current_layer_name <- SOM_layer_names[layer_index]
  current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "pairwise.coassignment.change"]
  current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
  if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
}
box()
if (show_assignment_margin_plot) {
  par(mar = c(0, 2, 0, 0), mgp = c(3, side_tick_label_gap, 0))
  boxplot(delta.mean.assignment.margin ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
  axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_names, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
  axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
  for (layer_index in seq_along(SOM_layer_names)) {
    current_layer_name <- SOM_layer_names[layer_index]
    current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "delta.mean.assignment.margin"]
    current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
    if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
  }
  box()
}
dev.off()


## Supplementary Figure S10E
plot_width_cm <- 9.48
plot_height_cm <- 5.4
axis.ticks.font.size <- 7.1
side_tick_label_gap <- 0.6
axis_bar_gap <- 1.82
font_family <- "Arial"
font_color <- "black"
plot_file_name <- "Supplementary_Figure_S10E.svg"

Viburnum_ancestry_plot <- as.matrix(Viburnum_SOM$ancestry_matrix)
Viburnum_dominant_cluster <- max.col(Viburnum_ancestry_plot, ties.method = "first")
Viburnum_assignment_strength <- apply(Viburnum_ancestry_plot, 1, max)
Viburnum_order <- order(Viburnum_dominant_cluster, -Viburnum_assignment_strength)
Viburnum_ancestry_plot <- Viburnum_ancestry_plot[Viburnum_order, , drop = FALSE]
svg_scaling_factor <- 96 / 72
cluster_colors <- viridis::viridis(ncol(Viburnum_ancestry_plot))
plotting_assignment_coefficients <- apply(cbind(0, Viburnum_ancestry_plot), 1, cumsum)
svg(file.path(figure_files_folder, plot_file_name),
    width = (plot_width_cm / 2.54) * svg_scaling_factor,
    height = (plot_height_cm / 2.54) * svg_scaling_factor,
    family = font_family)
base_font_size <- par("ps")
axis_ticks_relative_font_size <- (axis.ticks.font.size * svg_scaling_factor) / base_font_size
par(mar = c(2, 2, 1.5, 1.5),
    family = font_family,
    fg = font_color,
    col.axis = font_color,
    col.lab = font_color,
    col.main = font_color,
    bty = "n",
    cex = 1,
    cex.axis = 1,
    cex.lab = 1,
    cex.main = 1)
plot(0,
     xlim = c(-axis_bar_gap, nrow(Viburnum_ancestry_plot)),
     ylim = c(0, 1),
     type = "n",
     ylab = "",
     xlab = "",
     xaxt = "n",
     yaxt = "n",
     xaxs = "i",
     yaxs = "i",
     frame.plot = FALSE)
axis(side = 2,
     las = 3,
     mgp = c(3, side_tick_label_gap, 0),
     col = font_color,
     col.axis = font_color,
     cex.axis = axis_ticks_relative_font_size)
for (cluster_index in seq_len(ncol(Viburnum_ancestry_plot))) {
  for (individual_index in seq_len(nrow(Viburnum_ancestry_plot))) {
    polygon(x = c(individual_index - 1, individual_index, individual_index, individual_index - 1),
            y = c(plotting_assignment_coefficients[cluster_index, individual_index],
                  plotting_assignment_coefficients[cluster_index, individual_index],
                  plotting_assignment_coefficients[cluster_index + 1, individual_index],
                  plotting_assignment_coefficients[cluster_index + 1, individual_index]),
            col = cluster_colors[cluster_index],
            border = cluster_colors[cluster_index],
            lwd = 0.5)
  }
}
dev.off()


## Supplementary Figure S10F
plot_width_cm <- 6.27
plot_height_cm <- 8.79
row_gap <- 1.45
column_gap <- 1.45
bottom_tick_label_gap <- 0.6
top_margin_mm <- 2
left_margin_mm <- 2
right_margin_mm <- 2
figure_name <- "Supplementary_Figure_S10F.svg"
bars_threshold_N <- 20
importance_threshold <- 0.0001
layer_label_font_size <- 9
axis_ticks_font_size <- 7.1
overwrite_variable_importance <- FALSE
variable_importance_file <- file.path(intermediate_files_folder, "Viburnum_SOM_variable_importance.rds")

matrix_names <- Viburnum_SOM$input_data_names
first_codebook_list <- kohonen::getCodes(Viburnum_SOM$som_models[[1]])
if (!is.list(first_codebook_list)) first_codebook_list <- list(first_codebook_list)
number_of_layers <- length(first_codebook_list)
if (length(matrix_names) != number_of_layers) matrix_names <- paste0("layer", seq_len(number_of_layers))

if (file.exists(variable_importance_file) && !overwrite_variable_importance) {
  Viburnum_SOM_variable_importance <- readRDS(variable_importance_file)
} else {
  for (layer_index in seq_len(number_of_layers)) {
    if (is.null(colnames(first_codebook_list[[layer_index]]))) colnames(first_codebook_list[[layer_index]]) <- paste0("V", seq_len(ncol(first_codebook_list[[layer_index]])))
  }
  retained_replicate_indices <- seq_along(Viburnum_SOM$som_models)
  Viburnum_SOM_variable_importance <- vector("list", number_of_layers)
  names(Viburnum_SOM_variable_importance) <- matrix_names
  for (layer_index in seq_len(number_of_layers)) {
    Viburnum_SOM_variable_importance[[layer_index]] <- matrix(NA_real_, nrow = length(retained_replicate_indices), ncol = ncol(first_codebook_list[[layer_index]]), dimnames = list(paste0("R", retained_replicate_indices), colnames(first_codebook_list[[layer_index]])))
  }
  for (retained_replicate_position in seq_along(retained_replicate_indices)) {
    retained_replicate_index <- retained_replicate_indices[retained_replicate_position]
    som_model <- Viburnum_SOM$som_models[[retained_replicate_index]]
    neuron_cluster_vector <- Viburnum_SOM$som_clusters[[retained_replicate_index]]
    codebook_list <- kohonen::getCodes(som_model)
    if (!is.list(codebook_list)) codebook_list <- list(codebook_list)
    for (layer_index in seq_len(number_of_layers)) {
      if (is.null(colnames(codebook_list[[layer_index]]))) colnames(codebook_list[[layer_index]]) <- paste0("V", seq_len(ncol(codebook_list[[layer_index]])))
      codebook_matrix <- as.matrix(codebook_list[[layer_index]])
      valid_cluster_rows <- is.finite(neuron_cluster_vector) & !is.na(neuron_cluster_vector)
      neuron_sample_counts <- tabulate(som_model$unit.classif, nbins = length(neuron_cluster_vector))
      codebook_matrix <- codebook_matrix[valid_cluster_rows, , drop = FALSE]
      cluster_labels <- neuron_cluster_vector[valid_cluster_rows]
      neuron_sample_counts <- neuron_sample_counts[valid_cluster_rows]
      neuron_weights <- neuron_sample_counts + 1
      Viburnum_SOM_variable_importance[[layer_index]][retained_replicate_position, ] <- apply(codebook_matrix, 2, function(variable_values) {
        valid_rows <- is.finite(variable_values) & !is.na(variable_values) & is.finite(neuron_weights) & !is.na(neuron_weights)
        variable_values <- variable_values[valid_rows]
        variable_cluster_labels <- cluster_labels[valid_rows]
        variable_weights <- neuron_weights[valid_rows]
        if (length(variable_values) < 2 || length(unique(variable_cluster_labels)) < 2 || sum(variable_weights) <= 0) return(NA_real_)
        weighted_variable_mean <- sum(variable_weights * variable_values) / sum(variable_weights)
        total_sum_of_squares <- sum(variable_weights * (variable_values - weighted_variable_mean)^2)
        if (!is.finite(total_sum_of_squares) || total_sum_of_squares <= 0) return(0)
        cluster_means <- tapply(variable_weights * variable_values, variable_cluster_labels, sum) / tapply(variable_weights, variable_cluster_labels, sum)
        cluster_sizes <- tapply(variable_weights, variable_cluster_labels, sum)
        between_cluster_sum_of_squares <- sum(cluster_sizes * (cluster_means - weighted_variable_mean)^2)
        as.numeric(between_cluster_sum_of_squares / total_sum_of_squares)
      })
    }
  }
  saveRDS(Viburnum_SOM_variable_importance, variable_importance_file)
}

svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72), family = "Arial")
base_font_size <- par("ps")
layer_label_relative_font_size <- (layer_label_font_size * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (axis_ticks_font_size * (96 / 72)) / base_font_size
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
left_margin_lines <- (left_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
panel_layout <- if (number_of_layers == 2) c(2, 1) else if (number_of_layers <= 3) c(1, number_of_layers) else if (number_of_layers == 4) c(2, 2) else if (number_of_layers <= 6) c(2, 3) else if (number_of_layers <= 8) c(2, 4) else if (number_of_layers == 9) c(3, 3) else c(ceiling(number_of_layers / 3), 3)
par(mfrow = panel_layout, oma = c(0, left_margin_lines, top_margin_lines, right_margin_lines), mar = c(2.2, column_gap / 2, row_gap, column_gap / 2), mgp = c(3, bottom_tick_label_gap, 0), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1, bty = "o", family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black")
layer_colors <- setNames(viridis::turbo(length(matrix_names)), matrix_names)
for (layer_index in seq_along(Viburnum_SOM_variable_importance)) {
  variable_importance_matrix <- Viburnum_SOM_variable_importance[[layer_index]]
  if (is.null(variable_importance_matrix) || nrow(variable_importance_matrix) == 0 || ncol(variable_importance_matrix) == 0) {
    plot.new()
    next
  }
  median_metric_per_variable <- matrixStats::colMedians(variable_importance_matrix, na.rm = TRUE, useNames = TRUE)
  median_metric_per_variable[!is.finite(median_metric_per_variable)] <- NA_real_
  retained_variable_indices <- which(is.finite(median_metric_per_variable) & !is.na(median_metric_per_variable) & median_metric_per_variable > importance_threshold)
  if (length(retained_variable_indices) == 0) {
    plot.new()
    next
  }
  variable_importance_matrix <- variable_importance_matrix[, retained_variable_indices, drop = FALSE]
  retained_variable_medians <- median_metric_per_variable[colnames(variable_importance_matrix)]
  variable_sort_indices <- order(retained_variable_medians, decreasing = FALSE)
  variable_importance_matrix <- variable_importance_matrix[, variable_sort_indices, drop = FALSE]
  number_of_bars <- ncol(variable_importance_matrix)
  boxplot(variable_importance_matrix, horizontal = TRUE, las = 1, notch = FALSE, outline = FALSE, col = rep(layer_colors[matrix_names[layer_index]], number_of_bars), border = "black", axes = FALSE, whisklty = if (number_of_bars > bars_threshold_N) 0 else 1, staplelty = if (number_of_bars > bars_threshold_N) 0 else 1, names = rep("", number_of_bars))
  axis(1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size, col = "black", col.axis = "black", font.axis = 1)
  mtext(matrix_names[layer_index], side = 3, line = 0.3, cex = layer_label_relative_font_size, font = 1, family = "Arial", col = "black")
}
dev.off()

top5_variables_with_values <- lapply(Viburnum_SOM_variable_importance, function(variable_importance_matrix) {
  median_importance <- matrixStats::colMedians(variable_importance_matrix, na.rm = TRUE, useNames = TRUE)
  median_importance <- median_importance[is.finite(median_importance) & !is.na(median_importance) & median_importance > importance_threshold]
  head(sort(median_importance, decreasing = TRUE), 5)
})
top5_variables_with_values




#### Supplementary Figure S11 ##################################################

## Install packages
rm(list = ls()) #clear environment
CRAN_packages <- c("tibble", "dplyr", "vcfR", "kohonen", "mclust", "viridis", "matrixStats") #CRAN packages
Bioconductor_packages <- "SeqArray" #Bioconductor package
for(p in CRAN_packages) if(!requireNamespace(p, quietly = TRUE)) install.packages(p) #install missing CRAN packages
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") #install BiocManager
if(!requireNamespace(Bioconductor_packages, quietly = TRUE)) BiocManager::install(Bioconductor_packages, ask = FALSE, update = FALSE) #install missing Bioconductor package
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv") #import (NicheDiv::transform.skewed.variables function from my other package (under development)
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")


## Set directories
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
intermediate_files_folder <- "Empirical_examples/Intermediate_files"
figure_files_folder <- "Figure_files"


## Import and prepare data
Elysia_COI <- process.SNP.data.SOM(nexus.path = "./Empirical_examples/Krug_et_al_2026/Elysia_mtDNA_expanded.nex",
                                   missing.loci.cutoff.lenient = 0.7,
                                   missing.loci.cutoff.final = 0.5,
                                   missing.individuals.cutoff = 0.5)
Elysia_metadata <- read.csv("./Empirical_examples/Krug_et_al_2026/Elysia_metadata_updated.csv",
                            header = TRUE,
                            stringsAsFactors = FALSE,
                            check.names = FALSE)
rownames(Elysia_metadata) <- Elysia_metadata$New_sample_name
Elysia_spatial <- dplyr::select(Elysia_metadata, Latitude, Longitude, Seadepth)
rownames(Elysia_spatial) <- rownames(Elysia_metadata)
Elysia_environmental <- read.csv("./Empirical_examples/Krug_et_al_2026/Elysia_environmental.csv",
                                 header = TRUE,
                                 stringsAsFactors = FALSE,
                                 row.names = 1)
Elysia_environmental <- (NicheDiv::transform.skewed.variables(Elysia_environmental))$transformed #transform skewed variables
Elysia_environmental <- remove.lowCV.multicollinearity.SOM(Elysia_environmental, #remove highly correlated and low-variance variables
                                                           CV.threshold = 0.05,
                                                           cor.threshold = 0.9)
Elysia_host <- dplyr::select(Elysia_metadata, Host)
rownames(Elysia_host) <- rownames(Elysia_metadata)
Elysia_host <- make.cols.binary.SOM(Elysia_host, #convert Host to binary columns and keep original
                                    make.binary.cols = "Host",
                                    append.to.original = TRUE)
Elysia_development <- dplyr::select(Elysia_metadata, Developmental_mode)
rownames(Elysia_development) <- rownames(Elysia_metadata)
Elysia_development <- make.cols.binary.SOM(Elysia_development, #convert Developmental_mode to binary columns and keep original
                                           make.binary.cols = "Developmental_mode",
                                           append.to.original = TRUE)
Elysia_host_development <- cbind(Elysia_host, Elysia_development)
Elysia_host_development <- Elysia_host_development[, !duplicated(colnames(Elysia_host_development)), drop = FALSE]
Elysia_host_development <- Elysia_host_development[rownames(Elysia_metadata), , drop = FALSE]
Elysia_common_IDs <- Reduce(intersect, list(rownames(Elysia_COI), rownames(Elysia_spatial), rownames(Elysia_environmental), rownames(Elysia_metadata), rownames(Elysia_host_development)))
Elysia_COI <- Elysia_COI[Elysia_common_IDs, , drop = FALSE]
Elysia_spatial <- Elysia_spatial[Elysia_common_IDs, , drop = FALSE]
Elysia_environmental <- Elysia_environmental[Elysia_common_IDs, , drop = FALSE]
Elysia_metadata <- Elysia_metadata[Elysia_common_IDs, , drop = FALSE]
Elysia_host_development <- Elysia_host_development[Elysia_common_IDs, , drop = FALSE]
Elysia_COI <- Elysia_COI[rowSums(!is.na(Elysia_COI)) > 0, , drop = FALSE]
Elysia_spatial <- Elysia_spatial[rowSums(!is.na(Elysia_spatial)) > 0, , drop = FALSE]
Elysia_environmental <- Elysia_environmental[rowSums(!is.na(Elysia_environmental)) > 0, , drop = FALSE]
Elysia_metadata <- Elysia_metadata[rowSums(!is.na(Elysia_metadata)) > 0, , drop = FALSE]
Elysia_host_development <- Elysia_host_development[rowSums(!is.na(Elysia_host_development)) > 0, , drop = FALSE]
Elysia_common_IDs <- Reduce(intersect, list(rownames(Elysia_COI), rownames(Elysia_spatial), rownames(Elysia_environmental), rownames(Elysia_metadata), rownames(Elysia_host_development)))
Elysia_COI <- Elysia_COI[Elysia_common_IDs, , drop = FALSE]
Elysia_spatial <- Elysia_spatial[Elysia_common_IDs, , drop = FALSE]
Elysia_environmental <- Elysia_environmental[Elysia_common_IDs, , drop = FALSE]
Elysia_metadata <- Elysia_metadata[Elysia_common_IDs, , drop = FALSE]
Elysia_host_development <- Elysia_host_development[Elysia_common_IDs, , drop = FALSE]
Elysia_COI <- Elysia_COI[rownames(Elysia_metadata), , drop = FALSE]
Elysia_spatial <- Elysia_spatial[rownames(Elysia_metadata), , drop = FALSE]
Elysia_environmental <- Elysia_environmental[rownames(Elysia_metadata), , drop = FALSE]
Elysia_host_development <- Elysia_host_development[rownames(Elysia_metadata), , drop = FALSE]
nrow(Elysia_metadata) #number of shared samples: 276
length(unique(Elysia_metadata$Species_name)) #number of species: 7
Elysia_SOM_data <- list(mtDNA = Elysia_COI,
                        Host_development = Elysia_host_development,
                        Environmental = Elysia_environmental,
                        Spatial = Elysia_spatial)
Elysia_SOM_tr <- train.SOM(Elysia_SOM_data, #276 samples, 1.7min
                           max.NA.row = 0.5,
                           max.NA.col = 0.5,
                           save.SOM.results = TRUE,
                           save.SOM.results.name = file.path(intermediate_files_folder, "Elysia_SOM_tr.Rdata"))
Elysia_SOM <- clustering.SOM(Elysia_SOM_tr, #1.7min
                             max.k = 10,
                             clustering.method = "kmeans+BICelbow",
                             save.SOM.results.name = file.path(intermediate_files_folder, "Elysia_SOM_kmeansBICelbow.Rdata"))
Elysia_SOM$optim_k_summary #k4 68%, k2 29%


## Supplementary Figure S11A & Supplementary Figure S11B
plot_width_cm <- 15
plot_height_cm <- 4.335
panel_gap_cm <- 0.52
boundary_color_clusters <- "red"
boundary_line_width_clusters <- 3
point_color_clusters <- "white"
point_shape_clusters <- 19
point_size_clusters <- 0.9
cluster_shape_clusters <- "straight"
cluster_shape_neighbor_dist <- "straight"
legend_font_size <- 7
figure_name <- "Supplementary_Figure_S11AB.svg"

calc_unit_neighbor_dist <- function(som_model) {
  number_of_units <- nrow(som_model$grid$pts)
  codebook_distance_matrix <- as.matrix(kohonen::object.distances(som_model, type = "codes"))
  grid_distance_matrix <- as.matrix(kohonen::unit.distances(som_model$grid))
  neighbor_matrix <- abs(grid_distance_matrix - 1) <= 0.001
  codebook_distance_matrix[!neighbor_matrix] <- NA_real_
  unit_mean_neighbor_distances <- colMeans(codebook_distance_matrix, na.rm = TRUE)
  unit_mean_neighbor_distances[!is.finite(unit_mean_neighbor_distances)] <- NA_real_
  unit_mean_neighbor_distances
}
count_SOM_clusters <- function(cluster_vector) {
  cluster_vector <- as.integer(cluster_vector)
  cluster_vector <- cluster_vector[is.finite(cluster_vector) & !is.na(cluster_vector) & cluster_vector >= 1]
  if (length(cluster_vector) == 0) return(NA_integer_)
  length(unique(cluster_vector))
}
choose_representative_replicate <- function(som_models, som_clusters) {
  number_of_replicates <- length(som_clusters)
  if (number_of_replicates == 1) return(1L)
  sample_cluster_assignments <- vector("list", number_of_replicates)
  for (replicate_index in seq_len(number_of_replicates)) {
    unit_classif <- as.integer(som_models[[replicate_index]]$unit.classif)
    unit_cluster_labels <- as.integer(som_clusters[[replicate_index]])
    sample_cluster_assignments[[replicate_index]] <- unit_cluster_labels[unit_classif]
  }
  k_values <- vapply(som_clusters, count_SOM_clusters, integer(1))
  k_frequency <- table(k_values[!is.na(k_values)])
  modal_k_values <- as.integer(names(k_frequency)[k_frequency == max(k_frequency)])
  selected_k <- min(modal_k_values)
  candidate_replicates <- which(k_values == selected_k)
  if (length(candidate_replicates) == 1) return(candidate_replicates)
  pairwise_adjusted_rand_index <- matrix(NA_real_, nrow = length(candidate_replicates), ncol = length(candidate_replicates))
  diag(pairwise_adjusted_rand_index) <- NA_real_
  for (candidate_index_1 in seq_len(length(candidate_replicates) - 1)) {
    for (candidate_index_2 in seq.int(candidate_index_1 + 1, length(candidate_replicates))) {
      replicate_index_1 <- candidate_replicates[candidate_index_1]
      replicate_index_2 <- candidate_replicates[candidate_index_2]
      current_adjusted_rand_index <- mclust::adjustedRandIndex(sample_cluster_assignments[[replicate_index_1]], sample_cluster_assignments[[replicate_index_2]])
      pairwise_adjusted_rand_index[candidate_index_1, candidate_index_2] <- current_adjusted_rand_index
      pairwise_adjusted_rand_index[candidate_index_2, candidate_index_1] <- current_adjusted_rand_index
    }
  }
  mean_adjusted_rand_index <- rowMeans(pairwise_adjusted_rand_index, na.rm = TRUE)
  if (all(!is.finite(mean_adjusted_rand_index) | is.na(mean_adjusted_rand_index))) return(candidate_replicates[1])
  representative_candidate_index <- which.max(replace(mean_adjusted_rand_index, !is.finite(mean_adjusted_rand_index) | is.na(mean_adjusted_rand_index), -Inf))
  candidate_replicates[representative_candidate_index]
}
som_models_use <- Elysia_SOM$som_models
som_clusters_use <- Elysia_SOM$som_clusters
replicate_index <- choose_representative_replicate(som_models = som_models_use, som_clusters = som_clusters_use)
som_model <- som_models_use[[replicate_index]]
som_cluster <- as.integer(som_clusters_use[[replicate_index]])
neighbor_distances <- calc_unit_neighbor_dist(som_model)
som_cluster <- as.integer(factor(som_cluster, levels = sort(unique(som_cluster))))
number_of_clusters <- length(unique(som_cluster))
cluster_colors <- viridis::viridis(number_of_clusters)
S11B_cluster_colors <- cluster_colors
S11B_sample_clusters <- som_cluster[as.integer(som_model$unit.classif)]
names(S11B_sample_clusters) <- rownames(Elysia_SOM$ancestry_matrix)
SOM_cluster_plot_col <- cluster_colors[som_cluster]
device_width_in <- (plot_width_cm / 2.54) * (96 / 72)
device_height_in <- (plot_height_cm / 2.54) * (96 / 72)
measurement_file <- tempfile(fileext = ".svg")
svg(measurement_file, width = device_width_in, height = device_height_in, family = "Arial")
base_font_size <- par("ps")
legend_relative_font_size <- (legend_font_size * (96 / 72)) / base_font_size
line_height_in <- par("csi") * par("mex")
panel_gap_in <- (panel_gap_cm / 2.54) * (96 / 72)
neighbor_panel_margin_width_in <- (4.6 + 0.6) * line_height_in
cluster_panel_margin_width_in <- (0.6 + 0.6) * line_height_in
map_plot_width_in <- (device_width_in - panel_gap_in - neighbor_panel_margin_width_in - cluster_panel_margin_width_in) / 2
neighbor_panel_width_in <- map_plot_width_in + neighbor_panel_margin_width_in
neighbor_panel_end_initial <- neighbor_panel_width_in / device_width_in
cluster_panel_start_initial <- (neighbor_panel_width_in + panel_gap_in) / device_width_in
par(family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(fig = c(0, neighbor_panel_end_initial, 0, 1), new = FALSE)
plot(x = som_model, type = "property", property = neighbor_distances, main = "", shape = cluster_shape_neighbor_dist, cex = legend_relative_font_size, palette.name = function(n) rev(viridis::cividis(n)))
left_map_right <- grconvertX(max(som_model$grid$pts[, 1]) + 0.5, from = "user", to = "ndc")
par(fig = c(cluster_panel_start_initial, 1, 0, 1), new = TRUE)
plot(x = som_model, shape = cluster_shape_clusters, type = "mapping", bgcol = SOM_cluster_plot_col, main = "", pchs = point_shape_clusters, cex = point_size_clusters, col = point_color_clusters)
right_map_left <- grconvertX(min(som_model$grid$pts[, 1]) - 0.5, from = "user", to = "ndc")
dev.off()
unlink(measurement_file)
actual_gap_fraction <- right_map_left - left_map_right
desired_gap_fraction <- panel_gap_cm / plot_width_cm
panel_shift_fraction <- (actual_gap_fraction - desired_gap_fraction) / 2
neighbor_panel_start <- panel_shift_fraction
neighbor_panel_end <- neighbor_panel_end_initial + panel_shift_fraction
cluster_panel_start <- cluster_panel_start_initial - panel_shift_fraction
cluster_panel_end <- 1 - panel_shift_fraction
svg(file.path(figure_files_folder, figure_name), width = device_width_in, height = device_height_in, family = "Arial")
base_font_size <- par("ps")
legend_relative_font_size <- (legend_font_size * (96 / 72)) / base_font_size
par(family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(fig = c(neighbor_panel_start, neighbor_panel_end, 0, 1), new = FALSE)
plot(x = som_model, type = "property", property = neighbor_distances, main = "", shape = cluster_shape_neighbor_dist, cex = legend_relative_font_size, palette.name = function(n) rev(viridis::cividis(n)))
par(fig = c(cluster_panel_start, cluster_panel_end, 0, 1), new = TRUE)
plot(x = som_model, shape = cluster_shape_clusters, type = "mapping", bgcol = SOM_cluster_plot_col, main = "", pchs = point_shape_clusters, cex = point_size_clusters, col = point_color_clusters)
if (number_of_clusters > 1) kohonen::add.cluster.boundaries(som_model, som_cluster, lwd = boundary_line_width_clusters, col = boundary_color_clusters)
dev.off()


## Supplementary Figure S11C
plot_width_cm <- 5.68
plot_height_cm <- 9.33
panel_heights <- c(1, 1, 1.1)
panel_gap <- 0.85
side_tick_label_gap <- 0.65
bottom_tick_label_gap <- 0.5
figure_name <- "Supplementary_Figure_S11C.svg"

svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72))
max_k <- Elysia_SOM$max_k
optim_k_vals <- as.numeric(Elysia_SOM$optim_k_vals)
BIC_values <- Elysia_SOM$BIC_values[seq_len(max_k), , drop = FALSE]
support_values <- if (!is.null(Elysia_SOM$support_values)) Elysia_SOM$support_values[seq_len(max_k), , drop = FALSE] else BIC_values
k_colors <- viridis::magma(max_k)
base_font_size <- par("ps")
bottom_numbers_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
layout(matrix(1:3, ncol = 1), heights = panel_heights)
par(bty = "n", oma = c(0, 0, 0, 0), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mar = c(panel_gap, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
finite_k_rows <- apply(support_values, 1, function(x) any(is.finite(x) & !is.na(x)))
plotted_k_values <- seq_len(max_k)[finite_k_rows]
values_for_plot <- t(support_values[finite_k_rows, , drop = FALSE])
boxplot(values_for_plot, at = plotted_k_values, xlim = c(0.5, max_k + 0.5), outline = FALSE, notch = FALSE, axes = FALSE, ylab = "", main = "", whisklty = 1, staplelty = 1, col = k_colors[plotted_k_values])
axis(1, at = seq_len(max_k), labels = FALSE)
y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = y_axis_breaks, labels = round(y_axis_breaks, 1), las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
par(mar = c(panel_gap, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
delta_BIC_matrix <- apply(BIC_values, 2, function(x) { previous_BIC <- x[-length(x)]; current_BIC <- x[-1]; delta_BIC <- previous_BIC - current_BIC; delta_BIC[!is.finite(previous_BIC) | !is.finite(current_BIC)] <- NA_real_; delta_BIC })
if (is.null(dim(delta_BIC_matrix))) delta_BIC_matrix <- matrix(delta_BIC_matrix, ncol = 1)
finite_delta_rows <- apply(delta_BIC_matrix, 1, function(x) any(is.finite(x) & !is.na(x)))
plotted_delta_k_values <- seq.int(2, max_k)[finite_delta_rows]
delta_BIC_for_plot <- t(delta_BIC_matrix[finite_delta_rows, , drop = FALSE])
boxplot(delta_BIC_for_plot, at = plotted_delta_k_values, xlim = c(0.5, max_k + 0.5), outline = FALSE, notch = FALSE, axes = FALSE, ylab = "", main = "", whisklty = 1, staplelty = 1, col = k_colors[plotted_delta_k_values])
axis(1, at = seq_len(max_k), labels = FALSE)
y_axis_breaks <- seq(par("usr")[3], par("usr")[4], length.out = 4)
axis(2, at = y_axis_breaks, labels = round(y_axis_breaks, 1), las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
par(mar = c(2, 2, panel_gap, 1), mgp = c(3, side_tick_label_gap, 0))
k_frequency_values <- table(factor(optim_k_vals, levels = seq_len(max_k))) / length(optim_k_vals)
bar_midpoints <- barplot(k_frequency_values, ylim = c(0, 1), col = k_colors, axes = FALSE, axisnames = FALSE, ylab = "", main = "")
axis(1, at = bar_midpoints, labels = seq_len(max_k), mgp = c(3, bottom_tick_label_gap, 0), font.axis = 1, cex.axis = bottom_numbers_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
dev.off()


## Supplementary Figure S11D
plot_width_cm <- 9.1
plot_height_cm <- 5.4
panel_gap <- 1.55
side_tick_label_gap <- 0.65
bottom_tick_label_gap <- 0.75
bottom_margin_mm <- 34
top_margin_mm <- 2
right_margin_mm <- 2
figure_name <- "Supplementary_Figure_S11D.svg"
leave_one_layer_out_file <- file.path(intermediate_files_folder, "Elysia_SOM_lolo.Rdata")

if (!file.exists(leave_one_layer_out_file)) {
  plot.layer.importance.leaveoneout.SOM(Elysia_SOM,
                                        save.leave.one.layer.out.results = TRUE,
                                        save.leave.one.layer.out.results.name = leave_one_layer_out_file)
}
load(leave_one_layer_out_file)
layer_summary_table <- leave.one.layer.out.results$layer.summary
replicate_matched_results_table <- leave.one.layer.out.results$replicate.matched.results
successful_replicate_matched_results_table <- replicate_matched_results_table[is.na(replicate_matched_results_table$error) | replicate_matched_results_table$error == "", , drop = FALSE]
SOM_layer_names <- as.character(layer_summary_table$layer)
SOM_layer_labels <- SOM_layer_names
SOM_layer_labels[SOM_layer_labels == "Host_development"] <- "Host+development"
successful_replicate_matched_results_table$layer <- factor(successful_replicate_matched_results_table$layer, levels = SOM_layer_names)
layer_colors <- setNames(viridis::turbo(length(SOM_layer_names)), SOM_layer_names)
show_assignment_margin_plot <- any(is.finite(successful_replicate_matched_results_table$delta.mean.assignment.margin) & !is.na(successful_replicate_matched_results_table$delta.mean.assignment.margin))
svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72))
base_font_size <- par("ps")
axis_labels_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (7.1 * (96 / 72)) / base_font_size
bottom_margin_lines <- (bottom_margin_mm / 25.4) / (par("csi") * par("mex"))
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
if (show_assignment_margin_plot) {
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1))
} else {
  layout(matrix(1:2, nrow = 1), widths = c(1, 1))
}
par(bty = "o", oma = c(bottom_margin_lines, 0, top_margin_lines, right_margin_lines), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mar = c(0, 2, 0, panel_gap), mgp = c(3, side_tick_label_gap, 0))
boxplot(absolute.k.deviation ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_labels, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
for (layer_index in seq_along(SOM_layer_names)) {
  current_layer_name <- SOM_layer_names[layer_index]
  current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "absolute.k.deviation"]
  current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
  if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
}
box()
par(mar = c(0, 2, 0, if (show_assignment_margin_plot) panel_gap else 0), mgp = c(3, side_tick_label_gap, 0))
boxplot(pairwise.coassignment.change ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_labels, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
for (layer_index in seq_along(SOM_layer_names)) {
  current_layer_name <- SOM_layer_names[layer_index]
  current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "pairwise.coassignment.change"]
  current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
  if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
}
box()
if (show_assignment_margin_plot) {
  par(mar = c(0, 2, 0, 0), mgp = c(3, side_tick_label_gap, 0))
  boxplot(delta.mean.assignment.margin ~ layer, data = successful_replicate_matched_results_table, col = layer_colors[SOM_layer_names], outline = FALSE, axes = FALSE, ylab = "", xlab = "", main = "", whisklty = 1, staplelty = 1)
  axis(1, at = seq_along(SOM_layer_names), labels = SOM_layer_labels, las = 2, font.axis = 1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_labels_relative_font_size)
  axis(2, las = 3, mgp = c(3, side_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size)
  for (layer_index in seq_along(SOM_layer_names)) {
    current_layer_name <- SOM_layer_names[layer_index]
    current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, "delta.mean.assignment.margin"]
    current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
    if (length(current_layer_values) > 0) points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.12), current_layer_values, pch = 16, cex = 0.8, col = adjustcolor(layer_colors[current_layer_name], alpha.f = 0.65))
  }
  box()
}
dev.off()


## Supplementary Figure S11E
plot_width_cm <- 9.48
plot_height_cm <- 5.4
axis.ticks.font.size <- 7.1
side_tick_label_gap <- 0.6
axis_bar_gap <- 1.82
font_family <- "Arial"
font_color <- "black"
plot_file_name <- "Supplementary_Figure_S11E.svg"

Elysia_SOM_k4 <- clustering.SOM(Elysia_SOM_tr,
                                set.k = 4,
                                clustering.method = "kmeans+BICelbow",
                                save.SOM.results.name = file.path(intermediate_files_folder, "Elysia_SOM_kmeansBICelbow_k4.Rdata"))
Elysia_ancestry_plot <- as.matrix(Elysia_SOM_k4$ancestry_matrix)
if (ncol(Elysia_ancestry_plot) != length(S11B_cluster_colors)) stop("Number of clusters differs between S11B and S11E")
Elysia_k4_dominant_cluster <- max.col(Elysia_ancestry_plot, ties.method = "first")
names(Elysia_k4_dominant_cluster) <- rownames(Elysia_ancestry_plot)
shared_samples <- intersect(names(S11B_sample_clusters), names(Elysia_k4_dominant_cluster))
if (length(shared_samples) == 0) stop("No shared samples between S11B and S11E")
reference_clusters <- S11B_sample_clusters[shared_samples]
k4_clusters <- Elysia_k4_dominant_cluster[shared_samples]
number_of_clusters <- ncol(Elysia_ancestry_plot)
cluster_permutations <- expand.grid(rep(list(seq_len(number_of_clusters)), number_of_clusters))
cluster_permutations <- as.matrix(cluster_permutations[apply(cluster_permutations, 1, function(x) length(unique(x)) == number_of_clusters), , drop = FALSE])
permutation_scores <- apply(cluster_permutations, 1, function(cluster_mapping) sum(reference_clusters == cluster_mapping[k4_clusters], na.rm = TRUE))
best_cluster_mapping <- as.integer(cluster_permutations[which.max(permutation_scores), ])
Elysia_ancestry_plot <- Elysia_ancestry_plot[, order(best_cluster_mapping), drop = FALSE]
Elysia_dominant_cluster <- max.col(Elysia_ancestry_plot, ties.method = "first")
Elysia_assignment_strength <- apply(Elysia_ancestry_plot, 1, max)
Elysia_order <- order(Elysia_dominant_cluster, -Elysia_assignment_strength)
Elysia_ancestry_plot <- Elysia_ancestry_plot[Elysia_order, , drop = FALSE]
svg_scaling_factor <- 96 / 72
cluster_colors <- S11B_cluster_colors
plotting_assignment_coefficients <- apply(cbind(0, Elysia_ancestry_plot), 1, cumsum)
svg(file.path(figure_files_folder, plot_file_name),
    width = (plot_width_cm / 2.54) * svg_scaling_factor,
    height = (plot_height_cm / 2.54) * svg_scaling_factor,
    family = font_family)
base_font_size <- par("ps")
axis_ticks_relative_font_size <- (axis.ticks.font.size * svg_scaling_factor) / base_font_size
par(mar = c(2, 2, 1.5, 1.5),
    family = font_family,
    fg = font_color,
    col.axis = font_color,
    col.lab = font_color,
    col.main = font_color,
    bty = "n",
    cex = 1,
    cex.axis = 1,
    cex.lab = 1,
    cex.main = 1)
plot(0,
     xlim = c(-axis_bar_gap, nrow(Elysia_ancestry_plot)),
     ylim = c(0, 1),
     type = "n",
     ylab = "",
     xlab = "",
     xaxt = "n",
     yaxt = "n",
     xaxs = "i",
     yaxs = "i",
     frame.plot = FALSE)
axis(side = 2,
     las = 3,
     mgp = c(3, side_tick_label_gap, 0),
     col = font_color,
     col.axis = font_color,
     cex.axis = axis_ticks_relative_font_size)
for (cluster_index in seq_len(ncol(Elysia_ancestry_plot))) {
  for (individual_index in seq_len(nrow(Elysia_ancestry_plot))) {
    polygon(x = c(individual_index - 1, individual_index, individual_index, individual_index - 1),
            y = c(plotting_assignment_coefficients[cluster_index, individual_index],
                  plotting_assignment_coefficients[cluster_index, individual_index],
                  plotting_assignment_coefficients[cluster_index + 1, individual_index],
                  plotting_assignment_coefficients[cluster_index + 1, individual_index]),
            col = cluster_colors[cluster_index],
            border = cluster_colors[cluster_index],
            lwd = 0.5)
  }
}
dev.off()


## Supplementary Figure S11E
plot_width_cm <- 15.85
plot_height_cm <- 3
axis.ticks.font.size <- 7.1
side_tick_label_gap <- 0.6
axis_bar_gap <- 5.43
font_family <- "Arial"
font_color <- "black"
plot_file_name <- "Supplementary_Figure_S11E.svg"

Elysia_SOM_k4 <- clustering.SOM(Elysia_SOM_tr,
                                set.k = 4,
                                clustering.method = "kmeans+BICelbow",
                                save.SOM.results.name = file.path(intermediate_files_folder, "Elysia_SOM_kmeansBICelbow_k4.Rdata"))
Elysia_ancestry_plot <- as.matrix(Elysia_SOM_k4$ancestry_matrix)
if (ncol(Elysia_ancestry_plot) != length(S11B_cluster_colors)) stop("Number of clusters differs between S11B and S11E")
Elysia_k4_dominant_cluster <- max.col(Elysia_ancestry_plot, ties.method = "first")
names(Elysia_k4_dominant_cluster) <- rownames(Elysia_ancestry_plot)
shared_samples <- intersect(names(S11B_sample_clusters), names(Elysia_k4_dominant_cluster))
if (length(shared_samples) == 0) stop("No shared samples between S11B and S11E")
reference_clusters <- S11B_sample_clusters[shared_samples]
k4_clusters <- Elysia_k4_dominant_cluster[shared_samples]
number_of_clusters <- ncol(Elysia_ancestry_plot)
cluster_permutations <- expand.grid(rep(list(seq_len(number_of_clusters)), number_of_clusters))
cluster_permutations <- as.matrix(cluster_permutations[apply(cluster_permutations, 1, function(x) length(unique(x)) == number_of_clusters), , drop = FALSE])
permutation_scores <- apply(cluster_permutations, 1, function(cluster_mapping) sum(reference_clusters == cluster_mapping[k4_clusters], na.rm = TRUE))
best_cluster_mapping <- as.integer(cluster_permutations[which.max(permutation_scores), ])
Elysia_ancestry_plot <- Elysia_ancestry_plot[, order(best_cluster_mapping), drop = FALSE]
Elysia_dominant_cluster <- max.col(Elysia_ancestry_plot, ties.method = "first")
Elysia_assignment_strength <- apply(Elysia_ancestry_plot, 1, max)
Elysia_order <- order(Elysia_dominant_cluster, -Elysia_assignment_strength)
Elysia_species_ordering <- as.character(Elysia_metadata[rownames(Elysia_ancestry_plot), "Species_name"])
Elysia_species_ordering_clean <- gsub("_", " ", Elysia_species_ordering)
Elysia_zuleicae_group <- grepl("zuleicae", Elysia_species_ordering_clean, ignore.case = TRUE)
Elysia_cf_zuleicae <- Elysia_zuleicae_group & grepl("\\bcf\\.?\\b", Elysia_species_ordering_clean, ignore.case = TRUE)
Elysia_zuleicae <- Elysia_zuleicae_group & !Elysia_cf_zuleicae
if (any(Elysia_zuleicae) && any(Elysia_cf_zuleicae)) {
  Elysia_first_zuleicae_position <- min(which(Elysia_zuleicae_group[Elysia_order]))
  Elysia_zuleicae_indices <- Elysia_order[Elysia_zuleicae[Elysia_order]]
  Elysia_cf_zuleicae_indices <- Elysia_order[Elysia_cf_zuleicae[Elysia_order]]
  Elysia_zuleicae_block <- c(Elysia_zuleicae_indices, Elysia_cf_zuleicae_indices)
  Elysia_order <- Elysia_order[!Elysia_zuleicae_group[Elysia_order]]
  Elysia_order <- append(Elysia_order, Elysia_zuleicae_block, after = Elysia_first_zuleicae_position - 1)
}
Elysia_species_block_order <- unique(Elysia_species_ordering_clean[Elysia_order])
Elysia_patina_position <- which(Elysia_species_block_order == "Elysia patina")
Elysia_christinae_position <- which(Elysia_species_block_order == "Elysia christinae")
if (length(Elysia_patina_position) == 1 && length(Elysia_christinae_position) == 1) {
  Elysia_species_block_order[c(Elysia_patina_position, Elysia_christinae_position)] <- Elysia_species_block_order[c(Elysia_christinae_position, Elysia_patina_position)]
  Elysia_order <- unlist(lapply(Elysia_species_block_order, function(current_species) Elysia_order[Elysia_species_ordering_clean[Elysia_order] == current_species]), use.names = FALSE)
}
Elysia_ancestry_plot <- Elysia_ancestry_plot[Elysia_order, , drop = FALSE]
svg_scaling_factor <- 96 / 72
cluster_colors <- S11B_cluster_colors
plotting_assignment_coefficients <- apply(cbind(0, Elysia_ancestry_plot), 1, cumsum)
svg(file.path(figure_files_folder, plot_file_name),
    width = (plot_width_cm / 2.54) * svg_scaling_factor,
    height = (plot_height_cm / 2.54) * svg_scaling_factor,
    family = font_family)
base_font_size <- par("ps")
axis_ticks_relative_font_size <- (axis.ticks.font.size * svg_scaling_factor) / base_font_size
par(mar = c(2, 2, 1.5, 1.5),
    family = font_family,
    fg = font_color,
    col.axis = font_color,
    col.lab = font_color,
    col.main = font_color,
    bty = "n",
    cex = 1,
    cex.axis = 1,
    cex.lab = 1,
    cex.main = 1)
plot(0,
     xlim = c(-axis_bar_gap, nrow(Elysia_ancestry_plot)),
     ylim = c(0, 1),
     type = "n",
     ylab = "",
     xlab = "",
     xaxt = "n",
     yaxt = "n",
     xaxs = "i",
     yaxs = "i",
     frame.plot = FALSE)
axis(side = 2,
     las = 3,
     mgp = c(3, side_tick_label_gap, 0),
     col = font_color,
     col.axis = font_color,
     cex.axis = axis_ticks_relative_font_size)
for (cluster_index in seq_len(ncol(Elysia_ancestry_plot))) {
  for (individual_index in seq_len(nrow(Elysia_ancestry_plot))) {
    polygon(x = c(individual_index - 1, individual_index, individual_index, individual_index - 1),
            y = c(plotting_assignment_coefficients[cluster_index, individual_index],
                  plotting_assignment_coefficients[cluster_index, individual_index],
                  plotting_assignment_coefficients[cluster_index + 1, individual_index],
                  plotting_assignment_coefficients[cluster_index + 1, individual_index]),
            col = cluster_colors[cluster_index],
            border = cluster_colors[cluster_index],
            lwd = 0.5)
  }
}
dev.off()


## Supplementary Figure S11E_species
plot_width_cm <- 15.85
plot_height_cm <- 3
axis.ticks.font.size <- 7.1
side_tick_label_gap <- 0.6
axis_bar_gap <- 5.43
font_family <- "Arial"
font_color <- "black"
plot_file_name <- "Supplementary_Figure_S11E_species.svg"

Elysia_SOM_k4 <- clustering.SOM(Elysia_SOM_tr,
                                set.k = 4,
                                clustering.method = "kmeans+BICelbow",
                                save.SOM.results.name = file.path(intermediate_files_folder, "Elysia_SOM_kmeansBICelbow_k4.Rdata"))
Elysia_ancestry_plot <- as.matrix(Elysia_SOM_k4$ancestry_matrix)
if (ncol(Elysia_ancestry_plot) != length(S11B_cluster_colors)) stop("Number of clusters differs between S11B and S11E")
Elysia_k4_dominant_cluster <- max.col(Elysia_ancestry_plot, ties.method = "first")
names(Elysia_k4_dominant_cluster) <- rownames(Elysia_ancestry_plot)
shared_samples <- intersect(names(S11B_sample_clusters), names(Elysia_k4_dominant_cluster))
if (length(shared_samples) == 0) stop("No shared samples between S11B and S11E")
reference_clusters <- S11B_sample_clusters[shared_samples]
k4_clusters <- Elysia_k4_dominant_cluster[shared_samples]
number_of_clusters <- ncol(Elysia_ancestry_plot)
cluster_permutations <- expand.grid(rep(list(seq_len(number_of_clusters)), number_of_clusters))
cluster_permutations <- as.matrix(cluster_permutations[apply(cluster_permutations, 1, function(x) length(unique(x)) == number_of_clusters), , drop = FALSE])
permutation_scores <- apply(cluster_permutations, 1, function(cluster_mapping) sum(reference_clusters == cluster_mapping[k4_clusters], na.rm = TRUE))
best_cluster_mapping <- as.integer(cluster_permutations[which.max(permutation_scores), ])
Elysia_ancestry_plot <- Elysia_ancestry_plot[, order(best_cluster_mapping), drop = FALSE]
Elysia_dominant_cluster <- max.col(Elysia_ancestry_plot, ties.method = "first")
Elysia_assignment_strength <- apply(Elysia_ancestry_plot, 1, max)
Elysia_order <- order(Elysia_dominant_cluster, -Elysia_assignment_strength)
Elysia_species_ordering <- as.character(Elysia_metadata[rownames(Elysia_ancestry_plot), "Species_name"])
Elysia_species_ordering_clean <- gsub("_", " ", Elysia_species_ordering)
Elysia_zuleicae_group <- grepl("zuleicae", Elysia_species_ordering_clean, ignore.case = TRUE)
Elysia_cf_zuleicae <- Elysia_zuleicae_group & grepl("\\bcf\\.?\\b", Elysia_species_ordering_clean, ignore.case = TRUE)
Elysia_zuleicae <- Elysia_zuleicae_group & !Elysia_cf_zuleicae
if (any(Elysia_zuleicae) && any(Elysia_cf_zuleicae)) {
  Elysia_first_zuleicae_position <- min(which(Elysia_zuleicae_group[Elysia_order]))
  Elysia_zuleicae_indices <- Elysia_order[Elysia_zuleicae[Elysia_order]]
  Elysia_cf_zuleicae_indices <- Elysia_order[Elysia_cf_zuleicae[Elysia_order]]
  Elysia_zuleicae_block <- c(Elysia_zuleicae_indices, Elysia_cf_zuleicae_indices)
  Elysia_order <- Elysia_order[!Elysia_zuleicae_group[Elysia_order]]
  Elysia_order <- append(Elysia_order, Elysia_zuleicae_block, after = Elysia_first_zuleicae_position - 1)
}
species_color_order <- unique(Elysia_species_ordering[Elysia_order])
Elysia_species_block_order <- unique(Elysia_species_ordering_clean[Elysia_order])
Elysia_patina_position <- which(Elysia_species_block_order == "Elysia patina")
Elysia_christinae_position <- which(Elysia_species_block_order == "Elysia christinae")
if (length(Elysia_patina_position) == 1 && length(Elysia_christinae_position) == 1) {
  Elysia_species_block_order[c(Elysia_patina_position, Elysia_christinae_position)] <- Elysia_species_block_order[c(Elysia_christinae_position, Elysia_patina_position)]
  Elysia_order <- unlist(lapply(Elysia_species_block_order, function(current_species) Elysia_order[Elysia_species_ordering_clean[Elysia_order] == current_species]), use.names = FALSE)
}
Elysia_ancestry_plot <- Elysia_ancestry_plot[Elysia_order, , drop = FALSE]

Elysia_species <- as.character(Elysia_metadata[rownames(Elysia_ancestry_plot), "Species_name"])
if (any(is.na(Elysia_species) | Elysia_species == "")) stop("Species assignment missing for one or more Elysia samples")
species_names <- unique(Elysia_species)
species_colors <- setNames(viridis::viridis(length(species_color_order)), species_color_order)
bar_colors <- species_colors[Elysia_species]
species_legend_labels <- sub("^Elysia ", "E. ", gsub("_", " ", species_names))

svg_scaling_factor <- 96 / 72
svg(file.path(figure_files_folder, plot_file_name),
    width = (plot_width_cm / 2.54) * svg_scaling_factor,
    height = (plot_height_cm / 2.54) * svg_scaling_factor,
    family = font_family)
base_font_size <- par("ps")
axis_ticks_relative_font_size <- (axis.ticks.font.size * svg_scaling_factor) / base_font_size
legend_relative_font_size <- (legend_font_size * svg_scaling_factor) / base_font_size
par(mar = c(2, 2, 3, 1.5),
    family = font_family,
    fg = font_color,
    col.axis = font_color,
    col.lab = font_color,
    col.main = font_color,
    bty = "n",
    cex = 1,
    cex.axis = 1,
    cex.lab = 1,
    cex.main = 1)
plot(0,
     xlim = c(-axis_bar_gap, nrow(Elysia_ancestry_plot)),
     ylim = c(0, 1),
     type = "n",
     ylab = "",
     xlab = "",
     xaxt = "n",
     yaxt = "n",
     xaxs = "i",
     yaxs = "i",
     frame.plot = FALSE)
axis(side = 2,
     las = 3,
     mgp = c(3, side_tick_label_gap, 0),
     col = font_color,
     col.axis = font_color,
     cex.axis = axis_ticks_relative_font_size)
for (individual_index in seq_len(nrow(Elysia_ancestry_plot))) {
  polygon(x = c(individual_index - 1, individual_index, individual_index, individual_index - 1),
          y = c(0, 0, 1, 1),
          col = bar_colors[individual_index],
          border = bar_colors[individual_index],
          lwd = 0.5)
}
legend(x = mean(par("usr")[1:2]),
       y = 1.06,
       legend = species_legend_labels,
       fill = species_colors[species_names],
       border = species_colors[species_names],
       bty = "n",
       cex = legend_relative_font_size,
       text.col = font_color,
       ncol = ceiling(length(species_legend_labels) / 2),
       xjust = 0.5,
       yjust = 0,
       xpd = NA)
dev.off()


## Supplementary Figure S11G
plot_width_cm <- 15.79
plot_height_cm <- 3.01
row_gap <- 1.45
column_gap <- 1.45
bottom_tick_label_gap <- 0.6
top_margin_mm <- 2
left_margin_mm <- 2
right_margin_mm <- 2
figure_name <- "Supplementary_Figure_S11G.svg"
bars_threshold_N <- 20
importance_threshold <- 0.0001
layer_label_font_size <- 9
axis_ticks_font_size <- 7.1
overwrite_variable_importance <- FALSE
variable_importance_file <- file.path(intermediate_files_folder, "Elysia_SOM_variable_importance.rds")

matrix_names <- Elysia_SOM$input_data_names
matrix_labels <- matrix_names
matrix_labels[matrix_labels == "Host_development"] <- "Host+Development"
first_codebook_list <- kohonen::getCodes(Elysia_SOM$som_models[[1]])
if (!is.list(first_codebook_list)) first_codebook_list <- list(first_codebook_list)
number_of_layers <- length(first_codebook_list)
if (length(matrix_names) != number_of_layers) {
  matrix_names <- paste0("layer", seq_len(number_of_layers))
  matrix_labels <- matrix_names
}

if (file.exists(variable_importance_file) && !overwrite_variable_importance) {
  Elysia_SOM_variable_importance <- readRDS(variable_importance_file)
} else {
  for (layer_index in seq_len(number_of_layers)) {
    if (is.null(colnames(first_codebook_list[[layer_index]]))) colnames(first_codebook_list[[layer_index]]) <- paste0("V", seq_len(ncol(first_codebook_list[[layer_index]])))
  }
  retained_replicate_indices <- seq_along(Elysia_SOM$som_models)
  Elysia_SOM_variable_importance <- vector("list", number_of_layers)
  names(Elysia_SOM_variable_importance) <- matrix_names
  for (layer_index in seq_len(number_of_layers)) {
    Elysia_SOM_variable_importance[[layer_index]] <- matrix(NA_real_, nrow = length(retained_replicate_indices), ncol = ncol(first_codebook_list[[layer_index]]), dimnames = list(paste0("R", retained_replicate_indices), colnames(first_codebook_list[[layer_index]])))
  }
  for (retained_replicate_position in seq_along(retained_replicate_indices)) {
    retained_replicate_index <- retained_replicate_indices[retained_replicate_position]
    som_model <- Elysia_SOM$som_models[[retained_replicate_index]]
    neuron_cluster_vector <- Elysia_SOM$som_clusters[[retained_replicate_index]]
    codebook_list <- kohonen::getCodes(som_model)
    if (!is.list(codebook_list)) codebook_list <- list(codebook_list)
    for (layer_index in seq_len(number_of_layers)) {
      if (is.null(colnames(codebook_list[[layer_index]]))) colnames(codebook_list[[layer_index]]) <- paste0("V", seq_len(ncol(codebook_list[[layer_index]])))
      codebook_matrix <- as.matrix(codebook_list[[layer_index]])
      valid_cluster_rows <- is.finite(neuron_cluster_vector) & !is.na(neuron_cluster_vector)
      neuron_sample_counts <- tabulate(som_model$unit.classif, nbins = length(neuron_cluster_vector))
      codebook_matrix <- codebook_matrix[valid_cluster_rows, , drop = FALSE]
      cluster_labels <- neuron_cluster_vector[valid_cluster_rows]
      neuron_sample_counts <- neuron_sample_counts[valid_cluster_rows]
      neuron_weights <- neuron_sample_counts + 1
      Elysia_SOM_variable_importance[[layer_index]][retained_replicate_position, ] <- apply(codebook_matrix, 2, function(variable_values) {
        valid_rows <- is.finite(variable_values) & !is.na(variable_values) & is.finite(neuron_weights) & !is.na(neuron_weights)
        variable_values <- variable_values[valid_rows]
        variable_cluster_labels <- cluster_labels[valid_rows]
        variable_weights <- neuron_weights[valid_rows]
        if (length(variable_values) < 2 || length(unique(variable_cluster_labels)) < 2 || sum(variable_weights) <= 0) return(NA_real_)
        weighted_variable_mean <- sum(variable_weights * variable_values) / sum(variable_weights)
        total_sum_of_squares <- sum(variable_weights * (variable_values - weighted_variable_mean)^2)
        if (!is.finite(total_sum_of_squares) || total_sum_of_squares <= 0) return(0)
        cluster_means <- tapply(variable_weights * variable_values, variable_cluster_labels, sum) / tapply(variable_weights, variable_cluster_labels, sum)
        cluster_sizes <- tapply(variable_weights, variable_cluster_labels, sum)
        between_cluster_sum_of_squares <- sum(cluster_sizes * (cluster_means - weighted_variable_mean)^2)
        as.numeric(between_cluster_sum_of_squares / total_sum_of_squares)
      })
    }
  }
  saveRDS(Elysia_SOM_variable_importance, variable_importance_file)
}

svg(file.path(figure_files_folder, figure_name), width = (plot_width_cm / 2.54) * (96 / 72), height = (plot_height_cm / 2.54) * (96 / 72), family = "Arial")
base_font_size <- par("ps")
layer_label_relative_font_size <- (layer_label_font_size * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (axis_ticks_font_size * (96 / 72)) / base_font_size
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
left_margin_lines <- (left_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
panel_layout <- c(1, number_of_layers)
par(mfrow = panel_layout, oma = c(0, left_margin_lines, top_margin_lines, right_margin_lines), mar = c(2.2, column_gap / 2, row_gap, column_gap / 2), mgp = c(3, bottom_tick_label_gap, 0), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1, bty = "o", family = "Arial", fg = "black", col.axis = "black", col.lab = "black", col.main = "black")
layer_colors <- setNames(viridis::turbo(length(matrix_names)), matrix_names)
for (layer_index in seq_along(Elysia_SOM_variable_importance)) {
  variable_importance_matrix <- Elysia_SOM_variable_importance[[layer_index]]
  if (is.null(variable_importance_matrix) || nrow(variable_importance_matrix) == 0 || ncol(variable_importance_matrix) == 0) {
    plot.new()
    next
  }
  median_metric_per_variable <- matrixStats::colMedians(variable_importance_matrix, na.rm = TRUE, useNames = TRUE)
  median_metric_per_variable[!is.finite(median_metric_per_variable)] <- NA_real_
  retained_variable_indices <- which(is.finite(median_metric_per_variable) & !is.na(median_metric_per_variable) & median_metric_per_variable > importance_threshold)
  if (length(retained_variable_indices) == 0) {
    plot.new()
    next
  }
  variable_importance_matrix <- variable_importance_matrix[, retained_variable_indices, drop = FALSE]
  retained_variable_medians <- median_metric_per_variable[colnames(variable_importance_matrix)]
  variable_sort_indices <- order(retained_variable_medians, decreasing = FALSE)
  variable_importance_matrix <- variable_importance_matrix[, variable_sort_indices, drop = FALSE]
  number_of_bars <- ncol(variable_importance_matrix)
  boxplot(variable_importance_matrix, horizontal = TRUE, las = 1, notch = FALSE, outline = FALSE, col = rep(layer_colors[matrix_names[layer_index]], number_of_bars), border = "black", axes = FALSE, whisklty = if (number_of_bars > bars_threshold_N) 0 else 1, staplelty = if (number_of_bars > bars_threshold_N) 0 else 1, names = rep("", number_of_bars))
  axis(1, mgp = c(3, bottom_tick_label_gap, 0), cex.axis = axis_ticks_relative_font_size, col = "black", col.axis = "black", font.axis = 1)
  mtext(matrix_labels[layer_index], side = 3, line = 0.3, cex = layer_label_relative_font_size, font = 1, family = "Arial", col = "black")
}
dev.off()

top5_variables_with_values <- lapply(Elysia_SOM_variable_importance, function(variable_importance_matrix) {
  median_importance <- matrixStats::colMedians(variable_importance_matrix, na.rm = TRUE, useNames = TRUE)
  median_importance <- median_importance[is.finite(median_importance) & !is.na(median_importance) & median_importance > importance_threshold]
  head(sort(median_importance, decreasing = TRUE), 5)
})
top5_variables_with_values



#### Supplementary Figure S14 ##################################################


## Set environment
rm(list = ls()) #clear environment
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")


## Install and load required packages
required_packages <- c("viridis")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
  library(package_name, character.only = TRUE)
}


## Set directories
intermediate_files_folder <- "Empirical_examples/Intermediate_files"
figure_files_folder <- "Figure_files"


## Set plot parameters
plot_width_cm <- 16
plot_height_cm <- 18.79
row_gap <- 1.55
column_gap <- 2.8
bottom_tick_label_gap <- 0.50
side_tick_label_gap <- 0.45
bottom_margin_mm <- 4
top_margin_mm <- 2
left_margin_mm <- 4
right_margin_mm <- 2
figure_name <- "Supplementary_Figure_S14.svg"
panel_title_font_size <- 9
axis_label_font_size <- 7.1
axis_ticks_font_size <- 7.1


## Load SOM outputs
load.SOM.object <- function(file_path) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  load_environment <- new.env()
  loaded_object_names <- load(file_path, envir = load_environment)
  SOM_candidates <- loaded_object_names[vapply(loaded_object_names, function(object_name_candidate) {
    object_candidate <- get(object_name_candidate, envir = load_environment)
    is.list(object_candidate) && !is.null(object_candidate$distance_weights_matrix)
  }, logical(1))]
  if (length(SOM_candidates) == 0) stop("No SOM object with distance_weights_matrix found in ", file_path)
  if (length(SOM_candidates) > 1) stop("More than one SOM object with distance_weights_matrix found in ", file_path, ": ", paste(SOM_candidates, collapse = ", "))
  get(SOM_candidates[1], envir = load_environment)
}
SOM_plot_specs <- list(
  list(panel_title = expression(italic(Polygonia)),
       file_path = file.path(intermediate_files_folder, "Polygonia_SOM_kmeansBICelbow.Rdata")),
  list(panel_title = expression(italic("D. aeneus")),
       file_path = file.path(intermediate_files_folder, "Aeneus_SOM_kmeansBICelbow.Rdata")),
  list(panel_title = expression(italic("D. monticola")~"+"~italic("D. cheaha")),
       file_path = file.path(intermediate_files_folder, "Monticola71_SOM_kmeansBICelbow.Rdata")),
  list(panel_title = expression(italic("D. pascagoula")~"+"~italic("D. valentinei")),
       file_path = file.path(intermediate_files_folder, "Pascagoula_SOM_kmeansBICelbow.Rdata")),
  list(panel_title = expression(italic(Viburnum)),
       file_path = file.path(intermediate_files_folder, "Viburnum_SOM_kmeansBICelbow.Rdata")),
  list(panel_title = expression(italic(Elysia)),
       file_path = file.path(intermediate_files_folder, "Elysia_SOM_kmeansBICelbow.Rdata")),
  list(panel_title = expression(italic(Pocillopora)),
       file_path = file.path(intermediate_files_folder, "Pocillopora_SOM_kmeansBICelbow.Rdata")),
  list(panel_title = expression(italic(Microcebus)),
       file_path = file.path(intermediate_files_folder, "Microcebus_SOM_kmeansBICelbow.Rdata"))
)
SOM_outputs <- lapply(SOM_plot_specs, function(plot_spec) {
  load.SOM.object(file_path = plot_spec$file_path)
})


## Create additional functions
rename.layer.names <- function(layer_names) {
  layer_names <- as.character(layer_names)
  layer_names[layer_names == "Alleles"] <- "SNP"
  layer_names[layer_names == "COI"] <- "mtDNA"
  layer_names[layer_names == "Host_development"] <- "Host+development"
  layer_names[layer_names == "Development_Host"] <- "Host+development"
  layer_names[layer_names == "Morphology_binary"] <- "Morphology 2"
  layer_names[layer_names == "Morphology_2"] <- "Morphology 2"
  layer_names[layer_names == "Symbiosis_haplotypes"] <- "Symbiosis"
  layer_names[layer_names == "Symbionts"] <- "Symbiosis"
  layer_names[layer_names == "Host_haplotypes"] <- "Host haplotypes"
  layer_names
}
get.layer.names <- function(SOM.output) {
  if ("input_data_names" %in% names(SOM.output)) {
    layer_names <- SOM.output$input_data_names
  } else {
    layer_names <- colnames(SOM.output$distance_weights_matrix)
  }
  if (is.null(layer_names) || length(layer_names) != ncol(SOM.output$distance_weights_matrix) || any(is.na(layer_names))) {
    layer_names <- paste0("Layer", seq_len(ncol(SOM.output$distance_weights_matrix)))
  }
  make.unique(rename.layer.names(layer_names))
}
all_layer_names <- unique(unlist(lapply(SOM_outputs, get.layer.names)))
preferred_layer_name_order <- c("SNP",
                                "mtDNA",
                                "Morphology",
                                "Morphology 2",
                                "Environmental",
                                "Spatial",
                                "Watershed",
                                "Host",
                                "Host haplotypes",
                                "Host+development",
                                "Symbiosis",
                                "Microsatellites")
ordered_layer_names <- c(preferred_layer_name_order[preferred_layer_name_order %in% all_layer_names],
                         setdiff(all_layer_names, preferred_layer_name_order))
global_layer_colors <- setNames(viridis::turbo(length(ordered_layer_names)), ordered_layer_names)

plot.layer.distance.scale.panel <- function(SOM.output,
                                            panel_title,
                                            axis_labels_relative_font_size,
                                            axis_ticks_relative_font_size,
                                            panel_title_relative_font_size,
                                            global_layer_colors) {
  if (is.null(SOM.output$distance_weights_matrix)) stop("Plotting aborted: SOM.output does not contain distance_weights_matrix")
  distance_weights_matrix <- SOM.output$distance_weights_matrix
  if (!is.matrix(distance_weights_matrix)) distance_weights_matrix <- as.matrix(distance_weights_matrix)
  if (!is.numeric(distance_weights_matrix)) stop("Plotting aborted: distance_weights_matrix must be numeric")
  if (any(!is.finite(distance_weights_matrix) | is.na(distance_weights_matrix))) stop("Plotting aborted: distance_weights_matrix contains NA or non-finite values")
  if (any(distance_weights_matrix <= 0)) stop("Plotting aborted: distance_weights_matrix must contain only positive values")
  if (ncol(distance_weights_matrix) < 2) stop("Plotting aborted: at least two layers are required for plotting")
  layer_names <- get.layer.names(SOM.output)
  pairwise_distance_matrix <- 1 / distance_weights_matrix
  colnames(pairwise_distance_matrix) <- layer_names
  mean_pairwise_distances <- colMeans(pairwise_distance_matrix, na.rm = TRUE)
  if (any(!is.finite(mean_pairwise_distances) | is.na(mean_pairwise_distances))) stop("Plotting aborted: calculated mean pairwise distances contain NA or non-finite values")
  layer_order <- order(mean_pairwise_distances, decreasing = TRUE)
  ordered_mean_pairwise_distances <- rev(mean_pairwise_distances[layer_order])
  ordered_layer_names <- rev(layer_names[layer_order])
  bar_midpoints <- barplot(height = ordered_mean_pairwise_distances,
                           col = global_layer_colors[ordered_layer_names],
                           border = "black",
                           names.arg = rep("", length(ordered_layer_names)),
                           axes = FALSE,
                           axisnames = FALSE,
                           horiz = TRUE,
                           xlab = "",
                           main = "")
  axis(1,
       mgp = c(3, bottom_tick_label_gap, 0),
       cex.axis = axis_ticks_relative_font_size)
  axis(2,
       at = bar_midpoints,
       labels = ordered_layer_names,
       las = 1,
       tick = FALSE,
       mgp = c(3, side_tick_label_gap, 0),
       cex.axis = axis_labels_relative_font_size,
       font.axis = 1)
  box()
  mtext(panel_title,
        side = 3,
        line = 0.3,
        adj = 0.5,
        padj = 0,
        cex = panel_title_relative_font_size,
        font = 1)
}


## Plot all distance scale plots
if (!dir.exists(figure_files_folder)) dir.create(figure_files_folder, recursive = TRUE)
svg(file.path(figure_files_folder, figure_name),
    width = (plot_width_cm / 2.54) * (96 / 72),
    height = (plot_height_cm / 2.54) * (96 / 72),
    family = "Arial")
base_font_size <- par("ps")
panel_title_relative_font_size <- (panel_title_font_size * (96 / 72)) / base_font_size
axis_labels_relative_font_size <- (axis_label_font_size * (96 / 72)) / base_font_size
axis_ticks_relative_font_size <- (axis_ticks_font_size * (96 / 72)) / base_font_size
bottom_margin_lines <- (bottom_margin_mm / 25.4) / (par("csi") * par("mex"))
top_margin_lines <- (top_margin_mm / 25.4) / (par("csi") * par("mex"))
left_margin_lines <- (left_margin_mm / 25.4) / (par("csi") * par("mex"))
right_margin_lines <- (right_margin_mm / 25.4) / (par("csi") * par("mex"))
par(mfrow = c(4, 2),
    oma = c(bottom_margin_lines, left_margin_lines, top_margin_lines, right_margin_lines),
    mar = c(2.2, 5.5, row_gap, column_gap / 2),
    cex = 1,
    cex.axis = 1,
    cex.lab = 1,
    cex.main = 1,
    bty = "o",
    family = "Arial",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black")
for (plot_index in seq_along(SOM_outputs)) {
  plot.layer.distance.scale.panel(SOM.output = SOM_outputs[[plot_index]],
                                  panel_title = SOM_plot_specs[[plot_index]]$panel_title,
                                  axis_labels_relative_font_size = axis_labels_relative_font_size,
                                  axis_ticks_relative_font_size = axis_ticks_relative_font_size,
                                  panel_title_relative_font_size = panel_title_relative_font_size,
                                  global_layer_colors = global_layer_colors)
}
dev.off()

