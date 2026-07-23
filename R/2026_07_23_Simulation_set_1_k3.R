#### Set environment and install/load packages #################################

## Set environment
rm(list = ls()) #clear environment
#setwd("./")
setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")


## Install and load required R packages
required_packages <- c("clue",
                       "emmeans",
                       "glmmTMB",
                       "kohonen",
                       "MASS",
                       "mclust",
                       "nlme",
                       "sn",
                       "tidyverse")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
  library(package_name, character.only = TRUE) #install missing packages and load all required packages
}




#### Set main simulation parameters ############################################
simulation_k3_dir <- file.path("Simulations", "Simulation_set_1", "k3")
if (!dir.exists(simulation_k3_dir)) dir.create(simulation_k3_dir, recursive = TRUE)

N_individuals <- 60 #number of individuals
N_SNP_loci <- 1000 #number of SNP loci
N_morph_traits <- 15 #number of morphological traits
N_climate_variables <- 30 #number of climatic variables
N_hosts <- 12 #number of host categories

overwrite <- FALSE
N_simulations <- 60
N_steps_SOM <- 100
N_replicates_SOM <- 110
SOM_grid_size <- rep(floor(sqrt(4 * sqrt(N_individuals))), 2)
max_NA_row_SOM <- 1
max_NA_col_SOM <- 1
verbose_SOM <- FALSE
BIC_threshold_SOM <- 6
learning_rate_tuning <- FALSE

SNP_target_Fst_range <- c(0.10, 0.20) #target average Fst range across all loci (overall divergence) - range reflects moderate divergence often observed among conspecific populations (Hasbún et al. 2016 https://doi.org/10.1155/2016/3654093; Haenel et al. 2021 https://doi.org/10.1038/s41467-021-25039-y; Hall 2022 https://doi.org/10.3390/ani12091115)
SNP_differentiated_prop_range <- c(0.15, 0.4) #proportion of loci with between-cluster differentiation
SNP_differentiated_Fst_range <- c(0.85, 1) #Fst range for differentiated loci (D) with very strong structure
SNP_random_prop_range <- c(0.01, 0.06) #proportion of random loci (R) representing drift and sequencing/genotyping noise (Pompanon et al. 2005 https://doi.org/10.1038/nrg1707; Helyar et al. 2011 https://doi.org/10.1111/j.1755-0998.2010.02943.x; Mastretta-Yanes et al. 2015 https://doi.org/10.1111/1755-0998.12291)
morph_trait_distance_range <- c(0.8, 1.8) #derived distance among cluster means for differentiated morphological traits in standardized multivariate trait space (Leinonen et al. 2008 https://doi.org/10.1111/j.1420-9101.2007.01445.x; De Kort et al. 2013 https://doi.org/10.1007/s10682-012-9624-9; Siefert et al. 2015 https://doi.org/10.1111/ele.12508; Westerband et al. 2021 https://doi.org/10.1093/aob/mcab011; Opedal et al. 2023 https://doi.org/10.1073/pnas.2203228120
morph_trait_sd_range <- c(0.6, 1.4) #derived standard deviation range for morphological traits in standardized trait units (Siefert et al. 2015 https://doi.org/10.1111/ele.12508; Westerband et al. 2021 https://doi.org/10.1093/aob/mcab011
climate_variables_distance_range <- c(0.8, 1.8) #derived distance among cluster means for differentiated climate variables in standardized environmental space (Broennimann et al. 2012 https://doi.org/10.1111/j.1466-8238.2011.00698.x; Liu et al. 2020 https://doi.org/10.1073/pnas.2004289117; Bates & Bertelsmeier 2021 https://doi.org/10.1016/j.cub.2021.08.035
climate_variables_sd_range <- c(0.6, 1.4) #derived standard deviation range for climate variables in standardized environmental units (Liu et al. 2020 https://doi.org/10.1073/pnas.2004289117; Carscadden et al. 2020 https://doi.org/10.1086/710388; Bates & Bertelsmeier 2021 https://doi.org/10.1016/j.cub.2021.08.035
host_dominant_prop_range <- c(0.7, 0.96) #proportion of individuals per cluster assigned to that cluster's dominant host (Ramírez-Martínez & Tlapaya-Romero 2023 https://doi.org/10.1016/j.ijppaw.2023.05.001; Feder et al. 1994 https://doi.org/10.1073/pnas.91.17.7990; Wehmeyer et al. 2024 https://doi.org/10.1186/s13071-024-06439-7

missing_bad_individual_prop_mean <- 0.20
missing_bad_individual_prop_CI95 <- c(0.05, 0.35)
missing_bad_variable_prop_mean <- 0.20
missing_bad_variable_prop_CI95 <- c(0.05, 0.35)

missing_data_prop_default <- 0.2




#### Set plotting parameters ###################################################
benchmark_plot_raw_point_alpha <- 0.3
benchmark_plot_raw_point_size <- 2
benchmark_plot_estimate_point_size <- 3.5
benchmark_plot_estimate_linewidth <- 1.5
benchmark_plot_errorbar_linewidth <- 1.5
benchmark_plot_errorbar_width <- 0.3
benchmark_plot_categorical_jitter_width <- 0.03
benchmark_plot_SE_multiplier <- 1
benchmark_plot_base_size <- 12
benchmark_plot_axis_title_size <- 12
benchmark_plot_axis_text_size <- 10
benchmark_plot_x_axis_text_angle <- 45
benchmark_plot_x_axis_text_hjust <- 1
benchmark_plot_color_ARI <- "#abd9e9"
benchmark_plot_color_Acc <- "#74add1"
benchmark_plot_color_K_correct <- "#fdae61"
benchmark_plot_color_K_far_off <- "#d73027"
benchmark_plot_color_QE <- "#66c2a5"
benchmark_plot_color_TE <- "#fc8d62"
benchmark_plot_color_Time <- "#fc9192"




#### Create functions #############################################################

## Create function to simulate data
simulate.data <- function(N_clusters, missing_data_prop = missing_data_prop_default, sim_id = NULL) {
  
  # Basic checks
  stopifnot(is.numeric(N_clusters), length(N_clusters) == 1, !is.na(N_clusters), N_clusters >= 1, N_clusters %% 1 == 0)
  stopifnot(is.numeric(missing_data_prop), length(missing_data_prop) == 1, !is.na(missing_data_prop),  missing_data_prop >= 0, missing_data_prop <= 1)
  stopifnot(N_clusters - 1 <= N_morph_traits) #must have enough morphology traits to place clusters in K-1 dimensions
  stopifnot(N_clusters - 1 <= N_climate_variables) #must have enough climate variables to place clusters in K-1 dimensions
  stopifnot(N_hosts >= N_clusters) #each cluster must have a matching dominant host column
  
  # Create function to place clusters at equal distances in multivariate space using regular simplex
  generate.simplex.coordinates <- function(N.points, point_distance) {
    if (N.points == 1) return(matrix(0, nrow = 1, ncol = 1)) #for single cluster, return origin
    centering_matrix <- diag(N.points) - matrix(1 / N.points, nrow = N.points, ncol = N.points) #create centering matrix
    eigen_decomposition <- eigen(centering_matrix, symmetric = TRUE) #calculate eigenvectors defining centered simplex space
    coordinate_matrix <- eigen_decomposition$vectors[, seq_len(N.points - 1), drop = FALSE] #retain orthonormal basis for K - 1 simplex dimensions
    coordinate_matrix <- coordinate_matrix * point_distance / sqrt(2) #scale coordinates so every pairwise distance equals point_distance
    coordinate_matrix #return simplex coordinates
  }
  
  # Assign individuals to clusters as evenly as possible
  base_samples_per_cluster <- floor(N_individuals / N_clusters)
  remainder_samples <- N_individuals %% N_clusters
  samples_per_cluster <- rep(base_samples_per_cluster, N_clusters)
  if (remainder_samples > 0) samples_per_cluster[1:remainder_samples] <- samples_per_cluster[1:remainder_samples] + 1
  true_cluster_labels <- unlist(mapply(rep, seq_len(N_clusters), samples_per_cluster))
  individual_ids <- paste0("ID", seq_len(N_individuals)) #create individual IDs
  
  # SNP data
  # Simulate SNP genotypes using three locus classes:
  # 1. Differentiated loci (D): strong between-cluster structure with high Fst
  # 2. Neutral loci (N): same allele frequency distribution across all clusters
  # 3. Random loci (R): mostly unstructured loci with random genotype perturbation to mimic noise
  SNP_differentiated_prop <- runif(1, SNP_differentiated_prop_range[1], SNP_differentiated_prop_range[2]) #draw proportion of differentiated loci
  SNP_random_prop <- runif(1, SNP_random_prop_range[1], SNP_random_prop_range[2]) #draw proportion of random loci
  SNP_differentiated_n <- round(N_SNP_loci * SNP_differentiated_prop) #number of differentiated loci
  SNP_random_n <- round(N_SNP_loci * SNP_random_prop) #number of random loci
  SNP_neutral_n <- N_SNP_loci - SNP_differentiated_n - SNP_random_n #number of neutral loci
  SNP_locus_types <- c(rep("D", SNP_differentiated_n), rep("N", SNP_neutral_n), rep("R", SNP_random_n)) #assign locus classes
  SNP_locus_types <- sample(SNP_locus_types) #randomize locus order
  SNP_colnames <- sprintf("SNP_%s_%d", SNP_locus_types, as.integer(ave(SNP_locus_types, SNP_locus_types, FUN = seq_along))) #create informative SNP names
  SNP_matrix <- matrix(NA, nrow = N_individuals, ncol = N_SNP_loci, dimnames = list(individual_ids, SNP_colnames)) #initialize SNP genotype matrix
  for (locus_index in which(SNP_locus_types == "N")) { #simulate neutral loci shared across all clusters
    shared_allele_frequency <- rbeta(1, 0.5, 0.5) #draw allele frequency from U-shaped beta distribution
    SNP_matrix[, locus_index] <- rbinom(N_individuals, 2, shared_allele_frequency) #simulate diploid genotypes coded as 0/1/2
  }
  for (locus_index in which(SNP_locus_types == "R")) { #simulate random loci with additional genotype perturbation
    random_locus_allele_frequency <- rbeta(1, 0.5, 0.5) #draw allele frequency from U-shaped beta distribution
    genotype_values <- rbinom(N_individuals, 2, random_locus_allele_frequency) #simulate diploid genotypes
    N_randomized_individuals <- floor(N_individuals * runif(1, SNP_random_prop_range[1], SNP_random_prop_range[2])) #number of individuals to randomize at this locus
    if (N_randomized_individuals > 0) {
      randomized_individual_indices <- sample(N_individuals, N_randomized_individuals) #select individuals to randomize
      genotype_values[randomized_individual_indices] <- sample(0:2, N_randomized_individuals, replace = TRUE) #replace some genotypes to mimic noise
    }
    SNP_matrix[, locus_index] <- genotype_values
  }
  if (N_clusters > 1) {
    reached_target_Fst <- FALSE
    max_iterations <- 200 #maximum number of attempts to generate differentiated loci with target Fst structure
    for (iteration_index in seq_len(max_iterations)) {
      for (locus_index in which(SNP_locus_types == "D")) { #simulate differentiated loci with cluster-specific allele frequencies
        locus_target_Fst <- runif(1, SNP_differentiated_Fst_range[1], min(SNP_differentiated_Fst_range[2], 1)) #draw strong per-locus Fst value
        global_allele_frequency <- runif(1, 0.25, 0.75) #draw global allele frequency across clusters
        beta_shape1 <- global_allele_frequency * (1 - locus_target_Fst) / (locus_target_Fst * 0.8) #derive beta parameter a for cluster-specific allele frequencies
        beta_shape2 <- (1 - global_allele_frequency) * (1 - locus_target_Fst) / (locus_target_Fst * 0.8) #derive beta parameter b for cluster-specific allele frequencies
        cluster_specific_allele_frequencies <- rbeta(N_clusters, beta_shape1, beta_shape2) #draw cluster-specific allele frequencies
        for (cluster_index in seq_len(N_clusters)) {
          cluster_member_indices <- which(true_cluster_labels == cluster_index)
          SNP_matrix[cluster_member_indices, locus_index] <- rbinom(length(cluster_member_indices), 2, cluster_specific_allele_frequencies[cluster_index]) #simulate diploid genotypes per cluster
        }
      }
      pairwise_Fst_matrix <- matrix(NA,  nrow = N_clusters, ncol = N_clusters, dimnames = list(paste0("k", seq_len(N_clusters)), paste0("k", seq_len(N_clusters)))) #initialize pairwise Fst matrix
      for (cluster_i in seq_len(N_clusters - 1)) {
        for (cluster_j in (cluster_i + 1):N_clusters) {
          selected_individuals <- true_cluster_labels %in% c(cluster_i, cluster_j) #select individuals from the focal cluster pair
          per_locus_Fst_values <- sapply(seq_len(N_SNP_loci), function(locus_index) {
            cluster_allele_frequencies <- tapply(SNP_matrix[selected_individuals, locus_index] / 2, true_cluster_labels[selected_individuals], mean) #calculate allele frequency per cluster
            mean_allele_frequency <- mean(cluster_allele_frequencies) #calculate mean allele frequency across the two clusters
            if (mean_allele_frequency * (1 - mean_allele_frequency) > 0) {
              mean((cluster_allele_frequencies - mean_allele_frequency)^2) / (mean_allele_frequency * (1 - mean_allele_frequency))
            } else {
              0
            } #calculate simple per-locus Fst
          })
          pairwise_Fst_matrix[cluster_i, cluster_j] <- mean(per_locus_Fst_values, na.rm = TRUE) #store mean pairwise Fst across loci
          pairwise_Fst_matrix[cluster_j, cluster_i] <- pairwise_Fst_matrix[cluster_i, cluster_j] #mirror lower and upper triangle
        }
      }
      pairwise_Fst_values <- pairwise_Fst_matrix[lower.tri(pairwise_Fst_matrix)] #extract pairwise Fst values
      observed_mean_pairwise_Fst <- mean(pairwise_Fst_values, na.rm = TRUE) #mean pairwise Fst across cluster pairs
      Fst_tolerance <- 0.005 #tolerance around allowed Fst range for simulation
      Fst_lower_bound <- SNP_target_Fst_range[1] - Fst_tolerance
      Fst_upper_bound <- SNP_target_Fst_range[2] + Fst_tolerance
      if (!is.na(observed_mean_pairwise_Fst) &&
          observed_mean_pairwise_Fst >= Fst_lower_bound &&
          observed_mean_pairwise_Fst <= Fst_upper_bound) {
        reached_target_Fst <- TRUE
        break
      }
    }
    if (!reached_target_Fst) stop("After ", max_iterations, " iterations, observed mean pairwise Fst never fell in [", round(Fst_lower_bound, 3), ", ", round(Fst_upper_bound, 3), "]")
    observed_Fst <- mean(pairwise_Fst_matrix[lower.tri(pairwise_Fst_matrix)], na.rm = TRUE) #calculate observed average pairwise Fst
  } else {
    for (locus_index in which(SNP_locus_types == "D")) { #for k = 1, differentiated loci are simulated without cluster structure
      unstructured_allele_frequency <- runif(1, 0.15, 0.5)
      SNP_matrix[, locus_index] <- rbinom(N_individuals, 2, unstructured_allele_frequency)
    }
    observed_Fst <- 0
  }
  SNP_df <- as.data.frame(SNP_matrix)
  SNP_df$Cluster <- factor(true_cluster_labels, labels = paste0("k", seq_len(N_clusters))) #append true cluster labels for convenience
  rownames(SNP_df) <- individual_ids
  
  # Morphological variables
  # Simulate mixture of neutral traits and differentiated traits:
  # 1. Neutral traits share same normal distribution across all clusters
  # 2. Differentiated traits have cluster-specific means arranged on simplex
  morph_differentiated_n <- N_clusters - 1
  morph_neutral_n <- N_morph_traits - morph_differentiated_n
  stopifnot(morph_differentiated_n + morph_neutral_n == N_morph_traits)
  morph_trait_types <- c(rep("D", morph_differentiated_n), rep("N", morph_neutral_n)) #assign differentiated and neutral trait classes
  morph_trait_types <- sample(morph_trait_types) #randomize trait order
  morph_differentiated_col_indices <- which(morph_trait_types == "D")
  morph_colnames <- sprintf("Trait_%s_%d", morph_trait_types, as.integer(ave(morph_trait_types, morph_trait_types, FUN = seq_along))) #create trait names
  morph_matrix <- matrix(NA, nrow = N_individuals, ncol = N_morph_traits, dimnames = list(individual_ids, morph_colnames)) #initialize morphology matrix
  morph_trait_sd <- runif(1, morph_trait_sd_range[1], morph_trait_sd_range[2]) #draw shared standard deviation for morphology
  for (trait_index in which(morph_trait_types == "N")) {
    morph_trait_mean_neutral <- sample(1:100, 1) #draw a trait-specific mean
    morph_matrix[, trait_index] <- rnorm(N_individuals, mean = morph_trait_mean_neutral, sd = morph_trait_sd) #simulate neutral traits with trait-specific mean
  }
  if (morph_differentiated_n > 0) {
    morph_cluster_mean_matrix <- matrix(0, nrow = N_clusters, ncol = morph_differentiated_n) #initialize cluster means for differentiated traits
    if (N_clusters == 2) {
      morph_trait_distance <- runif(1, morph_trait_distance_range[1], morph_trait_distance_range[2]) #draw separation distance between the two clusters
      morph_cluster_mean_matrix[, 1] <- c(-morph_trait_distance / 2, morph_trait_distance / 2) #place two clusters symmetrically around zero
    } else if (morph_differentiated_n >= (N_clusters - 1)) {
      morph_trait_distance <- runif(1, morph_trait_distance_range[1], morph_trait_distance_range[2]) #draw simplex spacing distance
      simplex_mean_coordinates <- generate.simplex.coordinates(N_clusters, morph_trait_distance) #generate equally spaced cluster means
      morph_cluster_mean_matrix[, 1:(N_clusters - 1)] <- simplex_mean_coordinates
    }
    if (morph_differentiated_n > (N_clusters - 1)) morph_cluster_mean_matrix[, N_clusters:morph_differentiated_n] <- matrix(rnorm((morph_differentiated_n - (N_clusters - 1)) * N_clusters, 10, 2), nrow = N_clusters) #fill any extra differentiated traits with random cluster means
    for (differentiated_trait_index in seq_along(morph_differentiated_col_indices)) {
      trait_col_index <- morph_differentiated_col_indices[differentiated_trait_index]
      for (cluster_index in seq_len(N_clusters)) {
        cluster_member_indices <- which(true_cluster_labels == cluster_index)
        morph_matrix[cluster_member_indices, trait_col_index] <- rnorm(length(cluster_member_indices), mean = morph_cluster_mean_matrix[cluster_index, differentiated_trait_index], sd = morph_trait_sd) #simulate cluster-specific trait values
      }
    }
  }
  morph_df <- as.data.frame(morph_matrix)
  rownames(morph_df) <- individual_ids
  
  # Climate variables
  # Simulate environmental variables similar to morphology but using skew-normal distributions:
  # 1. Neutral variables share same mean across all clusters
  # 2. Differentiated variables have cluster-specific means arranged on simplex
  # 3. Variable-level skewness introduces realistic non-normality
  climate_differentiated_n <- N_clusters - 1
  climate_neutral_n <- N_climate_variables - climate_differentiated_n
  climate_variable_types <- c(rep("D", climate_differentiated_n), rep("N", climate_neutral_n))
  climate_variable_types <- sample(climate_variable_types) #randomize variable order
  climate_differentiated_col_indices <- which(climate_variable_types == "D")
  climate_neutral_col_indices <- which(climate_variable_types == "N")
  climate_colnames <- sprintf("Climate_%s_%d", climate_variable_types, as.integer(ave(climate_variable_types, climate_variable_types, FUN = seq_along))) #create climate variable names
  climate_matrix <- matrix(NA, nrow = N_individuals, ncol = N_climate_variables, dimnames = list(individual_ids, climate_colnames)) #initialize climate matrix
  climate_variable_sd <- runif(1, climate_variables_sd_range[1], climate_variables_sd_range[2]) #draw shared scale parameter
  climate_cluster_mean_matrix <- matrix(0, nrow = N_clusters, ncol = N_climate_variables) #initialize cluster means for all climate variables
  if (climate_differentiated_n > 0) {
    if (N_clusters == 2) {
      climate_variable_distance <- runif(1, climate_variables_distance_range[1], climate_variables_distance_range[2]) #draw mean separation for the two clusters
      climate_cluster_mean_matrix[, climate_differentiated_col_indices[1]] <- c(-climate_variable_distance / 2, climate_variable_distance / 2)
    } else if (climate_differentiated_n >= (N_clusters - 1)) {
      climate_variable_distance <- runif(1, climate_variables_distance_range[1], climate_variables_distance_range[2]) #draw simplex spacing distance
      simplex_mean_coordinates <- generate.simplex.coordinates(N_clusters, climate_variable_distance) #generate equally spaced cluster means
      climate_cluster_mean_matrix[, climate_differentiated_col_indices[1:(N_clusters - 1)]] <- simplex_mean_coordinates
    }
    if (climate_differentiated_n > (N_clusters - 1)) {
      extra_differentiated_indices <- N_clusters:climate_differentiated_n
      climate_cluster_mean_matrix[, climate_differentiated_col_indices[extra_differentiated_indices]] <- matrix(rnorm(length(extra_differentiated_indices) * N_clusters, 0, 2), nrow = N_clusters) #fill any extra differentiated variables with random means
    }
  }
  differentiated_variable_skews <- matrix(runif(N_clusters * climate_differentiated_n, -6, 6), nrow = N_clusters) #draw cluster-specific skew parameters for differentiated variables
  neutral_variable_skews <- runif(length(climate_neutral_col_indices), -6, 6) #draw variable-specific skew parameters for neutral variables
  for (variable_index in seq_len(N_climate_variables)) {
    if (variable_index %in% climate_differentiated_col_indices) {
      for (cluster_index in seq_len(N_clusters)) {
        cluster_member_indices <- which(true_cluster_labels == cluster_index)
        cluster_mean_value <- climate_cluster_mean_matrix[cluster_index, variable_index]
        cluster_skew_value <- differentiated_variable_skews[cluster_index, which(climate_differentiated_col_indices == variable_index)]
        climate_matrix[cluster_member_indices, variable_index] <- sn::rsn(length(cluster_member_indices), xi = cluster_mean_value, omega = climate_variable_sd, alpha = cluster_skew_value) #simulate skew-normal climate values for each cluster
      }
    } else {
      neutral_variable_position <- which(climate_neutral_col_indices == variable_index)
      climate_mean_neutral <- sample(1:100, 1)
      climate_matrix[, variable_index] <- sn::rsn(N_individuals, xi = climate_mean_neutral, omega = climate_variable_sd, alpha = neutral_variable_skews[neutral_variable_position]) #simulate skew-normal climate values shared across clusters
    }
  }
  climate_df <- as.data.frame(climate_matrix)
  rownames(climate_df) <- individual_ids
  
  # Host data
  # For each cluster:
  # 1. Assign dominant host to most individuals
  # 2. Assign remaining individuals to one or more secondary hosts
  # 3. Use gamma-distributed weights to create uneven secondary-host use
  host_ids <- paste0("Host_", seq_len(N_hosts))
  dominant_host_prop <- runif(1, host_dominant_prop_range[1], host_dominant_prop_range[2]) #draw dominant host proportion for this simulation
  host_matrix <- matrix(0, nrow = N_individuals, ncol = N_hosts, dimnames = list(individual_ids, host_ids)) #initialize binary host assignment matrix
  for (cluster_index in seq_len(N_clusters)) {
    cluster_member_indices <- which(true_cluster_labels == cluster_index)
    N_cluster_members <- length(cluster_member_indices)
    N_dominant_host_members <- round(dominant_host_prop * N_cluster_members)
    dominant_host_member_indices <- sample(cluster_member_indices, N_dominant_host_members) #select individuals assigned to the dominant host
    host_matrix[dominant_host_member_indices, cluster_index] <- 1 #assign dominant host
    remaining_member_indices <- setdiff(cluster_member_indices, dominant_host_member_indices) #identify individuals to assign to secondary hosts
    N_secondary_hosts <- rpois(1, lambda = 1) + 1 #draw number of secondary hosts from a zero-truncated Poisson distribution
    available_secondary_hosts <- setdiff(seq_len(N_hosts), cluster_index) #exclude the dominant host column
    N_available_secondary_hosts <- min(N_secondary_hosts, length(available_secondary_hosts))
    selected_secondary_hosts <- if (N_available_secondary_hosts > 0) sample(available_secondary_hosts, N_available_secondary_hosts, replace = FALSE) else integer(0) #sample actual secondary hosts
    if (length(selected_secondary_hosts) > 0 && length(remaining_member_indices) > 0) {
      if (length(selected_secondary_hosts) == 1) {
        host_matrix[remaining_member_indices, selected_secondary_hosts] <- 1 #assign all remaining individuals to the only available secondary host
      } else {
        raw_secondary_host_weights <- rgamma(length(selected_secondary_hosts), shape = 0.8, rate = 1) #draw skewed weights for secondary hosts
        if (all(!is.na(raw_secondary_host_weights)) && all(raw_secondary_host_weights > 0) && sum(raw_secondary_host_weights) > 0 && length(raw_secondary_host_weights) == length(selected_secondary_hosts)) {
          secondary_host_probabilities <- raw_secondary_host_weights / sum(raw_secondary_host_weights) #normalize weights to probabilities
          if (length(secondary_host_probabilities) == length(selected_secondary_hosts) && all(!is.na(secondary_host_probabilities)) && all(secondary_host_probabilities > 0) && length(remaining_member_indices) > 0) {
            assigned_secondary_hosts <- sample(selected_secondary_hosts, length(remaining_member_indices), prob = secondary_host_probabilities, replace = TRUE) #assign each remaining individual to a secondary host
            for (secondary_host_index in selected_secondary_hosts) host_matrix[remaining_member_indices[assigned_secondary_hosts == secondary_host_index], secondary_host_index] <- 1
          } else {
            fallback_secondary_hosts <- rep(selected_secondary_hosts, length.out = length(remaining_member_indices)) #fallback: recycle secondary hosts evenly
            host_matrix[cbind(remaining_member_indices, fallback_secondary_hosts)] <- 1
          }
        } else {
          fallback_secondary_hosts <- rep(selected_secondary_hosts, length.out = length(remaining_member_indices)) #fallback: recycle secondary hosts evenly
          host_matrix[cbind(remaining_member_indices, fallback_secondary_hosts)] <- 1
        }
      }
    } else if (length(remaining_member_indices) > 0) {
      host_matrix[remaining_member_indices, cluster_index] <- 1 #fallback if no secondary hosts are available
    }
  }
  
  # Calculate number of different secondary host categories used per cluster
  secondary_host_counts <- sapply(seq_len(N_clusters), function(cluster_index) {
    cluster_member_indices <- which(true_cluster_labels == cluster_index)
    dominant_host_column <- cluster_index
    cluster_host_assignments <- host_matrix[cluster_member_indices, , drop = FALSE]
    secondary_host_assignments <- cluster_host_assignments[, -dominant_host_column, drop = FALSE]
    sum(colSums(secondary_host_assignments) > 0) #count how many secondary hosts are used by at least one individual
  })
  secondary_host_counts <- mean(secondary_host_counts) #take mean number of secondary hosts across clusters
  
  # Simulation summary
  sim_stats <- data.frame(
    N_clusters = N_clusters,
    missing_data_prop = missing_data_prop,
    SNP_differentiated_prop = SNP_differentiated_prop,
    SNP_neutral_prop = SNP_neutral_n / N_SNP_loci,
    SNP_random_prop = SNP_random_prop,
    SNP_differentiated_n = SNP_differentiated_n,
    SNP_neutral_n = SNP_neutral_n,
    SNP_random_n = SNP_random_n,
    SNP_target_Fst_min = SNP_target_Fst_range[1],
    SNP_target_Fst_max = SNP_target_Fst_range[2],
    SNP_observed_Fst = observed_Fst,
    morph_differentiated_n = morph_differentiated_n,
    morph_neutral_n = morph_neutral_n,
    climate_differentiated_n = climate_differentiated_n,
    climate_neutral_n = climate_neutral_n,
    dominant_host_prop = dominant_host_prop,
    secondary_host_species_N = secondary_host_counts,
    stringsAsFactors = FALSE
  )
  
  # Create function to randomly introduce missing values
  add.NA <- function(input_matrix, missing_prop) {
    input_matrix <- as.matrix(input_matrix)
    non_missing_indices <- which(!is.na(input_matrix), arr.ind = TRUE)
    N_missing_cells <- floor(missing_prop * prod(dim(input_matrix)))
    if (N_missing_cells > 0) {
      selected_missing_indices <- non_missing_indices[sample(seq_len(nrow(non_missing_indices)), N_missing_cells), , drop = FALSE]
      for (missing_index in seq_len(nrow(selected_missing_indices))) {
        input_matrix[selected_missing_indices[missing_index, 1], selected_missing_indices[missing_index, 2]] <- NA
      }
    }
    input_matrix
  }
  
  # Add missing data
  SNP_df_missing <- as.data.frame(add.NA(SNP_matrix, missing_data_prop))
  morph_df_missing <- as.data.frame(add.NA(morph_matrix, missing_data_prop))
  climate_df_missing <- as.data.frame(add.NA(climate_matrix, missing_data_prop))
  host_df_missing <- as.data.frame(add.NA(host_matrix, missing_data_prop))
  rownames(SNP_df_missing) <- individual_ids
  rownames(morph_df_missing) <- individual_ids
  rownames(climate_df_missing) <- individual_ids
  rownames(host_df_missing) <- individual_ids
  
  # Return results
  list(SNP = SNP_df_missing,
       Morphology = morph_df_missing,
       Climate = climate_df_missing,
       Host = host_df_missing,
       sim_stats = sim_stats,
       Cluster = true_cluster_labels)
}


## Create function to add missing values to simulation using normally distributed selected bad individuals and selected bad variables
add.missing.data.to.simulation <- function(simulation_data,
                                           missing_data_prop,
                                           missing_seed = 1,
                                           bad_individual_prop_mean = missing_bad_individual_prop_mean,
                                           bad_individual_prop_CI95 = missing_bad_individual_prop_CI95,
                                           bad_variable_prop_mean = missing_bad_variable_prop_mean,
                                           bad_variable_prop_CI95 = missing_bad_variable_prop_CI95,
                                           max_sampling_attempts = 1000) {
  sample.normal.proportion.from.CI95 <- function(mean_value, CI95_range, max_attempts) {
    if (!is.numeric(mean_value) || length(mean_value) != 1 || is.na(mean_value)) stop("mean_value must be one non-NA numeric value")
    if (!is.numeric(CI95_range) || length(CI95_range) != 2 || any(is.na(CI95_range))) stop("CI95_range must contain two non-NA numeric values")
    if (CI95_range[1] >= CI95_range[2]) stop("CI95_range lower value must be smaller than upper value")
    if (mean_value <= 0 || mean_value >= 1) stop("mean_value must be between 0 and 1")
    if (CI95_range[1] <= 0 || CI95_range[2] >= 1) stop("CI95_range must be fully between 0 and 1")
    if (abs(mean_value - mean(CI95_range)) > 1e-8) stop("For a normal distribution, mean_value must be the midpoint of CI95_range")
    sd_value <- (CI95_range[2] - CI95_range[1]) / (2 * stats::qnorm(0.975))
    for (attempt_index in seq_len(max_attempts)) {
      sampled_value <- stats::rnorm(1, mean = mean_value, sd = sd_value)
      if (sampled_value > 0 && sampled_value < 1) return(sampled_value)
    }
    stop("Could not sample valid normally distributed missingness proportion between 0 and 1")
  }
  set.seed(missing_seed)
  bad_individual_prop <- sample.normal.proportion.from.CI95(mean_value = bad_individual_prop_mean, CI95_range = bad_individual_prop_CI95,  max_attempts = max_sampling_attempts)
  bad_variable_prop <- sample.normal.proportion.from.CI95(mean_value = bad_variable_prop_mean, CI95_range = bad_variable_prop_CI95,  max_attempts = max_sampling_attempts)
  N_individuals_current <- nrow(simulation_data$SNP)
  N_bad_individuals <- max(1, round(bad_individual_prop * N_individuals_current))
  bad_individual_ids <- rownames(simulation_data$SNP)[sample(seq_len(N_individuals_current), N_bad_individuals, replace = FALSE)]
  add.NA <- function(input_matrix, missing_prop, seed_offset = 0) {
    input_matrix <- as.matrix(input_matrix)
    bad_individual_indices <- which(rownames(input_matrix) %in% bad_individual_ids)
    set.seed(missing_seed + seed_offset)
    N_bad_variables <- max(1, round(bad_variable_prop * ncol(input_matrix)))
    N_bad_variables <- min(N_bad_variables, ncol(input_matrix))
    bad_variable_indices <- sample(seq_len(ncol(input_matrix)), N_bad_variables, replace = FALSE)
    good_variable_indices <- setdiff(seq_len(ncol(input_matrix)), bad_variable_indices)
    target_NA_cells_per_bad_individual <- round(missing_prop * ncol(input_matrix))
    target_NA_cells_per_bad_variable <- round(missing_prop * nrow(input_matrix))
    required_bad_variable_NA_per_bad_individual <- max(0, target_NA_cells_per_bad_individual - length(good_variable_indices))
    if (required_bad_variable_NA_per_bad_individual > 0) {
      for (bad_individual_index in bad_individual_indices) {
        available_variable_indices <- bad_variable_indices[!is.na(input_matrix[bad_individual_index, bad_variable_indices])]
        available_variable_NA_counts <- sapply(available_variable_indices, function(variable_index) sum(is.na(input_matrix[, variable_index])))
        available_variable_indices <- available_variable_indices[available_variable_NA_counts < target_NA_cells_per_bad_variable]
        available_variable_NA_counts <- available_variable_NA_counts[available_variable_NA_counts < target_NA_cells_per_bad_variable]
        if (length(available_variable_indices) < required_bad_variable_NA_per_bad_individual) stop("Not enough bad variables available to reserve overlap for bad individuals")
        set.seed(missing_seed + seed_offset + 500 + bad_individual_index)
        variable_order <- available_variable_indices[order(available_variable_NA_counts, runif(length(available_variable_indices)))]
        selected_variable_indices <- variable_order[seq_len(required_bad_variable_NA_per_bad_individual)]
        input_matrix[bad_individual_index, selected_variable_indices] <- NA
      }
    }
    if (target_NA_cells_per_bad_variable > 0) {
      for (bad_variable_index in bad_variable_indices) {
        current_NA_cells <- sum(is.na(input_matrix[, bad_variable_index]))
        N_missing_cells <- target_NA_cells_per_bad_variable - current_NA_cells
        if (N_missing_cells < 0) stop("A selected bad variable already has more NA values than requested by missing_data_prop")
        if (N_missing_cells > 0) {
          available_individual_indices <- which(!is.na(input_matrix[, bad_variable_index]))
          non_bad_available_individual_indices <- setdiff(available_individual_indices, bad_individual_indices)
          bad_available_individual_indices <- intersect(available_individual_indices, bad_individual_indices)
          bad_available_individual_NA_counts <- rowSums(is.na(input_matrix[bad_available_individual_indices, , drop = FALSE]))
          bad_available_individual_indices <- bad_available_individual_indices[bad_available_individual_NA_counts < target_NA_cells_per_bad_individual]
          bad_available_individual_NA_counts <- bad_available_individual_NA_counts[bad_available_individual_NA_counts < target_NA_cells_per_bad_individual]
          set.seed(missing_seed + seed_offset + 1000 + bad_variable_index)
          non_bad_individual_order <- non_bad_available_individual_indices[order(runif(length(non_bad_available_individual_indices)))]
          bad_individual_order <- bad_available_individual_indices[order(bad_available_individual_NA_counts, runif(length(bad_available_individual_indices)))]
          individual_order <- c(non_bad_individual_order, bad_individual_order)
          if (N_missing_cells > length(individual_order)) stop("Not enough available individuals to reach requested bad-variable NA proportion")
          selected_individual_indices <- individual_order[seq_len(N_missing_cells)]
          input_matrix[selected_individual_indices, bad_variable_index] <- NA
        }
      }
    }
    if (target_NA_cells_per_bad_individual > 0) {
      for (bad_individual_index in bad_individual_indices) {
        current_NA_cells <- sum(is.na(input_matrix[bad_individual_index, ]))
        N_missing_cells <- target_NA_cells_per_bad_individual - current_NA_cells
        if (N_missing_cells < 0) stop("A selected bad individual already has more NA values than requested by missing_data_prop")
        if (N_missing_cells > 0) {
          available_variable_indices <- good_variable_indices[!is.na(input_matrix[bad_individual_index, good_variable_indices])]
          if (N_missing_cells > length(available_variable_indices)) stop("Not enough good variables available to reach requested bad-individual NA proportion")
          set.seed(missing_seed + seed_offset + 2000 + bad_individual_index)
          variable_order <- available_variable_indices[order(runif(length(available_variable_indices)))]
          selected_variable_indices <- variable_order[seq_len(N_missing_cells)]
          input_matrix[bad_individual_index, selected_variable_indices] <- NA
        }
      }
    }
    final_bad_individual_NA <- unname(rowMeans(is.na(input_matrix[bad_individual_indices, , drop = FALSE])))
    final_bad_variable_NA <- unname(colMeans(is.na(input_matrix[, bad_variable_indices, drop = FALSE])))
    target_bad_individual_NA <- target_NA_cells_per_bad_individual / ncol(input_matrix)
    target_bad_variable_NA <- target_NA_cells_per_bad_variable / nrow(input_matrix)
    if (any(abs(final_bad_individual_NA - target_bad_individual_NA) > 1e-12)) stop("Final bad-individual NA proportion does not match requested missing_data_prop")
    if (any(abs(final_bad_variable_NA - target_bad_variable_NA) > 1e-12)) stop("Final bad-variable NA proportion does not match requested missing_data_prop")
    output_matrix <- as.data.frame(input_matrix)
    attr(output_matrix, "N_bad_variables") <- N_bad_variables
    output_matrix
  }
  simulation_data_missing <- simulation_data
  simulation_data_missing$SNP <- add.NA(simulation_data$SNP, missing_data_prop, seed_offset = 1)
  simulation_data_missing$Morphology <- add.NA(simulation_data$Morphology, missing_data_prop, seed_offset = 2)
  simulation_data_missing$Climate <- add.NA(simulation_data$Climate, missing_data_prop, seed_offset = 3)
  simulation_data_missing$Host <- add.NA(simulation_data$Host, missing_data_prop, seed_offset = 4)
  rownames(simulation_data_missing$SNP) <- rownames(simulation_data$SNP)
  rownames(simulation_data_missing$Morphology) <- rownames(simulation_data$Morphology)
  rownames(simulation_data_missing$Climate) <- rownames(simulation_data$Climate)
  rownames(simulation_data_missing$Host) <- rownames(simulation_data$Host)
  simulation_data_missing$sim_stats$missing_data_prop <- missing_data_prop
  simulation_data_missing$sim_stats$bad_individual_prop <- bad_individual_prop
  simulation_data_missing$sim_stats$bad_variable_prop <- bad_variable_prop
  simulation_data_missing$sim_stats$bad_individual_prop_mean <- bad_individual_prop_mean
  simulation_data_missing$sim_stats$bad_individual_prop_CI95_lower <- bad_individual_prop_CI95[1]
  simulation_data_missing$sim_stats$bad_individual_prop_CI95_upper <- bad_individual_prop_CI95[2]
  simulation_data_missing$sim_stats$bad_variable_prop_mean <- bad_variable_prop_mean
  simulation_data_missing$sim_stats$bad_variable_prop_CI95_lower <- bad_variable_prop_CI95[1]
  simulation_data_missing$sim_stats$bad_variable_prop_CI95_upper <- bad_variable_prop_CI95[2]
  simulation_data_missing$sim_stats$N_bad_individuals <- N_bad_individuals
  simulation_data_missing$sim_stats$N_bad_variables_SNP <- attr(simulation_data_missing$SNP, "N_bad_variables")
  simulation_data_missing$sim_stats$N_bad_variables_Morphology <- attr(simulation_data_missing$Morphology, "N_bad_variables")
  simulation_data_missing$sim_stats$N_bad_variables_Climate <- attr(simulation_data_missing$Climate, "N_bad_variables")
  simulation_data_missing$sim_stats$N_bad_variables_Host <- attr(simulation_data_missing$Host, "N_bad_variables")
  simulation_data_missing
}


## Create function to extract mean QE and TE across all layers for SOM model replicate
extract.QE.TE <- function(som.model) {
  codebook_layers <- kohonen::getCodes(som.model)
  if (!is.list(codebook_layers)) codebook_layers <- list(codebook_layers)
  codebook_layers <- lapply(codebook_layers, function(codebook_matrix) {
    codebook_matrix <- as.matrix(codebook_matrix)
    if (ncol(codebook_matrix) == 1) dim(codebook_matrix) <- c(nrow(codebook_matrix), 1) #ensure one-column codebooks remain two-dimensional
    codebook_matrix
  })
  data_layers <- som.model$data
  if (!is.list(data_layers)) data_layers <- list(as.matrix(data_layers))
  QE_per_layer <- numeric(length(data_layers)) #mean distance from each sample to its BMU
  TE_per_layer <- numeric(length(data_layers)) #topographic error rate per layer
  for (layer_index in seq_along(data_layers)) {
    layer_codebook <- codebook_layers[[layer_index]]
    layer_data <- data_layers[[layer_index]]
    sample_QE <- numeric(nrow(layer_data))
    for (sample_index in seq_len(nrow(layer_data))) {
      non_missing_dimensions <- which(!is.na(layer_data[sample_index, ]))
      if (length(non_missing_dimensions) > 0) {
        best_matching_unit_index <- som.model$unit.classif[sample_index]
        squared_differences <- (layer_data[sample_index, non_missing_dimensions] - layer_codebook[best_matching_unit_index, non_missing_dimensions])^2
        sample_QE[sample_index] <- sqrt(sum(squared_differences))
      } else {
        sample_QE[sample_index] <- NA
      }
    }
    QE_per_layer[layer_index] <- mean(sample_QE, na.rm = TRUE)
    unit_distance_matrix <- kohonen::unit.distances(som.model$grid)
    sample_TE <- logical(nrow(layer_data))
    for (sample_index in seq_len(nrow(layer_data))) {
      non_missing_dimensions <- which(!is.na(layer_data[sample_index, ]))
      if (length(non_missing_dimensions) < 2) {
        sample_TE[sample_index] <- NA
        next
      }
      difference_matrix <- t(layer_codebook[, non_missing_dimensions]) - layer_data[sample_index, non_missing_dimensions]
      squared_distances_to_units <- colSums(difference_matrix^2)
      ordered_unit_indices <- order(squared_distances_to_units)
      sample_TE[sample_index] <- unit_distance_matrix[ordered_unit_indices[1], ordered_unit_indices[2]] > 1
    }
    TE_per_layer[layer_index] <- mean(sample_TE, na.rm = TRUE)
  }
  list(QE = mean(QE_per_layer, na.rm = TRUE), TE = mean(TE_per_layer, na.rm = TRUE)
  )
}


## Create function to map predicted clusters to true clusters using Hungarian algorithm for maximum agreement
get.best.assignment <- function(true_labels, predicted_labels) {
  contingency_table <- table(true_labels, predicted_labels)
  optimal_assignment <- clue::solve_LSAP(contingency_table, maximum = TRUE)
  true_levels <- rownames(contingency_table)
  predicted_levels <- colnames(contingency_table)
  predicted_to_true_map <- setNames(true_levels, predicted_levels[optimal_assignment]) #map predicted cluster -> matched true cluster
  mapped_predicted_labels <- unname(predicted_to_true_map[as.character(predicted_labels)])
  mapped_predicted_labels
}


## Create function to calculate proportion correctly assigned after best matching
get.accuracy <- function(true_labels, predicted_labels) {
  if (length(unique(predicted_labels)) != length(unique(true_labels))) return(NA)
  contingency_table <- table(true_labels, predicted_labels)
  optimal_assignment <- clue::solve_LSAP(contingency_table, maximum = TRUE)
  true_levels <- rownames(contingency_table)
  predicted_levels <- colnames(contingency_table)
  predicted_to_true_map <- setNames(true_levels, predicted_levels[optimal_assignment]) #map predicted cluster -> matched true cluster
  mapped_predicted_labels <- unname(predicted_to_true_map[as.character(predicted_labels)])
  mean(as.character(true_labels) == as.character(mapped_predicted_labels))
}


## Create function to extract replicate-specific cluster labels from clustering output
get.cluster.labels.for.replicate <- function(cluster_assignment, replicate_index) {
  cluster_assignment <- as.matrix(cluster_assignment)
  if (ncol(cluster_assignment) >= replicate_index) {
    as.integer(cluster_assignment[, replicate_index]) #use replicate-specific cluster labels when available
  } else {
    as.integer(cluster_assignment[, 1]) #fallback to first clustering column if only one is available
  }
}


## Create function to add simulation summary statistics to benchmark results
add.sim.stats <- function(stats_df, sim_stats_row = NULL) {
  if (is.null(sim_stats_row)) return(stats_df)
  sim_stats_row <- as.data.frame(sim_stats_row, stringsAsFactors = FALSE)
  for (column_name in colnames(sim_stats_row)) {
    stats_df[[column_name]] <- sim_stats_row[[column_name]][1]
  }
  stats_df
}


## Create function to return placeholder rows when SOM/clustering fails
make.failed.stats.df <- function(group_column_name,
                                 group_value,
                                 sim_stats_row = NULL,
                                 N_rows = N_replicates_SOM) {
  stats_df <- data.frame(Time = rep(NA_real_, N_rows),
                         K_inferred = rep(NA_real_, N_rows),
                         ARI = rep(NA_real_, N_rows),
                         Acc = rep(NA_real_, N_rows),
                         K_correct = rep(NA, N_rows),
                         K_far_off = rep(NA, N_rows),
                         QE = rep(NA_real_, N_rows),
                         TE = rep(NA_real_, N_rows),
                         stringsAsFactors = FALSE)
  stats_df[[group_column_name]] <- rep(group_value, N_rows)
  stats_df <- stats_df[, c(group_column_name, setdiff(colnames(stats_df), group_column_name)), drop = FALSE]
  stats_df <- add.sim.stats(stats_df, sim_stats_row)
  stats_df
}


## Create function to bind simulation result data frames and retain original simulation indices
bind.results.with.sim <- function(result_list) {
  keep_indices <- which(vapply(result_list, function(result_object) {
    !is.null(result_object) && is.data.frame(result_object) && nrow(result_object) > 0
  }, logical(1)))
  if (length(keep_indices) == 0) return(NULL)
  for (simulation_index in keep_indices) {
    result_list[[simulation_index]]$sim <- simulation_index
  }
  do.call(rbind, result_list[keep_indices])
}


## Create function to simulate standard datasets with fixed seeds
simulate.standard.datasets <- function(N.simulations,
                                       N.clusters,
                                       missing_data_prop = 0.05,
                                       simulation_seeds = NULL,
                                       max_seed_tries_multiplier = 20) {
  if (is.null(simulation_seeds)) {
    set.seed(2)
    simulation_seeds <- sample(1:1e8, N.simulations * max_seed_tries_multiplier)
  }
  if (length(simulation_seeds) < N.simulations) stop("Length of simulation_seeds must be at least N.simulations")
  simulation_results <- vector("list", N.simulations)
  failed_simulations <- list()
  successful_simulation_index <- 1
  seed_index <- 1
  while (successful_simulation_index <= N.simulations) {
    if (seed_index > length(simulation_seeds)) stop("Ran out of simulation seeds before obtaining ", N.simulations, " successful simulations - increase simulation_seeds or relax the Fst constraint")
    set.seed(simulation_seeds[seed_index])
    simulation_output <- tryCatch({
      simulate.data(N_clusters = N.clusters, missing_data_prop = missing_data_prop, sim_id = successful_simulation_index)
    }, error = function(error_message) {error_message})
    if (inherits(simulation_output, "error")) {
      failed_simulations[[length(failed_simulations) + 1]] <- data.frame(sim = successful_simulation_index, seed = simulation_seeds[seed_index], error = conditionMessage(simulation_output), stringsAsFactors = FALSE)
      seed_index <- seed_index + 1
      next
    }
    simulation_results[[successful_simulation_index]] <- simulation_output
    message("Successfully simulated dataset ", successful_simulation_index, " of ", N.simulations)
    successful_simulation_index <- successful_simulation_index + 1
    seed_index <- seed_index + 1
  }
  
  attr(simulation_results, "failed_simulations") <- if (length(failed_simulations) > 0) {
    do.call(rbind, failed_simulations)
  } else {
    NULL
  }
  simulation_results
}


## Create function to run one SOM benchmark and return replicate-level statistics
run.SOM.benchmark <- function(input_data,
                              true_labels,
                              group_column_name,
                              group_value,
                              clustering_method,
                              training_neighborhoods,
                              max_k,
                              learning_rate_tuning_option = FALSE,
                              sim_stats_row = NULL,
                              verbose = FALSE) {
  if (is.null(rownames(input_data[[1]]))) stop("Input data must have row names so retained SOM samples can be matched to true labels")
  true_labels_original <- true_labels
  names(true_labels_original) <- rownames(input_data[[1]])
  true_K <- length(unique(true_labels_original))
  som_output <- NULL
  clustering_output <- NULL
  elapsed_time <- system.time({
    som_output <- train.SOM(input_data = input_data,
                            learning.rate.tuning = learning_rate_tuning_option,
                            training.neighborhoods = training_neighborhoods,
                            N.steps = N_steps_SOM,
                            grid.size = SOM_grid_size,
                            N.replicates = N_replicates_SOM,
                            max.NA.row = max_NA_row_SOM,
                            max.NA.col = max_NA_col_SOM,
                            save.SOM.results = FALSE,
                            verbose = verbose)
    clustering_output <- clustering.SOM(SOM.output = som_output,
                                        max.k = max_k,
                                        BIC.thresh = BIC_threshold_SOM,
                                        verbose = verbose,
                                        save.SOM.results = FALSE,
                                        clustering.method = clustering_method)
  })[["elapsed"]]
  if (is.null(som_output$som_models) || length(som_output$som_models) == 0) stop("SOM training returned no SOM models for ", group_column_name, " = ", group_value)
  if (is.null(clustering_output$cluster_assignment)) stop("clustering.SOM returned no cluster_assignment for ", group_column_name, " = ", group_value)
  som_models <- som_output$som_models
  total_time <- elapsed_time
  Time = time_per_simulation,
  stats_per_replicate <- lapply(seq_along(som_models), function(replicate_index) {
    som_model <- som_models[[replicate_index]]
    QE_TE <- extract.QE.TE(som_model)
    predicted_labels <- get.cluster.labels.for.replicate(clustering_output$cluster_assignment, replicate_index)
    if (any(is.na(predicted_labels))) stop("Predicted cluster labels contain NA for ", group_column_name, " = ", group_value, ", replicate ", replicate_index)
    retained_ids <- rownames(som_model$data[[1]])
    if (is.null(retained_ids)) stop("SOM model retained sample IDs are missing for ", group_column_name, " = ", group_value, ", replicate ", replicate_index)
    true_labels_retained <- true_labels_original[retained_ids]
    if (any(is.na(true_labels_retained))) stop("Some retained SOM sample IDs could not be matched to true labels for ", group_column_name, " = ", group_value, ", replicate ", replicate_index)
    true_labels_retained <- as.integer(as.factor(true_labels_retained))
    if (length(predicted_labels) != length(true_labels_retained)) stop("Predicted labels and retained true labels have different lengths for ", group_column_name, " = ", group_value, ", replicate ", replicate_index, ": predicted = ", length(predicted_labels), ", true = ", length(true_labels_retained))
    K_inferred <- length(unique(predicted_labels))
    ARI <- mclust::adjustedRandIndex(true_labels_retained, predicted_labels)
    Acc <- if (K_inferred == true_K) get.accuracy(true_labels_retained, predicted_labels) else NA_real_
    K_correct <- (K_inferred == true_K)
    K_far_off <- (abs(K_inferred - true_K) >= 2)
    stats_df <- data.frame(Time = time_per_simulation,
                           K_inferred = K_inferred,
                           ARI = ARI,
                           Acc = Acc,
                           K_correct = K_correct,
                           K_far_off = K_far_off,
                           QE = QE_TE$QE,
                           TE = QE_TE$TE,
                           stringsAsFactors = FALSE)
    stats_df[[group_column_name]] <- group_value
    stats_df <- stats_df[, c(group_column_name, setdiff(colnames(stats_df), group_column_name)), drop = FALSE]
    stats_df
  })
  if (length(stats_per_replicate) == 0) stop("No replicate-level statistics were produced for ", group_column_name, " = ", group_value)
  stats_df <- do.call(rbind, stats_per_replicate)
  stats_df <- add.sim.stats(stats_df, sim_stats_row)
  stats_df
}


## Create function to calculate means without returning NaN for all-NA vectors
safe.mean <- function(input_vector) {
  if (all(is.na(input_vector))) return(NA_real_)
  mean(input_vector, na.rm = TRUE)
}


## Create function to safely convert logical/numeric/character TRUE-FALSE columns
convert.to.logical <- function(input_vector) {
  if (is.logical(input_vector)) return(input_vector)
  if (is.numeric(input_vector)) return(ifelse(is.na(input_vector), NA, input_vector != 0))
  input_vector <- tolower(as.character(input_vector))
  output_vector <- rep(NA, length(input_vector))
  output_vector[input_vector %in% c("true", "t", "1", "yes", "y")] <- TRUE
  output_vector[input_vector %in% c("false", "f", "0", "no", "n")] <- FALSE
  output_vector
}


## Create common plot theme
create.benchmark.plot.theme <- function(rotate_x_axis_text = FALSE) {
  plot_theme <- theme_classic(base_size = benchmark_plot_base_size) +
    theme(axis.title = element_text(size = benchmark_plot_axis_title_size),
          axis.text = element_text(size = benchmark_plot_axis_text_size))
  if (rotate_x_axis_text) {
    plot_theme <- plot_theme +
      theme(axis.text.x = element_text(angle = benchmark_plot_x_axis_text_angle,
                                       hjust = benchmark_plot_x_axis_text_hjust))
  }
  plot_theme
}


## Create function to prepare simulation-level continuous data by averaging SOM technical replicates
prepare.continuous.benchmark.data <- function(input_data,
                                              condition_variable,
                                              response_variable) {
  model_data <- input_data[, c("sim", condition_variable, response_variable)]
  colnames(model_data)[colnames(model_data) == condition_variable] <- "condition"
  colnames(model_data)[colnames(model_data) == response_variable] <- "response"
  model_data <- model_data[!is.na(model_data$condition) & !is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data <- model_data %>%
    dplyr::group_by(sim, condition) %>%
    dplyr::summarise(response = safe.mean(response),
                     .groups = "drop")
  model_data <- model_data[!is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data
}


## Create function to prepare row-level binary data for K outcomes
prepare.binary.benchmark.data <- function(input_data,
                                          condition_variable,
                                          response_variable) {
  model_data <- input_data[, c("sim", condition_variable, response_variable)]
  colnames(model_data)[colnames(model_data) == condition_variable] <- "condition"
  colnames(model_data)[colnames(model_data) == response_variable] <- "response"
  model_data$response <- as.numeric(convert.to.logical(model_data$response))
  model_data <- model_data[!is.na(model_data$condition) & !is.na(model_data$response), , drop = FALSE]
  model_data$sim <- factor(model_data$sim)
  model_data
}



## Create function to extract response-scale emmeans columns
extract.response.emmeans.columns <- function(emmeans_output_df) {
  response_column <- if ("prob" %in% colnames(emmeans_output_df)) {
    "prob"
  } else if ("response" %in% colnames(emmeans_output_df)) {
    "response"
  } else {
    stop("Could not find response column in emmeans output")
  }
  emmeans_output_df$estimated_response <- emmeans_output_df[[response_column]]
  emmeans_output_df$lower_SE <- pmax(0, emmeans_output_df$estimated_response - benchmark_plot_SE_multiplier * emmeans_output_df$SE)
  emmeans_output_df$upper_SE <- pmin(1, emmeans_output_df$estimated_response + benchmark_plot_SE_multiplier * emmeans_output_df$SE)
  emmeans_output_df
}


## Create categorical Gaussian mixed-model plot
plot.categorical.gaussian.benchmark.model <- function(input_data,
                                                      condition_variable,
                                                      response_variable,
                                                      condition_levels,
                                                      x_axis_label,
                                                      y_axis_label,
                                                      estimate_color,
                                                      connect_estimates = FALSE) {
  model_data <- prepare.continuous.benchmark.data(input_data = input_data,
                                                  condition_variable = condition_variable,
                                                  response_variable = response_variable)
  model_data$condition_factor <- factor(as.character(model_data$condition),
                                        levels = as.character(condition_levels))
  model_data <- model_data[!is.na(model_data$condition_factor), , drop = FALSE]
  if (length(unique(model_data$sim)) < 2) stop("Gaussian mixed model requires at least two simulation replicates with non-NA values")
  model_output <- nlme::lme(response ~ condition_factor,
                            random = ~ 1 | sim,
                            data = model_data,
                            na.action = na.omit,
                            method = "REML")
  emmeans_output <- emmeans::emmeans(model_output, ~ condition_factor)
  emmeans_output_df <- as.data.frame(emmeans_output)
  emmeans_output_df$lower_SE <- emmeans_output_df$emmean - benchmark_plot_SE_multiplier * emmeans_output_df$SE
  emmeans_output_df$upper_SE <- emmeans_output_df$emmean + benchmark_plot_SE_multiplier * emmeans_output_df$SE
  model_plot <- ggplot() +
    geom_point(data = model_data,
               aes(x = condition_factor, y = response),
               alpha = benchmark_plot_raw_point_alpha,
               size = benchmark_plot_raw_point_size,
               position = position_jitter(width = benchmark_plot_categorical_jitter_width, height = 0)) +
    geom_errorbar(data = emmeans_output_df,
                  aes(x = condition_factor, ymin = lower_SE, ymax = upper_SE),
                  colour = estimate_color,
                  linewidth = benchmark_plot_errorbar_linewidth,
                  width = benchmark_plot_errorbar_width)
  if (connect_estimates) {
    model_plot <- model_plot +
      geom_line(data = emmeans_output_df,
                aes(x = condition_factor, y = emmean, group = 1),
                colour = estimate_color,
                linewidth = benchmark_plot_estimate_linewidth)
  }
  model_plot <- model_plot +
    geom_point(data = emmeans_output_df,
               aes(x = condition_factor, y = emmean),
               colour = estimate_color,
               size = benchmark_plot_estimate_point_size) +
    scale_x_discrete(drop = FALSE) +
    create.benchmark.plot.theme(rotate_x_axis_text = TRUE) +
    labs(x = x_axis_label,
         y = y_axis_label)
  list(model = model_output,
       emmeans = emmeans_output,
       emmeans_df = emmeans_output_df,
       plot = model_plot)
}


## Create categorical binary mixed-model plot
plot.categorical.binary.benchmark.model <- function(input_data,
                                                    condition_variable,
                                                    response_variable,
                                                    condition_levels,
                                                    x_axis_label,
                                                    y_axis_label,
                                                    estimate_color,
                                                    connect_estimates = FALSE) {
  model_data <- prepare.binary.benchmark.data(input_data = input_data,
                                              condition_variable = condition_variable,
                                              response_variable = response_variable)
  model_data$condition_factor <- factor(as.character(model_data$condition),
                                        levels = as.character(condition_levels))
  model_data <- model_data[!is.na(model_data$condition_factor), , drop = FALSE]
  model_data$sim_condition <- interaction(model_data$sim, model_data$condition_factor, drop = TRUE)
  model_data$sim_condition <- factor(model_data$sim_condition)
  if (length(unique(model_data$sim)) < 2) stop("Binary mixed model requires at least two simulation replicates with non-NA values")
  model_output <- glmmTMB::glmmTMB(response ~ condition_factor + (1 | sim) + (1 | sim_condition),
                                   family = binomial,
                                   data = model_data)
  emmeans_output <- emmeans::emmeans(model_output, ~ condition_factor, type = "response")
  emmeans_output_df <- as.data.frame(emmeans_output)
  emmeans_output_df <- extract.response.emmeans.columns(emmeans_output_df)
  raw_plot_data <- model_data %>%
    dplyr::group_by(sim, condition_factor) %>%
    dplyr::summarise(response = safe.mean(response),
                     .groups = "drop")
  model_plot <- ggplot() +
    geom_point(data = raw_plot_data,
               aes(x = condition_factor, y = response),
               alpha = benchmark_plot_raw_point_alpha,
               size = benchmark_plot_raw_point_size,
               position = position_jitter(width = benchmark_plot_categorical_jitter_width,
                                          height = 0)) +
    geom_errorbar(data = emmeans_output_df,
                  aes(x = condition_factor, ymin = lower_SE, ymax = upper_SE),
                  colour = estimate_color,
                  linewidth = benchmark_plot_errorbar_linewidth,
                  width = benchmark_plot_errorbar_width)
  if (connect_estimates) {
    model_plot <- model_plot +
      geom_line(data = emmeans_output_df,
                aes(x = condition_factor, y = estimated_response, group = 1),
                colour = estimate_color,
                linewidth = benchmark_plot_estimate_linewidth)
  }
  model_plot <- model_plot +
    geom_point(data = emmeans_output_df,
               aes(x = condition_factor, y = estimated_response),
               colour = estimate_color,
               size = benchmark_plot_estimate_point_size) +
    scale_x_discrete(drop = FALSE) +
    create.benchmark.plot.theme(rotate_x_axis_text = TRUE) +
    labs(x = x_axis_label,
         y = y_axis_label)
  list(model = model_output,
       emmeans = emmeans_output,
       emmeans_df = emmeans_output_df,
       plot = model_plot)
}


## Create function to create safe file names from parameter values
create.safe.file.name <- function(input_name) {
  safe_name <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(input_name))
  safe_name <- gsub("_+", "_", safe_name)
  safe_name <- gsub("^_|_$", "", safe_name)
  safe_name
}


## Create function to save currently completed simulation results within one test condition
save.current.condition.results <- function(stats_list, condition_partial_csv, condition_partial_rds) {
  saveRDS(stats_list, condition_partial_rds)
  condition_results <- bind.results.with.sim(stats_list)
  if (!is.null(condition_results)) {
    condition_results <- as.data.frame(condition_results, stringsAsFactors = FALSE)
    write.csv(condition_results, condition_partial_csv, row.names = FALSE)
  }
}


## Create function to save combined results after each completed test condition
save.partial.results <- function(result_list, combined_csv, combined_rds) {
  saveRDS(result_list, combined_rds)
  keep_results <- vapply(result_list, function(result_object) {
    !is.null(result_object) && is.data.frame(result_object) && nrow(result_object) > 0
  }, logical(1))
  if (any(keep_results)) {
    combined_results <- do.call(rbind, result_list[keep_results])
    combined_results <- as.data.frame(combined_results, stringsAsFactors = FALSE)
    write.csv(combined_results, combined_csv, row.names = FALSE)
  }
}




#### Create/load simulations ###################################################

## Set shared simulation parameters
N_clusters <- 3
training_neighborhoods_SOM <- "gaussian"
max_k_SOM <- 8


## Create intermediate output directory
intermediate_results_dir <- file.path(simulation_k3_dir, "Intermediate_results")
if (!dir.exists(intermediate_results_dir)) dir.create(intermediate_results_dir, recursive = TRUE)


## File paths for saving/loading shared complete simulations
sim_data_base_file <- file.path(simulation_k3_dir, "Sim_data_base_complete.rds")
sim_failed_base_csv <- file.path(simulation_k3_dir, "Sim_data_base_complete_failed_simulations.csv")


## Load saved complete base simulations or simulate once
if (file.exists(sim_data_base_file) && !overwrite) {
  simulation_results_base <- readRDS(sim_data_base_file)
  message("Loaded shared complete base simulation data from RDS")
} else {
  set.seed(1)
  simulation_seeds <- sample(1:1e8, N_simulations * 100)
  simulation_results_base <- simulate.standard.datasets(N.simulations = N_simulations,
                                                        N.clusters = N_clusters,
                                                        missing_data_prop = 0,
                                                        simulation_seeds = simulation_seeds)
  failed_simulations <- attr(simulation_results_base, "failed_simulations")
  if (!is.null(failed_simulations)) write.csv(failed_simulations, sim_failed_base_csv, row.names = FALSE)
  saveRDS(simulation_results_base, sim_data_base_file)
  message("Simulations completed and saved to file")
}




#### Test effect of clustering method ##########################################

## Set clustering methods to test
clustering_methods <- c(
  "kmeans+BICelbow",
  "kmeans+BICthreshold",
  "GMM+BICthreshold",
  "hierarchical+DB",
  "HDBSCAN",
  "OPTICS+Silhouette"
)


## Set parameters
N_clusters <- 3
training_neighborhoods_SOM <- "gaussian"
max_k_SOM <- 8


## File paths for saving/loading
sim_results_clustering_methods_csv <- file.path(simulation_k3_dir, "Sim_results_clustering_methods.csv")
sim_results_clustering_methods_rds <- file.path(intermediate_results_dir, "Sim_results_clustering_methods_partial.rds")
if (file.exists(sim_results_clustering_methods_csv)) file.remove(sim_results_clustering_methods_csv)
if (file.exists(sim_results_clustering_methods_rds)) file.remove(sim_results_clustering_methods_rds)


## Load saved results or run as needed
if (file.exists(sim_results_clustering_methods_csv) && !overwrite) {
  full_sim_stats_clustering_methods <- read.csv(sim_results_clustering_methods_csv)
  message("Loaded clustering methods simulation results from CSV, skipping run")
} else {
  simulation_results <- simulation_results_base
  all_results_clustering_methods <- list()
  for (clustering_method_name in clustering_methods) {
    cat("Running for clustering method:", clustering_method_name, "\n")
    condition_safe_name <- create.safe.file.name(clustering_method_name)
    condition_complete_csv <- file.path(intermediate_results_dir, paste0("Sim_results_clustering_method_", condition_safe_name, "_complete.csv"))
    condition_partial_csv <- file.path(intermediate_results_dir, paste0("Sim_results_clustering_method_", condition_safe_name, "_partial.csv"))
    condition_partial_rds <- file.path(intermediate_results_dir, paste0("Sim_results_clustering_method_", condition_safe_name, "_partial.rds"))
    if (file.exists(condition_complete_csv) && !overwrite) {
      stats_df_all <- read.csv(condition_complete_csv)
      message("Loaded completed clustering-method intermediate results for ", clustering_method_name)
    } else {
      stats_list <- vector("list", length(simulation_results))
      skipped_simulations <- list()
      condition_skipped_csv <- file.path(intermediate_results_dir, paste0("Sim_results_clustering_method_", condition_safe_name, "_skipped.csv"))
      for (simulation_index in seq_along(simulation_results)) {
        cat("Testing simulation", simulation_index, "of", length(simulation_results), "for clustering_method:", clustering_method_name, "\n")
        flush.console()
        simulation_data <- simulation_results[[simulation_index]]
        som_input <- list(SNP = simulation_data$SNP,
                          Morphology = simulation_data$Morphology,
                          Climate = simulation_data$Climate,
                          Host = simulation_data$Host)
        stats_list[[simulation_index]] <- tryCatch({
          run.SOM.benchmark(input_data = som_input,
                            true_labels = simulation_data$Cluster,
                            group_column_name = "clustering_method",
                            group_value = clustering_method_name,
                            clustering_method = clustering_method_name,
                            training_neighborhoods = training_neighborhoods_SOM,
                            max_k = max_k_SOM,
                            learning_rate_tuning_option = learning_rate_tuning,
                            sim_stats_row = simulation_data$sim_stats,
                            verbose = verbose_SOM)
        }, error = function(error_message) {
          skip_message <- paste0("Skipped simulation ", simulation_index, " for clustering_method = ", clustering_method_name, " because: ", conditionMessage(error_message))
          message(skip_message)
          skipped_simulations[[length(skipped_simulations) + 1]] <<- data.frame(clustering_method = clustering_method_name, sim = simulation_index, error_message = conditionMessage(error_message), stringsAsFactors = FALSE)
          NULL
        })
        save.current.condition.results(stats_list = stats_list, condition_partial_csv = condition_partial_csv, condition_partial_rds = condition_partial_rds)
        if (length(skipped_simulations) > 0) write.csv(do.call(rbind, skipped_simulations), condition_skipped_csv, row.names = FALSE)
      }
      stats_df_all <- bind.results.with.sim(stats_list)
      if (is.null(stats_df_all)) stop("All simulations failed for clustering_method = ", clustering_method_name)
      write.csv(stats_df_all, condition_complete_csv, row.names = FALSE)
    }
    all_results_clustering_methods[[clustering_method_name]] <- stats_df_all
    save.partial.results(result_list = all_results_clustering_methods,
                         combined_csv = sim_results_clustering_methods_csv,
                         combined_rds = sim_results_clustering_methods_rds)
  }
  all_results_clustering_methods <- all_results_clustering_methods[vapply(all_results_clustering_methods, function(result_object) !is.null(result_object), logical(1))]
  if (length(all_results_clustering_methods) == 0) stop("No clustering-method results to summarize")
  full_sim_stats_clustering_methods <- do.call(rbind, all_results_clustering_methods)
  write.csv(full_sim_stats_clustering_methods, sim_results_clustering_methods_csv, row.names = FALSE)
  message("Clustering methods simulation results saved to CSV")
}


## Summarize results
summary_sim_stats_clustering_methods <- full_sim_stats_clustering_methods %>%
  dplyr::group_by(clustering_method) %>%
  dplyr::summarise(mean_ARI = safe.mean(ARI),
                   mean_K = safe.mean(K_inferred),
                   mean_Acc = safe.mean(Acc),
                   mean_Time = safe.mean(Time[!duplicated(sim)]),
                   failed_rows = sum(is.na(Time) & is.na(K_inferred) & is.na(ARI) & is.na(QE) & is.na(TE)),
                   n = dplyr::n()) %>%
  as.data.frame()
print(summary_sim_stats_clustering_methods)


## Plot results
clustering_methods_ARI_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_clustering_methods,
                                                                               condition_variable = "clustering_method",
                                                                               response_variable = "ARI",
                                                                               condition_levels = clustering_methods,
                                                                               x_axis_label = "Clustering method",
                                                                               y_axis_label = "Adjusted Rand Index",
                                                                               estimate_color = benchmark_plot_color_ARI)
clustering_methods_Acc_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_clustering_methods,
                                                                               condition_variable = "clustering_method",
                                                                               response_variable = "Acc",
                                                                               condition_levels = clustering_methods,
                                                                               x_axis_label = "Clustering method",
                                                                               y_axis_label = "Assignment accuracy",
                                                                               estimate_color = benchmark_plot_color_Acc)
clustering_methods_K_correct_model_plot <- plot.categorical.binary.benchmark.model(input_data = full_sim_stats_clustering_methods,
                                                                                   condition_variable = "clustering_method",
                                                                                   response_variable = "K_correct",
                                                                                   condition_levels = clustering_methods,
                                                                                   x_axis_label = "Clustering method",
                                                                                   y_axis_label = "Probability K correct",
                                                                                   estimate_color = benchmark_plot_color_K_correct)
clustering_methods_K_far_off_model_plot <- plot.categorical.binary.benchmark.model(input_data = full_sim_stats_clustering_methods,
                                                                                   condition_variable = "clustering_method",
                                                                                   response_variable = "K_far_off",
                                                                                   condition_levels = clustering_methods,
                                                                                   x_axis_label = "Clustering method",
                                                                                   y_axis_label = "Probability K far off",
                                                                                   estimate_color = benchmark_plot_color_K_far_off)
clustering_methods_Time_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_clustering_methods,
                                                                                condition_variable = "clustering_method",
                                                                                response_variable = "Time",
                                                                                condition_levels = clustering_methods,
                                                                                x_axis_label = "Clustering method",
                                                                                y_axis_label = "Time (seconds)",
                                                                                estimate_color = benchmark_plot_color_Time)
summary(clustering_methods_ARI_model_plot$model)
summary(clustering_methods_Acc_model_plot$model)
summary(clustering_methods_K_correct_model_plot$model)
summary(clustering_methods_K_far_off_model_plot$model)
summary(clustering_methods_Time_model_plot$model)
clustering_methods_ARI_model_plot$plot
clustering_methods_Acc_model_plot$plot
clustering_methods_K_correct_model_plot$plot
clustering_methods_K_far_off_model_plot$plot
clustering_methods_Time_model_plot$plot




#### Test effect of N.steps ####################################################

## Set parameters                                                                 
N_steps_values <- c(20, 50, 100, 200, 500, 1000, 2000)
clustering_method <- "kmeans+BICelbow"
N_clusters <- 3
training_neighborhoods_SOM <- "gaussian"
max_k_SOM <- 8


## File paths for saving/loading
sim_results_N_steps_csv <- file.path(simulation_k3_dir, "Sim_results_N_steps.csv")
sim_results_N_steps_rds <- file.path(intermediate_results_dir, "Sim_results_N_steps_partial.rds")
if (file.exists(sim_results_N_steps_csv)) file.remove(sim_results_N_steps_csv)
if (file.exists(sim_results_N_steps_rds)) file.remove(sim_results_N_steps_rds)


## Load saved results or run as needed
if (file.exists(sim_results_N_steps_csv) && !overwrite) {
  full_sim_stats_N_steps <- read.csv(sim_results_N_steps_csv)
  message("Loaded N.steps simulation results from CSV, skipping run")
} else {
  simulation_results <- simulation_results_base
  all_results_N_steps <- list()
  original_N_steps_SOM <- N_steps_SOM
  for (N_steps_value in N_steps_values) {
    cat("Running for N.steps:", N_steps_value, "\n")
    N_steps_SOM <- N_steps_value
    condition_safe_name <- create.safe.file.name(N_steps_value)
    condition_complete_csv <- file.path(intermediate_results_dir, paste0("Sim_results_N_steps_", condition_safe_name, "_complete.csv"))
    condition_partial_csv <- file.path(intermediate_results_dir, paste0("Sim_results_N_steps_", condition_safe_name, "_partial.csv"))
    condition_partial_rds <- file.path(intermediate_results_dir, paste0("Sim_results_N_steps_", condition_safe_name, "_partial.rds"))
    if (file.exists(condition_complete_csv) && !overwrite) {
      stats_df_all <- read.csv(condition_complete_csv)
      message("Loaded completed N.steps intermediate results for ", N_steps_value)
    } else {
      stats_list <- vector("list", length(simulation_results))
      skipped_simulations <- list()
      condition_skipped_csv <- file.path(intermediate_results_dir, paste0("Sim_results_N_steps_", condition_safe_name, "_skipped.csv"))
      for (simulation_index in seq_along(simulation_results)) {
        cat("Testing simulation", simulation_index, "of", length(simulation_results), "for N_steps:", N_steps_value, "\n")
        flush.console()
        simulation_data <- simulation_results[[simulation_index]]
        som_input <- list(SNP = simulation_data$SNP,
                          Morphology = simulation_data$Morphology,
                          Climate = simulation_data$Climate,
                          Host = simulation_data$Host)
        stats_list[[simulation_index]] <- tryCatch({
          run.SOM.benchmark(input_data = som_input,
                            true_labels = simulation_data$Cluster,
                            group_column_name = "N_steps",
                            group_value = N_steps_value,
                            clustering_method = clustering_method,
                            training_neighborhoods = training_neighborhoods_SOM,
                            max_k = max_k_SOM,
                            learning_rate_tuning_option = learning_rate_tuning,
                            sim_stats_row = simulation_data$sim_stats,
                            verbose = verbose_SOM)
        }, error = function(error_message) {
          skip_message <- paste0("Skipped simulation ", simulation_index, " for N_steps = ", N_steps_value, " because: ", conditionMessage(error_message))
          message(skip_message)
          skipped_simulations[[length(skipped_simulations) + 1]] <<- data.frame(N_steps = N_steps_value, sim = simulation_index, error_message = conditionMessage(error_message), stringsAsFactors = FALSE)
          NULL
        })
        save.current.condition.results(stats_list = stats_list, condition_partial_csv = condition_partial_csv, condition_partial_rds = condition_partial_rds)
        if (length(skipped_simulations) > 0) write.csv(do.call(rbind, skipped_simulations), condition_skipped_csv, row.names = FALSE)
      }
      stats_df_all <- bind.results.with.sim(stats_list)
      if (is.null(stats_df_all)) stop("All simulations failed for N_steps = ", N_steps_value)
      write.csv(stats_df_all, condition_complete_csv, row.names = FALSE)
    }
    if (!is.null(stats_df_all)) all_results_N_steps[[as.character(N_steps_value)]] <- stats_df_all
    save.partial.results(result_list = all_results_N_steps,
                         combined_csv = sim_results_N_steps_csv,
                         combined_rds = sim_results_N_steps_rds)
  }
  N_steps_SOM <- original_N_steps_SOM
  all_results_N_steps <- all_results_N_steps[vapply(all_results_N_steps, function(result_object) !is.null(result_object), logical(1))]
  if (length(all_results_N_steps) == 0) stop("No results to summarize")
  full_sim_stats_N_steps <- do.call(rbind, all_results_N_steps)
  full_sim_stats_N_steps <- as.data.frame(full_sim_stats_N_steps, stringsAsFactors = FALSE)
  write.csv(full_sim_stats_N_steps, sim_results_N_steps_csv, row.names = FALSE)
  message("N.steps simulation results saved to CSV")
}


## Summarize results
summary_sim_stats_N_steps <- full_sim_stats_N_steps %>%
  dplyr::group_by(N_steps) %>%
  dplyr::summarise(mean_ARI = safe.mean(ARI),
                   mean_K = safe.mean(K_inferred),
                   mean_Acc = safe.mean(Acc),
                   mean_Time = safe.mean(Time[!duplicated(sim)]),
                   mean_QE = safe.mean(QE),
                   mean_TE = safe.mean(TE),
                   failed_rows = sum(is.na(Time) & is.na(K_inferred) & is.na(ARI) & is.na(QE) & is.na(TE)),
                   n = dplyr::n()) %>%
  as.data.frame()
print(summary_sim_stats_N_steps)


## Plot results
N_steps_ARI_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_N_steps,
                                                                    condition_variable = "N_steps",
                                                                    response_variable = "ARI",
                                                                    condition_levels = N_steps_values,
                                                                    x_axis_label = "N.steps",
                                                                    y_axis_label = "Adjusted Rand Index",
                                                                    estimate_color = benchmark_plot_color_ARI,
                                                                    connect_estimates = TRUE)
N_steps_Acc_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_N_steps,
                                                                    condition_variable = "N_steps",
                                                                    response_variable = "Acc",
                                                                    condition_levels = N_steps_values,
                                                                    x_axis_label = "N.steps",
                                                                    y_axis_label = "Assignment accuracy",
                                                                    estimate_color = benchmark_plot_color_Acc,
                                                                    connect_estimates = TRUE)
N_steps_K_correct_model_plot <- plot.categorical.binary.benchmark.model(input_data = full_sim_stats_N_steps,
                                                                        condition_variable = "N_steps",
                                                                        response_variable = "K_correct",
                                                                        condition_levels = N_steps_values,
                                                                        x_axis_label = "N.steps",
                                                                        y_axis_label = "Probability K correct",
                                                                        estimate_color = benchmark_plot_color_K_correct,
                                                                        connect_estimates = TRUE)
N_steps_K_far_off_model_plot <- plot.categorical.binary.benchmark.model(input_data = full_sim_stats_N_steps,
                                                                        condition_variable = "N_steps",
                                                                        response_variable = "K_far_off",
                                                                        condition_levels = N_steps_values,
                                                                        x_axis_label = "N.steps",
                                                                        y_axis_label = "Probability K far off",
                                                                        estimate_color = benchmark_plot_color_K_far_off,
                                                                        connect_estimates = TRUE)
N_steps_QE_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_N_steps,
                                                                   condition_variable = "N_steps",
                                                                   response_variable = "QE",
                                                                   condition_levels = N_steps_values,
                                                                   x_axis_label = "N.steps",
                                                                   y_axis_label = "Quantization error",
                                                                   estimate_color = benchmark_plot_color_QE,
                                                                   connect_estimates = TRUE)
N_steps_TE_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_N_steps,
                                                                   condition_variable = "N_steps",
                                                                   response_variable = "TE",
                                                                   condition_levels = N_steps_values,
                                                                   x_axis_label = "N.steps",
                                                                   y_axis_label = "Topographic error",
                                                                   estimate_color = benchmark_plot_color_TE,
                                                                   connect_estimates = TRUE)
N_steps_Time_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_N_steps,
                                                                     condition_variable = "N_steps",
                                                                     response_variable = "Time",
                                                                     condition_levels = N_steps_values,
                                                                     x_axis_label = "N.steps",
                                                                     y_axis_label = "Time (seconds)",
                                                                     estimate_color = benchmark_plot_color_Time,
                                                                     connect_estimates = TRUE)
summary(N_steps_ARI_model_plot$model)
summary(N_steps_Acc_model_plot$model)
summary(N_steps_K_correct_model_plot$model)
summary(N_steps_K_far_off_model_plot$model)
summary(N_steps_QE_model_plot$model)
summary(N_steps_TE_model_plot$model)
summary(N_steps_Time_model_plot$model)
N_steps_ARI_model_plot$plot
N_steps_Acc_model_plot$plot
N_steps_K_correct_model_plot$plot
N_steps_K_far_off_model_plot$plot
N_steps_QE_model_plot$plot
N_steps_TE_model_plot$plot
N_steps_Time_model_plot$plot




#### Test effect of missing data proportion ####################################

## Set parameters
clustering_method <- "kmeans+BICelbow"
missing_data_props <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
N_clusters <- 3
training_neighborhoods_SOM <- "gaussian"
max_k_SOM <- 8


## File paths for saving/loading
sim_data_NA_file <- file.path(simulation_k3_dir, "Sim_data_NA.rds")
sim_data_NA_dir <- file.path(intermediate_results_dir, "NA_threshold_data")
if (!dir.exists(sim_data_NA_dir)) dir.create(sim_data_NA_dir, recursive = TRUE)
sim_results_NA_csv <- file.path(simulation_k3_dir, "Sim_results_NA.csv")
sim_results_NA_rds <- file.path(intermediate_results_dir, "Sim_results_NA_partial.rds")


## Load saved results or create missing-data copies/run as needed
if (file.exists(sim_results_NA_csv) && !overwrite) {
  full_sim_stats_NA <- read.csv(sim_results_NA_csv)
  message("Loaded missing-data simulation results from CSV, skipping run")
} else {
  simulation_results_all <- list()
  for (missing_prop in missing_data_props) {
    missing_prop_safe_name <- create.safe.file.name(missing_prop)
    missing_prop_data_file <- file.path(sim_data_NA_dir, paste0("Sim_data_NA_", missing_prop_safe_name, ".rds"))
    if (file.exists(missing_prop_data_file) && !overwrite) {
      simulation_results_all[[as.character(missing_prop)]] <- readRDS(missing_prop_data_file)
      message("Loaded missing-data simulation data for missing_data_prop = ", missing_prop)
    } else {
      simulation_results_all[[as.character(missing_prop)]] <- lapply(seq_along(simulation_results_base), function(simulation_index) {
        add.missing.data.to.simulation(simulation_data = simulation_results_base[[simulation_index]],
                                       missing_data_prop = as.numeric(missing_prop),
                                       missing_seed = 100000 + simulation_index)
      })
      saveRDS(simulation_results_all[[as.character(missing_prop)]], missing_prop_data_file)
      message("Saved missing-data simulation data for missing_data_prop = ", missing_prop)
    }
  }
  saveRDS(simulation_results_all, sim_data_NA_file)
  all_results_NA <- list()
  for (missing_prop in names(simulation_results_all)) {
    simulation_results <- simulation_results_all[[missing_prop]]
    missing_prop_safe_name <- create.safe.file.name(missing_prop)
    condition_complete_csv <- file.path(intermediate_results_dir, paste0("Sim_results_NA_", missing_prop_safe_name, "_complete.csv"))
    condition_partial_csv <- file.path(intermediate_results_dir, paste0("Sim_results_NA_", missing_prop_safe_name, "_partial.csv"))
    condition_partial_rds <- file.path(intermediate_results_dir, paste0("Sim_results_NA_", missing_prop_safe_name, "_partial.rds"))
    if (file.exists(condition_complete_csv) && !overwrite) {
      stats_df_all <- read.csv(condition_complete_csv)
      message("Loaded completed missing-data intermediate results for missing_data_prop = ", missing_prop)
    } else {
      stats_list <- vector("list", length(simulation_results))
      skipped_simulations <- list()
      condition_skipped_csv <- file.path(intermediate_results_dir, paste0("Sim_results_NA_", missing_prop_safe_name, "_skipped.csv"))
      for (simulation_index in seq_along(simulation_results)) {
        cat("Testing simulation", simulation_index, "of", length(simulation_results), "for missing_data_prop:", missing_prop, "\n")
        flush.console()
        simulation_data <- simulation_results[[simulation_index]]
        som_input <- list(SNP = simulation_data$SNP,
                          Morphology = simulation_data$Morphology,
                          Climate = simulation_data$Climate,
                          Host = simulation_data$Host)
        stats_list[[simulation_index]] <- tryCatch({
          run.SOM.benchmark(input_data = som_input,
                            true_labels = simulation_data$Cluster,
                            group_column_name = "missing_data_prop",
                            group_value = as.numeric(missing_prop),
                            clustering_method = clustering_method,
                            training_neighborhoods = training_neighborhoods_SOM,
                            max_k = max_k_SOM,
                            learning_rate_tuning_option = learning_rate_tuning,
                            sim_stats_row = simulation_data$sim_stats,
                            verbose = verbose_SOM)
        }, error = function(error_message) {
          skip_message <- paste0("Skipped simulation ", simulation_index, " for missing_data_prop = ", missing_prop, " because: ", conditionMessage(error_message))
          message(skip_message)
          skipped_simulations[[length(skipped_simulations) + 1]] <<- data.frame(missing_data_prop = as.numeric(missing_prop), sim = simulation_index, error_message = conditionMessage(error_message), stringsAsFactors = FALSE)
          NULL
        })
        save.current.condition.results(stats_list = stats_list, condition_partial_csv = condition_partial_csv, condition_partial_rds = condition_partial_rds)
        if (length(skipped_simulations) > 0) write.csv(do.call(rbind, skipped_simulations), condition_skipped_csv, row.names = FALSE)
      }
      stats_df_all <- bind.results.with.sim(stats_list)
      if (is.null(stats_df_all)) stop("All simulations failed for missing_data_prop = ", missing_prop)
      write.csv(stats_df_all, condition_complete_csv, row.names = FALSE)
    }
    if (!is.null(stats_df_all)) all_results_NA[[as.character(missing_prop)]] <- stats_df_all
    save.partial.results(result_list = all_results_NA,
                         combined_csv = sim_results_NA_csv,
                         combined_rds = sim_results_NA_rds)
  }
  all_results_NA <- all_results_NA[vapply(all_results_NA, function(result_object) !is.null(result_object), logical(1))]
  if (length(all_results_NA) == 0) stop("No results to summarize")
  full_sim_stats_NA <- do.call(rbind, all_results_NA)
  full_sim_stats_NA <- as.data.frame(full_sim_stats_NA, stringsAsFactors = FALSE)
  write.csv(full_sim_stats_NA, sim_results_NA_csv, row.names = FALSE)
  message("Missing-data simulation results saved to CSV")
}


## Summarize results
summary_sim_stats_NA <- full_sim_stats_NA %>%
  dplyr::group_by(missing_data_prop) %>%
  dplyr::summarise(mean_ARI = safe.mean(ARI),
                   mean_K = safe.mean(K_inferred),
                   mean_Acc = safe.mean(Acc),
                   mean_Time = safe.mean(Time[!duplicated(sim)]),
                   mean_QE = safe.mean(QE),
                   mean_TE = safe.mean(TE),
                   failed_rows = sum(is.na(Time) & is.na(K_inferred) & is.na(ARI) & is.na(QE) & is.na(TE)),
                   K1_rows = sum(K_inferred == 1, na.rm = TRUE),
                   n = dplyr::n()) %>%
  as.data.frame()
print(summary_sim_stats_NA)


## Plot results
NA_ARI_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_NA,
                                                               condition_variable = "missing_data_prop",
                                                               response_variable = "ARI",
                                                               condition_levels = missing_data_props,
                                                               x_axis_label = "Missing data proportion",
                                                               y_axis_label = "Adjusted Rand Index",
                                                               estimate_color = benchmark_plot_color_ARI,
                                                               connect_estimates = TRUE)
NA_Acc_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_NA,
                                                               condition_variable = "missing_data_prop",
                                                               response_variable = "Acc",
                                                               condition_levels = missing_data_props,
                                                               x_axis_label = "Missing data proportion",
                                                               y_axis_label = "Assignment accuracy",
                                                               estimate_color = benchmark_plot_color_Acc,
                                                               connect_estimates = TRUE)
NA_K_correct_model_plot <- plot.categorical.binary.benchmark.model(input_data = full_sim_stats_NA,
                                                                   condition_variable = "missing_data_prop",
                                                                   response_variable = "K_correct",
                                                                   condition_levels = missing_data_props,
                                                                   x_axis_label = "Missing data proportion",
                                                                   y_axis_label = "Probability K correct",
                                                                   estimate_color = benchmark_plot_color_K_correct,
                                                                   connect_estimates = TRUE)
NA_K_far_off_model_plot <- plot.categorical.binary.benchmark.model(input_data = full_sim_stats_NA,
                                                                   condition_variable = "missing_data_prop",
                                                                   response_variable = "K_far_off",
                                                                   condition_levels = missing_data_props,
                                                                   x_axis_label = "Missing data proportion",
                                                                   y_axis_label = "Probability K far off",
                                                                   estimate_color = benchmark_plot_color_K_far_off,
                                                                   connect_estimates = TRUE)
NA_QE_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_NA,
                                                              condition_variable = "missing_data_prop",
                                                              response_variable = "QE",
                                                              condition_levels = missing_data_props,
                                                              x_axis_label = "Missing data proportion",
                                                              y_axis_label = "Quantization error",
                                                              estimate_color = benchmark_plot_color_QE,
                                                              connect_estimates = TRUE)
NA_TE_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_NA,
                                                              condition_variable = "missing_data_prop",
                                                              response_variable = "TE",
                                                              condition_levels = missing_data_props,
                                                              x_axis_label = "Missing data proportion",
                                                              y_axis_label = "Topographic error",
                                                              estimate_color = benchmark_plot_color_TE,
                                                              connect_estimates = TRUE)

summary(NA_ARI_model_plot$model)
summary(NA_Acc_model_plot$model)
summary(NA_K_correct_model_plot$model)
summary(NA_K_far_off_model_plot$model)
summary(NA_QE_model_plot$model)
summary(NA_TE_model_plot$model)
NA_ARI_model_plot$plot
NA_Acc_model_plot$plot
NA_K_correct_model_plot$plot
NA_K_far_off_model_plot$plot
NA_QE_model_plot$plot
NA_TE_model_plot$plot




#### Test effect of neighborhood function ######################################

## Set neighborhood functions to test
neighborhoods <- c("gaussian", "bubble")
clustering_method <- "kmeans+BICelbow"
N_clusters <- 3
max_k_SOM <- 8


## File paths for saving/loading
sim_results_neighborhoods_csv <- file.path(simulation_k3_dir, "Sim_results_neighborhoods.csv")
sim_results_neighborhoods_rds <- file.path(intermediate_results_dir, "Sim_results_neighborhoods_partial.rds")
if (file.exists(sim_results_neighborhoods_csv)) file.remove(sim_results_neighborhoods_csv)
if (file.exists(sim_results_neighborhoods_rds)) file.remove(sim_results_neighborhoods_rds)


## Load saved results or run as needed
if (file.exists(sim_results_neighborhoods_csv) && !overwrite) {
  full_sim_stats_neighborhood <- read.csv(sim_results_neighborhoods_csv)
  message("Loaded neighborhood-function simulation results from CSV, skipping run")
} else {
  simulation_results <- simulation_results_base
  all_results_neigh <- list()
  for (neighborhood_function_name in neighborhoods) {
    cat("Running for neighborhood function:", neighborhood_function_name, "\n")
    condition_safe_name <- create.safe.file.name(neighborhood_function_name)
    condition_complete_csv <- file.path(intermediate_results_dir, paste0("Sim_results_neighborhood_", condition_safe_name, "_complete.csv"))
    condition_partial_csv <- file.path(intermediate_results_dir, paste0("Sim_results_neighborhood_", condition_safe_name, "_partial.csv"))
    condition_partial_rds <- file.path(intermediate_results_dir, paste0("Sim_results_neighborhood_", condition_safe_name, "_partial.rds"))
    if (file.exists(condition_complete_csv) && !overwrite) {
      stats_df_all <- read.csv(condition_complete_csv)
      message("Loaded completed neighborhood-function intermediate results for ", neighborhood_function_name)
    } else {
      stats_list <- vector("list", length(simulation_results))
      skipped_simulations <- list()
      condition_skipped_csv <- file.path(intermediate_results_dir, paste0("Sim_results_neighborhood_", condition_safe_name, "_skipped.csv"))
      for (simulation_index in seq_along(simulation_results)) {
        cat("Testing simulation", simulation_index, "of", length(simulation_results), "for neighborhood_function:", neighborhood_function_name, "\n")
        flush.console()
        simulation_data <- simulation_results[[simulation_index]]
        som_input <- list(SNP = simulation_data$SNP,
                          Morphology = simulation_data$Morphology,
                          Climate = simulation_data$Climate,
                          Host = simulation_data$Host)
        stats_list[[simulation_index]] <- tryCatch({
          run.SOM.benchmark(input_data = som_input,
                            true_labels = simulation_data$Cluster,
                            group_column_name = "neighborhood_function",
                            group_value = neighborhood_function_name,
                            clustering_method = clustering_method,
                            training_neighborhoods = neighborhood_function_name,
                            max_k = max_k_SOM,
                            learning_rate_tuning_option = learning_rate_tuning,
                            sim_stats_row = simulation_data$sim_stats,
                            verbose = verbose_SOM)
        }, error = function(error_message) {
          skip_message <- paste0("Skipped simulation ", simulation_index, " for neighborhood_function = ", neighborhood_function_name, " because: ", conditionMessage(error_message))
          message(skip_message)
          skipped_simulations[[length(skipped_simulations) + 1]] <<- data.frame(neighborhood_function = neighborhood_function_name, sim = simulation_index, error_message = conditionMessage(error_message), stringsAsFactors = FALSE)
          NULL
        })
        save.current.condition.results(stats_list = stats_list, condition_partial_csv = condition_partial_csv, condition_partial_rds = condition_partial_rds)
        if (length(skipped_simulations) > 0) write.csv(do.call(rbind, skipped_simulations), condition_skipped_csv, row.names = FALSE)
      }
      stats_df_all <- bind.results.with.sim(stats_list)
      if (is.null(stats_df_all)) stop("All simulations failed for neighborhood_function = ", neighborhood_function_name)
      write.csv(stats_df_all, condition_complete_csv, row.names = FALSE)
    }
    if (!is.null(stats_df_all)) all_results_neigh[[as.character(neighborhood_function_name)]] <- stats_df_all
    save.partial.results(result_list = all_results_neigh,
                         combined_csv = sim_results_neighborhoods_csv,
                         combined_rds = sim_results_neighborhoods_rds)
  }
  all_results_neigh <- all_results_neigh[vapply(all_results_neigh, function(result_object) !is.null(result_object), logical(1))]
  if (length(all_results_neigh) == 0) stop("No results to summarize")
  full_sim_stats_neighborhood <- do.call(rbind, all_results_neigh)
  full_sim_stats_neighborhood <- as.data.frame(full_sim_stats_neighborhood, stringsAsFactors = FALSE)
  write.csv(full_sim_stats_neighborhood, sim_results_neighborhoods_csv, row.names = FALSE)
  message("Neighborhood-function simulation results saved to CSV")
}


## Summarize results
summary_sim_stats_neighborhood <- full_sim_stats_neighborhood %>%
  dplyr::group_by(neighborhood_function) %>%
  dplyr::summarise(mean_ARI = safe.mean(ARI),
                   mean_K = safe.mean(K_inferred),
                   mean_Acc = safe.mean(Acc),
                   mean_Time = safe.mean(Time[!duplicated(sim)]),
                   mean_QE = safe.mean(QE),
                   mean_TE = safe.mean(TE),
                   failed_rows = sum(is.na(Time) & is.na(K_inferred) & is.na(ARI) & is.na(QE) & is.na(TE)),
                   n = dplyr::n()) %>%
  as.data.frame()
print(summary_sim_stats_neighborhood)


## Plot results
neighborhood_ARI_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_neighborhood,
                                                                         condition_variable = "neighborhood_function",
                                                                         response_variable = "ARI",
                                                                         condition_levels = neighborhoods,
                                                                         x_axis_label = "Neighborhood function",
                                                                         y_axis_label = "Adjusted Rand Index",
                                                                         estimate_color = benchmark_plot_color_ARI)
neighborhood_Acc_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_neighborhood,
                                                                         condition_variable = "neighborhood_function",
                                                                         response_variable = "Acc",
                                                                         condition_levels = neighborhoods,
                                                                         x_axis_label = "Neighborhood function",
                                                                         y_axis_label = "Assignment accuracy",
                                                                         estimate_color = benchmark_plot_color_Acc)
neighborhood_K_correct_model_plot <- plot.categorical.binary.benchmark.model(input_data = full_sim_stats_neighborhood,
                                                                             condition_variable = "neighborhood_function",
                                                                             response_variable = "K_correct",
                                                                             condition_levels = neighborhoods,
                                                                             x_axis_label = "Neighborhood function",
                                                                             y_axis_label = "Probability K correct",
                                                                             estimate_color = benchmark_plot_color_K_correct)
neighborhood_K_far_off_model_plot <- plot.categorical.binary.benchmark.model(input_data = full_sim_stats_neighborhood,
                                                                             condition_variable = "neighborhood_function",
                                                                             response_variable = "K_far_off",
                                                                             condition_levels = neighborhoods,
                                                                             x_axis_label = "Neighborhood function",
                                                                             y_axis_label = "Probability K far off",
                                                                             estimate_color = benchmark_plot_color_K_far_off)
neighborhood_QE_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_neighborhood,
                                                                        condition_variable = "neighborhood_function",
                                                                        response_variable = "QE",
                                                                        condition_levels = neighborhoods,
                                                                        x_axis_label = "Neighborhood function",
                                                                        y_axis_label = "Quantization error",
                                                                        estimate_color = benchmark_plot_color_QE)
neighborhood_TE_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_neighborhood,
                                                                        condition_variable = "neighborhood_function",
                                                                        response_variable = "TE",
                                                                        condition_levels = neighborhoods,
                                                                        x_axis_label = "Neighborhood function",
                                                                        y_axis_label = "Topographic error",
                                                                        estimate_color = benchmark_plot_color_TE)
neighborhood_Time_model_plot <- plot.categorical.gaussian.benchmark.model(input_data = full_sim_stats_neighborhood,
                                                                          condition_variable = "neighborhood_function",
                                                                          response_variable = "Time",
                                                                          condition_levels = neighborhoods,
                                                                          x_axis_label = "Neighborhood function",
                                                                          y_axis_label = "Time (seconds)",
                                                                          estimate_color = benchmark_plot_color_Time)
summary(neighborhood_ARI_model_plot$model)
summary(neighborhood_Acc_model_plot$model)
summary(neighborhood_K_correct_model_plot$model)
summary(neighborhood_K_far_off_model_plot$model)
summary(neighborhood_QE_model_plot$model)
summary(neighborhood_TE_model_plot$model)
summary(neighborhood_Time_model_plot$model)
neighborhood_ARI_model_plot$plot
neighborhood_Acc_model_plot$plot
neighborhood_K_correct_model_plot$plot
neighborhood_K_far_off_model_plot$plot
neighborhood_QE_model_plot$plot
neighborhood_TE_model_plot$plot
neighborhood_Time_model_plot$plot
