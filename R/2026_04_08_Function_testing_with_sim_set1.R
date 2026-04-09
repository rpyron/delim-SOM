################################################################################
#### Set environment and install/load packages
################################################################################

## Set environment
rm(list = ls()) #clear environment
setwd("./")
#setwd("C:/Users/danie/Desktop/PhD research/SOM package")
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")


## Install and load required R packages
required_packages <- c("clue",
                       "MASS",
                       "sn",
                       "tidyverse")
for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
  library(package_name, character.only = TRUE) #install missing packages and load all required packages
}



################################################################################
#### Simulate data
################################################################################

## Set main simulation parameters
N_individuals <- 60 #number of individuals
N_SNP_loci <- 1000 #number of SNP loci
N_morph_traits <- 15 #number of morphological traits
N_climate_variables <- 30 #number of climatic variables
N_hosts <- 12 #number of host categories

SNP_target_Fst_range <- c(0.10, 0.20) #target average Fst range across all loci (overall divergence) - range reflects moderate divergence often observed among conspecific populations (Hasbún et al. 2016 https://doi.org/10.1155/2016/3654093; Haenel et al. 2021 https://doi.org/10.1038/s41467-021-25039-y; Hall 2022 https://doi.org/10.3390/ani12091115)
SNP_differentiated_prop_range <- c(0.15, 0.4) #proportion of loci with between-cluster differentiation
SNP_differentiated_Fst_range <- c(0.85, 1) #Fst range for differentiated loci (D) with very strong structure
SNP_random_prop_range <- c(0.01, 0.06) #proportion of random loci (R) representing drift and sequencing/genotyping noise (Pompanon et al. 2005 https://doi.org/10.1038/nrg1707; Helyar et al. 2011 https://doi.org/10.1111/j.1755-0998.2010.02943.x; Mastretta-Yanes et al. 2015 https://doi.org/10.1111/1755-0998.12291)

morph_trait_distance_range <- c(1.0, 1.8) #derived distance among cluster means for differentiated morphological traits in standardized multivariate trait space (Leinonen et al. 2008 https://doi.org/10.1111/j.1420-9101.2007.01445.x; De Kort et al. 2013 https://doi.org/10.1007/s10682-012-9624-9; Siefert et al. 2015 https://doi.org/10.1111/ele.12508; Westerband et al. 2021 https://doi.org/10.1093/aob/mcab011; Opedal et al. 2023 https://doi.org/10.1073/pnas.2203228120
morph_trait_sd_range <- c(0.7, 1.2) #derived standard deviation range for morphological traits in standardized trait units (Siefert et al. 2015 https://doi.org/10.1111/ele.12508; Westerband et al. 2021 https://doi.org/10.1093/aob/mcab011

climate_variables_distance_range <- c(1.0, 1.8) #derived distance among cluster means for differentiated climate variables in standardized environmental space (Broennimann et al. 2012 https://doi.org/10.1111/j.1466-8238.2011.00698.x; Liu et al. 2020 https://doi.org/10.1073/pnas.2004289117; Bates & Bertelsmeier 2021 https://doi.org/10.1016/j.cub.2021.08.035
climate_variables_sd_range <- c(0.7, 1.2) #derived standard deviation range for climate variables in standardized environmental units (Liu et al. 2020 https://doi.org/10.1073/pnas.2004289117; Carscadden et al. 2020 https://doi.org/10.1086/710388; Bates & Bertelsmeier 2021 https://doi.org/10.1016/j.cub.2021.08.035

host_dominant_prop_range <- c(0.7, 0.96) #proportion of individuals per cluster assigned to that cluster's dominant host (Ramírez-Martínez & Tlapaya-Romero 2023 https://doi.org/10.1016/j.ijppaw.2023.05.001; Feder et al. 1994 https://doi.org/10.1073/pnas.91.17.7990; Wehmeyer et al. 2024 https://doi.org/10.1186/s13071-024-06439-7


## Create function to simulate data
simulate.data <- function(N_clusters, missing_data_prop = 0.2, sim_id = NULL) {
  
  # Print current simulation replicate
  if (!is.null(sim_id)) message("Running simulation replicate: ", sim_id)
  
  # Basic checks
  stopifnot(is.numeric(N_clusters),
            length(N_clusters) == 1,
            !is.na(N_clusters),
            N_clusters >= 1,
            N_clusters %% 1 == 0)
  stopifnot(is.numeric(missing_data_prop),
            length(missing_data_prop) == 1,
            !is.na(missing_data_prop),
            missing_data_prop >= 0,
            missing_data_prop <= 1)
  stopifnot(N_clusters - 1 <= N_morph_traits) #must have enough morphology traits to place clusters in K-1 dimensions
  stopifnot(N_clusters - 1 <= N_climate_variables) #must have enough climate variables to place clusters in K-1 dimensions
  stopifnot(N_hosts >= N_clusters) #each cluster must have a matching dominant host column
  
  # Create function to place clusters at equal distances in multivariate space using a simplex
  generate.simplex.coordinates <- function(N.points, point_distance) {
    if (N.points == 1) return(matrix(0, 1, 1)) #for a single cluster, return the origin
    coordinate_matrix <- diag(rep(point_distance, N.points - 1)) #create diagonal matrix defining simplex corners
    coordinate_matrix <- rbind(rep(0, N.points - 1), coordinate_matrix) #prepend origin so total number of points equals number of clusters
    coordinate_matrix <- scale(coordinate_matrix, scale = FALSE) #center simplex coordinates around zero
    coordinate_matrix #return coordinate matrix
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
      pairwise_Fst_matrix <- matrix(NA,
                                    nrow = N_clusters,
                                    ncol = N_clusters,
                                    dimnames = list(paste0("k", seq_len(N_clusters)), paste0("k", seq_len(N_clusters)))) #initialize pairwise Fst matrix
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
    if (!reached_target_Fst) {
      stop("After ", max_iterations,
           " iterations, observed mean pairwise Fst never fell in [",
           round(Fst_lower_bound, 3), ", ",
           round(Fst_upper_bound, 3), "]")
    }
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
    
    if (morph_differentiated_n > (N_clusters - 1)) {
      morph_cluster_mean_matrix[, N_clusters:morph_differentiated_n] <- matrix(rnorm((morph_differentiated_n - (N_clusters - 1)) * N_clusters, 10, 2), nrow = N_clusters) #fill any extra differentiated traits with random cluster means
    }
    
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
  host_matrix <- matrix(0,
                        nrow = N_individuals,
                        ncol = N_hosts,
                        dimnames = list(individual_ids, host_ids)) #initialize binary host assignment matrix
  
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
  
  # Calculate the number of different secondary host categories used per cluster
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
  
  # Add function to randomly introduce missing values
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


## Create function to extract mean QE and TE across all layers for SOM model replicate
extract.QE.TE <- function(som.model) {
  
  # Retrieve codebook matrices for each layer and corresponding input data layers
  codebook_layers <- kohonen::getCodes(som.model)
  if (!is.list(codebook_layers)) codebook_layers <- list(codebook_layers)
  codebook_layers <- lapply(codebook_layers, function(codebook_matrix) {
    codebook_matrix <- as.matrix(codebook_matrix)
    if (ncol(codebook_matrix) == 1) dim(codebook_matrix) <- c(nrow(codebook_matrix), 1) #ensure one-column codebooks remain two-dimensional
    codebook_matrix
  })
  data_layers <- som.model$data
  if (!is.list(data_layers)) data_layers <- list(as.matrix(data_layers))
  
  # Preallocate vectors to store layer-specific QE and TE values
  QE_per_layer <- numeric(length(data_layers)) #mean distance from each sample to its BMU
  TE_per_layer <- numeric(length(data_layers)) #topographic error rate per layer
  
  # Loop over each layer
  for (layer_index in seq_along(data_layers)) {
    layer_codebook <- codebook_layers[[layer_index]]
    layer_data <- data_layers[[layer_index]]
    
    # Quantization Error (QE)
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
    
    # Topographic Error (TE)
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
  
  # Return overall averages across layers
  list(
    QE = mean(QE_per_layer, na.rm = TRUE),
    TE = mean(TE_per_layer, na.rm = TRUE)
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
  if (length(simulation_seeds) < N.simulations) {
    stop("Length of simulation_seeds must be at least N.simulations")
  }
  simulation_results <- vector("list", N.simulations)
  failed_simulations <- list()
  successful_simulation_index <- 1
  seed_index <- 1
  while (successful_simulation_index <= N.simulations) {
    if (seed_index > length(simulation_seeds)) {
      stop("Ran out of simulation seeds before obtaining ", N.simulations,
           " successful simulations. Increase simulation_seeds or relax the Fst constraint.")
    }
    set.seed(simulation_seeds[seed_index])
    simulation_output <- tryCatch({
      simulate.data(N_clusters = N.clusters,
                    missing_data_prop = missing_data_prop,
                    sim_id = successful_simulation_index)
    }, error = function(error_message) {
      error_message
    })
    if (inherits(simulation_output, "error")) {
      failed_simulations[[length(failed_simulations) + 1]] <- data.frame(sim = successful_simulation_index,
                                                                         seed = simulation_seeds[seed_index],
                                                                         error = conditionMessage(simulation_output),
                                                                         stringsAsFactors = FALSE)
      seed_index <- seed_index + 1
      next
    }
    simulation_results[[successful_simulation_index]] <- simulation_output
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
  
  true_labels <- as.integer(as.factor(true_labels))
  true_K <- length(unique(true_labels))
  benchmark_output <- tryCatch({
    som_output <- NULL
    clustering_output <- NULL
    elapsed_time <- system.time({
      som_output <- train.SOM(input_data = input_data,
                              learning.rate.tuning = learning_rate_tuning_option,
                              training.neighborhoods = training_neighborhoods,
                              N.steps = N_steps_SOM,
                              N.replicates = N_replicates_SOM,
                              max.NA.row = max_NA_row_SOM,
                              max.NA.col = max_NA_col_SOM,
                              verbose = verbose)
      
      clustering_output <- clustering.SOM(SOM.output = som_output,
                                          max.k = max_k,
                                          BIC.thresh = BIC_threshold_SOM,
                                          clustering.method = clustering_method)
    })[["elapsed"]]
    list(som_models = som_output$som_models,
         clustering = clustering_output,
         elapsed = elapsed_time)
  }, error = function(error_message) {
    message("SOM/clustering failed for ", group_column_name, " = ", group_value, ": ", conditionMessage(error_message))
    NULL
  })
  if (is.null(benchmark_output) ||
      is.null(benchmark_output$som_models) ||
      length(benchmark_output$som_models) == 0 ||
      is.null(benchmark_output$clustering) ||
      is.null(benchmark_output$clustering$cluster_assignment)) {
    return(make.failed.stats.df(group_column_name = group_column_name,
                                group_value = group_value,
                                sim_stats_row = sim_stats_row))
  }
  som_models <- benchmark_output$som_models
  clustering_output <- benchmark_output$clustering
  total_time <- benchmark_output$elapsed
  time_per_replicate <- total_time / length(som_models)
  stats_per_replicate <- lapply(seq_along(som_models), function(replicate_index) {
    som_model <- som_models[[replicate_index]]
    QE_TE <- extract.QE.TE(som_model)
    tryCatch({
      predicted_labels <- get.cluster.labels.for.replicate(clustering_output$cluster_assignment, replicate_index)
      K_inferred <- length(unique(predicted_labels))
      ARI <- mclust::adjustedRandIndex(true_labels, predicted_labels)
      Acc <- if (K_inferred == true_K) get.accuracy(true_labels, predicted_labels) else NA_real_
      K_correct <- (K_inferred == true_K)
      K_far_off <- (abs(K_inferred - true_K) >= 2)
      stats_df <- data.frame(Time = time_per_replicate,
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
    }, error = function(error_message) {
      stats_df <- data.frame(Time = time_per_replicate,
                             K_inferred = NA_real_,
                             ARI = NA_real_,
                             Acc = NA_real_,
                             K_correct = NA,
                             K_far_off = NA,
                             QE = QE_TE$QE,
                             TE = QE_TE$TE,
                             stringsAsFactors = FALSE)
      stats_df[[group_column_name]] <- group_value
      stats_df <- stats_df[, c(group_column_name, setdiff(colnames(stats_df), group_column_name)), drop = FALSE]
      stats_df
    })
  })
  if (length(stats_per_replicate) == 0) {
    return(make.failed.stats.df(group_column_name = group_column_name,
                                group_value = group_value,
                                sim_stats_row = sim_stats_row))
  }
  stats_df <- do.call(rbind, stats_per_replicate)
  stats_df <- add.sim.stats(stats_df, sim_stats_row)
  stats_df
}



################################################################################
#### Evaluate and check simulated data
################################################################################
N_simulations <- 2 #number of simulations


#### k = 1

## Simulate data
N_clusters <- 1 #number of clusters to simulate
simulation_results <- list()
while (length(simulation_results) < N_simulations) {
  current_replicate <- length(simulation_results) + 1
  current_simulation <- try(simulate.data(N_clusters = N_clusters, sim_id = current_replicate), silent = TRUE)
  if (!inherits(current_simulation, "try-error")) {
    simulation_results[[current_replicate]] <- current_simulation
  }
}


## Evaluate density distributions of morphological and climate variables across clusters
simulation_1 <- simulation_results[[1]]
morphology_data <- simulation_1$Morphology %>% mutate(Cluster = factor(simulation_1$Cluster))
climate_data <- simulation_1$Climate %>% mutate(Cluster = factor(simulation_1$Cluster))
morphology_long <- morphology_data %>% pivot_longer(-Cluster, names_to = "Trait", values_to = "Value")
climate_long <- climate_data %>% pivot_longer(-Cluster, names_to = "Var", values_to = "Value")

ggplot(morphology_long, aes(x = Value, color = Cluster, fill = Cluster)) +
  geom_density(alpha = 0.5, size = 0.3) +
  facet_wrap(~Trait, scales = "free", ncol = 5) +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Density of morphological traits by true cluster",
       x = "Trait value",
       y = "Density")

ggplot(climate_long, aes(x = Value, color = Cluster, fill = Cluster)) +
  geom_density(alpha = 0.5, size = 0.3) +
  facet_wrap(~Var, scales = "free", ncol = 5) +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Density of climatic variables by true cluster",
       x = "Variable value",
       y = "Density")

## Evaluate Fst range 
fst_summary_k1 <- do.call(rbind, lapply(simulation_results, function(sim_result) sim_result$sim_stats))
fst_summary_k1$SNP_observed_Fst



#### k = 2

## Simulate data
N_clusters <- 2 #number of clusters to simulate
simulation_results <- list()
while (length(simulation_results) < N_simulations) {
  current_replicate <- length(simulation_results) + 1
  current_simulation <- try(simulate.data(N_clusters = N_clusters, sim_id = current_replicate), silent = TRUE)
  if (!inherits(current_simulation, "try-error")) {
    simulation_results[[current_replicate]] <- current_simulation
  }
}


## Evaluate one simulated dataset
SNP_sim_data <- as.data.frame(simulation_results[[1]]$SNP)
head(SNP_sim_data)
Morphology_sim_data <- as.data.frame(simulation_results[[1]]$Morphology)
head(Morphology_sim_data)
Climate_sim_data <- as.data.frame(simulation_results[[1]]$Climate)
head(Climate_sim_data)
Host_sim_data <- as.data.frame(simulation_results[[1]]$Host)
head(Host_sim_data)


## Evaluate summary statistics across all simulations
simulation_stats <- lapply(simulation_results, function(simulation_object) simulation_object$sim_stats)
simulation_stats <- do.call(rbind, simulation_stats)
head(simulation_stats)
hist(simulation_stats$SNP_observed_Fst)
mean(simulation_stats$SNP_observed_Fst)


## Evaluate density distributions of morphological and climate variables across clusters
simulation_1 <- simulation_results[[1]]
morphology_data <- simulation_1$Morphology %>% mutate(Cluster = factor(simulation_1$Cluster))
climate_data <- simulation_1$Climate %>% mutate(Cluster = factor(simulation_1$Cluster))
morphology_long <- morphology_data %>% pivot_longer(-Cluster, names_to = "Trait", values_to = "Value")
climate_long <- climate_data %>% pivot_longer(-Cluster, names_to = "Var", values_to = "Value")
ggplot(morphology_long, aes(x = Value, color = Cluster, fill = Cluster)) +
  geom_density(alpha = 0.5, size = 0.3) +
  facet_wrap(~Trait, scales = "free", ncol = 5) +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Density of morphological traits by true cluster",
       x = "Trait value",
       y = "Density")
ggplot(climate_long, aes(x = Value, color = Cluster, fill = Cluster)) +
  geom_density(alpha = 0.5, size = 0.3) +
  facet_wrap(~Var, scales = "free", ncol = 5) +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Density of climatic variables by true cluster",
       x = "Variable value",
       y = "Density")


## Evaluate host data
Host_sim_data <- as.data.frame(simulation_results[[1]]$Host)
head(Host_sim_data)
true_cluster_labels <- factor(simulation_results[[1]]$Cluster, labels = paste0("k", 1:N_clusters))
Host_sim_data$Cluster <- true_cluster_labels
host_row_sums <- rowSums(Host_sim_data[, grep("^Host_", colnames(Host_sim_data))], na.rm = FALSE)
table(host_row_sums, useNA = "ifany")
host_use_by_cluster <- Host_sim_data %>%
  pivot_longer(cols = starts_with("Host_"), names_to = "Host", values_to = "Assignment") %>%
  group_by(Cluster, Host) %>%
  summarise(prop_assigned = mean(Assignment, na.rm = TRUE), .groups = "drop")
host_use_by_cluster
ggplot(host_use_by_cluster, aes(x = Host, y = prop_assigned, fill = Cluster)) +
  geom_col(position = "dodge") +
  theme_classic() +
  labs(title = "Host assignment proportions by true cluster",
       x = "Host",
       y = "Proportion assigned") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#### k = 6

## Simulate data
N_clusters <- 6 #number of clusters to simulate
simulation_results <- list()
while (length(simulation_results) < N_simulations) {
  current_replicate <- length(simulation_results) + 1
  current_simulation <- try(simulate.data(N_clusters = N_clusters, sim_id = current_replicate, missing_data_prop = 0.4), silent = TRUE)
  if (!inherits(current_simulation, "try-error")) {
    simulation_results[[current_replicate]] <- current_simulation
  }
}


## Evaluate density distributions of morphological and climate variables across clusters
simulation_2 <- simulation_results[[2]]
morphology_data <- simulation_2$Morphology %>% mutate(Cluster = factor(simulation_2$Cluster))
climate_data <- simulation_2$Climate %>% mutate(Cluster = factor(simulation_2$Cluster))
morphology_long <- morphology_data %>% pivot_longer(-Cluster, names_to = "Trait", values_to = "Value")
climate_long <- climate_data %>% pivot_longer(-Cluster, names_to = "Var", values_to = "Value")

ggplot(morphology_long, aes(x = Value, color = Cluster, fill = Cluster)) +
  geom_density(alpha = 0.5, size = 0.3) +
  facet_wrap(~Trait, scales = "free", ncol = 5) +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Density of morphological traits by true cluster",
       x = "Trait value",
       y = "Density")
ggplot(climate_long, aes(x = Value, color = Cluster, fill = Cluster)) +
  geom_density(alpha = 0.5, size = 0.3) +
  facet_wrap(~Var, scales = "free", ncol = 5) +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Density of climatic variables by true cluster",
       x = "Variable value",
       y = "Density")


## Evaluate summary statistics across all simulations
simulation_stats <- lapply(simulation_results, function(simulation_object) simulation_object$sim_stats)
simulation_stats <- do.call(rbind, simulation_stats)
head(simulation_stats)
hist(simulation_stats$SNP_observed_Fst)
mean(simulation_stats$SNP_observed_Fst)



################################################################################
#### Evaluate and check simulated data with SOM
################################################################################

## Simulate testing data
N_clusters <- 4 #number of clusters to simulate
N_simulations <- 3 #number of simulations
NA_proportion <- 0.2

simulation_results <- list()
while (length(simulation_results) < N_simulations) {
  current_replicate <- length(simulation_results) + 1
  current_simulation <- try(simulate.data(N_clusters = N_clusters, sim_id = current_replicate, missing_data_prop = NA_proportion), silent = TRUE)
  if (!inherits(current_simulation, "try-error")) {
    simulation_results[[current_replicate]] <- current_simulation
  }
}


## Use simluation set 1 for testing
simulation_data <- simulation_results[[1]]
simulation_data_SOM <- simulation_data[c("SNP", "Morphology", "Climate", "Host")]


## Train SOM
simulation_SOM <- train.SOM(input_data = simulation_data_SOM,
                            N.steps = 200,
                            N.replicates = 100,
                            parallel = T,
                            max.NA.row = 0.5,
                            max.NA.col = 0.6,
                            save.SOM.results = F)


## Cluster codebook vectors
simulation_SOM_1 <- clustering.SOM(SOM.output = simulation_SOM, clustering.method = "kmeans+BICelbow")
simulation_SOM_2 <- clustering.SOM(SOM.output = simulation_SOM, clustering.method = "HDBSCAN")


## Evaluate results
simulation_SOM_1$optim_k_summary
simulation_SOM_2$optim_k_summary
simulation_data$sim_stats$SNP_observed_Fst
plot.model.SOM(simulation_SOM_1, replicate.mode = "representative", set.k = 4)
plot.model.SOM(simulation_SOM_2)
plot.K.SOM(simulation_SOM_1)
plot.variable.importance.SOM(simulation_SOM_1)
plot.variable.importance.SOM(simulation_SOM_1, mode = "Map.variance")
plot.variable.importance.SOM(simulation_SOM_2)



## Use simluation set 2 for testing
simulation_data <- simulation_results[[2]]
simulation_data_SOM <- simulation_data[c("SNP", "Morphology", "Climate", "Host")]


## Train SOM
simulation_SOM <- train.SOM(input_data = simulation_data_SOM,
                            N.steps = 200,
                            N.replicates = 100,
                            parallel = T,
                            max.NA.row = 0.5,
                            max.NA.col = 0.6,
                            save.SOM.results = F)


## Cluster codebook vectors
simulation_SOM_1 <- clustering.SOM(SOM.output = simulation_SOM, clustering.method = "kmeans+BICelbow")
simulation_SOM_2 <- clustering.SOM(SOM.output = simulation_SOM, clustering.method = "HDBSCAN")


## Evaluate results
simulation_SOM_1$optim_k_summary
simulation_SOM_2$optim_k_summary
simulation_data$sim_stats$SNP_observed_Fst
plot.model.SOM(simulation_SOM_1)
plot.model.SOM(simulation_SOM_2)
plot.variable.importance.SOM(simulation_SOM_1)
plot.variable.importance.SOM(simulation_SOM_2)