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
#### Main simulation parameters
################################################################################
overwrite <- TRUE
N_simulations <- 3
N_steps_SOM <- 10
N_replicates_SOM <- 20
max_NA_row_SOM <- 0.3
max_NA_col_SOM <- 0.6
verbose_SOM <- FALSE
BIC_threshold_SOM <- 6
learning_rate_tuning <- FALSE
if (!dir.exists("Simulations")) dir.create("Simulations")



################################################################################
#### Test effect of clustering method using simulated data
################################################################################

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
max_k_SOM <- 10


## File paths for saving/loading
sim_data_clustering_methods_file <- "Simulations/Sim_data_clustering_methods.rds"
sim_failed_clustering_methods_csv <- "Simulations/Sim_data_clustering_methods_failed_simulations.csv"
sim_results_clustering_methods_csv <- "Simulations/Sim_results_clustering_methods.csv"


## Load saved results or simulate/run as needed
if (file.exists(sim_results_clustering_methods_csv) && !overwrite) {
  full_sim_stats_clustering_methods <- read.csv(sim_results_clustering_methods_csv)
  message("Loaded clustering methods simulation results from CSV, skipping run")
} else {
  if (file.exists(sim_data_clustering_methods_file) && !overwrite) {
    simulation_results <- readRDS(sim_data_clustering_methods_file)
    message("Loaded clustering methods simulation data from RDS")
  } else {
    set.seed(1)
    simulation_seeds <- sample(1:1e8, N_simulations * 1e7)
    simulation_results <- simulate.standard.datasets(N.simulations = N_simulations,
                                                     N.clusters = N_clusters,
                                                     simulation_seeds = simulation_seeds)
    failed_simulations <- attr(simulation_results, "failed_simulations")
    if (!is.null(failed_simulations)) {
      write.csv(failed_simulations,
                sim_failed_clustering_methods_csv,
                row.names = FALSE)
    }
    saveRDS(simulation_results, sim_data_clustering_methods_file)
    message("Clustering-method simulations completed and saved to file")
  }
  all_results_clustering_methods <- list()
  for (clustering_method_name in clustering_methods) {
    cat("Running for clustering method:", clustering_method_name, "\n")
    stats_list <- vector("list", length(simulation_results))
    for (simulation_index in seq_along(simulation_results)) {
      simulation_data <- simulation_results[[simulation_index]]
      som_input <- list(SNP = simulation_data$SNP,
                        Morphology = simulation_data$Morphology,
                        Climate = simulation_data$Climate)
      stats_list[[simulation_index]] <- run.SOM.benchmark(input_data = som_input,
                                                          true_labels = simulation_data$Cluster,
                                                          group_column_name = "clustering_method",
                                                          group_value = clustering_method_name,
                                                          clustering_method = clustering_method_name,
                                                          training_neighborhoods = training_neighborhoods_SOM,
                                                          max_k = max_k_SOM,
                                                          learning_rate_tuning_option = learning_rate_tuning,
                                                          sim_stats_row = simulation_data$sim_stats,
                                                          verbose = verbose_SOM)
    }
    stats_df_all <- bind.results.with.sim(stats_list)
    all_results_clustering_methods[[clustering_method_name]] <- stats_df_all
  }
  all_results_clustering_methods <- all_results_clustering_methods[vapply(all_results_clustering_methods, function(result_object) !is.null(result_object), logical(1))]
  if (length(all_results_clustering_methods) == 0) {
    stop("No clustering-method results to summarize")
  }
  
  
  ## Combine and save all results
  full_sim_stats_clustering_methods <- do.call(rbind, all_results_clustering_methods)
  write.csv(full_sim_stats_clustering_methods, sim_results_clustering_methods_csv, row.names = FALSE)
  message("Clustering methods simulation results saved to CSV")
}


## Summarize results
summary_sim_stats_clustering_methods <- full_sim_stats_clustering_methods %>%
  dplyr::group_by(clustering_method) %>%
  dplyr::summarise(mean_ARI = mean(ARI, na.rm = TRUE),
                   mean_K = mean(K_inferred, na.rm = TRUE),
                   mean_Acc = mean(Acc, na.rm = TRUE),
                   mean_Time = mean(Time, na.rm = TRUE),
                   n = dplyr::n()) %>%
  as.data.frame()
print(summary_sim_stats_clustering_methods)


## Plot results
ggplot(full_sim_stats_clustering_methods, aes(x = factor(clustering_method), y = ARI)) +
  geom_boxplot(fill = "#abd9e9", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Clustering method", y = "Adjusted Rand Index (ARI)")
ggplot(full_sim_stats_clustering_methods, aes(x = factor(clustering_method), y = Acc)) +
  geom_boxplot(fill = "#74add1", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Clustering method", y = "Assignment accuracy")
ggplot(full_sim_stats_clustering_methods, aes(x = factor(clustering_method), y = as.numeric(K_correct))) +
  geom_boxplot(fill = "#fdae61", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Clustering method", y = "K correct (1 = correct)")
ggplot(full_sim_stats_clustering_methods, aes(x = factor(clustering_method), y = as.numeric(K_far_off))) +
  geom_boxplot(fill = "#d73027", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Clustering method", y = "K far off (1 = >=2 away from true K)")
ggplot(full_sim_stats_clustering_methods, aes(x = factor(clustering_method), y = Time)) +
  geom_boxplot(fill = "#fc9192", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Clustering method", y = "Time (seconds)")



################################################################################
#### Test effect of N.steps using simulated data
################################################################################

## Set N.steps values to test
N_steps_values <- c(20, 50, 100, 200, 500, 1000, 2000)


## Set parameters
clustering_method <- "kmeans+BICelbow"


## File paths for saving/loading
sim_data_N_steps_file <- "Simulations/Sim_data_N_steps.rds"
sim_results_N_steps_csv <- "Simulations/Sim_results_N_steps.csv"


## Load saved results or simulate/run as needed
if (file.exists(sim_results_N_steps_csv) && !overwrite) {
  full_sim_stats_N_steps <- read.csv(sim_results_N_steps_csv)
  message("Loaded N.steps simulation results from CSV, skipping run")
} else {
  if (file.exists(sim_data_N_steps_file) && !overwrite) {
    simulation_results <- readRDS(sim_data_N_steps_file)
    message("Loaded N.steps simulation data from RDS")
  } else {
    set.seed(1)
    simulation_seeds <- sample(1:1e8, N_simulations * 1e7)
    simulation_results <- simulate.standard.datasets(N.simulations = N_simulations,
                                                     N.clusters = N_clusters,
                                                     simulation_seeds = simulation_seeds)
    saveRDS(simulation_results, sim_data_N_steps_file)
    message("N.steps simulations completed and saved to file")
  }
  all_results_N_steps <- list()
  original_N_steps_SOM <- N_steps_SOM
  for (N_steps_value in N_steps_values) {
    cat("Running for N.steps:", N_steps_value, "\n")
    N_steps_SOM <- N_steps_value
    stats_list <- vector("list", length(simulation_results))
    for (simulation_index in seq_along(simulation_results)) {
      simulation_data <- simulation_results[[simulation_index]]
      som_input <- list(SNP = simulation_data$SNP,
                        Morphology = simulation_data$Morphology,
                        Climate = simulation_data$Climate)
      stats_list[[simulation_index]] <- run.SOM.benchmark(input_data = som_input,
                                                          true_labels = simulation_data$Cluster,
                                                          group_column_name = "N_steps",
                                                          group_value = N_steps_value,
                                                          clustering_method = clustering_method,
                                                          training_neighborhoods = training_neighborhoods_SOM,
                                                          max_k = max_k_SOM,
                                                          learning_rate_tuning_option = learning_rate_tuning,
                                                          sim_stats_row = simulation_data$sim_stats,
                                                          verbose = verbose_SOM)
    }
    stats_df_all <- bind.results.with.sim(stats_list)
    if (!is.null(stats_df_all)) {
      all_results_N_steps[[as.character(N_steps_value)]] <- stats_df_all
    }
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
  dplyr::summarise(mean_ARI = mean(ARI, na.rm = TRUE),
                   mean_K = mean(K_inferred, na.rm = TRUE),
                   mean_Acc = mean(Acc, na.rm = TRUE),
                   mean_Time = mean(Time, na.rm = TRUE),
                   mean_QE = mean(QE, na.rm = TRUE),
                   mean_TE = mean(TE, na.rm = TRUE),
                   n = dplyr::n()) %>%
  as.data.frame()
print(summary_sim_stats_N_steps)


## Plot results
ggplot(full_sim_stats_N_steps, aes(x = factor(N_steps), y = ARI)) +
  geom_boxplot(fill = "#abd9e9", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "N.steps", y = "Adjusted Rand Index (ARI)")
ggplot(full_sim_stats_N_steps, aes(x = factor(N_steps), y = Acc)) +
  geom_boxplot(fill = "#74add1", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "N.steps", y = "Assignment accuracy")
ggplot(full_sim_stats_N_steps, aes(x = factor(N_steps), y = as.numeric(K_correct))) +
  geom_boxplot(fill = "#fdae61", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "N.steps", y = "K correct (1 = correct)")
ggplot(full_sim_stats_N_steps, aes(x = factor(N_steps), y = as.numeric(K_far_off))) +
  geom_boxplot(fill = "#d73027", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "N.steps", y = "K far off (1 = >=2 away from true K)")
ggplot(full_sim_stats_N_steps, aes(x = factor(N_steps), y = QE)) +
  geom_boxplot(fill = "#66c2a5", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "N.steps", y = "Mean quantization error (QE)")
ggplot(full_sim_stats_N_steps, aes(x = factor(N_steps), y = TE)) +
  geom_boxplot(fill = "#fc8d62", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "N.steps", y = "Mean topographic error (TE)")
ggplot(full_sim_stats_N_steps, aes(x = factor(N_steps), y = Time)) +
  geom_boxplot(fill = "#fc9192", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "N.steps", y = "Time (seconds)")



################################################################################
#### Test effect of missing data proportion using simulated data
################################################################################

## Set missing data proportions (NA proportions)
missing_data_props <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)


## Set parameters
clustering_method <- "kmeans+BICelbow"


## File paths for saving/loading
sim_data_NA_file <- "Simulations/Sim_data_NA.rds"
sim_results_NA_csv <- "Simulations/Sim_results_NA.csv"


## Load saved results or simulate/run as needed
if (file.exists(sim_results_NA_csv) && !overwrite) {
  full_sim_stats_NA <- read.csv(sim_results_NA_csv)
  message("Loaded missing-data simulation results from CSV, skipping run")
} else {
  
  if (file.exists(sim_data_NA_file) && !overwrite) {
    simulation_results_all <- readRDS(sim_data_NA_file)
    message("Loaded missing-data simulation data from RDS")
  } else {
    set.seed(1)
    simulation_seeds <- sample(1:1e8, N_simulations * 1e7)
    
    simulation_results_all <- list()
    for (missing_prop in missing_data_props) {
      simulation_results_all[[as.character(missing_prop)]] <- simulate.standard.datasets(N.simulations = N_simulations,
                                                                                         N.clusters = N_clusters,
                                                                                         missing_data_prop = as.numeric(missing_prop),
                                                                                         simulation_seeds = simulation_seeds)
    }
    saveRDS(simulation_results_all, sim_data_NA_file)
    message("Missing-data simulations completed and saved to file")
  }
  
  all_results_NA <- list()
  
  for (missing_prop in names(simulation_results_all)) {
    simulation_results <- simulation_results_all[[missing_prop]]
    stats_list <- vector("list", length(simulation_results))
    
    for (simulation_index in seq_along(simulation_results)) {
      simulation_data <- simulation_results[[simulation_index]]
      som_input <- list(SNP = simulation_data$SNP,
                        Morphology = simulation_data$Morphology,
                        Climate = simulation_data$Climate)
      
      stats_list[[simulation_index]] <- run.SOM.benchmark(input_data = som_input,
                                                          true_labels = simulation_data$Cluster,
                                                          group_column_name = "missing_data_prop",
                                                          group_value = as.numeric(missing_prop),
                                                          clustering_method = clustering_method,
                                                          training_neighborhoods = training_neighborhoods_SOM,
                                                          max_k = max_k_SOM,
                                                          learning_rate_tuning_option = learning_rate_tuning,
                                                          sim_stats_row = simulation_data$sim_stats,
                                                          verbose = verbose_SOM)
    }
    
    stats_df_all <- bind.results.with.sim(stats_list)
    if (!is.null(stats_df_all)) {
      all_results_NA[[as.character(missing_prop)]] <- stats_df_all
    }
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
  dplyr::summarise(mean_ARI = mean(ARI, na.rm = TRUE),
                   mean_K = mean(K_inferred, na.rm = TRUE),
                   mean_Acc = mean(Acc, na.rm = TRUE),
                   mean_Time = mean(Time, na.rm = TRUE),
                   mean_QE = mean(QE, na.rm = TRUE),
                   mean_TE = mean(TE, na.rm = TRUE),
                   n = dplyr::n()) %>%
  as.data.frame()
print(summary_sim_stats_NA)


## Plot results
ggplot(full_sim_stats_NA, aes(x = factor(missing_data_prop), y = ARI)) +
  geom_boxplot(fill = "#abd9e9", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Missing data proportion", y = "Adjusted Rand Index (ARI)")
ggplot(full_sim_stats_NA, aes(x = factor(missing_data_prop), y = Acc)) +
  geom_boxplot(fill = "#74add1", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Missing data proportion", y = "Assignment accuracy")
ggplot(full_sim_stats_NA, aes(x = factor(missing_data_prop), y = as.numeric(K_correct))) +
  geom_boxplot(fill = "#fdae61", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Missing data proportion", y = "K correct (1 = correct)")
ggplot(full_sim_stats_NA, aes(x = factor(missing_data_prop), y = as.numeric(K_far_off))) +
  geom_boxplot(fill = "#d73027", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Missing data proportion", y = "K far off (1 = >=2 away from true K)")
ggplot(full_sim_stats_NA, aes(x = factor(missing_data_prop), y = QE)) +
  geom_boxplot(fill = "#66c2a5", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Missing data proportion", y = "Mean quantization error (QE)")
ggplot(full_sim_stats_NA, aes(x = factor(missing_data_prop), y = TE)) +
  geom_boxplot(fill = "#fc8d62", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Missing data proportion", y = "Mean topographic error (TE)")
ggplot(full_sim_stats_NA, aes(x = factor(missing_data_prop), y = Time)) +
  geom_boxplot(fill = "#fc9192", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Missing data proportion", y = "Time (seconds)")



################################################################################
#### Test effect of neighborhood function using simulated data
################################################################################

## Set neighborhood functions to test
neighborhoods <- c("gaussian", "bubble")


## File paths for saving/loading
sim_data_neighborhoods_file <- "Simulations/Sim_data_neighborhoods.rds"
sim_results_neighborhoods_csv <- "Simulations/Sim_results_neighborhoods.csv"


## Load saved results or simulate/run as needed
if (file.exists(sim_results_neighborhoods_csv) && !overwrite) {
  full_sim_stats_neighborhood <- read.csv(sim_results_neighborhoods_csv)
  message("Loaded neighborhood-function simulation results from CSV, skipping run")
} else {
  
  if (file.exists(sim_data_neighborhoods_file) && !overwrite) {
    simulation_results <- readRDS(sim_data_neighborhoods_file)
    message("Loaded neighborhood-function simulation data from RDS")
  } else {
    set.seed(1)
    simulation_seeds <- sample(1:1e8, N_simulations * 1e7)
    simulation_results <- simulate.standard.datasets(N.simulations = N_simulations,
                                                     N.clusters = N_clusters,
                                                     simulation_seeds = simulation_seeds)
    saveRDS(simulation_results, sim_data_neighborhoods_file)
    message("Neighborhood-function simulations completed and saved to file")
  }
  all_results_neigh <- list()
  for (neighborhood_function_name in neighborhoods) {
    stats_list <- vector("list", length(simulation_results))
    for (simulation_index in seq_along(simulation_results)) {
      simulation_data <- simulation_results[[simulation_index]]
      som_input <- list(SNP = simulation_data$SNP,
                        Morphology = simulation_data$Morphology,
                        Climate = simulation_data$Climate)
      stats_list[[simulation_index]] <- run.SOM.benchmark(input_data = som_input,
                                                          true_labels = simulation_data$Cluster,
                                                          group_column_name = "neighborhood_function",
                                                          group_value = neighborhood_function_name,
                                                          clustering_method = clustering_method,
                                                          training_neighborhoods = neighborhood_function_name,
                                                          max_k = max_k_SOM,
                                                          learning_rate_tuning_option = learning_rate_tuning,
                                                          sim_stats_row = simulation_data$sim_stats,
                                                          verbose = verbose_SOM)
    }
    stats_df_all <- bind.results.with.sim(stats_list)
    if (!is.null(stats_df_all)) {
      all_results_neigh[[as.character(neighborhood_function_name)]] <- stats_df_all
    }
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
  dplyr::summarise(mean_ARI = mean(ARI, na.rm = TRUE),
                   mean_K = mean(K_inferred, na.rm = TRUE),
                   mean_Acc = mean(Acc, na.rm = TRUE),
                   mean_Time = mean(Time, na.rm = TRUE),
                   mean_QE = mean(QE, na.rm = TRUE),
                   mean_TE = mean(TE, na.rm = TRUE),
                   n = dplyr::n()) %>%
  as.data.frame()
print(summary_sim_stats_neighborhood)


## Plot results
ggplot(full_sim_stats_neighborhood, aes(x = factor(neighborhood_function), y = ARI)) +
  geom_boxplot(fill = "#abd9e9", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Neighborhood function", y = "Adjusted Rand Index (ARI)")
ggplot(full_sim_stats_neighborhood, aes(x = factor(neighborhood_function), y = Acc)) +
  geom_boxplot(fill = "#74add1", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Neighborhood function", y = "Assignment accuracy")
ggplot(full_sim_stats_neighborhood, aes(x = factor(neighborhood_function), y = as.numeric(K_correct))) +
  geom_boxplot(fill = "#fdae61", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Neighborhood function", y = "K correct (1 = correct)")
ggplot(full_sim_stats_neighborhood, aes(x = factor(neighborhood_function), y = as.numeric(K_far_off))) +
  geom_boxplot(fill = "#d73027", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Neighborhood function", y = "K far off (1 = >=2 away from true K)")
ggplot(full_sim_stats_neighborhood, aes(x = factor(neighborhood_function), y = QE)) +
  geom_boxplot(fill = "#66c2a5", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Neighborhood function", y = "Mean quantization error (QE)")
ggplot(full_sim_stats_neighborhood, aes(x = factor(neighborhood_function), y = TE)) +
  geom_boxplot(fill = "#fc8d62", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Neighborhood function", y = "Mean topographic error (TE)")
ggplot(full_sim_stats_neighborhood, aes(x = factor(neighborhood_function), y = Time)) +
  geom_boxplot(fill = "#fc9192", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Neighborhood function", y = "Time (seconds)")



################################################################################
#### Test effect of codebook learning rate tuning using simulated data
################################################################################

## Set learning rate tuning options to test
learning_rate_tuning_SOM <- c(TRUE, FALSE)


## File paths for saving/loading
sim_data_learning_rate_tuning_file <- "Simulations/Sim_data_learning_rate_tuning.rds"
sim_results_learning_rate_tuning_csv <- "Simulations/Sim_results_learning_rate_tuning.csv"


## Load saved results or simulate/run as needed
if (file.exists(sim_results_learning_rate_tuning_csv) && !overwrite) {
  full_sim_stats_learning_rate_tuning <- read.csv(sim_results_learning_rate_tuning_csv)
  message("Loaded learning-rate-tuning simulation results from CSV, skipping run")
} else {
  
  if (file.exists(sim_data_learning_rate_tuning_file) && !overwrite) {
    simulation_results <- readRDS(sim_data_learning_rate_tuning_file)
    message("Loaded learning-rate-tuning simulation data from RDS")
  } else {
    set.seed(1)
    simulation_seeds <- sample(1:1e8, N_simulations * 1e7)
    simulation_results <- simulate.standard.datasets(N.simulations = N_simulations,
                                                     N.clusters = N_clusters,
                                                     simulation_seeds = simulation_seeds)
    saveRDS(simulation_results, sim_data_learning_rate_tuning_file)
    message("Learning-rate-tuning simulations completed and saved to file")
  }
  all_results_learning_rate_tuning <- list()
  for (learning_rate_tuning_option in learning_rate_tuning_SOM) {
    stats_list <- vector("list", length(simulation_results))
    for (simulation_index in seq_along(simulation_results)) {
      simulation_data <- simulation_results[[simulation_index]]
      som_input <- list(SNP = simulation_data$SNP,
                        Morphology = simulation_data$Morphology,
                        Climate = simulation_data$Climate)
      stats_list[[simulation_index]] <- run.SOM.benchmark(input_data = som_input,
                                                          true_labels = simulation_data$Cluster,
                                                          group_column_name = "learning_rate_tuning_codebooks",
                                                          group_value = learning_rate_tuning_option,
                                                          clustering_method = clustering_method,
                                                          training_neighborhoods = training_neighborhoods_SOM,
                                                          max_k = max_k_SOM,
                                                          learning_rate_tuning_option = learning_rate_tuning_option,
                                                          sim_stats_row = simulation_data$sim_stats,
                                                          verbose = verbose_SOM)
    }
    stats_df_all <- bind.results.with.sim(stats_list)
    if (!is.null(stats_df_all)) {
      all_results_learning_rate_tuning[[as.character(learning_rate_tuning_option)]] <- stats_df_all
    }
  }
  all_results_learning_rate_tuning <- all_results_learning_rate_tuning[vapply(all_results_learning_rate_tuning, function(result_object) !is.null(result_object), logical(1))]
  if (length(all_results_learning_rate_tuning) == 0) stop("No results to summarize")
  full_sim_stats_learning_rate_tuning <- do.call(rbind, all_results_learning_rate_tuning)
  full_sim_stats_learning_rate_tuning <- as.data.frame(full_sim_stats_learning_rate_tuning, stringsAsFactors = FALSE)
  write.csv(full_sim_stats_learning_rate_tuning, sim_results_learning_rate_tuning_csv, row.names = FALSE)
  message("Learning-rate-tuning simulation results saved to CSV")
}


## Summarize results
summary_sim_stats_learning_rate_tuning <- full_sim_stats_learning_rate_tuning %>%
  dplyr::group_by(learning_rate_tuning_codebooks) %>%
  dplyr::summarise(mean_ARI = mean(ARI, na.rm = TRUE),
                   mean_K = mean(K_inferred, na.rm = TRUE),
                   mean_Acc = mean(Acc, na.rm = TRUE),
                   mean_Time = mean(Time, na.rm = TRUE),
                   mean_QE = mean(QE, na.rm = TRUE),
                   mean_TE = mean(TE, na.rm = TRUE),
                   n = dplyr::n()) %>%
  as.data.frame()
print(summary_sim_stats_learning_rate_tuning)


## Plot results
ggplot(full_sim_stats_learning_rate_tuning, aes(x = factor(learning_rate_tuning_codebooks), y = ARI)) +
  geom_boxplot(fill = "#abd9e9", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "learning rate tuning of codebooks", y = "Adjusted Rand Index (ARI)")
ggplot(full_sim_stats_learning_rate_tuning, aes(x = factor(learning_rate_tuning_codebooks), y = Acc)) +
  geom_boxplot(fill = "#74add1", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "learning rate tuning of codebooks", y = "Assignment accuracy")
ggplot(full_sim_stats_learning_rate_tuning, aes(x = factor(learning_rate_tuning_codebooks), y = as.numeric(K_correct))) +
  geom_boxplot(fill = "#fdae61", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "learning rate tuning of codebooks", y = "K correct (1 = correct)")
ggplot(full_sim_stats_learning_rate_tuning, aes(x = factor(learning_rate_tuning_codebooks), y = as.numeric(K_far_off))) +
  geom_boxplot(fill = "#d73027", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "learning rate tuning of codebooks", y = "K far off (1 = >=2 away from true K)")
ggplot(full_sim_stats_learning_rate_tuning, aes(x = factor(learning_rate_tuning_codebooks), y = QE)) +
  geom_boxplot(fill = "#66c2a5", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "learning rate tuning of codebooks", y = "Mean quantization error (QE)")
ggplot(full_sim_stats_learning_rate_tuning, aes(x = factor(learning_rate_tuning_codebooks), y = TE)) +
  geom_boxplot(fill = "#fc8d62", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "learning rate tuning of codebooks", y = "Mean topographic error (TE)")
ggplot(full_sim_stats_learning_rate_tuning, aes(x = factor(learning_rate_tuning_codebooks), y = Time)) +
  geom_boxplot(fill = "#fc9192", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "learning rate tuning of codebooks", y = "Time (seconds)")



################################################################################
#### Test effect of multicollinearity removal with varying correlation cutoffs
################################################################################

## Set correlation cutoffs to test
#correlation_cutoffs <- c(0.6, 0.7, 0.8, 0.9, 1.0)
correlation_cutoffs <- c(0.6, 0.8)

## Set parameters
N_traits <- 20
N_trait_modules <- 4
min_traits_per_module <- 2
within_module_corr_range <- c(0.5, 0.85)
between_module_corr_range <- c(-0.15, 0.15)
correlation_jitter_SD <- 0.04
cluster_covariance_jitter_SD <- 0.03
trait_SD_range <- c(0.7, 2.0)
near_invariant_trait_prop <- 0.1
near_invariant_SD_range <- c(0.03, 0.12)
cluster_size_concentration <- 3
dataset_mean_shift_SD <- 0.75
dataset_trait_SD_jitter_SD <- 0.12
trait_mean_jitter_SD <- 0.5
trait_type_levels <- c("gaussian", "right_skewed", "left_skewed", "bounded", "count_like")
trait_type_probabilities <- c(0.35, 0.2, 0.15, 0.15, 0.15)
missing_data_prop_corr <- 0.02
outlier_prop_corr <- 0.01
outlier_SD_multiplier <- 4


## File paths
sim_data_corr_file <- "Simulations/Sim_data_corr.rds"
sim_results_corr_csv <- "Simulations/Sim_results_corr.csv"


## Create function to generate unequal cluster sizes
generate.cluster.sizes <- function(N.individuals = N_individuals,
                                   N.clusters = N_clusters,
                                   concentration = cluster_size_concentration) {
  stopifnot(is.numeric(N.individuals),
            length(N.individuals) == 1,
            !is.na(N.individuals),
            N.individuals >= N.clusters,
            N.individuals %% 1 == 0)
  stopifnot(is.numeric(N.clusters),
            length(N.clusters) == 1,
            !is.na(N.clusters),
            N.clusters >= 1,
            N.clusters %% 1 == 0)
  stopifnot(is.numeric(concentration),
            length(concentration) == 1,
            !is.na(concentration),
            concentration > 0)
  cluster_sizes <- rep(1, N.clusters)
  N_remaining_individuals <- N.individuals - sum(cluster_sizes)
  if (N_remaining_individuals > 0) {
    raw_cluster_weights <- rgamma(N.clusters, shape = concentration, rate = 1) #draw uneven cluster weights so cluster sizes differ across simulations
    extra_cluster_sizes <- as.vector(rmultinom(1, size = N_remaining_individuals, prob = raw_cluster_weights))
    cluster_sizes <- cluster_sizes + extra_cluster_sizes
  }
  cluster_sizes
}


## Create function to generate uneven module membership across traits
generate.trait.modules <- function(N.traits = N_traits,
                                   N.modules = N_trait_modules,
                                   min.traits.per.module = min_traits_per_module) {
  stopifnot(is.numeric(N.traits),
            length(N.traits) == 1,
            !is.na(N.traits),
            N.traits >= 2,
            N.traits %% 1 == 0)
  stopifnot(is.numeric(N.modules),
            length(N.modules) == 1,
            !is.na(N.modules),
            N.modules >= 1,
            N.modules %% 1 == 0,
            N.modules <= N.traits)
  stopifnot(is.numeric(min.traits.per.module),
            length(min.traits.per.module) == 1,
            !is.na(min.traits.per.module),
            min.traits.per.module >= 1,
            min.traits.per.module %% 1 == 0,
            N.modules * min.traits.per.module <= N.traits)
  module_sizes <- rep(min.traits.per.module, N.modules)
  N_remaining_traits <- N.traits - sum(module_sizes)
  if (N_remaining_traits > 0) {
    raw_module_weights <- rgamma(N.modules, shape = 1.5, rate = 1) #draw uneven module weights so module sizes differ across simulations
    extra_module_sizes <- as.vector(rmultinom(1, size = N_remaining_traits, prob = raw_module_weights))
    module_sizes <- module_sizes + extra_module_sizes
  }
  trait_module_assignments <- rep(seq_len(N.modules), times = module_sizes)
  trait_module_assignments <- sample(trait_module_assignments) #randomize module membership across traits
  
  list(trait_module_assignments = trait_module_assignments,
       module_sizes = module_sizes)
}


## Create function to generate a base correlation matrix from shared module structure
generate.base.correlation.matrix <- function(trait_module_assignments,
                                             within.module.corr.range = within_module_corr_range,
                                             between.module.corr.range = between_module_corr_range,
                                             jitter.SD = correlation_jitter_SD) {
  N.traits <- length(trait_module_assignments)
  correlation_matrix <- diag(1, N.traits)
  for (trait_i in seq_len(N.traits - 1)) {
    for (trait_j in (trait_i + 1):N.traits) {
      if (trait_module_assignments[trait_i] == trait_module_assignments[trait_j]) {
        correlation_value <- runif(1, within.module.corr.range[1], within.module.corr.range[2]) #sample stronger within-module correlations
      } else {
        correlation_value <- runif(1, between.module.corr.range[1], between.module.corr.range[2]) #sample weaker between-module correlations
      }
      correlation_value <- correlation_value + rnorm(1, mean = 0, sd = jitter.SD) #add small random variation around the module-level correlation
      correlation_value <- min(max(correlation_value, -0.99), 0.99) #constrain correlations to a realistic range
      correlation_matrix[trait_i, trait_j] <- correlation_value
      correlation_matrix[trait_j, trait_i] <- correlation_value
    }
  }
  eigenvalues <- eigen(correlation_matrix, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenvalues <= 0)) {
    correlation_matrix <- correlation_matrix + diag(abs(min(eigenvalues)) + 1e-6, N.traits) #shift diagonal to ensure positive definiteness
    correlation_matrix <- cov2cor(correlation_matrix) #rescale matrix so diagonal values equal 1
  }
  correlation_matrix
}


## Create function to jitter a base correlation matrix so covariance structure differs among clusters
jitter.correlation.matrix <- function(base.correlation.matrix,
                                      jitter.SD = cluster_covariance_jitter_SD) {
  N.traits <- ncol(base.correlation.matrix)
  jitter_matrix <- matrix(0, nrow = N.traits, ncol = N.traits)
  jitter_matrix[lower.tri(jitter_matrix)] <- rnorm(N.traits * (N.traits - 1) / 2, mean = 0, sd = jitter.SD)
  jitter_matrix <- jitter_matrix + t(jitter_matrix)
  correlation_matrix <- base.correlation.matrix + jitter_matrix
  diag(correlation_matrix) <- 1
  correlation_matrix[lower.tri(correlation_matrix)] <- pmin(pmax(correlation_matrix[lower.tri(correlation_matrix)], -0.99), 0.99)
  correlation_matrix[upper.tri(correlation_matrix)] <- t(correlation_matrix)[upper.tri(correlation_matrix)]
  eigenvalues <- eigen(correlation_matrix, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenvalues <= 0)) {
    correlation_matrix <- correlation_matrix + diag(abs(min(eigenvalues)) + 1e-6, N.traits) #shift diagonal to ensure positive definiteness
    correlation_matrix <- cov2cor(correlation_matrix) #rescale matrix so diagonal values equal 1
  }
  correlation_matrix
}


## Create function to convert correlation matrices to covariance matrices with heterogeneous trait variances
build.covariance.matrix <- function(correlation_matrix,
                                    trait_SDs) {
  correlation_matrix * (trait_SDs %o% trait_SDs)
}


## Create function to transform latent Gaussian traits into a mixture of realistic variable types
apply.trait.transforms <- function(input_matrix,
                                   trait_types) {
  input_matrix <- as.matrix(input_matrix)
  output_matrix <- input_matrix
  for (trait_index in seq_len(ncol(input_matrix))) {
    trait_values <- input_matrix[, trait_index]
    standardized_values <- if (sd(trait_values, na.rm = TRUE) > 0) {
      as.numeric(scale(trait_values))
    } else {
      rep(0, length(trait_values))
    }
    if (trait_types[trait_index] == "gaussian") {
      output_matrix[, trait_index] <- trait_values #retain approximately Gaussian traits
    } else if (trait_types[trait_index] == "right_skewed") {
      output_matrix[, trait_index] <- exp(standardized_values / 1.5) #create positive right-skewed traits
    } else if (trait_types[trait_index] == "left_skewed") {
      transformed_values <- exp(standardized_values / 1.5)
      output_matrix[, trait_index] <- max(transformed_values) + 1 - transformed_values #create bounded positive left-skewed traits
    } else if (trait_types[trait_index] == "bounded") {
      output_matrix[, trait_index] <- plogis(standardized_values) #create continuous traits bounded between 0 and 1
    } else if (trait_types[trait_index] == "count_like") {
      output_matrix[, trait_index] <- round(exp(standardized_values / 1.5)) #create non-negative count-like traits
    }
  }
  
  output_matrix
}


## Create function to add outliers and missing data
add.matrix.artifacts <- function(input_matrix,
                                 missing.prop = missing_data_prop_corr,
                                 outlier.prop = outlier_prop_corr,
                                 outlier.SD.multiplier = outlier_SD_multiplier) {
  input_matrix <- as.matrix(input_matrix)
  if (outlier.prop > 0) {
    N_outlier_cells <- floor(outlier.prop * prod(dim(input_matrix)))
    if (N_outlier_cells > 0) {
      outlier_indices <- cbind(sample(seq_len(nrow(input_matrix)), N_outlier_cells, replace = TRUE),
                               sample(seq_len(ncol(input_matrix)), N_outlier_cells, replace = TRUE))
      for (outlier_index in seq_len(nrow(outlier_indices))) {
        row_index <- outlier_indices[outlier_index, 1]
        col_index <- outlier_indices[outlier_index, 2]
        column_SD <- sd(input_matrix[, col_index], na.rm = TRUE)
        if (is.finite(column_SD) && column_SD > 0) {
          input_matrix[row_index, col_index] <- input_matrix[row_index, col_index] + rnorm(1, mean = 0, sd = column_SD * outlier.SD.multiplier) #add occasional outliers scaled to the variance of each trait
        }
      }
    }
  }
  if (missing.prop > 0) {
    non_missing_indices <- which(!is.na(input_matrix), arr.ind = TRUE)
    N_missing_cells <- floor(missing.prop * prod(dim(input_matrix)))
    if (N_missing_cells > 0 && nrow(non_missing_indices) > 0) {
      selected_missing_indices <- non_missing_indices[sample(seq_len(nrow(non_missing_indices)), min(N_missing_cells, nrow(non_missing_indices))), , drop = FALSE]
      for (missing_index in seq_len(nrow(selected_missing_indices))) {
        input_matrix[selected_missing_indices[missing_index, 1], selected_missing_indices[missing_index, 2]] <- NA #add low levels of missing data after simulation
      }
    }
  }
  
  input_matrix
}


## Create function to simulate two datasets with shared latent biology and more realistic covariance structure
simulate.two.simple.datasets <- function(N.individuals = N_individuals,
                                         N.clusters = N_clusters,
                                         N.traits = N_traits,
                                         N.modules = N_trait_modules,
                                         min.traits.per.module = min_traits_per_module,
                                         within.module.corr.range = within_module_corr_range,
                                         between.module.corr.range = between_module_corr_range,
                                         jitter.SD = correlation_jitter_SD,
                                         cluster.covariance.jitter.SD = cluster_covariance_jitter_SD,
                                         trait.SD.range = trait_SD_range,
                                         near.invariant.trait.prop = near_invariant_trait_prop,
                                         near.invariant.SD.range = near_invariant_SD_range,
                                         concentration = cluster_size_concentration,
                                         dataset.mean.shift.SD = dataset_mean_shift_SD,
                                         dataset.trait.SD.jitter.SD = dataset_trait_SD_jitter_SD,
                                         trait.mean.jitter.SD = trait_mean_jitter_SD,
                                         trait.type.levels = trait_type_levels,
                                         trait.type.probabilities = trait_type_probabilities,
                                         missing.prop = missing_data_prop_corr,
                                         outlier.prop = outlier_prop_corr,
                                         outlier.SD.multiplier = outlier_SD_multiplier) {
  cluster_sizes <- generate.cluster.sizes(N.individuals = N.individuals,
                                          N.clusters = N.clusters,
                                          concentration = concentration)
  cluster_labels <- rep(seq_len(N.clusters), times = cluster_sizes)
  individual_ids <- paste0("ID", seq_len(N.individuals))
  trait_names_1 <- paste0("Dataset1_V", seq_len(N.traits))
  trait_names_2 <- paste0("Dataset2_V", seq_len(N.traits))
  module_output <- generate.trait.modules(N.traits = N.traits,
                                          N.modules = N.modules,
                                          min.traits.per.module = min.traits.per.module)
  trait_module_assignments <- module_output$trait_module_assignments
  base_correlation_matrix <- generate.base.correlation.matrix(trait_module_assignments = trait_module_assignments,
                                                              within.module.corr.range = within.module.corr.range,
                                                              between.module.corr.range = between.module.corr.range,
                                                              jitter.SD = jitter.SD)
  shared_trait_SDs <- runif(N.traits, trait.SD.range[1], trait.SD.range[2]) #draw shared latent trait standard deviations across both datasets
  N_near_invariant_traits <- floor(near.invariant.trait.prop * N.traits)
  if (N_near_invariant_traits > 0) {
    near_invariant_trait_indices <- sample(seq_len(N.traits), N_near_invariant_traits)
    shared_trait_SDs[near_invariant_trait_indices] <- runif(N_near_invariant_traits, near.invariant.SD.range[1], near.invariant.SD.range[2]) #force a subset of traits to be near-invariant
  }
  trait_types <- sample(trait.type.levels, size = N.traits, replace = TRUE, prob = trait.type.probabilities) #simulate a mixture of Gaussian, skewed, bounded, and count-like traits
  module_baseline_means <- rnorm(N.modules, mean = 15, sd = 4)
  cluster_module_effects <- matrix(rnorm(N.clusters * N.modules, mean = 0, sd = 3), nrow = N.clusters) #simulate shared cluster-level module effects across both datasets
  shared_module_mean_matrix <- sweep(cluster_module_effects, 2, module_baseline_means, FUN = "+")
  dataset_module_mean_matrix_1 <- shared_module_mean_matrix + matrix(rnorm(N.clusters * N.modules, mean = 0, sd = dataset.mean.shift.SD), nrow = N.clusters) #add dataset-specific deviations while retaining shared latent biology
  dataset_module_mean_matrix_2 <- shared_module_mean_matrix + matrix(rnorm(N.clusters * N.modules, mean = 0, sd = dataset.mean.shift.SD), nrow = N.clusters) #add dataset-specific deviations while retaining shared latent biology
  dataset_trait_SDs_1 <- pmax(shared_trait_SDs + rnorm(N.traits, mean = 0, sd = dataset.trait.SD.jitter.SD), near.invariant.SD.range[1] / 2) #allow dataset-specific trait-scale deviations
  dataset_trait_SDs_2 <- pmax(shared_trait_SDs + rnorm(N.traits, mean = 0, sd = dataset.trait.SD.jitter.SD), near.invariant.SD.range[1] / 2) #allow dataset-specific trait-scale deviations
  dataset_1 <- matrix(NA, nrow = N.individuals, ncol = N.traits, dimnames = list(individual_ids, trait_names_1))
  dataset_2 <- matrix(NA, nrow = N.individuals, ncol = N.traits, dimnames = list(individual_ids, trait_names_2))
  for (cluster_index in seq_len(N.clusters)) {
    cluster_member_indices <- which(cluster_labels == cluster_index)
    cluster_correlation_matrix_1 <- jitter.correlation.matrix(base.correlation.matrix = base_correlation_matrix,
                                                              jitter.SD = cluster.covariance.jitter.SD)
    cluster_correlation_matrix_2 <- jitter.correlation.matrix(base.correlation.matrix = base_correlation_matrix,
                                                              jitter.SD = cluster.covariance.jitter.SD)
    cluster_trait_SDs_1 <- pmax(dataset_trait_SDs_1 + rnorm(N.traits, mean = 0, sd = 0.05), near.invariant.SD.range[1] / 2) #allow covariance structure to vary among clusters within dataset 1
    cluster_trait_SDs_2 <- pmax(dataset_trait_SDs_2 + rnorm(N.traits, mean = 0, sd = 0.05), near.invariant.SD.range[1] / 2) #allow covariance structure to vary among clusters within dataset 2
    covariance_matrix_1 <- build.covariance.matrix(correlation_matrix = cluster_correlation_matrix_1,
                                                   trait_SDs = cluster_trait_SDs_1)
    covariance_matrix_2 <- build.covariance.matrix(correlation_matrix = cluster_correlation_matrix_2,
                                                   trait_SDs = cluster_trait_SDs_2)
    cluster_mean_vector_1 <- dataset_module_mean_matrix_1[cluster_index, trait_module_assignments] + rnorm(N.traits, mean = 0, sd = trait.mean.jitter.SD) #simulate trait means from shared module-level biology with small trait-specific deviations
    cluster_mean_vector_2 <- dataset_module_mean_matrix_2[cluster_index, trait_module_assignments] + rnorm(N.traits, mean = 0, sd = trait.mean.jitter.SD) #simulate trait means from shared module-level biology with small trait-specific deviations
    dataset_1[cluster_member_indices, ] <- MASS::mvrnorm(n = length(cluster_member_indices), mu = cluster_mean_vector_1, Sigma = covariance_matrix_1) #simulate first correlated dataset
    dataset_2[cluster_member_indices, ] <- MASS::mvrnorm(n = length(cluster_member_indices), mu = cluster_mean_vector_2, Sigma = covariance_matrix_2) #simulate second correlated dataset
  }
  dataset_1 <- apply.trait.transforms(input_matrix = dataset_1,
                                      trait_types = trait_types)
  dataset_2 <- apply.trait.transforms(input_matrix = dataset_2,
                                      trait_types = trait_types)
  dataset_1 <- add.matrix.artifacts(input_matrix = dataset_1,
                                    missing.prop = missing.prop,
                                    outlier.prop = outlier.prop,
                                    outlier.SD.multiplier = outlier.SD.multiplier)
  dataset_2 <- add.matrix.artifacts(input_matrix = dataset_2,
                                    missing.prop = missing.prop,
                                    outlier.prop = outlier.prop,
                                    outlier.SD.multiplier = outlier.SD.multiplier)
  list(Dataset1 = as.data.frame(dataset_1, stringsAsFactors = FALSE),
       Dataset2 = as.data.frame(dataset_2, stringsAsFactors = FALSE),
       Cluster = cluster_labels)
}


## Load saved results or simulate/run as needed
if (file.exists(sim_results_corr_csv) && !overwrite) {
  full_sim_stats_corr <- read.csv(sim_results_corr_csv)
  message("Loaded correlation-cutoff simulation results from CSV, skipping run")
} else {
  if (file.exists(sim_data_corr_file) && !overwrite) {
    simulation_results <- readRDS(sim_data_corr_file)
    message("Loaded correlation-cutoff simulation data from RDS")
  } else {
    set.seed(1)
    simulation_seeds <- sample(1:1e8, N_simulations * 1e7)
    simulation_results <- vector("list", N_simulations)
    for (simulation_index in seq_len(N_simulations)) {
      set.seed(simulation_seeds[simulation_index])
      simulation_results[[simulation_index]] <- simulate.two.simple.datasets()
    }
    saveRDS(simulation_results, sim_data_corr_file)
    message("Correlation-cutoff simulations completed and saved to file")
  }
  all_results_corr <- list()
  for (correlation_cutoff_value in correlation_cutoffs) {
    filtered_simulation_results <- lapply(simulation_results, function(simulation_data) {
      filtered_dataset_1 <- remove.lowCV.multicollinearity.SOM(simulation_data$Dataset1,
                                                               cor.threshold = correlation_cutoff_value,
                                                               verbose = FALSE)
      filtered_dataset_2 <- remove.lowCV.multicollinearity.SOM(simulation_data$Dataset2,
                                                               cor.threshold = correlation_cutoff_value,
                                                               verbose = FALSE)
      filtered_dataset_1 <- as.data.frame(filtered_dataset_1, stringsAsFactors = FALSE)
      filtered_dataset_2 <- as.data.frame(filtered_dataset_2, stringsAsFactors = FALSE)
      if (is.null(rownames(filtered_dataset_1))) rownames(filtered_dataset_1) <- rownames(simulation_data$Dataset1)
      if (is.null(rownames(filtered_dataset_2))) rownames(filtered_dataset_2) <- rownames(simulation_data$Dataset2)
      list(Dataset1 = filtered_dataset_1,
           Dataset2 = filtered_dataset_2,
           Cluster = simulation_data$Cluster)
    })
    stats_list <- vector("list", length(filtered_simulation_results))
    for (simulation_index in seq_along(filtered_simulation_results)) {
      simulation_data <- filtered_simulation_results[[simulation_index]]
      som_input <- list(Dataset1 = simulation_data$Dataset1,
                        Dataset2 = simulation_data$Dataset2)
      stats_list[[simulation_index]] <- run.SOM.benchmark(input_data = som_input,
                                                          true_labels = simulation_data$Cluster,
                                                          group_column_name = "correlation_cutoff",
                                                          group_value = correlation_cutoff_value,
                                                          clustering_method = clustering_method,
                                                          training_neighborhoods = training_neighborhoods_SOM,
                                                          max_k = max_k_SOM,
                                                          learning_rate_tuning_option = learning_rate_tuning,
                                                          sim_stats_row = NULL,
                                                          verbose = verbose_SOM)
    }
    stats_df_all <- bind.results.with.sim(stats_list)
    if (!is.null(stats_df_all)) {
      all_results_corr[[as.character(correlation_cutoff_value)]] <- stats_df_all
    }
  }
  all_results_corr <- all_results_corr[vapply(all_results_corr, function(result_object) !is.null(result_object), logical(1))]
  if (length(all_results_corr) == 0) stop("No results to summarize")
  full_sim_stats_corr <- do.call(rbind, all_results_corr)
  write.csv(full_sim_stats_corr, sim_results_corr_csv, row.names = FALSE)
  message("Correlation-cutoff simulation results saved to CSV")
}


## Summarize results
summary_sim_stats_corr <- full_sim_stats_corr %>%
  dplyr::group_by(correlation_cutoff) %>%
  dplyr::summarise(mean_ARI = mean(ARI, na.rm = TRUE),
                   mean_K = mean(K_inferred, na.rm = TRUE),
                   mean_Acc = mean(Acc, na.rm = TRUE),
                   mean_Time = mean(Time, na.rm = TRUE),
                   mean_QE = mean(QE, na.rm = TRUE),
                   mean_TE = mean(TE, na.rm = TRUE),
                   n = dplyr::n()) %>%
  as.data.frame()
print(summary_sim_stats_corr)


## Plot results
ggplot(full_sim_stats_corr, aes(x = factor(correlation_cutoff), y = ARI)) +
  geom_boxplot(fill = "#abd9e9", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Correlation cutoff", y = "Adjusted Rand Index (ARI)")
ggplot(full_sim_stats_corr, aes(x = factor(correlation_cutoff), y = Acc)) +
  geom_boxplot(fill = "#74add1", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Correlation cutoff", y = "Assignment accuracy")
ggplot(full_sim_stats_corr, aes(x = factor(correlation_cutoff), y = as.numeric(K_correct))) +
  geom_boxplot(fill = "#fdae61", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Correlation cutoff", y = "K correct (1 = correct)")
ggplot(full_sim_stats_corr, aes(x = factor(correlation_cutoff), y = as.numeric(K_far_off))) +
  geom_boxplot(fill = "#d73027", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Correlation cutoff", y = "K far off (1 = >=2 away from true K)")
ggplot(full_sim_stats_corr, aes(x = factor(correlation_cutoff), y = QE)) +
  geom_boxplot(fill = "#66c2a5", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Correlation cutoff", y = "Mean quantization error (QE)")
ggplot(full_sim_stats_corr, aes(x = factor(correlation_cutoff), y = TE)) +
  geom_boxplot(fill = "#fc8d62", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Correlation cutoff", y = "Mean topographic error (TE)")
ggplot(full_sim_stats_corr, aes(x = factor(correlation_cutoff), y = Time)) +
  geom_boxplot(fill = "#fc9192", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Correlation cutoff", y = "Time (seconds)")