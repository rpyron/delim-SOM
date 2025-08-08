################################################################################
#### Set environment and install/load packages
################################################################################

## Set environment
rm(list = ls()) #clear environment
setwd("C:/Users/danie/Desktop/PhD research/SOM package")
source("delim-SOM v2 functions.R")


## Install and load required R packages
required_packages <- c("clue",
  "MASS",
  "mclust",
  "sn",
  "tidyverse")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = T)) install.packages(pkg)
  library(pkg, character.only = T) #install missing and load packages
}



################################################################################
#### Simulate data
################################################################################

## Set main simulation parameters
N_individuals <- 150 #number of individuals
N_SNP_loci <- 1000 #number of SNP loci
N_morph_traits <- 20 #number of morphological traits
N_climate_variables <- 30 #number of climatic variables
N_hosts <- 12 #number of host categories

SNP_target_Fst_range <- c(0.12, 0.25) #target average Fst range across all loci (overall divergence) - range 0.12–0.25: reflecting moderate divergence typically observed among conspecific populations (Hereford 2009 Am Nat 173)
SNP_differentiated_prop_range <- c(0.10, 0.3) #proportion of loci with between-cluster differentiation
SNP_differentiated_Fst_range <- c(0.9, 1) #Fst range for differentiated loci (D) (very high divergence)
SNP_random_prop_range <- c(0.01, 0.05) #range for proportion of random (R) drift loci (simulating drift and sequencing/genotyping error) - based on e.g. von Thaden et al. 2017 Scientific Reports 7 and Pompanon et al 2005 Nat. Rev. Genet. 6

morph_trait_distance_range <- c(2, 4) #range for distance among cluster means for differentiated morphological traits - based on empirical estimates (e.g. Chung et al. 2023 https://onlinelibrary.wiley.com/doi/10.1002/ece3.9926; Macholán 2006  https://doi.org/10.1111/j.1469-7998.2006.00156.x)
morph_trait_sd_range <- c(0.5, 1.5) #SD range for morphological traits based on empirical estimates (e.g. Chung et al. 2023 https://onlinelibrary.wiley.com/doi/10.1002/ece3.9926; Macholán 2006 https://doi.org/10.1111/j.1469-7998.2006.00156.x)

climate_variables_distance_range <- c(1.5, 3) #range for among cluster means for differentiated climate variables (e.g. Broennimann et al 2012 https://doi.org/10.1111/j.1466-8238.2011.00698.x; Peterson et al 1999 https://doi.org/10.1126/science.285.5431.1265)
climate_variables_sd_range <- c(0.3, 1.3) #SD range for climate variables (e.g. Title & Bemmels 2018 https://doi.org/10.1111/ecog.02880; Peterson 2011 https://doi.org/10.1111/j.1365-2699.2010.02456.x)

host_dominant_prop_range <- c(0.7, 0.95) #range of proportion for number of individuals per cluster assigned to dominant host (reasonable range based on empirical estimates, e.g., Martínez & Tlapaya 2023 10.1016/j.ijppaw.2023.05.001; Wehmeyer et al 2024 https://doi.org/10.1186/s13071-024-06439-7)


## Create function to simulate data
simulate.data <- function(missing_data_prop = 0.1) {
  
  # Create function to separate clusters with equal distances in multivariate space (i.e., compute K points on (K−1)-dimensional simplex scaled by d)
  simplex_coords <- function(k, d) {
    if (k == 1) return(matrix(0, 1, 1)) #if k = 1, return origin
    mat <- diag(rep(d, k - 1)) #build (K−1) × (K−1) diagonal matrix with entries = d
    mat <- rbind(rep(0, k - 1), mat) #prepend zero row to get K points total
    mat <- scale(mat, scale = FALSE) #center each column to center simplex
    mat #return K × (K−1) coordinate matrix
  }
  
  # Assign individuals to clusters (as evenly as possible)
  base_n <- floor(N_individuals / N_clusters)
  remainder <- N_individuals %% N_clusters
  samples_per_cluster <- rep(base_n, N_clusters)
  if (remainder > 0) samples_per_cluster[1:remainder] <- samples_per_cluster[1:remainder] + 1
  cluster_vector <- unlist(mapply(rep, seq_len(N_clusters), samples_per_cluster))
  sample_ids <- paste0("ID", seq_len(N_individuals)) #set sample IDs
  
  # SNP data
  # Simulate SNP genotypes for all individuals using three types of loci:
  # 1. Differentiated loci (D): structured among clusters (high Fst, strong signal, low proportion) - each cluster has its own MAF for each locus drawn from beta distribution with parameters a, b so that mean pairwise Fst between clusters matches target Fst range
  # 2. Neutral loci (N): same distribution across all clusters (unstructured, high proportion)
  # 3. Random loci (R): random flipping of genotypes to simulate drift and genotyping/sequencing error (no structure, very low proportion)
  SNP_differentiated_prop <- runif(1, SNP_differentiated_prop_range[1], SNP_differentiated_prop_range[2]) #randomly draw proportion of differentiated loci
  SNP_random_prop <- runif(1, SNP_random_prop_range[1], SNP_random_prop_range[2]) #randomly draw proportion of random drift loci
  SNP_differentiated_n <- round(N_SNP_loci * SNP_differentiated_prop) #calculate number of differentiated loci
  SNP_random_n <- round(N_SNP_loci * SNP_random_prop) #calculate number of drift loci
  SNP_neutral_n <- N_SNP_loci - SNP_differentiated_n - SNP_random_n #calculate number of neutral loci
  
  SNP_locus_types <- c(rep("D", SNP_differentiated_n), rep("N", SNP_neutral_n), rep("R", SNP_random_n)) #assign locus types
  SNP_locus_types <- sample(SNP_locus_types) #shuffle
  SNP_colnames <- sprintf("SNP_%s_%d", SNP_locus_types, as.integer(ave(SNP_locus_types, SNP_locus_types, FUN = seq_along))) #build locus names
  
  target_Fst <- runif(1, SNP_target_Fst_range[1], SNP_target_Fst_range[2]) #randomly draw overall target Fst
  SNP_matrix <- matrix(NA, nrow = N_individuals, ncol = N_SNP_loci, dimnames = list(sample_ids, SNP_colnames)) #create matrix for simulated SNP genotypes (individuals x loci)
  
  for (i in which(SNP_locus_types == "N")) { #for neutral loci N (shared distribution)
    allele_freq <- rbeta(1, 0.5, 0.5) #draw MAF from beta distribution for U‑shaped spectrum
    SNP_matrix[, i] <- rbinom(N_individuals, 2, allele_freq) #binomially sample genotypes (0,1,2)
  }
  
  for (i in which(SNP_locus_types == "R")) { #for random drift loci R (not structured)
    drift_maf <- rbeta(1, 0.5, 0.5) #draw MAF from beta distribution for U‑shaped spectrum (no structure)
    geno <- rbinom(N_individuals, 2, drift_maf)
    n_flip <- floor(N_individuals * runif(1, SNP_random_prop_range[1], SNP_random_prop_range[2])) #fraction of samples to randomize
    if (n_flip > 0) { #randomly flip some genotypes to simulate drift/genotyping errors
      flip_idx <- sample(N_individuals, n_flip)
      geno[flip_idx] <- sample(0:2, n_flip, replace = TRUE)
    }
    SNP_matrix[, i] <- geno
  }
  
  if (N_clusters > 1) {
    converged <- FALSE
    max_iter  <- 500 #max number of iterations
    for (iter in seq_len(max_iter)) { #repeat until mean Fst across all loci matches target Fst range (max 5000 times)
      for (i in which(SNP_locus_types == "D")) { #for differentiated loci D (cluster-structured)
        locus_target_Fst <- runif(1, SNP_differentiated_Fst_range[1], SNP_differentiated_Fst_range[2]) #randomly draw high Fst for D loci
        global_maf <- runif(1, 0.2, 0.8) #global minor allele frequency across all clusters
        a <- global_maf * (1 - locus_target_Fst) / (locus_target_Fst * 0.8) #beta parameters: a/b shape how different allele frequencies can be among clusters
        b <- (1 - global_maf) * (1 - locus_target_Fst) / (locus_target_Fst * 0.8) 
        cluster_mafs <- rbeta(N_clusters, a, b) #draw cluster-specific allele frequencies from beta distribution
        for (g in seq_len(N_clusters)) {
          inds <- which(cluster_vector == g)
          SNP_matrix[inds, i] <- rbinom(length(inds), 2, cluster_mafs[g])
        }
      }
      pairwise_fst_matrix <- matrix(NA, #prepare empty K × K matrix
                                    nrow = N_clusters,
                                    ncol = N_clusters,
                                    dimnames = list(paste0("k", 1:N_clusters), paste0("k", 1:N_clusters)
                                    )
      ) 
      for (cluster_i in seq_len(N_clusters - 1)) { #iterate clusters i
        for (cluster_j in (cluster_i + 1):N_clusters) { #iterate clusters j>i
          sel <- cluster_vector %in% c(cluster_i, cluster_j) #select individuals in i or j
          fst_vals <- sapply(seq_len(N_SNP_loci), function(locus) { #per-locus FST
            freqs <- tapply(SNP_matrix[sel, locus] / 2, cluster_vector[sel], mean) #allele freq in each pop
            p_bar <- mean(freqs) #mean allele freq
            if(p_bar * (1 - p_bar) > 0) var(freqs) / (p_bar * (1 - p_bar)) else 0 #per-locus FST
          })
          pairwise_fst_matrix[cluster_i, cluster_j] <- mean(fst_vals, na.rm = TRUE) #average FST
          pairwise_fst_matrix[cluster_j, cluster_i] <- pairwise_fst_matrix[cluster_i, cluster_j] #mirror
        }
      }
      fst_vals <- pairwise_fst_matrix[lower.tri(pairwise_fst_matrix)] #extract lower triangle values
      lower_bound <- SNP_target_Fst_range[1] #target Fst lower bound
      upper_bound <- SNP_target_Fst_range[2] #target Fst upper bound
      tolerance <- 0.03 #allowed tolerance outside target range
      in_range <- (fst_vals >= (lower_bound - tolerance)) & (fst_vals <= (upper_bound + tolerance)) #check if within tolerance
      prop_in_range <- mean(in_range) #proportion in tolerance range
      if (prop_in_range >= 0.75) { #accept if ≥75% pairs in range
        converged <- TRUE
        break #exit iteration loop
      }
    }
    if (!converged) {
      stop("After ", max_iter,
           " iterations, observed Fst never fell in [",
           SNP_target_Fst_range[1], ", ",
           SNP_target_Fst_range[2], "]")
    }
    observed_Fst <- mean(pairwise_fst_matrix[lower.tri(pairwise_fst_matrix)], na.rm = TRUE) #calculate observed average pairwise FST
  } else { #for k = 1: treat all D‐loci as unstructured
    for (i in which(SNP_locus_types == "D")) {
      allele_freq <- runif(1, 0.15, 0.5)
      SNP_matrix[, i] <- rbinom(N_individuals, 2, allele_freq)
    }
    observed_Fst <- 0
  }
  
  SNP_df <- as.data.frame(SNP_matrix)
  SNP_df$Cluster <- factor(cluster_vector, labels = paste0("k", seq_len(N_clusters)))
  rownames(SNP_df) <- sample_ids
  
  # Morphological variables
  # Simulate mixture of neutral (shared) and differentiated (cluster-specific) traits:
  # 1) Neutral traits (N) are sampled from same normal distribution for all clusters
  # 2) Differentiated traits (D) have cluster-specific means arranged on simplex in trait space (each cluster is positioned at different corner in multivariate space with equal distance between clusters)
  morph_differentiated_n <- N_clusters - 1
  morph_neutral_n <- N_morph_traits - morph_differentiated_n
  stopifnot(morph_differentiated_n + morph_neutral_n == N_morph_traits)
  morph_trait_types <- c(rep("D", morph_differentiated_n), rep("N", morph_neutral_n))
  morph_trait_types <- sample(morph_trait_types)
  D_cols <- which(morph_trait_types == "D")
  morph_colnames <- sprintf("Trait_%s_%d", morph_trait_types, as.integer(ave(morph_trait_types, morph_trait_types, FUN = seq_along)))
  morph_matrix <- matrix(NA, nrow = N_individuals, ncol = N_morph_traits, dimnames = list(sample_ids, morph_colnames))
  
  morph_trait_sd <- runif(1, morph_trait_sd_range[1], morph_trait_sd_range[2]) #draw SD value randomly from range
  for (i in which(morph_trait_types == "N")) { #simulate neutral morphological traits
    morph_matrix[, i] <- rnorm(N_individuals, mean = 0, sd = morph_trait_sd)
  }
  
  if (morph_differentiated_n > 0) { #simulate differentiated morph traits (clusters have different means arranged on simplex)
    means_mat <- matrix(0, N_clusters, morph_differentiated_n)
    if (N_clusters == 2) {
      morph_trait_distance <- runif(1, morph_trait_distance_range[1], morph_trait_distance_range[2])
      means_mat[, 1] <- c(-morph_trait_distance/2, morph_trait_distance/2)
    } else if (morph_differentiated_n >= (N_clusters - 1)) {
      morph_trait_distance <- runif(1, morph_trait_distance_range[1], morph_trait_distance_range[2])
      simple_means <- simplex_coords(N_clusters, morph_trait_distance)
      means_mat[, 1:(N_clusters - 1)] <- simple_means
    }
    if (morph_differentiated_n > (N_clusters - 1)) {
      means_mat[, (N_clusters):morph_differentiated_n] <- matrix(rnorm((morph_differentiated_n - (N_clusters - 1)) * N_clusters, 10, 2), nrow = N_clusters)
    }
    D_cols <- which(morph_trait_types == "D")
    for (j in seq_along(D_cols)) {
      col_idx <- D_cols[j]
      for (g in seq_len(N_clusters)) {
        inds <- which(cluster_vector == g)
        morph_matrix[inds, col_idx] <- rnorm(length(inds), mean = means_mat[g, j], sd = morph_trait_sd)
      }
    }
  }
  
  morph_df <- as.data.frame(morph_matrix)
  rownames(morph_df) <- sample_ids
  
  # Climate variables
  # Simulate environmental variables similar to morphological data but:
  # All variables have skewed (asymmetric) normal distribution with to capture potential real-world non-normality
  climate_differentiated_n <- N_clusters - 1 #number of D variables
  climate_neutral_n <- N_climate_variables - climate_differentiated_n #number of N variables
  climate_variables_types <- c(rep("D", climate_differentiated_n), rep("N", climate_neutral_n)) #build types vector
  climate_variables_types <- sample(climate_variables_types) #shuffle types
  C_cols <- which(climate_variables_types == "D") #indices of D variables
  N_cols <- which(climate_variables_types == "N") #indices of N variables
  climate_colnames <- sprintf("Climate_%s_%d", climate_variables_types, as.integer(ave(climate_variables_types, climate_variables_types, FUN = seq_along))) #variable names
  climate_matrix <- matrix(NA, nrow = N_individuals, ncol = N_climate_variables, dimnames = list(sample_ids, climate_colnames)) #empty matrix
  climate_variables_sd <- runif(1, climate_variables_sd_range[1], climate_variables_sd_range[2]) #draw SD value randomly
  
  climate_means_mat <- matrix(0, nrow = N_clusters, ncol = N_climate_variables) #initialize all means to 0
  if (climate_differentiated_n > 0) { #set differentiated variable means for clusters
    if (N_clusters == 2) {
      climate_variables_distance <- runif(1, climate_variables_distance_range[1], climate_variables_distance_range[2]) #draw mean separation
      climate_means_mat[, C_cols[1]] <- c(-climate_variables_distance/2, climate_variables_distance/2) #means for k=2
    } else if (climate_differentiated_n >= (N_clusters - 1)) {
      climate_variables_distance <- runif(1, climate_variables_distance_range[1], climate_variables_distance_range[2]) #draw mean separation
      mat <- diag(rep(climate_variables_distance, N_clusters - 1)) #simplex
      mat <- rbind(rep(0, N_clusters - 1), mat) #prepend zero row
      simple_means <- scale(mat, scale = FALSE) #center
      climate_means_mat[, C_cols[1:(N_clusters - 1)]] <- simple_means #assign means for D variables
    }
    if (climate_differentiated_n > (N_clusters - 1)) { #if more D variables, random means for rest
      idx <- (N_clusters):climate_differentiated_n
      climate_means_mat[, C_cols[idx]] <- matrix(rnorm(length(idx) * N_clusters, 0, 2), nrow = N_clusters)
    }
  }
  
  cluster_skews_d <- matrix(runif(N_clusters * climate_differentiated_n, -6, 6), nrow = N_clusters) #random skew for D cluster-variables
  neutral_skews <- runif(length(N_cols), -6, 6) #random skew per neutral variable (same across all clusters)
  
  for (j in seq_len(N_climate_variables)) { #iterate variables
    if (j %in% C_cols) { #differentiated variables
      for (g in seq_len(N_clusters)) { #iterate clusters
        inds <- which(cluster_vector == g) #sample indices for this cluster
        mean_val <- climate_means_mat[g, j] #mean for this cluster-variable
        skew_val <- cluster_skews_d[g, which(C_cols == j)] #skewness for this cluster-variable
        climate_matrix[inds, j] <- sn::rsn(length(inds), xi = mean_val, omega = climate_variables_sd, alpha = skew_val) #simulate skewed normal
      }
    } else { #neutral variables (same distribution for all clusters, skew varies by variable)
      idx <- which(N_cols == j) #index of neutral variable in neutral variables list
      climate_matrix[, j] <- sn::rsn(N_individuals, xi = 0, omega = climate_variables_sd, alpha = neutral_skews[idx])
    }
  }
  
  climate_df <- as.data.frame(climate_matrix)
  rownames(climate_df) <- sample_ids
  
  # Host data:
  # For each cluster, draw dominant‐host proportion of its members and assign them to cluster’s primary host category
  # For remaining individuals, draw secondary-host count via zero-truncated Poisson and select those hosts, then allocate each leftover individual among them using weights drawn from Gamma distribution (to induce skew)
  host_ids <- paste0("Host_", seq_len(N_hosts))
  dominant_host_prop <- runif(1, host_dominant_prop_range[1], host_dominant_prop_range[2]) #draw dominant host proportion (per simulation)
  host_matrix <- matrix(0, #prepare empty matrix
                        nrow = N_individuals,
                        ncol = N_hosts,
                        dimnames = list(sample_ids, host_ids))
  
  for(g in seq_len(N_clusters)) { #iterate clusters
    inds_g <- which(cluster_vector == g) #indices of individuals in cluster g
    n_g <- length(inds_g) #cluster size
    n_main <- round(dominant_host_prop * n_g) #primary host count
    main_i <- sample(inds_g, n_main) #pick individuals for primary host
    host_matrix[main_i, g] <- 1 #assign primary host
    other_i <- setdiff(inds_g, main_i) #remaining individuals to assign
    if (N_clusters > 1) { #only run secondary host logic if more than one cluster
      k_sec <- rpois(1, lambda = 2) + 1 #draw number of secondary hosts from zero-truncated Poisson
      k_sec <- min(k_sec, N_hosts - 1) #cap at N_hosts-1
      avail_hosts <- setdiff(seq_len(N_hosts), g) #exclude primary host
      actual_k_sec <- min(k_sec, length(avail_hosts)) #number of possible secondary hosts
      sec_hosts <- if (actual_k_sec > 0) sample(avail_hosts, actual_k_sec, replace = FALSE) else integer(0) #sample secondary hosts
      if (length(sec_hosts) > 0 && length(other_i) > 0) { #check for non-empty assignments
        if (length(sec_hosts) == 1) { #if only one possible secondary host, assign all directly
          host_matrix[other_i, sec_hosts] <- 1
        } else {
          pvec_raw <- rgamma(length(sec_hosts), shape = 0.8, rate = 1) #gamma weights for secondary hosts
          if (all(!is.na(pvec_raw)) && all(pvec_raw > 0) && sum(pvec_raw) > 0 && length(pvec_raw) == length(sec_hosts)) { #probabilities must be valid and positive
            pvec <- pvec_raw / sum(pvec_raw) #normalize weights
            if (length(pvec) == length(sec_hosts) && all(!is.na(pvec)) && all(pvec > 0) && length(other_i) > 0) { #final check for valid probabilities
              host_choice <- sample(sec_hosts, length(other_i), prob = pvec, replace = TRUE) #assign each leftover individual to secondary host
              for (h in sec_hosts) host_matrix[other_i[host_choice == h], h] <- 1 #update host matrix for secondary hosts
            } else {
              fallback <- rep(sec_hosts, length.out = length(other_i)) #fallback: assign evenly
              host_matrix[cbind(other_i, fallback)] <- 1
            }
          } else {
            fallback <- rep(sec_hosts, length.out = length(other_i)) #fallback: assign evenly
            host_matrix[cbind(other_i, fallback)] <- 1
          }
        }
      }
    }
  }
  
  # Calculate number of different secondary host species per cluster
  secondary_host_counts <- sapply(seq_len(N_clusters), function(g) { #for each cluster g
    inds_g <- which(cluster_vector == g) #indices of individuals in cluster g
    primary_host_col <- g #primary host column index for cluster g
    hosts_in_cluster <- host_matrix[inds_g, , drop = FALSE] #host assignments for cluster g individuals
    secondary_hosts <- hosts_in_cluster[, -primary_host_col, drop = FALSE] #exclude primary host column
    sum(colSums(secondary_hosts) > 0) #count secondary host species with any individuals assigned
  })
  secondary_host_counts <- mean(secondary_host_counts) #take mean
  
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
    SNP_target_Fst = target_Fst,
    SNP_observed_Fst = observed_Fst,
    morph_differentiated_n = morph_differentiated_n,
    morph_neutral_n = morph_neutral_n,
    climate_differentiated_n = climate_differentiated_n,
    climate_neutral_n = climate_neutral_n,
    dominant_host_prop = dominant_host_prop,
    secondary_host_species_N = secondary_host_counts,
    stringsAsFactors = FALSE
  )
  
  # Add function to add missing data
  add.NA <- function(mat, prop) {
    mat <- as.matrix(mat)
    idx <- which(!is.na(mat), arr.ind = TRUE)
    n_miss <- floor(prop * prod(dim(mat)))
    if (n_miss > 0) {
      miss_idx <- idx[sample(seq_len(nrow(idx)), n_miss), , drop = FALSE]
      for (i in seq_len(nrow(miss_idx))) {
        mat[miss_idx[i, 1], miss_idx[i, 2]] <- NA
      }
    }
    mat
  }
  
  # Add missing data
  SNP_df_missing <- as.data.frame(add.NA(SNP_matrix, missing_data_prop))
  morph_df_missing <- as.data.frame(add.NA(morph_matrix, missing_data_prop))
  climate_df_missing <- as.data.frame(add.NA(climate_matrix, missing_data_prop))
  host_df_missing <- as.data.frame(add.NA(host_matrix, missing_data_prop))
  rownames(SNP_df_missing) <- sample_ids
  rownames(morph_df_missing) <- sample_ids
  rownames(climate_df_missing) <- sample_ids
  rownames(host_df_missing) <- sample_ids
  
  # Return results
  list(SNP = SNP_df_missing,
       Morphology = morph_df_missing,
       Climate = climate_df_missing,
       Host = host_df_missing,
       sim_stats = sim_stats,
       Cluster = cluster_vector)
}


## Create function to extract mean QE and TE across all layers for SOM model replicate
extract.QE.TE <- function(som.model) {
  
  # Retrieve list of codebook matrices (one per layer) and corresponding input data layers
  codebook.layers <- lapply(kohonen::getCodes(som.model), function(x) {
    m <- as.matrix(x)
    # if it collapsed to a single vector, enforce 2D
    if (ncol(m) == 1) dim(m) <- c(nrow(m), 1)
    m
  })
  data.layers <- som.model$data
  
  # Preallocate vectors to store per-layer QE and TE values
  QE.per.layer <- numeric(length(data.layers)) # mean distance to BMU
  TE.per.layer <- numeric(length(data.layers)) # topographic error rate
  
  # Loop over each layer
  for (layer.idx in seq_along(data.layers)) {
    layer.codebook <- codebook.layers[[layer.idx]]
    layer.data     <- data.layers[[layer.idx]]
    
    # Quantization Error (QE)
    sample.QE <- numeric(nrow(layer.data))
    for (s in seq_len(nrow(layer.data))) {
      dims <- which(!is.na(layer.data[s, ]))
      if (length(dims) > 0) {
        bmu.index <- som.model$unit.classif[s]
        sqdiffs   <- (layer.data[s, dims] - layer.codebook[bmu.index, dims])^2
        sample.QE[s] <- sqrt(sum(sqdiffs))
      } else {
        sample.QE[s] <- NA
      }
    }
    QE.per.layer[layer.idx] <- mean(sample.QE, na.rm = TRUE)
    
    # Topographic Error (TE)
    unit.dist.mat <- kohonen::unit.distances(som.model$grid)
    sample.TE     <- logical(nrow(layer.data))
    for (s in seq_len(nrow(layer.data))) {
      dims <- which(!is.na(layer.data[s, ]))
      if (length(dims) < 2) {
        sample.TE[s] <- NA
        next
      }
      diff.mat       <- t(layer.codebook[, dims]) - layer.data[s, dims]
      sq.dists       <- colSums(diff.mat^2)
      ordered.idx    <- order(sq.dists)
      sample.TE[s]   <- unit.dist.mat[ordered.idx[1], ordered.idx[2]] > 1
    }
    TE.per.layer[layer.idx] <- mean(sample.TE, na.rm = TRUE)
  }
  
  # Return overall averages
  list(
    QE = mean(QE.per.layer, na.rm = TRUE),
    TE = mean(TE.per.layer, na.rm = TRUE)
  )
}


## Create function to map predicted clusters to true clusters using Hungarian algorithm for max agreement
get.best.assignment <- function(true_labels, pred_labels) {
  t <- table(true_labels, pred_labels)
  assignment <- clue::solve_LSAP(t, maximum = TRUE)
  pred_mapped <- assignment[match(pred_labels, as.integer(colnames(t)))]
  return(as.integer(rownames(t))[pred_mapped])
}


## Create function to calculate proportion correctly assigned, given best matching
get.accuracy <- function(true_labels, pred_labels) {
  if (length(unique(pred_labels)) != length(unique(true_labels))) return(NA)
  tab <- table(true_labels, pred_labels)
  assignment <- clue::solve_LSAP(tab, maximum = TRUE)
  mapped_pred <- as.integer(rownames(tab))[assignment[match(pred_labels, as.integer(colnames(tab)))]]
  mean(true_labels == mapped_pred)
}



################################################################################
#### Evaluate and check simulated data
################################################################################

N_simulations <- 5 #number of simulations


#### k = 1  

## Simulate data
N_clusters <- 1 #number of K to simulate
sim_results <- list()
while (length(sim_results) < N_simulations) {
  sim_result <- try(simulate.data(), silent = TRUE)
  if (!inherits(sim_result, "try-error")) {
    sim_results[[length(sim_results) + 1]] <- sim_result
  }
}


## Evaluate density distributions of morphological and climate variables across clusters
sim1 <- sim_results[[1]]
morph <- sim1$Morph %>% mutate(Cluster = factor(sim1$Cluster))
climate <- sim1$Climate %>% mutate(Cluster = factor(sim1$Cluster))
morph_long <- morph %>% pivot_longer(-Cluster, names_to = "Trait",  values_to = "Value")
climate_long <- climate %>% pivot_longer(-Cluster, names_to = "Var", values_to = "Value")

ggplot(morph_long, aes(x = Value, color = Cluster, fill = Cluster)) +
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


#### k = 2  

## Simulate data
N_clusters <- 2 #number of K to simulate
sim_results <- list()
while (length(sim_results) < N_simulations) {
  sim_result <- try(simulate.data(), silent = TRUE)
  if (!inherits(sim_result, "try-error")) {
    sim_results[[length(sim_results) + 1]] <- sim_result
  }
}

## Evaluate one simulated dataset
SNP_sim_data <- as.data.frame(sim_results[[1]]$SNP)
head(SNP_sim_data)
Morph_sim_data <- as.data.frame(sim_results[[1]]$Morph)
head(Morph_sim_data)
Clim_sim_data <- as.data.frame(sim_results[[1]]$Climate)
head(Clim_sim_data)
Host_sim_data <- as.data.frame(sim_results[[1]]$Host)
head(Host_sim_data)


## Evaluate summary statistics across all simulations
simulation_stats <- lapply(sim_results, function(x) x$sim_stats)
simulation_stats <- do.call(rbind, simulation_stats)
head(simulation_stats)


## Evaluate density distributions of morphological and climate variables across clusters
sim1 <- sim_results[[1]]
morph <- sim1$Morph %>% mutate(Cluster = factor(sim1$Cluster))
climate <- sim1$Climate %>% mutate(Cluster = factor(sim1$Cluster))
morph_long <- morph %>% pivot_longer(-Cluster, names_to = "Trait",  values_to = "Value")
climate_long <- climate %>% pivot_longer(-Cluster, names_to = "Var", values_to = "Value")

ggplot(morph_long, aes(x = Value, color = Cluster, fill = Cluster)) +
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


#### k = 8 

## Simulate data
N_clusters <- 8 #number of K to simulate
sim_results <- list()
while (length(sim_results) < N_simulations) {
  sim_result <- try(simulate.data(), silent = TRUE)
  if (!inherits(sim_result, "try-error")) {
    sim_results[[length(sim_results) + 1]] <- sim_result
  }
}

## Evaluate density distributions of morphological and climate variables across clusters
sim1 <- sim_results[[2]]
morph <- sim1$Morph %>% mutate(Cluster = factor(sim1$Cluster))
climate <- sim1$Climate %>% mutate(Cluster = factor(sim1$Cluster))
morph_long <- morph %>% pivot_longer(-Cluster, names_to = "Trait",  values_to = "Value")
climate_long <- climate %>% pivot_longer(-Cluster, names_to = "Var", values_to = "Value")

ggplot(morph_long, aes(x = Value, color = Cluster, fill = Cluster)) +
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



################################################################################
#### Main simulation parameters
################################################################################

overwrite <- TRUE
N_simulations <- 3
N_steps_SOM <- 5
N_replicates_SOM <- 5
max_NA_row_SOM <- 1
max_NA_col_SOM <- 1
verbose_SOM <- FALSE
BIC_threshold_SOM <- 6
learning_rate_tuning <- FALSE
pca_codebooks_SOM <- TRUE
if (!dir.exists("Simulations")) dir.create("Simulations")



################################################################################
#### Test effect of clustering method using simulated data
################################################################################

## Set clustering methods to test
clustering_methods <- c("GMM+BIC",
                        "GMM+BICthresh",
                        "GMM+gapstat",
                        "HDBSCAN",
                        "hierarchical+gapstat",
                        "kmeans+BICjump",
                        "kmeans+BICthresh",
                        "kmeans+BICjump+hierarch",
                        "kmeans+BICthresh+hierarch",
                        "kmeans+BICjump+spectral",
                        "kmeans+BICthresh+spectral",
                        "kmeans+gapstat",
                        "PAM+gapstat")


## Set parameters
N_clusters <- 3
training_neighborhoods_SOM <- "gaussian"
max_k_SOM <- 10


## File paths for saving/loading
sim_data_clustering_methods_file <- "Simulations/Sim_data_clustering_methods.rds"
sim_results_clustering_methods_csv <- "Simulations/Sim_results_clustering_methods.csv"


## Load saved results or simulate/run as needed
if (file.exists(sim_results_clustering_methods_csv) && !overwrite) {
  full_sim_stats_clustering_methods <- read.csv(sim_results_clustering_methods_csv)
  message("Loaded clustering methods simulation results from CSV, skipping run")
} else {
  
  if (file.exists(sim_data_clustering_methods_file) && !overwrite) {
    sim_results <- readRDS(sim_data_clustering_methods_file)
    message("Loaded clustering methods simulation data from RDS")
  } else {
    
    # Prepare seeds and simulate data
    set.seed(1)
    seeds <- sample(1:1e8, N_simulations)
    sim_results <- vector("list", N_simulations)
    for(i in seq_len(N_simulations)) {
      set.seed(seeds[i])
      sim_results[[i]] <- simulate.data()
    }
    saveRDS(sim_results, sim_data_clustering_methods_file)
    message("Simulations completed and saved to file")
  }
  
  all_results_clustering_methods <- list()
  
  for(cl_method in clustering_methods) {
    cat("Running for clustering method:", cl_method, "\n")
    
    som_results <- lapply(sim_results, function(sim_data) {
      som_input <- list(SNP = sim_data$SNP,
                        Morphology = sim_data$Morphology,
                        Climate = sim_data$Climate)
      # Time entire SOM+clustering step
      elapsed_time <- system.time({
        som_out <- tryCatch(
          train.SOM(input_data = som_input,
                    learning.rate.tuning = learning_rate_tuning,
                    training.neighborhoods = training_neighborhoods_SOM,
                    N.steps = N_steps_SOM,
                    N.replicates = N_replicates_SOM,
                    max.NA.row = max_NA_row_SOM,
                    max.NA.col = max_NA_col_SOM,
                    verbose = verbose_SOM),
          error = function(e) {
            message("SOM training failed: ", conditionMessage(e))
            return(NULL)
          }
        )
        clust <- if (!is.null(som_out)) {
          tryCatch(
            clustering.SOM(SOM.output = som_out,
                           max.k = max_k_SOM,
                           BIC.thresh = BIC_threshold_SOM,
                           pca.codebooks = pca_codebooks_SOM,
                           clustering.method = cl_method),
            error = function(e) {
              message("Clustering failed for method ", cl_method, 
                      ": ", conditionMessage(e))
              return(NULL)
            }
          )
        } else {
          NULL
        }
        list(som_models = if (!is.null(som_out)) som_out$som_models else NULL,
             clustering = clust)
      })["elapsed"]
      attr(elapsed_time, "value") # attach full value (SOM+clustering output) to elapsed_time for next steps
      attr(elapsed_time, "elapsed") <- elapsed_time
      elapsed_time
    })
    
    stats_list <- list()
    
    for(i_sim in seq_along(som_results)) {
      total_time <- som_results[[i_sim]]
      som_models <- attr(total_time, "value")$som_models
      clust <- attr(total_time, "value")$clustering
      sim_time <- attr(total_time, "elapsed")
      true_labels <- as.integer(as.factor(sim_results[[i_sim]]$Cluster))
      
      # If SOM or clustering failed, create NA rows for all replicates
      if (is.null(som_models) || is.null(clust)) {
        stats_df <- data.frame(
          clustering_method = cl_method,
          Time = sim_time,
          K_inferred = NA,
          ARI = NA,
          Acc = NA,
          K_correct = NA,
          K_far_off = NA,
          QE = NA,
          TE = NA,
          stringsAsFactors = FALSE
        )
        sim_stats_row <- sim_results[[i_sim]]$sim_stats
        for (col in colnames(sim_stats_row)) {
          stats_df[[col]] <- sim_stats_row[[col]]
        }
        stats_list[[i_sim]] <- stats_df[rep(1, N_replicates_SOM), , drop = FALSE]
        next
      }
      
      # Otherwise, calculate statistics for each replicate
      stats_per_rep <- lapply(seq_along(som_models), function(ridx) {
        som_model <- som_models[[ridx]]
        QE_TE <- extract.QE.TE(som_model)
        
        # Try to extract predicted labels; if fails, fill with NA
        res <- tryCatch({
          pred_labels <- as.integer(clust$cluster_assignment[, ridx])
          K_inferred <- length(unique(pred_labels))
          ARI <- adjustedRandIndex(true_labels, pred_labels)
          Acc <- if (K_inferred == N_clusters) get.accuracy(true_labels, pred_labels) else NA
          K_correct <- (K_inferred == N_clusters)
          K_far_off <- (abs(K_inferred - N_clusters) >= 2)
          data.frame(clustering_method = cl_method,
                     Time = sim_time,
                     K_inferred = K_inferred,
                     ARI = ARI,
                     Acc = Acc,
                     K_correct = K_correct,
                     K_far_off = K_far_off,
                     QE = QE_TE$QE,
                     TE = QE_TE$TE,
                     stringsAsFactors = FALSE)
        }, error = function(e) {
          data.frame(clustering_method = cl_method,
                     Time = sim_time,
                     K_inferred = NA,
                     ARI = NA,
                     Acc = NA,
                     K_correct = NA,
                     K_far_off = NA,
                     QE = NA,
                     TE = NA,
                     stringsAsFactors = FALSE)
        })
        res
      })
      
      stats_df <- do.call(rbind, stats_per_rep)
      sim_stats_row <- sim_results[[i_sim]]$sim_stats
      for (col in colnames(sim_stats_row)) {
        stats_df[[col]] <- sim_stats_row[[col]]
      }
      stats_list[[i_sim]] <- stats_df
    }
    
    stats_df_all <- do.call(rbind, stats_list)
    stats_df_all$sim <- rep(seq_along(stats_list), each = nrow(stats_list[[1]]))
    all_results_clustering_methods[[cl_method]] <- stats_df_all
  }
  
  
  ## Combine and save all results
  full_sim_stats_clustering_methods <- do.call(rbind, all_results_clustering_methods)
  write.csv(full_sim_stats_clustering_methods, sim_results_clustering_methods_csv, row.names = FALSE)
  message("Clustering methods simulation results saved to CSV")
}


## Summarize results
summary_sim_stats_clustering_methods <- full_sim_stats_clustering_methods %>%
  group_by(clustering_method) %>%
  summarise(mean_ARI = mean(ARI, na.rm = TRUE),
            mean_K = mean(K_inferred, na.rm = TRUE),
            mean_Acc = mean(Acc, na.rm = TRUE),
            mean_Time = mean(Time, na.rm = TRUE),
            mean_QE = mean(QE, na.rm = TRUE),
            mean_TE = mean(TE, na.rm = TRUE),
            n = n()) %>%
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
  labs(x = "Clustering method", y = "K far off (1 = ≥2 away from true K)")

ggplot(full_sim_stats_clustering_methods, aes(x = factor(clustering_method), y = QE)) +
  geom_boxplot(fill = "#66c2a5", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Clustering method", y = "Mean quantization error (QE)")

ggplot(full_sim_stats_clustering_methods, aes(x = factor(clustering_method), y = TE)) +
  geom_boxplot(fill = "#fc8d62", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Clustering method", y = "Mean topographic error (TE)")

ggplot(full_sim_stats_clustering_methods, aes(x = factor(clustering_method), y = Time)) +
  geom_boxplot(fill = "#fc9192", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "Clustering method", y = "Time (seconds)")



################################################################################
#### Test effect of codebook PCA using simulated data
################################################################################

## Set PCA options to test
pca_SOM_options <- c(TRUE, FALSE)


## Set parameters
clustering_method <- "kmeans+BICjump+hierarch"
N_clusters <- 3
training_neighborhoods_SOM <- "gaussian"
max_k_SOM <- 10


## File paths for saving/loading
sim_data_pca_file <- "Simulations/Sim_data_pca.rds"
sim_results_pca_csv <- "Simulations/Sim_results_pca.csv"


## Load saved results or simulate/run as needed
if (file.exists(sim_results_pca_csv) && !overwrite) {
  full_sim_stats_pca <- read.csv(sim_results_pca_csv)
} else {
  
  if (file.exists(sim_data_pca_file) && !overwrite) {
    sim_results <- readRDS(sim_data_pca_file)
  } else {
    
    set.seed(1)
    seeds <- sample(1:1e8, N_simulations)
    sim_results <- vector("list", N_simulations)
    for(i in seq_len(N_simulations)) {
      set.seed(seeds[i])
      sim_results[[i]] <- simulate.data()
    }
    saveRDS(sim_results, sim_data_pca_file)
  }
  
  all_results_pca <- list()
  
  for(pca_opt in pca_SOM_options) {
    
    som_results <- vector("list", length(sim_results))
    
    for(i in seq_along(sim_results)) {
      
      sim_data <- sim_results[[i]]
      
      som_input <- list(
        SNP = sim_data$SNP,
        Morphology = sim_data$Morphology,
        Climate = sim_data$Climate
      )
      
      som_results[[i]] <- tryCatch({
        result <- NULL
        pt <- system.time({
          som_out <- train.SOM(
            input_data = som_input,
            learning.rate.tuning = learning_rate_tuning,
            training.neighborhoods = training_neighborhoods_SOM,
            N.steps = N_steps_SOM,
            N.replicates = N_replicates_SOM,
            max.NA.row = max_NA_row_SOM,
            max.NA.col = max_NA_col_SOM,
            verbose = FALSE
          )
          clust <- clustering.SOM(
            SOM.output = som_out,
            max.k = max_k_SOM,
            BIC.thresh = BIC_threshold_SOM,
            pca.codebooks = pca_opt,
            clustering.method = clustering_method
          )
          result <<- list(som_models = som_out$som_models, clustering = clust)
        })
        result$elapsed <- pt[["elapsed"]]
        result
      }, error = function(e) {
        NULL
      })
    }
    
    stats_list <- list()
    
    for(i_sim in seq_along(som_results)) {
      sr <- som_results[[i_sim]]
      if (is.null(sr)) next
      
      som_models <- sr$som_models
      clust <- sr$clustering
      total_time <- sr$elapsed
      true_labels <- as.integer(as.factor(sim_results[[i_sim]]$Cluster))
      
      if (is.null(som_models) || is.null(clust) || is.null(clust$cluster_assignment)) next
      
      stats_per_rep <- lapply(som_models, function(som_model) {
        QE_TE <- extract.QE.TE(som_model)
        pred_labels <- as.integer(clust$cluster_assignment[, 1])
        K_inferred <- length(unique(pred_labels))
        ARI <- adjustedRandIndex(true_labels, pred_labels)
        Acc <- if (K_inferred == N_clusters) get.accuracy(true_labels, pred_labels) else NA
        K_correct <- (K_inferred == N_clusters)
        K_far_off <- (abs(K_inferred - N_clusters) >= 2)
        
        data.frame(pca_codebooks = pca_opt,
                   Time = total_time,
                   K_inferred = K_inferred,
                   ARI = ARI,
                   Acc = Acc,
                   K_correct = K_correct,
                   K_far_off = K_far_off,
                   QE = QE_TE$QE,
                   TE = QE_TE$TE,
                   stringsAsFactors = FALSE)
      })
      
      stats_df <- do.call(rbind, stats_per_rep)
      sim_stats_row <- sim_results[[i_sim]]$sim_stats
      for (col in colnames(sim_stats_row)) {
        stats_df[[col]] <- sim_stats_row[[col]]
      }
      stats_list[[i_sim]] <- stats_df
    }
    
    non_empty <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, stats_list)
    if (length(non_empty) > 0) {
      stats_df_all <- do.call(rbind, non_empty)
      stats_df_all$sim <- rep(seq_along(non_empty), each = nrow(non_empty[[1]]))
      all_results_pca[[as.character(pca_opt)]] <- stats_df_all
    }
  }
  
  all_results_pca <- all_results_pca[sapply(all_results_pca, is.data.frame)]
  if (length(all_results_pca) == 0) stop("No results to summarize")
  full_sim_stats_pca <- do.call(rbind, all_results_pca)
  full_sim_stats_pca <- as.data.frame(full_sim_stats_pca, stringsAsFactors = FALSE)
  write.csv(full_sim_stats_pca, sim_results_pca_csv, row.names = FALSE)
}


## Summarize results
summary_sim_stats_pca <- full_sim_stats_pca %>%
  dplyr::group_by(pca_codebooks) %>%
  dplyr::summarise(mean_ARI = mean(ARI, na.rm = TRUE),
                   mean_K = mean(K_inferred, na.rm = TRUE),
                   mean_Acc = mean(Acc, na.rm = TRUE),
                   mean_Time = mean(Time, na.rm = TRUE),
                   mean_QE = mean(QE, na.rm = TRUE),
                   mean_TE = mean(TE, na.rm = TRUE),
                   n = dplyr::n()) %>%
  as.data.frame()
print(summary_sim_stats_pca)


## Plot results
ggplot(full_sim_stats_pca, aes(x = factor(pca_codebooks), y = ARI)) +
  geom_boxplot(fill = "#abd9e9", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "PCA codebooks", y = "Adjusted Rand Index (ARI)")

ggplot(full_sim_stats_pca, aes(x = factor(pca_codebooks), y = Acc)) +
  geom_boxplot(fill = "#74add1", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "PCA codebooks", y = "Assignment accuracy")

ggplot(full_sim_stats_pca, aes(x = factor(pca_codebooks), y = as.numeric(K_correct))) +
  geom_boxplot(fill = "#fdae61", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "PCA codebooks", y = "K correct (1 = correct)")

ggplot(full_sim_stats_pca, aes(x = factor(pca_codebooks), y = as.numeric(K_far_off))) +
  geom_boxplot(fill = "#d73027", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "PCA codebooks", y = "K far off (1 = ≥2 away from true K)")

ggplot(full_sim_stats_pca, aes(x = factor(pca_codebooks), y = QE)) +
  geom_boxplot(fill = "#66c2a5", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "PCA codebooks", y = "Mean quantization error (QE)")

ggplot(full_sim_stats_pca, aes(x = factor(pca_codebooks), y = TE)) +
  geom_boxplot(fill = "#fc8d62", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "PCA codebooks", y = "Mean topographic error (TE)")

ggplot(full_sim_stats_pca, aes(x = factor(pca_codebooks), y = Time)) +
  geom_boxplot(fill = "#fc9192", outlier.size = 0.7) +
  theme_classic() +
  labs(x = "PCA codebooks", y = "Time (seconds)")



################################################################################
#### Test effect of missing data proportion using simulated data
################################################################################


## Set missing data proportions (NA proportions) to test
missing_data_props <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)


## Set parameters
clustering_method <- "kmeans+BICjump+hierarch"
N_clusters <- 3
training_neighborhoods_SOM <- "gaussian"
max_k_SOM <- 10
pca_codebooks_SOM <- TRUE


## File paths for saving/loading
sim_data_NA_file <- "Simulations/Sim_data_NA.rds"
sim_results_NA_csv <- "Simulations/Sim_results_NA.csv"


## Load saved results or simulate/run as needed
if (file.exists(sim_results_NA_csv) && !overwrite) {
  full_sim_stats_NA <- read.csv(sim_results_NA_csv)
} else {
  
  if (file.exists(sim_data_NA_file) && !overwrite) {
    sim_results_all <- readRDS(sim_data_NA_file)
  } else {
    
    set.seed(1)
    seeds <- sample(1:1e8, N_simulations)
    
    sim_results_all <- list()
    for(na_prop in missing_data_props) {
      sim_results <- vector("list", N_simulations)
      for(i in seq_len(N_simulations)) {
        set.seed(seeds[i])
        sim_results[[i]] <- simulate.data(missing_data_prop = as.numeric(na_prop))
      }
      sim_results_all[[as.character(na_prop)]] <- sim_results
    }
    saveRDS(sim_results_all, sim_data_NA_file)
  }
  
  all_results_NA <- list()
  
  for(na_prop in names(sim_results_all)) {
    
    sim_results <- sim_results_all[[na_prop]]
    
    som_results <- vector("list", length(sim_results))
    
    for(i in seq_along(sim_results)) {
      
      sim_data <- sim_results[[i]]
      
      som_input <- list(
        SNP = sim_data$SNP,
        Morphology = sim_data$Morphology,
        Climate = sim_data$Climate
      )
      
      som_results[[i]] <- tryCatch({
        result <- NULL
        pt <- system.time({
          som_out <- train.SOM(
            input_data = som_input,
            learning.rate.tuning = learning_rate_tuning,
            training.neighborhoods = training_neighborhoods_SOM,
            N.steps = N_steps_SOM,
            N.replicates = N_replicates_SOM,
            max.NA.row = max_NA_row_SOM,
            max.NA.col = max_NA_col_SOM,
            verbose = FALSE
          )
          clust <- clustering.SOM(
            SOM.output = som_out,
            max.k = max_k_SOM,
            BIC.thresh = BIC_threshold_SOM,
            pca.codebooks = pca_codebooks_SOM,
            clustering.method = clustering_method
          )
          result <<- list(som_models = som_out$som_models, clustering = clust)
        })
        result$elapsed <- pt[["elapsed"]]
        result
      }, error = function(e) {
        NULL
      })
    }
    
    stats_list <- list()
    
    for(i_sim in seq_along(som_results)) {
      sr <- som_results[[i_sim]]
      if (is.null(sr)) next
      
      som_models <- sr$som_models
      clust <- sr$clustering
      total_time <- sr$elapsed
      true_labels <- as.integer(as.factor(sim_results[[i_sim]]$Cluster))
      
      if (is.null(som_models) || is.null(clust) || is.null(clust$cluster_assignment)) next
      
      stats_per_rep <- lapply(som_models, function(som_model) {
        QE_TE <- extract.QE.TE(som_model)
        pred_labels <- as.integer(clust$cluster_assignment[, 1])
        K_inferred <- length(unique(pred_labels))
        ARI <- adjustedRandIndex(true_labels, pred_labels)
        Acc <- if (K_inferred == N_clusters) get.accuracy(true_labels, pred_labels) else NA
        K_correct <- (K_inferred == N_clusters)
        K_far_off <- (abs(K_inferred - N_clusters) >= 2)
        
        data.frame(missing_data_prop = as.numeric(na_prop),
                   Time = total_time,
                   K_inferred = K_inferred,
                   ARI = ARI,
                   Acc = Acc,
                   K_correct = K_correct,
                   K_far_off = K_far_off,
                   QE = QE_TE$QE,
                   TE = QE_TE$TE,
                   stringsAsFactors = FALSE)
      })
      
      stats_df <- do.call(rbind, stats_per_rep)
      sim_stats_row <- sim_results[[i_sim]]$sim_stats
      for (col in colnames(sim_stats_row)) {
        stats_df[[col]] <- sim_stats_row[[col]]
      }
      stats_list[[i_sim]] <- stats_df
    }
    
    non_empty <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, stats_list)
    if (length(non_empty) > 0) {
      stats_df_all <- do.call(rbind, non_empty)
      stats_df_all$sim <- rep(seq_along(non_empty), each = nrow(non_empty[[1]]))
      all_results_NA[[as.character(na_prop)]] <- stats_df_all
    }
  }
  
  all_results_NA <- all_results_NA[sapply(all_results_NA, is.data.frame)]
  if (length(all_results_NA) == 0) stop("No results to summarize")
  full_sim_stats_NA <- do.call(rbind, all_results_NA)
  full_sim_stats_NA <- as.data.frame(full_sim_stats_NA, stringsAsFactors = FALSE)
  write.csv(full_sim_stats_NA, sim_results_NA_csv, row.names = FALSE)
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
  labs(x = "Missing data proportion", y = "K far off (1 = ≥2 away from true K)")

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


## Set parameters
clustering_method <- "kmeans+BICjump+hierarch"
N_clusters <- 3
max_k_SOM <- 10
pca_codebooks_SOM <- TRUE


## File paths for saving/loading
sim_data_neighborhoods_file <- "Simulations/Sim_data_neighborhoods.rds"
sim_results_neighborhoods_csv <- "Simulations/Sim_results_neighborhoods.csv"


## Load saved results or simulate/run as needed
if (file.exists(sim_results_neighborhoods_csv) && !overwrite) {
  full_sim_stats_neighborhood <- read.csv(sim_results_neighborhoods_csv)
} else {
  
  if (file.exists(sim_data_neighborhoods_file) && !overwrite) {
    sim_results <- readRDS(sim_data_neighborhoods_file)
  } else {
    
    set.seed(1)
    seeds <- sample(1:1e8, N_simulations)
    sim_results <- vector("list", N_simulations)
    for(i in seq_len(N_simulations)) {
      set.seed(seeds[i])
      sim_results[[i]] <- simulate.data()
    }
    saveRDS(sim_results, sim_data_neighborhoods_file)
  }
  
  all_results_neigh <- list()
  
  for(neighborhoods_fun in neighborhoods) {
    
    som_results <- vector("list", length(sim_results))
    
    for(i in seq_along(sim_results)) {
      
      sim_data <- sim_results[[i]]
      
      som_input <- list(
        SNP = sim_data$SNP,
        Morphology = sim_data$Morphology,
        Climate = sim_data$Climate
      )
      
      som_results[[i]] <- tryCatch({
        result <- NULL
        pt <- system.time({
          som_out <- train.SOM(
            input_data = som_input,
            learning.rate.tuning = learning_rate_tuning,
            training.neighborhoods = neighborhoods_fun,
            N.steps = N_steps_SOM,
            N.replicates = N_replicates_SOM,
            max.NA.row = max_NA_row_SOM,
            max.NA.col = max_NA_col_SOM,
            verbose = FALSE
          )
          clust <- clustering.SOM(
            SOM.output = som_out,
            max.k = max_k_SOM,
            BIC.thresh = BIC_threshold_SOM,
            pca.codebooks = pca_codebooks_SOM,
            clustering.method = clustering_method
          )
          result <<- list(som_models = som_out$som_models, clustering = clust)
        })
        result$elapsed <- pt[["elapsed"]]
        result
      }, error = function(e) {
        NULL
      })
    }
    
    stats_list <- list()
    
    for(i_sim in seq_along(som_results)) {
      sr <- som_results[[i_sim]]
      if (is.null(sr)) next
      
      som_models <- sr$som_models
      clust <- sr$clustering
      total_time <- sr$elapsed
      true_labels <- as.integer(as.factor(sim_results[[i_sim]]$Cluster))
      
      if (is.null(som_models) || is.null(clust) || is.null(clust$cluster_assignment)) next
      
      stats_per_rep <- lapply(som_models, function(som_model) {
        QE_TE <- extract.QE.TE(som_model)
        pred_labels <- as.integer(clust$cluster_assignment[, 1])
        K_inferred <- length(unique(pred_labels))
        ARI <- adjustedRandIndex(true_labels, pred_labels)
        Acc <- if (K_inferred == N_clusters) get.accuracy(true_labels, pred_labels) else NA
        K_correct <- (K_inferred == N_clusters)
        K_far_off <- (abs(K_inferred - N_clusters) >= 2)
        
        data.frame(neighborhood_function = neighborhoods_fun,
                   Time = total_time,
                   K_inferred = K_inferred,
                   ARI = ARI,
                   Acc = Acc,
                   K_correct = K_correct,
                   K_far_off = K_far_off,
                   QE = QE_TE$QE,
                   TE = QE_TE$TE,
                   stringsAsFactors = FALSE)
      })
      
      stats_df <- do.call(rbind, stats_per_rep)
      sim_stats_row <- sim_results[[i_sim]]$sim_stats
      for (col in colnames(sim_stats_row)) {
        stats_df[[col]] <- sim_stats_row[[col]]
      }
      stats_list[[i_sim]] <- stats_df
    }
    
    non_empty <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, stats_list)
    if (length(non_empty) > 0) {
      stats_df_all <- do.call(rbind, non_empty)
      stats_df_all$sim <- rep(seq_along(non_empty), each = nrow(non_empty[[1]]))
      all_results_neigh[[as.character(neighborhoods_fun)]] <- stats_df_all
    }
  }
  
  all_results_neigh <- all_results_neigh[sapply(all_results_neigh, is.data.frame)]
  if (length(all_results_neigh) == 0) stop("No results to summarize")
  full_sim_stats_neighborhood <- do.call(rbind, all_results_neigh)
  full_sim_stats_neighborhood <- as.data.frame(full_sim_stats_neighborhood, stringsAsFactors = FALSE)
  write.csv(full_sim_stats_neighborhood, sim_results_neighborhoods_csv, row.names = FALSE)
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
  labs(x = "Neighborhood function", y = "K far off (1 = ≥2 away from true K)")

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


## Set parameters
clustering_method <- "kmeans+BICjump+hierarch"
N_clusters <- 3
training_neighborhoods_SOM <- "gaussian"
max_k_SOM <- 10
pca_codebooks_SOM <- TRUE


## File paths for saving/loading
sim_data_learning_rate_tuning_file <- "Simulations/Sim_data_learning_rate_tuning.rds"
sim_results_learning_rate_tuning_csv <- "Simulations/Sim_results_learning_rate_tuning.csv"


## Load saved results or simulate/run as needed
if (file.exists(sim_results_learning_rate_tuning_csv) && !overwrite) {
  full_sim_stats_learning_rate_tuning <- read.csv(sim_results_learning_rate_tuning_csv)
} else {
  
  if (file.exists(sim_data_learning_rate_tuning_file) && !overwrite) {
    sim_results <- readRDS(sim_data_learning_rate_tuning_file)
  } else {
    set.seed(1)
    seeds <- sample(1:1e8, N_simulations)
    sim_results <- vector("list", N_simulations)
    for(i in seq_len(N_simulations)) {
      set.seed(seeds[i])
      sim_results[[i]] <- simulate.data()
    }
    saveRDS(sim_results, sim_data_learning_rate_tuning_file)
  }
  
  all_results_learning_rate_tuning <- list()
  
  for(learning_rate_tuning_opt in learning_rate_tuning_SOM) {
    
    som_results <- vector("list", length(sim_results))
    
    for(i in seq_along(sim_results)) {
      
      sim_data <- sim_results[[i]]
      
      som_input <- list(
        SNP = sim_data$SNP,
        Morphology = sim_data$Morphology,
        Climate = sim_data$Climate
      )
      
      som_results[[i]] <- tryCatch({
        result <- NULL
        pt <- system.time({
          som_out <- train.SOM(
            input_data = som_input,
            learning.rate.tuning = learning_rate_tuning_opt,
            training.neighborhoods = training_neighborhoods_SOM,
            N.steps = N_steps_SOM,
            N.replicates = N_replicates_SOM,
            max.NA.row = max_NA_row_SOM,
            max.NA.col = max_NA_col_SOM,
            verbose = FALSE
          )
          clust <- clustering.SOM(
            SOM.output = som_out,
            max.k = max_k_SOM,
            BIC.thresh = BIC_threshold_SOM,
            pca.codebooks = pca_codebooks_SOM,
            clustering.method = clustering_method
          )
          result <<- list(som_models = som_out$som_models, clustering = clust)
        })
        result$elapsed <- pt[["elapsed"]]
        result
      }, error = function(e) {
        NULL
      })
    }
    
    stats_list <- list()
    
    for(i_sim in seq_along(som_results)) {
      sr <- som_results[[i_sim]]
      if (is.null(sr)) next
      
      som_models <- sr$som_models
      clust <- sr$clustering
      total_time <- sr$elapsed
      true_labels <- as.integer(as.factor(sim_results[[i_sim]]$Cluster))
      
      if (is.null(som_models) || is.null(clust) || is.null(clust$cluster_assignment)) next
      
      stats_per_rep <- lapply(som_models, function(som_model) {
        QE_TE <- extract.QE.TE(som_model)
        pred_labels <- as.integer(clust$cluster_assignment[, 1])
        K_inferred <- length(unique(pred_labels))
        ARI <- adjustedRandIndex(true_labels, pred_labels)
        Acc <- if (K_inferred == N_clusters) get.accuracy(true_labels, pred_labels) else NA
        K_correct <- (K_inferred == N_clusters)
        K_far_off <- (abs(K_inferred - N_clusters) >= 2)
        
        data.frame(learning_rate_tuning_codebooks = learning_rate_tuning_opt,
                   Time = total_time,
                   K_inferred = K_inferred,
                   ARI = ARI,
                   Acc = Acc,
                   K_correct = K_correct,
                   K_far_off = K_far_off,
                   QE = QE_TE$QE,
                   TE = QE_TE$TE,
                   stringsAsFactors = FALSE)
      })
      
      stats_df <- do.call(rbind, stats_per_rep)
      sim_stats_row <- sim_results[[i_sim]]$sim_stats
      for (col in colnames(sim_stats_row)) {
        stats_df[[col]] <- sim_stats_row[[col]]
      }
      stats_list[[i_sim]] <- stats_df
    }
    
    non_empty <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, stats_list)
    if (length(non_empty) > 0) {
      stats_df_all <- do.call(rbind, non_empty)
      stats_df_all$sim <- rep(seq_along(non_empty), each = nrow(non_empty[[1]]))
      all_results_learning_rate_tuning[[as.character(learning_rate_tuning_opt)]] <- stats_df_all
    }
  }
  
  all_results_learning_rate_tuning <- all_results_learning_rate_tuning[sapply(all_results_learning_rate_tuning, is.data.frame)]
  if (length(all_results_learning_rate_tuning) == 0) stop("No results to summarize")
  full_sim_stats_learning_rate_tuning <- do.call(rbind, all_results_learning_rate_tuning)
  full_sim_stats_learning_rate_tuning <- as.data.frame(full_sim_stats_learning_rate_tuning, stringsAsFactors = FALSE)
  write.csv(full_sim_stats_learning_rate_tuning, sim_results_learning_rate_tuning_csv, row.names = FALSE)
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
  labs(x = "learning rate tuning of codebooks", y = "K far off (1 = ≥2 away from true K)")

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
correlation_cutoffs <- c(0.6, 0.7, 0.8, 0.9, 1.0)


## Set parameters
N_traits <- 100
clustering_method <- "kmeans+BICjump+hierarch"
correlation_mean <- 0.5
correlation_SD <- 0.25


## File paths
sim_data_corr_file <- "Simulations/Sim_data_corr.rds"
sim_results_corr_csv <- "Simulations/Sim_results_corr.csv"


## Create function to generate positive-definite correlation matrix
generate.correlation.matrix <- function(N.traits = N_traits,
                                        mean.corr = correlation_mean,
                                        SD.corr = correlation_SD) {
  n_off_diag <- N.traits * (N.traits - 1) / 2
  off_diag <- rnorm(n_off_diag, mean = mean.corr, sd = SD.corr)
  off_diag <- pmin(pmax(off_diag, 0), 1)
  corr_mat <- diag(1, N.traits)
  corr_mat[lower.tri(corr_mat)] <- off_diag
  corr_mat[upper.tri(corr_mat)] <- t(corr_mat)[upper.tri(corr_mat)]
  eigvals <- eigen(corr_mat, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigvals <= 0)) {
    corr_mat <- corr_mat + diag(abs(min(eigvals)) + 1e-6, N.traits)
  }
  return(corr_mat)
}


## Create function to simulate two datasets with clusters and correlation structure
simulate.two.simple.datasets <- function(N.individuals = N_individuals,
                                         N.clusters = N_clusters,
                                         N.traits = N_traits,
                                         mean.corr = correlation_mean,
                                         SD.corr = correlation_SD) {
  cluster.labels <- rep(seq_len(N.clusters), length.out = N.individuals)
  cov.mat1 <- generate.correlation.matrix(N.traits, mean.corr, SD.corr)
  cov.mat2 <- generate.correlation.matrix(N.traits, mean.corr, SD.corr)
  cluster.means1 <- matrix(runif(N.clusters * N.traits, 10, 30), nrow = N.clusters)
  cluster.means2 <- matrix(runif(N.clusters * N.traits, 10, 30), nrow = N.clusters)
  data1 <- matrix(NA, nrow = N.individuals, ncol = N.traits)
  data2 <- matrix(NA, nrow = N.individuals, ncol = N.traits)
  for (clust in seq_len(N.clusters)) {
    inds <- which(cluster.labels == clust)
    data1[inds, ] <- MASS::mvrnorm(n = length(inds), mu = cluster.means1[clust, ], Sigma = cov.mat1)
    data2[inds, ] <- MASS::mvrnorm(n = length(inds), mu = cluster.means2[clust, ], Sigma = cov.mat2)
  }
  list(Dataset1 = data1, Dataset2 = data2, Cluster = cluster.labels)
}


## Load saved results or simulate/run as needed
if (file.exists(sim_results_corr_csv) && !overwrite) {
  full_sim_stats_corr <- read.csv(sim_results_corr_csv)
} else {
  
  if (file.exists(sim_data_corr_file) && !overwrite) {
    sim_results <- readRDS(sim_data_corr_file)
  } else {
    set.seed(1)
    seeds <- sample(1:1e8, N_simulations)
    sim_results <- vector("list", N_simulations)
    for (i in seq_len(N_simulations)) {
      set.seed(seeds[i])
      sim_results[[i]] <- simulate.two.simple.datasets()
    }
    saveRDS(sim_results, sim_data_corr_file)
  }
  
  all_results_corr <- list()
  
  for (cutoff in correlation_cutoffs) {
    
    sim_results_filtered <- lapply(sim_results, function(sim_data) {
      list(
        Dataset1 = remove.lowCV.multicollinearity.SOM(sim_data$Dataset1,
                                                      cor.threshold = cutoff,
                                                      verbose = FALSE),
        Dataset2 = remove.lowCV.multicollinearity.SOM(sim_data$Dataset2,
                                                      cor.threshold = cutoff,
                                                      verbose = FALSE),
        Cluster = sim_data$Cluster
      )
    })
    
    stats_list <- vector("list", length(sim_results_filtered))
    
    for (i_sim in seq_along(sim_results_filtered)) {
      sim_data <- sim_results_filtered[[i_sim]]
      
      result <- tryCatch({
        res <- NULL
        pt <- system.time({
          som_out <- train.SOM(input_data = list(Dataset1 = sim_data$Dataset1,
                                                 Dataset2 = sim_data$Dataset2),
                               learning.rate.tuning = learning_rate_tuning,
                               training.neighborhoods = training_neighborhoods_SOM,
                               N.steps = N_steps_SOM,
                               N.replicates = N_replicates_SOM,
                               max.NA.row = max_NA_row_SOM,
                               max.NA.col = max_NA_col_SOM,
                               verbose = FALSE)
          clust <- clustering.SOM(SOM.output = som_out,
                                  max.k = max_k_SOM,
                                  BIC.thresh = BIC_threshold_SOM,
                                  pca.codebooks = pca_codebooks_SOM,
                                  clustering.method = clustering_method)
          res <<- list(som_models = som_out$som_models, clustering = clust)
        })
        list(res = res, elapsed = pt["elapsed"])
      }, error = function(e) {
        NULL
      })
      
      if (is.null(result)) next
      
      som_models <- result$res$som_models
      clust <- result$res$clustering
      total_time <- result$elapsed
      true_labels <- as.integer(factor(sim_data$Cluster))
      
      if (is.null(som_models) || is.null(clust) || is.null(clust$cluster_assignment)) next
      
      stats_per_rep <- lapply(som_models, function(som_model) {
        QE_TE <- extract.QE.TE(som_model)
        pred_labels <- as.integer(clust$cluster_assignment[, 1])
        K_inferred <- length(unique(pred_labels))
        ARI <- adjustedRandIndex(true_labels, pred_labels)
        Acc <- if (K_inferred == N_clusters) get.accuracy(true_labels, pred_labels) else NA
        K_correct <- (K_inferred == N_clusters)
        K_far_off <- (abs(K_inferred - N_clusters) >= 2)
        data.frame(correlation_cutoff = cutoff,
                   Time = total_time,
                   K_inferred = K_inferred,
                   ARI = ARI,
                   Acc = Acc,
                   K_correct = K_correct,
                   K_far_off = K_far_off,
                   QE = QE_TE$QE,
                   TE = QE_TE$TE,
                   stringsAsFactors = FALSE)
      })
      
      stats_list[[i_sim]] <- do.call(rbind, stats_per_rep)
    }
    
    non_empty <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, stats_list)
    if (length(non_empty) > 0) {
      stats_df_all <- do.call(rbind, non_empty)
      stats_df_all$sim <- rep(seq_along(non_empty), each = nrow(non_empty[[1]]))
      all_results_corr[[as.character(cutoff)]] <- stats_df_all
    }
  }
  
  all_results_corr <- all_results_corr[sapply(all_results_corr, is.data.frame)]
  full_sim_stats_corr <- do.call(rbind, all_results_corr)
  write.csv(full_sim_stats_corr, sim_results_corr_csv, row.names = FALSE)
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
  labs(x = "Correlation cutoff", y = "K far off (1 = ≥2 away from true K)")

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
