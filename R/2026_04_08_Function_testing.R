################################################################################
#### Set environment
################################################################################

rm(list = ls()) #clear environment
#setwd("C:/Users/danie/Desktop/PhD research/SOM package")
setwd("./")
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")



###############################################################################
#### Test SOM functions and arguments with simple simulated datasets
################################################################################

## Specify parameters for all datasets
set.seed(1) #set seed for reproducibility
n_individuals <- 90 #simulate dataset with N individuals (needs to be consistent across all datasets)
rownames_datasets <- "Individual" #name rows (needs to be consistent across all datasets)


## Allele data (Alleles)
n_Alleles <- 120 #simulate dataset with 120 alleles
allele_frequencies <- c('1' = 0.4, '0.5' = 0.25, '0' = 0.35) #set allele frequencies
Alleles <- data.frame(lapply(1:n_Alleles, function(x) sample(names(allele_frequencies), n_individuals, replace = T, prob = allele_frequencies)),
                      row.names = paste0(rownames_datasets, 1:n_individuals)) #generate Alleles based on these frequencies
Alleles[] <- lapply(Alleles, function(x) as.numeric(as.character(x))) #convert character values to numeric (0, 0.5, 1)
colnames(Alleles) <- paste("Allele", 1:n_Alleles) #rename columns as "Allele 1", "Allele 2", ..., "Allele n"


## Environmental data (ENV)
n_env <- 40 #simulate dataset with 40 environmental variables
ENV <- matrix(runif(n_individuals * n_env, min = 0, max = 100), nrow = n_individuals, ncol = n_env)
rownames(ENV) <- paste0(rownames_datasets, 1:n_individuals)
na_indices_env <- sample(1:(n_individuals * n_env), size = round(n_individuals * n_env * 0.3), replace = F) #introduce some rare NAs for realism
colnames(ENV) <- paste("BIO", 1:n_env) #rename columns as "BIO 1", "BIO 2", ..., "BIO n"
ENV[na_indices_env] <- NA
ENV <- as.data.frame(ENV)
ENV$non_numeric_col <- sample(c("A", "B", "C"), n_individuals, replace = TRUE) #add non-numeric columns
ENV$zero_var_col <- 0.5 #add single zero-variance column
ENV$all_NA_col <- NA #add single all-NA column


## Morphological data (MORPH)
n_morph <- 20 #simulate morphological data for 20 traits
MORPH <- matrix(rnorm(n_individuals * n_morph, mean = 5, sd = 2), nrow = n_individuals, ncol = n_morph) #simulate morphological data
rownames(MORPH) <- paste0(rownames_datasets, 1:n_individuals) #name rows
na_indices_morph <- sample(1:(n_individuals * n_morph), size = round(n_morph * n_individuals * 0.4), replace = F) #introduce some NAs
MORPH[na_indices_morph] <- NA
MORPH <- cbind(MORPH, Trait_constant = 1)
colnames(MORPH) <- c(paste("Trait", 1:n_morph), "Trait_constant")


## Simulate k3 data for each cluster (simulating k = 3)
n_clusters <- 3 #number of clusters
n_k3_test <- 50  #number of traits
clusters <- sample(1:n_clusters, n_individuals, replace = T) #assign each individual to cluster
k3_test <- matrix(NA, nrow = n_individuals, ncol = n_k3_test)
colnames(k3_test) <- paste("Trait", 1:n_k3_test)
rownames(k3_test) <- paste0("Individual", 1:n_individuals)

k3_test_means <- list( #define cluster-specific means for k3 traits
  cluster_1 = rnorm(n_k3_test, mean = 10, sd = 1.6),
  cluster_2 = rnorm(n_k3_test, mean = 15, sd = 1.6),
  cluster_3 = rnorm(n_k3_test, mean = 20, sd = 1.6))

for (i in 1:n_clusters) { #assign data for each cluster with some random variation
  cluster_indices <- which(clusters == i)
  k3_test[cluster_indices, ] <- matrix(rnorm(length(cluster_indices) * n_k3_test,
                                             mean = k3_test_means[[paste0("cluster_", i)]],
                                             sd = 3),
                                       nrow = length(cluster_indices),
                                       ncol = n_k3_test)}
na_indices_k3_test <- sample(1:(n_individuals * n_k3_test), #introduce some missing values
                             size = round(n_k3_test * n_individuals * 0.2),
                             replace = F)
k3_test[na_indices_k3_test] <- NA
k3_test <- as.data.frame(k3_test)
par(mfrow = c(3, 1))
hist(k3_test$`Trait 1`, breaks = 30, main = "Histogram of Trait 1") #plot histograms for trait 1
hist(k3_test$`Trait 2`, breaks = 30, main = "Histogram of Trait 2") #plot histograms for trait 2
hist(k3_test$`Trait 3`, breaks = 30, main = "Histogram of Trait 3") #plot histograms for trait 3
par(mfrow = c(1, 1))

## Test_D dataset
n <- 100
mat1 <- data.frame(
  A = sample(c(1:5, NA), n, TRUE),
  B = sample(c(1:5, NA), n, TRUE),
  C = sample(c(1:5, NA), n, TRUE),
  D = sample(c(1:5, NA), n, TRUE),
  row.names = paste0("S", 1:n))
mat2 <- data.frame(
  E = sample(c(1:5, NA), n, TRUE),
  F = sample(c(1:5, NA), n, TRUE),
  G = sample(c(1:5, NA), n, TRUE),
  row.names = paste0("S", 1:n))
TestD <- list(Layer1 = mat1, Layer2 = mat2)



## Simulate coordinate data
GER_Coordinates <- data.frame(Latitude = runif(n = n_individuals, min = 47.0, max = 53.0),
                              Longitude = runif(n = n_individuals, min = 7.5, max = 14.5))
rownames(GER_Coordinates) <- rownames(MORPH)

US_Coordinates <- data.frame(Latitude = runif(n = n_individuals, min = 29, max = 41),
                             Longitude = runif(n = n_individuals, min = -90, max = -80))
rownames(US_Coordinates) <- rownames(MORPH)


## Test train SOM function
SOM_single_Alleles <- train.SOM(Alleles, max.NA.row = 0.3, N.steps = 20, N.replicates = 3)
try(SOM_single_ENV <- train.SOM(ENV, max.NA.col = 0.05, max.NA.row = 0.4)) #will fail with error
try(SOM_single_ENV <- train.SOM(ENV, max.NA.col = 0.15, max.NA.row = 0.4)) #will fail with error
try(SOM_single_ENV <- train.SOM(ENV, max.NA.col = 0.3, max.NA.row = 0.2)) #will fail with error
SOM_single_ENV <- train.SOM(ENV, max.NA.col = 0.3, max.NA.row = 0.2, grid.multiplier = 3)
SOM_single_k3 <- train.SOM(k3_test, parallel = F, N.cores = 4, N.steps = 200, N.replicates = 50)
try(SOM_single_k3_2 <- train.SOM(k3_test, max.NA.col = 0.1, max.NA.row = 0.1)) #will fail with message
SOM_single_k3_2 <- train.SOM(k3_test, max.NA.col = 0.1, max.NA.row = 0.2)
try(SOM_single_MORPH <- train.SOM(MORPH, max.NA.row = 0.3)) #will fail with message
try(SOM_single_MORPH <- train.SOM(MORPH, max.NA.row = 0.3, grid.size = 4)) #will fail with message
try(SOM_single_MORPH <- train.SOM(MORPH, max.NA.row = 0.3, grid.size = c(5, 5))) #will fail with message
SOM_single_MORPH <- train.SOM(MORPH, max.NA.row = 0.3, grid.size = c(3, 2))
SOM_single_TestD <- train.SOM(TestD, max.NA.row = 0.3)
try(SOM_multi_MORPH_ENV <- train.SOM(list(MORPH, ENV), max.NA.row = 0.2)) #will fail with message
SOM_multi_MORPH_ENV <- train.SOM(list(MORPH, ENV), max.NA.row = 0.4, learning.rate.tuning = T)
try(SOM_multi_ENV_k3_Alleles <- train.SOM(list(ENV, MORPH, Alleles), max.NA.row = 0.2)) #will fail
try(SOM_multi_ENV_k3_Alleles <- train.SOM(list(ENV, MORPH, Alleles), max.NA.row = 0.5, max.NA.col = 0.1)) #will fail
SOM_multi_ENV_k3_Alleles <- train.SOM(list(ENV, MORPH, Alleles), N.replicates = 10,
                                      message.N.replicates = 5, grid.size = c(2, 2))
SOM_multi_ENV_k3_Alleles_2 <- train.SOM(list(ENV, MORPH, Alleles), save.SOM.results = T, message.N.replicates = 5, grid.multiplier = 2)


## Evaluate SOM object structures for one dataset after SOM training
head(SOM_multi_MORPH_ENV$distance_weights_matrix)
head(SOM_multi_MORPH_ENV$learning_values_list)
head(SOM_multi_MORPH_ENV$input_data_names)
head(SOM_multi_MORPH_ENV$N_steps)
head(SOM_multi_MORPH_ENV$N_replicates)
head(SOM_multi_MORPH_ENV$train.SOM.set.seed.N)
head(SOM_multi_MORPH_ENV$learning_rate_initial)
head(SOM_multi_MORPH_ENV$learning_rate_final)
head(SOM_multi_MORPH_ENV$codebook_vectors)
head(SOM_multi_MORPH_ENV$som_models)
head(SOM_multi_MORPH_ENV$som_models[[1]])
head(SOM_multi_MORPH_ENV$som_models[[1]]$dist.fcts)
head(SOM_multi_MORPH_ENV$som_models[[1]]$radius)
head(SOM_multi_MORPH_ENV$quantization_error)
head(SOM_multi_MORPH_ENV$topographic_error)
head(SOM_multi_MORPH_ENV$train.SOM.args)


## Test cluster SOM function
SOM_single_Alleles <- clustering.SOM(SOM_single_Alleles, max.k = 5, clustering.method = "kmeans+BICelbow")
try(SOM_single_ENV <- clustering.SOM(SOM_single_ENV, max.k = 20, clustering.method = "kmeans+BICelbow")) #will fail with error
SOM_single_ENV <- clustering.SOM(SOM_single_ENV, max.k = 2, clustering.method = "kmeans+BICelbow")
SOM_single_k3 <- clustering.SOM(SOM_single_k3, clustering.method = "kmeans+BICelbow", set.k = 3)
SOM_single_k3_2 <- clustering.SOM(SOM_single_k3_2, clustering.method = "HDBSCAN", max.k = 20)
SOM_single_TestD <- clustering.SOM(SOM_single_TestD, clustering.method = "kmeans+BICelbow", set.k = 2)
try(SOM_multi_ENV_k3_Alleles <- clustering.SOM(SOM_multi_ENV_k3_Alleles, clustering.method = "kmeans+BICelbow")) #will fail with error
SOM_multi_ENV_k3_Alleles <- clustering.SOM(SOM_multi_ENV_k3_Alleles, clustering.method = "hierarchical+DB", max.k = 3)
SOM_multi_MORPH_ENV <- clustering.SOM(SOM_multi_MORPH_ENV, clustering.method = "GMM+BICthreshold", set.k = 3)
SOM_multi_ENV_k3_Alleles_2 <- clustering.SOM(SOM_multi_ENV_k3_Alleles_2, clustering.method = "OPTICS+Silhouette", max.k = 3)


## Evaluate clustering result objects
head(SOM_single_Alleles$cluster_assignment)
SOM_single_Alleles$optim_k_vals
SOM_single_Alleles$optim_k_mean
SOM_single_Alleles$optim_k_summary
SOM_multi_ENV_k3_Alleles$optim_k_summary
SOM_multi_MORPH_ENV$optim_k_summary
SOM_single_k3$median_etasquared_variable_importance
SOM_single_k3$median_map_variance_variable_importance
SOM_single_k3$removed_replicates
SOM_single_k3$quantization.error.quantile
SOM_single_k3$clustering.SOM.args
SOM_single_k3$mean_assignment_margin
SOM_single_k3$mean_normalized_assignment_entropy


## Test learning plot
plot.learning.SOM(SOM_single_k3_2)
plot.learning.SOM(SOM_single_TestD)
plot.learning.SOM(SOM_single_Alleles, top.margin = 4.5)
plot.learning.SOM(SOM_single_k3, lines.alpha = 0.9, lines.thickness = 1)
plot.learning.SOM(SOM_multi_ENV_k3_Alleles, col.pal = viridis::viridis, lines.thickness = 1.5, lines.alpha = 0.2)
plot.learning.SOM(SOM_multi_MORPH_ENV)


## Test layer distance scale plot
try(plot.layer.distance.scale.SOM(SOM_single_Alleles)) #will fail with message (single layer)
plot.layer.distance.scale.SOM(SOM_multi_ENV_k3_Alleles)
plot.layer.distance.scale.SOM(SOM_multi_MORPH_ENV)


## Test K value plot
plot.K.SOM(SOM_single_Alleles)
plot.K.SOM(SOM_single_ENV) #k1-2 only example
plot.K.SOM(SOM_single_k3)
plot.K.SOM(SOM_single_k3_2) #HDBSCAN example
plot.K.SOM(SOM_single_TestD)
plot.K.SOM(SOM_multi_ENV_k3_Alleles, N.axis.labels.BIC.plot = 3) #DB example
plot.K.SOM(SOM_multi_ENV_k3_Alleles_2) #OPTICS+Silhouette example
plot.K.SOM(SOM_multi_MORPH_ENV, col.pal = viridisLite::cividis, round.axis.labels.BIC.plot = 1)


## Test model plot
plot.model.SOM(SOM_single_Alleles, replicate.mode = "average")
plot.model.SOM(SOM_single_Alleles, replicate.mode = "first")
plot.model.SOM(SOM_single_Alleles, replicate.mode = "representative")
plot.model.SOM(SOM_single_ENV)
plot.model.SOM(SOM_single_k3, point.col.clusters = "orange", cluster.shape.clusters = "round", 
               replicate.mode = "average")
plot.model.SOM(SOM_single_k3, replicate.mode = "representative")
plot.model.SOM(SOM_single_k3, replicate.mode = "first")
plot.model.SOM(SOM_single_TestD, col.pal.neighbor.dist = viridis::mako, col.pal.clusters = viridis::magma, point.col.clusters =  "orange")
plot.model.SOM(SOM_multi_ENV_k3_Alleles)
plot.model.SOM(SOM_multi_MORPH_ENV, replicate.mode = "average")
plot.model.SOM(SOM_multi_MORPH_ENV, replicate.mode = "representative")
plot.model.SOM(SOM_multi_MORPH_ENV, replicate.mode = "first")


## Test species coefficient plot
try(plot.structure.SOM(SOM_single_Alleles, bottom.margin = 4)) #will fail with message
plot.structure.SOM(SOM_single_TestD, Individual.labels.font.size = 0.5, bar.border.col = "white")
plot.structure.SOM(SOM_single_k3_2)
plot.structure.SOM(SOM_single_k3_2, sort.by.col = 2)



## Plot maps
plot.map.SOM(SOM_single_ENV, GER_Coordinates, lat.buffer.range = 2)
plot.map.SOM(SOM_single_Alleles, US_Coordinates)
try(plot.map.SOM(SOM_single_Alleles, US_Coordinates, legend.cluster.names = c("Species 1", "Species 2"))) #will fail with message
plot.map.SOM(SOM_single_Alleles, US_Coordinates, legend.cluster.names = c("Species 1: Species name"))


## Additional dataset with NAs and some non-matching rownames to test mapping
ancestry_mat <- matrix(runif(6 * 3), ncol = 3)
rownames(ancestry_mat) <- c("ind1", "ind2", "ind3", "ind4", "ind5", "ind6")
colnames(ancestry_mat) <- paste0("Cluster", 1:3)
Coordinates <- data.frame(
  Latitude = c(48.1, NA, 49.3, 47.5, 50.0),
  Longitude = c(11.6, 8.9, 10.2, 9.7, NA))
rownames(Coordinates) <- c("ind2", "ind3", "ind7", "ind6", "ind1") #ind7" does not match ancestry; "ind4", "ind5" missing
plot.map.SOM(SOM.output = list(ancestry_matrix = ancestry_mat, input_data_names = c("testdata")),
             Coordinates,
             north.arrow.position = c(0.04, 0.85),
             north.arrow.length = 0.4,
             north.arrow.N.position = 0.15)


## Test variable importance plot
try(plot.variable.importance.SOM(SOM_single_Alleles)) #will fail with error message
plot.variable.importance.SOM(SOM_single_Alleles, mode = "Map.variance") 
plot.variable.importance.SOM(SOM_single_TestD, mode = "Cluster.separation")
plot.variable.importance.SOM(SOM_single_TestD, mode = "Map.variance")
try(plot.variable.importance.SOM(SOM_single_ENV)) #will fail
plot.variable.importance.SOM(SOM_single_k3)
plot.variable.importance.SOM(SOM_single_k3_2)
plot.variable.importance.SOM(SOM_multi_ENV_k3_Alleles, top.margin = 2, bar.label.font.size = 0.4, 
                             bars.threshold.N = 200, mode = "Map.variance")
plot.variable.importance.SOM(SOM_multi_MORPH_ENV, top.margin = 6, bar.label.font.size = 0.5, 
                             mode = "Map.variance")


## Test layer importance layer plot
plot.layer.importance.varimp.SOM(SOM_single_Alleles) 
plot.layer.importance.varimp.SOM(SOM_single_TestD)
plot.layer.importance.varimp.SOM(SOM_single_TestD, add.boxplot.whiskers = F)
plot.layer.importance.varimp.SOM(SOM_single_ENV)
plot.layer.importance.varimp.SOM(SOM_single_k3)
plot.layer.importance.varimp.SOM(SOM_single_k3_2)
plot.layer.importance.varimp.SOM(SOM_multi_ENV_k3_Alleles)
plot.layer.importance.varimp.SOM(SOM_multi_ENV_k3_Alleles, col.pal = viridis::rocket)
plot.layer.importance.varimp.SOM(SOM_multi_MORPH_ENV)


## Test layer importance layer plot
try(plot.layer.importance.leaveoneout.SOM(SOM_single_Alleles)) #will fail because of single layer
plot.layer.importance.leaveoneout.SOM(SOM_single_TestD)
plot.layer.importance.leaveoneout.SOM(SOM_multi_ENV_k3_Alleles, 
                                      add.points = F, col.pal = viridis::rocket, save = T)
plot.layer.importance.leaveoneout.SOM(SOM_multi_MORPH_ENV)
