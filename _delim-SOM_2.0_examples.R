################################################################################
#### Set environment
################################################################################

## Install and load packages
CRAN_packages <- c("adegenet", 
                   "ape", 
                   "dplyr", 
                   "poppr", 
                   "readr", 
                   "sf", 
                   "stringr", 
                   "tibble", 
                   "vcfR") #CRAN packages
Bioconductor_packages <- "SeqArray" #Bioconductor package
for (p in CRAN_packages) if (!requireNamespace(p, quietly = TRUE)) install.packages(p) #install missing CRAN packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") #install BiocManager
if (!requireNamespace(Bioconductor_packages, quietly = TRUE)) BiocManager::install(Bioconductor_packages) #install missing Bioconductor package
lapply(c(CRAN_packages, Bioconductor_packages), library, character.only = TRUE) #load all packages

rm(list = ls()) #clear environment
setwd("C:/Users/danie/Desktop/PhD research/SOM package")
source("delim-SOM v2 functions.R")



################################################################################
#### Test SOM functions with small simulated datasets
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
n_k3_test <- 40  #number of traits
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


## Test_D dataset
n <- 10
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
SOM_single_Alleles <- train.SOM(Alleles, max.NA.row = 0.3)
SOM_single_ENV <- train.SOM(ENV, max.NA.col = 0.05, max.NA.row = 0.4) #will fail with error
SOM_single_ENV <- train.SOM(ENV, max.NA.col = 0.15, max.NA.row = 0.4) #will fail with error
SOM_single_ENV <- train.SOM(ENV, max.NA.col = 0.3, max.NA.row = 0.2) #will fail with error
SOM_single_ENV <- train.SOM(ENV, max.NA.col = 0.3, max.NA.row = 0.2, grid.multiplier = 3)
SOM_single_k3 <- train.SOM(k3_test, parallel = T, N_cores = 4, N.steps = 200, N.replicates = 50) #test parallel version
SOM_single_k3_2 <- train.SOM(k3_test, max.NA.col = 0.1, max.NA.row = 0.1) #will fail with message
SOM_single_k3_2 <- train.SOM(k3_test, max.NA.col = 0.1, max.NA.row = 0.2)
SOM_single_MORPH <- train.SOM(MORPH, max.NA.row = 0.3) #will fail with message
SOM_single_MORPH <- train.SOM(MORPH, max.NA.row = 0.3, grid.size = 4) #will fail with message
SOM_single_MORPH <- train.SOM(MORPH, max.NA.row = 0.3, grid.size = c(2, 3)) #will fail with message
SOM_single_TestD <- train.SOM(TestD, max.NA.row = 0.3) #will fail with message
SOM_single_TestD <- train.SOM(TestD, max.NA.row = 0.3, grid.multiplier = 2)
SOM_multi_MORPH_ENV <- train.SOM(list(MORPH, ENV), max.NA.row = 0.2) #will fail with message
SOM_multi_MORPH_ENV <- train.SOM(list(MORPH, ENV), max.NA.row = 0.3) #will fail with message
SOM_multi_MORPH_ENV <- train.SOM(list(MORPH, ENV), max.NA.row = 0.4)
SOM_multi_ENV_k3_Alleles <- train.SOM(list(ENV, MORPH, Alleles), max.NA.row = 0.2) #will fail
SOM_multi_ENV_k3_Alleles <- train.SOM(list(ENV, MORPH, Alleles), max.NA.row = 0.5, max.NA.col = 0.1) #will fail
SOM_multi_ENV_k3_Alleles <- train.SOM(list(ENV, MORPH, Alleles), learning.rate.tuning = F, message.N.replicates = 5, grid.size = c(2, 2))
SOM_multi_ENV_k3_Alleles_2 <- train.SOM(list(ENV, MORPH, Alleles), learning.rate.tuning = F, save.SOM.results = T, message.N.replicates = 5)


## Evaluate SOM object structures for one dataset after SOM training
head(SOM_multi_MORPH_ENV$distance_matrix)
head(SOM_multi_MORPH_ENV$learning_values_list)
head(SOM_multi_MORPH_ENV$input_data_names)
head(SOM_multi_MORPH_ENV$N_steps)
head(SOM_multi_MORPH_ENV$N_replicates)
head(SOM_multi_MORPH_ENV$learning_rate_initial)
head(SOM_multi_MORPH_ENV$learning_rate_final)
head(SOM_multi_MORPH_ENV$codebook_vectors)
head(SOM_multi_MORPH_ENV$median_variable_importance)
head(SOM_multi_MORPH_ENV$som_models)
head(SOM_multi_MORPH_ENV$som_models[[1]])


## Test cluster SOM function
SOM_single_Alleles <- clustering.SOM(SOM_single_Alleles, max.k = 5, clustering.method = "kmeans+BICjump+hierarch")
SOM_single_ENV <- clustering.SOM(SOM_single_ENV, max.k = 20, clustering.method = "kmeans+BICjump+hierarch") #will fail with error
SOM_single_ENV <- clustering.SOM(SOM_single_ENV, max.k = 8, clustering.method = "kmeans+BICjump+hierarch")
SOM_single_k3 <- clustering.SOM(SOM_single_k3, clustering.method = "kmeans+BICjump+hierarch", set.k = 5)
SOM_single_k3_2 <- clustering.SOM(SOM_single_k3_2, clustering.method = "kmeans+BICjump+hierarch", max.k = 20, pca.codebooks = F)
SOM_single_TestD <- clustering.SOM(SOM_single_TestD, clustering.method = "kmeans+BICjump+hierarch", max.k = 3, pca.codebooks = F)
SOM_multi_ENV_k3_Alleles <- clustering.SOM(SOM_multi_ENV_k3_Alleles, clustering.method = "kmeans+BICjump+hierarch") #will fail with error
SOM_multi_ENV_k3_Alleles <- clustering.SOM(SOM_multi_ENV_k3_Alleles, clustering.method = "kmeans+BICjump+hierarch", max.k = 2)
SOM_multi_MORPH_ENV <- clustering.SOM(SOM_multi_MORPH_ENV, clustering.method = "kmeans+BICjump+hierarch", set.k = 3)


## Evaluate clustering result objects
head(SOM_single_Alleles$cluster_assignment)
SOM_single_Alleles$optim_k_vals
SOM_single_Alleles$optim_k_mean
SOM_single_Alleles$optim_k_summary
SOM_multi_ENV_k3_Alleles$optim_k_summary
SOM_multi_MORPH_ENV$optim_k_summary


## Test learning plot
plot.Learning.SOM(SOM_single_k3_2)
plot.Learning.SOM(SOM_single_TestD)
plot.Learning.SOM(SOM_single_Alleles, margin.top = 4.5)
plot.Learning.SOM(SOM_single_k3, lines.alpha = 0.9, lines.thickness = 1)
plot.Learning.SOM(SOM_single_MORPH, save = F, title = "Training progress")
plot.Learning.SOM(SOM_multi_ENV_k3_Alleles, col.pal = mako, lines.thickness = 1.5, lines.alpha = 1)
plot.Learning.SOM(SOM_multi_MORPH_ENV)


## Test layer weights plot
plot.Layers.SOM(SOM_single_Alleles) #will fail with message
plot.Layers.SOM(SOM_multi_ENV_k3_Alleles)
plot.Layers.SOM(SOM_multi_MORPH_ENV)


## Test K value plot
plot.K.SOM(SOM_single_Alleles)
plot.K.SOM(SOM_single_ENV)
plot.K.SOM(SOM_single_k3)
plot.K.SOM(SOM_single_k3_2)
plot.K.SOM(SOM_single_TestD)
plot.K.SOM(SOM_multi_ENV_k3_Alleles, N.axis.labels.BIC.plot = 3)
plot.K.SOM(SOM_multi_ENV_k3_Alleles_2) #will fail with error
plot.K.SOM(SOM_multi_MORPH_ENV, col.pal = cividis, round.axis.labels.BIC.plot = 1)


## Test model plot
plot.Model.SOM(SOM_single_Alleles)
plot.Model.SOM(SOM_single_ENV, boundary.col.clusters = "orange", boundary.lwd.clusters = 5, save = T)
plot.Model.SOM(SOM_single_k3, point.col.clusters = "black", cluster.shape.clusters = "round")
plot.Model.SOM(SOM_single_TestD, col.pal.neighbor.dist = mako, col.pal.clusters = magma, point.col.clusters =  "orange")
plot.Model.SOM(SOM_multi_ENV_k3_Alleles)
plot.Model.SOM(SOM_multi_MORPH_ENV)


## Test species coefficient plot
plot.Structure.SOM(SOM_single_Alleles, margin.bottom = 4)
plot.Structure.SOM(SOM_single_ENV, Individual.labels.font.size = 0.7, bar.border.col = "white", bar.border.lwd = 1)
plot.Structure.SOM(SOM_single_TestD, Individual.labels.font.size = 0.9, bar.border.col = "black")
plot.Structure.SOM(SOM_multi_ENV_k3_Alleles)
plot.Structure.SOM(SOM_multi_MORPH_ENV)


## Plot maps
plot.Map.SOM(SOM_single_ENV, GER_Coordinates, lat.buffer.range = 2)
plot.Map.SOM(SOM_single_Alleles, US_Coordinates)
plot.Map.SOM(SOM_single_Alleles, US_Coordinates, legend.cluster.names = c("Species 1", "Species 2")) #will fail with message
plot.Map.SOM(SOM_single_Alleles, US_Coordinates, legend.cluster.names = c("Species 1", "Species 2", "Species 3", "Species 4"))


## Additional dataset with NAs and some non-matching rownames to test mapping
ancestry_mat <- matrix(runif(6 * 3), ncol = 3)
rownames(ancestry_mat) <- c("ind1", "ind2", "ind3", "ind4", "ind5", "ind6")
colnames(ancestry_mat) <- paste0("Cluster", 1:3)
Coordinates <- data.frame(
  Latitude = c(48.1, NA, 49.3, 47.5, 50.0),
  Longitude = c(11.6, 8.9, 10.2, 9.7, NA))
rownames(Coordinates) <- c("ind2", "ind3", "ind7", "ind6", "ind1") #ind7" does not match ancestry; "ind4", "ind5" missing
plot.Map.SOM(SOM.output = list(ancestry_matrix = ancestry_mat, input_data_names = c("testdata")), 
             Coordinates,
             north.arrow.position = c(0.04, 0.85),
             north.arrow.length = 0.4,
             north.arrow.N.position = 0.15)


## Test variable importance plot
plot.Variable.Importance.SOM(SOM_single_Alleles)
plot.Variable.Importance.SOM(SOM_single_TestD)
plot.Variable.Importance.SOM(SOM_single_ENV)
plot.Variable.Importance.SOM(SOM_single_k3)
plot.Variable.Importance.SOM(SOM_single_k3_2)
plot.Variable.Importance.SOM(SOM_multi_ENV_k3_Alleles, col.pal = inferno)
plot.Variable.Importance.SOM(SOM_multi_ENV_k3_Alleles, top.margin = 6, bar.label.font.size = 0.4, bars.threshold.N = 200)
plot.Variable.Importance.SOM(SOM_multi_MORPH_ENV, top.margin = 6, bar.label.font.size = 0.5)



################################################################################
#### Empirical example: original Monticola71 salamander data (Pyron et al. 2023)
################################################################################

## Import and process genetic SNP data
Monticola71original_genind <- read.structure(file = "Test data/seal_in_c90.str",
                                             n.ind = 71,
                                             n.loc = 7809,
                                             onerowperind = F,
                                             col.lab = 1,
                                             col.pop = 0,
                                             col.others = 0,
                                             row.marknames = 0,
                                             NA.char = -9) #read in str
Monticola71original_SNP <- process.SNP.data.SOM(genind.input = Monticola71original_genind,
                                                missing.loci.cutoff.lenient = 0.5,
                                                missing.loci.cutoff.final = 0.2,
                                                missing.individuals.cutoff = 0.35,
                                                singleton.loci.filter = TRUE,
                                                invariant.loci.filter = TRUE)
ncol(Monticola71original_SNP) #number of SNPs: 2290
nrow(Monticola71original_SNP) #number of samples: 66



## Train and cluster SOM
Monticola71original_SOM_data <- list(SNP = Monticola71original_SNP)
Monticola71original_SOM <- train.SOM(input_data = Monticola71original_SOM_data, 
                                     N.steps = 150,
                                     N.replicates = 100, 
                                     learning.rate.tuning = T, 
                                     parallel = T,
                                     save.SOM.results = T,
                                     overwrite.SOM.results = F,
                                     N_cores = 5)
Monticola71original_SOM <- clustering.SOM(SOM.output = Monticola71original_SOM, 
                                          clustering.method = "kmeans+BICjump+hierarch",
                                          max.k = 7)


## Evaluate and plot results
Monticola71original_SOM$optim_k_summary
plot.Learning.SOM(Monticola71original_SOM)
plot.Model.SOM(Monticola71original_SOM)
plot.Structure.SOM(Monticola71original_SOM, sort.by.col = 2)
plot.Variable.Importance.SOM(Monticola71original_SOM)
head(sort(Monticola71original_SOM$median_variable_importance[[1]], decreasing = T))



################################################################################
#### Empirical example: Desmognathus salamanders - Monticola71 (original GBS data from Pyron et al. 2023 with updated environmental variables)
################################################################################

## Import sample data
Monticola71_data <- read.csv(file = "Test data/monticola71.csv",
                             row.names = 1, 
                             header = TRUE,
                             colClasses = c(huc2 = "character",
                                            huc4 = "character",
                                            huc6 = "character",
                                            huc8 = "character",
                                            huc10 = "character",
                                            huc12 = "character"))


## Import and process genetic SNP data
Monticola71_SNP <- process.SNP.data.SOM(vcf.path = "Test data/Monticola71.vcf.gz", #filter loci and individuals and create SNP matrix dataframe
                                        missing.loci.cutoff.lenient = 0.4,
                                        missing.loci.cutoff.final = 0.2,
                                        missing.individuals.cutoff = 0.4)
Monticola71_snp_to_sample <- Monticola71_data$Sample[match(rownames(Monticola71_SNP), rownames(Monticola71_data))] #returns RAP**** names
rownames(Monticola71_SNP) <- Monticola71_snp_to_sample #rename SNP matrix to RAP codes
ncol(Monticola71_SNP) # number of loci: 2669
nrow(Monticola71_SNP) # number of samples: 64


## Create spatial dataset with coordinates and elevation
Monticola71_spatial <- Monticola71_data[, c("Lat", "Long", "Elev"), drop = FALSE] #extract numeric spatial data
Monticola71_spatial <- Monticola71_spatial %>% rename(Latitude = Lat, Longitude = Long, Elevation = Elev)
rownames(Monticola71_spatial) <- Monticola71_data$Sample #assign rownames
nrow(Monticola71_spatial) #number of samples: 71


## Create environmental dataset, including binary watershed variables (other variables extracted and processed by separate R script based on coordinates)
Monticola71_environmental <- read.csv("Test data/Monticola71_environmental.csv", row.names = 1, header = TRUE) #read CSV with rownames
Monticola71_environmental <- Monticola71_environmental %>% 
  dplyr::select(-Latitude, -Longitude, -Elevation) #remove spatial variables
Monticola71_environmental <- as.data.frame(lapply(Monticola71_environmental, as.numeric)) #ensure numeric
rownames(Monticola71_environmental) <- Monticola71_data$Sample #assign rownames
Monticola71_watershed_binary <- make.cols.binary.SOM(dataframe = Monticola71_data, #make watershed columns binary
                                                     make.binary.cols = c("lvl4", "huc2", "huc4", "huc6", "huc8", "huc10", "huc12")) #select watershed columns
Monticola71_watershed_binary <- Monticola71_watershed_binary[Monticola71_data$Sample, ] #ensure matched order
Monticola71_environmental <- cbind(Monticola71_environmental, Monticola71_watershed_binary) #append to environmental data
Monticola71_environmental <- remove.lowCV.multicollinearity.SOM(Monticola71_environmental, #remove highly correlated and low-variance variables
                                                                CV.threshold = 0.05,
                                                                cor.threshold = 0.9)
ncol(Monticola71_environmental) #number of variables: 160
nrow(Monticola71_environmental) #number of samples: 71


## Create morphological trait dataset
Monticola71_trait_names <- c("SVL", "TL", "AG", "CW", "FL", "HL", "SG", "TW", "TO", "FI", "HW", "ED", "IN", "ES", "ON", "IO", "IC")
Monticola71_trait_data <- Monticola71_data[, Monticola71_trait_names] #extract variables
Monticola71_log_traits <- log(Monticola71_trait_data) #log-transform all traits
Monticola71_filtered_log_traits <- remove.lowCV.multicollinearity.SOM(Monticola71_log_traits, #filter log-transformed traits by CV and correlation, excluding SVL from removal
                                                                      CV.threshold = 0.05,
                                                                      cor.threshold = 0.9,
                                                                      exclude.cols = "SVL")
Monticola71_SVL <- Monticola71_filtered_log_traits[, "SVL"] #extract SVL and residualize all others
Monticola71_residuals_mat <- sapply(colnames(Monticola71_filtered_log_traits)[colnames(Monticola71_filtered_log_traits) != "SVL"], function(trait) {resid(lm(Monticola71_filtered_log_traits[, trait] ~ Monticola71_SVL))}) #regress each trait on SVL and store residuals
rownames(Monticola71_filtered_log_traits) <- Monticola71_data$Sample #set rownames for log-transformed traits
rownames(Monticola71_residuals_mat) <- Monticola71_data$Sample #set rownames for residualized traits
Monticola71_morphology <- as.data.frame(cbind(SVL = Monticola71_SVL, Monticola71_residuals_mat)) #combine log(SVL) and residuals
ncol(Monticola71_morphology) #number of traits: 7
nrow(Monticola71_morphology) #number of samples: 71


## Train and cluster SOM
SOM_data_Monticola71 <- list(SNP = Monticola71_SNP,
                             Spatial = Monticola71_spatial,
                             Environmental = Monticola71_environmental,
                             Morphology = Monticola71_morphology)
SOM_Monticola71 <- train.SOM(input_data = SOM_data_Monticola71,
                             N.steps = 2,
                             N.replicates = 2,
                             parallel = F, 
                             N_cores = 5,
                             learning.rate.tuning = F,
                             max.NA.row = 0.2, 
                             max.NA.col = 0.15, 
                             save.SOM.results = F,
                             overwrite.SOM.results = F)
SOM_Monticola71 <- clustering.SOM(SOM.output = SOM_Monticola71, 
                                  clustering.method = "kmeans+BICjump+hierarch", 
                                  max.k = 10)


## Evaluate and plot results
plot.Learning.SOM(SOM_Monticola71)
plot.Layers.SOM(SOM_Monticola71)
plot.K.SOM(SOM_Monticola71)
plot.Model.SOM(SOM_Monticola71)
plot.Structure.SOM(SOM_Monticola71)
plot.Map.SOM(SOM.output = SOM_Monticola71, 
             Coordinates = Monticola71_spatial[, c("Latitude", "Longitude")], 
             USA.add.counties = T,
             scale.position = c(0.8, 0.05))
plot.Variable.Importance.SOM(SOM_Monticola71,
                             left.margin = 5.8,
                             bar.label.font.size = 0.4, 
                             importance.threshold = 0.1)

head(round(sort(SOM_Monticola71$median_variable_importance[[1]], decreasing = T), 2), 800)
head(round(sort(SOM_Monticola71$median_variable_importance[[2]], decreasing = T), 2), 3)
head(round(sort(SOM_Monticola71$median_variable_importance[[3]], decreasing = T), 2), 40)
head(round(sort(SOM_Monticola71$median_variable_importance[[4]], decreasing = T), 2), 10)

tail(round(sort(SOM_Monticola71$median_variable_importance[[3]], decreasing = T), 2), 100)




################################################################################
## Empirical example: aeneus56 - original GBS data from Pyron et al. 2024 with updated environmental variables
################################################################################

## Read in sample data
Aeneus_data <- read.csv(file = "Test data/aeneus56.csv",
                        row.names = 1,
                        header = T, colClasses = c(huc2 = "character",
                                                   huc4 = "character",
                                                   huc6 = "character",
                                                   huc8 = "character",
                                                   huc10 = "character",
                                                   huc12 = "character"))


## Import and process genetic SNP data
Aeneus_SNP <- process.SNP.data.SOM(vcf.path = "Test data/aeneus56.vcf.gz", #filter loci and individuals and create SNP matrix dataframe
                                   missing.loci.cutoff.lenient = 0.4,
                                   missing.loci.cutoff.final = 0.2,
                                   missing.individuals.cutoff = 0.3)
rownames(Aeneus_SNP) <- Aeneus_data$Sample[match(rownames(Aeneus_SNP), rownames(Aeneus_data))] #rename alleles
ncol(Aeneus_SNP) #number of loci: 1626
nrow(Aeneus_SNP) #number of samples: 47


## Create spatial dataset with coordinates and elevation
Aeneus_spatial <- Aeneus_data[, c("Lat", "Long", "Elev")] #extract coordinates and elevation
Aeneus_spatial <- Aeneus_spatial %>% rename(Latitude = Lat, Longitude = Long, Elevation = Elev) #rename variables
rownames(Aeneus_spatial) <- Aeneus_data$Sample #assign rownames
nrow(Aeneus_spatial) #number of samples: 56


## Create environmental dataset, including binary watershed variables (other variables extracted and processed by separate R script based on coordinates)
Aeneus_environmental <- read.csv("Test data/Aeneus_environmental.csv", row.names = 1, header = TRUE) #read CSV with rownames
Aeneus_environmental <- Aeneus_environmental %>% 
  dplyr::select(-Latitude, -Longitude, -Elevation) #remove spatial variables
Aeneus_environmental <- as.data.frame(lapply(Aeneus_environmental, as.numeric)) #ensure numeric
rownames(Aeneus_environmental) <- Aeneus_data$Sample #assign rownames
Aeneus_watershed_binary <- make.cols.binary.SOM(dataframe = Aeneus_data, #make watershed columns binary
                                                make.binary.cols = c("lvl4", "huc2", "huc4", "huc6", "huc8", "huc10", "huc12")) #select watershed columns
Aeneus_watershed_binary <- Aeneus_watershed_binary[Aeneus_data$Sample, ] #ensure matched order
Aeneus_environmental <- cbind(Aeneus_environmental, Aeneus_watershed_binary) #append to environmental data
Aeneus_environmental <- remove.lowCV.multicollinearity.SOM(df = Aeneus_environmental, #remove highly correlated and low-variance variables
                                                           CV.threshold = 0.05,
                                                           cor.threshold = 0.9)
ncol(Aeneus_environmental) #number of variables: 34
nrow(Aeneus_environmental) #number of samples: 56


# Create morphological trait dataset
Aeneus_trait_names <- c("SVL", "TL", "AG", "CW", "FL", "HL", "SG", "TW", "TO", "FI", "HW", "ED", "IN", "ES", "ON", "IO", "IC") #define trait columns
Aeneus_trait_data <- Aeneus_data[, Aeneus_trait_names] #extract variables
rownames(Aeneus_trait_data) <- Aeneus_data$Sample #assign sample IDs as rownames
Aeneus_trait_data <- Aeneus_trait_data[rowSums(!is.na(Aeneus_trait_data)) > 0, ] #remove samples with only NA values
Aeneus_log_traits <- log(Aeneus_trait_data) #log-transform all traits
Aeneus_filtered_log_traits <- remove.lowCV.multicollinearity.SOM( #filter log-transformed traits by CV and correlation, excluding SVL from removal
  Aeneus_log_traits,
  CV.threshold = 0.05,
  cor.threshold = 0.9,
  exclude.cols = "SVL")
rownames(Aeneus_filtered_log_traits) <- rownames(Aeneus_trait_data) #set rownames for filtered log traits
Aeneus_SVL <- Aeneus_filtered_log_traits[, "SVL"] #extract SVL and residualize all others
Aeneus_residuals_mat <- sapply(colnames(Aeneus_filtered_log_traits)[colnames(Aeneus_filtered_log_traits) != "SVL"],
                               function(trait) resid(lm(Aeneus_filtered_log_traits[, trait] ~ Aeneus_SVL, na.action = na.exclude))) #regress each trait on SVL and retain NA alignment
rownames(Aeneus_residuals_mat) <- rownames(Aeneus_filtered_log_traits) #assign rownames to residual matrix
Aeneus_morphology <- as.data.frame(cbind(SVL = Aeneus_SVL, Aeneus_residuals_mat)) #combine log(SVL) and residuals
Aeneus_morphology <- Aeneus_morphology[rownames(Aeneus_trait_data), ] #keep only samples with trait data
ncol(Aeneus_morphology) #number of traits: 11
nrow(Aeneus_morphology) #number of samples: 45


## Train and cluster SOM
SOM_Aeneus_SOM_data <- list(Alleles = Aeneus_SNP,
                            Spatial = Aeneus_spatial,
                            Environmental = Aeneus_environmental,
                            Morphology = Aeneus_morphology)
SOM_Aeneus <- train.SOM(input_data = SOM_Aeneus_SOM_data,
                        N.steps = 2,
                        N.replicates = 2,
                        parallel = F, 
                        N_cores = 5,
                        grid.multiplier = 4,
                        learning.rate.tuning = F,
                        max.NA.row = 0.2, 
                        max.NA.col = 0.15, 
                        save.SOM.results = F,
                        overwrite.SOM.results = F)
SOM_Aeneus <- clustering.SOM(SOM.output = SOM_Aeneus, 
                             clustering.method = "kmeans+BICjump+hierarch", 
                             max.k = 10)


## Evaluate and plot results
plot.Learning.SOM(SOM_Aeneus)
plot.Layers.SOM(SOM_Aeneus)
plot.K.SOM(SOM_Aeneus)
plot.Model.SOM(SOM_Aeneus)
plot.Structure.SOM(SOM_Aeneus)
plot.Map.SOM(SOM.output = SOM_Aeneus, 
             Coordinates = Aeneus_spatial[, c("Latitude", "Longitude")], 
             USA.add.counties = T,
             north.arrow.position = c(0.05, 0.9),
             north.arrow.length = 0.4,
             north.arrow.N.position = 0.15,
             scale.position = c(0.79, 0.05))
plot.Variable.Importance.SOM(SOM_Aeneus,
                             left.margin = 5.8,
                             bar.label.font.size = 0.4, 
                             importance.threshold = 0.1)

head(round(sort(SOM_Aeneus$median_variable_importance[[1]], decreasing = T), 2), 400)
head(round(sort(SOM_Aeneus$median_variable_importance[[2]], decreasing = T), 2), 3)
head(round(sort(SOM_Aeneus$median_variable_importance[[3]], decreasing = T), 2), 40)
head(round(sort(SOM_Aeneus$median_variable_importance[[4]], decreasing = T), 2), 10)

tail(round(sort(SOM_Aeneus$median_variable_importance[[3]], decreasing = T), 2), 100)



################################################################################
#### Pocillopora corals in Indo-Pacific (Oury et al 2023)
################################################################################

## Import and process genetic data
Pocillopora_vcf_file <- "Test data/Oury et al 2023/Pocillopora_361ADN_1559SNP.vcf" #VCF file path
Pocillopora_gds_file <- "Test data/Oury et al 2023/Pocillopora.gds" #GDS file path
seqVCF2GDS(Pocillopora_vcf_file, Pocillopora_gds_file, storage.option = "LZ4_RA.max") #convert VCF to GDS
Pocillopora_gds <- seqOpen(Pocillopora_gds_file) #open GDS file
Pocillopora_geno <- seqGetData(Pocillopora_gds, "genotype") #get genotype array
Pocillopora_SNP <- Pocillopora_geno[1, , ] + Pocillopora_geno[2, , ] #sum allele counts
rownames(Pocillopora_SNP) <- seqGetData(Pocillopora_gds, "sample.id") #assign rownames
colnames(Pocillopora_SNP) <- paste0("SNP", seq_len(ncol(Pocillopora_SNP))) #assign colnames
Pocillopora_SNP <- as.data.frame(Pocillopora_SNP) #convert to data frame
seqClose(Pocillopora_gds) #close GDS connection
Pocillopora_genotypes <- as.data.frame(matrix(ifelse(is.na(unlist(Pocillopora_SNP)), NA, ifelse(unlist(Pocillopora_SNP) == 0, "A/A", ifelse(unlist(Pocillopora_SNP) == 1, "A/B", ifelse(unlist(Pocillopora_SNP) == 2, "B/B", NA)))), nrow = nrow(Pocillopora_SNP), dimnames = dimnames(Pocillopora_SNP))) #convert dosages to genotype strings
Pocillopora_genind <- df2genind(Pocillopora_genotypes, sep = "/", ncode = 1, ploidy = 2) #convert to genind
Pocillopora <- process.SNP.data.SOM(genind.input = Pocillopora_genind, #filter loci and individuals and create SNP matrix dataframe
                                    missing.loci.cutoff.lenient = 0.4,
                                    missing.loci.cutoff.final = 0.2,
                                    missing.individuals.cutoff = 0.5)
ncol(Pocillopora_SNP) #number of loci: 1559
nrow(Pocillopora_SNP) #number of samples: 361


## Import and process morphology dataset
Pocillopora_morphology <- read_delim(file = "Test data/Oury et al 2023/Micromorphometry_Pocillopora_170ind.csv", #import csv
                                     delim = ";", 
                                     quote = "\"", 
                                     escape_double = TRUE,
                                     locale = locale(decimal_mark = ","), 
                                     trim_ws = TRUE)
Pocillopora_morphology <- as.data.frame(Pocillopora_morphology) #make dataframe
rownames(Pocillopora_morphology) <- Pocillopora_morphology$Sample_Name #set rownames
Pocillopora_morphology$GSH <- NULL #remove GSH column
Pocillopora_morphology$Sample_Name <- NULL #remove sample name column
Pocillopora_trait_names <- c( #shorten trait names
  "Max_calice_diameter_1", #(v1) Maximum calice diameter
  "Max_calice_diameter_2", #(v2) Perp. diameter to v1
  "Distance_corallite", #(v3) Center-to-center distance
  "Distance_denticles", #(v4) Denticle spacing
  "Height_septa", #(v5) Septa height/teeth
  "Max_columella_diameter_1", #(v6) Maximum columella diameter
  "Max_columella_diameter_2", #(v7) Perp. columella diameter
  "Shape_septa", #(v8) Septa shape
  "Shape_columella") #(v9) Columella shape
colnames(Pocillopora_morphology) <- Pocillopora_trait_names
Pocillopora_morphology <- make.cols.binary.SOM(dataframe = Pocillopora_morphology, #make columella and septa shapes binary
                                               make.binary.cols = c("Shape_columella", "Shape_septa"),
                                               append.to.original = TRUE)


## Import csv file with multiple traits and meta data
Pocillopora_multiple_traits <- read_delim(file = "Test data/Oury et al 2023/DB_Pocillopora_genomic_364ind.csv",
                                          delim = ";", 
                                          quote = "\"", 
                                          escape_double = TRUE,
                                          show_col_types = FALSE,
                                          locale = locale(decimal_mark = ","),
                                          trim_ws = TRUE) %>%
  as.data.frame()

Pocillopora_species_names <- Pocillopora_multiple_traits %>% 
  distinct(Sample_Name, .keep_all = TRUE) %>%
  column_to_rownames("Sample_Name") %>%
  dplyr::select(Locality, SSH, GSH, PSH)
Pocillopora_species_names[Pocillopora_species_names == "-"] <- NA
Pocillopora_species_names[Pocillopora_species_names == "?"] <- NA


## Microsatellite genotypes
Pocillopora_microsat_cols <- grep("\\.(1|2)$", names(Pocillopora_multiple_traits), value = TRUE)
Pocillopora_microsatellites <- Pocillopora_multiple_traits %>%
  dplyr::select(Sample_Name, all_of(Pocillopora_microsat_cols)) %>%
  mutate(across(-Sample_Name, ~ na_if(., "-"))) %>%
  mutate(across(-Sample_Name, ~ na_if(., "?"))) %>%
  mutate(na_count = rowSums(is.na(across(-Sample_Name)))) %>%
  group_by(Sample_Name) %>%
  slice_min(na_count, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-na_count) %>%
  column_to_rownames("Sample_Name") %>%
  mutate(across(everything(), as.numeric))


## Haplotype markers (ORF and PocHistone) for symbiosis as binary
Pocillopora_symbiosis_haplotypes <- Pocillopora_multiple_traits %>% #copy original
  mutate(across(c(ORF, PocHistone), ~ na_if(.x, "-"))) %>% #replace "-" with NA
  mutate(across(c(ORF, PocHistone), ~ na_if(.x, "?"))) %>% #replace "?" with NA
  dplyr::select(Sample_Name, ORF, PocHistone) %>% #select only relevant columns
  mutate(na_count = rowSums(is.na(across(c(ORF, PocHistone))))) %>% #count NAs
  group_by(Sample_Name) %>% slice_min(na_count, with_ties = FALSE) %>% ungroup() %>% #keep best per sample
  dplyr::select(-na_count) #drop NA count
Pocillopora_symbiosis_haplotypes[] <- lapply(Pocillopora_symbiosis_haplotypes, function(x) {
  if (is.character(x)) x <- as.factor(x) #convert character to factor
  if (is.factor(x)) x <- droplevels(x) #drop unused levels
  return(x)
})
Pocillopora_symbiosis_haplotypes <- as.data.frame(Pocillopora_symbiosis_haplotypes) #ensure data.frame
rownames(Pocillopora_symbiosis_haplotypes) <- Pocillopora_symbiosis_haplotypes$Sample_Name #set rownames
Pocillopora_symbiosis_haplotypes$Sample_Name <- NULL #remove Sample_Name
Pocillopora_symbiosis_haplotypes <- make.cols.binary.SOM( #convert ORF and PocHistone to binary variables
  dataframe = Pocillopora_symbiosis_haplotypes,
  make.binary.cols = c("ORF", "PocHistone"),
  append.to.original = TRUE)
colnames(Pocillopora_symbiosis_haplotypes) <- colnames(Pocillopora_symbiosis_haplotypes) %>%
  gsub("^ORF_NA$|^NA$", "ORF_missing", .) %>% #handle ORF missing
  gsub("^ORF_", "", .) %>% #remove ORF prefix
  gsub("^PocHistone_NA$", "Poc_missing", .) %>% #handle Poc missing
  gsub("^PocHistone_", "Poc", .) #rename PocHistone
Pocillopora_symbiosis_haplotypes <- Pocillopora_symbiosis_haplotypes %>% #drop missing flags
  dplyr::select(!matches("(^missing$|^Poc_missing$)"))


## Biogeographic region as binary
Pocillopora_biogeography <- Pocillopora_multiple_traits %>% 
  dplyr::select(Sample_Name, Province) %>% 
  mutate(Province = na_if(Province, "-")) %>% 
  mutate(na_count = is.na(Province)) %>% 
  group_by(Sample_Name) %>% 
  slice_min(na_count, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  dplyr::select(-na_count) %>% 
  column_to_rownames("Sample_Name")
Pocillopora_biogeography <- make.cols.binary.SOM(dataframe = Pocillopora_biogeography, #convert Province to binary matrix
                                                 make.binary.cols = "Province")


## Add morphotype as binary to morphology dataset
Pocillopora_morpho_map <- Pocillopora_multiple_traits %>% #extract Sample_Name and Morphotype
  dplyr::select(Sample_Name, Morphotype) %>% #select relevant columns
  mutate(na_count = is.na(Morphotype)) %>% #count NA entries
  group_by(Sample_Name) %>% #group by sample
  slice_min(na_count, n = 1, with_ties = FALSE) %>% #keep most complete entry
  ungroup() %>% #ungroup
  dplyr::select(-na_count) %>% #drop NA count
  column_to_rownames("Sample_Name") #set Sample_Name as rownames
Pocillopora_morphology$Morphotype <- Pocillopora_morpho_map[rownames(Pocillopora_morphology), "Morphotype"] #merge morphotype into morphology
morpho_lookup <- c(ac = "Acicular", #define label mapping
                   br = "Branching",
                   `br,da` = "Branching/Digitate",
                   da = "Digitate",
                   ef = "Effuse",
                   `ef,li` = "Effuse/Lineate",
                   fu = "Fused branches",
                   gr = "Granulate surface",
                   me = "Meandroid",
                   `me,gr` = "Meandroid/Granulate",
                   ve = "Verrucose",
                   `ve,da` = "Verrucose/Digitate",
                   `ve,ke` = "Verrucose/Keeled",
                   `ve,me` = "Verrucose/Meandroid",
                   `ve,me,da` = "Verrucose/Meandroid/Digitate",
                   velvety = "Velvety")
Pocillopora_morphology <- Pocillopora_morphology %>% #recode to full names
  mutate(Morphotype_full = recode(Morphotype, !!!morpho_lookup, .default = NA_character_)) %>% #recode morphotypes
  dplyr::select(-Morphotype) #drop original Morphotype column
simple_traits <- unique(unlist(str_split(na.omit(as.character(Pocillopora_morphology$Morphotype_full)), "/"))) #get individual morphotype components
for(tr in simple_traits){ #iterate over traits
  col <- paste0("Morphotype_", str_replace_all(tr, "\\s+", "_")) #define column name
  Pocillopora_morphology[[col]] <- ifelse(is.na(Pocillopora_morphology$Morphotype_full), NA_integer_, as.integer(str_detect(Pocillopora_morphology$Morphotype_full, tr))) #binary encode
}
Pocillopora_morphology$Morphotype_full <- NULL #drop full label column
Pocillopora_morphology <- remove.lowCV.multicollinearity.SOM(Pocillopora_morphology, #remove highly correlated and low variance variables
                                                             CV.threshold = 0.05, 
                                                             cor.threshold = 0.9)
ncol(Pocillopora_morphology) #number of traits: 22
nrow(Pocillopora_morphology) #number of samples: 175


## Train and cluster SOM
Pocillopora_SOM_data <- list(SNP = Pocillopora_SNP,
                             Microsatellites = Pocillopora_microsatellites,
                             morphology = Pocillopora_morphology,
                             Biogeography = Pocillopora_biogeography,
                             Symbiosis_haplotypes = Pocillopora_symbiosis_haplotypes)
Pocillopora_SOM <- train.SOM(input_data = Pocillopora_SOM_data, 
                             learning.rate.tuning = F,
                             max.NA.row = 0.4,
                             max.NA.col = 0.3,
                             parallel = F,
                             N_cores = 5,
                             N.steps = 5,
                             N.replicates = 5)
Pocillopora_SOM <- clustering.SOM(Pocillopora_SOM, 
                                  max.k = 20,
                                  clustering.method = "kmeans+BICjump+hierarch")


## Plot and evaluate results
Pocillopora_SOM$optim_k_summary
plot.Learning.SOM(Pocillopora_SOM)
plot.Layers.SOM(Pocillopora_SOM)
plot.Model.SOM(Pocillopora_SOM)
plot.Structure.SOM(Pocillopora_SOM, Individual.labels.font.size = 0.3)
plot.K.SOM(Pocillopora_SOM)

plot.Variable.Importance.SOM(SOM.output = Pocillopora_SOM, 
                             top.margin = 5,
                             left.margin = 6)
head(round(sort(Pocillopora_SOM$median_variable_importance[[1]], decreasing = T), 2), 300)
head(round(sort(Pocillopora_SOM$median_variable_importance[[2]], decreasing = T), 2), 20)
head(round(sort(Pocillopora_SOM$median_variable_importance[[3]], decreasing = T), 2), 40)
head(round(sort(Pocillopora_SOM$median_variable_importance[[4]], decreasing = T), 2), 10)
head(round(sort(Pocillopora_SOM$median_variable_importance[[5]], decreasing = T), 2), 20)

Pocillopora_SOM_ancestry_matrix <- as.data.frame(Pocillopora_SOM$ancestry_matrix)
Pocillopora_SOM_common_samples <- intersect(rownames(Pocillopora_SOM_ancestry_matrix), rownames(Pocillopora_species_names))
Pocillopora_SOM_ancestry_mat_sub <- Pocillopora_SOM_ancestry_matrix[Pocillopora_SOM_common_samples, , drop = FALSE]
Pocillopora_SOM_Pocillopora_SOM_species_names_sub <- Pocillopora_species_names[Pocillopora_SOM_common_samples, , drop = FALSE]
Pocillopora_SOM_ancestry_matrix <- cbind(Pocillopora_SOM_ancestry_mat_sub, Pocillopora_SOM_Pocillopora_SOM_species_names_sub) #compare ancestry matrix with species hypotheses
unique(sort(Pocillopora_SOM_ancestry_matrix$PSH))
length(unique(sort(Pocillopora_SOM_ancestry_matrix$PSH)))
table(Pocillopora_SOM_ancestry_matrix$PSH) #primary species hypotheses



################################################################################
#### Empirical example: Polygonia anglewing butterflies in Western Canada (Dupuis et al 2018)
################################################################################

## Import and process genetic SNP data
Polygonia_SNP <- process.SNP.data.SOM(vcf.path = "Test data/Dupuis et al 2018/Polygonia_961SNPs.vcf", #filter loci and individuals and create SNP matrix dataframe
                                      missing.loci.cutoff.lenient = 0.4,
                                      missing.loci.cutoff.final = 0.2,
                                      missing.individuals.cutoff = 0.15)
rownames(Polygonia_SNP) <- sub(".*?(\\d+)$", "\\1", rownames(Polygonia_SNP)) #only keep numeric identifier as rownames
ncol(Polygonia_SNP) #number of loci: 955
nrow(Polygonia_SNP) #number of samples: 22


## Import and filter COI data
Polygonia_COI <- process.SNP.data.SOM(nexus.path = "Test data/Dupuis et al 2018/Polygonia_COI.nex",
                                      missing.loci.cutoff.lenient = 0.4,
                                      missing.loci.cutoff.final = 0.2,
                                      missing.individuals.cutoff = 0.2)
Polygonia_COI_numeric_rownames <- sub(".*?(\\d+)$", "\\1", rownames(Polygonia_COI)) #extract numeric code from each rowname (e.g., "pf_8301" -> "8301")
Polygonia_COI <- Polygonia_COI[!duplicated(Polygonia_COI_numeric_rownames), , drop = FALSE] #keep only first occurrence for each numeric code (remove duplicates)
rownames(Polygonia_COI) <- Polygonia_COI_numeric_rownames[!duplicated(Polygonia_COI_numeric_rownames)] #set rownames to unique numeric codes
ncol(Polygonia_COI) #number of loci: 215
nrow(Polygonia_COI) #number of samples: 254


## Import and process RGB values
Polygonia_RGB <- read.delim("Test data/Dupuis et al 2018/Polygonia_RGB_characters.txt", stringsAsFactors = FALSE)
rownames(Polygonia_RGB) <- Polygonia_RGB$Species
Polygonia_RGB <- Polygonia_RGB %>% 
  dplyr::select(-Name, -Species) #remove columns


## Import and process wing character scores
Polygonia_wing_scores <- read.delim("Test data/Dupuis et al 2018/Polygonia_visually_scored.txt", stringsAsFactors = FALSE)
rownames(Polygonia_wing_scores) <- Polygonia_wing_scores$Name
Polygonia_wing_scores <- Polygonia_wing_scores %>% 
  dplyr::select(-Name, -Species) %>% #remove columns
  rename(Wing_character_1 = Ch1,
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
Wing_character_8_states <- model.matrix(~ Wing_character_8 - 1, data = Polygonia_wing_scores) #make wing character 8 binary (since it is not ordinal)
colnames(Wing_character_8_states) <- paste0("Wing_character_8_state_", 1:4)
Polygonia_wing_scores <- cbind(Polygonia_wing_scores, Wing_character_8_states)
Polygonia_wing_scores$Wing_character_8 <- NULL


## Import and process meta data with spatial data, species names and morphotype
Polygonia_metadata <- read.csv("Test data/Dupuis et al 2018/Polygonia_metadata.csv",header = T, sep = ";")
rownames(Polygonia_metadata) <- Polygonia_metadata$ID
nrow(Polygonia_metadata) #number of samples: 265

Polygonia_spatial <- Polygonia_metadata %>% 
  dplyr::select(Latitude, Longitude) #create dataframe with Lat and Long
Polygonia_spatial$Elevation <- NA #initialize elevation column with NA
Polygonia_spatial_sf <- st_as_sf(Polygonia_spatial[!is.na(Polygonia_spatial$Latitude) & !is.na(Polygonia_spatial$Longitude), ], coords = c("Longitude", "Latitude"), crs = 4326) #extract elevation
Polygonia_spatial$Elevation[!is.na(Polygonia_spatial$Latitude) & !is.na(Polygonia_spatial$Longitude)] <- 
  elevatr::get_elev_point(locations = Polygonia_spatial_sf, prj = st_crs(Polygonia_spatial_sf)$proj4string, src = "aws")$elevation

Polygonia_morphotype <- Polygonia_metadata %>% 
  dplyr::select(Morphotype)

Polygonia_metadata <- Polygonia_metadata %>% 
  dplyr::select(Species, ID)


## Merge morphometric dataframes and make morphotype binary
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
Polygonia_morphology <- remove.lowCV.multicollinearity.SOM(Polygonia_morphology, #remove highly correlated and low-variance variables
                                                           CV.threshold = 0.05,
                                                           cor.threshold = 0.9)
ncol(Polygonia_morphology) #number of traits: 33
nrow(Polygonia_morphology) #number of samples: 217


## Import and process environmental dataset (variables extracted and processed by separate R script based on coordinates)
Polygonia_environmental <- read.csv("Test data/Dupuis et al 2018/Polygonia_environmental.csv",
                                    row.names = 1, header = TRUE)
Polygonia_environmental <- Polygonia_environmental %>% 
  dplyr::select(-Latitude, -Longitude, -Elevation)
Polygonia_environmental_rownames <- rownames(Polygonia_environmental) #save rownames
Polygonia_environmental <- as.data.frame(lapply(Polygonia_environmental, as.numeric)) #ensure all columns are numeric
rownames(Polygonia_environmental) <- Polygonia_environmental_rownames #reassign saved row names
Polygonia_environmental <- remove.lowCV.multicollinearity.SOM(Polygonia_environmental, #remove highly correlated and low-variance variables
                                                              CV.threshold = 0.05,
                                                              cor.threshold = 0.9)
ncol(Polygonia_environmental) #number of variables: 98
nrow(Polygonia_environmental) #number of samples: 265


## Match datasets and remove all NA rows
for (Polygonia_shared_data in c("Polygonia_morphology", "Polygonia_SNP", "Polygonia_COI",
                                "Polygonia_spatial", "Polygonia_environmental", "Polygonia_metadata")) {
  Polygonia_shared_data_mat <- get(Polygonia_shared_data)
  rownames(Polygonia_shared_data_mat) <- make.unique(as.character(rownames(Polygonia_shared_data_mat)))
  assign(Polygonia_shared_data, Polygonia_shared_data_mat, envir = .GlobalEnv)}
Polygonia_common_IDs <- Reduce(intersect, list(
  rownames(Polygonia_morphology),
  rownames(Polygonia_SNP),
  rownames(Polygonia_COI),
  rownames(Polygonia_spatial),
  rownames(Polygonia_environmental),
  rownames(Polygonia_metadata)))
Polygonia_morphology2 <- Polygonia_morphology[Polygonia_common_IDs, , drop = FALSE]
Polygonia_SNP2 <- Polygonia_SNP[Polygonia_common_IDs, , drop = FALSE]
Polygonia_COI2 <- Polygonia_COI[Polygonia_common_IDs, , drop = FALSE]
Polygonia_spatial2 <- Polygonia_spatial[Polygonia_common_IDs, , drop = FALSE]
Polygonia_environmental2 <- Polygonia_environmental[Polygonia_common_IDs, , drop = FALSE]
Polygonia_metadata2 <- Polygonia_metadata[Polygonia_common_IDs, , drop = FALSE]

Polygonia_morphology <- Polygonia_morphology2[rowSums(!is.na(Polygonia_morphology2)) > 0, , drop = FALSE]
Polygonia_SNP <- Polygonia_SNP2[rowSums(!is.na(Polygonia_SNP2)) > 0, , drop = FALSE]
Polygonia_COI <- Polygonia_COI2[rowSums(!is.na(Polygonia_COI2)) > 0, , drop = FALSE]
Polygonia_spatial <- Polygonia_spatial2[rowSums(!is.na(Polygonia_spatial2)) > 0, , drop = FALSE]
Polygonia_environmental <- Polygonia_environmental2[rowSums(!is.na(Polygonia_environmental2)) > 0, , drop = FALSE]
Polygonia_metadata <- Polygonia_metadata2[rowSums(!is.na(Polygonia_metadata2)) > 0, , drop = FALSE]
nrow(Polygonia_metadata) #number of shared samples: 179


## Update rownames by adding species names to ID
Polygonia_species_vec <- Polygonia_metadata$Species[match(rownames(Polygonia_morphology), Polygonia_metadata$ID)]
Polygonia_new_rownames <- paste(rownames(Polygonia_metadata), Polygonia_species_vec, sep = "_")
rownames(Polygonia_morphology) <- Polygonia_new_rownames
rownames(Polygonia_SNP) <- Polygonia_new_rownames
rownames(Polygonia_COI) <- Polygonia_new_rownames
rownames(Polygonia_spatial) <- Polygonia_new_rownames
rownames(Polygonia_metadata) <- Polygonia_new_rownames
rownames(Polygonia_environmental) <- Polygonia_new_rownames


## Train and cluster SOM
Polygonia_all_data <- list(Morphology = Polygonia_morphology,
                           SNP = Polygonia_SNP,
                           COI = Polygonia_COI,
                           Environmental = Polygonia_environmental,
                           Spatial = Polygonia_spatial)
Polygonia_SOM <- train.SOM(input_data = Polygonia_all_data, 
                           max.NA.row = 0.15,
                           max.NA.col = 0.4,
                           parallel = F,
                           N_cores = 5,
                           learning.rate.tuning = F,
                           save.SOM.results = T,
                           N.steps = 3, 
                           N.replicates = 3)
Polygonia_SOM <- clustering.SOM(SOM.output = Polygonia_SOM,
                                max.k = 10,
                                clustering.method = "kmeans+BICjump+hierarch")


## Evaluate and plot results
Polygonia_SOM$optim_k_summary
plot.K.SOM(Polygonia_SOM)
plot.Layers.SOM(Polygonia_SOM)
plot.Model.SOM(Polygonia_SOM)
plot.Map.SOM(Polygonia_SOM, 
             Coordinates = Polygonia_spatial[, c(1:2)], 
             lat.buffer.range = 1, 
             lon.buffer.range = 2,
             north.arrow.position = c(0.01, 0.92), #position (x, y) of north arrow relative to map
             north.arrow.length = 0.7, #length of north arrow
             north.arrow.N.position = 0.3, #position of north arrow "N"
             north.arrow.N.size = 1, #size of north arrow "N"
             scale.position = c(0.78, 0.02)) #relative position (x, y) of scale
plot.Structure.SOM(Polygonia_SOM, Individual.labels.font.size = 0.25)
plot.Learning.SOM(Polygonia_SOM)
plot.Variable.Importance.SOM(Polygonia_SOM, left.margin = 8, bars.threshold.N = 50)
head(Polygonia_SOM$ancestry_matrix)
head(Polygonia_metadata)

head(round(sort(Polygonia_SOM$median_variable_importance[[1]], decreasing = T), 2), 10)
head(round(sort(Polygonia_SOM$median_variable_importance[[2]], decreasing = T), 2), 600)
head(round(sort(Polygonia_SOM$median_variable_importance[[3]], decreasing = T), 2), 30)
head(round(sort(Polygonia_SOM$median_variable_importance[[4]], decreasing = T), 2), 10)
head(round(sort(Polygonia_SOM$median_variable_importance[[5]], decreasing = T), 2), 10)

Polygonia_ancestry <- as.data.frame(Polygonia_SOM$ancestry_matrix)
Polygonia_ancestry$ID <- rownames(Polygonia_ancestry)
Polygonia_metadata$ID <- as.character(rownames(Polygonia_metadata))
Polygonia_ancestry <- merge(Polygonia_metadata, Polygonia_ancestry, by = "ID") #examine ancestry matrix and species labels
table(Polygonia_ancestry$Species)



  ################################################################################
  #### Viburnum shrubs dataset - Eastern North America (Spriggs et al 2018)
  ################################################################################
  
  ## Import and process genetic SNP data
  Viburnum_SNP <- process.SNP.data.SOM(vcf.path = "Test data/Spriggs et al 2018/nudum-c88-d6-min50.vcf",
                                       missing.loci.cutoff.lenient = 0.4,
                                       missing.loci.cutoff.final = 0.2, 
                                       missing.individuals.cutoff = 0.5)
  ncol(Viburnum_SNP) #number of SNPs: 2511
  nrow(Viburnum_SNP) #number of samples: 65
  
  
  
  ## Import and process morphological dataset
  Viburnum_morphology <- read.delim("Test data/Spriggs et al 2018/morphological_trait_data2.txt", stringsAsFactors = FALSE)
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
  Viburnum_morphology <- remove.lowCV.multicollinearity.SOM(Viburnum_morphology, #remove highly correlated and low-variance variables
                                                            CV.threshold = 0.05,
                                                            cor.threshold = 0.9)
  ncol(Viburnum_morphology) #number of traits: 9
  nrow(Viburnum_morphology) #number of samples: 145
  
  
  ## Import and process metadata
  Viburnum_metadata <- read.delim("Test data/Spriggs et al 2018/morphological_trait_data2.txt", stringsAsFactors = FALSE)
  Viburnum_metadata <- Viburnum_metadata[!duplicated(Viburnum_metadata$Individual), ] #remove duplicate IDs
  rownames(Viburnum_metadata) <- Viburnum_metadata$Individual #add rownames
  Viburnum_metadata <- Viburnum_metadata %>% #only keep State and County columns
    dplyr::select(State, County)
  nrow(Viburnum_metadata) #number of samples: 145
  
  
  ## Update rownames and remove outgroup species
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
  
  
  ## Match datasets and remove rows with all NA
  Viburnum_shared_ids <- Reduce(intersect, list(rownames(Viburnum_SNP), rownames(Viburnum_morphology), rownames(Viburnum_metadata)))
  Viburnum_SNP <- Viburnum_SNP[Viburnum_shared_ids, , drop = FALSE]
  Viburnum_morphology <- Viburnum_morphology[Viburnum_shared_ids, , drop = FALSE]
  Viburnum_metadata <- Viburnum_metadata[Viburnum_shared_ids, , drop = FALSE]
  nrow(Viburnum_metadata) #shared number of samples: 52
  
  
  ## Train and cluster SOM
  Viburnum_SOM_data <- list(Morphology = Viburnum_morphology, SNP = Viburnum_SNP)
  Viburnum_SOM <- train.SOM(Viburnum_SOM_data,
                            max.NA.row = 0.4,
                            max.NA.col = 0.33,
                            N.steps = 15, 
                            N.replicates = 10,
                            learning.rate.tuning = F, 
                            parallel = T,
                            save.SOM.results = T,
                            N_cores = 3)
  Viburnum_SOM <- clustering.SOM(Viburnum_SOM, 
                                 clustering.method = "kmeans+BICjump+hierarch", 
                                 max.k = 8)
  
  
  ## Evaluate and plot results
  Viburnum_SOM$optim_k_summary
  plot.Learning.SOM(Viburnum_SOM)
  plot.Layers.SOM(Viburnum_SOM)
  plot.K.SOM(Viburnum_SOM)
  plot.Model.SOM(Viburnum_SOM)
  plot.Variable.Importance.SOM(Viburnum_SOM, 
                               left.margin = 4,
                               bar.label.font.size = 0.4)
  plot.Structure.SOM(Viburnum_SOM, sort.by.col = 2, Individual.labels.font.size = 0.8)
  
  Viburnum_SOM_species_state_vector <- c(ELS002 = "cassinoides, NY", #create named character vector with species and state info based on Figure 4 in Spriggs et al. 2018
                                         ELS003 = "nitidum, FL",
                                         ELS050 = "nudum, NC",
                                         ELS054 = "nitidum, NC",
                                         ELS079 = "cassinoides, NH",
                                         ELS081 = "cassinoides, CT",
                                         ELS198 = "nitidum, FL",
                                         ELS251 = "nitidum, FL",
                                         ELS252 = "nitidum, FL",
                                         ELS265 = "nudum, SC",
                                         ELS309 = "cassinoides, NC",
                                         ELS322 = "nitidum, NC",
                                         ELS324 = "nitidum, NC",
                                         ELS328 = "nudum, NC",
                                         ELS340 = "cassinoides, NC",
                                         ELS384 = "cassinoides, ME",
                                         ELS410 = "cassinoides, NH",
                                         ELS413 = "cassinoides, MA",
                                         ELS423 = "cassinoides, MA",
                                         ELS537 = "cassinoides, MA",
                                         ELS542 = "cassinoides, MA",
                                         ELS564 = "cassinoides, NC",
                                         ELS586 = "nudum, AL",
                                         ELS607 = "nudum, SC",
                                         ELS612 = "nudum, SC",
                                         ELS619 = "nudum, NC",
                                         ELS621 = "nudum, NC",
                                         ELS630 = "nudum, VA",
                                         ELS631 = "nudum, VA",
                                         ELS645 = "nudum, DE",
                                         ELS659 = "nudum, GA",
                                         ELS660 = "nudum, GA",
                                         ELS661 = "nitidum, GA",
                                         ELS667 = "nudum, GA",
                                         ELS670 = "nudum, GA",
                                         ELS671 = "nudum, GA",
                                         ELS674 = "nudum, GA",
                                         ELS675 = "nudum, GA",
                                         ELS677 = "nitidum, GA",
                                         ELS678 = "nitidum, GA",
                                         ELS684 = "nudum, GA",
                                         ELS685 = "nudum, GA",
                                         ELS688 = "nitidum, GA")
  Viburnum_ancestry <- as.data.frame(Viburnum_SOM$ancestry_matrix)
  Viburnum_ancestry$Species_State <- Viburnum_SOM_species_state_vector[rownames(Viburnum_ancestry)]
  Viburnum_ancestry <- Viburnum_ancestry %>%
    tidyr::separate(Species_State, into = c("Species", "State"), sep = ", ", remove = TRUE)
  table(Viburnum_ancestry$Species)



################################################################################
#### Empirical example: Microcebus lemurs dataset - Madagaskar (van Elst et al 2024)
################################################################################

## Import and process genetic SNP data
Microcebus_SNP <- process.SNP.data.SOM(
  vcf.path = "Test data/van Elst et al 2024/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss0.05.thinned.vcf", #VCF file path
  missing.loci.cutoff.lenient = 0.4, #remove loci with >40% missing data (lenient)
  missing.loci.cutoff.final = 0.05, #remove loci with >5% missing data (stricter)
  missing.individuals.cutoff = 0.5) #remove individuals with >50% missing data
Microcebus_species_split <- stringr::str_split_fixed(rownames(Microcebus_SNP), "_", n = 2) #split rownames by underscore
rownames(Microcebus_SNP) <- Microcebus_species_split[, 2] #set rownames to just the ID part
ncol(Microcebus_SNP) #number of SNPs: 9385
nrow(Microcebus_SNP) #number of samples: 213


## Import and process multiple data dataset 2 containing range of data types
Microcebus_multiple_data2 <- read.csv("Test data/van Elst et al 2024/01_Microcebus_morphological_data.csv", 
                                      stringsAsFactors = FALSE, header = T, sep = ";")
Microcebus_multiple_data2 <- Microcebus_multiple_data2 %>% #only keep individuals that are Rad sequenced (have SNP data)
  dplyr::filter(RADSeq.available != "no" & !is.na(RADSeq.available))
rownames(Microcebus_multiple_data2) <- Microcebus_multiple_data2$Individual.ID


## Import and process multiple data dataset containing range of data types
Microcebus_multiple_data <- read.csv("Test data/van Elst et al 2024/data.csv", 
                                     stringsAsFactors = FALSE, header = T, sep = ";")
Microcebus_multiple_data <- Microcebus_multiple_data[!duplicated(Microcebus_multiple_data$Individual.ID), ] #remove duplicate IDs
Microcebus_multiple_data <- Microcebus_multiple_data[!is.na(Microcebus_multiple_data$Individual.ID) & Microcebus_multiple_data$Individual.ID != "", ] # drop rows where Individual.ID is NA or empty-string
rownames(Microcebus_multiple_data) <- Microcebus_multiple_data$Individual.ID
Microcebus_multiple_data <- Microcebus_multiple_data[rownames(Microcebus_multiple_data) %in% rownames(Microcebus_multiple_data2), , drop = FALSE] #only keep rownames that are present in morphological and SNP data


## Update rownames
extract.sample.id.and.species <- function(sample_name) { #create extractor function
  if (grepl("^carlesii_", sample_name)) {
    species_name <- "carlesii"
    sample_id <- sub("^carlesii_", "", sample_name)
  } else if (grepl("^rhytidophyllum_", sample_name)) {
    species_name <- "rhytidophyllum"
    sample_id <- sub("^rhytidophyllum_", "", sample_name)
  } else {
    parts <- strsplit(sample_name, "_")[[1]]
    if (length(parts) > 2 &&
        grepl("^[a-z]+$", parts[length(parts) - 1]) &&
        grepl("^[A-Z]+$", parts[length(parts)])) {
      species_name <- paste(parts[(length(parts) - 1):length(parts)], collapse = "_")
      sample_id <- paste(parts[1:(length(parts) - 2)], collapse = "_")
    } else {
      species_name <- parts[length(parts)]
      sample_id <- paste(parts[1:(length(parts) - 1)], collapse = "_")}}
  data.frame(Sample_ID = sample_id, Species_Name = species_name, stringsAsFactors = FALSE)}
Microcebus_id_species_df <- do.call(rbind, lapply(rownames(Microcebus_SNP), extract.sample.id.and.species)) #update rownames using function

Microcebus_duplicates <- duplicated(Microcebus_id_species_df$Sample_ID) #extract duplicates
if (any(Microcebus_duplicates)) {message("Dropping duplicates: ", paste(unique(Microcebus_id_species_df$Sample_ID[Microcebus_duplicates]), collapse = ", "))}
Microcebus_id_species_df <- Microcebus_id_species_df[!Microcebus_duplicates, , drop = FALSE] #drop duplicates
Microcebus_SNP <- Microcebus_SNP[!Microcebus_duplicates, , drop = FALSE] #drop duplicates
rownames(Microcebus_SNP) <- Microcebus_id_species_df$Sample_ID #assign rownames


## Update some rownames manually due to not-fully matching rownames
Microcebus_duplicates_original_ids <- rownames(Microcebus_SNP) #extract current rownames
Microcebus_duplicates_new_ids <- str_replace(Microcebus_duplicates_original_ids, "^(\\d+)y(\\d+)_(.+)$", "\\1-\\2\\3") #modify rownames
Microcebus_duplicates_new_ids <- str_replace(Microcebus_duplicates_new_ids, "^(.*?)y(\\d+)$", "\\1-\\2") #modify rownames
rownames(Microcebus_SNP) <- Microcebus_duplicates_new_ids #update rownames

Microcebus_id_map <- c( #define mapping of SNP IDs (left) and morphology IDs (right)
  "00016A8577"   = "00-016A-8577",
  "00016A875E"   = "00-016A-875E",
  "01-06hely"    = "M01-06hely",
  "01-13bom"     = "F01-13Bom",
  "02-08befa"    = "M02-08Befa",
  "04-06hely"    = "W04-06hely",
  "04-13hely"    = "F04-13hely",
  "05-02oka"     = "F05-02oka",
  "06-02mahi"    = "F06-02Mah",
  "06-05mamy"    = "F06-05mamy",
  "06-08hely"    = "M06-08hely",
  "07-02mahi"    = "F07-02Mah",
  "07-06habe"    = "M07-06habe",
  "07-08hely"    = "M07-08hely",
  "08-08hely"    = "M08-08hely",
  "10-07hely"    = "W10-07hely",
  "10-13jbb"     = "M10-13JBB",
  "11-02oka"     = "F11-02oka",
  "12-02oka"     = "M12-02oka",
  "15-02mahi"    = "F15-02Mah",
  "3-11hara"     = "F03-11hara",
  "55-04bibo"    = "F55-04bibo",
  "71-04tsin"    = "M71-04Tsin",
  "83-04injo"    = "M83-04injo",
  "90-04anji"    = "F90-04Anji",
  "A100"         = "A100_2014",
  "A3"           = "A032014",
  "A8"           = "A082014",
  "A12"          = "A-MM-A12",
  "A17"          = "A17_2014",
  "A23"          = "A-MM-A23",
  "A24"          = "A-MM-A24",
  "A45"          = "A45_2014",
  "A61"          = "A61_2014",
  "AMBH_D7"      = "AMBH_D07_2013",
  "ANALF_B22"    = "ANALF_B22_2013",
  "ANALF_B28"    = "ANALF_B28_2013",
  "ANALM_C56"    = "ANALM_C56_2012",
  "ANALV_A19"    = "ANALV_A19_2013",
  "ANKA_F137"    = "ANKA_F137_2012",
  "ANTSB_H53"    = "ANTSB_H53_2011",
  "ANTSR_L16"    = "ANTSR_L16_2011",
  "B13"          = "B-MM-B13",
  "B14"          = "B14_2014",
  "B24"          = "B-MM-B24",
  "B34"          = "B-MM-B34",
  "B77"          = "B77_2014",
  "BC3"          = "B-MF-BC3",
  "BC49"         = "BC49_2017",
  "BD1"          = "B-MF-BD1",
  "BE51"         = "BE51_2017",
  "BEN_G62"      = "BEN_G62_2011",
  "BEZAV_C30"    = "BEZAV_C30_2013",
  "BIN_C32"      = "BIN_C32_2010",
  "C23"          = "C23_2013",
  "C24"          = "C24_2013",
  "D72"          = "D72_2015",
  "E23"          = "E23_2014",
  "F06-02Mah"    = "F06-02Mah",
  "F07-02Mah"    = "F07-02Mah",
  "F15-02Mah"    = "F15-02Mah",
  "F28"          = "F28_2015",
  "F60"          = "F60_2015",
  "I01"          = "I01_2015",
  "JMR001"       = "JMR-001",
  "JMR002"       = "JMR-002",
  "M07-06habe"   = "M07-06habe",
  "M10-13JBB"    = "M10-13JBB",
  "M12-02oka"    = "M12-02oka",
  "M12-16mibe"   = "M12-16mibe",
  "M16-16Lok"    = "M16-16Lok",
  "M22-16Anji"   = "M22-16Anji",
  "MBB 019"      = "MBB 019",
  "MBB 020"      = "MBB 020",
  "SALAF_B66"    = "SALAF_B66_2013",
  "SALAF_B77"    = "SALAF_B77_2013",
  "SOL_K18"      = "SOL_K18_2011",
  "f08-13"       = "F08-13zana",
  "f09-13"       = "F09-13zana",
  "f13-15"       = "F13-15sely",
  "f142-17"      = "M142-17jabe",
  "f17-13"       = "F17-13zana")
Microcebus_original_ids_2 <- rownames(Microcebus_SNP) #apply mapping to rownames
Microcebus_new_ids_2 <- ifelse(Microcebus_original_ids_2 %in% names(Microcebus_id_map), Microcebus_id_map[Microcebus_original_ids_2], Microcebus_original_ids_2)
rownames(Microcebus_SNP) <- Microcebus_new_ids_2 #update rownames


## Match rownames and remove rows with only NA
Microcebus_multiple_data_all_rows <- union(rownames(Microcebus_multiple_data), rownames(Microcebus_multiple_data2))
Microcebus_multiple_data_combined <- data.frame(row.names = Microcebus_multiple_data_all_rows) #make empty “combined” with those rows
Microcebus_multiple_data_combined[rownames(Microcebus_multiple_data), names(Microcebus_multiple_data)] <- Microcebus_multiple_data
Microcebus_multiple_data_new_cols <- setdiff(names(Microcebus_multiple_data2), names(Microcebus_multiple_data))
Microcebus_multiple_data_combined[rownames(Microcebus_multiple_data2), Microcebus_multiple_data_new_cols] <- Microcebus_multiple_data2[Microcebus_multiple_data_new_cols]
Microcebus_multiple_data_combined <- Microcebus_multiple_data_combined[rownames(Microcebus_SNP), , drop = FALSE] #re‑order and pad to SNP samples
Microcebus_multiple_data_combined <- Microcebus_multiple_data_combined %>% filter(if_any(everything(), ~ !is.na(.))) #remove rows with all NA
Microcebus_SNP <- Microcebus_SNP[rownames(Microcebus_multiple_data_combined), , drop = FALSE] #subset SNP dataset to combined samples


## Add Reproductive_period_month column (1–12) for reproductively active animals (based on Supplementary Methods information)
Microcebus_multiple_data_combined <- Microcebus_multiple_data_combined %>%
  mutate(MonthNum = as.integer(str_extract(month, "^\\d{2}")), #extract numeric month (01–12)
         Reproductive_period_month = case_when(
           male_repro == "enlarged" ~ MonthNum, #male testes enlarged indicates reproductive activity
           female_repro %in% c("swollen", "open") ~ MonthNum, #swollen/open female indicates reproductive activity
           female_repro == "pregnant" ~ ((MonthNum - 2 - 1) %% 12) + 1, #pregnant female indicates reproductive activity 2 months earlier
           female_repro == "lactating" ~ ((MonthNum - 3 - 1) %% 12) + 1, #lactating female indicates reproductive activity 3 months earlier
           TRUE ~ NA_integer_)) %>% #else NA (not in reproductive period)
  dplyr::select(-MonthNum) #drop helper column


## Create morphology dataset (including reproductive traits)
Microcebus_morphology <- Microcebus_multiple_data_combined %>% 
  dplyr::select(ear.length, ear.width, head.length, head.width, interorbital.dist, 
                intraorbital.dist, snout.length, lower.leg.length, hind.foot.length,
                third.toe.length, body.length, tail.length, body.mass, tail.circumference,
                testis.width.total, testis.width.left, testis.width.right, 
                testis.length.left, testis.length.right, Reproductive_period_month) %>%
  rename(Interorbital_Distance = interorbital.dist, #modify column names
         Intraorbital_Distance = intraorbital.dist) %>%
  rename_with(~ str_replace_all(.x, "\\.", "_"), everything()) %>% #rename dots to underscores
  rename_with(~ str_to_title(.x), everything()) %>% #capitalize first letter
  mutate(across(where(is.character), ~ na_if(., "#NV"))) %>% #replace “#NV” with NA
  mutate(across(where(is.character), ~ sub(",", ".", .))) %>% #change comma to decimal point
  mutate(across(where(is.character), as.numeric)) #convert remaining chars to numeric
Microcebus_morphology$Testis_width_total <- as.character(Microcebus_morphology$Testis_width_total)
Microcebus_morphology_Testis_replacement <- c("0", "large", "small") #replace 0, "large" and "small" with NA
Microcebus_morphology$Testis_width_total[Microcebus_morphology$Testis_width_total %in% Microcebus_morphology_Testis_replacement] <- NA
Microcebus_morphology$Testis_width_total <- as.numeric(Microcebus_morphology$Testis_width_total)
Microcebus_morphology <- remove.lowCV.multicollinearity.SOM(Microcebus_morphology, #remove highly correlated and low-variance variables
                                                            CV.threshold = 0.05,
                                                            cor.threshold = 0.9)
ncol(Microcebus_morphology) #number of traits: 17
nrow(Microcebus_morphology) #number of samples: 72


## Create metadata dataset with species assignments
Microcebus_metadata <- Microcebus_multiple_data_combined %>% 
  dplyr::select(Individual_number, sex, clade, population, fine_sub_pops) %>%
  rename(Sex = sex, Clade = clade, Population = population, Subpopulations = fine_sub_pops)
Microcebus_species_map <- list( #define species and revised species as list based on supplementary table S17
  "00-016A-8577"    = c("M. ganzhorni",	"M. murinus"),
  "00-016A-875E"    = c("M. ganzhorni",	"M. murinus"),
  "A100_2014"       = c("M. sp. 1", "M. arnholdi"),
  "A17_2014"        = c("M. sp. 1", "M. arnholdi"),
  "A45_2014"        = c("M. sp. 1", "M. arnholdi"),
  "A61_2014"        = c("M. sp. 1", "M. arnholdi"),
  "A-MM-A12"        = c("M. jonahi", "M. jonahi"),
  "A-MM-A23"        = c("M. jonahi", "M. jonahi"),
  "A-MM-A24"        = c("M. jonahi", "M. jonahi"),
  "AMBH_D07_2013"   = c("M. sp. 1", "M. arnholdi"),
  "ANALF_B22_2013"  = c("M. tavaratra", "M. tavaratra"),
  "ANALF_B28_2013"  = c("M. tavaratra", "M. tavaratra"),
  "ANALM_C56_2012"  = c("M. tavaratra", "M. tavaratra"),
  "ANALV_A19_2013"  = c("M. sp. 1", "M. arnholdi"),
  "ANKA_F137_2012"  = c("M. tavaratra", "M. tavaratra"),
  "ANTSB_H53_2011"  = c("M. tavaratra", "M. tavaratra"),
  "ANTSR_L16_2011"  = c("M. tavaratra", "M. tavaratra"),
  "B14_2014"        = c("M. sp. 1", "M. arnholdi"),
  "B77_2014"        = c("M. sp. 1", "M. arnholdi"),
  "B-MM-B13"        = c("M. jonahi", "M. jonahi"),
  "B-MM-B24"        = c("M. jonahi", "M. jonahi"),
  "B-MM-B34"        = c("M. jonahi", "M. jonahi"),
  "B-MF-BC3"        = c("M. lehilahytsara", "M. lehilahytsara"),
  "B-MF-BD1"        = c("M. jonahi", "M. jonahi"),
  "BC49_2017"       = c("M. arnholdi", "M. arnholdi"),
  "BE51_2017"       = c("M. arnholdi", "M. arnholdi"),
  "BEN_G62_2011"    = c("M. tavaratra", "M. tavaratra"),
  "BEZAV_C30_2013"  = c("M. sp. 1", "M. arnholdi"),
  "BIN_C32_2010"    = c("M. tavaratra", "M. tavaratra"),
  "C23_2013"        = c("M. tavaratra", "M. tavaratra"),
  "C24_2013"        = c("M. lehilahytsara", "M. lehilahytsara"),
  "D72_2015"        = c("M. tavaratra", "M. tavaratra"),
  "E23_2014"        = c("M. sp. 1", "M. arnholdi"),
  "F01-13Bom"       = c("M. myoxinus", "M. myoxinus"),
  "F03-11hara"      = c("M. bongolavensis", "M. ravelobensis"),
  "F04-13hely"      = c("M. macarthurii", "M. macarthurii"),
  "F05-02oka"       = c("M. mamiratra", "M. mamiratra"),
  "F06-02Mah"       = c("M. sambiranensis", "M. sambiranensis"),
  "F06-05mamy"      = c("M. sambiranensis", "M. sambiranensis"),
  "F07-02Mah"       = c("M. sambiranensis", "M. sambiranensis"),
  "F08-13zana"      = c("M. bongolavensis", "M. ravelobensis"),
  "F09-13zana"      = c("M. bongolavensis", "M. ravelobensis"),
  "F13-15sely"      = c("M. bongolavensis", "M. ravelobensis"),
  "F15-02Mah"       = c("M. sambiranensis", "M. sambiranensis"),
  "F17-13zana"      = c("M. bongolavensis", "M. ravelobensis"),
  "F28_2015"        = c("M. sp. 1", "M. arnholdi"),
  "F55-04bibo"      = c("M. bongolavensis", "M. ravelobensis"),
  "F60_2015"        = c("M. sp. 1", "M. arnholdi"),
  "F90-04Anji"      = c("M. danfossi", "M. danfossi"),
  "I01_2015"        = c("M. mamiratra", "M. mamiratra"),
  "JMR-001"         = c("M. lehilahytsara", "M. lehilahytsara"),
  "JMR-002"         = c("M. lehilahytsara", "M. lehilahytsara"),
  "M01-06hely"      = c("M. macarthurii", "M. macarthurii"),
  "M02-08Befa"      = c("M. danfossi", "M. danfossi"),
  "M06-08hely"      = c("M. macarthurii", "M. macarthurii"),
  "M07-06habe"      = c("M. mittermeieri", "M. lehilahytsara"),
  "M07-08hely"      = c("M. macarthurii", "M. macarthurii"),
  "M08-08hely"      = c("M. macarthurii", "M. macarthurii"),
  "M10-13JBB"       = c("M. ravelobensis", "M. ravelobensis"),
  "M12-02oka"       = c("M. mamiratra", "M. mamiratra"),
  "M142-17jabe"     = c("M. murinus (north)", "M. murinus"),
  "M71-04Tsin"      = c("M. murinus (north)", "M. murinus"),
  "M83-04injo"      = c("M. danfossi", "M. danfossi"),
  "RMR115"          = c("M. boraha", "M. simmonsi"),
  "RMR116"          = c("M. boraha", "M. simmonsi"),
  "RMR124"          = c("M. boraha", "M. simmonsi"),
  "RMR129"          = c("M. boraha", "M. simmonsi"),
  "RMR131"          = c("M. marohita", "M. jollyae"),
  "SALAF_B66_2013"  = c("M. sp. 1", "M. arnholdi"),
  "SALAF_B77_2013"  = c("M. sp. 1", "M. arnholdi"),
  "SOL_K18_2011"    = c("M. tavaratra", "M. tavaratra"),
  "W04-06hely"      = c("M. macarthurii", "M. macarthurii"),
  "W10-07hely"      = c("M. mittermeieri", "M. lehilahytsara"))
Microcebus_species_mat <- t(sapply(rownames(Microcebus_metadata), function(id) {if (!is.null(Microcebus_species_map[[id]])) {Microcebus_species_map[[id]]} else {c(NA, NA)}})) #extract species and species_revised using mapping list
Microcebus_metadata$Species <- Microcebus_species_mat[, 1] #add Species
Microcebus_metadata$Species_revised <- Microcebus_species_mat[, 2] #add Species revised
length(unique(Microcebus_metadata$Species))
length(unique(Microcebus_metadata$Species_revised))
nrow(Microcebus_metadata) #number of samples: 72


## Import and process environmental dataset (variables extracted and processed by separate R script based on coordinates)
Microcebus_environmental <- read.csv("Test data/van Elst et al 2024/Microcebus_environmental.csv", row.names = 1, stringsAsFactors = FALSE) %>%
  rownames_to_column("ID")
Microcebus_Bio_variables <- Microcebus_multiple_data_combined %>% #extract Bio variables
  rownames_to_column("ID") %>%
  dplyr::select(ID, bio01:bio19)
Microcebus_environmental <- Microcebus_environmental %>% #add Bio variables to environmental dataset
  left_join(Microcebus_Bio_variables, by = "ID") %>%
  column_to_rownames("ID")
Microcebus_environmental_rownames <- rownames(Microcebus_environmental) #save rownames for later
Microcebus_environmental <- Microcebus_environmental %>% 
  select(-Latitude, -Longitude, -Elevation)
Microcebus_environmental <- as.data.frame(lapply(Microcebus_environmental, as.numeric)) #ensure all columns are numeric
rownames(Microcebus_environmental) <- Microcebus_environmental_rownames #keep rownames
Microcebus_environmental <- remove.lowCV.multicollinearity.SOM(Microcebus_environmental, #remove highly correlated and low-variance variables
                                                               CV.threshold = 0.05,
                                                               cor.threshold = 0.9)
ncol(Microcebus_environmental) #number of variables: 34
nrow(Microcebus_environmental) #number of samples: 75


## Create spatial dataset with Latitude, Longitude and Elevation
Microcebus_spatial <- Microcebus_multiple_data_combined %>% 
  dplyr::select(latitude, longitude) %>% #add Latitude and Longitude
  rename(Latitude = latitude, Longitude = longitude)
Microcebus_environmental_spatial <- read.csv("Test data/van Elst et al 2024/Microcebus_environmental.csv", row.names = 1, stringsAsFactors = FALSE)
Microcebus_environmental_spatial <- Microcebus_environmental_spatial %>% dplyr::select(Elevation)
Microcebus_spatial$Elevation <- Microcebus_environmental_spatial[rownames(Microcebus_spatial), "Elevation"]
Microcebus_spatial$Elevation <- Microcebus_environmental_spatial[rownames(Microcebus_spatial), "Elevation"]
nrow(Microcebus_spatial) #number of samples: 72


#### Match all datasets to rownames in morphology dataset
shared_ids <- rownames(Microcebus_morphology) #define target set
Microcebus_SNP <- Microcebus_SNP[shared_ids, , drop = FALSE] #match to morphology
Microcebus_environmental <- Microcebus_environmental[shared_ids, , drop = FALSE]
Microcebus_spatial <- Microcebus_spatial[shared_ids, , drop = FALSE]
Microcebus_metadata <- Microcebus_metadata[shared_ids, , drop = FALSE]


## Train and cluster SOM on full data
Microcebus_SOM_full_data <- list(SNP = Microcebus_SNP,
                                 Morphology = Microcebus_morphology,
                                 Environmental = Microcebus_environmental,
                                 Spatial = Microcebus_spatial)
Microcebus_SOM_full <- train.SOM(Microcebus_SOM_full_data, 
                                 max.NA.row = 0.4,
                                 max.NA.col = 0.4,
                                 N.steps = 10, 
                                 parallel = F, 
                                 N_cores = 5,
                                 save.SOM.results = T,
                                 overwrite.SOM.results = T,
                                 learning.rate.tuning = F,
                                 N.replicates = 5)
Microcebus_SOM_full <- clustering.SOM(Microcebus_SOM_full, 
                                      max.k = 10, 
                                      clustering.method = "kmeans+BICjump+hierarch")


## Evaluate and plot results for full data
Microcebus_SOM_full$optim_k_summary
plot.Learning.SOM(Microcebus_SOM_full)
plot.Layers.SOM(Microcebus_SOM_full)
plot.K.SOM(Microcebus_SOM_full)
plot.Model.SOM(Microcebus_SOM_full)
plot.Variable.Importance.SOM(Microcebus_SOM_full,
                             left.margin = 5.8,
                             bar.label.font.size = 0.4, 
                             importance.threshold = 0.1)
plot.Structure.SOM(Microcebus_SOM_full, sort.by.col = 1)
plot.Map.SOM(SOM.output = Microcebus_SOM_full, 
             Coordinates = Microcebus_spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 1, #add coordinates as buffer range around latitude coordinates
             lon.buffer.range = 1, #add coordinates as buffer range around longitude coordinates
             pie.size = 2.5, #pie chart size
             north.arrow.position = c(0.03, 0.86), #position (x, y) of north arrow relative to map
             north.arrow.length = 0.5, #length of north arrow
             north.arrow.N.position = 0.2, #position of north arrow "N"
             north.arrow.N.size = 1) #size of north arrow "N"


## Compare ancestry coefficients with prior species and proposed ("revised") species labels for full data
Microcebus_ancestry <- as.data.frame(Microcebus_SOM_full$ancestry_matrix)
Microcebus_ancestry$Species <- Microcebus_metadata$Species[match(rownames(Microcebus_SOM_full$ancestry_matrix), rownames(Microcebus_metadata))]
Microcebus_ancestry$Species_revised <- Microcebus_metadata$Species[match(rownames(Microcebus_SOM_full$ancestry_matrix), rownames(Microcebus_metadata))]
nrow(Microcebus_ancestry) #number of samples in analyzed SOM dataset
length(unique(Microcebus_ancestry$Species)) #number of species present in data
length(unique(Microcebus_ancestry$Species_revised)) #number of proposed species present in data
table(Microcebus_ancestry$Species)
table(Microcebus_ancestry$Species_revised)


## Create species subsets
Microcebus_species_subsets <- list(
  arnh_sp1        = c("M. arnholdi", "M. sp. 1"),
  jona_maca       = c("M. jonahi", "M. macarthurii"),
  ganz_muri       = c("M. ganzhorni", "M. murinus (north)"),
  rave_bong_danf  = c("M. ravelobensis", "M. bongolavensis", "M. danfossi"),
  lehi_mitte      = c("M. lehilahytsara", "M. mittermeieri"),
  bora_clade      = c("M. boraha", "M. tavaratra", "M. sp. 1", "M. arnholdi", 
                      "M. sambiranensis", "M. mamiratra", "M. mittermeieri", 
                      "M. lehilahytsara", "M. myoxinus"),
  tava_clade      = c("M. tavaratra", "M. sp. 1", "M. arnholdi", "M. sambiranensis", 
                      "M. mamiratra", "M. mittermeieri", "M. lehilahytsara", "M. myoxinus"))

Microcebus_SOM_subsets <- lapply(names(Microcebus_species_subsets), function(group_name) { #loop through each subset and create matching data.frames
  Microcebus_group_species <- Microcebus_species_subsets[[group_name]]
  Microcebus_ids <- rownames(Microcebus_metadata)[Microcebus_metadata$Species %in% Microcebus_group_species]
  Microcebus_ids <- intersect(Microcebus_ids, rownames(Microcebus_morphology)) #restrict to IDs that are also in morphology dataset
  list(SNP = Microcebus_SNP[Microcebus_ids, , drop = FALSE], #build subset list
       Morphology = Microcebus_morphology[Microcebus_ids, , drop = FALSE],
       Environmental = Microcebus_environmental[Microcebus_ids, , drop = FALSE],
       Spatial = Microcebus_spatial[Microcebus_ids, , drop = FALSE],
       Metadata = Microcebus_metadata[Microcebus_ids, , drop = FALSE])})
names(Microcebus_SOM_subsets) <- names(Microcebus_species_subsets) #assign names to list


## Train and cluster SOM and evaluate and plot results - arnh_sp1
Microcebus_SOM_subsets$arnh_sp1$Metadata
table(Microcebus_SOM_subsets$arnh_sp1$Metadata$Species)
nrow(Microcebus_SOM_subsets$arnh_sp1$Metadata)
SOM_arnh_sp1 <- train.SOM(
  input_data = Microcebus_SOM_subsets$arnh_sp1[c("SNP", "Morphology", "Environmental", "Spatial")],
  max.NA.row = 0.4,
  max.NA.col = 0.4,
  N.steps = 10,
  parallel = FALSE,
  N_cores = 5,
  grid.multiplier = 3,
  save.SOM.results = TRUE,
  overwrite.SOM.results = TRUE,
  learning.rate.tuning = FALSE,
  N.replicates = 5)
SOM_arnh_sp1 <- clustering.SOM(SOM_arnh_sp1, max.k = 6, clustering.method = "kmeans+BICjump+hierarch")
SOM_arnh_sp1$optim_k_summary
plot.Learning.SOM(SOM_arnh_sp1)
plot.Layers.SOM(SOM_arnh_sp1)
plot.K.SOM(SOM_arnh_sp1)
plot.Model.SOM(SOM_arnh_sp1)
plot.Variable.Importance.SOM(SOM_arnh_sp1,
                             left.margin = 5.8,
                             bar.label.font.size = 0.4, 
                             importance.threshold = 0.1)
plot.Structure.SOM(SOM_arnh_sp1, Individual.labels.font.size = 0.7)
plot.Map.SOM(SOM.output = SOM_arnh_sp1, Coordinates = Microcebus_SOM_subsets$arnh_sp1$Spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 1.5, lon.buffer.range = 1.5, pie.size = 2.5,
             north.arrow.position = c(0.05, 0.86), north.arrow.length = 0.3,
             north.arrow.N.position = 0.15, north.arrow.N.size = 1)
SOM_ancestry_arnh_sp1 <- as.data.frame(SOM_arnh_sp1$ancestry_matrix)
SOM_ancestry_arnh_sp1$Species <- Microcebus_SOM_subsets$arnh_sp1$Metadata$Species[match(rownames(SOM_ancestry_arnh_sp1), rownames(Microcebus_SOM_subsets$arnh_sp1$Metadata))]
SOM_ancestry_arnh_sp1$Species_revised <- SOM_ancestry_arnh_sp1$Species
table(SOM_ancestry_arnh_sp1$Species)
table(SOM_ancestry_arnh_sp1$Species_revised)


## Train and cluster SOM and evaluate and plot results - jona_maca
Microcebus_SOM_subsets$jona_maca$Metadata
table(Microcebus_SOM_subsets$jona_maca$Metadata$Species)
nrow(Microcebus_SOM_subsets$jona_maca$Metadata)
SOM_jona_maca <- train.SOM(
  input_data = Microcebus_SOM_subsets$jona_maca[c("SNP", "Morphology", "Environmental", "Spatial")],
  max.NA.row = 0.3,
  max.NA.col = 0.3,
  N.steps = 10,
  parallel = FALSE,
  N_cores = 5,
  grid.multiplier = 3,
  save.SOM.results = TRUE,
  overwrite.SOM.results = TRUE,
  learning.rate.tuning = FALSE,
  N.replicates = 5)
SOM_jona_maca <- clustering.SOM(SOM_jona_maca, max.k = 6, clustering.method = "kmeans+BICjump+hierarch")
SOM_jona_maca$optim_k_summary
plot.Learning.SOM(SOM_jona_maca)
plot.Layers.SOM(SOM_jona_maca)
plot.K.SOM(SOM_jona_maca)
plot.Model.SOM(SOM_jona_maca)
plot.Variable.Importance.SOM(SOM_jona_maca,
                             left.margin = 5.8,
                             bar.label.font.size = 0.4, 
                             importance.threshold = 0.1)
plot.Structure.SOM(SOM_jona_maca, Individual.labels.font.size = 0.7)
plot.Map.SOM(SOM.output = SOM_jona_maca, Coordinates = Microcebus_SOM_subsets$jona_maca$Spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 1.5, lon.buffer.range = 1.5, pie.size = 2.5,
             north.arrow.position = c(0.05, 0.86), north.arrow.length = 0.3,
             north.arrow.N.position = 0.15, north.arrow.N.size = 1)
SOM_ancestry_jona_maca <- as.data.frame(SOM_jona_maca$ancestry_matrix)
SOM_ancestry_jona_maca$Species <- Microcebus_SOM_subsets$jona_maca$Metadata$Species[match(rownames(SOM_ancestry_jona_maca), rownames(Microcebus_SOM_subsets$jona_maca$Metadata))]
SOM_ancestry_jona_maca$Species_revised <- SOM_ancestry_jona_maca$Species
table(SOM_ancestry_jona_maca$Species)
table(SOM_ancestry_jona_maca$Species_revised)


## Train and cluster SOM and evaluate and plot results - ganz_muri
Microcebus_SOM_subsets$ganz_muri$Metadata
table(Microcebus_SOM_subsets$ganz_muri$Metadata$Species)
nrow(Microcebus_SOM_subsets$ganz_muri$Metadata)
SOM_ganz_muri <- train.SOM(
  input_data = Microcebus_SOM_subsets$ganz_muri[c("SNP", "Environmental", "Spatial")],
  max.NA.row = 0.3,
  max.NA.col = 0.2,
  N.steps = 10,
  parallel = FALSE,
  N_cores = 5,
  grid.multiplier = 2,
  save.SOM.results = TRUE,
  overwrite.SOM.results = TRUE,
  learning.rate.tuning = FALSE,
  N.replicates = 5)
SOM_ganz_muri <- clustering.SOM(SOM_ganz_muri, max.k = 3, clustering.method = "kmeans+BICjump+hierarch")
SOM_ganz_muri$optim_k_summary
plot.Learning.SOM(SOM_ganz_muri)
plot.Layers.SOM(SOM_ganz_muri)
plot.K.SOM(SOM_ganz_muri)
plot.Model.SOM(SOM_ganz_muri)
plot.Variable.Importance.SOM(SOM_ganz_muri,
                             left.margin = 5.8,
                             bar.label.font.size = 0.4, 
                             importance.threshold = 0.1)
plot.Structure.SOM(SOM_ganz_muri, Individual.labels.font.size = 0.7)
plot.Map.SOM(SOM.output = SOM_ganz_muri, Coordinates = Microcebus_SOM_subsets$ganz_muri$Spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 2, lon.buffer.range = 3, pie.size = 2.5,
             north.arrow.position = c(0.06, 0.88), north.arrow.length = 0.8,
             north.arrow.N.position = 0.3, north.arrow.N.size = 1)
SOM_ancestry_ganz_muri <- as.data.frame(SOM_ganz_muri$ancestry_matrix)
SOM_ancestry_ganz_muri$Species <- Microcebus_SOM_subsets$ganz_muri$Metadata$Species[match(rownames(SOM_ancestry_ganz_muri), rownames(Microcebus_SOM_subsets$ganz_muri$Metadata))]
SOM_ancestry_ganz_muri$Species_revised <- SOM_ancestry_ganz_muri$Species
table(SOM_ancestry_ganz_muri$Species)
table(SOM_ancestry_ganz_muri$Species_revised)


## Train and cluster SOM and evaluate and plot results - rave_bong_danf
Microcebus_SOM_subsets$rave_bong_danf$Metadata
table(Microcebus_SOM_subsets$rave_bong_danf$Metadata$Species)
nrow(Microcebus_SOM_subsets$rave_bong_danf$Metadata)
SOM_rave_bong_danf <- train.SOM(
  input_data = Microcebus_SOM_subsets$rave_bong_danf[c("SNP", "Morphology", "Environmental", "Spatial")],
  max.NA.row = 0.3,
  max.NA.col = 0.4,
  N.steps = 10,
  parallel = FALSE,
  N_cores = 5,
  grid.multiplier = 3,
  save.SOM.results = TRUE,
  overwrite.SOM.results = TRUE,
  learning.rate.tuning = FALSE,
  N.replicates = 5)
SOM_rave_bong_danf <- clustering.SOM(SOM_rave_bong_danf, max.k = 8, clustering.method = "kmeans+BICjump+hierarch")
SOM_rave_bong_danf$optim_k_summary
plot.Learning.SOM(SOM_rave_bong_danf)
plot.Layers.SOM(SOM_rave_bong_danf)
plot.K.SOM(SOM_rave_bong_danf)
plot.Model.SOM(SOM_rave_bong_danf)
plot.Variable.Importance.SOM(SOM_rave_bong_danf,
                             left.margin = 5.8,
                             bar.label.font.size = 0.4, 
                             importance.threshold = 0.1)
plot.Structure.SOM(SOM_rave_bong_danf, Individual.labels.font.size = 0.7)
plot.Map.SOM(SOM.output = SOM_rave_bong_danf, Coordinates = Microcebus_SOM_subsets$rave_bong_danf$Spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 1.5, lon.buffer.range = 1.5, pie.size = 2.5,
             north.arrow.position = c(0.05, 0.86), north.arrow.length = 0.3,
             north.arrow.N.position = 0.15, north.arrow.N.size = 1)
SOM_ancestry_rave_bong_danf <- as.data.frame(SOM_rave_bong_danf$ancestry_matrix)
SOM_ancestry_rave_bong_danf$Species <- Microcebus_SOM_subsets$rave_bong_danf$Metadata$Species[match(rownames(SOM_ancestry_rave_bong_danf), rownames(Microcebus_SOM_subsets$rave_bong_danf$Metadata))]
SOM_ancestry_rave_bong_danf$Species_revised <- SOM_ancestry_rave_bong_danf$Species
table(SOM_ancestry_rave_bong_danf$Species)
table(SOM_ancestry_rave_bong_danf$Species_revised)


## Train and cluster SOM and evaluate and plot results - lehi_mitte
Microcebus_SOM_subsets$lehi_mitte$Metadata
table(Microcebus_SOM_subsets$lehi_mitte$Metadata$Species)
nrow(Microcebus_SOM_subsets$lehi_mitte$Metadata)
SOM_lehi_mitte <- train.SOM(
  input_data = Microcebus_SOM_subsets$lehi_mitte[c("SNP", "Morphology", "Environmental", "Spatial")],
  max.NA.row = 0.5,
  max.NA.col = 0.4,
  N.steps = 10,
  parallel = FALSE,
  N_cores = 5,
  grid.multiplier = 2,
  save.SOM.results = TRUE,
  overwrite.SOM.results = TRUE,
  learning.rate.tuning = FALSE,
  N.replicates = 5)
SOM_lehi_mitte <- clustering.SOM(SOM_lehi_mitte, max.k = 3, clustering.method = "kmeans+BICjump+hierarch")
SOM_lehi_mitte$optim_k_summary
plot.Learning.SOM(SOM_lehi_mitte)
plot.Layers.SOM(SOM_lehi_mitte)
plot.K.SOM(SOM_lehi_mitte)
plot.Model.SOM(SOM_lehi_mitte)
plot.Variable.Importance.SOM(SOM_lehi_mitte,
                             left.margin = 5.8,
                             bar.label.font.size = 0.4, 
                             importance.threshold = 0.1)
plot.Structure.SOM(SOM_lehi_mitte, Individual.labels.font.size = 0.7)
plot.Map.SOM(SOM.output = SOM_lehi_mitte, Coordinates = Microcebus_SOM_subsets$lehi_mitte$Spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 2, lon.buffer.range = 3, pie.size = 2.5,
             north.arrow.position = c(0.05, 0.88), north.arrow.length = 0.5,
             north.arrow.N.position = 0.2, north.arrow.N.size = 1)
SOM_ancestry_lehi_mitte <- as.data.frame(SOM_lehi_mitte$ancestry_matrix)
SOM_ancestry_lehi_mitte$Species <- Microcebus_SOM_subsets$lehi_mitte$Metadata$Species[match(rownames(SOM_ancestry_lehi_mitte), rownames(Microcebus_SOM_subsets$lehi_mitte$Metadata))]
SOM_ancestry_lehi_mitte$Species_revised <- SOM_ancestry_lehi_mitte$Species
table(SOM_ancestry_lehi_mitte$Species)
table(SOM_ancestry_lehi_mitte$Species_revised)


## Train and cluster SOM and evaluate and plot results - bora_clade
Microcebus_SOM_subsets$bora_clade$Metadata
table(Microcebus_SOM_subsets$bora_clade$Metadata$Species)
nrow(Microcebus_SOM_subsets$bora_clade$Metadata)
SOM_bora_clade <- train.SOM(
  input_data = Microcebus_SOM_subsets$bora_clade[c("SNP", "Morphology", "Environmental", "Spatial")],
  max.NA.row = 0.3,
  max.NA.col = 0.2,
  N.steps = 10,
  parallel = FALSE,
  N_cores = 5,
  grid.multiplier = 4,
  save.SOM.results = TRUE,
  overwrite.SOM.results = TRUE,
  learning.rate.tuning = FALSE,
  N.replicates = 5)
SOM_bora_clade <- clustering.SOM(SOM_bora_clade, max.k = 15, clustering.method = "kmeans+BICjump+hierarch")
SOM_bora_clade$optim_k_summary
plot.Learning.SOM(SOM_bora_clade)
plot.Layers.SOM(SOM_bora_clade)
plot.K.SOM(SOM_bora_clade)
plot.Model.SOM(SOM_bora_clade)
plot.Variable.Importance.SOM(SOM_bora_clade,
                             left.margin = 5.8,
                             bar.label.font.size = 0.4, 
                             importance.threshold = 0.1)
plot.Structure.SOM(SOM_bora_clade, Individual.labels.font.size = 0.6, margin.bottom = 6)
plot.Map.SOM(SOM.output = SOM_bora_clade, Coordinates = Microcebus_SOM_subsets$bora_clade$Spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 1.5, lon.buffer.range = 2, pie.size = 2.5,
             north.arrow.position = c(0.05, 0.88), north.arrow.length = 0.4,
             north.arrow.N.position = 0.15, north.arrow.N.size = 1)
SOM_ancestry_bora_clade <- as.data.frame(SOM_bora_clade$ancestry_matrix)
SOM_ancestry_bora_clade$Species <- Microcebus_SOM_subsets$bora_clade$Metadata$Species[match(rownames(SOM_ancestry_bora_clade), rownames(Microcebus_SOM_subsets$bora_clade$Metadata))]
SOM_ancestry_bora_clade$Species_revised <- SOM_ancestry_bora_clade$Species
table(SOM_ancestry_bora_clade$Species)
table(SOM_ancestry_bora_clade$Species_revised)


## Train and cluster SOM and evaluate and plot results - tava_clade
Microcebus_SOM_subsets$tava_clade$Metadata
table(Microcebus_SOM_subsets$tava_clade$Metadata$Species)
nrow(Microcebus_SOM_subsets$tava_clade$Metadata)
SOM_tava_clade <- train.SOM(
  input_data = Microcebus_SOM_subsets$tava_clade[c("SNP", "Morphology", "Environmental", "Spatial")],
  max.NA.row = 0.3,
  max.NA.col = 0.2,
  N.steps = 10,
  parallel = FALSE,
  N_cores = 5,
  grid.multiplier = 4,
  save.SOM.results = TRUE,
  overwrite.SOM.results = TRUE,
  learning.rate.tuning = FALSE,
  N.replicates = 5)
SOM_tava_clade <- clustering.SOM(SOM_tava_clade, max.k = 15, clustering.method = "kmeans+BICjump+hierarch")
SOM_tava_clade$optim_k_summary
plot.Learning.SOM(SOM_tava_clade)
plot.Layers.SOM(SOM_tava_clade)
plot.K.SOM(SOM_tava_clade)
plot.Model.SOM(SOM_tava_clade)
plot.Variable.Importance.SOM(SOM_tava_clade,
                             left.margin = 5.8,
                             bar.label.font.size = 0.4, 
                             importance.threshold = 0.1)
plot.Structure.SOM(SOM_tava_clade, Individual.labels.font.size = 0.7, margin.bottom = 7)
plot.Map.SOM(SOM.output = SOM_tava_clade, Coordinates = Microcebus_SOM_subsets$tava_clade$Spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 1.5, lon.buffer.range = 2, pie.size = 2.5,
             north.arrow.position = c(0.05, 0.88), north.arrow.length = 0.4,
             north.arrow.N.position = 0.15, north.arrow.N.size = 1)
SOM_ancestry_tava_clade <- as.data.frame(SOM_tava_clade$ancestry_matrix)
SOM_ancestry_tava_clade$Species <- Microcebus_SOM_subsets$tava_clade$Metadata$Species[match(rownames(SOM_ancestry_tava_clade), rownames(Microcebus_SOM_subsets$tava_clade$Metadata))]
SOM_ancestry_tava_clade$Species_revised <- SOM_ancestry_tava_clade$Species
table(SOM_ancestry_tava_clade$Species)
table(SOM_ancestry_tava_clade$Species_revised)
