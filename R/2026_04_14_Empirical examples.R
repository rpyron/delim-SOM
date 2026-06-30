#### Set environment ###########################################################

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
#setwd("C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package")
intermediate_files_folder <- "Empirical_examples/Intermediate_files"



#### Desmognathus dusky salamanders in Eastern US (Pyron 2023) #################

## https://doi.org/10.1016/j.ympev.2023.107939
## "Monticola71"
## GBS data
## Updated environmental data

## Import sample data
Monticola71_data <- read.csv(file = "../Empirical_examples/Pyron_2023/monticola71.csv",
                             row.names = 1,
                             header = TRUE,
                             colClasses = c(huc2 = "character",
                                            huc4 = "character",
                                            huc6 = "character",
                                            huc8 = "character",
                                            huc10 = "character",
                                            huc12 = "character"))


## Import and process genetic SNP data
Monticola71_SNP <- process.SNP.data.SOM(vcf.path = "../Empirical_examples/Pyron_2023/Monticola71.vcf.gz", #filter loci and individuals and create SNP matrix dataframe
                                        missing.loci.cutoff.lenient = 0.7,
                                        missing.loci.cutoff.final = 0.5,
                                        missing.individuals.cutoff = 0.6)
Monticola71_snp_to_sample <- Monticola71_data$Sample[match(rownames(Monticola71_SNP), rownames(Monticola71_data))] #returns RAP**** names
rownames(Monticola71_SNP) <- Monticola71_snp_to_sample #rename SNP matrix to RAP codes
ncol(Monticola71_SNP) #number of loci: 13031
nrow(Monticola71_SNP) #number of samples: 71


## Create spatial dataset with coordinates and elevation
Monticola71_spatial <- Monticola71_data[, c("Lat", "Long", "Elev"), drop = FALSE] #extract numeric spatial data
Monticola71_spatial <- dplyr::rename(Monticola71_spatial, Latitude = Lat, Longitude = Long, Elevation = Elev)
rownames(Monticola71_spatial) <- Monticola71_data$Sample #assign rownames
nrow(Monticola71_spatial) #number of samples: 71


## Create environmental dataset and binary watershed variables (other variables extracted and processed by separate R script based on coordinates)
Monticola71_environmental <- read.csv("../Empirical_examples/Pyron_2023/Monticola71_environmental.csv", header = TRUE) #read CSV
rownames(Monticola71_environmental) <- Monticola71_environmental$Sample
Monticola71_environmental <- Monticola71_environmental[, !names(Monticola71_environmental) %in% c("Sample", "ID")] #remove ID columns
Monticola71_environmental <- Monticola71_environmental[, !names(Monticola71_environmental) %in% c("Latitude", "Longitude", "Elevation")] #remove spatial variables
Monticola71_environmental[] <- lapply(Monticola71_environmental, as.numeric) #ensure numeric
rownames(Monticola71_data) <- Monticola71_data$Sample #assign rownames
Monticola71_watershed <- make.cols.binary.SOM(dataframe = Monticola71_data,
                                              make.binary.cols = c("lvl4", "huc2", "huc4", "huc6", "huc8", "huc10", "huc12"))
Monticola71_watershed <- Monticola71_watershed[rownames(Monticola71_data), , drop = FALSE]
Monticola71_environmental <- (NicheDiv::transform.skewed.variables(Monticola71_environmental))$transformed #transform skewed variables
Monticola71_environmental <- remove.lowCV.multicollinearity.SOM(Monticola71_environmental, #remove highly correlated and low-variance variables
                                                                CV.threshold = 0.05,
                                                                cor.threshold = 0.9)
ncol(Monticola71_environmental) #number of variables: 71
nrow(Monticola71_environmental) #number of samples: 71
ncol(Monticola71_watershed) #number of variables: 197
nrow(Monticola71_watershed) #number of samples: 71


## Create morphological trait dataset
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
ncol(Monticola71_morphology) #number of traits: 7
nrow(Monticola71_morphology) #number of samples: 71


## Train and cluster SOM
Monticola71_SOM_data <- list(SNP = Monticola71_SNP,
                             Spatial = Monticola71_spatial,
                             Environmental = Monticola71_environmental,
                             Watershed = Monticola71_watershed,
                             Morphology = Monticola71_morphology)
print(unname(round(system.time({
Monticola71_SOM_tr <- train.SOM(input_data = Monticola71_SOM_data, #71 samples, 14.7min
                                  save.SOM.results = TRUE,
                                  save.SOM.results.name = file.path(intermediate_files_folder, "Monticola71_SOM_tr.Rdata"),
                                  max.NA.row = 0.6,
                                  max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Monticola71_SOM_kmeansBICthreshold <- clustering.SOM(Monticola71_SOM_tr, #17.7min
                                                     clustering.method = "kmeans+BICthreshold", 
                                                     save.SOM.results.name = file.path(intermediate_files_folder, "Monticola71_SOM_kmeansBICthreshold.Rdata"))
})[3] / 60, 1)))
Monticola71_SOM_kmeansBICthreshold$optim_k_summary #k2 100%
print(unname(round(system.time({
Monticola71_SOM_HDBSCAN <- clustering.SOM(Monticola71_SOM_tr, #7.6min
                                          clustering.method = "HDBSCAN",
                                          save.SOM.results.name = file.path(intermediate_files_folder, "Monticola71_SOM_HDBSCAN.Rdata"))
})[3] / 60, 1)))
Monticola71_SOM_HDBSCAN$optim_k_summary #k2 83%, k3 13%
print(unname(round(system.time({
Monticola71_SOM_hierarchicalDB <- clustering.SOM(Monticola71_SOM_tr, #222.0min
                                                 clustering.method = "hierarchical+DB",
                                                 save.SOM.results.name = file.path(intermediate_files_folder, "Monticola71_SOM_hierarchicalDB.Rdata"))
})[3] / 60, 1)))
Monticola71_SOM_hierarchicalDB$optim_k_summary #k2 74%, k10 26%
print(unname(round(system.time({
Monticola71_SOM_GMMBICthreshold <- clustering.SOM(Monticola71_SOM_tr, #568.7min
                                                  clustering.method = "GMM+BICthreshold",
                                                  message.N.replicates = 1,
                                                  save.SOM.results.name = file.path(intermediate_files_folder, "Monticola71_SOM_GMMBICthreshold.Rdata"))
})[3] / 60, 1)))
Monticola71_SOM_GMMBICthreshold$optim_k_summary #k3 53%, k2 46%
print(unname(round(system.time({
Monticola71_SOM_OPTICSSilhouette <- clustering.SOM(Monticola71_SOM_tr, #3.1min
                                                   clustering.method = "OPTICS+Silhouette",
                                                   save.SOM.results.name = file.path(intermediate_files_folder, "Monticola71_SOM_OPTICSSilhouette.Rdata"))
})[3] / 60, 1)))
Monticola71_SOM_OPTICSSilhouette$optim_k_summary #k1 80%, k2 20%
print(unname(round(system.time({
Monticola71_SOM_kmeansBICelbow <- clustering.SOM(Monticola71_SOM_tr, #15.8min
                                                 clustering.method = "kmeans+BICelbow",
                                                 save.SOM.results.name = file.path(intermediate_files_folder, "Monticola71_SOM_kmeansBICelbow.Rdata"))
})[3] / 60, 1)))
Monticola71_SOM_kmeansBICelbow$optim_k_summary #k2 100%


## Evaluate and plot results
Monticola71_SOM <- Monticola71_SOM_kmeansBICelbow
plot.learning.SOM(Monticola71_SOM)
plot.layer.distance.scale.SOM(Monticola71_SOM)
plot.K.SOM(Monticola71_SOM)
plot.model.SOM(Monticola71_SOM, replicate.mode = "representative")
plot.model.SOM(Monticola71_SOM, replicate.mode = "first")
plot.structure.SOM(Monticola71_SOM)
plot.map.SOM(SOM.output = Monticola71_SOM,
             Coordinates = Monticola71_spatial[, c("Latitude", "Longitude")],
             USA.add.counties = TRUE,
             scale.position = c(0.78, 0.05))
plot.variable.importance.SOM(Monticola71_SOM, mode = "Cluster.separation", left.margin = 1.5)
plot.variable.importance.SOM(Monticola71_SOM, mode = "Map.variance", left.margin = 1.5)
plot.layer.importance.varimp.SOM(Monticola71_SOM, bottom.margin = 3.5)
plot.layer.importance.leaveoneout.SOM(Monticola71_SOM, 
                                      bottom.margin = 6.5,
                                      save.leave.one.layer.out.results = TRUE,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Monticola71_SOM_lolo.Rdata")) #this will take 10-20min (running 2 x N replicates for train and clustering SOM)


## Hierarchical analyses based on recovered clusters
Monticola71_clusters <- apply(Monticola71_SOM$ancestry_matrix, 1, which.max) #assign each sample to cluster with highest ancestry proportion
Monticola71_clusters <- paste0("cluster", Monticola71_clusters) #rename clusters
table(Monticola71_clusters)
Monticola71_cluster_samples <- split(rownames(Monticola71_SOM$ancestry_matrix), Monticola71_clusters)
Monticola71_cluster1_data <- lapply(Monticola71_SOM$input_data, function(x) x[Monticola71_cluster_samples$cluster1, , drop = FALSE]) #cluster 1 subset
Monticola71_cluster2_data <- lapply(Monticola71_SOM$input_data, function(x) x[Monticola71_cluster_samples$cluster2, , drop = FALSE]) #cluster 2 subset

Monticola71_SOM_tr_cluster1 <- train.SOM(Monticola71_cluster1_data, #48 samples
                                         grid.multiplier = 4,
                                         max.NA.row = 0.5,
                                         max.NA.col = 0.5)
Monticola71_SOM_cluster1 <- clustering.SOM(Monticola71_SOM_tr_cluster1,
                                           clustering.method = "kmeans+BICelbow",
                                           max.k = 5)
Monticola71_SOM_cluster1$optim_k_summary #k1 100% support
Monticola71_SOM_tr_cluster2 <- train.SOM(Monticola71_cluster2_data, #21 samples
                                         grid.multiplier = 4,
                                         max.NA.row = 0.5,
                                         max.NA.col = 0.5)
Monticola71_SOM_cluster2 <- clustering.SOM(Monticola71_SOM_tr_cluster2,
                                           clustering.method = "kmeans+BICelbow",
                                           max.k = 5)
Monticola71_SOM_cluster2$optim_k_summary #k1 100% support




#### Desmognathus salamanders in Alabama/Mississippi (Pyron et al. 2022) #######

## https://doi.org/10.11646/zootaxa.5133.1.3
## "Pascagoula"
## GBS data
## With updated environmental variables
## k2 example (Desmognathus valentinei and D. pascagoula sp. nov.)

## Read in sample data
Pascagoula_data <- read.csv(file = "../Empirical_examples/Pyron_et_al_2022/pascagoula22.csv",
                            row.names = 1,
                            header = TRUE, 
                            colClasses = c(huc2 = "character",
                                           huc4 = "character",
                                           huc6 = "character",
                                           huc8 = "character",
                                           huc10 = "character",
                                           huc12 = "character"))


## Import and process genetic SNP data
Pascagoula_SNP <- process.SNP.data.SOM(vcf.path = "../Empirical_examples/Pyron_et_al_2022/pascagoula22.vcf.gz", #filter loci and individuals and create SNP matrix dataframe
                                       missing.loci.cutoff.lenient = 0.7,
                                       missing.loci.cutoff.final = 0.5,
                                       missing.individuals.cutoff = 0.5)
rownames(Pascagoula_SNP) <- Pascagoula_data$Sample[match(rownames(Pascagoula_SNP), rownames(Pascagoula_data))] #rename alleles
ncol(Pascagoula_SNP) #number of loci: 3728
nrow(Pascagoula_SNP) #number of samples: 22


## Create spatial dataset with coordinates and elevation
Pascagoula_spatial <- Pascagoula_data[, c("Lat", "Long", "Elev")] #extract coordinates and elevation
Pascagoula_spatial <- dplyr::rename(Pascagoula_spatial, Latitude = Lat, Longitude = Long, Elevation = Elev)
rownames(Pascagoula_spatial) <- Pascagoula_data$Sample #assign rownames
nrow(Pascagoula_spatial) #number of samples: 22


## Create environmental dataset and binary watershed variables (other variables extracted and processed by separate R script based on coordinates)
Pascagoula_environmental <- read.csv("../Empirical_examples/Pyron_et_al_2022/Pascagoula22_environmental.csv", header = TRUE) #read CSV
rownames(Pascagoula_environmental) <- Pascagoula_environmental$Sample
Pascagoula_environmental <- Pascagoula_environmental[, !names(Pascagoula_environmental) %in% c("Sample", "ID")] #remove ID columns
Pascagoula_environmental <- Pascagoula_environmental[, !names(Pascagoula_environmental) %in% c("Latitude", "Longitude", "Elevation")] #remove spatial variables
Pascagoula_environmental[] <- lapply(Pascagoula_environmental, as.numeric) #ensure numeric
rownames(Pascagoula_data) <- Pascagoula_data$Sample #assign rownames
Pascagoula_watershed <- make.cols.binary.SOM(dataframe = Pascagoula_data,
                                             make.binary.cols = c("lvl4", "huc2", "huc4", "huc6", "huc8", "huc10", "huc12"))
Pascagoula_watershed <- Pascagoula_watershed[rownames(Pascagoula_data), , drop = FALSE]
Pascagoula_environmental <- (NicheDiv::transform.skewed.variables(Pascagoula_environmental))$transformed #transform skewed variables
Pascagoula_environmental <- remove.lowCV.multicollinearity.SOM(Pascagoula_environmental, #remove highly correlated and low-variance variables
                                                               CV.threshold = 0.05,
                                                               cor.threshold = 0.9)
ncol(Pascagoula_environmental) #number of variables: 59
nrow(Pascagoula_environmental) #number of samples: 22
ncol(Pascagoula_watershed) #number of variables: 42
nrow(Pascagoula_watershed) #number of samples: 22


## Create morphological trait dataset
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
ncol(Pascagoula_morphology) #number of traits: 9
nrow(Pascagoula_morphology) #number of samples: 22


## Train and cluster SOM
Pascagoula_SOM_data <- list(Alleles = Pascagoula_SNP,
                            Spatial = Pascagoula_spatial,
                            Environmental = Pascagoula_environmental,
                            Watershed = Pascagoula_watershed,
                            Morphology = Pascagoula_morphology)
print(unname(round(system.time({
Pascagoula_SOM_tr <- train.SOM(input_data = Pascagoula_SOM_data, #22 samples, 0.7min
                               max.NA.row = 0.5,
                               max.NA.col = 0.5,
                               save.SOM.results = TRUE,
                               save.SOM.results.name = file.path(intermediate_files_folder, "Pascagoula_SOM_tr.Rdata"),
                               grid.multiplier = 4)
})[3] / 60, 1)))

print(unname(round(system.time({
Pascagoula_SOM_kmeansBICthreshold <- clustering.SOM(Pascagoula_SOM_tr, #2.6min
                                                    clustering.method = "kmeans+BICthreshold",
                                                    save.SOM.results.name = file.path(intermediate_files_folder, "Pascagoula_SOM_kmeansBICthreshold.Rdata"))
})[3] / 60, 1)))
Pascagoula_SOM_kmeansBICthreshold$optim_k_summary #k2 99%
print(unname(round(system.time({
Pascagoula_SOM_HDBSCAN <- clustering.SOM(Pascagoula_SOM_tr, #2.0min
                                         clustering.method = "HDBSCAN",
                                         save.SOM.results.name = file.path(intermediate_files_folder, "Pascagoula_SOM_HDBSCAN.Rdata"))
})[3] / 60, 1)))
Pascagoula_SOM_HDBSCAN$optim_k_summary #k2 90%, k1 5%, k3 5%
print(unname(round(system.time({
Pascagoula_SOM_hierarchicalDB <- clustering.SOM(Pascagoula_SOM_tr, #29.3min
                                                clustering.method = "hierarchical+DB",
                                                save.SOM.results.name = file.path(intermediate_files_folder, "Pascagoula_SOM_hierarchicalDB.Rdata"))
})[3] / 60, 1)))
Pascagoula_SOM_hierarchicalDB$optim_k_summary #k10 94%, k9 6%
print(unname(round(system.time({
Pascagoula_SOM_GMMBICthreshold <- clustering.SOM(Pascagoula_SOM_tr, #25.7min
                                                 clustering.method = "GMM+BICthreshold",
                                                 save.SOM.results.name = file.path(intermediate_files_folder, "Pascagoula_SOM_GMMBICthreshold.Rdata"))
})[3] / 60, 1)))
Pascagoula_SOM_GMMBICthreshold$optim_k_summary #k4 35%, k3 24%, k1 21%, k5 16%
print(unname(round(system.time({
Pascagoula_SOM_OPTICSSilhouette <- clustering.SOM(Pascagoula_SOM_tr, #1.5min
                                                  clustering.method = "OPTICS+Silhouette",
                                                  save.SOM.results.name = file.path(intermediate_files_folder, "Pascagoula_SOM_OPTICSSilhouette.Rdata"))
})[3] / 60, 1)))
Pascagoula_SOM_OPTICSSilhouette$optim_k_summary #k1 52%, k2 46%
print(unname(round(system.time({
Pascagoula_SOM_kmeansBICelbow <- clustering.SOM(Pascagoula_SOM_tr, #2.6min
                                                clustering.method = "kmeans+BICelbow",
                                                save.SOM.results.name = file.path(intermediate_files_folder, "Pascagoula_SOM_kmeansBICelbow.Rdata"))
})[3] / 60, 1)))
Pascagoula_SOM_kmeansBICelbow$optim_k_summary #k2 99%


## Evaluate and plot results
Pascagoula_SOM <- Pascagoula_SOM_kmeansBICelbow
plot.learning.SOM(Pascagoula_SOM)
plot.layer.distance.scale.SOM(Pascagoula_SOM)
plot.K.SOM(Pascagoula_SOM)
plot.model.SOM(Pascagoula_SOM, replicate.mode = "representative")
plot.model.SOM(Pascagoula_SOM, replicate.mode = "first")
plot.structure.SOM(Pascagoula_SOM)
plot.map.SOM(SOM.output = Pascagoula_SOM,
             Coordinates = Pascagoula_spatial[, c("Latitude", "Longitude")],
             USA.add.counties = TRUE,
             north.arrow.position = c(0.05, 0.9),
             north.arrow.length = 0.35,
             north.arrow.N.position = 0.15,
             scale.position = c(0.76, 0.06))
plot.variable.importance.SOM(Pascagoula_SOM,
                             mode = "Cluster.separation", 
                             bars.threshold.N = 15)
plot.variable.importance.SOM(Pascagoula_SOM,
                             mode = "Map.variance",
                             bars.threshold.N = 15)
plot.layer.importance.varimp.SOM(Pascagoula_SOM, bottom.margin = 3.7)
plot.layer.importance.leaveoneout.SOM(Pascagoula_SOM, #this will take 20min (running 2 x N replicates for train and clustering SOM)
                                      save.leave.one.layer.out.results = TRUE,
                                      bottom.margin = 6.5,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Pascagoula_SOM_lolo.Rdata")) 


## Hierarchical analyses based on recovered clusters
Pascagoula_clusters <- apply(Pascagoula_SOM$ancestry_matrix, 1, which.max) #assign each sample to cluster with highest ancestry proportion
Pascagoula_clusters <- paste0("cluster", Pascagoula_clusters) #rename clusters
table(Pascagoula_clusters)
Pascagoula_cluster_samples <- split(rownames(Pascagoula_SOM$ancestry_matrix), Pascagoula_clusters)
Pascagoula_cluster1_data <- lapply(Pascagoula_SOM$input_data, function(x) x[Pascagoula_cluster_samples$cluster1, , drop = FALSE]) #cluster 1 subset
Pascagoula_cluster2_data <- lapply(Pascagoula_SOM$input_data, function(x) x[Pascagoula_cluster_samples$cluster2, , drop = FALSE]) #cluster 2 subset

Pascagoula_SOM_tr_cluster1 <- train.SOM(Pascagoula_cluster1_data, grid.multiplier = 3) #10 samples
Pascagoula_SOM_cluster1 <- clustering.SOM(Pascagoula_SOM_tr_cluster1, clustering.method = "kmeans+BICelbow", max.k = 5)
Pascagoula_SOM_cluster1$optim_k_summary #k1 100%

Pascagoula_SOM_tr_cluster2 <- train.SOM(Pascagoula_cluster2_data, grid.multiplier = 3) #12 samples
Pascagoula_SOM_cluster2 <- clustering.SOM(Pascagoula_SOM_tr_cluster2, clustering.method = "kmeans+BICelbow", max.k = 5)
Pascagoula_SOM_cluster2$optim_k_summary #k1 100%



#### Desmognathus seepage salamanders in southeastern US (Pyron et al. 2024) ####

## www.https://doi.org/10.1111/mec.17219
## "Aeneus56"
## GBS data
## With updated environmental data
## One species consisting of three structured lineages)

## Read in sample data
Aeneus_data <- read.csv(file = "../Empirical_examples/Pyron_et_al_2024/aeneus56.csv",
                        row.names = 1,
                        header = TRUE, 
                        colClasses = c(huc2 = "character",
                                       huc4 = "character",
                                       huc6 = "character",
                                       huc8 = "character",
                                       huc10 = "character",
                                       huc12 = "character"))


## Import and process genetic SNP data
Aeneus_SNP <- process.SNP.data.SOM(vcf.path = "../Empirical_examples/Pyron_et_al_2024/aeneus56.vcf.gz", #filter loci and individuals and create SNP matrix dataframe
                                   missing.loci.cutoff.lenient = 0.7,
                                   missing.loci.cutoff.final = 0.5,
                                   missing.individuals.cutoff = 0.5)
rownames(Aeneus_SNP) <- Aeneus_data$Sample[match(rownames(Aeneus_SNP), rownames(Aeneus_data))] #rename alleles
ncol(Aeneus_SNP) #number of loci: 7667 
nrow(Aeneus_SNP) #number of samples: 51


## Create spatial dataset with coordinates and elevation
Aeneus_spatial <- Aeneus_data[, c("Lat", "Long", "Elev")] #extract coordinates and elevation
Aeneus_spatial <- dplyr::rename(Aeneus_spatial, Latitude = Lat, Longitude = Long, Elevation = Elev) #rename variables
rownames(Aeneus_spatial) <- Aeneus_data$Sample #assign rownames
nrow(Aeneus_spatial) #number of samples: 56


## Create environmental dataset and binary watershed variables (other variables extracted and processed by separate R script based on coordinates)
Aeneus_environmental <- read.csv("../Empirical_examples/Pyron_et_al_2024/Aeneus56_environmental.csv", row.names = 1, header = TRUE) #read CSV with rownames
Aeneus_environmental <- Aeneus_environmental[, !names(Aeneus_environmental) %in% c("Latitude", "Longitude", "Elevation")] #remove spatial variables
Aeneus_environmental <- as.data.frame(lapply(Aeneus_environmental, as.numeric)) #ensure numeric
rownames(Aeneus_environmental) <- Aeneus_data$Sample #assign rownames
rownames(Aeneus_data) <- Aeneus_data$Sample #assign rownames
Aeneus_watershed <- make.cols.binary.SOM(dataframe = Aeneus_data,
                                         make.binary.cols = c("lvl4", "huc2", "huc4", "huc6", "huc8", "huc10", "huc12"))
Aeneus_watershed <- Aeneus_watershed[rownames(Aeneus_data), , drop = FALSE]
Aeneus_environmental <- (NicheDiv::transform.skewed.variables(Aeneus_environmental))$transformed #transform skewed variables
Aeneus_environmental <- remove.lowCV.multicollinearity.SOM(Aeneus_environmental, #remove highly correlated and low-variance variables
                                                           CV.threshold = 0.05,
                                                           cor.threshold = 0.9)
ncol(Aeneus_environmental) #number of variables: 53
nrow(Aeneus_environmental) #number of samples: 56
ncol(Aeneus_watershed) #number of variables: 98
nrow(Aeneus_watershed) #number of samples: 56


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
                               function(trait) stats::resid(stats::lm(Aeneus_filtered_log_traits[, trait] ~ Aeneus_SVL, na.action = stats::na.exclude))) #regress each trait on SVL and retain NA alignment
rownames(Aeneus_residuals_mat) <- rownames(Aeneus_filtered_log_traits) #assign rownames to residual matrix
Aeneus_morphology <- as.data.frame(cbind(SVL = Aeneus_SVL, Aeneus_residuals_mat)) #combine log(SVL) and residuals
Aeneus_morphology <- Aeneus_morphology[rownames(Aeneus_trait_data), ] #keep only samples with trait data
ncol(Aeneus_morphology) #number of traits: 13  
nrow(Aeneus_morphology) #number of samples: 45


## Train and cluster SOM
Aeneus_SOM_data <- list(Alleles = Aeneus_SNP,
                        Spatial = Aeneus_spatial,
                        Environmental = Aeneus_environmental,
                        Watershed = Aeneus_watershed,
                        Morphology = Aeneus_morphology)
print(unname(round(system.time({
Aeneus_SOM_tr <- train.SOM(input_data = Aeneus_SOM_data, #40 samples, 2.8min
                           save.SOM.results = TRUE,
                           save.SOM.results.name = file.path(intermediate_files_folder, "Aeneus_SOM_tr.Rdata"),
                           max.NA.row = 0.5,
                           max.NA.col = 0.5)
})[3] / 60, 1)))

print(unname(round(system.time({
Aeneus_SOM_kmeansBICthreshold <- clustering.SOM(Aeneus_SOM_tr, #4.2min
                                                clustering.method = "kmeans+BICthreshold",
                                                save.SOM.results.name = file.path(intermediate_files_folder, "Aeneus_SOM_kmeansBICthreshold.Rdata"))
})[3] / 60, 1)))
Aeneus_SOM_kmeansBICthreshold$optim_k_summary #k1 100%
print(unname(round(system.time({
Aeneus_SOM_HDBSCAN <- clustering.SOM(Aeneus_SOM_tr, #3.3min
                                     clustering.method = "HDBSCAN",
                                     save.SOM.results.name = file.path(intermediate_files_folder, "Aeneus_SOM_HDBSCAN.Rdata"))
})[3] / 60, 1)))
Aeneus_SOM_HDBSCAN$optim_k_summary #k2 64%, k1 23%, k3 13%
print(unname(round(system.time({
Aeneus_SOM_hierarchicalDB <- clustering.SOM(Aeneus_SOM_tr, #68.7min
                                            clustering.method = "hierarchical+DB",
                                            save.SOM.results.name = file.path(intermediate_files_folder, "Aeneus_SOM_hierarchicalDB.Rdata"))
})[3] / 60, 1)))
Aeneus_SOM_hierarchicalDB$optim_k_summary #k10 91%, k9 9%
print(unname(round(system.time({
Aeneus_SOM_GMMBICthreshold <- clustering.SOM(Aeneus_SOM_tr, #208.6min
                                             clustering.method = "GMM+BICthreshold",
                                             save.SOM.results.name = file.path(intermediate_files_folder, "Aeneus_SOM_GMMBICthreshold.Rdata"))
})[3] / 60, 1)))
Aeneus_SOM_GMMBICthreshold$optim_k_summary #k2 47%, k3 42%, k4 8%
print(unname(round(system.time({
Aeneus_SOM_OPTICSSilhouette <- clustering.SOM(Aeneus_SOM_tr, #1.6min
                                              clustering.method = "OPTICS+Silhouette",
                                              save.SOM.results.name = file.path(intermediate_files_folder, "Aeneus_SOM_OPTICSSilhouette.Rdata"))
})[3] / 60, 1)))
Aeneus_SOM_OPTICSSilhouette$optim_k_summary #k1 93%, k2 7%
print(unname(round(system.time({
Aeneus_SOM_kmeansBICelbow <- clustering.SOM(Aeneus_SOM_tr, #3.9min
                                            clustering.method = "kmeans+BICelbow",
                                            save.SOM.results.name = file.path(intermediate_files_folder, "Aeneus_SOM_kmeansBICelbow.Rdata"))
})[3] / 60, 1)))
Aeneus_SOM_kmeansBICelbow$optim_k_summary #k1 100%


## Evaluate and plot results
Aeneus_SOM <- Aeneus_SOM_kmeansBICelbow
plot.learning.SOM(Aeneus_SOM)
plot.layer.distance.scale.SOM(Aeneus_SOM)
plot.K.SOM(Aeneus_SOM)
plot.model.SOM(Aeneus_SOM, replicate.mode = "first")
plot.model.SOM(Aeneus_SOM, replicate.mode = "representative")
plot.map.SOM(Aeneus_SOM,
             Coordinates = Aeneus_spatial[, c("Latitude", "Longitude")],
             USA.add.counties = TRUE,
             north.arrow.position = c(0.05, 0.9),
             north.arrow.length = 0.4,
             north.arrow.N.position = 0.15,
             scale.position = c(0.79, 0.05))
plot.variable.importance.SOM(Aeneus_SOM,
                             mode = "Cluster.separation",
                             left.margin = 5.8,
                             bar.label.font.size = 0.4) 
plot.variable.importance.SOM(Aeneus_SOM,
                             mode = "Map.variance",
                             left.margin = 5.8,
                             bar.label.font.size = 0.4)
plot.layer.importance.varimp.SOM(Aeneus_SOM, bottom.margin = 6)
plot.layer.importance.leaveoneout.SOM(Aeneus_SOM, #this will take 10-20min (running 2 x N replicates for train and clustering SOM)
                                      save.leave.one.layer.out.results = TRUE,
                                      bottom.margin = 6.5,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Aeneus_SOM_lolo.Rdata"))




#### Pocillopora corals in Indo-Pacific (Oury et al. 2023) ######################

## https://doi.org/10.1016/j.ympev.2023.107803
## 364 colonies
## Target enrichment of 1,248 UCE and 1,385 exon loci
## Morphological data: 14 traits
## Four biogeographic region as binary (somewhat equivalent to spatial data)
## Two haplotype markers for symbiosis: 31 columns as binary
## 21 species hypotheses inferred by genomics
## 13 species strongly supported by all approaches (while six could represent either undescribed or nominal species that have been synonymised incorrectly)

## Import and process genetic data
Pocillopora_vcf_file <- "../Empirical_examples/Oury_et_al_2023/Pocillopora_361ADN_1559SNP.vcf" #VCF file path
Pocillopora_gds_file <- "../Empirical_examples/Oury_et_al_2023/Pocillopora.gds" #GDS file path
SeqArray::seqVCF2GDS(Pocillopora_vcf_file, Pocillopora_gds_file, storage.option = "LZ4_RA.max", verbose = FALSE) #convert VCF to GDS
Pocillopora_gds <- SeqArray::seqOpen(Pocillopora_gds_file) #open GDS file
Pocillopora_geno <- SeqArray::seqGetData(Pocillopora_gds, "genotype") #get genotype array
Pocillopora_SNP_raw <- Pocillopora_geno[1, , ] + Pocillopora_geno[2, , ] #sum allele counts
rownames(Pocillopora_SNP_raw) <- SeqArray::seqGetData(Pocillopora_gds, "sample.id") #assign rownames
colnames(Pocillopora_SNP_raw) <- paste0("SNP", seq_len(ncol(Pocillopora_SNP_raw))) #assign colnames
Pocillopora_SNP_raw <- as.data.frame(Pocillopora_SNP_raw) #convert to data frame
SeqArray::seqClose(Pocillopora_gds) #close GDS connection
Pocillopora_genotypes <- as.data.frame(matrix(ifelse(is.na(unlist(Pocillopora_SNP_raw)), NA, ifelse(unlist(Pocillopora_SNP_raw) == 0, "A/A", ifelse(unlist(Pocillopora_SNP_raw) == 1, "A/B", ifelse(unlist(Pocillopora_SNP_raw) == 2, "B/B", NA)))), nrow = nrow(Pocillopora_SNP_raw), dimnames = dimnames(Pocillopora_SNP_raw))) #convert dosages to genotype strings
Pocillopora_genind <- adegenet::df2genind(Pocillopora_genotypes, sep = "/", ncode = 1, ploidy = 2) #convert to genind
Pocillopora_SNP <- process.SNP.data.SOM(genind.input = Pocillopora_genind, #filter loci and individuals and create SNP matrix dataframe
                                        missing.loci.cutoff.lenient = 0.7,
                                        missing.loci.cutoff.final = 0.5,
                                        missing.individuals.cutoff = 0.5)
nrow(Pocillopora_SNP) #number of samples: 350
ncol(Pocillopora_SNP) #number of loci: 1559


## Import and process morphology dataset
Pocillopora_morphology <- readr::read_delim(file = "../Empirical_examples/Oury_et_al_2023/Micromorphometry_Pocillopora_170ind.csv", #import csv
                                            delim = ";",
                                            quote = "\"",
                                            escape_double = TRUE,
                                            locale = readr::locale(decimal_mark = ","),
                                            trim_ws = TRUE)
Pocillopora_morphology <- as.data.frame(Pocillopora_morphology) #make dataframe
rownames(Pocillopora_morphology) <- Pocillopora_morphology$Sample_Name #set rownames
Pocillopora_morphology$GSH <- NULL #remove Genomic_species_hypothesis column
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
Pocillopora_numeric_traits <- c(
  "Max_calice_diameter_1",
  "Max_calice_diameter_2",
  "Distance_corallite",
  "Distance_denticles",
  "Height_septa",
  "Max_columella_diameter_1",
  "Max_columella_diameter_2")
Pocillopora_morphology[, Pocillopora_numeric_traits] <- as.data.frame(setNames(lapply(Pocillopora_numeric_traits, function(trait_name) { #process numeric traits
  trait_values <- as.numeric(Pocillopora_morphology[, trait_name])
  if(any(trait_values <= 0, na.rm = TRUE)) { #leave traits with zero or negative values untransformed
    trait_values
  } else { #log transform strictly positive traits
    log(trait_values)
  }
}), Pocillopora_numeric_traits), row.names = rownames(Pocillopora_morphology)) #process numeric traits
head(Pocillopora_morphology)
Pocillopora_morphology_binary  <- make.cols.binary.SOM(dataframe = Pocillopora_morphology, #make columella and septa shapes binary
                                                       make.binary.cols = c("Shape_columella", "Shape_septa"),
                                                       append.to.original = FALSE)
Pocillopora_morphology <- Pocillopora_morphology[, Pocillopora_numeric_traits, drop = FALSE] #keep only continuous traits
Pocillopora_morphology <- remove.lowCV.multicollinearity.SOM(Pocillopora_morphology, CV.threshold = 0.05, cor.threshold = 0.9) #remove highly correlated and low variance variables
ncol(Pocillopora_morphology) #number of variables: 6
nrow(Pocillopora_morphology) #number of samples: 175


## Import csv file with multiple traits and meta data
Pocillopora_multiple_traits <- readr::read_delim(file = "../Empirical_examples/Oury_et_al_2023/DB_Pocillopora_genomic_364ind.csv",
                                                 delim = ";",
                                                 quote = "\"",
                                                 escape_double = TRUE,
                                                 show_col_types = FALSE,
                                                 locale = readr::locale(decimal_mark = ","),
                                                 trim_ws = TRUE)
Pocillopora_multiple_traits <- as.data.frame(Pocillopora_multiple_traits)
colnames(Pocillopora_multiple_traits)[colnames(Pocillopora_multiple_traits) == "GSH"] <- "Genomic_species_hypothesis"
colnames(Pocillopora_multiple_traits)[colnames(Pocillopora_multiple_traits) == "PSH"] <- "Primary_species_hypothesis"
colnames(Pocillopora_multiple_traits)[colnames(Pocillopora_multiple_traits) == "SSH"] <- "Secondary_species_hypothesis"
Pocillopora_species_names <- dplyr::distinct(Pocillopora_multiple_traits, Sample_Name, .keep_all = TRUE)
Pocillopora_species_names <- tibble::column_to_rownames(Pocillopora_species_names, "Sample_Name")
Pocillopora_species_names <- dplyr::select(Pocillopora_species_names, Locality, Secondary_species_hypothesis, Genomic_species_hypothesis, Primary_species_hypothesis)
Pocillopora_species_names[Pocillopora_species_names == "-"] <- NA
Pocillopora_species_names[Pocillopora_species_names == "?"] <- NA


## Extract and process microsatellite genotypes
Pocillopora_microsat_cols <- grep("\\.(1|2)$", names(Pocillopora_multiple_traits), value = TRUE)
Pocillopora_microsatellites <- dplyr::select(Pocillopora_multiple_traits, Sample_Name, dplyr::all_of(Pocillopora_microsat_cols))
Pocillopora_microsatellites <- dplyr::mutate(Pocillopora_microsatellites, dplyr::across(-Sample_Name, ~ dplyr::na_if(., "-")))
Pocillopora_microsatellites <- dplyr::mutate(Pocillopora_microsatellites, dplyr::across(-Sample_Name, ~ dplyr::na_if(., "?")))
Pocillopora_microsatellites$na_count <- rowSums(is.na(Pocillopora_microsatellites[, colnames(Pocillopora_microsatellites) != "Sample_Name"]))
Pocillopora_microsatellites <- dplyr::group_by(Pocillopora_microsatellites, Sample_Name)
Pocillopora_microsatellites <- dplyr::slice_min(Pocillopora_microsatellites, na_count, n = 1, with_ties = FALSE)
Pocillopora_microsatellites <- dplyr::ungroup(Pocillopora_microsatellites)
Pocillopora_microsatellites <- dplyr::select(Pocillopora_microsatellites, -na_count)
Pocillopora_microsatellites <- as.data.frame(Pocillopora_microsatellites)
rownames(Pocillopora_microsatellites) <- Pocillopora_microsatellites$Sample_Name
Pocillopora_microsatellites$Sample_Name <- NULL
Pocillopora_microsatellites <- dplyr::mutate(Pocillopora_microsatellites, dplyr::across(dplyr::everything(), as.numeric))
ncol(Pocillopora_microsatellites) #number of variables: 26
nrow(Pocillopora_microsatellites) #number of samples: 367


## Extract and process haplotype markers (ORF and PocHistone) for symbiosis as binary
Pocillopora_symbiosis_haplotypes <- Pocillopora_multiple_traits #copy original
Pocillopora_symbiosis_haplotypes <- dplyr::mutate(Pocillopora_symbiosis_haplotypes, dplyr::across(c(ORF, PocHistone), function(marker_values) dplyr::na_if(marker_values, "-"))) #replace "-" with NA
Pocillopora_symbiosis_haplotypes <- dplyr::mutate(Pocillopora_symbiosis_haplotypes, dplyr::across(c(ORF, PocHistone), function(marker_values) dplyr::na_if(marker_values, "?"))) #replace "?" with NA
Pocillopora_symbiosis_haplotypes <- dplyr::select(Pocillopora_symbiosis_haplotypes, Sample_Name, ORF, PocHistone) #select only relevant columns
Pocillopora_symbiosis_haplotypes$na_count <- rowSums(is.na(Pocillopora_symbiosis_haplotypes[, c("ORF", "PocHistone")])) #count NAs
Pocillopora_symbiosis_haplotypes <- dplyr::group_by(Pocillopora_symbiosis_haplotypes, Sample_Name)
Pocillopora_symbiosis_haplotypes <- dplyr::slice_min(Pocillopora_symbiosis_haplotypes, na_count, with_ties = FALSE) #keep best per sample
Pocillopora_symbiosis_haplotypes <- dplyr::ungroup(Pocillopora_symbiosis_haplotypes)
Pocillopora_symbiosis_haplotypes <- dplyr::select(Pocillopora_symbiosis_haplotypes, -na_count) #drop NA count
Pocillopora_symbiosis_haplotypes[] <- lapply(Pocillopora_symbiosis_haplotypes, function(column_values) {
  if (is.character(column_values)) column_values <- as.factor(column_values) #convert character to factor
  if (is.factor(column_values)) column_values <- droplevels(column_values) #drop unused levels
  return(column_values)
})
Pocillopora_symbiosis_haplotypes <- as.data.frame(Pocillopora_symbiosis_haplotypes) #ensure data.frame
rownames(Pocillopora_symbiosis_haplotypes) <- Pocillopora_symbiosis_haplotypes$Sample_Name #set rownames
Pocillopora_symbiosis_haplotypes$Sample_Name <- NULL #remove Sample_Name
Pocillopora_symbiosis_haplotypes <- make.cols.binary.SOM( #convert ORF and PocHistone to binary variables
  dataframe = Pocillopora_symbiosis_haplotypes,
  make.binary.cols = c("ORF", "PocHistone"),
  append.to.original = TRUE)
colnames(Pocillopora_symbiosis_haplotypes) <- gsub("^ORF_NA$", "ORF_missing", colnames(Pocillopora_symbiosis_haplotypes)) #handle ORF missing
colnames(Pocillopora_symbiosis_haplotypes) <- gsub("^NA$", "ORF_missing", colnames(Pocillopora_symbiosis_haplotypes)) #handle ORF missing if created without prefix
colnames(Pocillopora_symbiosis_haplotypes) <- gsub("^ORF_", "", colnames(Pocillopora_symbiosis_haplotypes)) #remove ORF prefix
colnames(Pocillopora_symbiosis_haplotypes) <- gsub("^PocHistone_NA$", "Poc_missing", colnames(Pocillopora_symbiosis_haplotypes)) #handle Poc missing
colnames(Pocillopora_symbiosis_haplotypes) <- gsub("^PocHistone_", "Poc", colnames(Pocillopora_symbiosis_haplotypes)) #rename PocHistone
Pocillopora_symbiosis_haplotypes <- dplyr::select(Pocillopora_symbiosis_haplotypes, !dplyr::matches("(^missing$|^Poc_missing$)")) #drop missing flags
ncol(Pocillopora_symbiosis_haplotypes) #number of variables: 31
nrow(Pocillopora_symbiosis_haplotypes) #number of samples: 367


## Extract and process biogeographic region as binary
Pocillopora_biogeography <- dplyr::select(Pocillopora_multiple_traits, Sample_Name, Province)
Pocillopora_biogeography <- dplyr::mutate(Pocillopora_biogeography, Province = dplyr::na_if(Province, "-"))
Pocillopora_biogeography$na_count <- is.na(Pocillopora_biogeography$Province)
Pocillopora_biogeography <- dplyr::group_by(Pocillopora_biogeography, Sample_Name)
Pocillopora_biogeography <- dplyr::slice_min(Pocillopora_biogeography, na_count, n = 1, with_ties = FALSE)
Pocillopora_biogeography <- dplyr::ungroup(Pocillopora_biogeography)
Pocillopora_biogeography <- dplyr::select(Pocillopora_biogeography, -na_count)
Pocillopora_biogeography <- as.data.frame(Pocillopora_biogeography)
rownames(Pocillopora_biogeography) <- Pocillopora_biogeography$Sample_Name
Pocillopora_biogeography$Sample_Name <- NULL
Pocillopora_biogeography <- make.cols.binary.SOM(dataframe = Pocillopora_biogeography, #convert Province to binary matrix
                                                 make.binary.cols = "Province")
ncol(Pocillopora_biogeography) #number of variables: 3
nrow(Pocillopora_biogeography) #number of samples: 367


## Create binary morphology dataset
Pocillopora_morpho_map <- dplyr::select(Pocillopora_multiple_traits, Sample_Name, Morphotype) #select relevant columns
Pocillopora_morpho_map$na_count <- is.na(Pocillopora_morpho_map$Morphotype) #count NA entries
Pocillopora_morpho_map <- dplyr::group_by(Pocillopora_morpho_map, Sample_Name) #group by sample
Pocillopora_morpho_map <- dplyr::slice_min(Pocillopora_morpho_map, na_count, n = 1, with_ties = FALSE) #keep most complete entry
Pocillopora_morpho_map <- dplyr::ungroup(Pocillopora_morpho_map) #ungroup
Pocillopora_morpho_map <- dplyr::select(Pocillopora_morpho_map, -na_count) #drop NA count
Pocillopora_morpho_map <- as.data.frame(Pocillopora_morpho_map) #convert to data.frame
rownames(Pocillopora_morpho_map) <- Pocillopora_morpho_map$Sample_Name #set Sample_Name as rownames
Pocillopora_morpho_map$Sample_Name <- NULL #remove Sample_Name column
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
Pocillopora_morphology_binary$Morphotype <- Pocillopora_morpho_map[rownames(Pocillopora_morphology_binary), "Morphotype"] #match morphotypes to binary morphology samples
Pocillopora_morphology_binary$Morphotype_full <- dplyr::recode(Pocillopora_morphology_binary$Morphotype, !!!morpho_lookup, .default = NA_character_) #recode morphotypes
Pocillopora_morphology_binary$Morphotype <- NULL #drop original Morphotype column
simple_traits <- unique(unlist(stringr::str_split(stats::na.omit(as.character(Pocillopora_morphology_binary$Morphotype_full)), "/"))) #get individual morphotype components
for(trait_name in simple_traits){ #iterate over traits
  column_name <- paste0("Morphotype_", stringr::str_replace_all(trait_name, "\\s+", "_")) #define column name
  Pocillopora_morphology_binary[[column_name]] <- ifelse(is.na(Pocillopora_morphology_binary$Morphotype_full), NA_integer_, as.integer(stringr::str_detect(Pocillopora_morphology_binary$Morphotype_full, trait_name))) #binary encode
}
Pocillopora_morphology_binary$Morphotype_full <- NULL #drop full label column
ncol(Pocillopora_morphology_binary) #number of traits: 19
nrow(Pocillopora_morphology_binary) #number of samples: 175


## Train and cluster SOM
Pocillopora_SOM_data <- list(SNP = Pocillopora_SNP,
                             Microsatellites = Pocillopora_microsatellites,
                             Morphology = Pocillopora_morphology,
                             Morphology_binary = Pocillopora_morphology_binary,
                             Biogeography = Pocillopora_biogeography,
                             Symbiosis_haplotypes = Pocillopora_symbiosis_haplotypes)
print(unname(round(system.time({
Pocillopora_SOM_tr <- train.SOM(input_data = Pocillopora_SOM_data, #80 samples
                                max.NA.row = 0.5,
                                max.NA.col = 0.5,
                                save.SOM.results = TRUE,
                                save.SOM.results.name = file.path(intermediate_files_folder, "Pocillopora_SOM_tr.Rdata"))
})[3] / 60, 1)))

print(unname(round(system.time({
Pocillopora_SOM_kmeansBICthreshold <- clustering.SOM(Pocillopora_SOM_tr, #
                                                     max.k = 35,
                                                     clustering.method = "kmeans+BICthreshold",
                                                     save.SOM.results.name = file.path(intermediate_files_folder, "Pocillopora_SOM_kmeansBICthreshold.Rdata"))
})[3] / 60, 1)))
Pocillopora_SOM_kmeansBICthreshold$optim_k_summary #k3 88%, k4 7%
print(unname(round(system.time({
Pocillopora_SOM_HDBSCAN <- clustering.SOM(Pocillopora_SOM_tr, #
                                          max.k = 35,
                                          clustering.method = "HDBSCAN",
                                          save.SOM.results.name = file.path(intermediate_files_folder, "Pocillopora_SOM_HDBSCAN.Rdata"))
})[3] / 60, 1)))
Pocillopora_SOM_HDBSCAN$optim_k_summary #k2 35%, k3 27%, k4 15%, k5 13%, k6 9%
print(unname(round(system.time({
Pocillopora_SOM_kmeansBICelbow <- clustering.SOM(Pocillopora_SOM_tr, #
                                                 max.k = 35,
                                                 clustering.method = "kmeans+BICelbow",
                                                 save.SOM.results.name = file.path(intermediate_files_folder, "Pocillopora_SOM_kmeansBICelbow.Rdata"))
})[3] / 60, 1)))
Pocillopora_SOM_kmeansBICelbow$optim_k_summary #k3 51%, k35 31%, k34 7%
print(unname(round(system.time({
Pocillopora_SOM_OPTICSSilhouette <- clustering.SOM(Pocillopora_SOM_tr, #
                                                   max.k = 35,
                                                   clustering.method = "OPTICS+Silhouette",
                                                   save.SOM.results.name = file.path(intermediate_files_folder, "Pocillopora_SOM_OPTICSSilhouette.Rdata"))
})[3] / 60, 1)))
Pocillopora_SOM_OPTICSSilhouette$optim_k_summary #k2 50%, k1 32%, k3 18%
print(unname(round(system.time({
Pocillopora_SOM_GMMBICthreshold <- clustering.SOM(Pocillopora_SOM_tr, #
                                                  max.k = 35,
                                                  clustering.method = "GMM+BICthreshold",
                                                  save.SOM.results.name = file.path(intermediate_files_folder, "Pocillopora_SOM_GMMBICthreshold.Rdata"))
})[3] / 60, 1)))
Pocillopora_SOM_GMMBICthreshold$optim_k_summary #k9 17%, k10 16%, k8 13%, k11 10%, k7 9%, k2 7%, k6 7% etc
print(unname(round(system.time({
Pocillopora_SOM_hierarchicalDB <- clustering.SOM(Pocillopora_SOM_tr, #
                                                 max.k = 35,
                                                 clustering.method = "hierarchical+DB",
                                                 save.SOM.results.name = file.path(intermediate_files_folder, "Pocillopora_SOM_hierarchicalDB.Rdata"))
})[3] / 60, 1)))
Pocillopora_SOM_hierarchicalDB$optim_k_summary #k35 99%



## Plot and evaluate results
Pocillopora_SOM <- Pocillopora_SOM_kmeansBICelbow
plot.learning.SOM(Pocillopora_SOM)
plot.layer.distance.scale.SOM(Pocillopora_SOM)
plot.K.SOM(Pocillopora_SOM)
plot.model.SOM(Pocillopora_SOM, replicate.mode = "first")
plot.model.SOM(Pocillopora_SOM, replicate.mode = "representative")
plot.model.SOM(Pocillopora_SOM, replicate.mode = "representative", set.k = 3)
plot.model.SOM(Pocillopora_SOM, replicate.mode = "representative", set.k = 35)
plot.structure.SOM(Pocillopora_SOM, Individual.labels.font.size = 0.7)
print(unname(round(system.time({
Pocillopora_SOM_kmeansBICelbow_k3 <- clustering.SOM(Pocillopora_SOM, set.k = 3,
                                                    clustering.method = "kmeans+BICelbow")
})[3] / 60, 1)))
plot.structure.SOM(Pocillopora_SOM_kmeansBICelbow_k3, Individual.labels.font.size = 0.7)
plot.variable.importance.SOM(Pocillopora_SOM, left.margin = 9.5)
plot.layer.importance.varimp.SOM(Pocillopora_SOM)
plot.layer.importance.leaveoneout.SOM(Pocillopora_SOM, 
                                      bottom.margin = 9,
                                      save.leave.one.layer.out.results = TRUE,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Pocillopora_SOM_lolo.Rdata"))


## Compare results with species assignments of study
Pocillopora_SOM_ancestry_matrix <- as.data.frame(Pocillopora_SOM$ancestry_matrix)
Pocillopora_SOM_common_samples <- intersect(rownames(Pocillopora_SOM_ancestry_matrix), rownames(Pocillopora_species_names))
Pocillopora_SOM_ancestry_mat_sub <- Pocillopora_SOM_ancestry_matrix[Pocillopora_SOM_common_samples, , drop = FALSE]
Pocillopora_SOM_species_names_sub <- Pocillopora_species_names[Pocillopora_SOM_common_samples, , drop = FALSE]
Pocillopora_SOM_ancestry_matrix <- cbind(Pocillopora_SOM_ancestry_mat_sub, Pocillopora_SOM_species_names_sub) #compare ancestry matrix with species hypotheses
unique(sort(Pocillopora_SOM_ancestry_matrix$Primary_species_hypothesis))
length(unique(sort(Pocillopora_SOM_ancestry_matrix$Primary_species_hypothesis))) #11
length(unique(sort(Pocillopora_SOM_ancestry_matrix$Secondary_species_hypothesis))) ## 19
length(unique(sort(Pocillopora_SOM_ancestry_matrix$Genomic_species_hypothesis))) #20
table(Pocillopora_SOM_ancestry_matrix$Primary_species_hypothesis) #primary species hypotheses


## Assign majority cluster to each sample
majority_proportion <- 0.7
cluster.cols <- grep("^Cluster_", colnames(Pocillopora_SOM_ancestry_matrix))
cluster.mat <- as.matrix(Pocillopora_SOM_ancestry_matrix[, cluster.cols, drop = FALSE])
max.prop <- apply(cluster.mat, 1, max, na.rm = TRUE) #get highest ancestry proportion per sample
max.cluster <- apply(cluster.mat, 1, function(x) colnames(cluster.mat)[which.max(x)]) #get cluster with highest ancestry proportion
n.max <- apply(cluster.mat, 1, function(x) sum(x == max(x, na.rm = TRUE))) #check whether maximum is unique
Pocillopora_SOM_ancestry_matrix$Majority_cluster <- ifelse(max.prop >= majority_proportion & n.max == 1, max.cluster, "Mixed") #assign majority cluster only if >=N% and unique
table(Pocillopora_SOM_ancestry_matrix$Majority_cluster)
Pocillopora_SOM_ancestry_matrix[, c("Majority_cluster", "Primary_species_hypothesis")]


## Use k3
Pocillopora_SOM_ancestry_matrix_k3 <- as.data.frame(Pocillopora_SOM_kmeansBICelbow_k3$ancestry_matrix)
Pocillopora_SOM_common_samples_k3 <- intersect(rownames(Pocillopora_SOM_ancestry_matrix_k3), rownames(Pocillopora_species_names))
Pocillopora_SOM_ancestry_mat_sub_k3 <- Pocillopora_SOM_ancestry_matrix_k3[Pocillopora_SOM_common_samples_k3, , drop = FALSE]
Pocillopora_SOM_species_names_sub_k3 <- Pocillopora_species_names[Pocillopora_SOM_common_samples_k3, , drop = FALSE]
Pocillopora_SOM_ancestry_matrix_k3 <- cbind(Pocillopora_SOM_ancestry_mat_sub_k3, Pocillopora_SOM_species_names_sub_k3) #compare ancestry matrix with species hypotheses
unique(sort(Pocillopora_SOM_ancestry_matrix_k3$Primary_species_hypothesis))

cluster.cols.k3 <- grep("^Cluster_", colnames(Pocillopora_SOM_ancestry_matrix_k3))
cluster.mat.k3 <- as.matrix(Pocillopora_SOM_ancestry_matrix_k3[, cluster.cols.k3, drop = FALSE])
max.prop.k3 <- apply(cluster.mat.k3, 1, max, na.rm = TRUE) #get highest ancestry proportion per sample
max.cluster.k3 <- apply(cluster.mat.k3, 1, function(x) colnames(cluster.mat.k3)[which.max(x)]) #get cluster with highest ancestry proportion
n.max.k3 <- apply(cluster.mat.k3, 1, function(x) sum(x == max(x, na.rm = TRUE))) #check whether maximum is unique
Pocillopora_SOM_ancestry_matrix_k3$Majority_cluster_k3 <- ifelse(max.prop.k3 >= majority_proportion & n.max.k3 == 1, max.cluster.k3, "Mixed") #assign majority cluster only if >=50% and unique
table(Pocillopora_SOM_ancestry_matrix_k3$Majority_cluster_k3)
Pocillopora_SOM_ancestry_matrix_k3$Majority_cluster_k3 <- factor(Pocillopora_SOM_ancestry_matrix_k3$Majority_cluster_k3, levels = c(sort(unique(Pocillopora_SOM_ancestry_matrix_k3$Majority_cluster_k3[Pocillopora_SOM_ancestry_matrix_k3$Majority_cluster_k3 != "Mixed"])), "Mixed"))
Pocillopora_SOM_ancestry_matrix_k3[order(Pocillopora_SOM_ancestry_matrix_k3$Majority_cluster_k3), 
                                   c("Majority_cluster_k3", "Primary_species_hypothesis", 
                                     "Secondary_species_hypothesis", "Genomic_species_hypothesis")]


## Get unique PSH, SSH and GSH names for each cluster
unique_names_per_cluster_k3 <- lapply(split(Pocillopora_SOM_ancestry_matrix_k3, Pocillopora_SOM_ancestry_matrix_k3$Majority_cluster_k3), function(x) {
  list(
    Primary_species_hypothesis = unique(sort(x$Primary_species_hypothesis[!is.na(x$Primary_species_hypothesis) & x$Primary_species_hypothesis != ""])),
    Secondary_species_hypothesis = unique(sort(x$Secondary_species_hypothesis[!is.na(x$Secondary_species_hypothesis) & x$Secondary_species_hypothesis != ""])),
    Genomic_species_hypothesis = unique(sort(x$Genomic_species_hypothesis[!is.na(x$Genomic_species_hypothesis) & x$Genomic_species_hypothesis != ""])),
    Genomic_species_hypothesis_n = table(x$Genomic_species_hypothesis[!is.na(x$Genomic_species_hypothesis) & x$Genomic_species_hypothesis != ""])
  )
})
unique_names_per_cluster_k3 <- unique_names_per_cluster_k3[sapply(unique_names_per_cluster_k3, function(x) length(x$Primary_species_hypothesis) > 0 | length(x$Secondary_species_hypothesis) > 0 | length(x$Genomic_species_hypothesis) > 0)]
names(unique_names_per_cluster_k3)
for(cluster in names(unique_names_per_cluster_k3)) {
  cat("\n", cluster, "\n", sep = "")
  cat("PSH:", paste(unique_names_per_cluster_k3[[cluster]]$Primary_species_hypothesis, collapse = ", "), "\n")
  cat("SSH:", paste(unique_names_per_cluster_k3[[cluster]]$Secondary_species_hypothesis, collapse = ", "), "\n")
  cat("GSH:", paste(unique_names_per_cluster_k3[[cluster]]$Genomic_species_hypothesis, collapse = ", "), "\n")
  cat("GSH_n:", paste(names(unique_names_per_cluster_k3[[cluster]]$Genomic_species_hypothesis_n), unique_names_per_cluster_k3[[cluster]]$Genomic_species_hypothesis_n, sep = "=", collapse = ", "), "\n")
}


## Hierarchical analyses based on recovered clusters
Pocillopora_clusters_k3 <- Pocillopora_SOM_ancestry_matrix_k3$Majority_cluster_k3 #use assigned majority clusters
table(Pocillopora_clusters_k3)
Pocillopora_cluster_samples_k3 <- split(rownames(Pocillopora_SOM_ancestry_matrix_k3), Pocillopora_clusters_k3)
Pocillopora_cluster1_data_k3 <- lapply(Pocillopora_SOM_kmeansBICelbow_k3$input_data, function(x) x[Pocillopora_cluster_samples_k3$Cluster_1, , drop = FALSE]) #cluster 1 subset
Pocillopora_cluster2_data_k3 <- lapply(Pocillopora_SOM_kmeansBICelbow_k3$input_data, function(x) x[Pocillopora_cluster_samples_k3$Cluster_2, , drop = FALSE]) #cluster 2 subset
Pocillopora_cluster3_data_k3 <- lapply(Pocillopora_SOM_kmeansBICelbow_k3$input_data, function(x) x[Pocillopora_cluster_samples_k3$Cluster_3, , drop = FALSE]) #cluster 3 subset
print(unname(round(system.time({
SOM_Pocillopora_cluster1_k3 <- train.SOM(Pocillopora_cluster1_data_k3, 
                                         max.NA.row = 0.5, max.NA.col = 0.5, grid.multiplier = 4) #25 samples
})[3] / 60, 1)))
print(unname(round(system.time({
SOM_Pocillopora_cluster1_k3 <- clustering.SOM(SOM_Pocillopora_cluster1_k3, clustering.method = "kmeans+BICelbow", max.k = 15)
})[3] / 60, 1)))
SOM_Pocillopora_cluster1_k3$optim_k_summary #k1 100%
print(unname(round(system.time({
SOM_Pocillopora_cluster2_k3 <- train.SOM(Pocillopora_cluster2_data_k3, 
                                         max.NA.row = 0.5, max.NA.col = 0.5, grid.multiplier = 4) #26 samples
})[3] / 60, 1)))
print(unname(round(system.time({
SOM_Pocillopora_cluster2_k3 <- clustering.SOM(SOM_Pocillopora_cluster2_k3, clustering.method = "kmeans+BICelbow", max.k = 15)
})[3] / 60, 1)))
SOM_Pocillopora_cluster2_k3$optim_k_summary #k1 48%, k5 20%, k15 17%
print(unname(round(system.time({
SOM_Pocillopora_cluster3_k3 <- train.SOM(Pocillopora_cluster3_data_k3, 
                                         max.NA.row = 0.5, max.NA.col = 0.55, grid.multiplier = 4) #25 samples
})[3] / 60, 1)))
print(unname(round(system.time({
SOM_Pocillopora_cluster3_k3 <- clustering.SOM(SOM_Pocillopora_cluster3_k3, clustering.method = "kmeans+BICelbow", max.k = 15)
})[3] / 60, 1)))
SOM_Pocillopora_cluster3_k3$optim_k_summary #k1 99%


## Add GSH to sample names in ancestry matrix
Pocillopora_SOM_kmeansBICelbow_k3_updated <- Pocillopora_SOM_kmeansBICelbow_k3
rownames(Pocillopora_SOM_kmeansBICelbow_k3_updated$ancestry_matrix) <- paste0(rownames(Pocillopora_SOM_kmeansBICelbow_k3_updated$ancestry_matrix), "_", Pocillopora_SOM_ancestry_matrix_k3[rownames(Pocillopora_SOM_kmeansBICelbow_k3_updated$ancestry_matrix), "Genomic_species_hypothesis"])
plot.structure.SOM(Pocillopora_SOM_kmeansBICelbow_k3_updated, bottom.margin = 8.5, Individual.labels.font.size = 0.65)




#### Polygonia anglewing butterflies in Western Canada (Dupuis et al. 2018) ####

## https://doi.org/10.1093/zoolinnean/zlx081
## 4 inferred species
## Mitochondrial COI
## GBS SNPs
## Spatial data
## Updated environmental data
## Categorical/binary and continuous morphological data (wing color, categorical color scores and morphotype)

## Import and process genetic SNP data
Polygonia_SNP <- process.SNP.data.SOM(vcf.path = "../Empirical_examples/Dupuis_et_al_2018/Polygonia_961SNPs.vcf", #filter loci and individuals and create SNP matrix dataframe
                                      missing.loci.cutoff.lenient = 0.7,
                                      missing.loci.cutoff.final = 0.5,
                                      missing.individuals.cutoff = 0.5)
rownames(Polygonia_SNP) <- sub(".*?(\\d+)$", "\\1", rownames(Polygonia_SNP)) #only keep numeric identifier as rownames
ncol(Polygonia_SNP) #number of loci: 961
nrow(Polygonia_SNP) #number of samples: 237


## Import and filter COI data
Polygonia_COI <- process.SNP.data.SOM(nexus.path = "../Empirical_examples/Dupuis_et_al_2018/Polygonia_COI.nex",
                                      missing.loci.cutoff.lenient = 0.7, 
                                      missing.loci.cutoff.final = 0.5,
                                      missing.individuals.cutoff = 0.5)
Polygonia_COI_numeric_rownames <- sub(".*?(\\d+)$", "\\1", rownames(Polygonia_COI)) #extract numeric code from each rowname (e.g., "pf_8301" -> "8301")
Polygonia_COI <- Polygonia_COI[!duplicated(Polygonia_COI_numeric_rownames), , drop = FALSE] #keep only first occurrence for each numeric code (remove duplicates)
rownames(Polygonia_COI) <- Polygonia_COI_numeric_rownames[!duplicated(Polygonia_COI_numeric_rownames)] #set rownames to unique numeric codes
ncol(Polygonia_COI) #number of loci: 213
nrow(Polygonia_COI) #number of samples: 255


## Import and process RGB values
Polygonia_RGB <- read.delim("../Empirical_examples/Dupuis_et_al_2018/Polygonia_RGB_characters.txt", stringsAsFactors = FALSE)
rownames(Polygonia_RGB) <- Polygonia_RGB$Species
Polygonia_RGB <- magrittr::`%>%`(Polygonia_RGB, dplyr::select(-Name, -Species)) #remove columns


## Import and process wing character scores
Polygonia_wing_scores <- read.delim("../Empirical_examples/Dupuis_et_al_2018/Polygonia_visually_scored.txt", stringsAsFactors = FALSE)
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


## Import and process meta data with spatial data, species names and morphotype
Polygonia_metadata <- read.csv("../Empirical_examples/Dupuis_et_al_2018/Polygonia_metadata.csv", header = TRUE, sep = ";")
rownames(Polygonia_metadata) <- Polygonia_metadata$ID
nrow(Polygonia_metadata) #number of samples: 265

Polygonia_spatial <- dplyr::select(Polygonia_metadata, Latitude, Longitude) #create dataframe with Lat and Long
Polygonia_spatial$Elevation <- NA #initialize elevation column with NA
Polygonia_spatial_sf <- sf::st_as_sf(Polygonia_spatial[!is.na(Polygonia_spatial$Latitude) & !is.na(Polygonia_spatial$Longitude), ], coords = c("Longitude", "Latitude"), crs = 4326) #extract elevation
Polygonia_spatial$Elevation[!is.na(Polygonia_spatial$Latitude) & !is.na(Polygonia_spatial$Longitude)] <-
  elevatr::get_elev_point(locations = Polygonia_spatial_sf, prj = sf::st_crs(Polygonia_spatial_sf)$proj4string, src = "aws")$elevation

Polygonia_morphotype <- dplyr::select(Polygonia_metadata, Morphotype)
Polygonia_metadata <- dplyr::select(Polygonia_metadata, Species, ID)


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
Polygonia_morphology$Morphotype
Polygonia_morphology <- make.cols.binary.SOM(Polygonia_morphology, #convert Morphotype to binary columns and remove original
                                             make.binary.cols = "Morphotype",
                                             append.to.original = TRUE)
non.continuous.cols <- grepl("^Wing_character_|^Morphotype_", colnames(Polygonia_morphology)) #identify non-continuous traits
Polygonia_morphology_categorical <- Polygonia_morphology[, non.continuous.cols, drop = FALSE] #non-continuous traits
Polygonia_morphology <- Polygonia_morphology[, !non.continuous.cols, drop = FALSE] #continuous traits
Polygonia_morphology <- remove.lowCV.multicollinearity.SOM(Polygonia_morphology, #remove highly correlated and low-variance variables
                                                           CV.threshold = 0.05,
                                                           cor.threshold = 0.9)
ncol(Polygonia_morphology) #number of variables: 15
nrow(Polygonia_morphology) #number of samples: 217

ncol(Polygonia_morphology_categorical) #number of variables 15
nrow(Polygonia_morphology_categorical) #number of samples: 217


## Import and process environmental dataset (variables extracted and processed by separate R script based on coordinates)
Polygonia_environmental <- read.csv("../Empirical_examples/Dupuis_et_al_2018/Polygonia_environmental.csv",
                                    row.names = 1, header = TRUE)
Polygonia_environmental <- dplyr::select(Polygonia_environmental, -Latitude, -Longitude, -Elevation)
Polygonia_environmental_rownames <- rownames(Polygonia_environmental) #save rownames
Polygonia_environmental <- as.data.frame(lapply(Polygonia_environmental, as.numeric)) #ensure all columns are numeric
rownames(Polygonia_environmental) <- Polygonia_environmental_rownames #reassign saved row names
Polygonia_environmental <- (NicheDiv::transform.skewed.variables(Polygonia_environmental))$transformed #transform skewed variables
Polygonia_environmental <- remove.lowCV.multicollinearity.SOM(Polygonia_environmental, #remove highly correlated and low-variance variables
                                                              CV.threshold = 0.05,
                                                              cor.threshold = 0.9)
ncol(Polygonia_environmental) #number of variables: 125
nrow(Polygonia_environmental) #number of samples: 265


## Match datasets and remove all NA rows
for (Polygonia_shared_data in c("Polygonia_morphology", "Polygonia_morphology_categorical", "Polygonia_SNP", "Polygonia_COI",
                                "Polygonia_spatial", "Polygonia_environmental", "Polygonia_metadata")) {
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
nrow(Polygonia_metadata) #number of shared samples: 200


## Update rownames by adding species names to ID
Polygonia_species_vec <- Polygonia_metadata[rownames(Polygonia_morphology), "Species"]
Polygonia_new_rownames <- paste(rownames(Polygonia_morphology), Polygonia_species_vec, sep = "_")
rownames(Polygonia_morphology) <- Polygonia_new_rownames
rownames(Polygonia_morphology_categorical) <- Polygonia_new_rownames
rownames(Polygonia_SNP) <- Polygonia_new_rownames
rownames(Polygonia_COI) <- Polygonia_new_rownames
rownames(Polygonia_spatial) <- Polygonia_new_rownames
rownames(Polygonia_metadata) <- Polygonia_new_rownames
rownames(Polygonia_environmental) <- Polygonia_new_rownames



## Train and cluster SOM
Polygonia_all_data <- list(Morphology = Polygonia_morphology,
                           Morphology_2 = Polygonia_morphology_categorical,
                           SNP = Polygonia_SNP,
                           COI = Polygonia_COI,
                           Environmental = Polygonia_environmental,
                           Spatial = Polygonia_spatial)
print(unname(round(system.time({
Polygonia_SOM_tr <- train.SOM(input_data = Polygonia_all_data, #200 samples
                              save.SOM.results = TRUE,
                              save.SOM.results.name = file.path(intermediate_files_folder, "Polygonia_SOM_tr.Rdata"),
                              max.NA.row = 0.5,
                              max.NA.col = 0.5)
})[3] / 60, 1)))

print(unname(round(system.time({
Polygonia_SOM_kmeansBICthreshold <- clustering.SOM(Polygonia_SOM_tr, #
                                                   clustering.method = "kmeans+BICthreshold",
                                                   save.SOM.results.name = file.path(intermediate_files_folder, "Polygonia_SOM_kmeansBICthreshold.Rdata"))
})[3] / 60, 1)))
Polygonia_SOM_kmeansBICthreshold$optim_k_summary #k5 46%, k6 36%, k7 15%
print(unname(round(system.time({
Polygonia_SOM_HDBSCAN <- clustering.SOM(Polygonia_SOM_tr, #
                                        clustering.method = "HDBSCAN",
                                        save.SOM.results.name = file.path(intermediate_files_folder, "Polygonia_SOM_HDBSCAN.Rdata"))
})[3] / 60, 1)))
Polygonia_SOM_HDBSCAN$optim_k_summary #k3 82%, k4 10%, k2 8%
print(unname(round(system.time({
Polygonia_SOM_hierarchicalDB <- clustering.SOM(Polygonia_SOM_tr, #
                                               clustering.method = "hierarchical+DB",
                                               save.SOM.results.name = file.path(intermediate_files_folder, "Polygonia_SOM_hierarchicalDB.Rdata"))
})[3] / 60, 1)))
Polygonia_SOM_hierarchicalDB$optim_k_summary #k3 76%, k4 10%, k5 6%
print(unname(round(system.time({
Polygonia_SOM_GMMBICthreshold <- clustering.SOM(Polygonia_SOM_tr, #ca. 15min
                                                clustering.method = "GMM+BICthreshold",
                                                save.SOM.results.name = file.path(intermediate_files_folder, "Polygonia_SOM_GMMBICthreshold.Rdata"))
})[3] / 60, 1)))
Polygonia_SOM_GMMBICthreshold$optim_k_summary #k4 32%, k2 28%, k5 21%, k6 10%
print(unname(round(system.time({
Polygonia_SOM_OPTICSSilhouette <- clustering.SOM(Polygonia_SOM_tr, #ca 5min
                                                 clustering.method = "OPTICS+Silhouette",
                                                 save.SOM.results.name = file.path(intermediate_files_folder, "Polygonia_SOM_OPTICSSilhouette.Rdata"))
})[3] / 60, 1)))
Polygonia_SOM_OPTICSSilhouette$optim_k_summary #k3 94%, k4 6%
print(unname(round(system.time({
Polygonia_SOM_kmeansBICelbow <- clustering.SOM(Polygonia_SOM_tr, #ca 3min
                                               clustering.method = "kmeans+BICelbow",
                                               save.SOM.results.name = file.path(intermediate_files_folder, "Polygonia_SOM_kmeansBICelbow.Rdata"))
})[3] / 60, 1)))
Polygonia_SOM_kmeansBICelbow$optim_k_summary #k3 83%, k4 16%


## Evaluate and plot results
Polygonia_SOM <- Polygonia_SOM_kmeansBICelbow
plot.K.SOM(Polygonia_SOM)
plot.learning.SOM(Polygonia_SOM)
plot.layer.distance.scale.SOM(Polygonia_SOM)
plot.model.SOM(Polygonia_SOM, replicate.mode = "representative")
plot.model.SOM(Polygonia_SOM, replicate.mode = "representative", set.k = 3)
plot.model.SOM(Polygonia_SOM, replicate.mode = "representative", set.k = 4)
plot.model.SOM(Polygonia_SOM, replicate.mode = "first")
plot.map.SOM(Polygonia_SOM,
             Coordinates = Polygonia_spatial[, c(1:2)],
             lat.buffer.range = 1,
             lon.buffer.range = 2,
             north.arrow.position = c(0.01, 0.92), #position (x, y) of north arrow relative to map
             north.arrow.length = 0.7, #length of north arrow
             north.arrow.N.position = 0.3, #position of north arrow "N"
             north.arrow.N.size = 1, #size of north arrow "N"
             scale.position = c(0.78, 0.02)) #relative position (x, y) of scale
plot.structure.SOM(Polygonia_SOM, Individual.labels.font.size = 0.25)


## Evaluate layer importance
plot.layer.importance.varimp.SOM(Polygonia_SOM, bottom.margin = 6.5)
plot.layer.importance.leaveoneout.SOM(Polygonia_SOM, 
                                      bottom.margin = 9,
                                      save.leave.one.layer.out.results = TRUE,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Poligonia_SOM_lolo.Rdata"))


## Evaluate variable importance
plot.variable.importance.SOM(Polygonia_SOM, mode = "Cluster.separation")
plot.variable.importance.SOM(Polygonia_SOM, mode = "Map.variance")
sort(roundPolygonia_SOM$median_etasquared_variable_importance[[1]], decreasing = T)
head(round(sort(Polygonia_SOM$median_etasquared_variable_importance[[1]], decreasing = T), 2), 10)
head(round(sort(Polygonia_SOM$median_etasquared_variable_importance[[2]], decreasing = T), 2), 20)
head(round(sort(Polygonia_SOM$median_etasquared_variable_importance[[3]], decreasing = T), 2), 500)
head(round(sort(Polygonia_SOM$median_etasquared_variable_importance[[4]], decreasing = T), 2), 50)
head(round(sort(Polygonia_SOM$median_etasquared_variable_importance[[5]], decreasing = T), 2), 10)

head(round(sort(Polygonia_SOM$median_map_variance_variable_importance[[1]], decreasing = T), 2), 10)
head(round(sort(Polygonia_SOM$median_map_variance_variable_importance[[2]], decreasing = T), 2), 20)
head(round(sort(Polygonia_SOM$median_map_variance_variable_importance[[3]], decreasing = T), 2), 500)
head(round(sort(Polygonia_SOM$median_map_variance_variable_importance[[4]], decreasing = T), 2), 50)
head(round(sort(Polygonia_SOM$median_map_variance_variable_importance[[5]], decreasing = T), 2), 10)

compare.top.overlap <- function(etasq.list, mapvar.list, top.n) { #check for overlap in both metrices
  if(length(etasq.list) != length(mapvar.list)) stop("'etasq.list' and 'mapvar.list' must have the same number of layers")
  if(length(top.n) == 1) top.n <- rep(top.n, length(etasq.list))
  if(length(top.n) != length(etasq.list)) stop(paste0("'top.n' must have length 1 or one value per layer. Found ", length(top.n), " value(s) for ", length(etasq.list), " layer(s)."))
  lapply(seq_along(etasq.list), function(i) intersect(names(sort(etasq.list[[i]], decreasing = T))[1:min(top.n[i], length(etasq.list[[i]]))], names(sort(mapvar.list[[i]], decreasing = T))[1:min(top.n[i], length(mapvar.list[[i]]))]))
}
top.overlap <- compare.top.overlap(Polygonia_SOM$median_etasquared_variable_importance, Polygonia_SOM$median_map_variance_variable_importance, c(5, 5, 50, 20, 5, 1))
for(i in seq_along(top.overlap)) cat("Layer ", i, " overlap: ", if(length(top.overlap[[i]]) == 0) "none" else paste(top.overlap[[i]], collapse = ", "), "\n", sep = "")


## Evaluate species IDs from study
Polygonia_ancestry <- as.data.frame(Polygonia_SOM$ancestry_matrix)
Polygonia_ancestry$ID <- rownames(Polygonia_ancestry)
Polygonia_metadata$ID <- as.character(rownames(Polygonia_metadata))
Polygonia_ancestry <- merge(Polygonia_metadata, Polygonia_ancestry, by = "ID") #examine ancestry matrix and species labels
table(Polygonia_ancestry$Species)



## Hierarchical analyses based on recovered clusters
Polygonia_clusters <- apply(Polygonia_SOM$ancestry_matrix, 1, which.max) #assign each sample to cluster with highest ancestry proportion
Polygonia_clusters <- paste0("cluster", Polygonia_clusters) #rename clusters
table(Polygonia_clusters)
Polygonia_cluster_samples <- split(rownames(Polygonia_SOM$ancestry_matrix), Polygonia_clusters)
Polygonia_cluster1_data <- lapply(Polygonia_SOM$input_data, function(x) x[Polygonia_cluster_samples$cluster1, , drop = FALSE]) #cluster 1 subset
Polygonia_cluster2_data <- lapply(Polygonia_SOM$input_data, function(x) x[Polygonia_cluster_samples$cluster2, , drop = FALSE]) #cluster 2 subset
Polygonia_cluster3_data <- lapply(Polygonia_SOM$input_data, function(x) x[Polygonia_cluster_samples$cluster3, , drop = FALSE]) #cluster 3 subset

print(unname(round(system.time({
Polygonia_SOM_tr_cluster1 <- train.SOM(Polygonia_cluster1_data, #75 samples
                                       max.NA.row = 0.5,
                                       max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Polygonia_SOM_cluster1 <- clustering.SOM(Polygonia_SOM_tr_cluster1,
                                         clustering.method = "kmeans+BICelbow")
})[3] / 60, 1)))
Polygonia_SOM_cluster1$optim_k_summary #k1 100%

print(unname(round(system.time({
Polygonia_SOM_tr_cluster2 <- train.SOM(Polygonia_cluster2_data,
                                       max.NA.row = 0.5,
                                       max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Polygonia_SOM_cluster2 <- clustering.SOM(Polygonia_SOM_tr_cluster2, #39 samples
                                         clustering.method = "kmeans+BICelbow")
})[3] / 60, 1)))
Polygonia_SOM_cluster2$optim_k_summary #k1 100%

print(unname(round(system.time({
Polygonia_SOM_tr_cluster3 <- train.SOM(Polygonia_cluster3_data, #72 samples
                                       max.NA.row = 0.5,
                                       max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Polygonia_SOM_cluster3 <- clustering.SOM(Polygonia_SOM_tr_cluster3,
                                         clustering.method = "kmeans+BICelbow")
})[3] / 60, 1)))
Polygonia_SOM_cluster3$optim_k_summary #k2 100%


plot.model.SOM(Polygonia_SOM_cluster3, replicate.mode = "representative")
plot.model.SOM(Polygonia_SOM_cluster3, replicate.mode = "first")
plot.structure.SOM(Polygonia_SOM_cluster3, bottom.margin = 8)
plot.K.SOM(Polygonia_SOM_cluster3)
plot.map.SOM(SOM.output = Polygonia_SOM_cluster3,
             Coordinates = Polygonia_spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 5, #add coordinates as buffer range around latitude coordinates
             lon.buffer.range = 5, #add coordinates as buffer range around longitude coordinates
             pie.size = 2.5, #pie chart size
             north.arrow.position = c(0.04, 0.9), #position (x, y) of north arrow relative to map
             north.arrow.length = 1, #length of north arrow
             north.arrow.N.position = 0.3, #position of north arrow "N"
             north.arrow.N.size = 1) #size of north arrow "N"
plot.variable.importance.SOM(Polygonia_SOM_cluster3,
                             mode = "Cluster.separation", 
                             left.margin = 6.5,
                             bar.label.font.size = 0.4)
plot.variable.importance.SOM(Polygonia_SOM_cluster3, 
                             mode = "Map.variance", 
                             left.margin = 6.5,
                             bar.label.font.size = 0.4)
plot.layer.importance.varimp.SOM(Polygonia_SOM_cluster3, bottom.margin = 6)
plot.layer.importance.leaveoneout.SOM(Polygonia_SOM_cluster3, 
                                      save.leave.one.layer.out.results = TRUE,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Polygonia_SOM_cluster3_lolo.Rdata")) #this will take 10-20min (running 2 x N replicates for train and clustering SOM)

Polygonia_ancestry_SOM_cluster3 <- as.data.frame(Polygonia_SOM_cluster3$ancestry_matrix)
Polygonia_ancestry_SOM_cluster3$Species <- Polygonia_metadata$Species[match(rownames(Polygonia_SOM_cluster3$ancestry_matrix), rownames(Polygonia_metadata))]
Polygonia_ancestry_SOM_cluster3$Species_revised <- Polygonia_metadata$Species[match(rownames(Polygonia_SOM_cluster3$ancestry_matrix), rownames(Polygonia_metadata))]
length(unique(Polygonia_ancestry_SOM_cluster3$Species)) #number of species present in data
length(unique(Polygonia_ancestry_SOM_cluster3$Species_revised)) #number of proposed species present in data
table(Polygonia_ancestry_SOM_cluster3$Species)
table(Polygonia_ancestry_SOM_cluster3$Species_revised)


## Calculate pairwise Weir and Cockerham Fst among species
SNP.ids <- rownames(Polygonia_SNP)
metadata.ids <- as.character(Polygonia_metadata$ID)
overlapping.ids <- intersect(SNP.ids, metadata.ids)
if (length(overlapping.ids) < 2) stop("Fewer than two overlapping IDs were found between Polygonia_SNP and Polygonia_metadata")
Polygonia_SNP <- Polygonia_SNP[overlapping.ids, , drop = FALSE]
Polygonia_metadata <- Polygonia_metadata[match(overlapping.ids, metadata.ids), , drop = FALSE]
rownames(Polygonia_metadata) <- overlapping.ids
Polygonia_metadata$Species <- as.character(Polygonia_metadata$Species)
Polygonia_metadata$Species[Polygonia_metadata$Species %in% c("Polygonia gracilis gracilis", "Polygonia gracilis zephyrus")] <- "Polygonia gracilis"
species.rows <- !is.na(Polygonia_metadata$Species)
species.rows <- species.rows & !Polygonia_metadata$Species %in% c("Polygonia comma", "Polygonia interrogationis")
Polygonia_SNP <- Polygonia_SNP[species.rows, , drop = FALSE]
Polygonia_metadata <- Polygonia_metadata[species.rows, , drop = FALSE]
SNP.matrix <- as.matrix(Polygonia_SNP)
storage.mode(SNP.matrix) <- "numeric"
SNP.values <- unique(as.vector(SNP.matrix))
SNP.values <- SNP.values[!is.na(SNP.values)]
if (!all(SNP.values %in% c(0, 1, 2))) stop("Polygonia_SNP must only contain 0, 1, 2, or NA values")
SNP.genotypes <- as.data.frame(apply(SNP.matrix, 2, function(x) {
  genotype.vector <- rep(NA_character_, length(x))
  genotype.vector[x == 0] <- "A/A"
  genotype.vector[x == 1] <- "A/B"
  genotype.vector[x == 2] <- "B/B"
  return(genotype.vector)
}))
rownames(SNP.genotypes) <- rownames(Polygonia_SNP)
Polygonia_genind <- adegenet::df2genind(SNP.genotypes, sep = "/", ncode = 1, ploidy = 2, pop = Polygonia_metadata$Species)
calculate.overall.Fst <- function(genind.object) {
  genotype.table <- genind.object@tab #extract allele-count table
  locus.factor <- genind.object@loc.fac #extract locus factor for allele columns
  population.vector <- as.character(adegenet::pop(genind.object)) #extract population assignments
  if (is.null(genotype.table) || is.null(locus.factor)) stop("genind.object does not contain genotype table or locus factor")
  if (length(population.vector) != nrow(genotype.table)) stop("Length mismatch between populations and genotype table rows")
  population.names <- sort(unique(population.vector[!is.na(population.vector)]))
  if (length(population.names) != 2) stop("Manual Weir and Cockerham Fst calculation expects exactly two populations")
  locus.names <- unique(as.character(locus.factor))
  a.vector <- rep(NA_real_, length(locus.names))
  b.vector <- rep(NA_real_, length(locus.names))
  c.vector <- rep(NA_real_, length(locus.names))
  names(a.vector) <- locus.names
  names(b.vector) <- locus.names
  names(c.vector) <- locus.names
  for (locus.index in seq_along(locus.names)) {
    current.locus <- locus.names[locus.index]
    locus.columns <- which(as.character(locus.factor) == current.locus)
    if (length(locus.columns) != 2) next
    locus.table <- genotype.table[, locus.columns, drop = FALSE]
    population.sample.sizes <- rep(NA_real_, length(population.names))
    population.alt.allele.frequencies <- rep(NA_real_, length(population.names))
    population.observed.heterozygosities <- rep(NA_real_, length(population.names))
    for (population.index in seq_along(population.names)) {
      population.rows <- which(population.vector == population.names[population.index])
      population.locus.table <- locus.table[population.rows, , drop = FALSE]
      complete.rows <- rowSums(is.na(population.locus.table)) == 0
      population.locus.table <- population.locus.table[complete.rows, , drop = FALSE]
      if (nrow(population.locus.table) == 0) next
      allele.copy.counts <- rowSums(population.locus.table)
      diploid.rows <- allele.copy.counts == 2
      population.locus.table <- population.locus.table[diploid.rows, , drop = FALSE]
      if (nrow(population.locus.table) == 0) next
      population.sample.sizes[population.index] <- nrow(population.locus.table)
      population.alt.allele.frequencies[population.index] <- sum(population.locus.table[, 2]) / (2 * nrow(population.locus.table))
      population.observed.heterozygosities[population.index] <- mean(population.locus.table[, 1] == 1 & population.locus.table[, 2] == 1)
    }
    if (any(is.na(population.sample.sizes)) || any(population.sample.sizes <= 1)) next
    n.populations <- length(population.names)
    total.sample.size <- sum(population.sample.sizes)
    n.bar <- mean(population.sample.sizes)
    n.c <- (total.sample.size - sum(population.sample.sizes^2) / total.sample.size) / (n.populations - 1)
    p.bar <- sum(population.sample.sizes * population.alt.allele.frequencies) / total.sample.size
    s.squared <- sum(population.sample.sizes * (population.alt.allele.frequencies - p.bar)^2) / ((n.populations - 1) * n.bar)
    h.bar <- sum(population.sample.sizes * population.observed.heterozygosities) / total.sample.size
    if (n.c <= 0 || p.bar <= 0 || p.bar >= 1) next
    a.vector[locus.index] <- (n.bar / n.c) * (s.squared - (1 / (n.bar - 1)) * ((p.bar * (1 - p.bar)) - ((n.populations - 1) / n.populations) * s.squared - 0.25 * h.bar))
    b.vector[locus.index] <- (n.bar / (n.bar - 1)) * ((p.bar * (1 - p.bar)) - ((n.populations - 1) / n.populations) * s.squared - ((2 * n.bar - 1) / (4 * n.bar)) * h.bar)
    c.vector[locus.index] <- 0.5 * h.bar
  }
  denominator <- sum(a.vector + b.vector + c.vector, na.rm = TRUE)
  if (is.na(denominator) || denominator == 0) return(NA_real_)
  overall.fst <- sum(a.vector, na.rm = TRUE) / denominator
  return(overall.fst)
}
calculate.pairwise.Fst <- function(genind.object) {
  population.vector <- as.character(adegenet::pop(genind.object)) #extract population assignments
  if (length(population.vector) != adegenet::nInd(genind.object)) stop("Length mismatch between populations and genind individuals")
  population.names <- sort(unique(population.vector[!is.na(population.vector)]))
  if (length(population.names) < 2) stop("At least two populations are required")
  population.pairs <- utils::combn(population.names, 2, simplify = FALSE)
  fst.matrix <- matrix(NA_real_, nrow = length(population.names), ncol = length(population.names))
  rownames(fst.matrix) <- population.names
  colnames(fst.matrix) <- population.names
  diag(fst.matrix) <- 0
  fst.table <- data.frame(Species_1 = character(), Species_2 = character(), Fst = numeric(), stringsAsFactors = FALSE)
  for (pair.index in seq_along(population.pairs)) {
    current.pair <- population.pairs[[pair.index]]
    pair.rows <- population.vector %in% current.pair
    pair.genind <- genind.object[pair.rows, ]
    pair.genind@pop <- factor(population.vector[pair.rows], levels = current.pair)
    current.fst <- calculate.overall.Fst(pair.genind)
    fst.matrix[current.pair[1], current.pair[2]] <- current.fst
    fst.matrix[current.pair[2], current.pair[1]] <- current.fst
    fst.table <- rbind(fst.table, data.frame(Species_1 = current.pair[1], Species_2 = current.pair[2], Fst = current.fst, stringsAsFactors = FALSE))
  }
  return(list(Fst_matrix = fst.matrix, Fst_table = fst.table))
}
Polygonia_pairwise_Fst <- calculate.pairwise.Fst(Polygonia_genind)
Polygonia_pairwise_Fst$Fst_matrix
Polygonia_pairwise_Fst$Fst_table




#### Viburnum woody shrubs in Eastern North America (Spriggs et al. 2018) ######

## https://doi.org/10.1093/sysbio/syy084
## SNP data based on RAD seq: 2-3k SNPs
## Morphology: 9 traits
## Two previous species (two clades), revised to three species
## XX shared samples of all three species
## Their niche modelling and coordinates data were not available for individuals


## Import and process genetic SNP data
Viburnum_SNP <- process.SNP.data.SOM(vcf.path = "../Empirical_examples/Spriggs_et_al_2018/nudum-c88-d6-min50.vcf.gz",
                                     missing.loci.cutoff.lenient = 0.7,
                                     missing.loci.cutoff.final = 0.5,
                                     missing.individuals.cutoff = 0.5)
ncol(Viburnum_SNP) #number of SNPs: 42159
nrow(Viburnum_SNP) #number of samples: 65



## Import and process morphological dataset
Viburnum_morphology <- read.delim("../Empirical_examples/Spriggs_et_al_2018/morphological_trait_data2.txt", stringsAsFactors = FALSE)
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
ncol(Viburnum_morphology) #number of traits: 8
nrow(Viburnum_morphology) #number of samples: 145


## Import and process metadata
Viburnum_metadata <- read.delim("../Empirical_examples/Spriggs_et_al_2018/morphological_trait_data2.txt", stringsAsFactors = FALSE)
Viburnum_metadata <- Viburnum_metadata[!duplicated(Viburnum_metadata$Individual), ] #remove duplicate IDs
rownames(Viburnum_metadata) <- Viburnum_metadata$Individual #add rownames
Viburnum_metadata <- Viburnum_metadata[, c("State", "County"), drop = FALSE] #only keep State and County columns
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
Viburnum_SOM_data <- list(Morphology = Viburnum_morphology, 
                          SNP = Viburnum_SNP)
print(unname(round(system.time({
Viburnum_SOM_tr <- train.SOM(Viburnum_SOM_data, #46 samples
                             max.NA.row = 0.5,
                             max.NA.col = 0.5,
                             save.SOM.results.name = file.path(intermediate_files_folder, "Viburnum_SOM_tr.Rdata"),
                             save.SOM.results = T)
})[3] / 60, 1)))

print(unname(round(system.time({
Viburnum_SOM_kmeansBICthreshold <- clustering.SOM(Viburnum_SOM_tr, #takes ca 2min!
                                                  clustering.method = "kmeans+BICthreshold",
                                                  save.SOM.results.name = file.path(intermediate_files_folder, "Viburnum_SOM_kmeansBICthreshold.Rdata"))
})[3] / 60, 1)))
Viburnum_SOM_kmeansBICthreshold$optim_k_summary #k2 100%
print(unname(round(system.time({
Viburnum_SOM_HDBSCAN <- clustering.SOM(Viburnum_SOM_tr, #takes ca 2min!
                                       clustering.method = "HDBSCAN",
                                       save.SOM.results.name = file.path(intermediate_files_folder, "Viburnum_SOM_HDBSCAN.Rdata"))
})[3] / 60, 1)))
Viburnum_SOM_HDBSCAN$optim_k_summary #k2 97%
print(unname(round(system.time({
Viburnum_SOM_hierarchicalDB <- clustering.SOM(Viburnum_SOM_tr, #takes ca 45min!
                                              clustering.method = "hierarchical+DB",
                                              save.SOM.results.name = file.path(intermediate_files_folder, "Viburnum_SOM_hierarchicalDB.Rdata"))
})[3] / 60, 1)))
Viburnum_SOM_hierarchicalDB$optim_k_summary #k2 100%
print(unname(round(system.time({
Viburnum_SOM_GMMBICthreshold <- clustering.SOM(Viburnum_SOM_tr, #ca. 15min
                                               clustering.method = "GMM+BICthreshold",
                                               save.SOM.results.name = file.path(intermediate_files_folder, "Viburnum_SOM_GMMBICthreshold.Rdata"))
})[3] / 60, 1)))
Viburnum_SOM_GMMBICthreshold$optim_k_summary #k2 52%, k3 44%
print(unname(round(system.time({
Viburnum_SOM_OPTICSSilhouette <- clustering.SOM(Viburnum_SOM_tr, #ca 5min
                                                clustering.method = "OPTICS+Silhouette",
                                                save.SOM.results.name = file.path(intermediate_files_folder, "Viburnum_SOM_OPTICSSilhouette.Rdata"))
})[3] / 60, 1)))
Viburnum_SOM_OPTICSSilhouette$optim_k_summary #k2 100%
print(unname(round(system.time({
Viburnum_SOM_kmeansBICelbow <- clustering.SOM(Viburnum_SOM_tr, #ca 3min
                                              clustering.method = "kmeans+BICelbow",
                                              save.SOM.results.name = file.path(intermediate_files_folder, "Viburnum_SOM_kmeansBICelbow.Rdata"))
})[3] / 60, 1)))
Viburnum_SOM_kmeansBICelbow$optim_k_summary #k2 100%



## Evaluate and plot results
Viburnum_SOM <- Viburnum_SOM_kmeansBICelbow
plot.learning.SOM(Viburnum_SOM)
plot.layer.distance.scale.SOM(Viburnum_SOM)
plot.K.SOM(Viburnum_SOM)
plot.model.SOM(Viburnum_SOM)
plot.model.SOM(Viburnum_SOM, replicate.mode = "first")
plot.variable.importance.SOM(Viburnum_SOM,mode = "Cluster.separation",
                             left.margin = 5.5,
                             bar.label.font.size = 0.4)
plot.variable.importance.SOM(Viburnum_SOM,mode = "Map.variance",
                             left.margin = 5.5,
                             bar.label.font.size = 0.4)
plot.structure.SOM(Viburnum_SOM, Individual.labels.font.size = 0.8)
plot.layer.importance.leaveoneout.SOM(Viburnum_SOM, 
                                      bottom.margin = 6.5,
                                      save.leave.one.layer.out.results = TRUE,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Viburnum_SOM_lolo.Rdata"))
plot.layer.importance.varimp.SOM(Viburnum_SOM, bottom.margin = 5.5)

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
Viburnum_ancestry <- tidyr::separate(Viburnum_ancestry,
                                     Species_State,
                                     into = c("Species", "State"),
                                     sep = ", ",
                                     remove = TRUE)
table(Viburnum_ancestry$Species)


## Hierarchical analyses based on recovered clusters
Viburnum_clusters <- apply(Viburnum_SOM$ancestry_matrix, 1, which.max) #assign each sample to cluster with highest ancestry proportion
Viburnum_clusters <- paste0("cluster", Viburnum_clusters) #rename clusters
table(Viburnum_clusters)
Viburnum_cluster_samples <- split(rownames(Viburnum_SOM$ancestry_matrix), Viburnum_clusters)
Viburnum_cluster1_data <- lapply(Viburnum_SOM$input_data, function(x) x[Viburnum_cluster_samples$cluster1, , drop = FALSE]) #cluster 1 subset
Viburnum_cluster2_data <- lapply(Viburnum_SOM$input_data, function(x) x[Viburnum_cluster_samples$cluster2, , drop = FALSE]) #cluster 2 subset

print(unname(round(system.time({
Viburnum_SOM_tr_cluster1 <- train.SOM(Viburnum_cluster1_data, #25 samples
                                      grid.multiplier = 4,
                                      max.NA.row = 0.5,
                                      max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Viburnum_SOM_cluster1 <- clustering.SOM(Viburnum_SOM_tr_cluster1,
                                        clustering.method = "kmeans+BICelbow")
})[3] / 60, 1)))
Viburnum_SOM_cluster1$optim_k_summary #k1 100%

print(unname(round(system.time({
Viburnum_SOM_tr_cluster2 <- train.SOM(Viburnum_cluster2_data, #21 samples
                                      grid.multiplier = 4,
                                      max.NA.row = 0.5,
                                      max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Viburnum_SOM_cluster2 <- clustering.SOM(Viburnum_SOM_tr_cluster2,
                                        clustering.method = "kmeans+BICelbow")
})[3] / 60, 1)))
Viburnum_SOM_cluster2$optim_k_summary #k1 100%




#### Microcebus lemurs dataset - Madagaskar (van Elst et al 2024) #######
library(dplyr)

## Import and process genetic SNP data
Microcebus_SNP <- process.SNP.data.SOM(
  vcf.path = "../Empirical_examples/van_Elst_et_al_2024/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss0.05.thinned.vcf.gz", #VCF file path
  missing.loci.cutoff.lenient = 0.7,
  missing.loci.cutoff.final = 0.5,
  missing.individuals.cutoff = 0.5)
Microcebus_species_split <- stringr::str_split_fixed(rownames(Microcebus_SNP), "_", n = 2) #split rownames by underscore
rownames(Microcebus_SNP) <- Microcebus_species_split[, 2] #set rownames to just the ID part
ncol(Microcebus_SNP) #number of SNPs: 9429
nrow(Microcebus_SNP) #number of samples: 213


## Import and process multiple data dataset 2 containing range of data types
Microcebus_multiple_data2 <- utils::read.csv("../Empirical_examples/van_Elst_et_al_2024/01_Microcebus_morphological_data.csv", 
                                             stringsAsFactors = FALSE, header = TRUE, sep = ";")
Microcebus_multiple_data2 <- Microcebus_multiple_data2 %>% #only keep individuals that are Rad sequenced (have SNP data)
  dplyr::filter(RADSeq.available != "no" & !is.na(RADSeq.available))
rownames(Microcebus_multiple_data2) <- Microcebus_multiple_data2$Individual.ID


## Import and process multiple data dataset containing range of data types
Microcebus_multiple_data <- utils::read.csv("../Empirical_examples/van_Elst_et_al_2024/data.csv", 
                                            stringsAsFactors = FALSE, header = TRUE, sep = ";")
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
Microcebus_duplicates_new_ids <- stringr::str_replace(Microcebus_duplicates_original_ids, "^(\\d+)y(\\d+)_(.+)$", "\\1-\\2\\3") #modify rownames
Microcebus_duplicates_new_ids <- stringr::str_replace(Microcebus_duplicates_new_ids, "^(.*?)y(\\d+)$", "\\1-\\2") #modify rownames
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
Microcebus_multiple_data_combined <- Microcebus_multiple_data_combined[rownames(Microcebus_SNP), , drop = FALSE] #re-order and pad to SNP samples
Microcebus_multiple_data_combined <- Microcebus_multiple_data_combined %>% dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.))) #remove rows with all NA
Microcebus_SNP <- Microcebus_SNP[rownames(Microcebus_multiple_data_combined), , drop = FALSE] #subset SNP dataset to combined samples


## Add Reproductive_period_month column (1–12) for reproductively active animals (based on Supplementary Methods information)
Microcebus_multiple_data_combined <- Microcebus_multiple_data_combined %>%
  dplyr::mutate(MonthNum = as.integer(stringr::str_extract(month, "^\\d{2}")), #extract numeric month (01–12)
                Reproductive_period_month = dplyr::case_when(
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
  dplyr::rename(Interorbital_Distance = interorbital.dist, #modify column names
                Intraorbital_Distance = intraorbital.dist) %>%
  dplyr::rename_with(~ stringr::str_replace_all(.x, "\\.", "_"), dplyr::everything()) %>% #rename dots to underscores
  dplyr::rename_with(~ stringr::str_to_title(.x), dplyr::everything()) %>% #capitalize first letter
  dplyr::mutate(dplyr::across(dplyr::where(is.character), ~ dplyr::na_if(., "#NV")),
                Testis_width_total = dplyr::case_when(Testis_width_total %in% c("0", "large", "small") ~ NA_character_, TRUE ~ Testis_width_total)) %>% #replace “#NV”, 0, "large" and "small" with NA before numeric conversion
  dplyr::mutate(dplyr::across(dplyr::where(is.character), ~ sub(",", ".", .))) %>% #change comma to decimal point
  dplyr::mutate(dplyr::across(dplyr::where(is.character), as.numeric)) #convert remaining chars to numeric
Microcebus_log_transform_traits <- c(
  "Ear_length",
  "Ear_width",
  "Head_length",
  "Head_width",
  "Interorbital_distance",
  "Intraorbital_distance",
  "Snout_length",
  "Lower_leg_length",
  "Hind_foot_length",
  "Third_toe_length",
  "Body_length",
  "Tail_length",
  "Body_mass",
  "Tail_circumference")
Microcebus_morphology[, Microcebus_log_transform_traits] <- lapply(Microcebus_morphology[, Microcebus_log_transform_traits, drop = FALSE], function(trait_values) { #log transform positive somatic size traits
  trait_values <- as.numeric(trait_values)
  if(any(trait_values <= 0, na.rm = TRUE)) trait_values
  else log(trait_values)
}) #log transform positive somatic size traits
Microcebus_row_names <- rownames(Microcebus_morphology) #store original rownames
Microcebus_body_size <- Microcebus_morphology[, "Body_mass"] #extract log Body_mass and residualize other log-transformed somatic traits
Microcebus_somatic_shape_trait_names <- Microcebus_log_transform_traits[Microcebus_log_transform_traits != "Body_mass"] #select log-transformed somatic traits to residualize
Microcebus_somatic_residuals_mat <- sapply(Microcebus_somatic_shape_trait_names, function(trait_name) {stats::resid(stats::lm(Microcebus_morphology[, trait_name] ~ Microcebus_body_size, na.action = stats::na.exclude))}) #regress each log-transformed somatic trait on log Body_length and store residuals
Microcebus_reproductive_trait_names <- c("Testis_width_total",
                                         "Testis_width_left",
                                         "Testis_width_right",
                                         "Testis_length_left",
                                         "Testis_length_right",
                                         "Reproductive_period_month") #traits not included in body-size correction
Microcebus_morphology <- as.data.frame(cbind(Body_length = Microcebus_body_size,
                                             Microcebus_somatic_residuals_mat,
                                             Microcebus_morphology[, Microcebus_reproductive_trait_names, drop = FALSE])) #combine log Body_length, size-corrected somatic traits and raw reproductive traits
rownames(Microcebus_morphology) <- Microcebus_row_names #restore original rownames
Microcebus_morphology <- remove.lowCV.multicollinearity.SOM(Microcebus_morphology, #remove highly correlated and low-variance variables
                                                            CV.threshold = 0.05,
                                                            cor.threshold = 0.9)
ncol(Microcebus_morphology) #number of  traits: 17
nrow(Microcebus_morphology) #number of samples: 73


## Create metadata dataset with species assignments
Microcebus_metadata <- Microcebus_multiple_data_combined %>% 
  dplyr::select(Individual_number, sex, clade, population, fine_sub_pops) %>%
  dplyr::rename(Sex = sex, Clade = clade, Population = population, Subpopulations = fine_sub_pops)
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
nrow(Microcebus_metadata) #number of samples: 73


## Import and process environmental dataset (variables extracted and processed by separate R script based on coordinates)
Microcebus_environmental <- utils::read.csv("../Empirical_examples/van_Elst_et_al_2024/Microcebus_environmental.csv", row.names = 1, stringsAsFactors = FALSE)
Microcebus_environmental_rownames <- Microcebus_environmental$Individual.ID #save IDs for later
Microcebus_environmental <- Microcebus_environmental %>% 
  dplyr::select(-Latitude, -Longitude, -Elevation, -Individual.ID)
Microcebus_environmental <- as.data.frame(lapply(Microcebus_environmental, as.numeric)) #ensure all columns are numeric
rownames(Microcebus_environmental) <- Microcebus_environmental_rownames #keep rownames
Microcebus_environmental <- (NicheDiv::transform.skewed.variables(Microcebus_environmental))$transformed #transform skewed variables
Microcebus_environmental <- remove.lowCV.multicollinearity.SOM(Microcebus_environmental, #remove highly correlated and low-variance variables
                                                               CV.threshold = 0.05,
                                                               cor.threshold = 0.9)
ncol(Microcebus_environmental) #number of variables: 56
nrow(Microcebus_environmental) #number of samples: 73


## Create spatial dataset with Latitude, Longitude and Elevation
Microcebus_spatial <- Microcebus_multiple_data_combined %>% 
  dplyr::select(latitude, longitude) %>% #add Latitude and Longitude
  dplyr::rename(Latitude = latitude, Longitude = longitude)
Microcebus_environmental_spatial <- utils::read.csv("../Empirical_examples/van_Elst_et_al_2024/Microcebus_environmental.csv", row.names = 1, stringsAsFactors = FALSE)
Microcebus_environmental_spatial <- Microcebus_environmental_spatial %>% dplyr::select(Elevation)
Microcebus_spatial$Elevation <- Microcebus_environmental_spatial[rownames(Microcebus_spatial), "Elevation"]
Microcebus_spatial$Elevation <- Microcebus_environmental_spatial[rownames(Microcebus_spatial), "Elevation"]
nrow(Microcebus_spatial) #number of samples: 73


## Match all datasets to rownames in morphology dataset
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
print(unname(round(system.time({
Microcebus_SOM_tr <- train.SOM(Microcebus_SOM_full_data, #?? samples
                               max.NA.row = 0.5,
                               max.NA.col = 0.5,
                               save.SOM.results.name = file.path(intermediate_files_folder, "Microcebus_SOM_tr.Rdata"),
                               save.SOM.results = T)
})[3] / 60, 1)))

print(unname(round(system.time({
Microcebus_SOM_kmeansBICelbow <- clustering.SOM(Microcebus_SOM_tr, max.k = 20,
                                                clustering.method = "kmeans+BICelbow",
                                                save.SOM.results.name = file.path(intermediate_files_folder, "Microcebus_SOM_kmeansBICelbow.Rdata"))
})[3] / 60, 1)))
Microcebus_SOM_kmeansBICelbow$optim_k_summary #k3 93%
print(unname(round(system.time({
Microcebus_SOM_kmeansBICthreshold <- clustering.SOM(Microcebus_SOM_tr, max.k = 20, #takes ca 2min!
                                                    clustering.method = "kmeans+BICthreshold",
                                                    save.SOM.results.name = file.path(intermediate_files_folder, "Microcebus_SOM_kmeansBICthreshold.Rdata"))
})[3] / 60, 1)))
Microcebus_SOM_kmeansBICthreshold$optim_k_summary #k3 99%
print(unname(round(system.time({
Microcebus_SOM_HDBSCAN <- clustering.SOM(Microcebus_SOM_tr, max.k = 20, #takes ca 2min!
                                         clustering.method = "HDBSCAN",
                                         save.SOM.results.name = file.path(intermediate_files_folder, "Microcebus_SOM_HDBSCAN.Rdata"))
})[3] / 60, 1)))
Microcebus_SOM_HDBSCAN$optim_k_summary #k3 31%, k2 29%, k4 22%, k5 10%
print(unname(round(system.time({
Microcebus_SOM_hierarchicalDB <- clustering.SOM(Microcebus_SOM_tr, max.k = 20, #takes ca 45min!
                                                clustering.method = "hierarchical+DB",
                                                save.SOM.results.name = file.path(intermediate_files_folder, "Microcebus_SOM_hierarchicalDB.Rdata"))
})[3] / 60, 1)))
Microcebus_SOM_hierarchicalDB$optim_k_summary #k20 81%, k19 10%
print(unname(round(system.time({
Microcebus_SOM_GMMBICthreshold <- clustering.SOM(Microcebus_SOM_tr, max.k = 20, #ca. 20min
                                                 clustering.method = "GMM+BICthreshold",
                                                 save.SOM.results.name = file.path(intermediate_files_folder, "Microcebus_SOM_GMMBICthreshold.Rdata"))
})[3] / 60, 1)))
Microcebus_SOM_GMMBICthreshold$optim_k_summary #k2 22%, k7 14%, k10 10%, k3 9%, k5 8%, k8 8%, k12 7%
print(unname(round(system.time({
Microcebus_SOM_OPTICSSilhouette <- clustering.SOM(Microcebus_SOM_tr, max.k = 20, #ca 5min
                                                  clustering.method = "OPTICS+Silhouette",
                                                  save.SOM.results.name = file.path(intermediate_files_folder, "Microcebus_SOM_OPTICSSilhouette.Rdata"))
})[3] / 60, 1)))
Microcebus_SOM_OPTICSSilhouette$optim_k_summary #k3 49%, k2 38%, k4 7%




## Evaluate and plot results for full data
Microcebus_SOM <- Microcebus_SOM_kmeansBICelbow
plot.learning.SOM(Microcebus_SOM)
plot.layer.weights.SOM(Microcebus_SOM_full)
plot.K.SOM(Microcebus_SOM)
plot.model.SOM(Microcebus_SOM, replicate.mode = "representative")
plot.model.SOM(Microcebus_SOM, replicate.mode = "first")
plot.variable.importance.SOM(Microcebus_SOM, 
                             mode = "Cluster.separation",
                             left.margin = 5.8,
                             bar.label.font.size = 0.4)
plot.variable.importance.SOM(Microcebus_SOM, 
                             mode = "Map.variance",
                             left.margin = 5.8,
                             bar.label.font.size = 0.4)
plot.structure.SOM(Microcebus_SOM, sort.by.col = 2)
plot.map.SOM(SOM.output = Microcebus_SOM, 
             Coordinates = Microcebus_spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 1, #add coordinates as buffer range around latitude coordinates
             lon.buffer.range = 1, #add coordinates as buffer range around longitude coordinates
             pie.size = 2.5, #pie chart size
             north.arrow.position = c(0.03, 0.86), #position (x, y) of north arrow relative to map
             north.arrow.length = 0.5, #length of north arrow
             north.arrow.N.position = 0.2, #position of north arrow "N"
             north.arrow.N.size = 1) #size of north arrow "N"

plot.layer.distance.scale.SOM(Microcebus_SOM)
plot.layer.importance.varimp.SOM(Microcebus_SOM, bottom.margin = 6.5)
plot.layer.importance.leaveoneout.SOM(Microcebus_SOM, 
                                      bottom.margin = 6.5,
                                      save.leave.one.layer.out.results = TRUE,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Microcebus_SOM_lolo.Rdata"))



## Compare ancestry coefficients with prior species and proposed ("revised") species labels for full data
Microcebus_ancestry <- as.data.frame(Microcebus_SOM$ancestry_matrix)
Microcebus_ancestry$Species <- Microcebus_metadata$Species[match(rownames(Microcebus_SOM$ancestry_matrix), rownames(Microcebus_metadata))]
Microcebus_ancestry$Species_revised <- Microcebus_metadata$Species[match(rownames(Microcebus_SOM$ancestry_matrix), rownames(Microcebus_metadata))]
length(unique(Microcebus_ancestry$Species)) #number of species present in data
length(unique(Microcebus_ancestry$Species_revised)) #number of proposed species present in data
table(Microcebus_ancestry$Species)
table(Microcebus_ancestry$Species_revised)


## Compare ancestry coefficients with prior species and proposed ("revised") species labels for full data
Microcebus_ancestry <- as.data.frame(Microcebus_SOM$ancestry_matrix)
Microcebus_major_cluster_index <- apply(Microcebus_SOM$ancestry_matrix, 1, which.max) #assign each sample to cluster with highest ancestry proportion
Microcebus_major_cluster_prop <- apply(Microcebus_SOM$ancestry_matrix, 1, max) #extract highest ancestry proportion per sample
Microcebus_ancestry$Major_cluster <- ifelse(Microcebus_major_cluster_prop >= 0.7,
                                            paste0("cluster", Microcebus_major_cluster_index),
                                            "mixed") #assign mixed if max ancestry proportion is below threshold
Microcebus_ancestry$Major_cluster_prop <- Microcebus_major_cluster_prop #store highest ancestry proportion
Microcebus_ancestry$Species <- Microcebus_metadata$Species[match(rownames(Microcebus_SOM$ancestry_matrix), rownames(Microcebus_metadata))]
Microcebus_ancestry$Species_revised <- Microcebus_metadata$Species_revised[match(rownames(Microcebus_SOM$ancestry_matrix), rownames(Microcebus_metadata))]
length(unique(Microcebus_ancestry$Species)) #number of species present in data
length(unique(Microcebus_ancestry$Species_revised)) #number of proposed species present in data
length(unique(Microcebus_ancestry$Major_cluster)) #number of major cluster categories present in data
table(Microcebus_ancestry$Species)
table(Microcebus_ancestry$Species_revised)
table(Microcebus_ancestry$Species, Microcebus_ancestry$Major_cluster)
table(Microcebus_ancestry$Major_cluster, Microcebus_ancestry$Species_revised)


## Hierarchical analyses based on recovered clusters
Microcebus_clusters <- apply(Microcebus_SOM$ancestry_matrix, 1, which.max) #assign each sample to cluster with highest ancestry proportion
Microcebus_clusters <- paste0("cluster", Microcebus_clusters) #rename clusters
table(Microcebus_clusters)
Microcebus_cluster_samples <- split(rownames(Microcebus_SOM$ancestry_matrix), Microcebus_clusters)
Microcebus_cluster1_data <- lapply(Microcebus_SOM$input_data, function(x) x[Microcebus_cluster_samples$cluster1, , drop = FALSE]) #cluster 1 subset
Microcebus_cluster2_data <- lapply(Microcebus_SOM$input_data, function(x) x[Microcebus_cluster_samples$cluster2, , drop = FALSE]) #cluster 2 subset
Microcebus_cluster3_data <- lapply(Microcebus_SOM$input_data, function(x) x[Microcebus_cluster_samples$cluster3, , drop = FALSE]) #cluster 3 subset

print(unname(round(system.time({
Microcebus_SOM_tr_cluster1 <- train.SOM(Microcebus_cluster1_data, #?? samples
                                        grid.multiplier = 3,
                                        max.NA.row = 0.5,
                                        max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Microcebus_SOM_cluster1 <- clustering.SOM(Microcebus_SOM_tr_cluster1
                                          ,max.k = 5,
                                          clustering.method = "kmeans+BICelbow")
})[3] / 60, 1)))
Microcebus_SOM_cluster1$optim_k_summary #k1 99%

print(unname(round(system.time({
Microcebus_SOM_tr_cluster2 <- train.SOM(Microcebus_cluster2_data, #? samples
                                        grid.multiplier = 5,
                                        max.NA.row = 0.5,
                                        max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Microcebus_SOM_cluster2 <- clustering.SOM(Microcebus_SOM_tr_cluster2,
                                          clustering.method = "kmeans+BICelbow",
                                          max.k = 10)
})[3] / 60, 1)))
Microcebus_SOM_cluster2$optim_k_summary #k1 98%

print(unname(round(system.time({
Microcebus_SOM_tr_cluster3 <- train.SOM(Microcebus_cluster3_data, #?? samples
                                        grid.multiplier = 3,
                                        max.NA.row = 0.5,
                                        max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Microcebus_SOM_cluster3 <- clustering.SOM(Microcebus_SOM_tr_cluster3,
                                          clustering.method = "kmeans+BICelbow",
                                          max.k = 10)
})[3] / 60, 1)))
Microcebus_SOM_cluster3$optim_k_summary #k2 89%, k10 11%


plot.model.SOM(Microcebus_SOM_cluster3, replicate.mode = "representative")
plot.model.SOM(Microcebus_SOM_cluster3, replicate.mode = "first")
plot.structure.SOM(Microcebus_SOM_cluster3)
plot.K.SOM(Microcebus_SOM_cluster3)
plot.map.SOM(SOM.output = Microcebus_SOM_cluster3,
             Coordinates = Microcebus_spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 0.1, #add coordinates as buffer range around latitude coordinates
             lon.buffer.range = 0.3, #add coordinates as buffer range around longitude coordinates
             pie.size = 2.5, #pie chart size
             north.arrow.position = c(0.04, 0.85), #position (x, y) of north arrow relative to map
             north.arrow.length = 0.1, #length of north arrow
             north.arrow.N.position = 0.05, #position of north arrow "N"
             north.arrow.N.size = 1) #size of north arrow "N"
plot.variable.importance.SOM(Microcebus_SOM_cluster3,
                             mode = "Cluster.separation", 
                             left.margin = 6.5,
                             bar.label.font.size = 0.4)
plot.variable.importance.SOM(Microcebus_SOM_cluster3, 
                             mode = "Map.variance", 
                             left.margin = 6.5,
                             bar.label.font.size = 0.4)
plot.layer.importance.varimp.SOM(Microcebus_SOM_cluster3, bottom.margin = 6)
plot.layer.importance.leaveoneout.SOM(Microcebus_SOM_cluster3, 
                                      save.leave.one.layer.out.results = TRUE,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Microcebus_SOM_cluster3_lolo.Rdata")) #this will take 10-20min (running 2 x N replicates for train and clustering SOM)

Microcebus_ancestry_SOM_cluster3 <- as.data.frame(Microcebus_SOM_cluster3$ancestry_matrix)
Microcebus_ancestry_SOM_cluster3$Species <- Microcebus_metadata$Species[match(rownames(Microcebus_SOM_cluster3$ancestry_matrix), rownames(Microcebus_metadata))]
Microcebus_ancestry_SOM_cluster3$Species_revised <- Microcebus_metadata$Species[match(rownames(Microcebus_SOM_cluster3$ancestry_matrix), rownames(Microcebus_metadata))]
length(unique(Microcebus_ancestry_SOM_cluster3$Species)) #number of species present in data
length(unique(Microcebus_ancestry_SOM_cluster3$Species_revised)) #number of proposed species present in data
table(Microcebus_ancestry_SOM_cluster3$Species)
table(Microcebus_ancestry_SOM_cluster3$Species_revised)




#### Elysia sea slugs from the Western Atlantic (Krug et al. 2026) #############

## Import and filter mitochondrial DNA data
Elysia_COI <- process.SNP.data.SOM(nexus.path = "../Empirical_examples/Krug_et_al_2026/Elysia_mtDNA_expanded.nex",
                                   missing.loci.cutoff.lenient = 0.7,
                                   missing.loci.cutoff.final = 0.5,
                                   missing.individuals.cutoff = 0.5)
ncol(Elysia_COI) #number of loci: 160
nrow(Elysia_COI) #number of samples: 282


## Import meta data
Elysia_metadata <- read.csv("../Empirical_examples/Krug_et_al_2026/Elysia_metadata_updated.csv",
                            header = TRUE,
                            stringsAsFactors = FALSE,
                            check.names = FALSE)
rownames(Elysia_metadata) <- Elysia_metadata$New_sample_name
nrow(Elysia_metadata) #number of samples: 282


## Extract spatial data (Coordinates and seadepth; latter obtained via marmap::getNOAA.bathy at 0.25 arcmin resolution)
Elysia_spatial <- dplyr::select(Elysia_metadata, Latitude, Longitude, Seadepth) #create dataframe with Lat and Long
rownames(Elysia_spatial) <- rownames(Elysia_metadata)
nrow(Elysia_spatial) #number of samples: 282


## Import environmental data (obtained via geodata::bio_oracle)
Elysia_environmental <- read.csv("../Empirical_examples/Krug_et_al_2026/Elysia_environmental.csv",
                                 header = TRUE,
                                 stringsAsFactors = FALSE,
                                 row.names = 1)
Elysia_environmental <- (NicheDiv::transform.skewed.variables(Elysia_environmental))$transformed #transform skewed variables
Elysia_environmental <- remove.lowCV.multicollinearity.SOM(Elysia_environmental, #remove highly correlated and low-variance variables
                                                           CV.threshold = 0.05,
                                                           cor.threshold = 0.9)
nrow(Elysia_environmental) #number of samples: 282
ncol(Elysia_environmental) #number of variables: 5


## Extract host data
Elysia_host <- dplyr::select(Elysia_metadata, Host)
rownames(Elysia_host) <- rownames(Elysia_metadata)
Elysia_host <- make.cols.binary.SOM(Elysia_host, #convert Host to binary columns and keep original
                                    make.binary.cols = "Host",
                                    append.to.original = TRUE)


## Extract developmental type data
Elysia_development <- dplyr::select(Elysia_metadata, Developmental_mode)
rownames(Elysia_development) <- rownames(Elysia_metadata)
Elysia_development <- make.cols.binary.SOM(Elysia_development, #convert Developmental_mode to binary columns and keep original
                                           make.binary.cols = "Developmental_mode",
                                           append.to.original = TRUE)


## Combine host and developmental data
Elysia_host_development <- cbind(Elysia_host, Elysia_development)
Elysia_host_development <- Elysia_host_development[, !duplicated(colnames(Elysia_host_development)), drop = FALSE]
Elysia_host_development <- Elysia_host_development[rownames(Elysia_metadata), , drop = FALSE]
nrow(Elysia_host_development) #number of samples: 282
ncol(Elysia_host_development) #number of variables: 7


## Match datasets and remove rows that are entirely NA within each dataset
Elysia_common_IDs <- Reduce(intersect, list(
  rownames(Elysia_COI),
  rownames(Elysia_spatial),
  rownames(Elysia_environmental),
  rownames(Elysia_metadata),
  rownames(Elysia_host_development)
))

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
Elysia_common_IDs <- Reduce(intersect, list( #re-match after NA filtering
  rownames(Elysia_COI),
  rownames(Elysia_spatial),
  rownames(Elysia_environmental),
  rownames(Elysia_metadata),
  rownames(Elysia_host_development)
))
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


## Train and cluster SOM
Elysia_all_data <- list(mtDNA = Elysia_COI,
                        Host_development = Elysia_host_development,
                        Environmental = Elysia_environmental,
                        Spatial = Elysia_spatial)
print(unname(round(system.time({
Elysia_SOM_tr <- train.SOM(input_data = Elysia_all_data, #?? samples
                           max.NA.row = 0.5,
                           max.NA.col = 0.5,
                           save.SOM.results = TRUE,
                           save.SOM.results.name = file.path(intermediate_files_folder, "Elysia_SOM_tr.Rdata"))
})[3] / 60, 1)))

print(unname(round(system.time({
Elysia_SOM_kmeansBICelbow <- clustering.SOM(Elysia_SOM_tr, max.k = 10,
                                            clustering.method = "kmeans+BICelbow",
                                            save.SOM.results.name = file.path(intermediate_files_folder, "Elysia_SOM_kmeansBICelbow.Rdata"))
})[3] / 60, 1)))
Elysia_SOM_kmeansBICelbow$optim_k_summary #k4 96%
print(unname(round(system.time({
Elysia_SOM_kmeansBICthreshold <- clustering.SOM(Elysia_SOM_tr, max.k = 10, #takes ca 2min!
                                                clustering.method = "kmeans+BICthreshold",
                                                save.SOM.results.name = file.path(intermediate_files_folder, "Elysia_SOM_kmeansBICthreshold.Rdata"))
})[3] / 60, 1)))
Elysia_SOM_kmeansBICthreshold$optim_k_summary #k8 43%, k7 34%, k9 10%, k6 8%
print(unname(round(system.time({
Elysia_SOM_HDBSCAN <- clustering.SOM(Elysia_SOM_tr, max.k = 10, #takes ca 2min!
                                     clustering.method = "HDBSCAN",
                                     save.SOM.results.name = file.path(intermediate_files_folder, "Elysia_SOM_HDBSCAN.Rdata"))
})[3] / 60, 1)))
Elysia_SOM_HDBSCAN$optim_k_summary #k3 39%, k4 39%, k5 12%, k6 7%
print(unname(round(system.time({
Elysia_SOM_hierarchicalDB <- clustering.SOM(Elysia_SOM_tr, max.k = 10, #takes ca 45min!
                                            clustering.method = "hierarchical+DB",
                                            save.SOM.results.name = file.path(intermediate_files_folder, "Elysia_SOM_hierarchicalDB.Rdata"))
})[3] / 60, 1)))
Elysia_SOM_hierarchicalDB$optim_k_summary #k5 39%, k4 38%, k6 8%, k15 8%
print(unname(round(system.time({
Elysia_SOM_GMMBICthreshold <- clustering.SOM(Elysia_SOM_tr, max.k = 10, #ca. 15min
                                             clustering.method = "GMM+BICthreshold",
                                             save.SOM.results.name = file.path(intermediate_files_folder, "Elysia_SOM_GMMBICthreshold.Rdata"))
})[3] / 60, 1)))
Elysia_SOM_GMMBICthreshold$optim_k_summary #k7 23%, k9 21%, k8 16%, k6 7%
print(unname(round(system.time({
Elysia_SOM_OPTICSSilhouette <- clustering.SOM(Elysia_SOM_tr, max.k = 10, #ca 5min
                                              clustering.method = "OPTICS+Silhouette",
                                              save.SOM.results.name = file.path(intermediate_files_folder, "Elysia_SOM_OPTICSSilhouette.Rdata"))
})[3] / 60, 1)))
Elysia_SOM_OPTICSSilhouette$optim_k_summary #k4 84%, k3 10%

Elysia_SOM <- Elysia_SOM_kmeansBICelbow
plot.structure.SOM(Elysia_SOM, bottom.margin = 7, Individual.labels.font.size = 0.2)
plot.learning.SOM(Elysia_SOM)
plot.layer.distance.scale.SOM(Elysia_SOM)
plot.K.SOM(Elysia_SOM)
plot.model.SOM(Elysia_SOM, replicate.mode = "representative")
plot.model.SOM(Elysia_SOM, replicate.mode = "first")
plot.variable.importance.SOM(Elysia_SOM, left.margin = 12)
plot.variable.importance.SOM(Elysia_SOM, left.margin = 12, mode = "Map.variance")
plot.map.SOM(SOM.output = Elysia_SOM,
             Coordinates = Elysia_spatial[, c("Latitude", "Longitude")],
             USA.add.counties = F,
             north.arrow.length = 1.5,
             north.arrow.position = c(0.0001, 0.95),
             north.arrow.N.position = 0.6,
             lon.buffer.range = 3,
             scale.position = c(0.6, 0.001))
Elysia_SOM$clustering.SOM.args
plot.layer.importance.leaveoneout.SOM(Elysia_SOM, 
                                      bottom.margin = 6.5,
                                      save.leave.one.layer.out.results = TRUE,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Elysia_SOM_SOM_lolo.Rdata"))
plot.layer.importance.varimp.SOM(Elysia_SOM, bottom.margin = 7.5)


Elysia_species_names <- dplyr::select(Elysia_metadata, Species = Species_name)
rownames(Elysia_species_names) <- rownames(Elysia_metadata)
Elysia_SOM_ancestry_matrix <- as.data.frame(Elysia_SOM$ancestry_matrix)
Elysia_SOM_common_samples <- intersect(rownames(Elysia_SOM_ancestry_matrix), rownames(Elysia_species_names))
Elysia_SOM_ancestry_mat_sub <- Elysia_SOM_ancestry_matrix[Elysia_SOM_common_samples, , drop = FALSE]
Elysia_SOM_species_names_sub <- Elysia_species_names[Elysia_SOM_common_samples, , drop = FALSE]
Elysia_SOM_ancestry_matrix <- cbind(Elysia_SOM_ancestry_mat_sub, Elysia_SOM_species_names_sub) #compare ancestry matrix with species hypotheses
head(Elysia_SOM_ancestry_matrix)


## Assign majority cluster only if ancestry ≥ ancestry_limit
ancestry_limit <- 0.8
cluster.cols <- grep("^Cluster_", colnames(Elysia_SOM_ancestry_matrix), value = TRUE)
cluster.matrix <- Elysia_SOM_ancestry_matrix[, cluster.cols, drop = FALSE]
max.cluster.idx <- max.col(cluster.matrix, ties.method = "first")
max.cluster.val <- apply(cluster.matrix, 1, max)
Elysia_SOM_ancestry_matrix$Majority_cluster <- "Mix"
Elysia_SOM_ancestry_matrix$Majority_cluster[max.cluster.val >= ancestry_limit] <- cluster.cols[max.cluster.idx[max.cluster.val >= ancestry_limit]]
table(Elysia_SOM_ancestry_matrix$Majority_cluster, useNA = "ifany")


## Species × majority cluster table
table(Elysia_SOM_ancestry_matrix$Species, Elysia_SOM_ancestry_matrix$Majority_cluster, useNA = "ifany")


## Hierarchical analyses based on recovered clusters
Elysia_clusters <- apply(Elysia_SOM$ancestry_matrix, 1, which.max) #assign each sample to cluster with highest ancestry proportion
Elysia_clusters <- paste0("cluster", Elysia_clusters) #rename clusters
table(Elysia_clusters)
Elysia_cluster_samples <- split(rownames(Elysia_SOM$ancestry_matrix), Elysia_clusters)
Elysia_cluster1_data <- lapply(Elysia_SOM$input_data, function(x) x[Elysia_cluster_samples$cluster1, , drop = FALSE]) #cluster 1 subset
Elysia_cluster2_data <- lapply(Elysia_SOM$input_data, function(x) x[Elysia_cluster_samples$cluster2, , drop = FALSE]) #cluster 2 subset
Elysia_cluster3_data <- lapply(Elysia_SOM$input_data, function(x) x[Elysia_cluster_samples$cluster3, , drop = FALSE]) #cluster 3 subset
Elysia_cluster4_data <- lapply(Elysia_SOM$input_data, function(x) x[Elysia_cluster_samples$cluster4, , drop = FALSE]) #cluster 4 subset

print(unname(round(system.time({
Elysia_SOM_tr_cluster1 <- train.SOM(Elysia_cluster1_data, #? samples
                                    grid.multiplier = 4,
                                    max.NA.row = 0.5,
                                    max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Elysia_SOM_cluster1 <- clustering.SOM(Elysia_SOM_tr_cluster1,
                                      clustering.method = "kmeans+BICelbow",
                                      max.k = 5)
})[3] / 60, 1)))
Elysia_SOM_cluster1$optim_k_summary #k3 84%, k4 11%

print(unname(round(system.time({
Elysia_SOM_tr_cluster2 <- train.SOM(Elysia_cluster2_data, #63 samples
                                    grid.multiplier = 5,
                                    max.NA.row = 0.5,
                                    max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Elysia_SOM_cluster2 <- clustering.SOM(Elysia_SOM_tr_cluster2,
                                      clustering.method = "kmeans+BICelbow",
                                      max.k = 10)
})[3] / 60, 1)))
Elysia_SOM_cluster2$optim_k_summary #k4 68%, k5 27%

print(unname(round(system.time({
Elysia_SOM_tr_cluster3 <- train.SOM(Elysia_cluster3_data[names(Elysia_cluster3_data) != "Host_development"],
                                    grid.multiplier = 5, #50 samples
                                    max.NA.row = 0.5,
                                    max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Elysia_SOM_cluster3 <- clustering.SOM(Elysia_SOM_tr_cluster3,
                                      clustering.method = "kmeans+BICelbow",
                                      max.k = 10)
})[3] / 60, 1)))
Elysia_SOM_cluster3$optim_k_summary #k1 85%, k3 6%

print(unname(round(system.time({
Elysia_SOM_tr_cluster4 <- train.SOM(Elysia_cluster4_data[names(Elysia_cluster4_data) != "Host_development"],
                                    grid.multiplier = 5, #131 samples
                                    max.NA.row = 0.5,
                                    max.NA.col = 0.5)
})[3] / 60, 1)))
print(unname(round(system.time({
Elysia_SOM_cluster4 <- clustering.SOM(Elysia_SOM_tr_cluster4,
                                      clustering.method = "kmeans+BICelbow",
                                      max.k = 25)
})[3] / 60, 1)))
Elysia_SOM_cluster4$optim_k_summary #k2 65%, k4 13%


## Cluster 1
plot.model.SOM(Elysia_SOM_cluster1, replicate.mode = "representative")
plot.model.SOM(Elysia_SOM_cluster1, replicate.mode = "first")
plot.structure.SOM(Elysia_SOM_cluster1)
plot.K.SOM(Elysia_SOM_cluster1)
plot.map.SOM(SOM.output = Elysia_SOM_cluster1,
             Coordinates = Elysia_spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 0.1, #add coordinates as buffer range around latitude coordinates
             lon.buffer.range = 0.3, #add coordinates as buffer range around longitude coordinates
             pie.size = 2.5, #pie chart size
             north.arrow.position = c(0.04, 0.85), #position (x, y) of north arrow relative to map
             north.arrow.length = 0.1, #length of north arrow
             north.arrow.N.position = 0.05, #position of north arrow "N"
             north.arrow.N.size = 1) #size of north arrow "N"
plot.variable.importance.SOM(Elysia_SOM_cluster1,
                             mode = "Cluster.separation", 
                             left.margin = 10,
                             bar.label.font.size = 0.6)
plot.variable.importance.SOM(Elysia_SOM_cluster1, 
                             mode = "Map.variance", 
                             left.margin = 10,
                             bar.label.font.size = 0.6)
plot.layer.importance.varimp.SOM(Elysia_SOM_cluster1, bottom.margin = 8)
plot.layer.importance.leaveoneout.SOM(Elysia_SOM_cluster1, 
                                      save.leave.one.layer.out.results = TRUE,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Elysia_SOM_cluster1_lolo.Rdata")) #this will take 10-20min (running 2 x N replicates for train and clustering SOM)

Elysia_ancestry_SOM_cluster1 <- as.data.frame(Elysia_SOM_cluster1$ancestry_matrix)
Elysia_cluster1_major_cluster_index <- apply(Elysia_SOM_cluster1$ancestry_matrix, 1, which.max) #assign each sample to cluster with highest ancestry proportion
Elysia_cluster1_major_cluster_prop <- apply(Elysia_SOM_cluster1$ancestry_matrix, 1, max) #extract highest ancestry proportion per sample
Elysia_ancestry_SOM_cluster1$Major_cluster <- ifelse(Elysia_cluster1_major_cluster_prop >= 0.8,
                                                     paste0("cluster", Elysia_cluster1_major_cluster_index),
                                                     "mixed") #assign mixed if max ancestry proportion is below threshold
Elysia_ancestry_SOM_cluster1$Major_cluster_prop <- Elysia_cluster1_major_cluster_prop #store highest ancestry proportion
Elysia_ancestry_SOM_cluster1$Species <- Elysia_metadata$Species_name[match(rownames(Elysia_SOM_cluster1$ancestry_matrix), rownames(Elysia_metadata))]
length(unique(Elysia_ancestry_SOM_cluster1$Species)) #number of species present in data
length(unique(Elysia_ancestry_SOM_cluster1$Major_cluster)) #number of major cluster categories present in data
table(Elysia_ancestry_SOM_cluster1$Species)
table(Elysia_ancestry_SOM_cluster1$Major_cluster, Elysia_ancestry_SOM_cluster1$Species)


## Cluster 2
plot.model.SOM(Elysia_SOM_cluster2, replicate.mode = "representative")
plot.model.SOM(Elysia_SOM_cluster2, replicate.mode = "first")
plot.structure.SOM(Elysia_SOM_cluster2, bottom.margin = 9)
plot.K.SOM(Elysia_SOM_cluster2)
plot.map.SOM(SOM.output = Elysia_SOM_cluster2,
             Coordinates = Elysia_spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 0.1, #add coordinates as buffer range around latitude coordinates
             lon.buffer.range = 0.3, #add coordinates as buffer range around longitude coordinates
             pie.size = 2.5, #pie chart size
             north.arrow.position = c(0.04, 0.85), #position (x, y) of north arrow relative to map
             north.arrow.length = 0.1, #length of north arrow
             north.arrow.N.position = 0.05, #position of north arrow "N"
             north.arrow.N.size = 1) #size of north arrow "N"
plot.variable.importance.SOM(Elysia_SOM_cluster2,
                             mode = "Cluster.separation", 
                             left.margin = 10,
                             bar.label.font.size = 0.6)
plot.variable.importance.SOM(Elysia_SOM_cluster2, 
                             mode = "Map.variance", 
                             left.margin = 10,
                             bar.label.font.size = 0.6)
plot.layer.importance.varimp.SOM(Elysia_SOM_cluster2, bottom.margin = 8)
plot.layer.importance.leaveoneout.SOM(Elysia_SOM_cluster2, 
                                      save.leave.one.layer.out.results = TRUE,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Elysia_SOM_cluster2_lolo.Rdata")) #this will take 10-20min (running 2 x N replicates for train and clustering SOM)

Elysia_ancestry_SOM_cluster2 <- as.data.frame(Elysia_SOM_cluster2$ancestry_matrix)
Elysia_cluster2_major_cluster_index <- apply(Elysia_SOM_cluster2$ancestry_matrix, 1, which.max) #assign each sample to cluster with highest ancestry proportion
Elysia_cluster2_major_cluster_prop <- apply(Elysia_SOM_cluster2$ancestry_matrix, 1, max) #extract highest ancestry proportion per sample
Elysia_ancestry_SOM_cluster2$Major_cluster <- ifelse(Elysia_cluster2_major_cluster_prop >= 0.6,
                                                     paste0("cluster", Elysia_cluster2_major_cluster_index),
                                                     "mixed") #assign mixed if max ancestry proportion is below threshold
Elysia_ancestry_SOM_cluster2$Major_cluster_prop <- Elysia_cluster2_major_cluster_prop #store highest ancestry proportion
Elysia_ancestry_SOM_cluster2$Species <- Elysia_metadata$Species_name[match(rownames(Elysia_SOM_cluster2$ancestry_matrix), rownames(Elysia_metadata))]
length(unique(Elysia_ancestry_SOM_cluster2$Species)) #number of species present in data
length(unique(Elysia_ancestry_SOM_cluster2$Major_cluster)) #number of major cluster categories present in data
table(Elysia_ancestry_SOM_cluster2$Species)
table(Elysia_ancestry_SOM_cluster2$Major_cluster, Elysia_ancestry_SOM_cluster2$Species)


## Cluster 4
print(unname(round(system.time({
Elysia_SOM_cluster4_k2_k2 <- clustering.SOM(Elysia_SOM_tr_cluster4,
                                            set.k = 2,
                                            clustering.method = "kmeans+BICelbow")
})[3] / 60, 1)))
plot.model.SOM(Elysia_SOM_cluster4_k2, replicate.mode = "representative")
plot.model.SOM(Elysia_SOM_cluster4_k2, replicate.mode = "first")
plot.structure.SOM(Elysia_SOM_cluster4_k2, bottom.margin = 9)
plot.K.SOM(Elysia_SOM_cluster4_k2)
plot.map.SOM(SOM.output = Elysia_SOM_cluster4_k2,
             Coordinates = Elysia_spatial[, c("Latitude", "Longitude")],
             lat.buffer.range = 0.1, #add coordinates as buffer range around latitude coordinates
             lon.buffer.range = 0.3, #add coordinates as buffer range around longitude coordinates
             pie.size = 2.5, #pie chart size
             north.arrow.position = c(0.04, 0.85), #position (x, y) of north arrow relative to map
             north.arrow.length = 0.1, #length of north arrow
             north.arrow.N.position = 0.05, #position of north arrow "N"
             north.arrow.N.size = 1) #size of north arrow "N"
plot.variable.importance.SOM(Elysia_SOM_cluster4_k2,
                             mode = "Cluster.separation", 
                             left.margin = 10,
                             bar.label.font.size = 0.6)
plot.variable.importance.SOM(Elysia_SOM_cluster4_k2, 
                             mode = "Map.variance", 
                             left.margin = 10,
                             bar.label.font.size = 0.6)
plot.layer.importance.varimp.SOM(Elysia_SOM_cluster4_k2, bottom.margin = 8)
plot.layer.importance.leaveoneout.SOM(Elysia_SOM_cluster4_k2, 
                                      save.leave.one.layer.out.results = TRUE,
                                      save.leave.one.layer.out.results.name = file.path(intermediate_files_folder, "Elysia_SOM_cluster4_k2_lolo.Rdata")) #this will take 10-20min (running 2 x N replicates for train and clustering SOM)

Elysia_ancestry_SOM_cluster4 <- as.data.frame(Elysia_SOM_cluster4_k2$ancestry_matrix)
Elysia_cluster4_major_cluster_index <- apply(Elysia_SOM_cluster4_k2$ancestry_matrix, 1, which.max) #assign each sample to cluster with highest ancestry proportion
Elysia_cluster4_major_cluster_prop <- apply(Elysia_SOM_cluster4_k2$ancestry_matrix, 1, max) #extract highest ancestry proportion per sample
Elysia_ancestry_SOM_cluster4$Major_cluster <- ifelse(Elysia_cluster4_major_cluster_prop >= 0.6,
                                                     paste0("cluster", Elysia_cluster4_major_cluster_index),
                                                     "mixed") #assign mixed if max ancestry proportion is below threshold
Elysia_ancestry_SOM_cluster4$Major_cluster_prop <- Elysia_cluster4_major_cluster_prop #store highest ancestry proportion
Elysia_ancestry_SOM_cluster4$Species <- Elysia_metadata$Species_name[match(rownames(Elysia_SOM_cluster4_k2$ancestry_matrix), rownames(Elysia_metadata))]
length(unique(Elysia_ancestry_SOM_cluster4$Species)) #number of species present in data
length(unique(Elysia_ancestry_SOM_cluster4$Major_cluster)) #number of major cluster categories present in data
table(Elysia_ancestry_SOM_cluster4$Species)
table(Elysia_ancestry_SOM_cluster4$Major_cluster, Elysia_ancestry_SOM_cluster4$Species)
