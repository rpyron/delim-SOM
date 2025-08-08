#####################################################################################################################
#### Empirical dataset: Monticola - original GBS data from Pyron et al. 2023 with updated environmental variables ###
#####################################################################################################################
rm(list = ls()) #clear environment
setwd("./")
source("./R/delim-SOM_2.0_functions.R")

## Read in sample data
Monticola_data <- read.csv(file = "Data/monticola71/monticola71.csv",
                            row.names = 1,
                            header = T, colClasses = c(huc2 = "character",
                                                       huc4 = "character",
                                                       huc6 = "character",
                                                       huc8 = "character",
                                                       huc10 = "character",
                                                       huc12 = "character"))


## Import and process genetic SNP data
Monticola_SNP <- process.SNP.data.SOM(vcf.path = "Data/monticola71/monticola71.vcf.gz",
                                   missing.loci.cutoff.lenient = 0.2,
                                   missing.loci.cutoff.final = 0.2,
                                   missing.individuals.cutoff = 0.6)
rownames(Monticola_SNP) <- Monticola_data$Sample[match(rownames(Monticola_SNP), rownames(Monticola_data))] #rename alleles
head(Monticola_SNP[,1:10])#check the SNPs
dim(Monticola_SNP) #number of loci: 2669, number of samples: 71


## Create spatial dataset with coordinates, elevation, and HUC
Monticola_spatial <- Monticola_data[, c("Latitude",
                                          "Longitude",
                                          "Elevation")] #extract coordinates and elevation
Monticola_spatial$Elevation <- log(Monticola_spatial$Elevation) #log Elevation
rownames(Monticola_spatial) <- Monticola_data$Sample #assign rownames
head(Monticola_spatial)


## Create environmental dataset, one-hot encoding lvl4
env_vars <- grep("^(wc|current|huc|lvl)", names(Monticola_data), value = TRUE)
Monticola_environmental <- Monticola_data[env_vars]
rownames(Monticola_environmental) <- Monticola_data$Sample
Monticola_environmental <- make.cols.binary.SOM(dataframe = Monticola_environmental, #make watershed columns binary
                                             make.binary.cols = c("lvl4", "huc2", "huc4", "huc6", "huc8", "huc10", "huc12"),
                                             remove.original.cols = TRUE, #remove.original.cols - if TRUE, remove original categorical columns after processing
                                             append.to.original = TRUE #append.to.original - if TRUE, append to input; if FALSE, return only binary indicators
                                                ) #select watershed columns
Monticola_environmental <- remove.lowCV.multicollinearity.SOM(Monticola_environmental, #remove highly correlated and low-variance variables
                                                           CV.threshold = 0.05,
                                                           cor.threshold = 0.7)
head(Monticola_environmental)


## Create morphological trait dataset
trait_names <- c("SVL", "TL", "AG", "CW", "FL", "HL", "SG", "TW", "TO", "FI", "HW", "ED", "IN", "ES", "ON", "IO", "IC")
log_traits <- log(Monticola_data[, trait_names])
SVL <- log_traits[, "SVL"]
residuals_mat <- sapply(log_traits[, trait_names[-1]], function(y) resid(lm(y ~ SVL)))
Monticola_morphology <- cbind(SVL = SVL, residuals_mat)
rownames(Monticola_morphology) <- Monticola_data$Sample
head(Monticola_morphology)


## Train and cluster SOM
SOM_Monticola_SOM_data <- list(Alleles = Monticola_SNP,
                                Spatial = Monticola_spatial,
                                Environmental = Monticola_environmental,
                                Morphology = Monticola_morphology)
SOM_Monticola <- train.SOM(input_data = SOM_Monticola_SOM_data,
                            N.steps = 100, #number of training iterations S for SOM
                            N.replicates = 30, #number of SOM runs R
                            parallel = TRUE, #run SOM training in parallel 
                            N_cores = 4, #number of cores for training SOM in parallel
                            grid.size = NULL, #grid size - specify as c(x, y) if desired
                            grid.multiplier = 5, #tuning parameter defining how “fine” or “coarse” SOM grid is relative to number of samples (recommended: 5)
                            learning.rate.initial = 0.5, #initial learning rate for SOM training
                            learning.rate.final = 0.1, #final learning rate for SOM training
                            learning.rate.tuning = FALSE, #set learning.rate tuning to choose best initial and final learning rate values
                            max.NA.row = 0.4, #maximum fraction of missing values allowed per row (sample) in input data to prevent row to be removed
                            max.NA.col = 0.7, #maximum fraction of missing values allowed per column (variable) in input data to prevent row to be removed
                            training.neighborhoods = "gaussian", #neighborhood function used for SOM training (options: "gaussian" or "bubble")
                            save.SOM.results = FALSE, #whether to save SOM results to file
                            save.SOM.results.name = NULL, #file name for saving SOM results (if NULL, default name based on input_data is generated)
                            overwrite.SOM.results = TRUE, #if FALSE, existing results are loaded instead of re-running SOM
                            verbose = TRUE, #show messages
                            message.N.replicates = 1 #frequency of progress messages during training (message is printed every message.N.replicates iterations)
                            )
SOM_Monticola_clusters <- clustering.SOM(SOM.output = SOM_Monticola, 
                                          max.k = 7, #maximum of considered clusters K + 1
                                          set.k = NULL, #set to test single value of K
                                          clustering.method = "kmeans+BICjump+hierarch", 
                                          pca.codebooks = FALSE, #run PCA on codebook vectors
                                          pca.codebooks.var.threshold = 0.99, #variance threshold for codebook PCA to choose how many PCs to retain
                                          BIC.thresh = 6 #BIC threshold for selecting K > 1 - we suggest using Raftery (1995) ranges: 2, 6, or 10 for weak, medium or strong support
                                          )


## Evaluate and plot results
plot.Learning.SOM(SOM_Monticola_clusters)
plot.Layers.SOM(SOM_Monticola_clusters)
plot.K.SOM(SOM_Monticola_clusters)
plot.Model.SOM(SOM_Monticola_clusters)
plot.Structure.SOM(SOM_Monticola_clusters)
plot.Map.SOM(SOM.output = SOM_Monticola_clusters, 
             Coordinates = Monticola_spatial[, c("Latitude", "Longitude")], 
             USA.add.counties = T,
             north.arrow.position = c(0.05, 0.9),
             north.arrow.length = 0.4,
             north.arrow.N.position = 0.15,
             scale.position = c(0.79, 0.05))
plot.Variable.Importance.SOM(SOM_Monticola_clusters,
                             left.margin = 5.8,
                             bar.label.font.size = 0.4, 
                             importance.threshold = 0.1)

#Variable importance across layers
head(round(sort(SOM_Monticola_clusters$median_variable_importance[[1]], decreasing = T), 2), 400)
round(sort(SOM_Monticola_clusters$median_variable_importance[[2]], decreasing = T), 2)
round(sort(SOM_Monticola_clusters$median_variable_importance[[3]], decreasing = T), 2)
round(sort(SOM_Monticola_clusters$median_variable_importance[[4]], decreasing = T), 2)

