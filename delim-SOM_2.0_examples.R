################################################################################
#### SOM-based species delimitation (Pyron et al. 2022, 2023)
################################################################################

## Use single-layer (single data matrix) or multi-layer (multiple data matrices) Kohonen Self-Organizing Maps (SOMs or SuperSOMs) to infer clusters and delimit species using unsupervised machine learning
## SOMs are used to reduce high-dimensional data into two-dimensional grid to which individuals are assigned based on their clustering of similar values across input features
## Number of clusters are chosen by k-means clustering of grid cells into proximate units based on BIC from the weighted sum of squares (WSS) of neighbor distances between adjacent cells
## Uses kohonen package (Wehrens & Buydens 2007)
## Contribution of each layer to final model output is recorded along with clustering assignment of each individual over multiple learning replicates
## Results are similar to STRUCTURE analysis including admixture estimates but represent unified delimitation model that incorporates various dimensions of ecological and evolutionary divergence for integrative taxonomy
## If only allelic data are used, then assignment probabilities approximate individual ancestry coefficients
## If multiple layers are used, assignment probabilities represent "species coefficients"
## Method is flexible and can take almost any data type as matrix or dataframe
## Each matrix needs same number of rows (representing individuals or populations) with similar rownames but each can have any number of columns (representing variables such as allele 1 or trait 1)
## Examples of input matrices: allele frequencies (biallelic SNP matrix), morphological trait data, behavioral data, spatial data (latitude, longitude, elevation), climatic data, host or habitat data, physiological data
## Data are normalized to range between 0 and 1 (if not done so already)

rm(list = ls()) #clear environment
setwd("./");source("./R/kohonen_2.0_code.R")

#### Simulate three datasets with N individuals: allele data, environmental data and morphological data

## Specify parameters for all datasets

set.seed(1) #set seed for reproducibility
n_individuals <- 80 #simulate dataset with N individuals (needs to be consistent across all datasets)
rownames_datasets <- "Individual" #name rows (needs to be consistent across all datasets)


## Allele data (Alleles)
n_Alleles <- 120 #simulate dataset with 120 alleles
allele_frequencies <- c('1' = 0.4, '0.5' = 0.25, '0' = 0.35) #set allele frequencies
Alleles <- data.frame(lapply(1:n_Alleles, function(x) sample(names(allele_frequencies), n_individuals, replace = T, prob = allele_frequencies)),
                      row.names = paste0(rownames_datasets, 1:n_individuals)) #generate Alleles based on these frequencies
Alleles[] <- lapply(Alleles, function(x) as.numeric(as.character(x))) #convert character values to numeric (0, 0.5, 1)
colnames(Alleles) <- paste("Allele", 1:n_Alleles) #rename columns as "Allele 1", "Allele 2", ..., "Allele n"
head(Alleles);hist(unlist(Alleles))


## Environmental data (ENV)
n_env <- 38 #simulate dataset with 38 environmental variables
ENV <- matrix(runif(n_individuals * n_env, min = 0, max = 100), nrow = n_individuals, ncol = n_env)
rownames(ENV) <- paste0(rownames_datasets, 1:n_individuals)
na_indices_env <- sample(1:(n_individuals * n_env), size = round(n_individuals * n_env * 0.1), replace = F) #introduce some rare NAs for realism
colnames(ENV) <- paste("BIO", 1:n_env) #rename columns as "BIO 1", "BIO 2", ..., "BIO n"
ENV[na_indices_env] <- NA;ENV <- as.data.frame(ENV)
head(ENV);pairs(ENV[,1:4])


## Morphological data (MORPH)
n_morph <- 21 #simulate morphological data for 21 traits
MORPH <- matrix(rnorm(n_individuals * n_morph, mean = 5, sd = 2), nrow = n_individuals, ncol = n_morph) #simulate morphological data
rownames(MORPH) <- paste0(rownames_datasets, 1:n_individuals) #name rows
na_indices_morph <- sample(1:(n_individuals * n_morph), size = round(n_morph * n_individuals * 0.15), replace = F) #introduce some rare NAs for realism
MORPH[na_indices_morph] <- NA
colnames(MORPH) <- paste("Trait", 1:n_morph) #rename columns as "Trait 1", "Trait 2", ..., "Trait n"
head(MORPH);pairs(MORPH[,1:4])


# Original monticola71 data from Pyron et al. 2023
seal_in_c90 <- read.structure("./data/seal_in_c90/seal_in_c90.str",
                    n.ind = 71,
                    n.loc = 7809,
                    onerowperind = F,
                    col.lab = 1,
                    col.pop = 0,
                    col.others = 0,
                    row.marknames = 0,
                    NA.char = -9)

# Filter
seal_in_c90 <- seal_in_c90[loc=which(seal_in_c90$loc.n.all==2),drop=TRUE]           # biallelic
seal_in_c90 <- seal_in_c90[loc=-which(minorAllele(seal_in_c90)<1/nInd(seal_in_c90))] # singletons
seal_in_c90 = missingno(seal_in_c90, type = "loci", cutoff = 0.20)                 # loci
seal_in_c90 = missingno(seal_in_c90, type = "geno", cutoff = 0.60)                 # individuals
seal_in_c90

# Make alleles object
alleles <- makefreq(seal_in_c90)


## Set number of clusters and individuals
n_clusters <- 3  #number of clusters
n_k3_test <- 5  #number of traits


## Assign each individual to cluster
clusters <- sample(1:n_clusters, n_individuals, replace = T)


## Simulate k3 data for each cluster
k3_test <- matrix(NA, nrow = n_individuals, ncol = n_k3_test)
colnames(k3_test) <- paste("Trait", 1:n_k3_test)
rownames(k3_test) <- paste0("ID", 1:n_individuals)

# Define cluster-specific means for k3 traits
k3_test_means <- list(
  cluster_1 = rnorm(n_k3_test, mean = 10, sd = 4),  
  cluster_2 = rnorm(n_k3_test, mean = 40, sd = 4),  
  cluster_3 = rnorm(n_k3_test, mean = 70, sd = 4)
)

# Assign data for each cluster with some random variation
for (i in 1:n_clusters) {
  cluster_indices <- which(clusters == i)
  k3_test[cluster_indices, ] <- matrix(rnorm(length(cluster_indices) * n_k3_test, 
                                             mean = k3_test_means[[paste0("cluster_", i)]], 
                                             sd = 3), 
                                       nrow = length(cluster_indices), 
                                       ncol = n_k3_test)
}


## Introduce some missing values
na_indices_k3_test <- sample(1:(n_individuals * n_k3_test), 
                             size = round(n_k3_test * n_individuals * 0.1), 
                             replace = F)
k3_test[na_indices_k3_test] <- NA
k3_test <- as.data.frame(k3_test)


## Plot histograms for trait 1
hist(k3_test$`Trait 1`, 
     breaks = 30, 
     main = "Histogram of Trait 1", 
     xlab = "Trait 1")


## Simulate coordinate data
GER_Coordinates <- data.frame(
  Individual = rownames(MORPH),
  Latitude = runif(n = n_individuals, min = 49.0, max = 53.0),
  Longitude = runif(n = n_individuals, min = 7.5, max = 14.5)
)

US_Coordinates <- data.frame(
  Individual = rownames(MORPH),
  Latitude = runif(n = n_individuals, min = 30, max = 40),
  Longitude = runif(n = n_individuals, min = -90, max = -80)
)



#### Test runs of functions

## Run SOMS for various datasets
SOM_single <- run.SOM(alleles)
SOM_single2 <- run.SOM(MORPH)
SOM_single3 <- run.SOM(ENV)
SOM_single4 <- run.SOM(k3_test)
SOM_multi <- run.SOM(list(MORPH, ENV))
SOM_multi2 <- run.SOM(list(MORPH, ENV, Alleles))

## Plot their learning
plot.Learning.SOM(SOM_single)
plot.Learning.SOM(SOM_single2)
plot.Learning.SOM(SOM_single3)
plot.Learning.SOM(SOM_single4)
plot.Learning.SOM(SOM_multi)
plot.Learning.SOM(SOM_multi2)

## Layer weights
plot.Layers.SOM(SOM_single)
plot.Layers.SOM(SOM_single2)
plot.Layers.SOM(SOM_single3)
plot.Layers.SOM(SOM_single4)
plot.Layers.SOM(SOM_multi)
plot.Layers.SOM(SOM_multi2)

## K values
plot.K.SOM(SOM_single)
plot.K.SOM(SOM_single2)
plot.K.SOM(SOM_single3)
plot.K.SOM(SOM_single4)
plot.K.SOM(SOM_multi)
plot.K.SOM(SOM_multi2)

## Plot models
plot.Model.SOM(SOM_single)
plot.Model.SOM(SOM_single2)
plot.Model.SOM(SOM_single3)
plot.Model.SOM(SOM_single4)
plot.Model.SOM(SOM_multi)
plot.Model.SOM(SOM_multi2)

## Plot species coefficients
plot.Structure.SOM(SOM_single)
plot.Structure.SOM(SOM_single2)
plot.Structure.SOM(SOM_single3)
plot.Structure.SOM(SOM_single4)
plot.Structure.SOM(SOM_multi)
plot.Structure.SOM(SOM_multi2)

## Plot maps
plot.Map.SOM(SOM_single3, GER_Coordinates)
plot.Map.SOM(SOM_single2, US_Coordinates)

## Plot variable importance
plot.Variable.Importance.SOM(SOM_single)
plot.Variable.Importance.SOM(SOM_single2)#Something is messed up here
plot.Variable.Importance.SOM(SOM_single3)#Something is messed up here
plot.Variable.Importance.SOM(SOM_single4)#Something is messed up here
plot.Variable.Importance.SOM(SOM_multi)
plot.Variable.Importance.SOM(SOM_multi2)


##################################################
### Empirical testing datasets - Desmognathus  ###
### These datasets are newly generated updates ###
### of previously published GBS assemblies,    ###
### with space, climate, and traits added      ###
##################################################

## monticola71 - original GBS data from Pyron et al. 2023 with updated climatic variables
# Read in sample data
monticola71data <- read.csv("./data/monticola71/monticola71.csv",
                            row.names=1,header=T,colClasses = c(huc2 = "character",
                                                                huc4 = "character",
                                                                huc6 = "character",
                                                                huc8 = "character",
                                                                huc10 = "character",
                                                                huc12 = "character"))

# Read in VCF
monticola71vcf <- read.vcfR("./data/monticola71/monticola71.vcf.gz")

# Convert to genind
genind_obj <- vcfR2genind(monticola71vcf);genind_obj

# Filter
genind_obj <- genind_obj[loc=which(genind_obj$loc.n.all==2),drop=TRUE]           # biallelic
genind_obj <- genind_obj[loc=-which(minorAllele(genind_obj)<1/nInd(genind_obj))] # singletons
genind_obj = missingno(genind_obj, type = "loci", cutoff = 0.20)                 # loci
genind_obj = missingno(genind_obj, type = "geno", cutoff = 0.60)                 # individuals
genind_obj

# Get allele frequency matrix
monticola71alleles <- makefreq(genind_obj, missing = NA)  # Or missing = 0/1 depending on how you want to treat NAs

# Rename and reorder alleles
rownames(monticola71alleles) <- monticola71data$Sample[match(rownames(monticola71alleles),rownames(monticola71data))]
monticola71alleles <- monticola71alleles[monticola71data$Sample,]

# Separate layers
  #Space
  vars <- c("Lat", "Long", "Elev", "lvl4", "huc2", "huc4", "huc6", "huc8", "huc10", "huc12")  # Select the columns
  factor_vars <- c("lvl4", "huc2", "huc4", "huc6", "huc8", "huc10", "huc12")  # One-hot encode the factor variables (lvl4 and all huc*)
  monticola71data[factor_vars] <- lapply(monticola71data[factor_vars], factor)  # Convert to factors (if not already)
  oh_encoded <- model.matrix(~ . - 1, data = monticola71data[factor_vars])  # Use model.matrix to one-hot encode (removes intercept column)
  numeric_data <- monticola71data[, c("Lat", "Long", "Elev")]  # Combine with numeric spatial data
  monticola71space <- cbind(numeric_data, oh_encoded) # Combine variables
  rownames(monticola71space) <- monticola71data$Sample  # Set row names
  monticola71space <- as.matrix(monticola71space)  # Optional: convert to matrix if needed (e.g., for clustering)

  #Climate
  climate_vars <- grep("^(wc|current)", names(monticola71data), value = TRUE)  # Select only columns starting with "wc" or "current"
  monticola71climate <- as.matrix(monticola71data[, climate_vars])  # Subset and convert to matrix
  rownames(monticola71climate) <- monticola71data$Sample  # Set row names

  #Traits
  trait_names <- c("SVL", "TL", "AG", "CW", "FL", "HL", "SG", "TW", "TO", "FI", "HW", "ED", "IN", "ES", "ON", "IO", "IC")   # Step 1: Subset the traits from the data
  trait_data <- monticola71data[, trait_names] # Extract the variables
  log_traits <- log(trait_data)  # Step 2: Log-transform all traits
  SVL <- log_traits[, "SVL"]  # Step 3: Extract SVL and residualize all others
  residuals_mat <- sapply(trait_names[-1], function(trait) {resid(lm(log_traits[, trait] ~ SVL))})  # Regress each trait on SVL and store residuals
  rownames(residuals_mat) <- monticola71data$Sample  # Step 4: Set row and column names
  monticola71traits <- cbind(SVL = SVL, residuals_mat)  # Step 5: Combine log(SVL) and residuals

# Remove multicollinearity
monticola71xyz <- monticola71space[,rownames(corSelect(monticola71space, sp.cols=NULL, var.cols=1:40, select="VIF", cor.th=0.7)$remaining.multicollinearity)]
monticola71clim <- monticola71climate[,rownames(corSelect(monticola71climate, sp.cols=NULL, var.cols=1:36, select="VIF", cor.th=0.7)$remaining.multicollinearity)]
monticola71pheno <- monticola71traits[,rownames(corSelect(monticola71traits, sp.cols=NULL, var.cols=1:17, select="VIF", cor.th=0.7)$remaining.multicollinearity)]
  
#SOM
SOM_monticola71 <- run.SOM(list(monticola71alleles,
                                monticola71space,
                                monticola71climate,
                                monticola71traits), #can be one matrix/dataframe or multiple matrices/dataframes provided as list()
                           N.steps = 100, #number of training iterations for SOM
                           N.replicates = 100, #number of SOM runs
                           max.k = 5, #maximum of considered clusters K + 1
                           set.k = NULL, #used to test a single value of K
                           BIC.thresh = 2, #BIC threshold for selecting K>1 - we suggest using Raftery (1995) ranges: 2, 6, or 10 for weak, medium, or strong support
                           learning.rate.initial = 0.6, #initial learning rate for SOM training
                           learning.rate.final = 0.2, #final learning rate for SOM training
                           max.NA = 0.9, #maximum fraction of missing values allowed per row in input data to prevent row to be removed
                           random.starts.kmeans = 25, #number of random starts for k-means clustering
                           training.neighborhoods = "gaussian", #neighborhood function used for SOM training SOM (options; "gaussian" or "bubble")
                           save.SOM.results = F, #whether to save SOM results to file
                           save.SOM.results.name = NULL, #file name for saving SOM results (if NULL, default name based on input_data is generated)
                           overwrite.SOM.results = T, #if FALSE, existing results are loaded instead of re-running SOM
                           message.N.replicates = 1) #frequency of progress messages during training (message is printed every message.N.replicates iterations)
  
#Plot results
plot.Learning.SOM(SOM_monticola71)
plot.Layers.SOM(SOM_monticola71)
plot.K.SOM(SOM_monticola71)
plot.Model.SOM(SOM_monticola71)
plot.Structure.SOM(SOM_monticola71)
plot.Map.SOM(SOM_monticola71, data.frame(Individual=monticola71data$Sample,
                                         Latitude=monticola71data$Lat,
                                         Longitude=monticola71data$Long))
plot.Variable.Importance.SOM(SOM_monticola71)
