################################################################################
#### Set environment
################################################################################

rm(list = ls()) #clear environment
setwd("./")
#setwd("C:/Users/danie/Desktop/PhD research/SOM package")
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")



################################################################################
#### SOM speed test
################################################################################

## Set SOM parameters
N.replicates <- 110
N.steps <- 150
grid.size <- c(4, 3)

## Create function to generate example allele datasets
generate.SOM.example.data <- function(n.samples = 50,
                                      n.variables,
                                      missing.fraction = 0,
                                      seed = 1)
{
  set.seed(seed)
    x <- matrix(rbinom(n.samples * n.variables, size = 2, prob = 0.5),
              nrow = n.samples,
              ncol = n.variables)
  row.names(x) <- paste0("sample_", seq_len(n.samples))
  colnames(x) <- paste0("var_", seq_len(n.variables))
  if(missing.fraction > 0)
  {
    na.idx <- sample(seq_len(length(x)),
                     size = round(length(x) * missing.fraction),
                     replace = FALSE)
    x[na.idx] <- NA
  }
  return(as.data.frame(x))
}


## Generate testing data with different numbers of variables
data_10 <- generate.SOM.example.data(n.variables = 10)
data_100 <- generate.SOM.example.data(n.variables = 100)
data_1000 <- generate.SOM.example.data(n.variables = 1000)
data_10000 <- generate.SOM.example.data(n.variables = 10000)
data_100000 <- generate.SOM.example.data(n.variables = 100000)


## Create function to train SOM using the original DNA.SOM structure with clustering steps removed
DNA.SOM.training.only <- function(alleles,
                                  grid.size,
                                  m,
                                  n)
{
  som_grid <- kohonen::somgrid(xdim = grid.size[1],
                      ydim = grid.size[2],
                      topo = "hexagonal")
  l_mat <- data.frame(row.names = 1:m) #to hold learning values
  for(j in 1:n) {
    som_model <- kohonen::som(as.matrix(alleles),
                     grid = som_grid,
                     maxNA.fraction = 0.9,
                     alpha = c(0.5, 0.1),
                     rlen = m)
    l_mat[, j] <- som_model$changes
  }
  return(list(l_mat = l_mat,
              som_model = som_model))
}


## Test: original implementation
time_10_original <- system.time(test_10 <- DNA.SOM.training.only(data_10,
                                                        grid.size = grid.size,
                                                        m = N.steps,
                                                        n = N.replicates))
time_100_original <- system.time(test_100 <- DNA.SOM.training.only(data_100,
                                                                 grid.size = grid.size,
                                                                 m = N.steps,
                                                                 n = N.replicates))
time_1000_original <- system.time(test_1000 <- DNA.SOM.training.only(data_1000,
                                                                   grid.size = grid.size,
                                                                   m = N.steps,
                                                                   n = N.replicates))
time_10000_original <- system.time(test_10000 <- DNA.SOM.training.only(data_10000,
                                                                    grid.size = grid.size,
                                                                    m = N.steps,
                                                                    n = N.replicates))
time_100000_original <- system.time(test_100000 <- DNA.SOM.training.only(data_100000,
                                                                       grid.size = grid.size,
                                                                       m = N.steps,
                                                                       n = N.replicates))


print(time_10_original[[3]])
print(time_100_original[[3]])
print(time_1000_original[[3]])
print(time_10000_original[[3]])
print(time_100000_original[[3]])


## Test: new implementation in parallel with 4 cores
time_10_new_parallel <- system.time(test_10 <- train.SOM(input_data = data_10,
                                                         verbose = FALSE,
                                                                grid.size = grid.size,
                                                                 N.steps = N.steps,
                                                         parallel = TRUE, 
                                                         max.NA.row = 0.9, 
                                                         max.NA.col = 0.9,
                                                         N.cores = 4,
                                                                 N.replicates = N.replicates))
time_100_new_parallel <- system.time(test_100 <- train.SOM(input_data = data_100,
                                                         verbose = FALSE,
                                                         grid.size = grid.size,
                                                         N.steps = N.steps,
                                                         parallel = TRUE, 
                                                         max.NA.row = 0.9, 
                                                         max.NA.col = 0.9,
                                                         N.cores = 4,
                                                         N.replicates = N.replicates))
time_1000_new_parallel <- system.time(test_1000 <- train.SOM(input_data = data_1000,
                                                           verbose = FALSE,
                                                           grid.size = grid.size,
                                                           N.steps = N.steps,
                                                           parallel = TRUE, 
                                                           max.NA.row = 0.9, 
                                                           max.NA.col = 0.9,
                                                           N.cores = 4,
                                                           N.replicates = N.replicates))
time_10000_new_parallel <- system.time(test_10000 <- train.SOM(input_data = data_10000,
                                                             verbose = FALSE,
                                                             grid.size = grid.size,
                                                             N.steps = N.steps,
                                                             parallel = TRUE, 
                                                             max.NA.row = 0.9, 
                                                             max.NA.col = 0.9,
                                                             N.cores = 4,
                                                             N.replicates = N.replicates))
time_100000_new_parallel <- system.time(test_100000 <- train.SOM(input_data = data_100000,
                                                               verbose = FALSE,
                                                               grid.size = grid.size,
                                                               N.steps = N.steps,
                                                               parallel = TRUE, 
                                                               max.NA.row = 0.9, 
                                                               max.NA.col = 0.9,
                                                               N.cores = 4,
                                                               N.replicates = N.replicates))




time_10_new_parallel[[3]]
time_100_new_parallel[[3]]
time_1000_new_parallel[[3]]
time_10000_new_parallel[[3]]
time_100000_new_parallel[[3]]


## Test: new implementation not in parallel
time_10_not_parallel <- system.time(test_10 <- train.SOM(input_data = data_10,
                                                         verbose = FALSE,
                                                         grid.size = grid.size,
                                                         N.steps = N.steps,
                                                         parallel = FALSE, 
                                                         max.NA.row = 0.9, 
                                                         max.NA.col = 0.9,
                                                         N.cores = 4,
                                                         N.replicates = N.replicates))
time_100_not_parallel <- system.time(test_100 <- train.SOM(input_data = data_100,
                                                           verbose = FALSE,
                                                           grid.size = grid.size,
                                                           N.steps = N.steps,
                                                           parallel = FALSE, 
                                                           max.NA.row = 0.9, 
                                                           max.NA.col = 0.9,
                                                           N.cores = 4,
                                                           N.replicates = N.replicates))
time_1000_not_parallel <- system.time(test_1000 <- train.SOM(input_data = data_1000,
                                                             verbose = FALSE,
                                                             grid.size = grid.size,
                                                             N.steps = N.steps,
                                                             parallel = FALSE, 
                                                             max.NA.row = 0.9, 
                                                             max.NA.col = 0.9,
                                                             N.cores = 4,
                                                             N.replicates = N.replicates))
time_10000_not_parallel <- system.time(test_10000 <- train.SOM(input_data = data_10000,
                                                               verbose = FALSE,
                                                               grid.size = grid.size,
                                                               N.steps = N.steps,
                                                               parallel = FALSE, 
                                                               max.NA.row = 0.9, 
                                                               max.NA.col = 0.9,
                                                               N.cores = 4,
                                                               N.replicates = N.replicates))
time_100000_not_parallel <- system.time(test_100000 <- train.SOM(input_data = data_100000,
                                                                 verbose = FALSE,
                                                                 grid.size = grid.size,
                                                                 N.steps = N.steps,
                                                                 parallel = FALSE, 
                                                                 max.NA.row = 0.9, 
                                                                 max.NA.col = 0.9,
                                                                 N.cores = 4,
                                                                 N.replicates = N.replicates))




time_10_not_parallel[[3]]
time_100_not_parallel[[3]]
time_1000_not_parallel[[3]]
time_10000_not_parallel[[3]]
time_100000_not_parallel[[3]]
