################################################################################
#### Set environment
################################################################################

rm(list = ls()) #clear environment
setwd("C:/Users/danie/Desktop/PhD research/SOM package")
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")


#################################################################################
#### Create function to simulate simple test data
################################################################################

## Create function to simulate two layers with user-defined numbers of clusters and mean cluster distances
simulate.simple.layers.SOM <- function(N.samples = 40, #number of samples
                                       N.variables.layer1 = 20, #number of variables in layer 1
                                       N.variables.layer2 = 20, #number of variables in layer 2
                                       layer1.K = 3, #number of clusters in layer 1
                                       layer2.K = 3, #number of clusters in layer 2
                                       layer1.mean.cluster.distance = 2, #pairwise Euclidean distance between cluster centers in layer 1
                                       layer2.mean.cluster.distance = 2, #pairwise Euclidean distance between cluster centers in layer 2
                                       layer1.within.cluster.sd = 1, #within-cluster sd in layer 1
                                       layer2.within.cluster.sd = 1, #within-cluster sd in layer 2
                                       prop.differentiated.variables.layer1 = 0.8, #proportion of variables with cluster mean differences in layer 1
                                       prop.differentiated.variables.layer2 = 0.8, #proportion of variables with cluster mean differences in layer 2
                                       between.layer.structure = c("independent",
                                                                   "layer2_nested_within_layer1",
                                                                   "layer1_nested_within_layer2"), #relationship between layer-specific sample memberships
                                       seed = 1 #random seed
) {
  if (!is.numeric(N.samples) || length(N.samples) != 1 || N.samples < 2 || N.samples %% 1 != 0) {
    stop("'N.samples' must be a single integer >= 2.")
  }
  if (!is.numeric(N.variables.layer1) || length(N.variables.layer1) != 1 || N.variables.layer1 < 1 || N.variables.layer1 %% 1 != 0) {
    stop("'N.variables.layer1' must be a single integer >= 1.")
  }
  if (!is.numeric(N.variables.layer2) || length(N.variables.layer2) != 1 || N.variables.layer2 < 1 || N.variables.layer2 %% 1 != 0) {
    stop("'N.variables.layer2' must be a single integer >= 1.")
  }
  if (!is.numeric(layer1.K) || length(layer1.K) != 1 || layer1.K < 1 || layer1.K %% 1 != 0) {
    stop("'layer1.K' must be a single integer >= 1.")
  }
  if (!is.numeric(layer2.K) || length(layer2.K) != 1 || layer2.K < 1 || layer2.K %% 1 != 0) {
    stop("'layer2.K' must be a single integer >= 1.")
  }
  if (layer1.K > N.samples) {
    stop("'layer1.K' cannot be larger than 'N.samples'.")
  }
  if (layer2.K > N.samples) {
    stop("'layer2.K' cannot be larger than 'N.samples'.")
  }
  if (!is.numeric(layer1.mean.cluster.distance) || length(layer1.mean.cluster.distance) != 1 || layer1.mean.cluster.distance < 0) {
    stop("'layer1.mean.cluster.distance' must be a single numeric value >= 0.")
  }
  if (!is.numeric(layer2.mean.cluster.distance) || length(layer2.mean.cluster.distance) != 1 || layer2.mean.cluster.distance < 0) {
    stop("'layer2.mean.cluster.distance' must be a single numeric value >= 0.")
  }
  if (!is.numeric(layer1.within.cluster.sd) || length(layer1.within.cluster.sd) != 1 || layer1.within.cluster.sd < 0) {
    stop("'layer1.within.cluster.sd' must be a single numeric value >= 0.")
  }
  if (!is.numeric(layer2.within.cluster.sd) || length(layer2.within.cluster.sd) != 1 || layer2.within.cluster.sd < 0) {
    stop("'layer2.within.cluster.sd' must be a single numeric value >= 0.")
  }
  if (!is.numeric(prop.differentiated.variables.layer1) || length(prop.differentiated.variables.layer1) != 1 || prop.differentiated.variables.layer1 <= 0 || prop.differentiated.variables.layer1 > 1) {
    stop("'prop.differentiated.variables.layer1' must be a single numeric value > 0 and <= 1.")
  }
  if (!is.numeric(prop.differentiated.variables.layer2) || length(prop.differentiated.variables.layer2) != 1 || prop.differentiated.variables.layer2 <= 0 || prop.differentiated.variables.layer2 > 1) {
    stop("'prop.differentiated.variables.layer2' must be a single numeric value > 0 and <= 1.")
  }
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("'seed' must be a single numeric value.")
  }
  between.layer.structure <- match.arg(between.layer.structure)
  if (between.layer.structure == "layer2_nested_within_layer1" && layer2.K < layer1.K) {
    stop("For 'layer2_nested_within_layer1', 'layer2.K' must be >= 'layer1.K'.")
  }
  if (between.layer.structure == "layer1_nested_within_layer2" && layer1.K < layer2.K) {
    stop("For 'layer1_nested_within_layer2', 'layer1.K' must be >= 'layer2.K'.")
  }
  base::set.seed(seed)
  sample.names <- paste0("Sample_", seq_len(N.samples))
  variable.names.layer1 <- paste0("Layer1_Var_", seq_len(N.variables.layer1))
  variable.names.layer2 <- paste0("Layer2_Var_", seq_len(N.variables.layer2))
  
  ## Create function to generate balanced group labels for one layer
  create.balanced.group.labels <- function(N.samples,
                                           N.clusters,
                                           sample.names,
                                           layer.name) {
    cluster.sizes <- rep(floor(N.samples / N.clusters), N.clusters)
    remainder.samples <- N.samples %% N.clusters
    if (remainder.samples > 0) {
      cluster.sizes[seq_len(remainder.samples)] <- cluster.sizes[seq_len(remainder.samples)] + 1
    }
    group.labels <- rep(paste0(layer.name, "_Group_", seq_len(N.clusters)), times = cluster.sizes)
    group.labels <- sample(group.labels, size = length(group.labels), replace = FALSE) #randomize sample membership to avoid order-driven overlap between layers
    group.labels <- factor(group.labels, levels = paste0(layer.name, "_Group_", seq_len(N.clusters)))
    names(group.labels) <- sample.names
    return(group.labels)
  }
  
  ## Create function to generate nested group labels for one layer within another
  create.nested.group.labels <- function(parent.labels,
                                         child.K,
                                         child.layer.name) {
    parent.indices <- split(seq_along(parent.labels), parent.labels)
    parent.K <- length(parent.indices)
    if (child.K < parent.K) {
      stop("'child.K' must be >= number of parent groups for strict nesting.")
    }
    child.groups.per.parent <- rep(floor(child.K / parent.K), parent.K)
    remainder.child.groups <- child.K %% parent.K
    if (remainder.child.groups > 0) {
      child.groups.per.parent[seq_len(remainder.child.groups)] <- child.groups.per.parent[seq_len(remainder.child.groups)] + 1
    }
    child.group.ids <- split(seq_len(child.K), rep(seq_len(parent.K), times = child.groups.per.parent))
    child.labels <- character(length(parent.labels))
    for (parent.index in seq_along(parent.indices)) {
      current.indices <- parent.indices[[parent.index]]
      current.child.ids <- child.group.ids[[parent.index]]
      if (length(current.indices) < length(current.child.ids)) {
        stop(paste0("Cannot nest ", child.K, " child groups within the parent groups because at least one parent group has fewer samples than assigned child groups."))
      }
      current.child.sizes <- rep(floor(length(current.indices) / length(current.child.ids)), length(current.child.ids))
      remainder.samples <- length(current.indices) %% length(current.child.ids)
      if (remainder.samples > 0) {
        current.child.sizes[seq_len(remainder.samples)] <- current.child.sizes[seq_len(remainder.samples)] + 1
      }
      current.child.labels <- rep(paste0(child.layer.name, "_Group_", current.child.ids), times = current.child.sizes)
      current.child.labels <- sample(current.child.labels, size = length(current.child.labels), replace = FALSE) #randomize assignment within each parent group
      child.labels[current.indices] <- current.child.labels
    }
    child.labels <- factor(child.labels, levels = paste0(child.layer.name, "_Group_", seq_len(child.K)))
    names(child.labels) <- names(parent.labels)
    return(child.labels)
  }
  
  create.regular.simplex.coordinates <- function(N.clusters) {
    if (N.clusters == 1) {
      return(matrix(0, nrow = 1, ncol = 1))
    }
    identity.matrix <- diag(N.clusters)
    centering.matrix <- matrix(1 / N.clusters, nrow = N.clusters, ncol = N.clusters)
    simplex.matrix <- identity.matrix - centering.matrix
    eigen.results <- eigen(crossprod(simplex.matrix), symmetric = TRUE)
    non.zero.eigenvalue.indices <- which(eigen.results$values > 1e-10)
    simplex.coordinates <- simplex.matrix %*% eigen.results$vectors[, non.zero.eigenvalue.indices, drop = FALSE] %*% diag(1 / sqrt(eigen.results$values[non.zero.eigenvalue.indices]), nrow = length(non.zero.eigenvalue.indices))
    return(simplex.coordinates)
  }
  
  generate.cluster.centers <- function(N.clusters,
                                       N.variables,
                                       mean.cluster.distance,
                                       prop.differentiated.variables,
                                       layer.name) {
    if (N.clusters == 1) {
      cluster.center.matrix <- matrix(0, nrow = 1, ncol = N.variables)
      pairwise.cluster.distances <- numeric(0)
      differentiated.variable.indices <- integer(0)
      return(list(cluster.center.matrix = cluster.center.matrix,
                  pairwise.cluster.distances = pairwise.cluster.distances,
                  differentiated.variable.indices = differentiated.variable.indices))
    }
    if (mean.cluster.distance == 0) {
      cluster.center.matrix <- matrix(0, nrow = N.clusters, ncol = N.variables)
      pairwise.cluster.distances <- as.vector(stats::dist(cluster.center.matrix))
      differentiated.variable.indices <- integer(0)
      return(list(cluster.center.matrix = cluster.center.matrix,
                  pairwise.cluster.distances = pairwise.cluster.distances,
                  differentiated.variable.indices = differentiated.variable.indices))
    }
    N.differentiated.variables <- max(1, round(N.variables * prop.differentiated.variables))
    if (N.differentiated.variables < (N.clusters - 1)) {
      stop(paste0("For ", layer.name, ", the number of differentiated variables must be >= N.clusters - 1 to generate exact equal pairwise cluster distances."))
    }
    differentiated.variable.indices <- sort(sample(seq_len(N.variables), size = N.differentiated.variables, replace = FALSE))
    simplex.coordinates <- create.regular.simplex.coordinates(N.clusters = N.clusters)
    current.pairwise.distance <- as.numeric(stats::dist(simplex.coordinates))[1]
    scaled.simplex.coordinates <- simplex.coordinates * ((mean.cluster.distance * sqrt(N.differentiated.variables)) / current.pairwise.distance)
    orthonormal.loading.matrix <- qr.Q(qr(matrix(stats::rnorm(N.differentiated.variables * (N.clusters - 1)),
                                                 nrow = N.differentiated.variables,
                                                 ncol = N.clusters - 1)))
    differentiated.cluster.centers <- scaled.simplex.coordinates %*% t(orthonormal.loading.matrix)
    cluster.center.matrix <- matrix(0, nrow = N.clusters, ncol = N.variables)
    cluster.center.matrix[, differentiated.variable.indices] <- differentiated.cluster.centers
    pairwise.cluster.distances <- as.vector(stats::dist(cluster.center.matrix))
    return(list(cluster.center.matrix = cluster.center.matrix,
                pairwise.cluster.distances = pairwise.cluster.distances,
                differentiated.variable.indices = differentiated.variable.indices))
  }
  
  simulate.single.layer <- function(group.labels,
                                    N.variables,
                                    mean.cluster.distance,
                                    within.cluster.sd,
                                    prop.differentiated.variables,
                                    variable.names,
                                    layer.name) {
    N.samples <- length(group.labels)
    N.clusters <- nlevels(group.labels)
    cluster.sizes <- as.integer(table(group.labels))
    names(cluster.sizes) <- levels(group.labels)
    cluster.center.results <- generate.cluster.centers(N.clusters = N.clusters,
                                                       N.variables = N.variables,
                                                       mean.cluster.distance = mean.cluster.distance,
                                                       prop.differentiated.variables = prop.differentiated.variables,
                                                       layer.name = layer.name)
    cluster.center.matrix <- cluster.center.results$cluster.center.matrix
    pairwise.cluster.distances <- cluster.center.results$pairwise.cluster.distances
    differentiated.variable.indices <- cluster.center.results$differentiated.variable.indices
    layer.matrix <- matrix(NA_real_, nrow = N.samples, ncol = N.variables)
    for (sample.index in seq_len(N.samples)) {
      current.cluster.index <- as.integer(group.labels[sample.index])
      current.cluster.center <- cluster.center.matrix[current.cluster.index, ]
      layer.matrix[sample.index, ] <- stats::rnorm(N.variables,
                                                   mean = current.cluster.center,
                                                   sd = within.cluster.sd)
    }
    rownames(layer.matrix) <- names(group.labels)
    colnames(layer.matrix) <- variable.names
    return(list(layer.matrix = layer.matrix,
                group.labels = group.labels,
                cluster.center.matrix = cluster.center.matrix,
                cluster.sizes = cluster.sizes,
                pairwise.cluster.distances = pairwise.cluster.distances,
                differentiated.variable.indices = differentiated.variable.indices,
                realized.mean.cluster.distance = ifelse(length(pairwise.cluster.distances) > 0, mean(pairwise.cluster.distances), NA_real_)))
  }
  
  ## Create group labels for both layers with explicit between-layer structure
  if (between.layer.structure == "independent") {
    layer1.group.labels <- create.balanced.group.labels(N.samples = N.samples,
                                                        N.clusters = layer1.K,
                                                        sample.names = sample.names,
                                                        layer.name = "Layer1")
    layer2.group.labels <- create.balanced.group.labels(N.samples = N.samples,
                                                        N.clusters = layer2.K,
                                                        sample.names = sample.names,
                                                        layer.name = "Layer2")
  }
  
  if (between.layer.structure == "layer2_nested_within_layer1") {
    layer1.group.labels <- create.balanced.group.labels(N.samples = N.samples,
                                                        N.clusters = layer1.K,
                                                        sample.names = sample.names,
                                                        layer.name = "Layer1")
    layer2.group.labels <- create.nested.group.labels(parent.labels = layer1.group.labels,
                                                      child.K = layer2.K,
                                                      child.layer.name = "Layer2")
  }
  
  if (between.layer.structure == "layer1_nested_within_layer2") {
    layer2.group.labels <- create.balanced.group.labels(N.samples = N.samples,
                                                        N.clusters = layer2.K,
                                                        sample.names = sample.names,
                                                        layer.name = "Layer2")
    layer1.group.labels <- create.nested.group.labels(parent.labels = layer2.group.labels,
                                                      child.K = layer1.K,
                                                      child.layer.name = "Layer1")
  }
  
  layer1.results <- simulate.single.layer(group.labels = layer1.group.labels,
                                          N.variables = N.variables.layer1,
                                          mean.cluster.distance = layer1.mean.cluster.distance,
                                          within.cluster.sd = layer1.within.cluster.sd,
                                          prop.differentiated.variables = prop.differentiated.variables.layer1,
                                          variable.names = variable.names.layer1,
                                          layer.name = "Layer1")
  
  layer2.results <- simulate.single.layer(group.labels = layer2.group.labels,
                                          N.variables = N.variables.layer2,
                                          mean.cluster.distance = layer2.mean.cluster.distance,
                                          within.cluster.sd = layer2.within.cluster.sd,
                                          prop.differentiated.variables = prop.differentiated.variables.layer2,
                                          variable.names = variable.names.layer2,
                                          layer.name = "Layer2")
  
  return(list(layer1 = layer1.results$layer.matrix,
              layer2 = layer2.results$layer.matrix,
              layer1.group.labels = layer1.results$group.labels,
              layer2.group.labels = layer2.results$group.labels,
              layer1.cluster.centers = layer1.results$cluster.center.matrix,
              layer2.cluster.centers = layer2.results$cluster.center.matrix,
              layer1.cluster.sizes = layer1.results$cluster.sizes,
              layer2.cluster.sizes = layer2.results$cluster.sizes,
              layer1.pairwise.cluster.distances = layer1.results$pairwise.cluster.distances,
              layer2.pairwise.cluster.distances = layer2.results$pairwise.cluster.distances,
              layer1.differentiated.variable.indices = layer1.results$differentiated.variable.indices,
              layer2.differentiated.variable.indices = layer2.results$differentiated.variable.indices,
              layer1.realized.mean.cluster.distance = layer1.results$realized.mean.cluster.distance,
              layer2.realized.mean.cluster.distance = layer2.results$realized.mean.cluster.distance,
              between.layer.structure = between.layer.structure))
}



################################################################################
#### k3 test
################################################################################

## No layer difference -> result: no difference in layer importance as expected
sim.data <- simulate.simple.layers.SOM(N.variables.layer1 = 30,
                                       N.variables.layer2 = 30,
                                       layer1.mean.cluster.distance = 3,
                                       layer2.mean.cluster.distance = 3,
                                       layer1.within.cluster.sd = 1,
                                       layer2.within.cluster.sd = 1)
Layer1 <- sim.data$layer1
Layer2 <- sim.data$layer2
test <- train.SOM(input_data = list(Layer1, Layer2),
                  parallel = F, verbose = F)
test <- clustering.SOM(SOM.output = test, clustering.method = "kmeans+BICelbow")
test$optim_k_summary
plot.layer.importance.varimp.SOM(test)
plot.layer.distance.scale.SOM(test)
plot.layer.importance.leaveoneout.SOM(test)


## Layer 2 with slightly more mean separation -> results: layer 2 slightly higher importance (as expected)
sim.data <- simulate.simple.layers.SOM(N.variables.layer1 = 30,
                                       N.variables.layer2 = 30,
                                       layer1.mean.cluster.distance = 2,
                                       layer2.mean.cluster.distance = 3,
                                       layer1.within.cluster.sd = 1,
                                       layer2.within.cluster.sd = 1)
Layer1 <- sim.data$layer1
Layer2 <- sim.data$layer2
test <- train.SOM(input_data = list(Layer1, Layer2),
                  parallel = F, verbose = F)
test <- clustering.SOM(SOM.output = test, clustering.method = "kmeans+BICelbow")
test$optim_k_summary
plot.layer.importance.varimp.SOM(test)
plot.layer.distance.scale.SOM(test)
plot.layer.importance.leaveoneout.SOM(test)



## Same mean distance but higher SD in layer 2 -> results: higher importance of layer 1 as expected
sim.data <- simulate.simple.layers.SOM(N.variables.layer1 = 30,
                                       N.variables.layer2 = 30,
                                       layer1.mean.cluster.distance = 2,
                                       layer2.mean.cluster.distance = 2,
                                       layer1.within.cluster.sd = 1,
                                       layer2.within.cluster.sd = 1.6)
Layer1 <- sim.data$layer1
Layer2 <- sim.data$layer2
test <- train.SOM(input_data = list(Layer1, Layer2),
                  parallel = F, verbose = F)
test <- clustering.SOM(SOM.output = test, clustering.method = "kmeans+BICelbow")
test$optim_k_summary
plot.layer.importance.varimp.SOM(test)
plot.layer.distance.scale.SOM(test)
plot.layer.importance.leaveoneout.SOM(test)



## Same mean distance and SD but different N -> results: same layer importance
sim.data <- simulate.simple.layers.SOM(N.variables.layer1 = 30,
                                       N.variables.layer2 = 60,
                                       layer1.mean.cluster.distance = 2.5,
                                       layer2.mean.cluster.distance = 2.5,
                                       layer1.within.cluster.sd = 1,
                                       layer2.within.cluster.sd = 1)
Layer1 <- sim.data$layer1
Layer2 <- sim.data$layer2
test <- train.SOM(input_data = list(Layer1, Layer2),
                  parallel = F, verbose = F)
test <- clustering.SOM(SOM.output = test, clustering.method = "kmeans+BICelbow")
test$optim_k_summary
plot.layer.importance.varimp.SOM(test)
plot.layer.distance.scale.SOM(test)
plot.layer.importance.leaveoneout.SOM(test)


## No mean distance for Layer 2 -> results: much higher importance of layer 1 as expected
sim.data <- simulate.simple.layers.SOM(N.variables.layer1 = 30,
                                       N.variables.layer2 = 30,
                                       layer1.mean.cluster.distance = 2.5,
                                       layer2.mean.cluster.distance = 0,
                                       layer1.within.cluster.sd = 1,
                                       layer2.within.cluster.sd = 1)
Layer1 <- sim.data$layer1
Layer2 <- sim.data$layer2
test <- train.SOM(input_data = list(Layer1, Layer2),
                  parallel = F, verbose = F)
test <- clustering.SOM(SOM.output = test, clustering.method = "kmeans+BICelbow")
test$optim_k_summary
plot.layer.importance.varimp.SOM(test)
plot.layer.distance.scale.SOM(test)
plot.layer.importance.leaveoneout.SOM(test)


## No mean distance for Layer 2 and low distance for Layer 1 -> results: much higher importance of layer 1 as expected but k1 support
sim.data <- simulate.simple.layers.SOM(N.variables.layer1 = 30,
                                       N.variables.layer2 = 30,
                                       layer1.mean.cluster.distance = 1,
                                       layer2.mean.cluster.distance = 0,
                                       layer1.within.cluster.sd = 1,
                                       layer2.within.cluster.sd = 1)
Layer1 <- sim.data$layer1
Layer2 <- sim.data$layer2
test <- train.SOM(input_data = list(Layer1, Layer2),
                  parallel = F, verbose = F)
test <- clustering.SOM(SOM.output = test, clustering.method = "kmeans+BICelbow")
test$optim_k_summary
plot.layer.importance.varimp.SOM(test)
plot.layer.distance.scale.SOM(test)
plot.layer.importance.leaveoneout.SOM(test)


## No mean distance for Layer 2 but 10x more N -> results: higher importance of layer 1 as expected but no support for k3 anymore
sim.data <- simulate.simple.layers.SOM(N.variables.layer1 = 30,
                                       N.variables.layer2 = 300,
                                       layer1.mean.cluster.distance = 2,
                                       layer2.mean.cluster.distance = 0,
                                       layer1.within.cluster.sd = 1,
                                       layer2.within.cluster.sd = 1)
Layer1 <- sim.data$layer1
Layer2 <- sim.data$layer2
test <- train.SOM(input_data = list(Layer1, Layer2),
                  parallel = F, verbose = F)
test <- clustering.SOM(SOM.output = test, clustering.method = "kmeans+BICelbow")
test$optim_k_summary
plot.layer.importance.varimp.SOM(test)
plot.layer.distance.scale.SOM(test)
plot.layer.importance.leaveoneout.SOM(test)


## No mean distance for Layer 2 but 10x more N -> results: higher importance of layer 1 as expected but no support for k3 anymore
sim.data <- simulate.simple.layers.SOM(N.variables.layer1 = 30,
                                       N.variables.layer2 = 300,
                                       layer1.mean.cluster.distance = 5,
                                       layer2.mean.cluster.distance = 0,
                                       layer1.within.cluster.sd = 1,
                                       layer2.within.cluster.sd = 1)
Layer1 <- sim.data$layer1
Layer2 <- sim.data$layer2
test <- train.SOM(input_data = list(Layer1, Layer2),
                  parallel = F, verbose = F)
test <- clustering.SOM(SOM.output = test, clustering.method = "kmeans+BICelbow")
test$optim_k_summary
plot.layer.importance.varimp.SOM(test)
plot.layer.distance.scale.SOM(test)
plot.layer.importance.leaveoneout.SOM(test)



################################################################################
#### k2 test
################################################################################

## No mean distance for Layer 2 but 10x more N -> results: higher importance of layer 1 as expected but support for k1 and k2
sim.data <- simulate.simple.layers.SOM(N.variables.layer1 = 30,
                                       N.variables.layer2 = 300,
                                       layer1.mean.cluster.distance = 6,
                                       layer2.mean.cluster.distance = 0,
                                       layer1.K = 2,
                                       layer2.K = 2,
                                       layer1.within.cluster.sd = 1,
                                       layer2.within.cluster.sd = 1)
Layer1 <- sim.data$layer1
Layer2 <- sim.data$layer2
test <- train.SOM(input_data = list(Layer1, Layer2),
                  parallel = F, verbose = F)
test <- clustering.SOM(SOM.output = test, clustering.method = "kmeans+BICelbow")
test$optim_k_summary
plot.layer.importance.varimp.SOM(test)
plot.layer.distance.scale.SOM(test)
plot.layer.importance.leaveoneout.SOM(test)

Layer2.long <- data.frame(sample = rownames(Layer2), cluster = sim.data$layer2.group.labels, Layer2, check.names = F)
Layer2.long <- reshape2::melt(Layer2.long, id.vars = c("sample", "cluster"), variable.name = "variable", value.name = "value")
ggplot2::ggplot(Layer2.long[Layer2.long$variable %in% colnames(Layer2)[1:12], ], ggplot2::aes(x = value, color = cluster, fill = cluster)) + ggplot2::geom_density(alpha = 0.25) + ggplot2::facet_wrap(~variable, scales = "free", ncol = 4)



################################################################################
#### k2+k3 test
################################################################################

## Same mean + SD but k3 for layer 2 -> results: same layer importance as expected and support for k2 and k3
sim.data <- simulate.simple.layers.SOM(N.variables.layer1 = 50,
                                       N.variables.layer2 = 50,
                                       layer1.mean.cluster.distance = 5,
                                       layer2.mean.cluster.distance = 5,
                                       layer1.K = 2,
                                       layer2.K = 3,
                                       layer1.within.cluster.sd = 1,
                                       layer2.within.cluster.sd = 1,
                                       between.layer.structure = "layer2_nested_within_layer1")
Layer1 <- sim.data$layer1
Layer2 <- sim.data$layer2
test <- train.SOM(input_data = list(Layer1, Layer2),
                  parallel = F, verbose = F)
test <- clustering.SOM(SOM.output = test, clustering.method = "kmeans+BICelbow")
test$optim_k_summary
plot.layer.importance.varimp.SOM(test)
plot.layer.distance.scale.SOM(test)
plot.layer.importance.leaveoneout.SOM(test)


## Higher mean for L1 + SD but k3 for layer 2 -> results: same layer importance as expected and support for k2 and k3
sim.data <- simulate.simple.layers.SOM(N.variables.layer1 = 50,
                                       N.variables.layer2 = 50,
                                       layer1.mean.cluster.distance = 5,
                                       layer2.mean.cluster.distance = 3,
                                       layer1.K = 2,
                                       layer2.K = 3,
                                       layer1.within.cluster.sd = 1,
                                       layer2.within.cluster.sd = 1,
                                       between.layer.structure = "layer2_nested_within_layer1")
Layer1 <- sim.data$layer1
Layer2 <- sim.data$layer2
test <- train.SOM(input_data = list(Layer1, Layer2),
                  parallel = F, verbose = F)
test <- clustering.SOM(SOM.output = test, clustering.method = "kmeans+BICelbow")
test$optim_k_summary
plot.layer.importance.varimp.SOM(test)
plot.layer.distance.scale.SOM(test)
plot.layer.importance.leaveoneout.SOM(test)
