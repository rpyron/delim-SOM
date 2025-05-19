###################################################################################
#### SOM-based species delimitation (Pyron et al. 2022, 2023) v.2.0 from Daniel ###
###################################################################################

#### Install and load required R packages
required_packages <- c("adegenet", 
                       "combinat", 
                       "conStruct",
                       "kohonen", 
                       "lsr", 
                       "maps", 
                       "poppr",
                       "scales", 
                       "viridis",
                       "vcfR",
                       "fuzzySim")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = T)) install.packages(pkg)
  library(pkg, character.only = T)} #install missing and load packages

#### Functions

## Function to run single-layer SOM (one matrix) or multi-layer Super-SOM (multiple matrices)
run.SOM <- function(input_data, #can be one matrix/dataframe or multiple matrices/dataframes provided as list()
                    N.steps = 100, #number of training iterations for SOM
                    N.replicates = 20, #number of SOM runs
                    grid.size = NULL, #grid size - specify as c(x,y) if desired
                    max.k = 6, #maximum of considered clusters K + 1
                    set.k = NULL, #used to test a single value of K
                    BIC.thresh = 2, #BIC threshold for selecting K>1 - we suggest using Raftery (1995) ranges: 2, 6, or 10 for weak, medium, strong support
                    learning.rate.initial = 0.6, #initial learning rate for SOM training
                    learning.rate.final = 0.2, #final learning rate for SOM training
                    max.NA = 0.9, #maximum fraction of missing values allowed per row in input data to prevent row to be removed
                    random.starts.kmeans = 25, #number of random starts for k-means clustering
                    training.neighborhoods = "gaussian", #neighborhood function used for SOM training SOM (options; "gaussian" or "bubble")
                    save.SOM.results = F, #whether to save SOM results to file
                    save.SOM.results.name = NULL, #file name for saving SOM results (if NULL, default name based on input_data is generated)
                    overwrite.SOM.results = T, #if FALSE, existing results are loaded instead of re-running SOM
                    message.N.replicates = 1) #frequency of progress messages during training (message is printed every message.N.replicates iterations)
{
  
  ## Set save.SOM.results.name and extract input_data_names ...
  
  # ... for list with multiple data sets
  if (is.list(input_data) && length(input_data) > 1 && !is.data.frame(input_data)) { 
    input_data_names <- sapply(substitute(input_data)[-1], deparse) #extract names of each dataset
    if (is.null(save.SOM.results.name)) { #only assign default name if save.SOM.results.name is NULL
      save.SOM.results.name <- paste0("SOM_results_", paste(input_data_names, collapse = "_"), ".Rdata")
    }
    
    # ... for list with one data set
  } else if (is.list(input_data) && length(input_data) == 1) {
    input_data_names <- deparse(substitute(input_data)) #extract names of dataset
    input_data_names <- gsub("^list\\((.+)\\)$", "\\1", input_data_names) #remove "list"
    input_data_names <- gsub("\"", "", input_data_names) #remove quotes
    if (is.null(save.SOM.results.name)) { #only assign default name if save.SOM.results.name is NULL
      save.SOM.results.name <- paste0("SOM_results_", input_data_names, ".Rdata")
    }
    
    # ... for non-list object with one dataset
  } else {
    input_data_names <- deparse(substitute(input_data)) #extract names of dataset
    if (is.null(save.SOM.results.name)) { #only assign default name if save.SOM.results.name is NULL
      save.SOM.results.name <- paste0("SOM_results_", input_data_names, ".Rdata")
    }
  }
  
  ## Check and transform data if necessary ...
  
  # ... for multiple datasets
  if (is.list(input_data) && length(input_data) > 1 && !is.data.frame(input_data)) { 
    
    # Check if row names are consistent across all matrices
    if (!all(sapply(1:(length(input_data) - 1), function(i) 
      identical(rownames(input_data[[i]]), rownames(input_data[[i + 1]]))))) {
      stop("Row names are not consistent across matrices - provide datasets with matching row names")
    }
    
    # For each dataset, ...
    input_data <- lapply(seq_along(input_data), function(index) {
      mat <- input_data[[index]]
      
      # Convert to matrix if not already
      if (!is.matrix(mat)) {
        mat <- as.matrix(mat)
        message(paste(input_data_names[[index]], "dataset was converted to matrix"))
      }
      
      # Normalize matrix if not already normalized
      if (any(mat < 0, na.rm = T) || any(mat > 1, na.rm = T)) {
        #mat <- (mat - min(mat, na.rm = T)) / (max(mat, na.rm = T) - min(mat, na.rm = T))
        #message(paste(input_data_names[[index]], "dataset was normalized to range [0, 1]"))
        mat <- normalize_columnwise(mat)
        message(paste(input_data_names[[index]], "dataset was normalized to range [0, 1] by column"))
      }
      
      # Check matrices
      if (all(is.na(mat))) {
        stop("All values are missing in input dataset")
      }
      if (max(mat, na.rm = T) == min(mat, na.rm = T)) {
        stop("Data normalization failed due to problematic data")
      }
      
      
      return(mat)
    })
    
    # ... for single dataset
  } else {
    
    # If input_data is list with single dataframe, convert it to list with matrix
    if (is.list(input_data) && is.data.frame(input_data[[1]])) {
      input_data[[1]] <- as.matrix(input_data[[1]]) #convert data frame to matrix
      message(input_data_names, paste0(" dataset was converted to matrix"))
    }
    
    # Convert to matrix if not already
    if (is.data.frame(input_data)) {
      input_data <- as.matrix(input_data) #convert data frame to matrix
      message(input_data_names, paste0(" dataset was converted to matrix"))
    }
    
    # If input_data is not a list, convert it to list
    if (!is.list(input_data)) {
      input_data <- list(input_data)
    }
    
    # Normalize input_data matrix if it's not already normalized
    if (any(input_data[[1]] < 0, na.rm = T) || any(input_data[[1]] > 1, na.rm = T)) {
      #mat <- (mat - min(mat, na.rm = T)) / (max(mat, na.rm = T) - min(mat, na.rm = T))
      #message(paste(input_data_names[[index]], "dataset was normalized to range [0, 1]"))
      mat <- input_data[[1]]
      mat <- normalize_columnwise(mat)
      message(paste(input_data_names[[1]], "dataset was normalized to range [0, 1] by column"))
    }
  }
  
  # If overwrite.SOM.results is FALSE and file already exists, return saved results
  if (!overwrite.SOM.results && file.exists(save.SOM.results.name)) {
    message("SOM results already exist - loading results from file and skipping SOM run")
    load(save.SOM.results.name)
    
    # Check if structure of SOM_results is correct
    required_fields <- c("cluster_assignment", 
                         "distance_matrix", 
                         "learning_values_list", 
                         "BIC_values", 
                         "optim_k")
    if (!all(required_fields %in% names(SOM_results))) {
      stop("Loaded SOM results do not contain expected objects such as cluster_assignment, distance_matrix, etc - please check saved file or rerun SOM")
    }
    
    return(SOM_results)
  }
  
  # Create SOM output grid
  if (!is.null(grid.size)) {
    # Use user-specified grid dimensions
    SOM_output_grid <- somgrid(xdim = grid.size[1],
                               ydim = grid.size[2],
                               topo = "hexagonal",
                               neighbourhood.fct = training.neighborhoods)
  } else {
    # Use default calculation based on number of samples
    som_dim <- round(sqrt(5 * sqrt(length(rownames(input_data[[1]])))))
    SOM_output_grid <- somgrid(xdim = som_dim,
                               ydim = som_dim,
                               topo = "hexagonal",
                               neighbourhood.fct = training.neighborhoods)
  }
  
  # Function to run SOM
  replicate_som <- function(j) {
    
    # Initialize results for replicate
    cluster_assignment <- matrix(numeric(nrow(input_data[[1]])), 
                                 nrow = nrow(input_data[[1]]),
                                 ncol = N.replicates)
    d_vec <- numeric(length(input_data))
    w_vec <- numeric(max.k)
    
    # Print message every N replicates (N specified by message.N.replicates)
    if (j %% message.N.replicates == 0) {
      message(paste("Running replicate", j, "of", N.replicates))
    }
    
    # Run SOM model
    som_model <- supersom(data = input_data, 
                          grid = SOM_output_grid, 
                          maxNA.fraction = max.NA, 
                          alpha = c(learning.rate.initial, 
                                    learning.rate.final), 
                          rlen = N.steps)
    
    # Store learning values for each matrix
    learning_values_list <- lapply(seq_along(input_data), function(i) som_model$changes[, i])
    
    # Store distance weights
    d_vec <- som_model$distance.weights 
    
    # Calculate within-cluster sum of squares (WSS) for each cluster
    if (length(input_data) > 1) {
      mydata <- do.call(cbind, lapply(1:length(input_data), function(i) getCodes(som_model)[[i]]))
    } else {mydata <- getCodes(som_model)}
    rownames(mydata) <- paste0("G", seq_len(nrow(mydata)))
    
    # Get wss across k values
    wss <- numeric(max.k)
    wss[1] <- (nrow(mydata) - 1) * sum(apply(mydata, 2, var))
    wss[2:max.k] <- sapply(2:max.k, function(i) 
      sum(kmeans(mydata, 
                 centers = i, 
                 nstart = random.starts.kmeans, 
                 iter.max = 1e5)$withinss))
    
    # Calculate BIC for each cluster using WSS
    w_vec <- dim(mydata)[1] * log(wss / dim(mydata)[1]) + log(dim(mydata)[1]) * (1:max.k)
    
    # Determine optimal number of clusters using hierarchical clustering based on BIC values
    #if (which(w_vec == min(w_vec)) == 1) {
    if (!is.null(set.k)) {
      # User-specified number of clusters
      num_clusters <- set.k
    } else if (abs(diff(w_vec))[1] < BIC.thresh) {
      # If the first difference in w_vec is small, set only one cluster
      num_clusters <- 1
    } else {
      # Automatically determine number of clusters via hierarchical clustering
      temp <- cutree(hclust(dist(diff(w_vec)), method = "ward.D"), k = 2)
      goodgrp <- which.min(tapply(diff(w_vec), temp, mean))  # cluster with smallest mean difference
      num_clusters <- max(which(temp == goodgrp)) + 1
    }
    som_cluster <- cutree(hclust(dist(mydata)), num_clusters) #use hierarchical clustering to cluster codebook vectors
    cluster_gridcell_assignments <- som_cluster[som_model$unit.classif] #get vector with cluster value for each original data sample
    cluster_gridcell_assignments <- as.matrix(cluster_gridcell_assignments)
    cluster_gridcell_assignments <- cbind(Cluster = as.numeric(cluster_gridcell_assignments), Gridcell = rownames(cluster_gridcell_assignments))
    rownames(cluster_gridcell_assignments) <- rownames(input_data[[1]])
    cluster_assignment <- as.matrix(som_cluster[som_model$unit.classif]) #add assignment to classifier matrix
    
    # Return results
    return(list(cluster_assignment = cluster_assignment, 
                d_vec = d_vec, 
                learning_values_list = learning_values_list, 
                w_vec = w_vec, 
                som_model = som_model,
                som_cluster = som_cluster,
                num_clusters = num_clusters,
                cluster_gridcell_assignments = cluster_gridcell_assignments
    ))
  }
  
  # Apply SOM training and collect results for all replicates
  results <- lapply(1:N.replicates, replicate_som)
  
  # Combine results from all replicates
  cluster_assignment <- do.call(cbind, lapply(results, `[[`, "cluster_assignment"))
  rownames(cluster_assignment) <- rownames(as.matrix(input_data[[1]]))
  colnames(cluster_assignment) <- paste0("R", seq_len(ncol(cluster_assignment)))
  
  distance_matrix <- do.call(rbind, lapply(results, `[[`, "d_vec"))
  rownames(distance_matrix) <- paste0("R", seq_len(nrow(distance_matrix)))
  colnames(distance_matrix) <- as.character(input_data_names)
  
  learning_values_list <- lapply(1:length(input_data), function(i) {
    learning_values_combined <- do.call(cbind, lapply(results, function(res) res$learning_values_list[[i]]))
    rownames(learning_values_combined) <- paste0("S", seq_len(N.steps))
    colnames(learning_values_combined) <- paste0("R", seq_len(ncol(learning_values_combined)))
    return(learning_values_combined)
  })
  
  BIC_values <- do.call(cbind, lapply(results, `[[`, "w_vec"))
  rownames(BIC_values) <- paste0("k", seq_len(max.k))
  colnames(BIC_values) <- paste0("R", seq_len(ncol(BIC_values)))
  
  optim_k <- sapply(results, `[[`, "num_clusters")
  optim_k <- t(as.matrix(optim_k))
  rownames(optim_k) <- "optim_k"
  colnames(optim_k) <- paste0("R", seq_len(N.replicates))
  optim_k_mean <- mean(optim_k, na.rm = T)
  if (all(is.na(optim_k))) stop("All optimal K values are NA - check input data")
  
  
  som_models <- lapply(results, `[[`, "som_model")
  som_clusters <- lapply(results, `[[`, "som_cluster")
  cluster_gridcell_assignments <- lapply(results, `[[`, "cluster_gridcell_assignments")
  
  # Replace NA values with 0.5 in input_data
  processed_input_data <- input_data[[1]]
  if (any(is.na(processed_input_data))) {
    processed_input_data[which(is.na(processed_input_data))] <- 0.5
  }
  
  # Preprocess input data and generate cluster labels
  cluster_labels <- rbind.data.frame(lapply(1:max.k, function(x) {
    kmeans(processed_input_data , x)$cluster #run k-means for each K and store cluster labels
  }))
  rownames(cluster_labels) <- rownames(input_data) #set row names for cluster labels
  colnames(cluster_labels) <- paste("k", 1:max.k, sep = '') #set column names
  
  # Extract SOM cluster assignments and filter replicates based on maximum K
  all_k <- apply(cluster_assignment, 2, max) #calculate maximum K for each replicate (column)
  assignment_matrix <- cluster_assignment[, which(optim_k <= max.k)] #filter to keep replicates with K <= max_k
  
  # Relabel across replicates
  relabeled_assignment_matrix <- assignment_matrix #create copy of filtered matrix for relabeling
  for (i in 1:N.replicates) { 
    run.k <- max(relabeled_assignment_matrix[, i]) #get number of clusters (K) for current replicate
    run.labels <- as.numeric(cluster_labels[, run.k]) #get corresponding labels for current K
    refactor <- data.frame(row.names = rownames(cluster_labels)) #create dataframe to store relabelings
    label.perm <- permn(1:run.k) #generate all permutations of labels
    for (j in 1:length(label.perm)) {
      refactor[, j] <- relabeled_assignment_matrix[, i] #apply label permutations and calculate best matching relabeling
    }
    for (k in 1:length(label.perm)) {
      refactor[, k] <- as.numeric(permuteLevels(factor(refactor[, k]), perm = label.perm[[k]])) #relabel using permutations
    }
    assignment_matrix[, i] <- refactor[, which.min(apply(refactor, 2, function(x) {
      sum(abs(x - run.labels))}))] #choose relabeling that minimizes distance to original labels
  } 
  
  # Calculate ancestry_matrix summarizing cluster memberships and return it
  ancestry_matrix <- t(apply(assignment_matrix, 1, FUN = function(x) {
    table(factor(unlist(x), levels = 1:max(all_k))) #create table for cluster assignments and normalize it by number of replicates
  } 
  ) / dim(relabeled_assignment_matrix)[2]) 
  colnames(ancestry_matrix) <- paste0("Cluster_", seq_len(ncol(ancestry_matrix)))
  
  # Save results
  call <- match.call()
  SOM_results <- list(cluster_assignment = cluster_assignment, 
                      distance_matrix = distance_matrix, 
                      learning_values_list = learning_values_list, 
                      BIC_values = BIC_values,
                      BIC_thresh = BIC.thresh,
                      ancestry_matrix = ancestry_matrix,
                      optim_k = optim_k,
                      optim_k_mean = optim_k_mean, 
                      input_data_names = as.character(input_data_names),
                      N_steps = N.steps,
                      N_replicates = N.replicates,
                      max_k = max.k,
                      set_k = set.k,
                      learning_rate_initial = learning.rate.initial,
                      learning_rate_final = learning.rate.final,
                      som_models = som_models,
                      som_clusters = som_clusters,
                      cluster_gridcell_assignments = cluster_gridcell_assignments)
  
  # Save results
  if (save.SOM.results) {
    
    # Check if directory exists
    dir_path <- dirname(save.SOM.results.name) #extract directory path
    if (!dir.exists(dir_path)) { 
      dir.create(dir_path, recursive = T) #create directory if it doesn't exist
      message(paste("The specified directory", dir_path, "did not exist and was created"))
    }
    
    #save results
    save(SOM_results, file = save.SOM.results.name)
    
    if (save.SOM.results && !overwrite.SOM.results) {
      message("SOM results saved as ", save.SOM.results.name)
    }
    
    if (save.SOM.results && overwrite.SOM.results) {
      message("SOM results overwritten as ", save.SOM.results.name)
    }
  }
  
  return(SOM_results)
}

## Function to plot learning progress for each SOM matrix
plot.Learning.SOM <- function(SOM.output, 
                              col.pal = viridis::turbo, #set color palette
                              lines.alpha = 0.4, #transparency for plot lines
                              lines.thickness = 0.8, #thickness for plot lines
                              save = F, #option to save plot
                              overwrite = T, #option to overwrite plot if it already exists
                              resolution = 300, #plot resolution in dpi
                              type = "svg", #options: "svg", "png", "jpg"
                              file.name = NULL, #set plot file.name (if NULL, default plot file.name is used)
                              width = 20, #plot width in cm
                              height = 15, #plot heigth in cm
                              margin.bottom = 5,
                              margin.left = 5,
                              margin.top = 3,
                              margin.right = 2.5,
                              title = NULL, #main title of plot (if NULL, default title name is used)
                              legend.position = "topright", #position of legend in plot
                              legend.lines.thickness = 3, #thickness for lines in legend
                              x.axis.label = "Training steps", 
                              y.axis.label = "Learning rate change") {
  
  # Check if learning_values or learning_values_list exists and is either matrix or list
  if ("learning_values" %in% names(SOM.output)) {
    SOM.output$learning_values_list <- list(SOM.output$learning_values) #convert to list for single-layer case
  } else if ("learning_values_list" %in% names(SOM.output)) { #for multi-layer case, no changes needed
  } else {
    stop("Neither learning_values nor learning_values_list found in SOM.output")
  }
  
  # Set default file.name for saving
  if (is.null(file.name)) { 
    file.name <- paste0("SOM_learning_plot_", paste(SOM.output$input_data_names, collapse = "_"), ".", type)
  }
  
  # Check if file already exists and overwrite option is set to FALSE
  if (!overwrite && file.exists(file.name)) {
    message(sprintf("File '%s' already exists. Skipping plot saving.", file.name))
    return(invisible(NULL)) #skip saving if file exists and overwrite is FALSE
  }
  
  # Ensure that learning_values_list is a list
  if (!is.list(SOM.output$learning_values_list)) {
    stop("learning_values_list of SOM.output is not a list") #check if learning_values_list is a list
  }
  
  # Convert data.frames to matrices if necessary
  SOM.output$learning_values_list <- lapply(SOM.output$learning_values_list, function(x) {
    if (is.data.frame(x)) {
      return(as.matrix(x))
    }
    return(x)
  })
  
  # Ensure all elements in learning_values_list are matrices
  if (!all(sapply(SOM.output$learning_values_list, is.matrix))) {
    stop("learning_values_list of SOM.output contains non-matrix elements after conversion")
  }
  
  # Extract matrix names (check for multi-layer or single-layer)
  if ("input_data_names" %in% names(SOM.output)) {
    matrix_names <- SOM.output$input_data_names
  } else {
    stop("Matrix names not found in provided SOM.output")
  }
  
  # Check if matrix names match number of layers
  if (length(matrix_names) != length(SOM.output$learning_values_list)) {
    stop("Number of matrix names does not match number of matrices")
  }
  
  # Determine global ylim
  global_ylim <- range(unlist(lapply(SOM.output$learning_values_list, function(mat) 
    range(mat, na.rm = T))), na.rm = T) * c(0.93, 1.07)
  
  # Set legend and main plot title based on number of layers
  if (length(SOM.output$learning_values_list) == 1) {
    main_title <- "Training progress across layer"
    legend_title <- "Layer"
  } else {
    main_title <- "Training progress across layers"
    legend_title <- "Layers"
  }
  
  # Use provided title if specified, otherwise default
  main_title <- ifelse(is.null(title), main_title, title)
  
  # Save plot if requested
  if (save) {
    if (type == "svg") {
      svg(file.name, 
          width = width / 2.54, 
          height = height / 2.54)
    } else if (type == "png") {
      png(file.name, 
          width = width, 
          height = height, 
          units = "cm", 
          res = resolution)
    } else if (type == "jpg") {
      jpeg(file.name, 
           width = width, 
           height = height, 
           units = "cm", 
           res = resolution)
    } else {
      stop("Invalid type - choose 'svg', 'png' or 'jpg'")
    }
  }
  
  # Prepare base plot
  first_matrix <- SOM.output$learning_values_list[[1]]
  par(mfrow = c(1, 1), 
      mar = c(margin.bottom, margin.left, margin.top, margin.right))
  plot(NULL, 
       xlim = c(1, nrow(first_matrix)), 
       ylim = global_ylim, 
       xlab = x.axis.label, 
       ylab = y.axis.label, 
       main = main_title,
       axes = T)
  
  # Plot each matrix
  layer_colors <- col.pal(length(SOM.output$learning_values_list))
  for (i in seq_along(SOM.output$learning_values_list)) {
    mat <- SOM.output$learning_values_list[[i]]
    for (j in 1:ncol(mat)) 
      lines(mat[, j], 
            col = alpha(layer_colors[i], lines.alpha), 
            lwd = lines.thickness)
  }
  
  # Add legend
  legend(legend.position, 
         legend = matrix_names, 
         col = layer_colors, 
         lty = 1, 
         lwd = legend.lines.thickness,
         title = legend_title)
  
  # Close graphics device
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}

## Function to plot layer weights
plot.Layers.SOM <- function(SOM.output, 
                            col.pal = viridis::turbo, #color palette
                            save = F, #option to save plot
                            overwrite = T, #option to overwrite plot if it already exists
                            type = "svg", #plot type options: "svg", "png", "jpg"
                            file.name = NULL, #set file.name (if NULL, default plot file.name is used)
                            resolution = 300, #plot resolution in dpi
                            width = 20, #plot width in cm  
                            height = 15, #plot width in cm
                            margin.bottom = 3,
                            margin.left = 5,
                            margin.top = 2.5,
                            margin.right = 1.5,
                            title = "Layer weights", #plot title name
                            y_axis_label = "Relative weights") { #set y axis label
  
  # Check if 'distance_matrix' exists in SOM.output
  if (!"distance_matrix" %in% names(SOM.output)) {
    message("distance_matrix could not be found in SOM output")
    return(invisible(NULL)) #skip plotting
  }
  
  # If single-layer SOM is present, skip plotting
  if (length(SOM.output$input_data_names) == 1) {
    message("Single-layer SOM detected - multiple layers are needed to plot and compare layer weights")
    return(invisible(NULL)) #skip plotting
  }
  
  # Prepare data
  d.mat <- SOM.output$distance_matrix #extract distance matrix
  layer.names <- SOM.output$input_data_names #extract layer names
  
  # Validate that names match number of layers
  if (length(layer.names) != ncol(d.mat)) {
    message("Mismatch between layer names and distance matrix columns - using generic layer names")
    layer.names <- paste0("Layer ", seq_len(ncol(d.mat)))
  }
  
  # Calculate relative layer weights
  raw.weights <- sqrt(1 / ifelse(colMeans(d.mat, na.rm = T) > 0, colMeans(d.mat, na.rm = T), NA)) #compute weights
  names(raw.weights) <- layer.names #assign names to weights
  sorted.weights <- sort(raw.weights[!is.na(raw.weights)], decreasing = T) #sort weights
  
  # Set default file name for saving based on layer.names
  if (is.null(file.name)) { 
    file.name <- paste0("SOM_layers_plot_", paste(layer.names, collapse = "_"), ".", type)
  }
  
  # Check for file existence if overwrite is FALSE
  if (!overwrite && file.exists(file.name)) {
    message(file.name, " already exists - skipping plot saving")
    return(invisible(NULL)) #skip saving
  }
  
  # Define color palette for layers
  layer.cols <- setNames(col.pal(length(layer.names)), layer.names) #assign colors
  
  # Save plot if requested
  if (save) {
    
    # Set plot type
    if (type == "svg") {
      svg(file.name, 
          width = width / 2.54, 
          height = height / 2.54)
    } else if (type == "png") {
      png(file.name, 
          width = width, 
          height = height, 
          units = "cm", 
          res = resolution)
    } else if (type == "jpg") {
      jpeg(file.name, 
           width = width, 
           height = height, 
           units = "cm", 
           res = resolution)
    } else {
      stop("Invalid type - choose 'svg', 'png', or 'jpg'")
    }
  }
  
  # Set layout and margins for plotting
  par(mfrow = c(1, 1), 
      mar = c(margin.bottom, 
              margin.left, 
              margin.top, 
              margin.right))
  
  # Plot barplot
  barplot(sorted.weights, 
          main = title, 
          col = layer.cols[names(sorted.weights)], 
          names.arg = names(sorted.weights), 
          ylab = y_axis_label)
  
  # Close graphics device
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}

## Function to evaluate K-values - PREVIOUS VERSION ONLY PLOTTED IF SAVE=TRUE
plot.K.SOM <- function(SOM.output,
                       col.pal = viridis::magma, # color palette
                       save = FALSE,             # option to save plot
                       overwrite = TRUE,         # option to overwrite plot
                       type = "svg",             # plot type: "svg", "png", "jpg"
                       file.name = NULL,         # custom file name
                       resolution = 300,         # plot resolution in dpi
                       width = 10,               # plot width in cm
                       height = 15,              # plot height in cm
                       margin.bottom = 3,
                       margin.left = 5, 
                       margin.top = 2,
                       margin.right = 1,
                       title = "Number of clusters (k)") {
  
  # --- Input validation ---
  if (is.null(SOM.output$BIC_values)) {
    message("BIC_values could not be found in SOM output")
    return(invisible(NULL))
  } else if (is.null(SOM.output$N_replicates)) {
    message("N_replicates could not be found in SOM output")
    return(invisible(NULL)) 
  } else if (is.null(SOM.output$optim_k)) {
    message("optim_k could not be found in SOM output")
    return(invisible(NULL)) 
  } else if (is.null(SOM.output$max_k)) {
    message("max_k could not be found in SOM output")
    return(invisible(NULL)) 
  }
  
  # --- Generate color palette ---
  k.cols <- col.pal(SOM.output$max_k-1) 
  
  # --- Save logic: open device if needed ---
  if (save) {
    if (is.null(file.name)) {
      file.name <- paste0("SOM_K_evaluation_plot_", paste(SOM.output$input_data_names, collapse = "_"), ".", type)
    }
    
    if (!overwrite && file.exists(file.name)) {
      message(file.name, " already exists - skipping plot saving")
      return(invisible(NULL))
    }
    
    if (type == "svg") {
      svg(file.name, width = width / 2.54, height = height / 2.54)
    } else if (type == "png") {
      png(file.name, width = width, height = height, res = resolution, units = "cm")
    } else if (type == "jpg") {
      jpeg(file.name, width = width, height = height, res = resolution, units = "cm")
    } else {
      stop("Unsupported file type. Choose 'svg', 'png', or 'jpg'.")
    }
  }
  
  # --- Plotting ---
  par(mfrow = c(3, 1), 
      mar = c(margin.bottom, margin.left, margin.top, margin.right))
  
  # 1. Boxplot of BIC values
  boxplot(t(SOM.output$BIC_values)[,1:SOM.output$max_k-1], 
          outline = FALSE, 
          notch = FALSE, 
          axes = FALSE, 
          ylab = "BIC", 
          ylim = range(unlist(SOM.output$BIC_values)), 
          col = k.cols)
  axis(1, at = 1:(SOM.output$max_k-1), labels = NA)#1:(SOM.output$max_k-1))
  axis(2, at = round(range(unlist(SOM.output$BIC_values))), las = 3)
  title(title, line = 0)
  
  # 2. Boxplot of delta BIC
  d_wss <- apply(SOM.output$BIC_values, 2, function(x) diff(diff(x)))
  rownames(d_wss) <- 2:(SOM.output$max_k-1);d_wss <- rbind(NA,d_wss)
  boxplot(t(d_wss), 
          outline = FALSE, 
          notch = FALSE, 
          axes = FALSE, 
          ylab = "Delta BIC", 
          col = k.cols)
  abline(h = 0, lty = 2, col = "black")
  axis(1, at = 1:(SOM.output$max_k-1), labels = NA)#1:(SOM.output$max_k-1))
  axis(2, at = sort(c(0, round(range(unlist(na.omit(d_wss)))))), las = 3)
  
  # 3. Barplot of sampling frequency
  barplot_data <- table(factor(SOM.output$optim_k, levels = 1:(SOM.output$max_k-1))) / SOM.output$N_replicates
  bar_positions <- barplot(barplot_data, 
                           ylab = "Sampling frequency", 
                           ylim = c(0, 1), 
                           col = k.cols,
                           axes = FALSE)
  #axis(1, at = bar_positions, labels = 1:(SOM.output$max_k-1))
  axis(2, las = 3)
  
  # --- Close device if saving ---
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}

## Function to plot model results as SOM grids (showing sample assignment to cells, cell distances and boundaries between cell clusters)
plot.Model.SOM <- function(SOM.output,
                           col.pal.neighbor.dist = viridis::cividis, #color palette of neighbor distance plot (top)
                           col.pal.clusters = viridis::viridis, #color palette of cluster plot (bottom)
                           boundary.col.clusters = "red", #color of cluster boundaries (in bottom plot)
                           boundary.lwd.clusters = 3, #linewidth of cluster boundaries (in bottom plot)
                           point.col.clusters = "white", #color of sample points (in bottom plot)
                           point.shape.clusters = 19, #shape of sample points (in bottom plot)
                           point.size.clusters = 0.8, #size of sample points (in bottom plot)
                           cluster.shape.clusters = "straight", #shape ("straight" or "round") of cluster cells (in bottom plot)
                           cluster.shape.neighbor.dist = "straight", #shape ("straight" or "round") of cluster cells (in top plot)
                           save = F, #option to save plot
                           overwrite = T, #option to overwrite plot if already present
                           type = "svg", #options: "svg", "png", "jpg"
                           file.name = NULL, #set plot file name (if NULL, default name is used)
                           width = 10, #plot width in cm
                           height = 15, #plot height in cm
                           resolution = 300, #plot resolution in dpi
                           margin.top = 1,
                           margin.bottom = 0,
                           margin.left = 0,
                           margin.right = 0,
                           shift.plot.clusters = 0.14, #shift bottom plot slightly to the right to align with gridraster of top plot
                           title.clusters = "SOM clusters", #title of cluster plot (bottom)
                           title.neighbor.dist = "SOM neighbor distances") { #title of neighbor distances plot
  
  # Extract cluster numbers 
  if (is.null(SOM.output$som_models) || is.null(SOM.output$som_clusters)) {
    stop("SOM.output is missing required elements: 'som_models' or 'som_clusters'")
  } else {
    som_model <- SOM.output$som_models[[1]]
    som_cluster <- SOM.output$som_clusters[[1]]
  }
  
  
  # Save file
  if (save) {
    
    # Set default file name
    if (is.null(file.name)) {
      file.name <- paste0("SOM_model_plot_", paste(SOM.output$input_data_names, collapse = "_"), ".", type)
    }
    
    # Set overwrite option
    if (file.exists(file.name) && !overwrite) {
      stop(paste("File", file.name, "already exists - set overwrite = T to overwrite"))
    }
    
    # Set plot format
    if (type == "svg") {
      svg(file.name, 
          width = width / 2.54, 
          height = height / 2.54)
    } else if (type == "png") {
      png(file.name, 
          width = width, 
          height = height, 
          res = resolution, 
          units = "cm")
    } else if (type == "jpg") {
      jpeg(file.name, 
           width = width,
           height = height, 
           res = resolution, 
           units = "cm")
    } else {
      stop("Unsupported file type for plotting - choose from 'svg', 'png' or 'jpg'")
    }
  }
  
  # Set plotting area
  par(mfrow = c(2, 1), 
      mar = c(margin.bottom, 
              margin.left, 
              margin.top, 
              margin.right))
  
  # Plot SOM neighbor distances (top plot)
  plot(x = som_model,
       type = "dist.neighbours",
       main = title.neighbor.dist,
       shape = cluster.shape.neighbor.dist,
       palette.name = function(n) rev(col.pal.neighbor.dist(n)))
  
  # Set color palette for bottom plot
  k.cols <- col.pal.clusters(max(som_cluster))
  SOM_cluster_plot_col <- rep(NA, length(som_cluster))
  for (i in seq_len(max(som_cluster))) {
    SOM_cluster_plot_col[som_cluster == i] <- k.cols[i]
  }
  
  # Ensure valid shift.plot.clusters value (for bottom plot)
  if (shift.plot.clusters >= 0.5 | shift.plot.clusters < 0) {
    message("Invalid shift.plot.clusters value (needs to be 0-0.5) - default of 0 is used")
    shift.plot.clusters <- 0
  }
  
  # Adjust position of bottom plot to align with grid of top plot
  par(fig = c(shift.plot.clusters, 1, 0, 0.5), new = T) #shift bottom plot
  
  # Plot SOM clusters (bottom plot)
  plot(x = som_model,
       shape = cluster.shape.clusters,
       type = "mapping",
       bgcol = SOM_cluster_plot_col,
       main = title.clusters,
       pch = point.shape.clusters,
       cex = point.size.clusters,
       col = point.col.clusters)
  
  # Add boundaries of SOM clusters
  if (max(som_cluster) > 1) {
    add.cluster.boundaries(som_model,
                           som_cluster,
                           lwd = boundary.lwd.clusters,
                           col = boundary.col.clusters)
  }
  
  # Close graphics device
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}

## Function to plot Structure-like barplots and save ancestry proportions as csv file
plot.Structure.SOM <- function(SOM.output, 
                               col.pal = viridis::viridis, #color palette
                               save = F, #save plot
                               overwrite = T, #overwrite plot if already present
                               type = "svg", #plot file type (choose: png, jpg or svg)
                               file.name = NULL, #plot file name (if NULL, default name is created)
                               width = 10, #plot width in cm
                               height = 15, #plot height in cm
                               resolution = 300, #plot resolution in dpi
                               margin.top = 2,
                               margin.bottom = 5,
                               margin.left = 4, 
                               margin.right = 2,
                               sort.by.col = 1, #specify integer giving column index of ancestry matrix for ordering rows of ancestry matrix (if NULL, hierarchical ordering is performed)
                               linkage.method = "single") { #agglomeration method used for hierarchical clustering (see hclust function)
  
  # Ensure ancestry_matrix is valid
  if (is.null(SOM.output$ancestry_matrix) || !is.matrix(SOM.output$ancestry_matrix)) {
    stop("ancestry_matrix of SOM.output not valid")
  }
  
  # If there's only one cluster, skip plotting
  if (ncol(SOM.output$ancestry_matrix) == 1) { 
    stop("Only one cluster detected - skipping structure plot")
  }
  
  # Order rows of ancestry_matrix
  if (!is.null(sort.by.col)) {
    if (sort.by.col > ncol(SOM.output$ancestry_matrix)) {
      stop(paste0("sort.by.col exceeds number of columns in ancestry_matrix - select number from 1-", ncol(SOM.output$ancestry_matrix), " or NULL to prevent sorting"))
    }
    ancestry_proportions <- SOM.output$ancestry_matrix[order(SOM.output$ancestry_matrix[, sort.by.col]), ]
  } else {
    ancestry_proportions <- SOM.output$ancestry_matrix #no sorting when set as NULL
  }
  
  # Perform hierarchical clustering on distance matrix
  cluster_order <- hclust(dist(ancestry_proportions), 
                          method = linkage.method)$order
  SOM_ancestry_proportions <- ancestry_proportions[cluster_order, ]
  
  # Generate layer colors
  layer_colors <- col.pal(ncol(SOM_ancestry_proportions)) 
  
  # Save file
  if (save) {
    
    # Assign file.name if it is NULL
    if (is.null(file.name)) {
      file.name <- paste0("SOM_structure_plot_", paste(SOM.output$input_data_names, collapse = "_"), ".", type)
    }
    
    # Set overwrite option
    if (file.exists(file.name) && !overwrite) {
      stop(paste(file.name, "already exists - set overwrite = TRUE to overwrite"))
    }
    
    # Set plot format type
    if (type == "svg") {
      svg(file.name, 
          width = width / 2.54, 
          height = height / 2.54)
    } else if (type == "png") {
      png(file.name, 
          width = width, 
          height = height, 
          res = resolution, 
          units = "cm")
    } else if (type == "jpg") {
      jpeg(file.name, 
           width = width, 
           height = height, 
           res = resolution, 
           units = "cm")
    } else {
      stop("Unsupported file type - choose from 'svg', 'png', or 'jpg'")
    }
  }
  
  # Set plot layout
  par(mfrow = c(1, 1))
  
  # Generate structure plot
  make.structure.plot(admix.proportions = SOM_ancestry_proportions, 
                      sample.names = rownames(SOM_ancestry_proportions), 
                      mar = c(margin.bottom, 
                              margin.left, 
                              margin.top, 
                              margin.right), 
                      layer.colors = layer_colors,
                      sort.by = sort.by.col)
  
  # Close graphics device
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}

## Function to plot sample map with cluster assignment for each individual
plot.Map.SOM <- function(SOM.output,
                         Coordinates, #coordinates as "Latitude" and "Longitude" columns in dataframe/matrix
                         lat.buffer.range = 2, #add coordinates as buffer range around latitude coordinates
                         lon.buffer.range = 2, #add coordinates as buffer range around longitude coordinates
                         save = F, #whether to save plot or not
                         overwrite = T, #whether to overwrite if file exists
                         file.name = NULL, #plot file name (NULL = default file name)
                         type = "svg", #plot format type (options: "png", "svg", "jpg")
                         width = 15, #plot width in cm
                         height = 20, #plot height in cm
                         resolution = 300, #plot resolution in dpi
                         pie.size = 1.9, #pie chart size
                         pie.col.pal = viridis::viridis, #color palette of pie charts
                         USA.add.states = T, #option to add US states (only works if range includes USA)
                         USA.add.counties = F, #option to add US counties (only works if range includes USA)
                         USA.state.lwd = 0.5, #linewidth of US state borders (only works if range includes USA)
                         USA.county.lwd = 0.3, #linewidth of US county borders (only works if range includes USA)
                         north.arrow.position = c(0.03, 0.88), #position (x, y) of north arrow relative to map
                         north.arrow.length = 0.7, #length of north arrow
                         north.arrow.lwd = 2, #linewidth of north arrow
                         north.arrow.N.position = 0.3, #position of north arrow "N"
                         north.arrow.N.size = 1, #size of north arrow "N"
                         scale.position = c(0.03, 0.05), #relative position (x, y) of scale
                         scale.size = 0.16, #size of scale
                         scale.font.size = 0.54, #font size of scale text
                         legend.position = "topright", #position of legend
                         legend.cluster.names = NULL, #names of clusters in legend (if NULL, default is used, otherwise use vector with length of number of clusters)
                         legend.font.size = 1, #font size of legend text
                         legend.box = T, #create white box around legend
                         legend.symbol.size = 1.5) { #size of legend symbols
  
  # Ensure ancestry_matrix is valid
  if (is.null(SOM.output$ancestry_matrix) || !is.matrix(SOM.output$ancestry_matrix)) {
    stop("ancestry_matrix of SOM.output is not valid")
  }
  
  # Check if Coordinates is data frame or matrix
  if (!is.data.frame(Coordinates) && !is.matrix(Coordinates)) {
    stop("Coordinates must be a data frame or matrix")
  }
  
  # Convert matrix to data frame if necessary
  if (is.matrix(Coordinates)) {
    Coordinates <- as.data.frame(Coordinates)
  }
  
  # Check if Coordinates contain Latitude and Longitude
  if (!all(c("Latitude", "Longitude") %in% names(Coordinates))) {
    stop("Coordinates must contain 'Latitude' and 'Longitude' columns")
  }
  
  # Check for NA values in Latitude and Longitude
  if (any(is.na(Coordinates$Latitude)) || any(is.na(Coordinates$Longitude))) {
    stop("Coordinates must not contain NA values in 'Latitude' or 'Longitude' columns")
  }
  
  # Prepare ancestry matrix
  q_matrix = as.data.frame(SOM.output$ancestry_matrix) #convert ancestry_matrix to dataframe
  q_matrix$Ind = rownames(SOM.output$ancestry_matrix) #add sample names
  
  # Check if number of rows in ancestry matrix matches number of coordinates
  if (nrow(q_matrix) != nrow(Coordinates)) {
    stop("Number of samples in ancestry_matrix does not match number of samples in Coordinates")
  }
  
  # Define color palette for pie charts
  k.cols = pie.col.pal(ncol(SOM.output$ancestry_matrix))
  
  # Define map boundaries
  lon_min = min(Coordinates$Longitude) - lon.buffer.range
  lon_max = max(Coordinates$Longitude) + lon.buffer.range
  lat_min = min(Coordinates$Latitude) - lat.buffer.range
  lat_max = max(Coordinates$Latitude) + lat.buffer.range
  
  # Create plot
  if (dev.cur() > 1) {
    dev.off() #close current device if open to avoid unwanted graphic distortions and other effects
  }
  
  # Set plot saving
  if (save) {
    
    # Set file name for saving
    if (is.null(file.name)) {
      file.name <- paste0("SOM_map_plot_", paste(SOM.output$input_data_names, collapse = "_"), ".", type)
    }
    
    # Set overwriting 
    if (file.exists(file.name) && !overwrite) {
      stop(paste(file.name, "already exists - set overwrite = T to overwrite"))
    }
    
    # Set file type
    if (type == "svg") {
      svg(file.name, 
          width = width / 2.54, 
          height = height / 2.54)
    } else if (type == "png") {
      png(file.name, 
          width = width, 
          height = height, 
          res = resolution, 
          units = "cm")
    } else if (type == "jpg") {
      jpeg(file.name, 
           width = width, 
           height = height, 
           res = resolution, 
           units = "cm")
    } else {
      stop("Unsupported plot file type - choose from 'svg', 'png', or 'jpg'")
    }
  }
  
  # Set layout and margins
  par(mfrow = c(1, 1),
      mar = c(1, 1, 1, 1))
  
  # Create map
  map("world", 
      fill = T, 
      col = "lightgrey", 
      xlim = c(lon_min, lon_max), 
      ylim = c(lat_min, lat_max))
  map.axes()
  
  # Add US counties if requested
  if (USA.add.counties) {
    map("county", 
        add = T, 
        col = "grey", 
        lwd = USA.county.lwd)
  }
  
  # Add US states if requested
  if (USA.add.states) {
    map("state",
        add = T, 
        col = "black", 
        lwd = USA.state.lwd)
  }
  
  # Add pie charts to map
  for (i in 1:nrow(q_matrix)) {
    coords = matrix(c(Coordinates$Longitude[i], Coordinates$Latitude[i]), ncol = 2, byrow = T)
    make.admix.pie.plot(admix.proportions = as.matrix(q_matrix[i, -ncol(q_matrix), drop = F]),
                        coords = coords,
                        layer.colors = k.cols,
                        radii = pie.size,
                        add = T)
  }
  
  # Define legend labels
  if (is.null(legend.cluster.names)) {
    legend.labels <- paste("Cluster", 1:length(k.cols)) #set default labels
  } else {
    if (length(legend.cluster.names) != length(k.cols)) {
      print("Length of legend.cluster.names must match number of clusters - default labels are used")
      legend.labels <- paste("Cluster", 1:length(k.cols)) #default labels
    }
    legend.labels <- legend.cluster.names #use custom labels
  }
  
  # Set legend box
  if (legend.box) {
    legend_box <- "o"
  } else {
    legend_box <- "n"  
  }
  
  # Add legend
  legend(legend.position, 
         legend = legend.labels,
         pch = 21,
         cex = legend.font.size,
         pt.cex = legend.symbol.size,
         pt.bg = k.cols,
         bty = legend_box)
  
  # Add scale
  scale_position_x = scale.position[1] * (lon_max - lon_min) + lon_min
  scale_position_y = scale.position[2] * (lat_max - lat_min) + lat_min
  map.scale(x = scale_position_x,
            y = scale_position_y,
            cex = scale.font.size,
            relwidth = scale.size,
            ratio = F)
  
  # Add north arrow
  north_arrow_x = north.arrow.position[1] * (lon_max - lon_min) + lon_min
  north_arrow_y = north.arrow.position[2] * (lat_max - lat_min) + lat_min
  arrows(x0 = north_arrow_x, 
         y0 = north_arrow_y, 
         x1 = north_arrow_x, 
         y1 = north_arrow_y + north.arrow.length, 
         length = 0.13, 
         col = "black", 
         lwd = north.arrow.lwd)
  
  # Add North "N" above north arrow
  text(x = north_arrow_x, 
       y = north_arrow_y + north.arrow.length + north.arrow.N.position, #adjust position for "N"
       labels = "N", 
       cex = north.arrow.N.size, 
       col = "black")
  
  # Close graphics device
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}

## Function to plot variable importance for each SOM layer (based on Codebook Vectors/Neuron Weights)
plot.Variable.Importance.SOM <- function(SOM.output, 
                                         col.pal = viridis::turbo, #color palette
                                         bars.threshold.N = 50, #threshold for leaving out bar labels
                                         save = F, #option to save plot
                                         overwrite = T, #option to overwrite plot if it alredy exists
                                         type = "svg", #plot file type (options: "png", "svg", "jpg")
                                         width = 20, #plot width in cm
                                         height = 15, #plot height in cm
                                         resolution = 300, #plot resolution in dpi
                                         file.name = NULL, #plot file name (if NULL, default file name is used)
                                         bottom.margin = 3,
                                         left.margin = 5,
                                         top.margin = 3,
                                         right.margin = 2,
                                         title.font.size = 1.2, #font size of title
                                         matrix.label.font.size = 1, #font size of matrix label
                                         bar.label.font.size = 0.65, #font size of bar labels
                                         importance.threshold = 0.001) { #threshold for showing variable importance
  
  # Extract SOM codes from SOM.output
  if (is.null(SOM.output$som_models)) {
    stop("som_models not found in SOM.output")
  } else {
    codes <- SOM.output$som_models[[1]]$codes
    m.codes <- lapply(codes, function(x) apply(x, 2, median)) #calculate median for each column
  }
  
  # Extract matrix names from SOM.output
  if (is.null(SOM.output$input_data_names)) {
    stop("input_data_names not found in SOM.output")
  } else {
    matrix_names <- SOM.output$input_data_names
  }
  
  # Set plot saving
  if (save) {
    
    # Set default plotting name
    if (is.null(file.name)) {
      file.name <- paste0("SOM_variable_importance_plot_", paste(matrix_names, collapse = "_"), ".", type)
    }
    
    # Check overwrite
    if (file.exists(file.name) && !overwrite) {
      stop(paste(file.name, "already exists - set overwrite = T to overwrite"))
    }
    
    # Set plot type
    if (type == "svg") {
      svg(file.name, 
          width = width / 2.54, 
          height = height / 2.54)
    } else if (type == "png") {
      png(file.name,
          width = width, 
          height = height, 
          res = resolution, 
          units = "cm")
    } else if (type == "jpg") {
      jpeg(file.name, 
           width = width, 
           height = height, 
           res = resolution, 
           units = "cm")
    } else {
      stop("Unsupported file type - choose from 'svg', 'png', or 'jpg'")
    }
  }
  
  # Set plot margins
  par(mar = c(bottom.margin, 
              left.margin, 
              top.margin, 
              right.margin))
  
  # Set plot layout based on number of SOM layers
  if (length(m.codes) == 1) {
    par(mfrow = c(1, 1))
  } else if (length(m.codes) == 2) {
    par(mfrow = c(1, 2))
  } else if (length(m.codes) == 3) {
    par(mfrow = c(1, 3))
  } else if (length(m.codes) == 4) {
    par(mfrow = c(2, 2))
  } else if (length(m.codes) == 5) {
    par(mfrow = c(2, 3))
  } else if (length(m.codes) == 6) {
    par(mfrow = c(2, 3))
  } else if (length(m.codes) == 7) {
    par(mfrow = c(2, 4))
  } else if (length(m.codes) == 8) {
    par(mfrow = c(2, 4))
  } else if (length(m.codes) == 9) {
    par(mfrow = c(3, 3))
  } else {
    num_rows <- ceiling(length(m.codes) / 3)
    par(mfrow = c(num_rows, 3))
  }
  
  # Iterate over each matrix and generate plots
  for (i in seq_along(m.codes)) {
    
    # Extract sorted data and corresponding labels
    sorted_data <- sort(m.codes[[i]][which(m.codes[[i]] > importance.threshold)])
    y_labels <- names(sorted_data)
    num_bars <- length(sorted_data)
    if (num_bars == 0) {
      warning(paste("No variables exceed importance.threshold of", importance.threshold, "for matrix", matrix_names[i], "- specify lower value for importance.threshold"))
      next #skip plot if no variables meet threshold
    }
    if (importance.threshold != 0.001) {
      message(paste("Only showing variables above importance threshold of", importance.threshold))
    }
    
    # Dynamically suppress labels when too many bars are present
    if (num_bars > bars.threshold.N) {
      y_labels <- rep("", num_bars) #suppress y-axis labels
      message(paste("Matrix", matrix_names[[i]], "has too many bars so that Y-axis labels are suppressed - increase bars.threshold.N to show labels"))
    }
    
    # Define color palette
    layer.cols <- col.pal(length(codes)) 
    
    # Create plot
    barplot(sorted_data, 
            horiz = T, 
            las = 1, 
            col = layer.cols[i], 
            xlim = c(0, 1), 
            names.arg = y_labels, 
            cex.names = bar.label.font.size)
    
    # Add input data name in bottom right corner of plot
    text(x = 0.95, 
         y = 0, 
         labels = matrix_names[i], 
         adj = c(1, 0), 
         cex = matrix.label.font.size,
         col = "black")
  }
  
  # Reset layout for title
  layout(matrix(c(1), 
                nrow = 1, 
                ncol = 1))
  par(mar = c(0, 0, 3, 0))
  
  # Add main title
  if (length(m.codes) > 1) {
    title(main = "Variable importance of SOM layers", 
          cex.main = title.font.size)
  } else {
    title(main = "Variable importance", 
          cex.main = title.font.size
    )
  }
  
  # Close graphics device
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
  
  SOM_results <- list(cluster_assignment = SOM.output$cluster_assignment, 
                      distance_matrix = SOM.output$distance_matrix, 
                      learning_values_list = SOM.output$learning_values_list, 
                      BIC_values = SOM.output$BIC_values,
                      ancestry_matrix = SOM.output$ancestry_matrix,
                      optim_k = SOM.output$optim_k,
                      optim_k_mean = SOM.output$optim_k_mean, 
                      input_data_names = as.character(SOM.output$input_data_names),
                      N_steps = SOM.output$N.steps,
                      N_replicates = SOM.output$N.replicates,
                      max_k = SOM.output$max.k,
                      learning_rate_initial = SOM.output$learning.rate.initial,
                      learning_rate_final = SOM.output$learning.rate.final,
                      som_models = SOM.output$som_models,
                      som_clusters = SOM.output$som_clusters,
                      cluster_gridcell_assignments = SOM.output$cluster_gridcell_assignments,
                      variable_importance = m.codes)
  
  return(SOM_results)
} 

## Function to normalize matrices by column, if needed
normalize_columnwise <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    col <- mat[, j]
    if (any(col < 0, na.rm = TRUE) || any(col > 1, na.rm = TRUE)) {
      min_val <- min(col, na.rm = TRUE)
      max_val <- max(col, na.rm = TRUE)
      if (max_val != min_val) {
        mat[, j] <- (col - min_val) / (max_val - min_val)
      } else {
        mat[, j] <- 0  # or NA, depending on context
      }
    }
  }
  mat
}
