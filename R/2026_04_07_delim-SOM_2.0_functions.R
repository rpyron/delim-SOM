#### Install required R packages

CRAN_packages <- c(
  "adegenet",    #genetic data manipulation
  "ape",         #phylogenetics
  "caroline",    #pie charts
  "clue",        #Hungarian algorithm / solve_LSAP
  "cluster",     #silhouette
  "clusterCrit", #Davies-Bouldin index
  "dbscan",      #HDBSCAN clustering
  "doRNG",       #reproducible RNG
  "doParallel",  #parallel backend
  "FactoMineR",  #multiple factor analysis
  "foreach",     #parallel processing
  "kohonen",     #SOM / supersom
  "maps",        #mapping
  "mclust",      #Gaussian mixture models
  "poppr",       #population genetics
  "vcfR",        #VCF file handling
  "viridis"      #color palettes
)


## Install missing CRAN packages
for (pkg in CRAN_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}



#### Functions

## Function to train single-layer SOM (one matrix) or multi-layer Super-SOM (multiple matrices)
train.SOM <- function(input_data, #one matrix/dataframe or multiple matrices/dataframes provided as list()
                      N.steps = 200, #number of training iterations S for SOM
                      N.replicates = 110, #number of SOM runs R
                      parallel = TRUE, #whether to run SOM training in parallel 
                      N.cores = 3, #number of cores for training SOM in parallel (if parallel = TRUE)
                      grid.size = NULL, #grid size - specify as c(x, y) if desired (if grid.size = NULL, size and shape is automatically determined)
                      grid.multiplier = 5, #tuning parameter defining how “fine” or “coarse” SOM grid is relative to number of samples (recommended: 5)
                      learning.rate.initial = 0.5, #initial learning rate for SOM training
                      learning.rate.final = 0.1, #final learning rate for SOM training
                      learning.rate.tuning = FALSE, #whether to perform learning.rate tuning to choose best initial and final learning rate values
                      layer.distance.functions = NULL, #distance function(s) for each data layer (if NULL, inferred automatically: continuous = sumofsquares, binary = tanimoto, count = manhattan)
                      manual.layer.weights = NULL, #whether to specify manual weights for multiple layers (NULL = equal weights)
                      max.NA.row = 0.3, #maximum fraction of missing values allowed per row (sample) in input data to prevent row to be removed
                      max.NA.col = 0.6, #maximum fraction of missing values allowed per column (variable) in input data to prevent row to be removed
                      training.neighborhoods = "gaussian", #neighborhood function used for SOM training (options: "gaussian" or "bubble")
                      save.SOM.results = FALSE, #whether to save SOM results to file
                      save.SOM.results.name = NULL, #file name for saving SOM results (if NULL, default name based on input_data is generated; if save.SOM.results = TRUE)
                      overwrite.SOM.results = FALSE, #if FALSE, existing results are loaded instead of re-running SOM
                      verbose = TRUE, #whether to show messages
                      message.N.replicates = 20, #frequency of progress messages during training (message is printed every message.N.replicates iterations)
                      set.seed.N = 1 #set seed for reproducibility
) {
  
  # Set messages
  messager <- function(...) if (isTRUE(verbose)) message(...)
  
  # Start processing input data
  messager("PROCESSING INPUT DATA ...")
  
  # Validate specified input_data
  if (is.null(input_data)) {
    stop("Data processing aborted: input_data cannot be NULL")
  }
  if (!(is.matrix(input_data) || is.data.frame(input_data) || is.list(input_data))) {
    stop("Data processing aborted: input_data must be a matrix, data.frame or list of such objects")
  }
  if (is.list(input_data) && length(input_data) == 0) {
    stop("Data processing aborted: input_data is an empty list")
  }
  
  # Validate specified N.steps
  if (!is.numeric(N.steps) || length(N.steps) != 1 || is.na(N.steps) || N.steps < 1 || (N.steps %% 1 != 0)) {
    stop("Data processing aborted: N.steps must be a single positive integer (>= 1)")
  }
  if (N.steps < 60) {
    messager("Warning: N.steps is low (", N.steps, ") - SOM training may be unstable (recommended: 60–200)")
  }
  if (N.steps > 200) {
    messager("Warning: N.steps is high (", N.steps, ") - computation will be slow (recommended: 60–200)")
  }
  
  # Validate specified N.replicates
  if (!is.numeric(N.replicates) || length(N.replicates) != 1 || is.na(N.replicates) || N.replicates < 1 || (N.replicates %% 1 != 0)) {
    stop("Data processing aborted: N.replicates must be a single positive integer (>= 1)")
  }
  if (N.replicates < 30) {
    messager("Warning: N.replicates is low (", N.replicates, ") - results may be unreliable (recommended: 30–100)")
  }
  if (N.replicates > 200) {
    messager("Warning: N.replicates is high (", N.replicates, ") - computation will be slow (recommended: 50–150)")
  }
  
  # Validate specified parallel
  if (!is.logical(parallel) || length(parallel) != 1 || is.na(parallel)) {
    stop("Data processing aborted: parallel must be TRUE or FALSE")
  }
  
  # Validate specified N.cores
  if (parallel) {
    if (!is.numeric(N.cores) || length(N.cores) != 1 || is.na(N.cores) || N.cores < 1 || (N.cores %% 1 != 0)) {
      stop("Data processing aborted: N.cores must be a single positive integer (>= 1)")
    }
    max_cores <- parallel::detectCores(logical = FALSE)
    if (is.na(max_cores) || max_cores < 1) max_cores <- parallel::detectCores(logical = TRUE)
    if (is.na(max_cores) || max_cores < 1) max_cores <- 1
    if (N.cores > max_cores) {
      messager(sprintf("Requested N.cores (%d) exceeds available cores (%d) - using %d cores", N.cores, max_cores, max_cores))
      N.cores <- max_cores
    }
  }
  
  # Validate specified grid.size
  if (!is.null(grid.size)) {
    if (!is.numeric(grid.size) || length(grid.size) != 2 || any(is.na(grid.size)) || any(grid.size <= 0) || any(grid.size %% 1 != 0)) {
      stop("Input aborted: grid.size must be NULL or numeric vector of length 2 with positive integers (e.g., c(5, 5))")
    }
  }
  
  # Validate specified grid.multiplier
  if (!is.numeric(grid.multiplier) || length(grid.multiplier) != 1 || grid.multiplier < 1) {
    stop("Data processing aborted: grid.multiplier must be a single numeric value (recommended: 5)")
  }
  
  # Validate specified learning.rate.tuning
  if (!is.logical(learning.rate.tuning) || length(learning.rate.tuning) != 1 || is.na(learning.rate.tuning)) {
    stop("Data processing aborted: learning.rate.tuning must be TRUE or FALSE")
  }
  
  # Validate specified learning rate parameters if learning rate tuning is not done
  if (!learning.rate.tuning) {
    if (!is.numeric(learning.rate.initial) || length(learning.rate.initial) != 1 || learning.rate.initial <= 0 || learning.rate.initial > 1) {
      stop("Data processing aborted: learning.rate.initial must be a single numeric value between 0 and 1 (e.g. 0.6)")
    }
    if (!is.numeric(learning.rate.final) || length(learning.rate.final) != 1 || is.na(learning.rate.final) || learning.rate.final < 0 || learning.rate.final > 1) {
      stop("Data processing aborted: learning.rate.final must be a single numeric value between 0 and 1 (e.g. 0.1)")
    }
    if (learning.rate.final > learning.rate.initial) {
      stop("Data processing aborted: learning.rate.final must be smaller than learning.rate.initial")
    }
  }
  
  # Validate specified layer.distance.functions
  if (!is.null(layer.distance.functions)) {
    n_layers_expected <- 1
    if (is.list(input_data) && length(input_data) > 1 && !is.data.frame(input_data)) n_layers_expected <- length(input_data)
    if (!(is.character(layer.distance.functions) || is.list(layer.distance.functions))) {
      stop("Data processing aborted: layer.distance.functions must be NULL, a character vector, or a list")
    }
    if (is.character(layer.distance.functions)) {
      if (any(is.na(layer.distance.functions)) || any(trimws(layer.distance.functions) == "")) {
        stop("Data processing aborted: layer.distance.functions contains NA or empty strings")
      }
      if (!(length(layer.distance.functions) %in% c(1, n_layers_expected))) {
        stop(sprintf("Data processing aborted: layer.distance.functions must have length 1 or %d (number of layers in input_data)", n_layers_expected))
      }
    }
    if (is.list(layer.distance.functions)) {
      if (!(length(layer.distance.functions) %in% c(1, n_layers_expected))) {
        stop(sprintf("Data processing aborted: layer.distance.functions must have length 1 or %d (number of layers in input_data)", n_layers_expected))
      }
      bad_types <- vapply(layer.distance.functions, function(x) !(is.character(x) || is.function(x)), logical(1))
      if (any(bad_types)) stop("Data processing aborted: each element of layer.distance.functions must be a character string or a function")
      bad_chr <- vapply(layer.distance.functions, function(x) is.character(x) && (length(x) != 1 || is.na(x) || trimws(x) == ""), logical(1))
      if (any(bad_chr)) stop("Data processing aborted: layer.distance.functions list contains invalid character entries (must be single non-empty strings)")
    }
  }
  
  # Validate specified manual.layer.weights
  if (!is.null(manual.layer.weights)) {
    if (!is.numeric(manual.layer.weights) || any(is.na(manual.layer.weights)) || any(manual.layer.weights <= 0)) {
      stop("Data processing aborted: manual.layer.weights must be NULL or a numeric vector of positive values")
    }
  }
  
  # Validate specified NA‐max.NA.row
  if (!is.numeric(max.NA.row) || length(max.NA.row) != 1 || max.NA.row < 0 || max.NA.row > 1) {
    stop("Data processing aborted: max.NA.row must be a single numeric value between 0 and 1 (recommended: 0.3)")
  }
  
  # Validate specified NA‐max.NA.col
  if (!is.numeric(max.NA.col) || length(max.NA.col) != 1 ||
      max.NA.col < 0 || max.NA.col > 1) {
    stop("Data processing aborted: max.NA.col must be a single numeric value between 0 and 1 (recommended: 0.7)")
  }
  
  # Validate specified training.neighborhoods
  if (!is.character(training.neighborhoods) || length(training.neighborhoods) != 1 || !(training.neighborhoods %in% c("gaussian", "bubble"))) {
    stop("Data processing aborted: training.neighborhoods must be 'gaussian' or 'bubble'")
  }
  
  # Validate specified save.SOM.results
  if (!is.logical(save.SOM.results) || length(save.SOM.results) != 1 || is.na(save.SOM.results)) {
    stop("Data processing aborted: save.SOM.results must be TRUE or FALSE")
  }
  
  # Validate specified save.SOM.results.name
  if (save.SOM.results && !is.null(save.SOM.results.name)) {
    if (!is.character(save.SOM.results.name) || length(save.SOM.results.name) != 1 || is.na(save.SOM.results.name) || trimws(save.SOM.results.name) == "") {
      stop("Data processing aborted: save.SOM.results.name must be non-empty character string (file path) if provided")
    }
    valid_ext <- tolower(tools::file_ext(save.SOM.results.name)) #extract extension
    if (valid_ext != "rdata") {
      stop("Data processing aborted: save.SOM.results.name must end with '.Rdata'") #abort if not .Rdata
    }
  }
  
  # Validate specified overwrite.SOM.results
  if (!is.logical(overwrite.SOM.results) || length(overwrite.SOM.results) != 1 || is.na(overwrite.SOM.results)) {
    stop("Data processing aborted: overwrite.SOM.results must be TRUE or FALSE")
  }
  
  # Validate specified verbose
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("Data processing aborted: verbose must be TRUE or FALSE")
  }
  
  # Validate message.N.replicates
  if (!is.numeric(message.N.replicates) || length(message.N.replicates) != 1 || message.N.replicates < 1 || (message.N.replicates %% 1 != 0)) {
    stop("Data processing aborted: message.N.replicates must be a single positive integer (>= 1)")
  }
  
  # Validate set.seed.N
  if (!is.numeric(set.seed.N) || length(set.seed.N) != 1 || set.seed.N < 1 || (set.seed.N %% 1 != 0)) {
    stop("Data processing aborted: set.seed.N must be a single positive integer (>= 1)")
  }
  
  
  # Extract input_data_names and set save.SOM.results.name for saving ...
  
  # ... for list with multiple data sets
  if (is.list(input_data) && length(input_data) > 1 && !is.data.frame(input_data)) {
    if (!is.null(names(input_data)) && all(names(input_data) != "")) { #use names in list if present
      input_data_names <- names(input_data)
    } else { #fallback: deparse list(...) call
      list_names <- match.call()$input_data
      input_data_names <- sapply(as.list(list_names)[-1], deparse)
    }
    if (is.null(save.SOM.results.name)) { #assign default saving name if save.SOM.results.name is NULL
      save.SOM.results.name <- paste0("SOM_results_", paste(input_data_names, collapse = "_"), ".Rdata")
    }
    
    # ... for list with one data set
  } else if (is.list(input_data) && length(input_data) == 1) {
    if (!is.null(names(input_data)) && names(input_data)[1] != "") {
      input_data_names <- names(input_data)[1]
    } else {
      input_data_names <- deparse(substitute(input_data)) #extract names of dataset
      input_data_names <- gsub("^list\\((.+)\\)$", "\\1", input_data_names) #remove "list" from name
      input_data_names <- gsub("\"", "", input_data_names) #remove quotes from name
    }
    if (is.null(save.SOM.results.name)) { #assign default name for saving if save.SOM.results.name is NULL
      save.SOM.results.name <- paste0("SOM_results_", input_data_names, ".Rdata")
    }
    
    # ... for non-list object with one dataset
  } else {
    input_data_names <- deparse(substitute(input_data)) #extract names of dataset
    if (is.null(save.SOM.results.name)) { #assign default name for savings if save.SOM.results.name is NULL
      save.SOM.results.name <- paste0("SOM_results_", input_data_names, ".Rdata")
    }
  }
  
  # If overwrite.SOM.results is FALSE and file already exists, return saved results
  if (!overwrite.SOM.results && file.exists(save.SOM.results.name)) {
    messager("SOM results already exist - loading results from file and skipping SOM run")
    load(save.SOM.results.name)
    
    required_fields <- c("distance_weights_matrix", "learning_values_list")
    if (!all(required_fields %in% names(SOM_results))) { #check if structure of SOM_results is correct
      stop("Data processing aborted: could not load SOM results (results do not contain expected objects 'distance_weights_matrix' and 'learning_values_list') - check saved file or rerun SOM")
    }
    return(SOM_results)
  }
  
  
  # Check and transform data if necessary ...
  
  # ... for multiple datasets
  if (is.list(input_data) && length(input_data) > 1 && !is.data.frame(input_data)) { 
    
    # Ensure each layer has rownames
    if (any(vapply(input_data, function(x) is.null(rownames(x)), logical(1)))) {
      stop("Data processing aborted: all provided data layers must have matching rownames")
    }
    
    # Extract shared samples across all matrices (before filtering)
    all_samples <- unique(unlist(lapply(input_data, rownames)))
    common_samples <- Reduce(intersect, lapply(input_data, rownames))
    not_shared <- setdiff(all_samples, common_samples)
    if (length(not_shared) > 0) { #show message if samples are removed
      n_all <- length(all_samples)
      n_removed <- length(not_shared)
      if (n_removed <= 30) {
        not_shared_rows <- paste(not_shared, collapse = ", ")
        messager(sprintf(
          "Removed %d of %d rows (samples) due to non-matching rownames: %s",
          n_removed, n_all, not_shared_rows
        ))
      } else {
        messager(sprintf(
          "Removed %d of %d rows (samples) due to non-matching rownames",
          n_removed, n_all
        ))
      }
    }
    if (length(common_samples) == 0) { #stop with message if no shared samples remain
      stop("Data processing aborted: no matching rownames (samples) across layers - check input data")
    }
    if (length(common_samples) == 1) { #stop with message if only one shared samples remain
      stop("Data processing aborted: only one row (sample) matches rownames across layers - check input data")
    }
    input_data <- lapply(input_data, function(mat) mat[common_samples, , drop = FALSE])
    
    # Filter by max.NA.row
    for (i in seq_along(input_data)) {
      mat <- input_data[[i]]
      na_frac <- rowMeans(is.na(mat))
      removed <- rownames(mat)[na_frac > max.NA.row]
      if (length(removed) > 0) { #show message if samples are removed
        if (length(removed) <= 30) {
          messager(sprintf(
            "Removed %d of %d rows (samples) in dataset %s due to more than %.0f%% NA in data (max.NA.row = %.2f): %s",
            length(removed), nrow(mat), input_data_names[i], max.NA.row * 100, max.NA.row, paste(removed, collapse = ", ")
          ))
        } else {
          messager(sprintf(
            "Removed %d of %d rows (samples) in dataset %s due to more than %.0f%% NA in data (max.NA.row = %.2f)",
            length(removed), nrow(mat), input_data_names[i], max.NA.row * 100, max.NA.row
          ))
        }
      }
      mat <- mat[na_frac <= max.NA.row, , drop = FALSE]
      if (nrow(mat) == 0) {
        stop(sprintf(
          "Data processing aborted: no rows (samples) remain in dataset %s after applying max.NA.row = %.2f - check input data or increase max.NA.row",
          input_data_names[i], max.NA.row
        ))
      }
      if (nrow(mat) == 1) {
        stop(sprintf(
          "Data processing aborted: only one row (sample) remains in dataset %s after applying max.NA.row = %.2f - check input data or increase max.NA.row",
          input_data_names[i], max.NA.row
        ))
      }
      input_data[[i]] <- mat
    }
    
    # Intersect again to obtain final shared samples
    all_samples <- unique(unlist(lapply(input_data, rownames)))
    common_samples <- Reduce(intersect, lapply(input_data, rownames))
    not_shared <- setdiff(all_samples, common_samples)
    if (length(common_samples) == 0) { #stop with message if no shared samples remain
      stop("Data processing aborted: no shared samples remain after NA filtering - check input data or increase max.NA.row")
    }
    if (length(common_samples) == 1) { #stop with message if only one shared sample remain
      stop("Data processing aborted: only one shared sample remains after NA filtering - check input data or increase max.NA.row")
    }
    input_data <- lapply(input_data, function(mat) mat[common_samples, , drop = FALSE])
    dataset_names <- input_data_names
    processed_data <- list()
    processed_names <- character(0)
    
    # For each dataset, ...
    for (i in seq_along(input_data)) {
      
      # Extract data and data names
      name <- dataset_names[i]
      mat <- as.data.frame(input_data[[i]], stringsAsFactors = FALSE)
      
      # Remove non-numeric columns
      non_numeric_cols <- which(vapply(mat, function(x) is.factor(x) || is.character(x), logical(1)))
      if (length(non_numeric_cols) > 0) { #print message if any columns were removed
        n_removed <- length(non_numeric_cols)
        n_total <- ncol(mat)
        if (n_removed <= 30) {
          messager(sprintf(
            "Removed %d of %d columns (variables) in dataset %s due to non-numeric type: %s",
            n_removed, n_total, name, paste(names(mat)[non_numeric_cols], collapse = ", ")
          ))
        } else {
          messager(sprintf(
            "Removed %d of %d columns (variables) in dataset %s due to non-numeric type",
            n_removed, n_total, name
          ))
        }
        mat <- mat[, -non_numeric_cols, drop = FALSE]
      }
      if (ncol(mat) == 0) stop(sprintf(
        "Data processing aborted: no columns (variables) remain in dataset %s after removing all non-numeric columns - check input data",
        name
      ))
      if (ncol(mat) == 1) {
        messager(sprintf(
          "Warning: dataset %s has only one column (variable) remaining after removing all non-numeric columns - proceeding because multiple layers are present",
          name
        ))
      }
      
      # Remove columns filled with only NA
      n_cols <- ncol(mat)
      all_na_cols <- vapply(mat, function(x) all(is.na(x)), logical(1)) #extract columns with only NA
      if (any(all_na_cols)) { #print message if any columns were removed
        n_removed <- sum(all_na_cols)
        if (n_removed <= 30) {
          messager(sprintf(
            "Removed %d of %d columns (variables) in dataset %s due to all NA: %s",
            n_removed, n_cols, name, paste(colnames(mat)[all_na_cols], collapse = ", ")
          ))
        } else {
          messager(sprintf(
            "Removed %d of %d columns (variables) in dataset %s due to all NA",
            n_removed, n_cols, name
          ))
        }
        mat <- mat[, !all_na_cols, drop = FALSE]
        if (ncol(mat) == 0) stop(sprintf(
          "Data processing aborted: no columns (variables) remain in dataset %s after removing columns with all NA - check input data",
          name
        ))
      }
      
      # Remove variables (columns) with > max.NA.col missing
      col_na_frac <- colMeans(is.na(mat))
      dropped_cols <- names(col_na_frac)[col_na_frac > max.NA.col]
      if (length(dropped_cols)) { #print message if any columns were removed
        n0 <- ncol(mat)
        n_dropped <- length(dropped_cols)
        if (n_dropped <= 30) {
          messager(sprintf(
            "Removed %d of %d columns (variables) in dataset %s due to more than %.0f%% NA in data (max.NA.col = %.2f): %s",
            n_dropped, n0, name, max.NA.col * 100, max.NA.col, paste(dropped_cols, collapse = ", ")
          ))
        } else {
          messager(sprintf(
            "Removed %d of %d columns (variables) in dataset %s due to more than %.0f%% NA in data (max.NA.col = %.2f)",
            n_dropped, n0, name, max.NA.col * 100, max.NA.col
          ))
        }
        mat <- mat[, !(colnames(mat) %in% dropped_cols), drop = FALSE]
        if (ncol(mat) == 0) {
          stop(sprintf(
            "Data processing aborted: no columns (variables) remain in dataset %s after applying max.NA.col = %.2f - check input data or increase max.NA.col",
            input_data_names[i], max.NA.col
          ))
        }
        if (ncol(mat) == 1) {
          messager(sprintf(
            "Warning: dataset %s has only one column (variable) remaining after applying max.NA.col = %.2f - proceeding because multiple layers are present",
            name, max.NA.col
          ))
        }
      }
      
      # Remove variables (columns) with zero variance
      zero_var_cols <- vapply(mat, function(x) {
        variance <- stats::var(x, na.rm = TRUE)
        is.finite(variance) && variance == 0
      }, logical(1))
      if (any(zero_var_cols)) {
        removed_cols <- names(which(zero_var_cols))
        n_removed <- length(removed_cols)
        total_cols <- ncol(mat)
        if (n_removed <= 30) {
          messager(sprintf(
            "Removed %d of %d columns (variables) in dataset %s due to zero variance: %s",
            n_removed, total_cols, name, paste(removed_cols, collapse = ", ")
          ))
        } else {
          messager(sprintf(
            "Removed %d of %d columns (variables) in dataset %s due to zero variance",
            n_removed, total_cols, name
          ))
        }
        mat <- mat[, !zero_var_cols, drop = FALSE]
        if (ncol(mat) == 0) stop(sprintf(
          "Data processing aborted: no columns (variables) remain in dataset %s after removing columns with zero variance - check input data",
          name
        ))
        if (ncol(mat) == 1) {
          messager(sprintf(
            "Warning: dataset %s has only one column (variable) remaining after removing columns with zero variance - proceeding because multiple layers are present",
            name
          ))
        }
      }
      
      # Remove all-NA rows after column filtering
      all_na_rows <- rowSums(is.na(mat)) == ncol(mat)
      if (any(all_na_rows)) { #print message if any rows are removed
        removed_names <- rownames(mat)[all_na_rows]
        n_removed <- sum(all_na_rows)
        n_rows <- nrow(mat)
        if (n_removed <= 30) {
          messager(sprintf(
            "Removed %d of %d rows (samples) in dataset %s due to all values being NA after column filtering: %s",
            n_removed, n_rows, name, paste(removed_names, collapse = ", ")
          ))
        } else {
          messager(sprintf(
            "Removed %d of %d rows (samples) in dataset %s due to all values being NA after column filtering",
            n_removed, n_rows, name
          ))
        }
        mat <- mat[!all_na_rows, , drop = FALSE]
        if (nrow(mat) == 0) {
          stop(sprintf(
            "Data processing aborted: no rows (samples) remain in dataset %s after removing all-NA rows after column filtering - check input data",
            name
          ))
        }
        if (nrow(mat) == 1) {
          stop(sprintf(
            "Data processing aborted: only one row (sample) remains in dataset %s after removing all-NA rows after column filtering - check input data",
            name
          ))
        }
      }
      
      # Add processed matrix to list for output (after filtering)
      processed_data[[length(processed_data) + 1]] <- mat
      processed_names <- c(processed_names, name)
    }
    
    # Intersect and match samples across processed matrices
    all_samples <- unique(unlist(lapply(processed_data, rownames)))
    common_samples <- Reduce(intersect, lapply(processed_data, rownames))
    if (length(common_samples) == 0) {
      stop("Data processing aborted: no shared samples remain after all-NA row filtering - check input data or increase max.NA.row")
    }
    if (length(common_samples) == 1) {
      stop("Data processing aborted: only one shared sample remains after final all-NA row filtering - check input data or increase max.NA.row")
    }
    processed_data <- lapply(processed_data, function(mat) mat[common_samples, , drop = FALSE])
    
    # Store pre-normalization layers for type detection
    type_detection_layers <- processed_data
    
    # Normalize columns for each matrix to 0-1 range
    for (i in seq_along(processed_data)) {
      input_matrix <- as.matrix(processed_data[[i]])
      input_row_names <- rownames(input_matrix)
      input_col_names <- colnames(input_matrix)
      normalized_matrix <- apply(input_matrix, 2, function(column_values) {
        column_range <- range(column_values, na.rm = TRUE)
        range_width <- diff(column_range)
        if (!all(is.finite(column_range)) || range_width == 0) return(rep(NA_real_, length(column_values)))
        (column_values - column_range[1]) / range_width
      })
      if (is.null(dim(normalized_matrix))) { #happens when ncol(input_matrix) == 1
        normalized_matrix <- matrix(normalized_matrix, ncol = 1, dimnames = list(input_row_names, input_col_names))
      } else {
        rownames(normalized_matrix) <- input_row_names
        colnames(normalized_matrix) <- input_col_names
      }
      normalized_matrix[is.nan(normalized_matrix) | is.infinite(normalized_matrix)] <- NA_real_
      if (ncol(normalized_matrix) == 0) stop(sprintf("Data processing aborted: dataset %s has no columns (variables) remaining after normalization - check input data", processed_names[i])) #check if all columns are gone after normalization
      if (nrow(normalized_matrix) == 0) stop(sprintf("Data processing aborted: dataset %s has no row (sample) remaining after normalization - check input data", processed_names[i])) #check if all rows are gone after normalization
      if (nrow(normalized_matrix) == 1) stop(sprintf("Data processing aborted: dataset %s has only one row (sample) remaining after normalization - check input data", processed_names[i])) #check if only one row remains after normalization
      
      # Remove any columns that are all NA after normalization
      mat <- as.data.frame(normalized_matrix, stringsAsFactors = FALSE)
      mat <- as.data.frame(mat, stringsAsFactors = FALSE)
      all_na_cols <- vapply(mat, function(x) all(is.na(x)), logical(1))
      if (any(all_na_cols)) { #print message if any columns were removed
        n_removed <- sum(all_na_cols)
        removed_names <- colnames(mat)[all_na_cols]
        if (n_removed <= 30) {
          messager(sprintf("Removed %d columns (variables) in dataset %s after normalization because all values are NA: %s",
                           n_removed, processed_names[i], paste(removed_names, collapse = ", ")))
        } else {
          messager(sprintf("Removed %d columns (variables) in dataset %s after normalization because all values are NA",
                           n_removed, processed_names[i]))
        }
        mat <- mat[, !all_na_cols, drop = FALSE]
        if (is.null(dim(mat)) || ncol(mat) == 0) stop("Data processing aborted: no columns (variables) remain after removing all-NA columns - check input data or increase max.NA.col")
        if (ncol(mat) == 1) {
          messager(sprintf(
            "Warning: dataset %s has only one column (variable) remaining after removing all-NA columns post-normalization - proceeding because multiple layers are present",
            processed_names[i]
          ))
        }
      }
      processed_data[[i]] <- as.matrix(mat)
    }
    
    # Restore names
    names(processed_data) <- processed_names
    input_data <- processed_data
    input_data_names <- processed_names
    
    # ... for single dataset
  } else {
    
    # Unwrap list if needed
    if (is.list(input_data) && length(input_data) == 1) {
      mat <- input_data[[1]] 
    } else {
      mat <- input_data
    }
    mat <- as.data.frame(mat, stringsAsFactors = FALSE) #convert to dataframe
    
    # Remove non-numeric columns
    non_numeric_cols <- which(sapply(mat, function(x) is.factor(x) || is.character(x)))
    if (length(non_numeric_cols) > 0) { #print message if any columns were removed
      n_removed <- length(non_numeric_cols)
      n_total <- ncol(mat) #number of columns before removal
      if (n_removed <= 30) {
        messager(sprintf(
          "Removed %d of %d columns (variables) in dataset %s due to non-numeric type: %s",
          n_removed, n_total, input_data_names[[1]], paste(names(mat)[non_numeric_cols], collapse = ", ")
        ))
      } else {
        messager(sprintf(
          "Removed %d of %d columns (variables) in dataset %s due to non-numeric type",
          n_removed, n_total, input_data_names[[1]]
        ))
      }
      mat <- mat[, -non_numeric_cols, drop = FALSE] #drop non-numeric columns
    }
    if (ncol(mat) == 0) stop("Data processing aborted: no columns (variables) remain after removing all non-numeric columns - check input data")
    if (ncol(mat) == 1) stop(sprintf(
      "Data processing aborted: dataset %s has only one column (variable) remaining after removing all non-numeric columns - check input data",
      input_data_names[[1]]
    ))
    
    # Remove rows with > max.NA.row missing
    na_props <- rowMeans(is.na(mat)) #fraction of NA per row
    bad_samples <- rownames(mat)[na_props > max.NA.row] #extract rows exceeding threshold
    if (length(bad_samples) > 0) { #print message if rows (samples) are removed
      if (length(bad_samples) <= 30) {
        messager(sprintf(
          "Removed %d of %d rows (samples) in dataset %s due to more than %.0f%% NA in data (max.NA.row = %.2f): %s",
          length(bad_samples), nrow(mat), input_data_names[[1]], max.NA.row * 100, max.NA.row, paste(bad_samples, collapse = ", ")
        ))
      } else {
        messager(sprintf(
          "Removed %d of %d rows (samples) in dataset %s due to more than %.0f%% NA in data (max.NA.row = %.2f)",
          length(bad_samples), nrow(mat), input_data_names[[1]], max.NA.row * 100, max.NA.row
        ))
      }
      mat <- mat[!rownames(mat) %in% bad_samples, , drop = FALSE] #remove rows
    }
    if (nrow(mat) == 0) {
      stop(sprintf(
        "Data processing aborted: no rows (samples) remain in dataset %s after applying max.NA.row = %.2f - check input data or increase max.NA.row",
        input_data_names[[1]], max.NA.row
      ))
    }
    if (nrow(mat) == 1) {
      stop(sprintf(
        "Data processing aborted: only one row (sample) remains in dataset %s after applying max.NA.row = %.2f - check input data or increase max.NA.row",
        input_data_names[[1]], max.NA.row
      ))
    }
    
    # Remove columns with all NA
    all_na_cols <- vapply(mat, function(x) all(is.na(x)), logical(1)) #identify columns with all NA
    if (any(all_na_cols)) { #print message if columns are removed
      n_removed <- sum(all_na_cols)
      n_total <- ncol(mat)
      removed_names <- colnames(mat)[all_na_cols]
      if (n_removed <= 30) {
        messager(sprintf("Removed %d of %d columns (variables) in dataset %s due to all NA: %s",
                         n_removed, n_total, input_data_names[[1]], paste(removed_names, collapse = ", ")))
      } else {
        messager(sprintf("Removed %d of %d columns (variables) in dataset %s due to all NA",
                         n_removed, n_total, input_data_names[[1]]))
      }
      mat <- mat[, !all_na_cols, drop = FALSE] #remove columns
    }
    if (ncol(mat) == 0) stop(paste0("Data processing aborted: dataset ", input_data_names[[1]], " has no columns (variables) remaining after removing all-NA columns - check input data")) #stop if no columns remain
    if (ncol(mat) == 1) stop(paste0("Data processing aborted: dataset ", input_data_names[[1]], " has only one column (variable) remaining after removing all-NA columns - check input data")) #stop if only one column remain
    
    # Remove columns with zero variance
    zero_var_cols <- vapply(mat, function(x) {
      variance <- stats::var(x, na.rm = TRUE)
      is.finite(variance) && variance == 0
    }, logical(1)) #extract zero-variance columns
    if (any(zero_var_cols)) { #print message if columns are removed
      removed_cols <- names(which(zero_var_cols))
      n_removed <- length(removed_cols)
      if (n_removed <= 30) {
        messager(sprintf(
          "Removed %d of %d columns (variables) in dataset %s due to zero variance: %s",
          n_removed, ncol(mat), input_data_names[[1]], paste(removed_cols, collapse = ", ")
        ))
      } else {
        messager(sprintf(
          "Removed %d of %d columns (variables) in dataset %s due to zero variance",
          n_removed, ncol(mat), input_data_names[[1]]
        ))
      }
      mat <- mat[, !zero_var_cols, drop = FALSE] #remove columns
    }
    if (ncol(mat) == 0) stop(paste0("Data processing aborted: dataset ", input_data_names[[1]], " has no columns (variables) remaining after removing columns with zero variance - check input data")) #stop if no columns (variables) remain
    if (ncol(mat) == 1) stop(paste0("Data processing aborted: dataset ", input_data_names[[1]], " has only one column (variable) remaining after removing columns with zero variance - check input data")) #stop if only one column (variable) remains
    
    # Remove columns with > max.NA.col missing
    col_na_frac <- colMeans(is.na(mat)) #fraction of NA per column
    dropped_cols <- names(col_na_frac)[col_na_frac > max.NA.col] #extract columns exceeding threshold
    if (length(dropped_cols)) { #print message if columns are removed
      n_dropped <- length(dropped_cols)
      if (n_dropped <= 30) {
        messager(sprintf(
          "Removed %d of %d columns (variables) due to more than %.0f%% NA in data (max.NA.col = %.2f): %s",
          n_dropped, ncol(mat), max.NA.col * 100, max.NA.col, paste(dropped_cols, collapse = ", ")
        ))
      } else {
        messager(sprintf(
          "Removed %d of %d columns (variables) due to more than %.0f%% NA in data (max.NA.col = %.2f)",
          n_dropped, ncol(mat), max.NA.col * 100, max.NA.col
        ))
      }
      mat <- mat[, !(colnames(mat) %in% dropped_cols), drop = FALSE] #remove columns
    }
    if (ncol(mat) == 0) {
      stop(sprintf(
        "Data processing aborted: no columns (variables) remain in dataset %s after applying max.NA.col = %.2f - check input data or increase max.NA.col",
        input_data_names[1], max.NA.col
      ))
    }
    if (ncol(mat) == 1) {
      stop(sprintf(
        "Data processing aborted: only one column (variable) remains in dataset %s after applying max.NA.col = %.2f - check input data or increase max.NA.col",
        input_data_names[1], max.NA.col
      ))
    }
    
    # Store pre-normalization layer for type detection
    type_detection_layers <- list(as.matrix(mat))
    
    # Normalize each column to 0–1 range
    mat <- apply(mat, 2, function(x) {
      col_range <- range(x, na.rm = TRUE)
      if (diff(col_range) == 0) { #avoid division by zero for constant columns
        return(rep(NA_real_, length(x)))
      }
      return((x - col_range[1]) / diff(col_range)) #apply normalization to each column
    })
    if (ncol(mat) == 0) stop(sprintf("Data processing aborted: dataset %s has no columns (variables) remaining after normalization - check input data", input_data_names[1])) #check if all columns are gone after normalization
    if (is.null(dim(mat))) stop(sprintf("Data processing aborted: dataset %s has only one column (variable) remaining after normalization - check input data", input_data_names[1])) #check if only one column remains after normalization
    if (nrow(mat) == 0) stop(sprintf("Data processing aborted: dataset %s has no row (sample) remaining after normalization - check input data", input_data_names[1])) #check if all rows are gone after normalization
    if (nrow(mat) == 1) stop(sprintf("Data processing aborted: dataset %s has only one row (sample) remaining after normalization - check input data", input_data_names[1])) #check if only one row remains after normalization
    mat <- as.data.frame(mat, stringsAsFactors = FALSE) #convert back to dataframe
    
    # Remove columns with all NA
    all_na_cols <- vapply(mat, function(x) all(is.na(x)), logical(1))
    if (any(all_na_cols)) { #print message if any columns were removed
      n_removed <- sum(all_na_cols)
      n_total <- ncol(mat)
      removed_names <- colnames(mat)[all_na_cols]
      if (n_removed <= 30) {
        messager(sprintf(
          "Removed %d of %d columns (variables) after normalization because all values are NA: %s",
          n_removed, n_total, paste(removed_names, collapse = ", ")
        ))
      } else {
        messager(sprintf(
          "Removed %d of %d columns (variables) after normalization because all values are NA",
          n_removed, n_total
        ))
      }
      mat <- mat[, !all_na_cols, drop = FALSE] #remove columns with all NA
    }
    if (ncol(mat) == 0) stop("Data processing aborted: no columns (variables) remain after removing all-NA columns post-normalization - check input data or increase max.NA.col") #print message and stop if no columns remain
    if (ncol(mat) == 1) stop("Data processing aborted: only one column (variable) remains after removing all-NA columns post-normalization - check input data or increase max.NA.col") #print message and stop if only one column remain
    
    # Remove any rows that are all NA after column filtering
    all_na_rows <- rowSums(is.na(mat)) == ncol(mat) #extract rows with all NA
    if (any(all_na_rows)) { #print message if rows with all NA are removed
      removed_names <- rownames(mat)[all_na_rows]
      n_removed <- sum(all_na_rows)
      n_total <- nrow(mat) #total number of rows before removal
      if (n_removed <= 30) {
        messager(sprintf(
          "Removed %d of %d rows (samples) with all NA after column filtering: %s",
          n_removed, n_total, paste(removed_names, collapse = ", ")
        ))
      } else {
        messager(sprintf(
          "Removed %d of %d rows (samples) with all NA after column filtering",
          n_removed, n_total
        ))
      }
      mat <- mat[!all_na_rows, , drop = FALSE] #remove rows
    }
    if (nrow(mat) == 0) stop("Data processing aborted: all rows (samples) removed after final all-NA row removal - check input data or increase max.NA.row") #stop if no rows remain
    if (nrow(mat) == 1) stop("Data processing aborted: only one row (sample) remains after final all-NA row removal - check input data or increase max.NA.row") #stop if only one row remain
    
    # Wrap matrix as list
    input_data <- list(as.matrix(mat)) #final output is list of matrix
    
  }
  
  # Infer layer types, set default distance functions (if not provided) and layer weights
  infer_layer_numeric_type <- function(layer_matrix) {
    layer_matrix <- as.matrix(layer_matrix) #force matrix
    observed_values <- stats::na.omit(as.vector(layer_matrix)) #remove all non-NA values
    if (length(observed_values) == 0) return("continuous") #fallback if layer is all NA
    unique_values <- sort(unique(observed_values)) #extract unique values in layer
    if (length(unique_values) <= 2 && all(unique_values %in% c(0, 1))) return("binary") #binary 0/1
    is_integer_like <- all(abs(observed_values - round(observed_values)) < 1e-8) #check integer-like
    min_value <- min(observed_values)
    if (is_integer_like && length(unique_values) == 3 && all(unique_values == c(0, 1, 2))) return("count") #SNP dosage 0/1/2
    if (is_integer_like && min_value >= 0) return("count") #count data
    return("continuous") #continuous data
  }
  layer_names <- names(input_data) #use list names if present
  if (is.null(layer_names) || any(layer_names == "")) layer_names <- as.character(input_data_names) #fallback to derived names
  if (is.null(manual.layer.weights)) { #set layer weights (manual.layer.weights)
    manual.layer.weights <- rep(1, length(input_data)) #equal weights
  } else {
    if (length(manual.layer.weights) == 1 && length(input_data) > 1) {
      manual.layer.weights <- rep(manual.layer.weights, length(input_data))
    }
    if (length(manual.layer.weights) != length(input_data)) {
      stop(sprintf(
        "Data processing aborted: manual.layer.weights has length %d but expected %d (number of layers)",
        length(manual.layer.weights), length(input_data)
      ))
    }
  }
  names(manual.layer.weights) <- layer_names
  layer_numeric_types <- vapply(type_detection_layers, infer_layer_numeric_type, character(1)) #binary / count / continuous
  names(layer_numeric_types) <- layer_names
  if (is.null(layer.distance.functions)) {
    layer.distance.functions <- rep("sumofsquares", length(input_data)) #default for numeric continuous layers
    names(layer.distance.functions) <- layer_names
    layer.distance.functions[layer_numeric_types == "binary"] <- "tanimoto"
    layer.distance.functions[layer_numeric_types == "binary" & 
                               grepl("SNP|COI|COII|nexus|genetic|DNA|mtDNA|cytb|VCF|ND1|genotype", layer_names, ignore.case = TRUE) & 
                               !grepl("host|plant|soil|parasite|pollinator|diet|phenology|habitat|behavior|ecology", layer_names, ignore.case = TRUE)] <- "manhattan"
    layer.distance.functions[layer_numeric_types == "count"] <- "manhattan"
    messager(sprintf("Layer types: %s", paste(sprintf("%s = %s", names(layer_numeric_types), layer_numeric_types), collapse = ", ")))
  } else {
    if (is.character(layer.distance.functions) && length(layer.distance.functions) == 1 && length(input_data) > 1) {
      layer.distance.functions <- rep(layer.distance.functions, length(input_data))
    }
    if (is.list(layer.distance.functions) && length(layer.distance.functions) == 1 && length(input_data) > 1) {
      layer.distance.functions <- rep(layer.distance.functions, length(input_data))
    }
    if (length(layer.distance.functions) != length(input_data)) {
      stop(sprintf("Data processing aborted: layer.distance.functions has length %d but expected %d (number of layers after filtering)", length(layer.distance.functions), length(input_data)))
    }
    names(layer.distance.functions) <- layer_names
  }
  
  # Report number of rows used as SOM input
  messager(sprintf("Data processing completed - using %d samples for SOM", nrow(input_data[[1]])))
  
  # Create SOM output grid
  messager("")
  messager("CREATING SOM GRID ...")
  
  # Use user-specified grid dimensions
  if (!is.null(grid.size)) {
    n_units <- prod(grid.size)
    if (n_units > nrow(input_data[[1]])) { #if number of grid cells exceeds number of samples
      stop(sprintf("Aborted SOM grid creation: requested SOM grid size %d x %d (%d units) exceeds number of samples (%d) - specify smaller grid.size", 
                   grid.size[1], grid.size[2], n_units, nrow(input_data[[1]])))
    }
    if (grid.size[1] < 2 || grid.size[2] < 2) { #if specified grid.size is too small
      stop(sprintf("Aborted SOM grid creation: custom SOM grid %d x %d is smaller than practical minimum 2 x 2 - use at least 2 x 2 for meaningful mapping", 
                   grid.size[1], grid.size[2]))
    }
    if (grid.size[1] > 50 || grid.size[2] > 50) { #if specified grid.size is very large
      messager(sprintf("Custom SOM grid %d x %d is very large! Consider smaller grid.size to avoid long training times", 
                       grid.size[1], grid.size[2]))
    }
    SOM_output_grid <- kohonen::somgrid(xdim = grid.size[1], #create SOM grid
                                        ydim = grid.size[2],
                                        topo = "hexagonal",
                                        neighbourhood.fct = training.neighborhoods)
    
    # Determine grid size and shape automatically based on number of samples and aspect ratio
  } else {
    n_samples <- nrow(input_data[[1]])
    n_units <- round(grid.multiplier * sqrt(n_samples))
    
    # Multi-layer input: use aspect ratio based on MFA (multiple factor analysis)
    if (is.list(input_data) && length(input_data) > 1 && !is.data.frame(input_data)) {
      data_list <- lapply(input_data, function(x) as.data.frame(x)) #convert each layer to data frame
      data_list <- lapply(data_list, function(df) {
        df[] <- lapply(df, function(col) { #set Inf to NA in each column
          col[is.infinite(col) | is.nan(col)] <- NA
          col
        }) 
        return(df)
      })
      MFA_data <- tryCatch({
        combined <- do.call(cbind, data_list) #combine all layers by column
        na_row_mask <- rowSums(is.na(combined)) < ncol(combined) #identify not-all-NA rows
        if (sum(na_row_mask) < 5) stop("Not enough non-NA rows for MFA - specify grid size manually") #require at least 5 rows
        combined <- combined[na_row_mask, , drop = FALSE] #keep rows with data in at least one column
        combined
      }, error = function(e) {
        messager("Aspect ratio calculation failed (not enough samples or too much NA) - using square/square-like default aspect ratio for SOM grid") #warn if failed
        return(NULL)
      })
      group_sizes <- sapply(data_list, ncol) #extract number of columns for each layer (MFA group sizes)
      if (!is.null(MFA_data) && nrow(MFA_data) >= 5 && ncol(MFA_data) >= 2) { #proceed if enough rows and columns are present
        MFA_results <- tryCatch({
          FactoMineR::MFA(MFA_data, #run MFA
                          group = group_sizes,
                          type = rep("s", length(group_sizes)),
                          ncp = 2,
                          graph = FALSE)
        }, error = function(e) {
          messager("Aspect ratio calculation failed (not enough samples or too much NA) - using square/square-like default aspect ratio for SOM grid")
          return(NULL)
        })
        if (!is.null(MFA_results)) {
          MFA_eigenvalues <- MFA_results$eig[1:2, 1] #extract first two eigenvalues from results
          if (any(is.na(MFA_eigenvalues)) || any(MFA_eigenvalues <= 0)) { #if eigenvalues invalid
            MFA_eigenvalues <- c(1, 1)
            messager("Aspect ratio calculation failed (not enough samples or too much NA) - using square/square-like default aspect ratio for SOM grid")
          }
          SOM_grid_aspect_ratio <- sqrt(MFA_eigenvalues[1] / MFA_eigenvalues[2]) #calculate aspect ratio
        } else {
          SOM_grid_aspect_ratio <- 1 #default to 1x1 if MFA fails
        }
      } else {
        SOM_grid_aspect_ratio <- 1 #default to 1x1 if not enough usable data
      }
      
      # Single-layer input: calculate SOM grid aspect ratio from input covariance matrix eigenvalues or if it fails using mean-imputation on columns with less than 25% missingness
    } else {
      mat <- as.matrix(input_data[[1]]) #convert input data to matrix
      mat[is.infinite(mat)] <- NA #set infinite values to NA
      col_good <- which(colSums(!is.na(mat)) >= 2) #keep columns with at least 2 non-NA values
      mat1 <- mat[, col_good, drop = FALSE] #subset to those columns
      row_good <- which(rowSums(!is.na(mat1)) >= 2) #keep rows with at least 2 non-NA
      mat1 <- mat1[row_good, , drop = FALSE] #subset to those rows
      mat_cc <- mat1[complete.cases(mat1), , drop = FALSE] #subset to complete-case rows
      
      if (nrow(mat_cc) >= 2 && ncol(mat_cc) >= 2) { #enough complete rows/cols for covariance
        covariance_mat <- stats::cov(mat_cc) #calculate covariance matrix
        eigenvalues <- tryCatch(eigen(covariance_mat)$values, error = function(e) c(1, 1)) #get eigenvalues or fallback
        if (length(eigenvalues) >= 2 && all(is.finite(eigenvalues)) && all(eigenvalues > 0)) {
          SOM_grid_aspect_ratio <- sqrt(eigenvalues[1] / eigenvalues[2]) #extract aspect ratio from eigenvalues
        } else {
          SOM_grid_aspect_ratio <- NA #flag to try next step
        }
      } else {
        SOM_grid_aspect_ratio <- NA #flag to try next step
      }
      
      if (is.na(SOM_grid_aspect_ratio)) { #mean-imputation step: only try if missing data is not excessive (<25% missing overall)
        col_missing_prop <- colMeans(is.na(mat)) #fraction of missing data per column
        col_good2 <- which(col_missing_prop < 0.25 & colSums(!is.na(mat)) >= 2) #columns with <25% missing and at least 2 non-NA
        row_good2 <- which(rowSums(!is.na(mat[, col_good2, drop = FALSE])) >= 2) #rows with at least 2 non-NA in these columns
        mat2 <- mat[row_good2, col_good2, drop = FALSE] #subset matrix to good rows and columns
        if (nrow(mat2) >= 2 && ncol(mat2) >= 2) { #proceed if enough data
          for (j in seq_len(ncol(mat2))) { #iterate over columns for imputation
            idx <- which(is.na(mat2[, j])) #find missing values in column
            if (length(idx) > 0 && sum(!is.na(mat2[, j])) > 0) mat2[idx, j] <- mean(mat2[, j], na.rm = TRUE) #impute mean if at least one value
          }
          covariance_mat2 <- stats::cov(mat2) #calculate covariance matrix
          eigenvalues2 <- tryCatch(eigen(covariance_mat2)$values, error = function(e) c(1, 1)) #extract eigenvalues, fallback if error
          if (length(eigenvalues2) >= 2 && all(is.finite(eigenvalues2)) && all(eigenvalues2 > 0)) {
            SOM_grid_aspect_ratio <- sqrt(eigenvalues2[1] / eigenvalues2[2]) #compute aspect ratio
          } else {
            SOM_grid_aspect_ratio <- 1 #fallback to square grid
            messager("Aspect ratio calculation failed (not enough samples or too much NA) - using square/square-like default aspect ratio for SOM grid") #warn if it fails
          }
        } else {
          SOM_grid_aspect_ratio <- 1 #fallback to square grid if not enough usable data
          messager("Aspect ratio calculation failed (not enough samples or too much NA) - using square/square-like default aspect ratio for SOM grid") #warn if it fails
        }
      }
    }
    
    # Compute ydim and xdim for closest integer fit to n_units and SOM_grid_aspect_ratio
    ydim <- max(2, round(sqrt(n_units / SOM_grid_aspect_ratio)))
    xdim <- max(2, round(n_units / ydim))
    SOM_grid_units <- xdim * ydim
    if (xdim > 50 || ydim > 50) { #large grid size
      messager(sprintf("Computed SOM grid size %d x %d is very large - consider reducing grid.multiplier or set smaller grid.size manually to avoid long training times", 
                       xdim, ydim))
    }
    if (SOM_grid_units > n_samples) { #grid cell number less than number of samples
      stop(sprintf("Aborted SOM grid creation: SOM grid size %d x %d (%d units) exceeds number of samples (%d) - decrease grid.multiplier or set grid.size manually", 
                   xdim, ydim, SOM_grid_units, n_samples))
    }
    if (xdim < 2 || ydim < 2) { #small grid size
      stop(sprintf("Aborted SOM grid creation: SOM grid size %d x %d is smaller than 2 x 2 - increase grid.multiplier or grid.size", 
                   xdim, ydim))
    }
    SOM_output_grid <- kohonen::somgrid(xdim = xdim, #create SOM grid
                                        ydim = ydim,
                                        topo = "hexagonal",
                                        neighbourhood.fct = training.neighborhoods)
  }
  
  # Report grid dimensions
  SOM_grid_x <- SOM_output_grid$xdim
  SOM_grid_y <- SOM_output_grid$ydim
  SOM_grid_cells <- SOM_grid_x * SOM_grid_y
  messager(sprintf(
    "Created %d x %d SOM grid (%d cells)", #message showing grid size and total cells
    SOM_grid_x, SOM_grid_y, SOM_grid_cells
  ))
  
  # Set neighborhood radius schedule
  nhbrdist <- kohonen::unit.distances(SOM_output_grid)
  radius.schedule <- c(as.numeric(stats::quantile(nhbrdist, 2/3, na.rm = TRUE)), 0)
  
  # Perform learning rate tuning
  if (learning.rate.tuning) {
    messager("")
    messager("TUNING LEARNING RATES ...")
    
    # Set tuning values and number of runs
    learning_rate_initial_tuning_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    learning_rate_final_tuning_values <- c(0.001, 0.01, 0.1, 0.3, 0.5)
    learning_rate_N_steps <- 50
    learning_rate_N_replicates <- 10
    tuning_values <- expand.grid(learning_rate_initial = learning_rate_initial_tuning_values,
                                 learning_rate_final = learning_rate_final_tuning_values,
                                 stringsAsFactors = FALSE)
    tuning_values <- tuning_values[tuning_values$learning_rate_final < tuning_values$learning_rate_initial, , drop = FALSE] #keep only valid combinations: final learning rate < initial learning rate
    tuning_values$qe_mean <- NA_real_
    tuning_values$qe_sd <- NA_real_
    
    # Run learning_rate_N_replicates replicates and record mean quantization error QE (average distance from each data point to its closest SOM node) each time
    for (i in seq_len(nrow(tuning_values))) {
      learning_rate_initial <- tuning_values$learning_rate_initial[i]
      learning_rate_final <- tuning_values$learning_rate_final[i]
      base::set.seed(i + set.seed.N)
      qes <- replicate(learning_rate_N_replicates, {
        som_tuning <- kohonen::supersom(data = input_data,
                                        grid = SOM_output_grid,
                                        maxNA.fraction = max.NA.row,
                                        alpha = c(learning_rate_initial, learning_rate_final),
                                        rlen = learning_rate_N_steps,
                                        dist.fcts = layer.distance.functions,
                                        mode = "online",
                                        keep.data = TRUE, 
                                        whatmap = NULL,
                                        radius = radius.schedule,
                                        user.weights = manual.layer.weights)
        mean(som_tuning$distances, na.rm = TRUE)
      })
      tuning_values$qe_mean[i] <- mean(qes, na.rm = TRUE)
      tuning_values$qe_sd[i] <- stats::sd(qes, na.rm = TRUE)
    }
    tuning_values_sorted <- tuning_values[order(tuning_values$qe_mean), ]
    tuning_values_best <- tuning_values_sorted[1, ]
    messager(sprintf("Best tuning parameters: initial = %s, final = %s",
                     tuning_values_best$learning_rate_initial, tuning_values_best$learning_rate_final
    ))
    learning.rate.initial <- tuning_values_best$learning_rate_initial
    learning.rate.final <- tuning_values_best$learning_rate_final
  }
  
  # Create function to calculate topographic error
  calculate.topographic.error <- function(som_model) {
    codes <- kohonen::getCodes(som_model)
    if (!is.list(codes)) codes <- list(codes)
    codebook_matrix <- do.call(cbind, lapply(codes, as.matrix))
    sample_matrix <- do.call(cbind, lapply(som_model$data, as.matrix))
    unit_distance_matrix <- kohonen::unit.distances(som_model$grid)
    adjacency_matrix <- unit_distance_matrix <= 1
    diag(adjacency_matrix) <- FALSE
    topographic_error_vector <- apply(sample_matrix, 1, function(sample_values) {
      distance_to_codes <- apply(codebook_matrix, 1, function(code_values) {
        sqrt(sum((sample_values - code_values)^2, na.rm = TRUE))
      })
      best_two_units <- order(distance_to_codes)[1:2]
      if (adjacency_matrix[best_two_units[1], best_two_units[2]]) 0 else 1
    })
    mean(topographic_error_vector, na.rm = TRUE)
  }
  
  # Create function to run SOM
  messager("")
  messager("TRAINING SOM ...")
  replicate_som <- function(j) {
    
    # Initialize results for replicate
    d_vec <- numeric(length(input_data))
    
    # Print message every N replicates (N specified by message.N.replicates)
    if (j %% message.N.replicates == 0) {
      messager(paste("Running replicate:", j, "of", N.replicates))
    }
    
    # Run SOM model
    som_model <- kohonen::supersom(data = input_data, 
                                   grid = SOM_output_grid, 
                                   maxNA.fraction = max.NA.row, 
                                   alpha = c(learning.rate.initial, learning.rate.final), 
                                   rlen = N.steps,
                                   mode = "online",
                                   whatmap = NULL,
                                   keep.data = TRUE, 
                                   radius = radius.schedule,
                                   dist.fcts = layer.distance.functions,
                                   user.weights = manual.layer.weights)
    
    # Store learning values for each matrix
    learning_values_list <- lapply(seq_along(input_data), function(i) som_model$changes[, i])
    
    # Store distance weights
    d_vec <- som_model$distance.weights 
    
    # Store quantization error
    quantization_error <- mean(som_model$distances, na.rm = TRUE)
    
    # Store topographic error
    topographic_error <- calculate.topographic.error(som_model)
    
    # Return results
    return(list(
      d_vec = d_vec, 
      learning_values_list = learning_values_list, 
      som_model = som_model,
      quantization_error = quantization_error,
      topographic_error = topographic_error
    ))
  }
  
  # Perform SOM training and collect results for all replicates
  
  # Run in parallel
  if (parallel) { 
    messager(sprintf("Running SOM in parallel with %d cores", N.cores))
    required_packages_parallel <- c(
      "kohonen",
      "foreach",
      "doParallel",
      "doRNG"
    )
    parallel_cluster <- parallel::makeCluster(N.cores) #start PSOCK cluster
    on.exit(parallel::stopCluster(parallel_cluster), add = TRUE)
    parallel::clusterExport( #export helper and all needed globals to each worker
      parallel_cluster,
      varlist = c("messager",
                  "replicate_som",
                  "calculate.topographic.error",
                  "input_data",
                  "SOM_output_grid",
                  "N.replicates",
                  "N.steps",
                  "max.NA.row",
                  "learning.rate.initial",
                  "learning.rate.final",
                  "layer.distance.functions",
                  "manual.layer.weights",
                  "radius.schedule",
                  "message.N.replicates"),
      envir = environment())
    doParallel::registerDoParallel(parallel_cluster) #register cluster for foreach
    doRNG::registerDoRNG(seed = set.seed.N) #set seed
    results <- foreach::`%dopar%`(
      foreach::foreach(
        j = seq_len(N.replicates),
        .packages = required_packages_parallel
      ),
      {
        tryCatch(
          replicate_som(j),
          error = function(e) {
            list(
              d_vec = rep(NA_real_, length(input_data)),
              learning_values_list = lapply(seq_along(input_data), function(k) rep(NA_real_, N.steps)),
              som_model = NULL,
              quantization_error = NA_real_,
              topographic_error = NA_real_,
              error = conditionMessage(e)
            )
          }
        )
      }
    )
  } else {
    
    # Run SOM normally (non-parallel)
    base::set.seed(set.seed.N)
    results <- tryCatch(
      lapply(seq_len(N.replicates), function(j) {
        base::set.seed(j + set.seed.N) # set unique seed for each replicate
        replicate_som(j)
      }),
      error = function(e) { #print error message if SOM training fails
        stop("SOM training aborted: try reducing grid.size or increasing max.NA.col and max.NA.row or check input data")
      }
    )
    if (is.null(results)) return(invisible(NULL))
  }
  
  # Drop failed replicates and preserve original replicate IDs
  original_replicate_ids <- paste0("R", seq_len(length(results)))
  names(results) <- original_replicate_ids
  failed <- vapply(results, function(x) !is.null(x$error) || is.null(x$som_model), logical(1))
  failed_replicate_ids <- names(results)[failed]
  if (any(failed)) {
    messager(sprintf("Dropped %d replicate(s) that failed during training", sum(failed)))
    results <- results[!failed]
  }
  if (length(results) < 1) stop("SOM training aborted: all replicates failed - check input data or relax NA thresholds / grid size")
  retained_replicate_ids <- names(results)
  
  # Combine results from all replicates
  distance_weights_matrix <- do.call(rbind, lapply(results, `[[`, "d_vec"))
  rownames(distance_weights_matrix) <- retained_replicate_ids
  colnames(distance_weights_matrix) <- as.character(input_data_names)
  
  learning_values_list <- lapply(1:length(input_data), function(i) {
    learning_values_combined <- do.call(cbind, lapply(results, function(res) res$learning_values_list[[i]]))
    rownames(learning_values_combined) <- paste0("S", seq_len(N.steps))
    colnames(learning_values_combined) <- retained_replicate_ids
    return(learning_values_combined)
  })
  
  quantization_error <- vapply(results, `[[`, numeric(1), "quantization_error")
  names(quantization_error) <- retained_replicate_ids
  
  topographic_error <- vapply(results, `[[`, numeric(1), "topographic_error")
  names(topographic_error) <- retained_replicate_ids
  
  som_models <- lapply(results, `[[`, "som_model")
  names(som_models) <- retained_replicate_ids
  
  # Compute codebook vectors for each layer across replicates
  codes0 <- kohonen::getCodes(som_models[[1]]) #get codes from first model
  if (is.list(codes0)) { #supersom multi-layer returns list
    n_layers <- length(codes0)
  } else { #single-layer may return matrix
    n_layers <- 1
  }
  all_layer_codes <- vector("list", n_layers)
  for (l in seq_len(n_layers)) {
    all_layer_codes[[l]] <- do.call(
      rbind,
      lapply(som_models, function(m) {
        codes_m <- kohonen::getCodes(m) #extract codes
        if (!is.list(codes_m)) return(codes_m) #single-layer: matrix
        codes_m[[l]] #multi-layer: list element
      })
    )
  }
  
  # Save results
  SOM_results <- list(
    distance_weights_matrix = distance_weights_matrix, 
    learning_values_list = learning_values_list, 
    input_data_names = as.character(input_data_names),
    N_steps = N.steps,
    N_replicates = length(som_models),
    N_replicates_requested = N.replicates,
    N_replicates_failed = length(failed_replicate_ids),
    retained_replicate_ids = retained_replicate_ids,
    failed_replicate_ids = failed_replicate_ids,
    learning_rate_initial = learning.rate.initial,
    learning_rate_final = learning.rate.final,
    codebook_vectors = all_layer_codes,
    som_models = som_models,
    layer_numeric_types = layer_numeric_types,
    layer.distance.functions = layer.distance.functions,
    quantization_error = quantization_error,
    topographic_error = topographic_error,
    input_data = input_data,
    train.SOM.set.seed.N = set.seed.N,
    train.SOM.args = list(
      N.steps = N.steps,
      N.replicates = length(som_models),
      N.replicates.requested = N.replicates,
      N.replicates.failed = length(failed_replicate_ids),
      parallel = parallel,
      N.cores = N.cores,
      grid.size = grid.size,
      grid.multiplier = grid.multiplier,
      learning.rate.initial = learning.rate.initial,
      learning.rate.final = learning.rate.final,
      learning.rate.tuning = learning.rate.tuning,
      layer.distance.functions = layer.distance.functions,
      manual.layer.weights = manual.layer.weights,
      max.NA.row = max.NA.row,
      max.NA.col = max.NA.col,
      training.neighborhoods = training.neighborhoods,
      save.SOM.results = save.SOM.results,
      save.SOM.results.name = save.SOM.results.name,
      overwrite.SOM.results = overwrite.SOM.results,
      verbose = verbose,
      message.N.replicates = message.N.replicates,
      set.seed.N = set.seed.N
    )
  )
  
  # Save results
  if (save.SOM.results) {
    
    # Check if directory exists
    dir_path <- dirname(save.SOM.results.name) #extract directory path
    if (!dir.exists(dir_path)) { 
      dir.create(dir_path, recursive = T) #create directory if it does not exist
      messager(paste("Specified directory", dir_path, "did not exist and was created"))
    }
    
    # Save results
    save(SOM_results, file = save.SOM.results.name)
    if (save.SOM.results && !overwrite.SOM.results) {
      messager("SOM results saved as ", save.SOM.results.name)
    }
    if (save.SOM.results && overwrite.SOM.results) {
      messager("SOM results overwritten as ", save.SOM.results.name)
    }
  }
  messager("")
  messager("FINISHED SUCCESSFULLY")
  return(SOM_results)
}


## Function to cluster SOM codebook vectors
clustering.SOM <- function(SOM.output,
                           max.k = 10, #maximum of considered clusters K + 1
                           set.k = NULL, #set to test single value of K
                           clustering.method, #set clustering method
                           BIC.thresh = 6, #BIC threshold for selecting K > 1 - we suggest using Raftery (1995) ranges: 2, 6, or 10 for weak, medium or strong support
                           quantization.error.quantile = 0.95, #remove mappings with QE above this quantile; NULL = no filtering
                           topographic.error.quantile = 0.95, #remove mappings with TE above this quantile; NULL = no filtering
                           set.seed.N = 1
) {
  
  
  # Validate specified SOM.output
  if (!is.list(SOM.output) || is.null(SOM.output$som_models) || length(SOM.output$som_models) < 1) {
    stop("Aborted SOM clustering: SOM.output must be list from train.SOM() with non-empty $som_models")
  }
  
  # Validate specified max.k
  if (!is.numeric(max.k) || length(max.k) != 1 || is.na(max.k) || max.k < 2 || (max.k %% 1 != 0)) {
    stop("Aborted SOM clustering: max.k must be a single integer >= 2")
  }
  
  # Validate specified set.k
  if (!is.null(set.k) && (!is.numeric(set.k) || length(set.k) != 1 || is.na(set.k) || set.k < 1 || (set.k %% 1 != 0))) {
    stop("Aborted SOM clustering: set.k must be NULL or single positive integer >= 1")
  }
  if (!is.null(set.k) && set.k > max.k) {
    stop("Aborted SOM clustering: set.k must be <= max.k")
  }
  
  # Validate specified clustering method
  valid.methods <- c(
    "kmeans+BICelbow",
    "kmeans+BICthreshold",
    "GMM+BICthreshold",
    "hierarchical+DB",
    "HDBSCAN",
    "OPTICS+Silhouette"
  )
  if (!clustering.method %in% valid.methods) {
    stop("Aborted SOM clustering: Invalid clustering.method - must be one of ", paste(valid.methods, collapse = ", "))
  }
  
  # Validate specified BIC.thresh
  if (!is.numeric(BIC.thresh) || length(BIC.thresh) != 1 || is.na(BIC.thresh) || BIC.thresh <= 0) {
    stop("Aborted SOM clustering: BIC.thresh must be a single positive numeric value (e.g., 2, 6, or 10 for low, moderate or strong support, respectively)")
  }
  
  # Validate specified quantization.error.quantile
  if (!is.null(quantization.error.quantile)) {
    if (!is.numeric(quantization.error.quantile) || length(quantization.error.quantile) != 1 || is.na(quantization.error.quantile) ||
        quantization.error.quantile <= 0 || quantization.error.quantile >= 1) {
      stop("Aborted SOM clustering: quantization.error.quantile must be NULL or a single numeric value in (0, 1)")
    }
  }
  
  # Validate specified topographic.error.quantile
  if (!is.null(topographic.error.quantile)) {
    if (!is.numeric(topographic.error.quantile) || length(topographic.error.quantile) != 1 || is.na(topographic.error.quantile) ||
        topographic.error.quantile <= 0 || topographic.error.quantile >= 1) {
      stop("Aborted SOM clustering: topographic.error.quantile must be NULL or a single numeric value in (0, 1)")
    }
  }
  
  # Create function to find nearest neighbors (following FNN get.knnx function)
  get.knnx.custom <- function(reference_data, query_data, k = 1) {
    if (k != 1) stop("get.knnx.custom currently supports only k = 1")
    reference_data <- as.matrix(reference_data) #ensure matrix input
    query_data <- as.matrix(query_data) #ensure matrix input
    if (!is.numeric(reference_data) || !is.numeric(query_data)) stop("Data non-numeric") #require numeric matrices
    if (anyNA(reference_data) || anyNA(query_data)) stop("Data include NAs") #NA values are not allowed
    if (ncol(reference_data) != ncol(query_data)) stop("Number of columns must be same!") #dimensions must match
    if (nrow(reference_data) == 0) stop("reference_data must contain at least one row") #need at least one candidate neighbor
    if (nrow(query_data) == 0) { #return empty result if there are no query points
      return(list(
        nn.index = matrix(integer(0), ncol = 1),
        nn.dist = matrix(numeric(0), ncol = 1)
      ))
    }
    nearest_neighbor_indices <- integer(nrow(query_data)) #store row index of nearest reference point for each query point
    nearest_neighbor_distances <- numeric(nrow(query_data)) #store distance to nearest reference point for each query point
    for (query_row_index in seq_len(nrow(query_data))) {
      query_point <- query_data[query_row_index, ] #extract current query point
      euclidean_distances <- sqrt(rowSums((sweep(reference_data, 2, query_point, FUN = "-"))^2)) #calculate Euclidean distance to all reference points
      nearest_reference_index <- which.min(euclidean_distances) #find closest reference point
      nearest_neighbor_indices[query_row_index] <- nearest_reference_index #average nearest neighbor index
      nearest_neighbor_distances[query_row_index] <- euclidean_distances[nearest_reference_index] #save nearest neighbor distance
    }
    return(list(
      nn.index = matrix(nearest_neighbor_indices, ncol = 1), #match FNN output structure
      nn.dist = matrix(nearest_neighbor_distances, ncol = 1) #match FNN output structure
    ))
  }
  
  # Create function to calculate layer-wise sample-to-unit distances
  compute.layer.sample.to.unit.distance.SOM <- function(sample_matrix,
                                                        codebook_matrix,
                                                        distance_function = "sumofsquares"
  ) {
    
    # Validate specified matrices
    sample_matrix <- as.matrix(sample_matrix)
    codebook_matrix <- as.matrix(codebook_matrix)
    storage.mode(sample_matrix) <- "numeric"
    storage.mode(codebook_matrix) <- "numeric"
    
    if (nrow(sample_matrix) == 0 || ncol(sample_matrix) == 0) {
      stop("Soft assignment calculation aborted: sample_matrix is empty")
    }
    if (nrow(codebook_matrix) == 0 || ncol(codebook_matrix) == 0) {
      stop("Soft assignment calculation aborted: codebook_matrix is empty")
    }
    if (ncol(sample_matrix) != ncol(codebook_matrix)) {
      stop("Soft assignment calculation aborted: sample_matrix and codebook_matrix must have identical numbers of columns")
    }
    
    # Calculate distances
    distance_matrix <- matrix(NA_real_,
                              nrow = nrow(sample_matrix),
                              ncol = nrow(codebook_matrix),
                              dimnames = list(rownames(sample_matrix),
                                              rownames(codebook_matrix)))
    
    for (unit_index in seq_len(nrow(codebook_matrix))) {
      current_codebook_vector <- codebook_matrix[unit_index, ]
      current_distance_vector <- apply(sample_matrix, 1, function(sample_vector) {
        valid_values <- is.finite(sample_vector) & !is.na(sample_vector) & is.finite(current_codebook_vector) & !is.na(current_codebook_vector)
        if (!any(valid_values)) return(NA_real_)
        sample_vector <- sample_vector[valid_values]
        current_codebook_vector_valid <- current_codebook_vector[valid_values]
        
        if (distance_function == "sumofsquares") {
          return(sqrt(sum((sample_vector - current_codebook_vector_valid) ^ 2)))
        }
        if (distance_function == "manhattan") {
          return(sum(abs(sample_vector - current_codebook_vector_valid)))
        }
        if (distance_function == "tanimoto") {
          sample_binary <- as.numeric(sample_vector > 0.5)
          codebook_binary <- as.numeric(current_codebook_vector_valid > 0.5)
          shared_one <- sum(sample_binary == 1 & codebook_binary == 1)
          any_one <- sum(sample_binary == 1 | codebook_binary == 1)
          if (any_one == 0) return(0)
          return(1 - (shared_one / any_one))
        }
        
        stop("Soft assignment calculation aborted: unsupported distance function in compute.layer.sample.to.unit.distance.SOM")
      })
      distance_matrix[, unit_index] <- current_distance_vector
    }
    
    return(distance_matrix)
  }
  
  
  # Create function to calculate integrated sample-to-unit distances across SOM layers
  compute.sample.to.unit.distance.matrix.SOM <- function(som_model,
                                                         layer.distance.functions,
                                                         layer.weights = NULL
  ) {
    
    # Validate specified som_model
    if (is.null(som_model) || is.null(som_model$data)) {
      stop("Soft assignment calculation aborted: som_model does not contain data")
    }
    
    # Extract sample data and codebook vectors
    sample_data_list <- som_model$data
    if (!is.list(sample_data_list)) {
      sample_data_list <- list(sample_data_list)
    }
    codebook_list <- kohonen::getCodes(som_model)
    if (!is.list(codebook_list)) {
      codebook_list <- list(codebook_list)
    }
    
    if (length(sample_data_list) != length(codebook_list)) {
      stop("Soft assignment calculation aborted: number of SOM data layers does not match number of codebook layers")
    }
    
    # Validate specified layer.distance.functions
    if (is.null(layer.distance.functions)) {
      layer.distance.functions <- rep("sumofsquares", length(sample_data_list))
    }
    if (length(layer.distance.functions) == 1 && length(sample_data_list) > 1) {
      layer.distance.functions <- rep(layer.distance.functions, length(sample_data_list))
    }
    if (length(layer.distance.functions) != length(sample_data_list)) {
      stop("Soft assignment calculation aborted: layer.distance.functions does not match number of SOM layers")
    }
    
    # Set layer weights
    if (is.null(layer.weights)) {
      if (!is.null(som_model$distance.weights)) {
        layer.weights <- as.numeric(som_model$distance.weights)
      } else {
        layer.weights <- rep(1, length(sample_data_list))
      }
    }
    if (length(layer.weights) == 1 && length(sample_data_list) > 1) {
      layer.weights <- rep(layer.weights, length(sample_data_list))
    }
    if (length(layer.weights) != length(sample_data_list)) {
      stop("Soft assignment calculation aborted: layer.weights does not match number of SOM layers")
    }
    
    # Calculate layer-specific distance matrices
    layer_distance_matrices <- vector("list", length(sample_data_list))
    for (layer_index in seq_along(sample_data_list)) {
      sample_matrix <- as.matrix(sample_data_list[[layer_index]])
      codebook_matrix <- as.matrix(codebook_list[[layer_index]])
      rownames(codebook_matrix) <- seq_len(nrow(codebook_matrix))
      layer_distance_matrices[[layer_index]] <- compute.layer.sample.to.unit.distance.SOM(
        sample_matrix = sample_matrix,
        codebook_matrix = codebook_matrix,
        distance_function = layer.distance.functions[layer_index]
      )
    }
    
    # Combine layer-specific distances
    integrated_distance_matrix <- matrix(0,
                                         nrow = nrow(layer_distance_matrices[[1]]),
                                         ncol = ncol(layer_distance_matrices[[1]]),
                                         dimnames = dimnames(layer_distance_matrices[[1]]))
    valid_weight_matrix <- matrix(0,
                                  nrow = nrow(integrated_distance_matrix),
                                  ncol = ncol(integrated_distance_matrix))
    
    for (layer_index in seq_along(layer_distance_matrices)) {
      current_distance_matrix <- layer_distance_matrices[[layer_index]]
      current_weight <- layer.weights[layer_index]
      valid_entries <- is.finite(current_distance_matrix) & !is.na(current_distance_matrix)
      integrated_distance_matrix[valid_entries] <- integrated_distance_matrix[valid_entries] + (current_distance_matrix[valid_entries] * current_weight)
      valid_weight_matrix[valid_entries] <- valid_weight_matrix[valid_entries] + current_weight
    }
    
    integrated_distance_matrix[valid_weight_matrix > 0] <- integrated_distance_matrix[valid_weight_matrix > 0] / valid_weight_matrix[valid_weight_matrix > 0]
    integrated_distance_matrix[valid_weight_matrix == 0] <- NA_real_
    
    return(integrated_distance_matrix)
  }
  
  
  # Create function to calculate sample-level soft cluster probabilities from sample-to-unit distances
  compute.replicate.ancestry.matrix.SOM <- function(sample_to_unit_distance_matrix,
                                                    unit_cluster_labels,
                                                    temperature = NULL
  ) {
    
    # Validate specified sample_to_unit_distance_matrix
    if (is.null(sample_to_unit_distance_matrix)) {
      stop("Replicate ancestry calculation aborted: sample_to_unit_distance_matrix is NULL")
    }
    sample_to_unit_distance_matrix <- as.matrix(sample_to_unit_distance_matrix)
    storage.mode(sample_to_unit_distance_matrix) <- "numeric"
    if (nrow(sample_to_unit_distance_matrix) == 0 || ncol(sample_to_unit_distance_matrix) == 0) {
      stop("Replicate ancestry calculation aborted: sample_to_unit_distance_matrix is empty")
    }
    
    # Validate specified unit_cluster_labels
    if (is.null(unit_cluster_labels) || length(unit_cluster_labels) != ncol(sample_to_unit_distance_matrix)) {
      stop("Replicate ancestry calculation aborted: unit_cluster_labels must have one label per SOM unit")
    }
    unit_cluster_labels <- as.character(unit_cluster_labels)
    unique_cluster_labels <- sort(unique(unit_cluster_labels))
    
    # Return one-column ancestry matrix if only one cluster is present
    if (length(unique_cluster_labels) == 1) {
      replicate_ancestry_matrix <- matrix(1,
                                          nrow = nrow(sample_to_unit_distance_matrix),
                                          ncol = 1,
                                          dimnames = list(rownames(sample_to_unit_distance_matrix),
                                                          unique_cluster_labels))
      return(replicate_ancestry_matrix)
    }
    
    # Choose temperature if not specified
    if (is.null(temperature)) {
      finite_distances <- as.numeric(sample_to_unit_distance_matrix)
      finite_distances <- finite_distances[is.finite(finite_distances) & !is.na(finite_distances)]
      if (length(finite_distances) == 0) {
        stop("Replicate ancestry calculation aborted: no finite sample-to-unit distances found")
      }
      temperature <- stats::median(finite_distances)
      if (!is.finite(temperature) || temperature <= 0) {
        temperature <- 1
      }
    }
    
    # Convert distances to soft unit weights
    unit_weight_matrix <- exp(-(sample_to_unit_distance_matrix ^ 2) / (2 * temperature ^ 2))
    unit_weight_row_sums <- rowSums(unit_weight_matrix, na.rm = TRUE)
    valid_rows <- is.finite(unit_weight_row_sums) & unit_weight_row_sums > 0
    if (!any(valid_rows)) {
      stop("Replicate ancestry calculation aborted: all unit weight row sums are invalid")
    }
    unit_weight_matrix[valid_rows, ] <- unit_weight_matrix[valid_rows, , drop = FALSE] / unit_weight_row_sums[valid_rows]
    
    # Aggregate unit weights to cluster probabilities
    replicate_ancestry_matrix <- sapply(unique_cluster_labels, function(current_cluster_label) {
      current_cluster_unit_indices <- which(unit_cluster_labels == current_cluster_label)
      rowSums(unit_weight_matrix[, current_cluster_unit_indices, drop = FALSE], na.rm = TRUE)
    })
    replicate_ancestry_matrix <- as.matrix(replicate_ancestry_matrix)
    
    # Restore row and column names
    rownames(replicate_ancestry_matrix) <- rownames(sample_to_unit_distance_matrix)
    colnames(replicate_ancestry_matrix) <- unique_cluster_labels
    
    # Renormalize cluster probabilities
    ancestry_row_sums <- rowSums(replicate_ancestry_matrix, na.rm = TRUE)
    valid_rows <- is.finite(ancestry_row_sums) & ancestry_row_sums > 0
    replicate_ancestry_matrix[valid_rows, ] <- replicate_ancestry_matrix[valid_rows, , drop = FALSE] / ancestry_row_sums[valid_rows]
    
    return(replicate_ancestry_matrix)
  }
  
  # Create function to calculate mean assignment margin
  calculate.mean.assignment.margin.SOM <- function(assignment_probability_matrix) {
    
    # Return missing if matrix is unavailable
    if (is.null(assignment_probability_matrix)) return(NA_real_)
    
    # k = 1 provides no meaningful separation metric
    if (ncol(assignment_probability_matrix) <= 1) {
      return(NA_real_)
    }
    
    # Calculate mean assignment margin
    row_sorted_probabilities <- t(apply(assignment_probability_matrix, 1, sort, decreasing = TRUE))
    mean_assignment_margin <- mean(row_sorted_probabilities[, 1] - row_sorted_probabilities[, 2], na.rm = TRUE)
    
    return(mean_assignment_margin)
  }
  
  
  # Create function to calculate mean normalized assignment entropy
  calculate.mean.normalized.assignment.entropy.SOM <- function(assignment_probability_matrix) {
    
    # Return missing if matrix is unavailable
    if (is.null(assignment_probability_matrix)) return(NA_real_)
    
    # k = 1 provides no meaningful separation metric
    if (ncol(assignment_probability_matrix) <= 1) {
      return(NA_real_)
    }
    
    # Calculate mean normalized assignment entropy
    safe_assignment_probability_matrix <- pmax(assignment_probability_matrix, .Machine$double.eps)
    row_entropies <- -rowSums(safe_assignment_probability_matrix * log(safe_assignment_probability_matrix), na.rm = TRUE)
    normalized_row_entropies <- row_entropies / log(ncol(assignment_probability_matrix))
    mean_normalized_assignment_entropy <- mean(normalized_row_entropies, na.rm = TRUE)
    
    # Return mean normalized assignment entropy
    return(mean_normalized_assignment_entropy)
  }
  
  # Extract number of replicates
  N.replicates <- length(SOM.output$som_models)
  
  # Filter poorly fitting SOM mappings based on quantization error and/or topographic error
  replicate_ids <- names(SOM.output$som_models)
  if (is.null(replicate_ids) || any(replicate_ids == "")) {
    replicate_ids <- paste0("R", seq_len(length(SOM.output$som_models)))
  }
  retained_replicates <- replicate_ids
  removed_replicates <- character(0)
  if (!is.null(quantization.error.quantile) || !is.null(topographic.error.quantile)) {
    replicates_to_remove_qe <- integer(0)
    replicates_to_remove_te <- integer(0)
    if (!is.null(quantization.error.quantile)) {
      if (is.null(SOM.output$quantization_error)) {
        stop("Aborted SOM clustering: SOM.output does not contain quantization_error - rerun train.SOM")
      }
      quantization_error <- SOM.output$quantization_error
      if (length(quantization_error) != length(SOM.output$som_models)) {
        stop("Aborted SOM clustering: length of quantization_error does not match number of som_models - rerun train.SOM")
      }
      if (any(!is.finite(quantization_error) | is.na(quantization_error))) {
        stop("Aborted SOM clustering: quantization_error contains NA or non-finite values")
      }
      quantization_error_cutoff <- stats::quantile(quantization_error,
                                                   probs = quantization.error.quantile,
                                                   na.rm = TRUE,
                                                   names = FALSE)
      replicates_to_remove_qe <- which(quantization_error > quantization_error_cutoff)
    }
    
    if (!is.null(topographic.error.quantile)) {
      if (is.null(SOM.output$topographic_error)) {
        stop("Aborted SOM clustering: SOM.output does not contain topographic_error - rerun train.SOM")
      }
      topographic_error <- SOM.output$topographic_error
      if (length(topographic_error) != length(SOM.output$som_models)) {
        stop("Aborted SOM clustering: length of topographic_error does not match number of som_models - rerun train.SOM")
      }
      if (any(!is.finite(topographic_error) | is.na(topographic_error))) {
        stop("Aborted SOM clustering: topographic_error contains NA or non-finite values")
      }
      topographic_error_cutoff <- stats::quantile(topographic_error,
                                                  probs = topographic.error.quantile,
                                                  na.rm = TRUE,
                                                  names = FALSE)
      replicates_to_remove_te <- which(topographic_error > topographic_error_cutoff)
    }
    removed_index <- sort(unique(c(replicates_to_remove_qe, replicates_to_remove_te)))
    removed_replicates <- replicate_ids[removed_index]
    retained_replicates <- setdiff(replicate_ids, removed_replicates)
    if (length(retained_replicates) < 2) {
      stop("Aborted SOM clustering: filtering removed too many SOM replicates - relax quantile thresholds")
    }
    retained_index <- match(retained_replicates, replicate_ids)
    SOM.output$som_models <- SOM.output$som_models[retained_index]
    if (!is.null(SOM.output$quantization_error)) {
      SOM.output$quantization_error <- SOM.output$quantization_error[retained_index]
    }
    if (!is.null(SOM.output$topographic_error)) {
      SOM.output$topographic_error <- SOM.output$topographic_error[retained_index]
    }
    if (!is.null(SOM.output$distance_weights_matrix)) {
      SOM.output$distance_weights_matrix <- SOM.output$distance_weights_matrix[retained_index, , drop = FALSE]
    }
    if (!is.null(SOM.output$learning_values_list)) {
      SOM.output$learning_values_list <- lapply(SOM.output$learning_values_list, function(x) {
        x[, retained_index, drop = FALSE]
      })
    }
  }
  
  # Update number of replicates after optional filtering
  N.replicates <- length(SOM.output$som_models)
  
  # Create function to cluster SOM models
  replicate_clust <- function(j) {
    
    # Set seed
    base::set.seed(j + set.seed.N)
    
    # Extract SOM models
    som_model <- SOM.output$som_models[[j]] 
    
    # Extract  SOM codebook vectors
    codes <- kohonen::getCodes(som_model)
    if (!is.list(codes)) codes <- list(codes)
    som_codes <- do.call(cbind, codes)
    rownames(som_codes) <- paste0("G", seq_len(nrow(som_codes)))
    
    support_vec <- rep(NA_real_, max.k) #method-specific support values across k
    support_label <- NULL #name of support metric
    support_higher_is_better <- NA #whether larger values indicate better support
    
    # Ensure som_codes is finite
    som_codes <- as.matrix(som_codes)
    invalid_codebook_entries <- !is.finite(som_codes) | is.na(som_codes)
    if (any(invalid_codebook_entries)) {
      for (variable_index in seq_len(ncol(som_codes))) {
        column_values <- som_codes[, variable_index]
        column_mean_value <- mean(column_values[is.finite(column_values)], na.rm = TRUE)
        if (!is.finite(column_mean_value)) column_mean_value <- 0.5
        column_values[!is.finite(column_values) | is.na(column_values)] <- column_mean_value
        som_codes[, variable_index] <- column_values
      }
    }
    
    # Ensure max.k is equal to or smaller than number of available codebook vectors (rows of som_codes)
    n_codes <- nrow(som_codes)
    if (!is.null(set.k) && set.k >= n_codes) {
      stop(sprintf(
        "Aborted SOM clustering: set.k = %d exceeds available codebook rows of %d - reduce set.k to ≤ %d",
        set.k, n_codes, n_codes - 1
      ))
    }
    if (max.k >= n_codes) {
      stop(sprintf(
        "Aborted SOM clustering: max.k = %d exceeds available codebook rows of %d - reduce max.k to ≤ %d",
        max.k, n_codes, n_codes - 1
      ))
    }
    
    # Create function to perform kmeans clustering and calculate within‐cluster sum of squares (wss) for each cluster (sum of squared Euclidean distances of SOM units to cluster center)
    calculate.wss <- function(som_codes, max.k) {
      wss <- numeric(max.k) #wss vector for k = 1 ... max.k
      fits <- vector("list", max.k) #store kmeans fits
      wss[1] <- (nrow(som_codes) - 1) * sum(apply(som_codes, 2, stats::var)) #calculate wss for k = 1 (= total sum of squared distances to overall mean)
      if(max.k >= 2){
        wss[2:max.k] <- sapply(2:max.k, function(i) { #calculate wss for k = 2 ... max.k via kmeans
          km <- stats::kmeans(som_codes, centers = i, nstart = 30, iter.max = 1e5)
          fits[[i]] <<- km
          sum(km$withinss)
        }) 
      }
      list(wss = wss, fits = fits)
    }
    
    # Create function to calculate BIC
    calculate.wssBIC <- function (wss, som_codes) {
      if (any(is.na(wss))) stop("wss contains NA - cannot calculate BIC")
      if (any(wss < 0)) stop("wss contains negative values - cannot calculate BIC")
      N <- nrow(som_codes) #number of output cells
      k_vals <- seq_along(wss)
      eps <- .Machine$double.eps #avoid log(0)
      BIC_vec <- N * log((wss + eps) / N) + log(N) * k_vals #BIC for each k
      if (any(is.na(BIC_vec) | is.infinite(BIC_vec))) stop("BIC calculation produced NA or Inf")
      BIC_vec
    }
    
    # Create function to determine optimal number of clusters by selecting smallest k where BIC improvement falls below BIC threshold, otherwise, choose k with minimum BIC
    select.k.BICthresh <- function(BIC_vec, BIC.thresh, set.k = NULL) {
      if (!is.null(set.k)) return(set.k) #if user-specified number of clusters is given, return it
      if (length(BIC_vec) < 2) return(1) #if there is only one or no BIC value, set default to k = 1
      if (all(is.na(BIC_vec))) stop("All BIC values are NA - cannot determine optimal number of clusters")
      if (is.na(BIC_vec[1]) || is.na(BIC_vec[2])) stop("The first or second BIC value is NA - cannot evaluate k = 1 vs k = 2")
      som_N_clusters <- NA #store optimal k
      for (k in 2:length(BIC_vec)) {
        if (!is.na(BIC_vec[k]) && !is.na(BIC_vec[k - 1]) && ((BIC_vec[k - 1] - BIC_vec[k]) < BIC.thresh)) { #if improvement is not at least as large as threshold, select previous k
          som_N_clusters <- k - 1
          break
        }
      }
      if (is.na(som_N_clusters)) {
        som_N_clusters <- which.min(replace(BIC_vec, is.na(BIC_vec), Inf)) #if all differences exceed threshold, pick k with lowest BIC
      }
      som_N_clusters
    }
    
    # Create function to determine optimal number of clusters using BIC elbow rule
    select.k.BICelbow <- function(BIC_vec, BIC.thresh, set.k = NULL) {
      if (!is.null(set.k)) return(set.k) #user-specified k
      if (length(BIC_vec) < 2) return(1) # if there is only one or no BIC value, return k = 1
      if (all(is.na(BIC_vec))) stop("All BIC values are NA - cannot determine optimal number of clusters")
      if (is.na(BIC_vec[1]) || is.na(BIC_vec[2])) stop("The first or second BIC value is NA - cannot evaluate k = 1 vs k = 2")
      if (((BIC_vec[1] - BIC_vec[2]) < BIC.thresh)) { #computes ΔBIC between k = 1 and k = 2 and if that drop is smaller than BIC threshold, there is no real improvement from adding second cluster, so pick k = 1
        som_N_clusters <- 1 #k = 1
      } else {
        if (length(BIC_vec) <= 2) { #if only 2 BIC values, pick k = 2
          som_N_clusters <- 2
        } else {
          delta_BIC_vec <- BIC_vec[-length(BIC_vec)] - BIC_vec[-1] #calculate BIC drop (improvement) from k-1 -> k
          valid_delta <- which(is.finite(delta_BIC_vec) & !is.na(delta_BIC_vec)) #keep only finite ΔBIC
          if (length(valid_delta) < 2) { #not enough info for elbow clustering
            som_N_clusters <- which.min(replace(BIC_vec, !is.finite(BIC_vec) | is.na(BIC_vec), Inf)) #fallback to best available BIC
            return(som_N_clusters)
          }
          delta_BIC_vec_valid <- delta_BIC_vec[valid_delta] #subset ΔBIC
          deltaBIC_clusters_valid <- stats::cutree(stats::hclust(stats::dist(delta_BIC_vec_valid), method = "ward.D2"), k = 2) #cluster valid ΔBIC values into two groups using hierarchical clustering
          deltaBIC_clusters <- rep(NA_integer_, length(delta_BIC_vec)) #expand back to full length
          deltaBIC_clusters[valid_delta] <- deltaBIC_clusters_valid #fill valid positions
          valid_groups <- which(!is.na(deltaBIC_clusters) & is.finite(delta_BIC_vec) & !is.na(delta_BIC_vec)) #valid grouped ΔBIC entries
          best_group <- unname(which.max(tapply(delta_BIC_vec[valid_groups], deltaBIC_clusters[valid_groups], mean, na.rm = TRUE))) #identify group with largest mean BIC drop ("real improvements", large drops)
          best_group_indices <- which(deltaBIC_clusters == best_group & is.finite(delta_BIC_vec) & !is.na(delta_BIC_vec)) #extract indices in ΔBIC vector that are "large drops"
          if (length(best_group_indices) > 0) { #if "large drop" group exists, pick k with lowest BIC within that group
            k_best_group <- best_group_indices + 1 #ΔBIC index i corresponds to drop into k = i + 1
            som_N_clusters <- k_best_group[which.min(replace(BIC_vec[k_best_group], is.na(BIC_vec[k_best_group]), Inf))] #pick k corresponding to lowest BIC in best group
          } else {
            som_N_clusters <- which.min(replace(BIC_vec, is.na(BIC_vec), Inf)) #pick k corresponding to global minimum BIC if best group is empty
            message("Warning: No large-drop group (i.e., the elbow) was detected - picking k corresponding to global minimum BIC")
          }
        }
      }
      som_N_clusters
    }
    
    # Create function to determine optimal number of clusters using DB elbow rule (for k >= 2 only)
    select.k.DBelbow <- function(DB_vec, set.k = NULL) {
      if (!is.null(set.k)) return(set.k) #return user-specified k
      if (length(DB_vec) < 3) return(2) #return k = 2 if fewer than k = 2..3 exist
      DB_vec_k2plus <- DB_vec[2:length(DB_vec)] #extract DB values for k >= 2
      if (all(is.na(DB_vec_k2plus))) stop("All DB values for k >= 2 are NA - cannot determine optimal number of clusters")
      if (length(DB_vec_k2plus) <= 1) return(2) #return k = 2 if only one k >= 2 value exists
      if (length(DB_vec_k2plus) <= 2) return(which.min(replace(DB_vec_k2plus, is.na(DB_vec_k2plus), Inf)) + 1) #select k with minimum DB if only k = 2..3 exist
      delta_DB_vec <- (DB_vec_k2plus[-length(DB_vec_k2plus)] - DB_vec_k2plus[-1]) / abs(DB_vec_k2plus[-length(DB_vec_k2plus)]) #calculate relative DB improvement from k-1 -> k for k >= 3
      deltaDB_clusters <- stats::cutree(stats::hclust(stats::dist(delta_DB_vec), method = "ward.D2"), k = 2) #cluster relative ΔDB values into two groups using hierarchical clustering
      best_group <- unname(which.max(tapply(delta_DB_vec, deltaDB_clusters, mean, na.rm = TRUE))) #identify group with largest mean relative ΔDB ("real improvements", large drops)
      best_group_indices <- which(deltaDB_clusters == best_group) #extract indices in relative ΔDB vector that are "large drops"
      if (length(best_group_indices) > 0) { #select k within best group if group exists
        k_best_group <- best_group_indices + 1 #map ΔDB index i to k index within DB_vec_k2plus
        som_N_clusters <- (k_best_group[which.min(replace(DB_vec_k2plus[k_best_group], is.na(DB_vec_k2plus[k_best_group]), Inf))] + 1) #convert back to original k by adding 1
      } else {
        som_N_clusters <- which.min(replace(DB_vec_k2plus, is.na(DB_vec_k2plus), Inf)) + 1 #select k with global minimum DB among k >= 2
        message("Warning: No large-drop group (i.e., the elbow) was detected - picking k corresponding to global minimum DB")
      }
      som_N_clusters
    }
    
    # Perform clustering of SOM codebook vectors based on specified method
    
    # Clustering method: kmeans + BICelbow
    if (clustering.method == "kmeans+BICelbow") {
      kmeans_results <- calculate.wss(som_codes, max.k) #perform kmeans and calculate within‐cluster sum of squares (wss)
      BIC_vec <- calculate.wssBIC(kmeans_results$wss, som_codes) #calculate BIC based on wss
      som_N_clusters <- select.k.BICelbow(BIC_vec, BIC.thresh, set.k) #determine optimal number of clusters using BIC elbow rule
      if (som_N_clusters == 1) { #if optimal k = 1
        som_cluster <- rep(1, nrow(som_codes)) #assign all units to a single cluster
      } else { 
        som_cluster <- kmeans_results$fits[[som_N_clusters]]$cluster #assign clusters from kmeans using selected k
      }
    }
    
    # Clustering method: kmeans + BICthresh
    if (clustering.method == "kmeans+BICthreshold") {
      kmeans_results <- calculate.wss(som_codes, max.k) #perform kmeans and calculate within‐cluster sum of squares (wss)
      BIC_vec <- calculate.wssBIC(kmeans_results$wss, som_codes) #calculate BIC based on wss
      som_N_clusters <- select.k.BICthresh(BIC_vec, BIC.thresh, set.k) #determine optimal number of clusters using BIC threshold rule (selecting smallest k where BIC improvement falls below BIC threshold, otherwise, choose k with minimum BIC)
      if (som_N_clusters == 1) { #if optimal k = 1
        som_cluster <- rep(1L, nrow(som_codes)) #assign all units to a single cluster
      } else { 
        som_cluster <- kmeans_results$fits[[som_N_clusters]]$cluster #assign clusters from kmeans using selected k
      }
    }
    
    # Clustering method: GMM (Gaussian mixture model) + BIC threshold
    if (clustering.method == "GMM+BICthreshold") {
      mclustBIC <- getFromNamespace("mclustBIC", "mclust")
      modelNames_1 <- c("EII", "VII") #stable set 1
      modelNames_2 <- c("EEI", "VVI") #stable set 2
      modelNames_3 <- c("VEI", "EVI") #mid set 1
      modelNames_4 <- c("EEV", "VEE") #mid set 2
      modelNames_5 <- c("EVE", "VEV") #mid set 3
      modelNames_6 <- c("VVE", "EVV") #expand set 1
      modelNames_7 <- c("EEE", "VVV") #expand set 2
      modelNames_stable <- list(modelNames_1, modelNames_2) #stable tier (preferred)
      modelNames_mid <- list(modelNames_3, modelNames_4, modelNames_5) #mid tier (if stable fails)
      modelNames_expand <- list(modelNames_6, modelNames_7) #expand tier (if stable+mid fails)
      modelNames_all <- unique(c(modelNames_1, modelNames_2, modelNames_3, modelNames_4, modelNames_5, modelNames_6, modelNames_7)) #all candidate covariance models
      mclust_BIC_matrix <- matrix(NA_real_, nrow = max.k, ncol = length(modelNames_all)) #initialize matrix of BIC values (rows = G, cols = modelNames)
      rownames(mclust_BIC_matrix) <- as.character(seq_len(max.k)) #set k rownames
      colnames(mclust_BIC_matrix) <- modelNames_all #set covariance model names
      BIC_vec <- rep(NA_real_, max.k) #initialize vector to store best BIC value per k
      available_k_values <- seq_len(max.k) #available k values (G values)
      for (k_value in available_k_values) { #extract best BIC values per k
        found_finite_BIC_for_k <- FALSE #track whether any finite BIC was obtained for this k
        best_BIC_value_for_k <- NA_real_ #store best BIC for this k (maximized)
        best_model_name_for_k <- NA_character_ #store covariance model name for best BIC
        for (modelNames_set in modelNames_stable) { #try stable tier first
          mclust_BIC_set <- try(mclustBIC(som_codes, G = k_value, modelNames = modelNames_set, verbose = FALSE), silent = TRUE) #fit GMM BIC for this k and model set
          if (inherits(mclust_BIC_set, "try-error")) next
          mclust_BIC_set <- as.matrix(mclust_BIC_set) #ensure matrix
          if (nrow(mclust_BIC_set) > 1) {
            row_idx_set <- which(suppressWarnings(as.integer(rownames(mclust_BIC_set))) == k_value) #match row if multiple rows returned
            if (length(row_idx_set) == 1) mclust_BIC_set <- mclust_BIC_set[row_idx_set, , drop = FALSE]
          }
          mclust_BIC_set <- mclust_BIC_set[1, , drop = FALSE] #force 1-row matrix
          mclust_BIC_set[!is.finite(mclust_BIC_set) | is.na(mclust_BIC_set)] <- NA_real_ #replace non-finite BIC with NA
          for (mn in colnames(mclust_BIC_set)) { #store BIC values in the full matrix
            mclust_BIC_matrix[as.character(k_value), mn] <- as.numeric(mclust_BIC_set[1, mn])
          }
          row_vals <- as.numeric(mclust_BIC_set[1, ]) #extract BIC values for this k across covariance models
          names(row_vals) <- colnames(mclust_BIC_set)
          row_vals <- row_vals[is.finite(row_vals) & !is.na(row_vals)] #keep only finite values
          if (length(row_vals) == 0) next
          found_finite_BIC_for_k <- TRUE
          if (is.na(best_BIC_value_for_k) || max(row_vals) > best_BIC_value_for_k) {
            best_BIC_value_for_k <- max(row_vals) #store highest (best) BIC across covariance models for this k (mclust BIC is maximized)
            best_model_name_for_k <- names(which.max(row_vals)) #store covariance model name for best BIC
          }
        }
        if (!found_finite_BIC_for_k) { #try mid tier only if stable tier failed completely
          for (modelNames_set in modelNames_mid) { #try mid tier
            mclust_BIC_set <- try(mclustBIC(som_codes, G = k_value, modelNames = modelNames_set), silent = TRUE) #fit GMM BIC for this k and model set
            if (inherits(mclust_BIC_set, "try-error")) next
            mclust_BIC_set <- as.matrix(mclust_BIC_set) #ensure matrix
            if (nrow(mclust_BIC_set) > 1) {
              row_idx_set <- which(suppressWarnings(as.integer(rownames(mclust_BIC_set))) == k_value) #match row if multiple rows returned
              if (length(row_idx_set) == 1) mclust_BIC_set <- mclust_BIC_set[row_idx_set, , drop = FALSE]
            }
            mclust_BIC_set <- mclust_BIC_set[1, , drop = FALSE] #force 1-row matrix
            mclust_BIC_set[!is.finite(mclust_BIC_set) | is.na(mclust_BIC_set)] <- NA_real_ #replace non-finite BIC with NA
            for (mn in colnames(mclust_BIC_set)) { #store BIC values in the full matrix
              mclust_BIC_matrix[as.character(k_value), mn] <- as.numeric(mclust_BIC_set[1, mn])
            }
            row_vals <- as.numeric(mclust_BIC_set[1, ]) #extract BIC values for this k across covariance models
            names(row_vals) <- colnames(mclust_BIC_set)
            row_vals <- row_vals[is.finite(row_vals) & !is.na(row_vals)] #keep only finite values
            if (length(row_vals) == 0) next
            found_finite_BIC_for_k <- TRUE
            if (is.na(best_BIC_value_for_k) || max(row_vals) > best_BIC_value_for_k) {
              best_BIC_value_for_k <- max(row_vals) #store highest (best) BIC across covariance models for this k (mclust BIC is maximized)
              best_model_name_for_k <- names(which.max(row_vals)) #store covariance model name for best BIC
            }
          }
        }
        if (!found_finite_BIC_for_k) { #try expand tier only if stable+mid tiers failed completely
          for (modelNames_set in modelNames_expand) { #try expand tier
            mclust_BIC_set <- try(mclustBIC(som_codes, G = k_value, modelNames = modelNames_set), silent = TRUE) #fit GMM BIC for this k and model set
            if (inherits(mclust_BIC_set, "try-error")) next
            mclust_BIC_set <- as.matrix(mclust_BIC_set) #ensure matrix
            if (nrow(mclust_BIC_set) > 1) {
              row_idx_set <- which(suppressWarnings(as.integer(rownames(mclust_BIC_set))) == k_value) #match row if multiple rows returned
              if (length(row_idx_set) == 1) mclust_BIC_set <- mclust_BIC_set[row_idx_set, , drop = FALSE]
            }
            mclust_BIC_set <- mclust_BIC_set[1, , drop = FALSE] #force 1-row matrix
            mclust_BIC_set[!is.finite(mclust_BIC_set) | is.na(mclust_BIC_set)] <- NA_real_ #replace non-finite BIC with NA
            for (mn in colnames(mclust_BIC_set)) { #store BIC values in the full matrix
              mclust_BIC_matrix[as.character(k_value), mn] <- as.numeric(mclust_BIC_set[1, mn])
            }
            row_vals <- as.numeric(mclust_BIC_set[1, ]) #extract BIC values for this k across covariance models
            names(row_vals) <- colnames(mclust_BIC_set)
            row_vals <- row_vals[is.finite(row_vals) & !is.na(row_vals)] #keep only finite values
            if (length(row_vals) == 0) next
            found_finite_BIC_for_k <- TRUE
            if (is.na(best_BIC_value_for_k) || max(row_vals) > best_BIC_value_for_k) {
              best_BIC_value_for_k <- max(row_vals) #store highest (best) BIC across covariance models for this k (mclust BIC is maximized)
              best_model_name_for_k <- names(which.max(row_vals)) #store covariance model name for best BIC
            }
          }
        }
        if (!is.na(best_BIC_value_for_k) && is.finite(best_BIC_value_for_k)) {
          BIC_vec[k_value] <- best_BIC_value_for_k #store best BIC for this k
        }
      }
      if (all(is.na(BIC_vec))) stop("All GMM BIC values are NA - cannot determine optimal number of clusters")
      BIC_vec <- -BIC_vec #invert sign so threshold selector (which assumes lower = better) works with mclust BIC (higher = better)
      som_N_clusters <- select.k.BICthresh(BIC_vec, BIC.thresh, set.k) #select k using BIC threshold rule (first k where improvement falls below threshold, otherwise global optimum)
      row_match <- which(available_k_values == som_N_clusters) #match selected k to BIC matrix row
      if (length(row_match) != 1) stop("Selected GMM k not found in mclust BIC matrix - check mclust output")
      best_model_name <- colnames(mclust_BIC_matrix)[which.max(replace(mclust_BIC_matrix[row_match, ], !is.finite(mclust_BIC_matrix[row_match, ]) | is.na(mclust_BIC_matrix[row_match, ]), -Inf))] #identify covariance model with highest BIC at selected k
      if (som_N_clusters == 1) { #if optimal k = 1
        som_cluster <- rep(1L, nrow(som_codes)) #assign all units to a single cluster
      } else {
        som_cluster <- mclust::Mclust(som_codes, G = som_N_clusters, verbose = FALSE, modelNames = best_model_name)$classification #refit GMM at selected k and best covariance model to obtain cluster assignments
      }
    }
    
    # Clustering method: hierarchical clustering + Davies-Bouldin pruning (allow k = 1 via permutation test on DB for k = 2)
    if (clustering.method == "hierarchical+DB") {
      davies_bouldin_values <- rep(NA_real_, max.k) #store DB values across k
      dist_som_codes <- stats::dist(som_codes) #compute pairwise Euclidean distances
      hierarchical_clust_som_codes <- stats::hclust(dist_som_codes, method = "ward.D2") #perform hierarchical clustering
      if (!is.null(set.k)) { #check if user specified k
        som_N_clusters <- as.integer(set.k) #set number of clusters to user-specified k
        if (som_N_clusters <= 1) { #check if user specified k <= 1
          som_N_clusters <- 1L #set number of clusters to 1
          som_cluster <- rep(1L, nrow(som_codes)) #assign all units to single cluster
        } else {
          som_cluster <- stats::cutree(hierarchical_clust_som_codes, som_N_clusters) #extract cluster assignments for user-specified k
        }
        BIC_vec <- rep(NA_real_, max.k) #initialize BIC vector as NA
      } else {
        cluster_assignments_k2 <- stats::cutree(hierarchical_clust_som_codes, 2) #extract cluster assignments for k = 2
        davies_bouldin_observed <- tryCatch({clusterCrit::intCriteria(as.matrix(som_codes), cluster_assignments_k2, "Davies_Bouldin")$davies_bouldin}, error = function(e) NA_real_) #calculate observed Davies-Bouldin index for k = 2
        if (is.na(davies_bouldin_observed) || !is.finite(davies_bouldin_observed)) { #check if observed DB is invalid
          som_N_clusters <- 1L #set number of clusters to 1
          som_cluster <- rep(1L, nrow(som_codes)) #assign all units to single cluster
          BIC_vec <- rep(NA_real_, max.k) #initialize BIC vector as NA
        } else {
          n_perm <- 500 #set number of permutations
          davies_bouldin_null <- numeric(n_perm) #initialize null DB vector
          for (i in seq_len(n_perm)) { #iterate over permutations
            permuted_som_codes <- apply(som_codes, 2, sample) #permute values within each column
            permuted_dist_som_codes <- stats::dist(permuted_som_codes) #compute distances for permuted data
            permuted_hierarchical_clust <- stats::hclust(permuted_dist_som_codes, method = "ward.D2") #perform hierarchical clustering on permuted data
            permuted_cluster_assignments_k2 <- stats::cutree(permuted_hierarchical_clust, 2) #extract permuted cluster assignments for k = 2
            davies_bouldin_null[i] <- tryCatch({clusterCrit::intCriteria(as.matrix(permuted_som_codes), permuted_cluster_assignments_k2, "Davies_Bouldin")$davies_bouldin}, error = function(e) NA_real_) #calculate null Davies-Bouldin index for k = 2
          }
          davies_bouldin_null <- davies_bouldin_null[!is.na(davies_bouldin_null) & is.finite(davies_bouldin_null)] #remove invalid null DB values
          pval <- if (length(davies_bouldin_null) == 0) 1 else mean(davies_bouldin_null <= davies_bouldin_observed) #calculate permutation p-value
          if (pval >= 0.05) { #check if observed DB is not significantly smaller than null
            som_N_clusters <- 1L #set number of clusters to 1
            som_cluster <- rep(1L, nrow(som_codes)) #assign all units to single cluster
            BIC_vec <- rep(NA_real_, max.k) #initialize BIC vector as NA
          } else {
            davies_bouldin_values <- rep(NA_real_, max.k) #initialize Davies-Bouldin vector for all k
            davies_bouldin_values[1] <- Inf #assign infinite DB to k = 1
            davies_bouldin_values[2] <- davies_bouldin_observed #store observed DB for k = 2
            if (max.k >= 3) { #check if additional k values exist
              for (k in 3:max.k) { #iterate over candidate k values
                cluster_assignments_k <- stats::cutree(hierarchical_clust_som_codes, k) #extract cluster assignments for k
                davies_bouldin_values[k] <- tryCatch({clusterCrit::intCriteria(as.matrix(som_codes), cluster_assignments_k, "Davies_Bouldin")$davies_bouldin}, error = function(e) NA_real_) #calculate Davies-Bouldin index for k
              }
            }
            if (all(is.na(davies_bouldin_values[2:max.k]) | is.infinite(davies_bouldin_values[2:max.k]))) { #check if all Davies-Bouldin values failed
              som_N_clusters <- 1L #set number of clusters to 1
              som_cluster <- rep(1L, nrow(som_codes)) #assign all units to single cluster
            } else {
              som_N_clusters <- which.min(replace(davies_bouldin_values[2:max.k], is.na(davies_bouldin_values[2:max.k]), Inf)) + 1L #select k with minimum Davies-Bouldin index
              som_cluster <- stats::cutree(hierarchical_clust_som_codes, som_N_clusters) #extract final cluster assignments
            }
            BIC_vec <- rep(NA_real_, max.k) #initialize BIC vector as NA
          }
        }
      }
    }  
    
    # Clustering method: HDBSCAN
    if (clustering.method == "HDBSCAN") {
      dist_som_codes <- stats::dist(som_codes) #compute pairwise Euclidean distances
      max_minPts <- max(2L, floor(nrow(som_codes) / 1.5)) #set max minPts
      minPts_min <- max(3L, floor(0.05 * nrow(som_codes))) #at least 3, or 5% of units
      minPts_min <- min(minPts_min, max_minPts) #ensure valid range
      minPts_vals <- seq(minPts_min, max_minPts) #grid of minPts values
      hdbscan_model_results <- data.frame(minPts = integer(), n_clusters = integer(), mean_mem = numeric()) #initialize result table
      for (m in minPts_vals) { #for each minPts value ...
        hdbscan_model <- dbscan::hdbscan(som_codes, minPts = m) #run HDBSCAN
        hdbscan_assignment <- hdbscan_model$cluster #extract clusters
        hdbscan_N_clusters <- length(unique(hdbscan_assignment[hdbscan_assignment != 0])) #count real clusters (exclude 0 = noise)
        hdbscan_membership_prob <- mean(hdbscan_model$membership_prob[hdbscan_assignment != 0], na.rm = TRUE) #record mean membership prob
        hdbscan_model_results <- rbind(hdbscan_model_results, #append results to result dataframe
                                       data.frame(minPts = m,
                                                  stringsAsFactors = FALSE,
                                                  n_clusters = hdbscan_N_clusters,
                                                  mean_mem = hdbscan_membership_prob))
      }
      valid_HDBSCAN <- which(hdbscan_model_results$n_clusters > 0 & !is.na(hdbscan_model_results$mean_mem)) #filter valid results
      if (length(valid_HDBSCAN) == 0) { #if no valid runs
        som_cluster <- rep(1L, nrow(som_codes)) #assign all points to one cluster
        som_N_clusters <- 1L
        BIC_vec <- rep(NA_real_, max.k)
      } else { #run best model and evaluate reassignment strategies
        hdbscan_model_results <- hdbscan_model_results[valid_HDBSCAN, , drop = FALSE] #subset to valid runs
        if (!is.null(set.k)) { #check if user specified k
          som_N_clusters <- as.integer(set.k) #set number of clusters to user-specified k
          if (som_N_clusters <= 1) { #check if user specified k <= 1
            som_N_clusters <- 1L #set number of clusters to 1
            som_cluster <- rep(1L, nrow(som_codes)) #assign all points to one cluster
            BIC_vec <- rep(NA_real_, max.k)
          } else {
            valid_HDBSCAN_k <- which(hdbscan_model_results$n_clusters == som_N_clusters & !is.na(hdbscan_model_results$mean_mem)) #filter valid runs matching set.k
            if (length(valid_HDBSCAN_k) == 0) { #if no valid runs match requested k
              som_cluster <- rep(1L, nrow(som_codes)) #assign all points to one cluster
              som_N_clusters <- 1L
              BIC_vec <- rep(NA_real_, max.k)
            } else {
              best_minPts_row <- valid_HDBSCAN_k[which.max(hdbscan_model_results$mean_mem[valid_HDBSCAN_k])] #select best minPts among matching-k runs
              best_minPts <- hdbscan_model_results$minPts[best_minPts_row]
              best_hdbscan_model <- dbscan::hdbscan(som_codes, minPts = best_minPts) #run HDBSCAN with selected minPts
              base_clusters <- best_hdbscan_model$cluster #extract clusters (0 = noise)
              if (all(base_clusters == 0)) { #all noise
                som_cluster <- rep(1L, nrow(som_codes)) #assign all points to one cluster
                som_N_clusters <- 1L
              } else {
                noise_indices <- which(base_clusters == 0) #noise points (cluster = 0)
                core_indices <- which(base_clusters != 0) #core points (real clusters)
                if (length(noise_indices) > 0 && length(core_indices) > 0) { #assign noise points to nearest core cluster
                  hdbscan_nn_result <- get.knnx.custom(reference_data = som_codes[core_indices, , drop = FALSE], query_data = som_codes[noise_indices, , drop = FALSE], k = 1)
                  nearest_clusters <- base_clusters[core_indices][hdbscan_nn_result$nn.index]
                  base_clusters[noise_indices] <- nearest_clusters
                }
                som_cluster <- as.integer(factor(base_clusters)) #relabel clusters sequentially
                som_N_clusters <- length(unique(som_cluster)) #extract number of clusters
              }
              BIC_vec <- rep(NA_real_, max.k) #fill BIC with NA (not used for HDBSCAN)
            }
          }
        } else {
          best_minPts_row <- which.max(hdbscan_model_results$mean_mem) #select best minPts by mean membership probability
          best_minPts <- hdbscan_model_results$minPts[best_minPts_row]
          best_hdbscan_model <- dbscan::hdbscan(som_codes, minPts = best_minPts) #run HDBSCAN with selected minPts
          base_clusters <- best_hdbscan_model$cluster #extract clusters (0 = noise)
          if (all(base_clusters == 0)) { #all noise
            som_cluster <- rep(1L, nrow(som_codes)) #assign all points to one cluster
            som_N_clusters <- 1L
            BIC_vec <- rep(NA_real_, max.k)
          } else {
            noise_indices <- which(base_clusters == 0) #noise points (cluster = 0)
            core_indices <- which(base_clusters != 0) #core points (real clusters)
            if (length(noise_indices) == 0 || length(core_indices) == 0) { #no reassignment needed
              som_cluster <- as.integer(factor(base_clusters)) #relabel clusters sequentially
              som_N_clusters <- length(unique(som_cluster)) #extract number of clusters
            } else {
              hdbscan_nn_result <- get.knnx.custom(reference_data = som_codes[core_indices, , drop = FALSE], query_data = som_codes[noise_indices, , drop = FALSE], k = 1)
              nearest_clusters <- base_clusters[core_indices][hdbscan_nn_result$nn.index]
              clusters_nearest <- base_clusters
              clusters_nearest[noise_indices] <- nearest_clusters
              clusters_nearest <- as.integer(factor(clusters_nearest)) #relabel clusters sequentially
              k_nearest <- length(unique(clusters_nearest)) #extract number of clusters
              silhouette_nearest <- NA_real_
              if (k_nearest >= 2) {
                silhouette_nearest <- tryCatch({mean(cluster::silhouette(clusters_nearest, dist_som_codes)[, 3])}, error = function(e) NA_real_) #calculate silhouette
              }
              dist_thresh <- stats::quantile(hdbscan_nn_result$nn.dist, 0.95) #distance threshold for singleton assignment
              clusters_singleton <- base_clusters
              next_cluster <- max(base_clusters)
              for (i in seq_along(noise_indices)) {
                ni <- noise_indices[i]
                dist_to_core <- hdbscan_nn_result$nn.dist[i]
                if (dist_to_core > dist_thresh) {
                  next_cluster <- next_cluster + 1
                  clusters_singleton[ni] <- next_cluster
                } else {
                  clusters_singleton[ni] <- base_clusters[core_indices][hdbscan_nn_result$nn.index[i]]
                }
              }
              clusters_singleton <- as.integer(factor(clusters_singleton)) #relabel clusters sequentially
              k_singleton <- length(unique(clusters_singleton)) #extract number of clusters
              silhouette_singleton <- NA_real_
              if (k_singleton >= 2) {
                silhouette_singleton <- tryCatch({mean(cluster::silhouette(clusters_singleton, dist_som_codes)[, 3])}, error = function(e) NA_real_) #calculate silhouette
              }
              if (is.na(silhouette_singleton) && is.na(silhouette_nearest)) { #compare both strategies and select best based on silhouette score (or cluster count fallback)
                if (k_singleton > k_nearest) {
                  som_cluster <- clusters_singleton
                } else {
                  som_cluster <- clusters_nearest
                }
              } else if (is.na(silhouette_singleton) || silhouette_nearest >= silhouette_singleton) {
                som_cluster <- clusters_nearest
              } else {
                som_cluster <- clusters_singleton
              }
              som_N_clusters <- length(unique(som_cluster)) #extract number of clusters
            }
            BIC_vec <- rep(NA_real_, max.k) #fill BIC with NA (not used for HDBSCAN)
          }
        }
      }
    }    
    
    # Clustering method: OPTICS + Silhouette score
    if (clustering.method == "OPTICS+Silhouette") {
      dist_som_codes <- stats::dist(som_codes) #compute pairwise Euclidean distances
      silhouette_values <- rep(NA_real_, max.k) #best silhouette support per k
      if (!is.null(set.k)) { #check if user specified k
        som_N_clusters <- as.integer(set.k) #set number of clusters to user-specified k
        if (som_N_clusters <= 1) { #check if user specified k <= 1
          som_N_clusters <- 1L #set number of clusters to 1
          som_cluster <- rep(1L, nrow(som_codes)) #assign all units to single cluster
          BIC_vec <- rep(NA_real_, max.k) #initialize BIC vector as NA
        } else {
          n_codes <- nrow(som_codes) #number of SOM codebook vectors
          minPts_start <- max(3L, floor(0.10 * n_codes)) #minimum minPts
          max_minPts <- min(15L, n_codes - 1L, floor(n_codes / 2)) #maximum minPts
          if (minPts_start > max_minPts) minPts_start <- max_minPts
          minPts_vals <- seq(minPts_start, max_minPts) #grid of minPts values
          xi_vals <- c(0.05, 0.10, 0.20, 0.30, 0.40) #Xi grid (include slightly larger values to make extraction less strict)
          optics_model_results <- data.frame(minPts = integer(),
                                             xi = numeric(),
                                             n_clusters = integer(),
                                             mean_silhouette = numeric(),
                                             result_index = integer(),
                                             stringsAsFactors = FALSE) #initialize result table
          optics_cluster_solutions <- list() #store cluster solutions
          result_counter <- 0L #result counter
          for (minPts_value in minPts_vals) { #iterate over minPts values
            optics_model <- try(dbscan::optics(som_codes, minPts = minPts_value, eps = NULL), silent = TRUE) #run OPTICS
            if (inherits(optics_model, "try-error")) next
            for (xi_value in xi_vals) { #iterate over Xi values
              xi_model <- try(withCallingHandlers(
                dbscan::extractXi(optics_model, xi = xi_value),
                warning = function(w) {
                  if (grepl("No clusters were found", conditionMessage(w))) invokeRestart("muffleWarning")
                }
              ), silent = TRUE)
              if (inherits(xi_model, "try-error") || is.null(xi_model$cluster)) next
              cluster_assignments <- as.integer(xi_model$cluster) #0 indicates noise/unassigned
              if (all(cluster_assignments == 0L)) next #no usable clustering for this parameter combination
              noise_indices <- which(cluster_assignments == 0L) #identify noise points
              core_indices <- which(cluster_assignments != 0L) #identify core points
              if (length(noise_indices) > 0 && length(core_indices) > 0) { #assign noise points to nearest core cluster
                optics_nn_result <- get.knnx.custom(reference_data = som_codes[core_indices, , drop = FALSE], query_data = som_codes[noise_indices, , drop = FALSE], k = 1)
                cluster_assignments[noise_indices] <- cluster_assignments[core_indices][optics_nn_result$nn.index]
              }
              cluster_assignments_relabelled <- as.integer(factor(cluster_assignments)) #relabel clusters sequentially
              number_of_clusters <- length(unique(cluster_assignments_relabelled)) #number of clusters
              if (number_of_clusters != som_N_clusters) next #enforce user-specified k
              silhouette_value <- tryCatch({mean(cluster::silhouette(cluster_assignments_relabelled, dist_som_codes)[, 3])}, error = function(e) NA_real_) #calculate silhouette
              if (!is.na(silhouette_value) && is.finite(silhouette_value) && number_of_clusters <= max.k) {
                if (is.na(silhouette_values[number_of_clusters]) || silhouette_value > silhouette_values[number_of_clusters]) {
                  silhouette_values[number_of_clusters] <- silhouette_value
                }
              }
              if (is.na(silhouette_value)) next
              result_counter <- result_counter + 1L #increment result counter
              optics_model_results <- rbind(optics_model_results,
                                            data.frame(minPts = minPts_value,
                                                       xi = xi_value,
                                                       stringsAsFactors = FALSE,
                                                       n_clusters = number_of_clusters,
                                                       mean_silhouette = silhouette_value,
                                                       result_index = result_counter))
              optics_cluster_solutions[[result_counter]] <- cluster_assignments_relabelled #store cluster solution
            }
          }
          if (nrow(optics_model_results) == 0) { #no parameter combination produced requested k
            som_N_clusters <- 1L #set number of clusters to 1
            som_cluster <- rep(1L, nrow(som_codes)) #assign all units to single cluster
          } else {
            best_row <- which.max(optics_model_results$mean_silhouette) #select best silhouette
            best_silhouette <- optics_model_results$mean_silhouette[best_row] #extract best silhouette
            best_cluster_solution <- optics_cluster_solutions[[optics_model_results$result_index[best_row]]] #extract best clustering
            if (is.na(best_silhouette) || best_silhouette <= 0.25) { #collapse weak structure to k = 1
              som_N_clusters <- 1L
              som_cluster <- rep(1L, nrow(som_codes))
            } else {
              som_cluster <- best_cluster_solution #assign best clustering
            }
          }
          BIC_vec <- rep(NA_real_, max.k) #initialize BIC vector as NA
        }
      } else {
        n_codes <- nrow(som_codes) #number of SOM codebook vectors
        minPts_start <- max(5L, floor(0.10 * n_codes)) #minimum minPts
        max_minPts <- max(2L, min(15L, n_codes - 1L, floor(n_codes / 2))) #maximum minPts
        minPts_vals <- seq(minPts_start, max_minPts) #grid of minPts values
        xi_vals <- c(0.05, 0.10, 0.20, 0.30, 0.40) #Xi grid
        optics_model_results <- data.frame(minPts = integer(),
                                           xi = numeric(),
                                           n_clusters = integer(),
                                           mean_silhouette = numeric(),
                                           result_index = integer(),
                                           stringsAsFactors = FALSE) #initialize result table
        optics_cluster_solutions <- list() #store cluster solutions
        result_counter <- 0L #result counter
        for (minPts_value in minPts_vals) { #iterate over minPts values
          optics_model <- try(dbscan::optics(som_codes, minPts = minPts_value, eps = NULL), silent = TRUE) #run OPTICS
          if (inherits(optics_model, "try-error")) next
          for (xi_value in xi_vals) { #iterate over Xi values
            xi_model <- try(withCallingHandlers(
              dbscan::extractXi(optics_model, xi = xi_value),
              warning = function(w) {
                if (grepl("No clusters were found", conditionMessage(w))) invokeRestart("muffleWarning")
              }
            ), silent = TRUE)
            if (inherits(xi_model, "try-error") || is.null(xi_model$cluster)) next
            cluster_assignments <- as.integer(xi_model$cluster) #0 indicates noise/unassigned
            if (all(cluster_assignments == 0L)) next #no usable clustering for this parameter combination
            noise_indices <- which(cluster_assignments == 0L) #identify noise points
            core_indices <- which(cluster_assignments != 0L) #identify core points
            if (length(noise_indices) > 0 && length(core_indices) > 0) { #assign noise points to nearest core cluster
              optics_nn_result <- get.knnx.custom(reference_data = som_codes[core_indices, , drop = FALSE], query_data = som_codes[noise_indices, , drop = FALSE], k = 1)
              cluster_assignments[noise_indices] <- cluster_assignments[core_indices][optics_nn_result$nn.index]
            }
            cluster_assignments_relabelled <- as.integer(factor(cluster_assignments)) #relabel clusters sequentially
            number_of_clusters <- length(unique(cluster_assignments_relabelled)) #number of clusters
            if (number_of_clusters < 2) next #require at least two clusters
            silhouette_value <- tryCatch({mean(cluster::silhouette(cluster_assignments_relabelled, dist_som_codes)[, 3])}, error = function(e) NA_real_) #calculate silhouette
            if (is.na(silhouette_value)) next
            result_counter <- result_counter + 1L #increment result counter
            optics_model_results <- rbind(optics_model_results,
                                          data.frame(minPts = minPts_value,
                                                     xi = xi_value,
                                                     stringsAsFactors = FALSE,
                                                     n_clusters = number_of_clusters,
                                                     mean_silhouette = silhouette_value,
                                                     result_index = result_counter))
            optics_cluster_solutions[[result_counter]] <- cluster_assignments_relabelled #store cluster solution
          }
        }
        if (nrow(optics_model_results) == 0) { #no valid clustering found
          som_cluster <- rep(1L, nrow(som_codes))
          som_N_clusters <- 1L
        } else {
          best_row <- which.max(optics_model_results$mean_silhouette) #select best silhouette
          best_silhouette <- optics_model_results$mean_silhouette[best_row] #extract best silhouette
          best_cluster_solution <- optics_cluster_solutions[[optics_model_results$result_index[best_row]]] #extract best clustering
          best_number_of_clusters <- optics_model_results$n_clusters[best_row] #extract best k
          if (is.na(best_silhouette) || best_silhouette <= 0.15) { #collapse weak structure to k = 1
            som_cluster <- rep(1L, nrow(som_codes))
            som_N_clusters <- 1L
          } else {
            som_cluster <- best_cluster_solution
            som_N_clusters <- best_number_of_clusters
          }
        }
        BIC_vec <- rep(NA_real_, max.k) #initialize BIC vector as NA
      }
    }
    
    # Extract SOM training samples
    som_training_samples <- rownames(som_model$data[[1]])
    
    # Create cluster_gridcell_assignments
    cluster_gridcell_assignments <- data.frame(Cluster = som_cluster[som_model$unit.classif],
                                               Gridcell = som_model$unit.classif,
                                               row.names = som_training_samples,
                                               stringsAsFactors = FALSE
    )
    
    # Create cluster_assignment
    cluster_assignment <- matrix(som_cluster[som_model$unit.classif],
                                 ncol = 1,
                                 dimnames = list(som_training_samples, NULL)
    )
    
    # Create replicate-specific soft ancestry matrix
    sample_to_unit_distance_matrix <- compute.sample.to.unit.distance.matrix.SOM(
      som_model = som_model,
      layer.distance.functions = SOM.output$layer.distance.functions,
      layer.weights = som_model$distance.weights
    )
    replicate_ancestry_matrix <- compute.replicate.ancestry.matrix.SOM(
      sample_to_unit_distance_matrix = sample_to_unit_distance_matrix,
      unit_cluster_labels = som_cluster
    )
    
    # Store generic support values for plotting
    if (exists("BIC_vec", inherits = FALSE) && any(is.finite(BIC_vec))) {
      support_vec <- BIC_vec
      support_label <- "BIC"
      support_higher_is_better <- FALSE
    }
    if (exists("davies_bouldin_values", inherits = FALSE) && any(is.finite(davies_bouldin_values))) {
      support_vec <- davies_bouldin_values
      support_label <- "Davies-Bouldin index"
      support_higher_is_better <- FALSE
    }
    if (exists("silhouette_values", inherits = FALSE) && any(is.finite(silhouette_values))) {
      support_vec <- silhouette_values
      support_label <- "Mean silhouette"
      support_higher_is_better <- TRUE
    }
    
    # Return results
    return(list(cluster_assignment = cluster_assignment, 
                BIC_vec = BIC_vec,
                support_vec = support_vec,
                support_label = support_label,
                support_higher_is_better = support_higher_is_better,
                som_model = som_model,
                som_cluster = som_cluster,
                som_N_clusters = som_N_clusters,
                cluster_gridcell_assignments = cluster_gridcell_assignments,
                replicate_ancestry_matrix = replicate_ancestry_matrix
    ))
  }
  
  # Collect results for all replicates
  results <- lapply(seq_len(N.replicates), replicate_clust)
  if (is.null(results)) return(invisible(NULL))
  
  # Combine results from all replicates
  cluster_assignment <- do.call(cbind, lapply(results, `[[`, "cluster_assignment"))
  sample_names <- rownames(results[[1]]$cluster_assignment)
  rownames(cluster_assignment) <- sample_names
  colnames(cluster_assignment) <- paste0("R", seq_len(ncol(cluster_assignment)))
  
  BIC_values <- do.call(cbind, lapply(results, `[[`, "BIC_vec"))
  rownames(BIC_values) <- paste0("k", seq_len(max.k))
  colnames(BIC_values) <- paste0("R", seq_len(ncol(BIC_values)))
  
  support_values <- do.call(cbind, lapply(results, `[[`, "support_vec"))
  rownames(support_values) <- paste0("k", seq_len(max.k))
  colnames(support_values) <- paste0("R", seq_len(ncol(support_values)))
  
  support_labels <- unique(unlist(lapply(results, function(x) x$support_label)))
  support_labels <- support_labels[!is.na(support_labels) & nzchar(support_labels)]
  support_label <- if (length(support_labels) == 0) NULL else support_labels[1]
  
  support_directions <- unique(unlist(lapply(results, function(x) x$support_higher_is_better)))
  support_directions <- support_directions[!is.na(support_directions)]
  support_higher_is_better <- if (length(support_directions) == 0) NA else support_directions[1]
  
  optim_k_vals <- sapply(results, `[[`, "som_N_clusters")
  optim_k_vals <- t(as.matrix(optim_k_vals))
  rownames(optim_k_vals) <- "optim_k_vals"
  colnames(optim_k_vals) <- paste0("R", seq_len(N.replicates))
  optim_k_mean <- mean(optim_k_vals, na.rm = T)
  if (all(is.na(optim_k_vals))) {
    stop("Aborted SOM clustering: all optimal K values are NA - check input data") 
  }
  optim_k_vector <- as.numeric(optim_k_vals)
  optim_k_vector <- optim_k_vector[is.finite(optim_k_vector)]
  k_levels <- sort(unique(optim_k_vector))
  optim_k_vals_counts <- as.numeric(table(factor(optim_k_vector, levels = k_levels)))
  optim_k_vals_props <- optim_k_vals_counts / sum(optim_k_vals_counts)
  optim_k_summary <- cbind(Count = optim_k_vals_counts,
                           Proportion = round(optim_k_vals_props, 2)
  )
  rownames(optim_k_summary) <- paste0("k", k_levels)
  
  som_models <- lapply(results, `[[`, "som_model")
  som_clusters <- lapply(results, `[[`, "som_cluster")
  cluster_gridcell_assignments <- lapply(results, `[[`, "cluster_gridcell_assignments")
  replicate_ancestry_matrices <- lapply(results, `[[`, "replicate_ancestry_matrix")
  
  mean_assignment_margin <- vapply(replicate_ancestry_matrices,
                                   calculate.mean.assignment.margin.SOM,
                                   numeric(1))
  names(mean_assignment_margin) <- paste0("R", seq_along(mean_assignment_margin))
  
  mean_normalized_assignment_entropy <- vapply(replicate_ancestry_matrices,
                                               calculate.mean.normalized.assignment.entropy.SOM,
                                               numeric(1))
  names(mean_normalized_assignment_entropy) <- paste0("R", seq_along(mean_normalized_assignment_entropy))
  
  # Replace NA values with 0.5 for labeling
  processed_input_data <- results[[1]]$som_model$data[[1]]
  processed_input_data[is.na(processed_input_data)] <- 0.5
  
  # Preprocess input data and generate cluster labels
  base::set.seed(set.seed.N) #set seed for reproducibility
  cluster_labels <- do.call(cbind, lapply(1:max.k, function(number_of_clusters) stats::kmeans(processed_input_data, centers = number_of_clusters, nstart = 30, iter.max = 1e5)$cluster)) #generate reference cluster labels (k = 1..max.k) for Hungarian relabeling
  rownames(cluster_labels) <- rownames(processed_input_data) #set row names for cluster labels
  colnames(cluster_labels) <- paste("k", 1:max.k, sep = '') #set column names
  
  # Extract SOM cluster assignments and filter replicates based on maximum K
  all_k <- apply(cluster_assignment, 2, max, na.rm = TRUE) #calculate maximum K for each replicate (column)
  if (length(which(optim_k_vals <= max.k)) == 0) stop("Aborted SOM clustering: no replicates have k ≤ max.k - increase max.k or check input data")
  assignment_matrix <- cluster_assignment[, which(optim_k_vals <= max.k), drop = FALSE] #filter to keep replicates with K <= max_k
  
  # Relabel across replicates with Hungarian algorithm
  for (replicate_index in seq_len(ncol(assignment_matrix))) {
    replicate_k <- max(assignment_matrix[, replicate_index], na.rm = TRUE)
    if (is.na(replicate_k) || replicate_k < 2) next
    reference_sample_cluster_labels <- cluster_labels[, replicate_k] #reference labels for this K
    replicate_sample_cluster_labels <- assignment_matrix[, replicate_index] #labels from this replicate
    cluster_overlap_count_matrix <- table(factor(reference_sample_cluster_labels, levels = 1:replicate_k), base::factor(replicate_sample_cluster_labels, levels = 1:replicate_k)) #create table (rows: ref clusters, cols: replicate clusters; with levels 1:replicate_k)
    if (nrow(cluster_overlap_count_matrix) != ncol(cluster_overlap_count_matrix)) { #pad overlap table to square if needed for Hungarian algorithm
      new_size <- max(nrow(cluster_overlap_count_matrix), ncol(cluster_overlap_count_matrix)) #determine new size for square table
      padded_cluster_overlap_count_matrix <- matrix(0, nrow = new_size, ncol = new_size) #create square matrix of zeros
      padded_cluster_overlap_count_matrix[1:nrow(cluster_overlap_count_matrix), 1:ncol(cluster_overlap_count_matrix)] <- cluster_overlap_count_matrix #copy original table into upper-left
      cluster_overlap_count_matrix <- padded_cluster_overlap_count_matrix #update table to padded table
    }
    assignment_cost_matrix <- max(cluster_overlap_count_matrix) - cluster_overlap_count_matrix #build cost matrix (max overlap = minimal cost)
    reference_cluster_to_replicate_cluster_map <- as.integer(clue::solve_LSAP(assignment_cost_matrix)) #returns replicate cluster (col) chosen for each reference cluster (row)
    replicate_cluster_to_reference_cluster_map <- integer(replicate_k) #invert mapping so we can relabel replicate -> reference
    replicate_cluster_to_reference_cluster_map[reference_cluster_to_replicate_cluster_map] <- seq_len(replicate_k) #replicate_cluster_to_reference_cluster_map[replicate_cluster] = reference_cluster
    assignment_matrix[, replicate_index] <- as.integer(replicate_cluster_to_reference_cluster_map[replicate_sample_cluster_labels]) #assign new cluster labels to samples
  }
  
  # Build ancestry from the re-labelled matrix
  k.max <- max(assignment_matrix, na.rm = TRUE)
  prop_list <- lapply(seq_len(nrow(assignment_matrix)), function(i) {
    prop.table(table(factor(assignment_matrix[i, ], levels = seq_len(k.max))))
  })
  ancestry_matrix <- do.call(rbind, prop_list)
  colnames(ancestry_matrix) <- paste0("Cluster_", seq_len(ncol(ancestry_matrix)))
  rownames(ancestry_matrix) <- rownames(assignment_matrix)
  
  # Calculate median map variance (neuron-weighted) and median eta squared (cluster separation effect size) for each variable in each layer
  codebook_list_1 <- kohonen::getCodes(som_models[[1]]) #extract codebook list to get number of layers
  if (!is.list(codebook_list_1)) codebook_list_1 <- list(codebook_list_1) #ensure list
  n_layers <- length(codebook_list_1) #number of SOM layers
  calculate.map.variance.per.variable <- function(codebook_matrix, som_model, baseline_weight = 1) { #calculate variance across map 
    neuron_sample_counts <- tabulate(som_model$unit.classif, nbins = nrow(codebook_matrix)) #number of samples per neuron
    neuron_weights <- neuron_sample_counts + baseline_weight #neuron weight = baseline + sample support (empty neurons still count)
    apply(codebook_matrix, 2, function(variable_values) {
      valid_rows <- is.finite(variable_values) & !is.na(variable_values) & is.finite(neuron_weights) & !is.na(neuron_weights) #identify valid values
      variable_values <- variable_values[valid_rows] #subset variable values
      weights <- neuron_weights[valid_rows] #subset weights
      if (length(variable_values) < 2) return(NA_real_) #require at least 2 observations
      if (sum(weights) <= 0) return(NA_real_) #require positive total weight
      grand_mean <- sum(weights * variable_values) / sum(weights) #weighted mean
      sum(weights * (variable_values - grand_mean)^2) / sum(weights) #weighted variance
    })
  }
  calculate.etasquared.per.variable <- function(codebook_matrix, neuron_cluster_vector, som_model, baseline_weight = 1) { #create function to calculate eta squared
    if (length(neuron_cluster_vector) != nrow(codebook_matrix)) return(rep(NA_real_, ncol(codebook_matrix))) #return NA if mismatch
    neuron_cluster_vector <- as.integer(neuron_cluster_vector) #ensure cluster vector is integer
    valid_cluster_rows <- is.finite(neuron_cluster_vector) & !is.na(neuron_cluster_vector) #identify valid cluster rows
    n_units <- length(neuron_cluster_vector) #store original number of neurons (before filtering)
    neuron_sample_counts <- tabulate(som_model$unit.classif, nbins = n_units) #number of samples per neuron (all neurons)
    codebook_matrix <- codebook_matrix[valid_cluster_rows, , drop = FALSE] #subset codebook to valid rows
    neuron_cluster_vector <- neuron_cluster_vector[valid_cluster_rows] #subset cluster vector to valid rows
    neuron_sample_counts <- neuron_sample_counts[valid_cluster_rows] #subset to valid rows
    neuron_weights <- neuron_sample_counts + baseline_weight #neuron weight = baseline + sample support (empty neurons still count)
    apply(codebook_matrix, 2, function(variable_values) {
      valid_variable_rows <- is.finite(variable_values) & !is.na(variable_values) & is.finite(neuron_weights) & !is.na(neuron_weights) #identify valid values
      variable_values <- variable_values[valid_variable_rows] #subset variable values
      cluster_labels <- neuron_cluster_vector[valid_variable_rows] #subset cluster labels
      weights <- neuron_weights[valid_variable_rows] #subset weights
      if (length(variable_values) < 2) return(NA_real_) #require at least 2 observations
      if (length(unique(cluster_labels)) < 2) return(NA_real_) #require at least 2 clusters
      if (sum(weights) <= 0) return(NA_real_) #require positive total weight
      grand_mean <- sum(weights * variable_values) / sum(weights) #weighted grand mean
      total_sum_of_squares <- sum(weights * (variable_values - grand_mean)^2) #weighted total sum of squares
      if (!is.finite(total_sum_of_squares) || total_sum_of_squares <= 0) return(0) #handle degenerate cases
      cluster_means <- tapply(weights * variable_values, cluster_labels, sum) / tapply(weights, cluster_labels, sum) #weighted cluster means
      cluster_sizes <- tapply(weights, cluster_labels, sum) #cluster "size" = sum of (baseline + hits) across neurons
      between_cluster_sum_of_squares <- sum(cluster_sizes * (cluster_means - grand_mean)^2) #weighted between-cluster sum of squares
      as.numeric(between_cluster_sum_of_squares / total_sum_of_squares) #eta squared
    })
  }
  median_map_variance_variable_importance <- vector("list", n_layers) #store median map variance per layer
  median_etasquared_variable_importance <- vector("list", n_layers) #store median eta^2 per layer
  n_clusters_per_rep <- vapply(som_clusters, function(x) length(unique(x[is.finite(x) & !is.na(x)])), integer(1)) #clusters per replicate
  for (layer_index in seq_len(n_layers)) {
    mapvar_across_replicates <- list() #store map variance per replicate
    etasquared_across_replicates <- list() #store eta^2 per replicate
    for (replicate_index in seq_along(som_models)) {
      som_model_current <- som_models[[replicate_index]] #current model
      codebook_list_current <- kohonen::getCodes(som_model_current) #extract codes
      if (!is.list(codebook_list_current)) codebook_list_current <- list(codebook_list_current) #ensure list
      if (layer_index > length(codebook_list_current)) next #skip if replicate has fewer layers
      codebook_matrix_current <- codebook_list_current[[layer_index]] #select layer
      if (is.null(colnames(codebook_matrix_current))) colnames(codebook_matrix_current) <- paste0("V", seq_len(ncol(codebook_matrix_current))) #ensure variable names exist
      mapvar_current <- calculate.map.variance.per.variable(codebook_matrix_current, som_model_current) #compute map variance
      names(mapvar_current) <- colnames(codebook_matrix_current) #assign names
      mapvar_across_replicates[[length(mapvar_across_replicates) + 1]] <- mapvar_current #store
      neuron_cluster_labels_current <- som_clusters[[replicate_index]] #current neuron clusters
      if (!is.finite(n_clusters_per_rep[replicate_index]) || n_clusters_per_rep[replicate_index] < 2L) {
        etasquared_current <- rep(NA_real_, ncol(codebook_matrix_current)) #k = 1 replicate
        names(etasquared_current) <- colnames(codebook_matrix_current) #assign names
      } else {
        etasquared_current <- calculate.etasquared.per.variable(codebook_matrix_current, neuron_cluster_labels_current, som_model_current) #compute eta^2
        names(etasquared_current) <- colnames(codebook_matrix_current) #assign names
      }
      etasquared_across_replicates[[length(etasquared_across_replicates) + 1]] <- etasquared_current #store
    }
    if (length(mapvar_across_replicates) == 0) {
      median_map_variance_variable_importance[[layer_index]] <- NULL
    } else {
      mapvar_matrix <- do.call(rbind, mapvar_across_replicates) #replicate x variable
      median_mapvar_values <- apply(mapvar_matrix, 2, stats::median, na.rm = TRUE) #median across replicates
      if (all(is.na(median_mapvar_values))) {
        median_map_variance_variable_importance[[layer_index]] <- NULL
      } else {
        median_map_variance_variable_importance[[layer_index]] <- median_mapvar_values
      }
    }
    if (length(etasquared_across_replicates) == 0) {
      median_etasquared_variable_importance[[layer_index]] <- NULL
    } else {
      etasquared_matrix <- do.call(rbind, etasquared_across_replicates) #replicate x variable
      median_etasquared_values <- apply(etasquared_matrix, 2, stats::median, na.rm = TRUE) #median across replicates
      if (all(is.na(median_etasquared_values))) {
        median_etasquared_variable_importance[[layer_index]] <- NULL
      } else {
        median_etasquared_variable_importance[[layer_index]] <- median_etasquared_values
      }
    }
  }
  if (all(n_clusters_per_rep < 2L)) {
    median_etasquared_variable_importance <- NULL
    warning("Eta squared effect size (variable importance) could not be computed because all replicates produced k = 1")
  }
  if (all(vapply(median_map_variance_variable_importance, is.null, logical(1)))) median_map_variance_variable_importance <- NULL #collapse to NULL if empty
  if (!is.null(median_etasquared_variable_importance) && all(vapply(median_etasquared_variable_importance, is.null, logical(1)))) median_etasquared_variable_importance <- NULL #collapse to NULL if empty  
  
  # Save results
  SOM_results <- SOM.output
  SOM_results$cluster_assignment <- cluster_assignment
  SOM_results$BIC_values <- BIC_values
  SOM_results$support_values <- support_values
  SOM_results$support_label <- support_label
  SOM_results$support_higher_is_better <- support_higher_is_better
  SOM_results$ancestry_matrix <- ancestry_matrix
  SOM_results$replicate_ancestry_matrices <- replicate_ancestry_matrices
  SOM_results$mean_assignment_margin <- mean_assignment_margin
  SOM_results$mean_normalized_assignment_entropy <- mean_normalized_assignment_entropy
  SOM_results$mean_mean_assignment_margin <- if (all(is.na(mean_assignment_margin))) NA_real_ else mean(mean_assignment_margin, na.rm = TRUE)
  SOM_results$mean_mean_normalized_assignment_entropy <- if (all(is.na(mean_normalized_assignment_entropy))) NA_real_ else mean(mean_normalized_assignment_entropy, na.rm = TRUE)
  SOM_results$clustering.SOM.set.seed.N <- set.seed.N
  SOM_results$clustering.SOM.args <- list(
    max.k = max.k,
    set.k = set.k,
    clustering.method = clustering.method,
    BIC.thresh = BIC.thresh,
    quantization.error.quantile = quantization.error.quantile,
    topographic.error.quantile = topographic.error.quantile,
    set.seed.N = set.seed.N
  )
  SOM_results$optim_k_vals <- optim_k_vals
  SOM_results$optim_k_mean <- optim_k_mean
  SOM_results$optim_k_summary <- optim_k_summary
  SOM_results$max_k <- max.k
  SOM_results$set_k <- set.k
  SOM_results$som_clusters <- som_clusters
  SOM_results$cluster_gridcell_assignments <- cluster_gridcell_assignments
  if (is.null(median_etasquared_variable_importance)) {
    SOM_results$median_etasquared_variable_importance <- list() #eta^2 unavailable (all replicates produced k = 1)
  } else {
    SOM_results$median_etasquared_variable_importance <- median_etasquared_variable_importance #median eta^2 per layer
  }
  if (is.null(median_map_variance_variable_importance)) {
    SOM_results$median_map_variance_variable_importance <- list() #map variance unavailable
  } else {
    SOM_results$median_map_variance_variable_importance <- median_map_variance_variable_importance #median map variance per layer
  }
  SOM_results$quantization.error.quantile <- quantization.error.quantile
  SOM_results$topographic.error.quantile <- topographic.error.quantile
  SOM_results$retained_replicates <- retained_replicates
  SOM_results$removed_replicates <- removed_replicates
  SOM_results$N_replicates_retained <- length(retained_replicates)
  if (!is.null(SOM.output$quantization_error)) SOM_results$quantization_error_retained <- SOM.output$quantization_error
  if (!is.null(SOM.output$topographic_error)) SOM_results$topographic_error_retained <- SOM.output$topographic_error
  
  # Return results
  return(SOM_results)
}


## Function to plot Structure-like barplots
plot.structure.SOM <- function(SOM.output,
                               col.pal = viridis::viridis, #color palette
                               save = F, #save plot
                               overwrite = T, #overwrite plot if already present (only if saving plot)
                               plot.type = "svg", #plot file type (choose: png, jpg or svg; only if saving plot)
                               file.name = NULL, #plot file name (if NULL, default name is created; only if saving plot)
                               width = 10, #plot width in cm (only if saving plot)
                               height = 15, #plot height in cm (only if saving plot)
                               resolution = 300, #plot resolution in dpi (only if saving plot)
                               bottom.margin = 5, #bottom margin
                               left.margin = 5, #left margin
                               top.margin = 1.5, #top margin
                               right.margin = 1.5, #right margin
                               Individual.labels.font.size = 0.45, #font size of individual labels on axis
                               Y.axis.title = "Cluster assignment coefficient", #set y axis title
                               sort.by.col = 1, #specify integer giving column index of ancestry matrix for ordering rows of ancestry matrix (if NULL, hierarchical ordering is performed)
                               linkage.method = "single", #agglomeration method used for hierarchical clustering (see hclust function)
                               bar.border.col = NULL, #color of separation lines (e.g. "black", "gray30"); NULL = no border
                               bar.border.lwd = 1 #line width of separation lines (ignored if color is NULL)
) {
  
  # Reset plotting parameters
  old_dev <- dev.cur()
  old_plotting_parameters <- par(no.readonly = TRUE)
  on.exit({
    if (dev.cur() == old_dev) par(old_plotting_parameters)
  }, add = TRUE)
  
  # Validate SOM.output
  if (is.null(SOM.output$ancestry_matrix) || !is.matrix(SOM.output$ancestry_matrix)) {
    stop("Plotting aborted: ancestry_matrix of SOM.output not valid - check SOM.output or rerun run.SOM")
  }
  if (nrow(SOM.output$ancestry_matrix) < 2) {
    stop("Plotting aborted: ancestry_matrix must have at least 2 rows (individuals)")
  }
  if (ncol(SOM.output$ancestry_matrix) < 1) {
    stop("Plotting aborted: ancestry_matrix must have at least one column (clusters)")
  }
  if (all(is.na(SOM.output$ancestry_matrix))) {
    stop("Plotting aborted: ancestry_matrix has all NA values")
  }
  
  # Validate specified color palette
  viridis_palettes <- list(
    viridis::viridis,
    viridis::magma,
    viridis::plasma,
    viridis::inferno,
    viridis::cividis,
    viridis::rocket,
    viridis::mako,
    viridis::turbo
  )
  if (!any(vapply(viridis_palettes, identical, logical(1), col.pal))) {
    stop("Plotting aborted: col.pal must viridis palette - viridis, magma, plasma, inferno, cividis, rocket, mako or turbo")
  }
  
  # Validate specified save
  if (!is.logical(save) || length(save) != 1) {
    stop("Plotting aborted: save must be TRUE or FALSE")
  }
  
  # Validate specified overwrite
  if (save) {
    if (!is.logical(overwrite) || length(overwrite) != 1) {
      stop("Plotting aborted: overwrite must be TRUE or FALSE")
    }
  }
  
  # Validate specified plot.type
  if (save) {
    allowed_plot.types <- c("svg", "png", "jpg")
    if (!is.character(plot.type) || length(plot.type) != 1 || is.na(plot.type) || !(plot.type %in% allowed_plot.types)) {
      stop("Plotting aborted: plot.type must be one of 'svg', 'png', or 'jpg'")
    }
  }
  
  # Validate specified file.name
  if (save) {
    if (!is.null(file.name) && (!is.character(file.name) || length(file.name) != 1 || is.na(file.name))) {
      stop("Plotting aborted: file.name must be NULL or single character string")
    }
  }
  
  # Validate specified width and height (reasonable values: 4–50 cm)
  if (save) {
    if (!is.numeric(width) || length(width) != 1 || is.na(width) || width <= 0) {
      stop("Plotting aborted: width must be a single positive number (cm)")
    }
    if (width < 4) {
      message("Warning: width is very small (", width, " cm) – plot may be hard to read")
    }
    if (width > 50) {
      message("Warning: width is very large (", width, " cm) – plot may be unwieldy")
    }
    if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
      stop("Plotting aborted: height must be a single positive number (cm)")
    }
    if (height < 4) {
      message("Warning: height is very small (", height, " cm) – plot may be hard to read")
    }
    if (height > 50) {
      message("Warning: height is very large (", height, " cm) – plot may be unwieldy")
    }
  }
  
  # Validate specified resolution (reasonable range: 72–1200 dpi)
  if (save) {
    if (!is.numeric(resolution) || length(resolution) != 1 || is.na(resolution) || resolution < 72) {
      stop("Plotting aborted: resolution must be a single number ≥ 72 (dpi)")
    }
    if (resolution > 1200) {
      message("Warning: resolution is very high (", resolution, " dpi) – file may be huge")
    }
  }
  
  # Validate specified margins (reasonable range: 0–10)
  margin.list <- c(bottom.margin, left.margin, top.margin, right.margin)
  margin.names <- c("bottom.margin", "left.margin", "top.margin", "right.margin")
  for (i in seq_along(margin.list)) {
    if (!is.numeric(margin.list[i]) || length(margin.list[i]) != 1 || is.na(margin.list[i])) {
      stop("Plotting aborted: ", margin.names[i], " must be a single numeric value")
    }
    if (margin.list[i] < 0) {
      stop("Plotting aborted: ", margin.names[i], " must be ≥ 0")
    }
    if (margin.list[i] > 10) {
      message("Warning: ", margin.names[i], " is large (", margin.list[i], ") – plot area may shrink")
    }
  }
  
  # Validate specified Individual.labels.font.size (reasonable: 0.2–4)
  if (!is.numeric(Individual.labels.font.size) || length(Individual.labels.font.size) != 1 || is.na(Individual.labels.font.size) || Individual.labels.font.size <= 0) {
    stop("Plotting aborted: Individual.labels.font.size must be a single positive number")
  }
  if (Individual.labels.font.size < 0.2) {
    message("Warning: Individual.labels.font.size is very small (", Individual.labels.font.size, ") – labels may not be readable")
  }
  if (Individual.labels.font.size > 4) {
    message("Warning: Individual.labels.font.size is large (", Individual.labels.font.size, ") – labels may overlap")
  }
  
  # Validate specified sort.by.col
  if (!is.null(sort.by.col)) {
    if (!is.numeric(sort.by.col) || length(sort.by.col) != 1 ||
        is.na(sort.by.col) || sort.by.col < 1 || sort.by.col > ncol(SOM.output$ancestry_matrix) ||
        (sort.by.col %% 1 != 0)) {
      stop(paste0("Plotting aborted: sort.by.col must be integer between 1 and ", ncol(SOM.output$ancestry_matrix), " or NULL"))
    }
  }
  
  # Validate specified linkage.method
  allowed_linkage <- c("single", "complete", "average", "ward.D", "ward.D2", "mcquitty", "median", "centroid")
  if (!linkage.method %in% allowed_linkage) {
    stop("Plotting aborted: linkage.method must be one of: ", paste(allowed_linkage, collapse = ", "))
  }
  
  # Validate bar.border.col
  if (!is.null(bar.border.col)) {
    if (!is.character(bar.border.col) || length(bar.border.col) != 1 || is.na(bar.border.col)) {
      stop("Plotting aborted: bar.border.col must be NULL or single character string (e.g. 'white' or 'black')")
    }
  }
  
  # Validate bar.border.lwd
  if (!is.null(bar.border.lwd)) {
    if (!is.numeric(bar.border.lwd) || length(bar.border.lwd) != 1 || is.na(bar.border.lwd) || bar.border.lwd <= 0) {
      stop("Plotting aborted: bar.border.lwd must be a single positive number")
    }
    if (bar.border.lwd > 5) {
      message("Warning: bar.border.lwd is large (", bar.border.lwd, ") — lines may obscure adjacent bars")
    }
  }
  
  # If there is only one cluster, skip plotting
  if (ncol(SOM.output$ancestry_matrix) == 1) {
    stop("Only one cluster detected - skipping Structure-like plot")
  }
  
  # Order rows of ancestry_matrix
  if (!is.null(sort.by.col)) {
    if (sort.by.col > ncol(SOM.output$ancestry_matrix)) {
      stop(paste0("Plotting aborted: sort.by.col exceeds number of columns in ancestry_matrix - select number from 1 - ", ncol(SOM.output$ancestry_matrix), " or NULL to prevent sorting"))
    }
    ancestry_proportions <- SOM.output$ancestry_matrix[order(SOM.output$ancestry_matrix[, sort.by.col]), ]
  } else {
    ancestry_proportions <- SOM.output$ancestry_matrix #no sorting when set as NULL
  }
  
  # Perform hierarchical clustering on distance matrix
  cluster_order <- stats::hclust(dist(ancestry_proportions), method = linkage.method)$order
  SOM_ancestry_proportions <- ancestry_proportions[cluster_order, ]
  
  # Generate layer colors
  layer_colors <- col.pal(ncol(SOM_ancestry_proportions))
  
  # Save file
  if (save) {
    if (is.null(file.name)) {
      file.name <- paste0("SOM_structure_plot_", paste(SOM.output$input_data_names, collapse = "_"), ".", plot.type)
    }
    if (file.exists(file.name) && !overwrite) {
      stop("Plotting aborted: ", paste(file.name), "already exists - set overwrite = TRUE to overwrite")
    }
    if (plot.type == "svg") {
      svg(file.name, 
          width = width / 2.54, 
          height = height / 2.54)
    } else if (plot.type == "png") {
      png(file.name,
          width = width, 
          height = height, 
          res = resolution, 
          units = "cm")
    } else if (plot.type == "jpg") {
      jpeg(file.name, 
           width = width,
           height = height, 
           res = resolution, 
           units = "cm")
    }
  }
  
  # Set plot layout
  par(mfrow = c(1, 1),
      mar = c(5.1, 4.1, 4.1, 2.1),
      oma = c(0, 0, 0, 0))
  
  # Create function to plot Structure-like plot (based on conStruct)
  plot.Structure <- function(admix.proportions, 
                             mar = c(2, 4, 2, 2), 
                             sample.order = NULL,
                             layer.order = NULL, 
                             sample.names = NULL, 
                             sort.by = NULL,
                             Y.axis.title = "Ancestry",
                             layer.colors = NULL,
                             bar.border.col = NULL,
                             bar.border.lwd = 0.5)
  {
    K <- ncol(admix.proportions) #number of clusters (layers)
    N <- nrow(admix.proportions) #number of individuals (samples)
    par(mar = mar) #set plot margins
    
    if (is.null(layer.order)) {
      layer.order <- seq(1:K) #default cluster order (1 to K)
    }
    if (is.null(sample.order)) {
      sample.order <- seq(1:N) #default sample order (1 to N)
    }
    use.colors <- layer.colors[1:K][layer.order]
    plot(0,
         xlim = c(0, N),
         ylim = c(0, 1),
         type = "n",
         ylab = Y.axis.title,
         xlab = "",
         xaxt = "n") #create empty plot with custom axes
    
    plotting.admix.props <- apply(cbind(0, admix.proportions[, layer.order]), 1, cumsum) #compute cumulative admixture proportions for polygons
    
    # Draw bars with no internal segment borders
    for (i in 1:K) {
      for (j in 1:N) {
        polygon(
          x = c(j - 1, j, j, j - 1),
          y = c(plotting.admix.props[i, sample.order[j]],
                plotting.admix.props[i, sample.order[j]],
                plotting.admix.props[i + 1, sample.order[j]],
                plotting.admix.props[i + 1, sample.order[j]]),
          col = use.colors[i],
          border = NA # disable within-bar borders
        )
      }
    }
    
    # Optional: draw vertical lines between bars (samples)
    if (!is.null(bar.border.col)) {
      for (j in 2:N) {
        segments(x0 = j - 1, y0 = 0, x1 = j - 1, y1 = 1,
                 col = bar.border.col, lwd = bar.border.lwd)
      }
    }
    
    
    if (!is.null(sample.names)) {
      axis(side = 1,
           at = seq(1:N) - 0.5,
           labels = sample.names[sample.order],
           cex.axis = Individual.labels.font.size,
           las = 2) #add sample names to x-axis
    }
    
    return(invisible(NULL)) #return invisible
  }
  
  # Generate structure plot
  plot.Structure(admix.proportions = SOM_ancestry_proportions,
                 sample.names = rownames(SOM_ancestry_proportions),
                 mar = c(bottom.margin,
                         left.margin,
                         top.margin,
                         right.margin),
                 layer.colors = layer_colors,
                 Y.axis.title = Y.axis.title,
                 bar.border.col = bar.border.col,
                 bar.border.lwd = bar.border.lwd,
                 sort.by = sort.by.col)
  
  # Close graphics device
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}


## Function to plot learning progress for each SOM matrix
plot.learning.SOM <- function(SOM.output, 
                              col.pal = viridis::turbo, #set color palette
                              save = F, #option to save plot
                              overwrite = T, #option to overwrite plot if it already exists (only if saving plot)
                              plot.type = "svg", #options: "svg", "png", "jpg" (only if saving plot)
                              file.name = NULL, #set plot file.name (if NULL, default plot file.name is used; only if saving plot)
                              width = 20, #plot width in cm (only if saving plot)
                              height = 15, #plot height in cm (only if saving plot)
                              resolution = 300, #plot resolution in dpi (only if saving plot)
                              bottom.margin = 5, #bottom margin
                              left.margin = 5, #left margin
                              top.margin = 3, #top margin
                              right.margin = 2.5, #right margin
                              lines.alpha = 0.4, #transparency for plot lines
                              lines.thickness = 0.9, #thickness for plot lines
                              title = NULL, #main title of plot (if NULL, default title name is used)
                              legend.position = "topright", #position of legend in plot
                              legend.lines.thickness = 3, #thickness for lines in legend
                              x.axis.label = "Training steps", #x axis label
                              y.axis.label = "Learning rate change" #y axis label
) {
  
  # Reset plotting parameters
  old_dev <- dev.cur()
  old_plotting_parameters <- par(no.readonly = TRUE)
  on.exit({
    if (dev.cur() == old_dev) par(old_plotting_parameters)
  }, add = TRUE)
  
  # Validate learning_values_list
  if (!is.list(SOM.output$learning_values_list)) {
    stop("Plotting aborted: learning_values_list must be list of matrices")
  }
  
  # Validate specified color palette
  viridis_palettes <- list(
    viridis::viridis,
    viridis::magma,
    viridis::plasma,
    viridis::inferno,
    viridis::cividis,
    viridis::rocket,
    viridis::mako,
    viridis::turbo
  )
  if (!any(vapply(viridis_palettes, identical, logical(1), col.pal))) {
    stop("Plotting aborted: col.pal must viridis palette - viridis, magma, plasma, inferno, cividis, rocket, mako or turbo")
  }
  
  # Validate specified save
  if (!is.logical(save) || length(save) != 1 || is.na(save)) {
    stop("Plotting aborted: save must be TRUE or FALSE")
  }
  
  # Validate specified overwrite
  if (save) {
    if (!is.logical(overwrite) || length(overwrite) != 1 || is.na(overwrite)) {
      stop("Plotting aborted: overwrite must be TRUE or FALSE")
    }
  }
  
  # Validate specified plot.type
  if (save) {
    allowed_plot.types <- c("svg", "png", "jpg")
    if (!is.character(plot.type) || length(plot.type) != 1 || is.na(plot.type) || !(plot.type %in% allowed_plot.types)) {
      stop("Plotting aborted: plot.type must be one of 'svg', 'png', or 'jpg'")
    }
  }
  
  # Validate specified file.name
  if (save) {
    if (!is.null(file.name) && (!is.character(file.name) || length(file.name) != 1 || is.na(file.name))) {
      stop("Plotting aborted: file.name must be NULL or single character string")
    }
  }
  
  # Validate specified width and height (reasonable values: 4–50 cm)
  if (save) {
    if (!is.numeric(width) || length(width) != 1 || is.na(width) || width <= 0) {
      stop("Plotting aborted: width must be a single positive number (cm)")
    }
    if (width < 4) {
      message("Warning: width is very small (", width, " cm) – plot may be hard to read")
    }
    if (width > 50) {
      message("Warning: width is very large (", width, " cm) – plot may be unwieldy")
    }
    if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
      stop("Plotting aborted: height must be a single positive number (cm)")
    }
    if (height < 4) {
      message("Warning: height is very small (", height, " cm) – plot may be hard to read")
    }
    if (height > 50) {
      message("Warning: height is very large (", height, " cm) – plot may be unwieldy")
    }
  }
  
  # Validate specified resolution (reasonable range: 72–1200 dpi)
  if (save) {
    if (!is.numeric(resolution) || length(resolution) != 1 || is.na(resolution) || resolution < 72) {
      stop("Plotting aborted: resolution must be a single number ≥ 72 (dpi)")
    }
    if (resolution > 1200) {
      message("Warning: resolution is very high (", resolution, " dpi) – file may be huge")
    }
  }
  
  # Validate specified margins (reasonable range: 0–10)
  margin.list <- c(bottom.margin, left.margin, top.margin, right.margin)
  margin.names <- c("bottom.margin", "left.margin", "top.margin", "right.margin")
  for (i in seq_along(margin.list)) {
    if (!is.numeric(margin.list[i]) || length(margin.list[i]) != 1 || is.na(margin.list[i])) {
      stop("Plotting aborted: ", margin.names[i], " must be a single numeric value")
    }
    if (margin.list[i] < 0) {
      stop("Plotting aborted: ", margin.names[i], " must be ≥ 0")
    }
    if (margin.list[i] > 10) {
      message("Warning: ", margin.names[i], " is large (", margin.list[i], ") – plot area may shrink")
    }
  }
  
  # Validate specified lines.alpha (must be 0–1)
  if (!is.numeric(lines.alpha) || length(lines.alpha) != 1 || is.na(lines.alpha) ||
      lines.alpha < 0 || lines.alpha > 1) {
    stop("Plotting aborted: lines.alpha must be numeric value between 0 and 1")
  }
  
  # Validate specified lines.thickness (must be positive)
  if (!is.numeric(lines.thickness) || length(lines.thickness) != 1 || is.na(lines.thickness) ||
      lines.thickness <= 0) {
    stop("Plotting aborted: lines.thickness must be a single positive numeric value")
  }
  
  # Validate specified title (must be NULL or character)
  if (!is.null(title) && (!is.character(title) || length(title) != 1 || is.na(title))) {
    stop("Plotting aborted: title must be NULL or single character string")
  }
  
  # Validate specified legend.position
  allowed.legend.positions <- c("topright", "topleft", "bottomright", "bottomleft", 
                                "right", "left", "top", "bottom", "center")
  if (!is.character(legend.position) || length(legend.position) != 1 || is.na(legend.position) ||
      !(legend.position %in% allowed.legend.positions)) {
    stop(paste0("Plotting aborted: legend.position must be one of ", 
                paste(allowed.legend.positions, collapse = ", ")))
  }
  
  # Validate specified legend.lines.thickness (must be positive)
  if (!is.numeric(legend.lines.thickness) || length(legend.lines.thickness) != 1 || is.na(legend.lines.thickness) ||
      legend.lines.thickness <= 0) {
    stop("Plotting aborted: legend.lines.thickness must be a single positive numeric value")
  }
  
  # Validate specified x.axis.label and y.axis.label (must be character)
  if (!is.character(x.axis.label) || length(x.axis.label) != 1 || is.na(x.axis.label)) {
    stop("Plotting aborted: x.axis.label must be a single character string")
  }
  if (!is.character(y.axis.label) || length(y.axis.label) != 1 || is.na(y.axis.label)) {
    stop("Plotting aborted: y.axis.label must be a single character string")
  }
  
  # Check if learning_values or learning_values_list exists and is either matrix or list
  if ("learning_values" %in% names(SOM.output)) {
    SOM.output$learning_values_list <- list(SOM.output$learning_values) #convert to list for single-layer case
  } else if ("learning_values_list" %in% names(SOM.output)) { #for multi-layer case, no changes needed
  } else {
    stop("Plotting aborted: neither learning_values nor learning_values_list found in SOM.output - recheck or rerun run.SOM")
  }
  
  # Set default file.name for saving
  if (save) {
    if (is.null(file.name)) { 
      file.name <- paste0("SOM_learning_plot_", paste(SOM.output$input_data_names, collapse = "_"), ".", plot.type)
    }
  }
  
  # Check if file already exists and overwrite option is set to FALSE
  if (save) {
    if (!overwrite && file.exists(file.name)) {
      stop(sprintf("Plotting aborted: file '%s' already exists - skipping plot saving", file.name))
    }
  }
  
  # Convert data.frames to matrices if necessary
  SOM.output$learning_values_list <- lapply(SOM.output$learning_values_list, function(x) {
    if (is.data.frame(x)) {
      return(as.matrix(x))
    }
    return(x)
  })
  
  # Extract matrix names (check for multi-layer or single-layer)
  if ("input_data_names" %in% names(SOM.output)) {
    matrix_names <- SOM.output$input_data_names
  } else {
    stop("Plotting aborted: matrix names (input_data_names) not found in provided SOM.output")
  }
  
  # Check if matrix names match number of layers
  if (length(matrix_names) != length(SOM.output$learning_values_list)) {
    stop("Plotting aborted: number of matrix names (input_data_names) does not match number of matrices")
  }
  
  # Determine global ylim
  global_ylim <- range(unlist(lapply(SOM.output$learning_values_list, function(mat) 
    range(mat, na.rm = T))), na.rm = T) * c(0.93, 1.07)
  
  # Set legend and main plot title based on number of layers
  if (length(SOM.output$learning_values_list) == 1) {
    main_title <- "Training progress across input layer"
    legend_title <- "Layer"
  } else {
    main_title <- "Training progress across input layers"
    legend_title <- "Layers"
  }
  
  # Use provided title if specified, otherwise use default
  main_title <- ifelse(is.null(title), yes = main_title, no = title)
  
  # Save plot if requested
  if (save) {
    if (plot.type == "svg") {
      svg(file.name, 
          width = width / 2.54, 
          height = height / 2.54)
    } else if (plot.type == "png") {
      png(file.name, 
          width = width, 
          height = height, 
          units = "cm", 
          res = resolution)
    } else if (plot.type == "jpg") {
      jpeg(file.name, 
           width = width, 
           height = height, 
           units = "cm", 
           res = resolution)
    }
  }
  
  # Prepare base plot
  first_matrix <- SOM.output$learning_values_list[[1]]
  par(mfrow = c(1, 1), 
      mar = c(bottom.margin, left.margin, top.margin, right.margin))
  plot(NULL, 
       xlim = c(1, nrow(first_matrix)), 
       ylim = global_ylim, 
       xlab = x.axis.label, 
       ylab = y.axis.label, 
       main = main_title,
       axes = T)
  
  # Plot each matrix
  layer_colors <- setNames(col.pal(length(matrix_names)), matrix_names)
  for (i in seq_along(SOM.output$learning_values_list)) {
    mat <- SOM.output$learning_values_list[[i]]
    current_layer_name <- matrix_names[i]
    for (j in 1:ncol(mat)) 
      lines(mat[, j], 
            col = grDevices::adjustcolor(layer_colors[current_layer_name], alpha.f = lines.alpha),
            lwd = lines.thickness)
  }
  
  # Add legend
  legend(legend.position, 
         legend = matrix_names, 
         col = layer_colors[matrix_names], 
         lty = 1, 
         lwd = legend.lines.thickness,
         title = legend_title)
  
  # Close graphics device
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}


## Function to plot mean pairwise distance across layers
plot.layer.distance.scale.SOM <- function(SOM.output,
                                          col.pal = viridis::turbo, #color palette
                                          save = F, #option to save plot
                                          overwrite = T, #option to overwrite plot if it already exists (only if saving plot)
                                          plot.type = "svg", #plot type options: "svg", "png", "jpg" (only if saving plot)
                                          file.name = NULL, #set file.name for saving (if NULL, default plot file.name is used; only if saving plot)
                                          width = 20, #plot width in cm (only if plot is saved)
                                          height = 15, #plot height in cm (only if plot is saved)
                                          resolution = 300, #plot resolution in dpi (only if plot is saved)
                                          bottom.margin = 3, #bottom margin
                                          left.margin = 5, #left margin
                                          top.margin = 3.5, #top margin
                                          right.margin = 0.5, #right margin
                                          title = "Layer distance scale across layers", #plot title name
                                          y_axis_label = "Average pairwise distance" #set y axis label
) {
  
  # Reset plotting parameters
  old_dev <- dev.cur()
  old_plotting_parameters <- par(no.readonly = TRUE)
  on.exit({
    if (dev.cur() != old_dev) dev.off()
    par(old_plotting_parameters)
  })
  
  # Validate SOM.output
  if (is.null(SOM.output$distance_weights_matrix)) {
    stop("Plotting aborted: SOM.output does not contain distance_weights_matrix")
  }
  
  # Extract replicate-wise distance weights matrix
  d.mat <- SOM.output$distance_weights_matrix
  if (!is.matrix(d.mat)) d.mat <- as.matrix(d.mat)
  
  # Require multilayer input
  if (ncol(d.mat) < 2) {
    stop("Plotting aborted: layer plot requires at least two layers")
  }
  
  # Extract matrix names
  if ("input_data_names" %in% names(SOM.output)) {
    matrix_names <- SOM.output$input_data_names
  } else {
    matrix_names <- colnames(d.mat)
  }
  
  if (is.null(matrix_names) || length(matrix_names) != ncol(d.mat)) {
    matrix_names <- paste0("Layer", seq_len(ncol(d.mat)))
  }
  
  # Convert distance weights back to mean pairwise distances
  mean_pairwise_distance_matrix <- 1 / d.mat
  colnames(mean_pairwise_distance_matrix) <- matrix_names
  
  # Calculate mean pairwise distance across replicates
  mean_pairwise_distance <- colMeans(mean_pairwise_distance_matrix, na.rm = TRUE)
  
  # Order layers by descending mean pairwise distance
  order_idx <- order(mean_pairwise_distance, decreasing = TRUE)
  mean_pairwise_distance <- mean_pairwise_distance[order_idx]
  matrix_names <- matrix_names[order_idx]
  
  # Set default file name
  if (is.null(file.name)) {
    file.name <- paste0("plot_layer_distance_scale.", plot.type)
  }
  
  # Check overwrite settings
  if (save && file.exists(file.name) && !overwrite) {
    stop(paste("Plotting aborted:", file.name, "already exists and overwrite = FALSE"))
  }
  
  # Save plot if requested
  if (save) {
    if (plot.type == "svg") {
      svg(file.name,
          width = width / 2.54,
          height = height / 2.54)
    } else if (plot.type == "png") {
      png(file.name,
          width = width,
          height = height,
          units = "cm",
          res = resolution)
    } else if (plot.type == "jpg") {
      jpeg(file.name,
           width = width,
           height = height,
           units = "cm",
           res = resolution)
    } else {
      stop("Plotting aborted: plot.type must be 'svg', 'png', or 'jpg'")
    }
  }
  
  # Prepare plot
  par(mfrow = c(1, 1),
      mar = c(bottom.margin, left.margin, top.margin, right.margin))
  layer_colors <- setNames(col.pal(length(SOM.output$input_data_names)), SOM.output$input_data_names)
  
  # Create barplot
  barplot(height = mean_pairwise_distance,
          col = layer_colors[matrix_names],
          border = "black",
          main = title,
          ylab = y_axis_label,
          names.arg = matrix_names)
  
  # Close graphics device
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}


## Function to evaluate K-values
plot.K.SOM <- function(SOM.output,
                       col.pal = viridis::magma, #set color palette (default magma)
                       save = FALSE, #set option to save plot (default: FALSE)
                       overwrite = TRUE, #option to overwrite plot (only if saving plot)
                       plot.type = "svg", #plot type: "svg", "png", "jpg" (only if saving plot)
                       file.name = NULL, #custom file name (only if saving plot)
                       width = 10, #plot width in cm (only if saving plot)
                       height = 15, #plot height in cm (only if saving plot)
                       resolution = 300, #plot resolution in dpi (only if saving plot)
                       bottom.margin = 4, #bottom margin
                       left.margin = 5.5, #left margin
                       top.margin = 2.5, #top margin
                       right.margin = 2.5, #right margin
                       N.axis.labels.BIC.plot = 4, #number of axis labels for support plot
                       N.axis.labels.deltaBIC.plot = 4, #number of axis labels for delta-BIC plot
                       round.axis.labels.BIC.plot = 0, #rounding axis labels of support plot to X digits
                       round.axis.labels.deltaBIC.plot = 0, #rounding axis labels of delta-BIC plot to X digits
                       title = "Number of clusters (k)") #change plot title 
{
  
  # Reset plotting parameters
  old_dev <- dev.cur()
  old_plotting_parameters <- par(no.readonly = TRUE)
  on.exit({
    if (dev.cur() == old_dev) par(old_plotting_parameters)
  }, add = TRUE)
  
  # Validate SOM.output
  if (is.null(SOM.output$N_replicates)) {
    stop("Plotting aborted: N_replicates could not be found in SOM output - check if clustering.SOM was run")
  } else if (is.null(SOM.output$optim_k_vals)) {
    stop("Plotting aborted: optim_k_vals could not be found in SOM output - check if clustering.SOM was run")
  } else if (is.null(SOM.output$max_k)) {
    stop("Plotting aborted: max_k could not be found in SOM output - check if clustering.SOM was run")
  }
  
  # Extract clustering method
  clustering.method <- if (!is.null(SOM.output$clustering.SOM.args$clustering.method)) {
    SOM.output$clustering.SOM.args$clustering.method
  } else {
    "unknown"
  }
  
  # Set default file name
  if (is.null(file.name)) {
    file.name <- paste0("K_plot.", plot.type)
  }
  
  # Open graphics device if saving
  if (save) {
    if (file.exists(file.name) && !overwrite) {
      stop(paste("Plotting aborted: file already exists:", file.name))
    }
    if (plot.type == "svg") {
      svg(filename = file.name, width = width / 2.54, height = height / 2.54)
    } else if (plot.type == "png") {
      png(filename = file.name, width = width, height = height, units = "cm", res = resolution)
    } else if (plot.type == "jpg") {
      jpeg(filename = file.name, width = width, height = height, units = "cm", res = resolution)
    } else {
      stop("Plotting aborted: plot.type must be 'svg', 'png', or 'jpg'")
    }
  }
  
  # Determine support values to plot
  plot_support_values <- NULL
  support_label <- NULL
  support_higher_is_better <- NA
  
  if (!is.null(SOM.output$support_values) && is.matrix(SOM.output$support_values) &&
      any(is.finite(SOM.output$support_values), na.rm = TRUE)) {
    plot_support_values <- SOM.output$support_values
    support_label <- SOM.output$support_label
    support_higher_is_better <- SOM.output$support_higher_is_better
  } else if (!is.null(SOM.output$BIC_values) && is.matrix(SOM.output$BIC_values) &&
             any(is.finite(SOM.output$BIC_values), na.rm = TRUE)) {
    plot_support_values <- SOM.output$BIC_values
    support_label <- "BIC"
    support_higher_is_better <- FALSE
  }
  
  support_available <- !is.null(plot_support_values)
  support_is_BIC <- isTRUE(identical(support_label, "BIC"))
  
  # Set colors
  k.cols <- col.pal(SOM.output$max_k)
  
  # Create helper for support axis breaks
  make.axis.breaks <- function(values, n.breaks = 4, digits = 0) {
    finite_values <- values[is.finite(values) & !is.na(values)]
    if (length(finite_values) == 0) return(NULL)
    axis_breaks <- pretty(range(finite_values), n = n.breaks)
    round(axis_breaks, digits)
  }
  
  # Plot support values if available
  plot.support.panel <- function(values_matrix, ylab_text, axis_digits = 0) {
    support_vals <- as.numeric(values_matrix)
    support_breaks <- make.axis.breaks(support_vals, n.breaks = N.axis.labels.BIC.plot, digits = axis_digits)
    
    boxplot(t(values_matrix)[, 1:(SOM.output$max_k), drop = FALSE],
            outline = FALSE,
            notch = FALSE,
            axes = FALSE,
            ylab = ylab_text,
            whisklty = 1,
            staplelty = 1,
            col = k.cols)
    
    axis(1, at = 1:(SOM.output$max_k), labels = 1:(SOM.output$max_k))
    if (!is.null(support_breaks)) axis(2, at = support_breaks, labels = support_breaks, las = 3)
    title(title, line = 0)
  }
  
  # Plot frequency panel
  plot.frequency.panel <- function() {
    optim_k_vector <- as.numeric(SOM.output$optim_k_vals)
    optim_k_vector <- optim_k_vector[is.finite(optim_k_vector) & !is.na(optim_k_vector)]
    barplot_data <- table(factor(optim_k_vector, levels = 1:(SOM.output$max_k))) / SOM.output$N_replicates
    
    bar_midpoints <- barplot(barplot_data,
                             ylab = "Frequency of selected K",
                             ylim = c(0, 1),
                             col = k.cols,
                             axes = FALSE,
                             axisnames = FALSE)
    axis(1, at = bar_midpoints, labels = seq_len(SOM.output$max_k))
    axis(2, las = 3)
  }
  
  # Plot delta-BIC panel if available
  plot.deltaBIC.panel <- function() {
    d_wss_raw <- apply(SOM.output$BIC_values, 2, function(x) {
      if (sum(is.finite(x) & !is.na(x)) < 3) {
        rep(NA_real_, max(0, length(x) - 2))
      } else {
        diff(diff(x))
      }
    })
    
    if (length(d_wss_raw) == 0 || all(is.na(d_wss_raw)) || all(!is.finite(d_wss_raw))) {
      return(FALSE)
    }
    
    if (is.null(dim(d_wss_raw))) {
      d_wss <- matrix(d_wss_raw,
                      nrow = 1,
                      dimnames = list(2:(SOM.output$max_k - 1),
                                      colnames(SOM.output$BIC_values)))
    } else {
      d_wss <- d_wss_raw
      rownames(d_wss) <- 2:(SOM.output$max_k - 1)
    }
    
    delta_vals <- as.numeric(d_wss)
    delta_breaks <- make.axis.breaks(delta_vals,
                                     n.breaks = N.axis.labels.deltaBIC.plot,
                                     digits = round.axis.labels.deltaBIC.plot)
    
    boxplot(t(d_wss),
            outline = FALSE,
            notch = FALSE,
            axes = FALSE,
            ylab = expression(Delta*BIC),
            whisklty = 1,
            staplelty = 1,
            col = k.cols[2:(SOM.output$max_k - 1)])
    
    axis(1, at = 1:nrow(d_wss), labels = rownames(d_wss))
    if (!is.null(delta_breaks)) axis(2, at = delta_breaks, labels = delta_breaks, las = 3)
    
    return(TRUE)
  }
  
  # Report unavailable panels via output messages
  if (!support_available) {
    message(paste0("plot.K.SOM: no finite support values available for clustering.method = '",
                   clustering.method,
                   "' - support panel will be omitted"))
  }
  
  if (support_available && !support_is_BIC) {
    message(paste0("plot.K.SOM: clustering.method = '",
                   clustering.method,
                   "' is not BIC-based - delta-BIC panel will be omitted"))
  }
  
  if (support_is_BIC && SOM.output$max_k <= 2) {
    message("plot.K.SOM: delta-BIC panel requires max.k >= 3 - delta-BIC panel will be omitted")
  }
  
  # Decide layout and plot
  if (support_available && support_is_BIC && SOM.output$max_k > 2) {
    par(mfrow = c(3, 1),
        bty = "n",
        mar = c(bottom.margin, left.margin, top.margin, right.margin))
    plot.support.panel(plot_support_values, ylab_text = "BIC", axis_digits = round.axis.labels.BIC.plot)
    deltaBIC_plotted <- plot.deltaBIC.panel()
    if (!deltaBIC_plotted) {
      message("plot.K.SOM: insufficient finite BIC values for delta-BIC panel - plotting support panel and K-frequency only")
      par(mfrow = c(2, 1),
          bty = "n",
          mar = c(bottom.margin, left.margin, top.margin, right.margin))
      plot.support.panel(plot_support_values, ylab_text = "BIC", axis_digits = round.axis.labels.BIC.plot)
      plot.frequency.panel()
    } else {
      plot.frequency.panel()
    }
  } else if (support_available) {
    par(mfrow = c(2, 1),
        bty = "n",
        mar = c(bottom.margin, left.margin, top.margin, right.margin))
    plot.support.panel(plot_support_values, ylab_text = support_label, axis_digits = round.axis.labels.BIC.plot)
    plot.frequency.panel()
  } else {
    par(mfrow = c(1, 1),
        bty = "n",
        mar = c(bottom.margin, left.margin, top.margin, right.margin))
    plot.frequency.panel()
    title(title, line = 0)
  }
  
  # Close graphics device
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}


## Function to plot model results as SOM grids (showing sample assignment to cells, cell distances and boundaries between cell clusters)
plot.model.SOM <- function(SOM.output,
                           col.pal.neighbor.dist = viridis::cividis, #color palette of neighbor distance plot (top)
                           col.pal.clusters = viridis::viridis, #color palette of cluster plot (bottom)
                           replicate.mode = "representative", #options: "first", "representative", "average"
                           set.k = NULL, #set to only use replicates with this k (NULL = no restriction)
                           save = F, #option to save plot
                           overwrite = T, #option to overwrite plot if already present (only if saving plot)
                           plot.type = "svg", #plot type options: "svg", "png", "jpg" (only if saving plot)
                           file.name = NULL, #set plot file name (if NULL, default name is used; only if saving plot)
                           width = 10, #plot width in cm (only if saving plot)
                           height = 15, #plot height in cm (only if saving plot)
                           resolution = 300, #plot resolution in dpi (only if saving plot)
                           bottom.margin = 0, #bottom margin
                           left.margin = 0, #left margin
                           top.margin = 1, #top margin
                           right.margin = 0, #right margin
                           boundary.col.clusters = "red", #color of cluster boundaries (in bottom plot)
                           boundary.lwd.clusters = 3, #line width of cluster boundaries (in bottom plot)
                           point.col.clusters = "white", #color of sample points (in bottom plot)
                           point.shape.clusters = 19, #shape of sample points (in bottom plot)
                           point.size.clusters = 0.9, #size of sample points (in bottom plot)
                           cluster.shape.clusters = "straight", #shape ("straight" or "round") of cluster cells (in bottom plot)
                           cluster.shape.neighbor.dist = "straight", #shape ("straight" or "round") of cluster cells (in top plot)
                           shift.plot.clusters = 0.099, #shift bottom plot slightly to the right to align with gridraster of top plot
                           title.clusters = "SOM clusters", #title of cluster plot (bottom)
                           title.neighbor.dist = "SOM neighbor distances" #title of neighbor distances plot
) {
  
  # Create function to calculate unit neighbor distances from SOM codebook vectors (ALL layers combined)
  calc.unit.neighbor.dist <- function(som_model) {
    codes <- kohonen::getCodes(som_model) #extract codes (matrix or list)
    if (is.null(codes)) stop("som_model has no codes")
    if (!is.list(codes)) codes <- list(codes)
    codebook_vectors <- do.call(cbind, codes) #combine all layers
    if (is.null(codebook_vectors) || nrow(codebook_vectors) < 1) stop("som_model codes are empty after combining layers")
    if (is.null(som_model$grid) || is.null(som_model$grid$pts)) stop("som_model has no grid points")
    grid_points <- som_model$grid$pts
    if (nrow(codebook_vectors) != nrow(grid_points)) stop("codes and grid points mismatch")
    
    # Replace invalid entries in combined codebook (NA/Inf) with column means
    codebook_vectors <- as.matrix(codebook_vectors)
    invalid <- !is.finite(codebook_vectors) | is.na(codebook_vectors)
    if (any(invalid)) {
      for (j in seq_len(ncol(codebook_vectors))) {
        col_vals <- codebook_vectors[, j]
        col_mean <- mean(col_vals[is.finite(col_vals)], na.rm = TRUE)
        if (!is.finite(col_mean)) col_mean <- 0.5
        col_vals[!is.finite(col_vals) | is.na(col_vals)] <- col_mean
        codebook_vectors[, j] <- col_vals
      }
    }
    
    grid_point_distance_matrix <- as.matrix(stats::dist(grid_points))
    grid_neighbor_step_distance <- min(grid_point_distance_matrix[grid_point_distance_matrix > 0])
    neighbor_index_pairs <- which(grid_point_distance_matrix > 0 & grid_point_distance_matrix <= (grid_neighbor_step_distance + 1e-12), arr.ind = TRUE)
    unit_mean_neighbor_distances <- rep(NA_real_, nrow(codebook_vectors))
    for (unit_index in seq_len(nrow(codebook_vectors))) {
      neighboring_unit_indices <- neighbor_index_pairs[neighbor_index_pairs[, 1] == unit_index, 2]
      if (length(neighboring_unit_indices) == 0) next
      unit_mean_neighbor_distances[unit_index] <- mean(sqrt(rowSums((codebook_vectors[neighboring_unit_indices, , drop = FALSE] - matrix(codebook_vectors[unit_index, ], nrow = length(neighboring_unit_indices), ncol = ncol(codebook_vectors), byrow = TRUE))^2)))
    }
    unit_mean_neighbor_distances
  }
  
  # Create function to extract combined codebook vectors from SOM (ALL layers combined)
  get.combined.codebook <- function(som_model) {
    codes <- kohonen::getCodes(som_model) #extract codes (matrix or list)
    if (is.null(codes)) stop("som_model has no codes")
    if (!is.list(codes)) codes <- list(codes)
    codebook_vectors <- do.call(cbind, codes) #combine all layers
    if (is.null(codebook_vectors) || nrow(codebook_vectors) < 1) stop("som_model codes are empty after combining layers")
    
    # Replace invalid entries in combined codebook (NA/Inf) with column means
    codebook_vectors <- as.matrix(codebook_vectors)
    invalid <- !is.finite(codebook_vectors) | is.na(codebook_vectors)
    if (any(invalid)) {
      for (j in seq_len(ncol(codebook_vectors))) {
        col_vals <- codebook_vectors[, j]
        col_mean <- mean(col_vals[is.finite(col_vals)], na.rm = TRUE)
        if (!is.finite(col_mean)) col_mean <- 0.5
        col_vals[!is.finite(col_vals) | is.na(col_vals)] <- col_mean
        codebook_vectors[, j] <- col_vals
      }
    }
    codebook_vectors
  }
  
  # Create function to align neuron positions across replicates using Hungarian algorithm on ALL-layer codebooks
  align.unit.positions <- function(reference_som_model, som_model_to_align) {
    reference_codebook_vectors <- get.combined.codebook(reference_som_model)
    aligned_codebook_vectors <- get.combined.codebook(som_model_to_align)
    if (nrow(reference_codebook_vectors) != nrow(aligned_codebook_vectors)) stop("Unit alignment aborted: codebook row mismatch across replicates")
    number_of_units <- nrow(reference_codebook_vectors)
    cost_matrix <- matrix(NA_real_, nrow = number_of_units, ncol = number_of_units)
    for (reference_unit_index in seq_len(number_of_units)) {
      reference_vector <- reference_codebook_vectors[reference_unit_index, ]
      diffs <- aligned_codebook_vectors - matrix(reference_vector, nrow = number_of_units, ncol = ncol(aligned_codebook_vectors), byrow = TRUE)
      cost_matrix[reference_unit_index, ] <- rowSums(diffs^2)
    }
    reference_to_aligned_map <- as.integer(clue::solve_LSAP(cost_matrix))
    aligned_to_reference_map <- rep(NA_integer_, number_of_units)
    aligned_to_reference_map[reference_to_aligned_map] <- seq_len(number_of_units)
    list(reference_to_aligned_map = reference_to_aligned_map,
         aligned_to_reference_map = aligned_to_reference_map)
  }
  
  # Create function to add boundaries of SOM clusters (works when kohonen::add.cluster.boundaries is unavailable)
  add.cluster.boundaries <- function(som_model, som_cluster, lwd = 3, col = "red") {
    if (is.null(som_model$grid) || is.null(som_model$grid$pts)) return(invisible(NULL))
    grid_points <- som_model$grid$pts
    if (length(som_cluster) != nrow(grid_points)) stop("Boundary plotting aborted: som_cluster length does not match number of grid units")
    grid_point_distance_matrix <- as.matrix(stats::dist(grid_points))
    grid_neighbor_step_distance <- min(grid_point_distance_matrix[grid_point_distance_matrix > 0])
    neighbor_index_pairs <- which(grid_point_distance_matrix > 0 & grid_point_distance_matrix <= (grid_neighbor_step_distance + 1e-12), arr.ind = TRUE)
    neighbor_index_pairs <- neighbor_index_pairs[neighbor_index_pairs[, 1] < neighbor_index_pairs[, 2], , drop = FALSE]
    if (nrow(neighbor_index_pairs) == 0) return(invisible(NULL))
    for (pair_index in seq_len(nrow(neighbor_index_pairs))) {
      unit_a <- neighbor_index_pairs[pair_index, 1]
      unit_b <- neighbor_index_pairs[pair_index, 2]
      if (som_cluster[unit_a] == som_cluster[unit_b]) next
      xa <- grid_points[unit_a, 1]
      ya <- grid_points[unit_a, 2]
      xb <- grid_points[unit_b, 1]
      yb <- grid_points[unit_b, 2]
      mx <- (xa + xb) / 2
      my <- (ya + yb) / 2
      vx <- xb - xa
      vy <- yb - ya
      d <- sqrt(vx^2 + vy^2)
      if (!is.finite(d) || d == 0) next
      px <- -vy / d
      py <- vx / d
      segment_length <- (d / sqrt(3)) * 0.98
      hx <- px * (segment_length / 2)
      hy <- py * (segment_length / 2)
      graphics::segments(mx - hx, my - hy, mx + hx, my + hy, col = col, lwd = lwd)
    }
    invisible(NULL)
  }
  
  # Create function to align cluster labels using Hungarian algorithm (deterministic relabeling)
  align.cluster.labels <- function(reference_cluster_labels, cluster_labels_to_align) {
    reference_cluster_labels <- as.integer(reference_cluster_labels)
    cluster_labels_to_align <- as.integer(cluster_labels_to_align)
    if (length(reference_cluster_labels) != length(cluster_labels_to_align)) stop("reference_cluster_labels and cluster_labels_to_align must be same length")
    if (anyNA(reference_cluster_labels) || anyNA(cluster_labels_to_align)) stop("NA labels not supported in alignment")
    number_of_reference_clusters <- max(reference_cluster_labels)
    number_of_clusters_to_align <- max(cluster_labels_to_align)
    if (!is.finite(number_of_reference_clusters) || !is.finite(number_of_clusters_to_align)) stop("Alignment aborted: invalid cluster label values")
    number_of_clusters <- max(number_of_reference_clusters, number_of_clusters_to_align)
    contingency_table <- table(factor(reference_cluster_labels, levels = seq_len(number_of_clusters)), factor(cluster_labels_to_align, levels = seq_len(number_of_clusters)))
    cost_matrix <- max(contingency_table) - contingency_table
    optimal_cluster_assignment <- clue::solve_LSAP(cost_matrix)
    relabel_mapping <- rep(NA_integer_, number_of_clusters)
    relabel_mapping[as.integer(optimal_cluster_assignment)] <- seq_len(number_of_clusters)
    relabel_mapping[cluster_labels_to_align]
  }
  
  # Create function to select representative SOM replicate based on mean pairwise Adjusted Rand Index (ARI)
  choose.representative.replicate <- function(som_clusters) {
    number_of_replicates <- length(som_clusters)
    if (number_of_replicates < 2) return(1)
    k_vals <- vapply(som_clusters, function(x) suppressWarnings(max(as.integer(x), na.rm = TRUE)), integer(1))
    mean_adjusted_rand_index_per_replicate <- rep(NA_real_, number_of_replicates)
    for (candidate_replicate_index in seq_len(number_of_replicates)) {
      candidate_k <- k_vals[candidate_replicate_index]
      if (!is.finite(candidate_k) || is.na(candidate_k)) next
      comparison_indices <- which(is.finite(k_vals) & !is.na(k_vals) & k_vals == candidate_k)
      comparison_indices <- comparison_indices[comparison_indices != candidate_replicate_index]
      if (length(comparison_indices) == 0) {
        mean_adjusted_rand_index_per_replicate[candidate_replicate_index] <- -Inf
        next
      }
      adjusted_rand_index_values <- rep(NA_real_, length(comparison_indices))
      for (i in seq_along(comparison_indices)) {
        comparison_replicate_index <- comparison_indices[i]
        candidate_replicate_unit_clusters <- as.integer(som_clusters[[candidate_replicate_index]])
        comparison_replicate_unit_clusters <- as.integer(som_clusters[[comparison_replicate_index]])
        if (length(candidate_replicate_unit_clusters) != length(comparison_replicate_unit_clusters)) stop("All som_clusters must have same length (same number of units)")
        comparison_replicate_unit_clusters_aligned <- align.cluster.labels(candidate_replicate_unit_clusters, comparison_replicate_unit_clusters)
        adjusted_rand_index_values[i] <- mclust::adjustedRandIndex(candidate_replicate_unit_clusters, comparison_replicate_unit_clusters_aligned)
      }
      mean_adjusted_rand_index_per_replicate[candidate_replicate_index] <- mean(adjusted_rand_index_values, na.rm = TRUE)
    }
    if (all(!is.finite(mean_adjusted_rand_index_per_replicate) | is.na(mean_adjusted_rand_index_per_replicate))) return(1)
    which.max(replace(mean_adjusted_rand_index_per_replicate, is.na(mean_adjusted_rand_index_per_replicate), -Inf))
  }
  
  # Create function to compute consensus SOM cluster labels across replicates with deterministic tie-breaking
  consensus.som.clusters <- function(som_clusters, reference_replicate_index = NULL) {
    number_of_replicates <- length(som_clusters)
    if (number_of_replicates == 0) stop("som_clusters is empty")
    if (number_of_replicates == 1) return(as.integer(som_clusters[[1]]))
    if (is.null(reference_replicate_index)) reference_replicate_index <- choose.representative.replicate(som_clusters)
    reference_unit_cluster_labels <- as.integer(som_clusters[[reference_replicate_index]])
    aligned_unit_cluster_labels_matrix <- matrix(NA_integer_, nrow = length(reference_unit_cluster_labels), ncol = number_of_replicates)
    aligned_unit_cluster_labels_matrix[, reference_replicate_index] <- reference_unit_cluster_labels
    for (replicate_index in seq_len(number_of_replicates)) {
      if (replicate_index == reference_replicate_index) next
      aligned_unit_cluster_labels_matrix[, replicate_index] <- align.cluster.labels(reference_unit_cluster_labels, as.integer(som_clusters[[replicate_index]]))
    }
    consensus_unit_cluster_labels <- rep(NA_integer_, nrow(aligned_unit_cluster_labels_matrix))
    for (unit_index in seq_len(nrow(aligned_unit_cluster_labels_matrix))) {
      unit_cluster_label_counts <- table(aligned_unit_cluster_labels_matrix[unit_index, ])
      most_frequent_cluster_labels <- as.integer(names(unit_cluster_label_counts)[unit_cluster_label_counts == max(unit_cluster_label_counts)])
      if (length(most_frequent_cluster_labels) == 1) {
        consensus_unit_cluster_labels[unit_index] <- most_frequent_cluster_labels
      } else if (reference_unit_cluster_labels[unit_index] %in% most_frequent_cluster_labels) {
        consensus_unit_cluster_labels[unit_index] <- reference_unit_cluster_labels[unit_index]
      } else {
        consensus_unit_cluster_labels[unit_index] <- min(most_frequent_cluster_labels)
      }
    }
    consensus_unit_cluster_labels
  }
  
  # Create function to compute consensus unit.classif across replicates with deterministic tie-breaking
  consensus.unit.classif <- function(unit_classif_matrix, preferred_unit_classif) {
    consensus_vec <- rep(NA_integer_, nrow(unit_classif_matrix))
    for (i in seq_len(nrow(unit_classif_matrix))) {
      vals <- as.integer(unit_classif_matrix[i, ])
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) next
      tab <- table(vals)
      top_vals <- as.integer(names(tab)[tab == max(tab)])
      if (length(top_vals) == 1) {
        consensus_vec[i] <- top_vals
      } else if (!is.na(preferred_unit_classif[i]) && preferred_unit_classif[i] %in% top_vals) {
        consensus_vec[i] <- preferred_unit_classif[i]
      } else {
        consensus_vec[i] <- min(top_vals)
      }
    }
    consensus_vec
  }
  
  # Reset plotting parameters
  old_dev <- dev.cur()
  old_plotting_parameters <- par(no.readonly = TRUE)
  on.exit({
    if (dev.cur() == old_dev) par(old_plotting_parameters)
  }, add = TRUE)
  
  # Validate SOM.output
  if (is.null(SOM.output$som_models) || is.null(SOM.output$som_clusters)) {
    stop("Plotting aborted: SOM.output is missing 'som_models' or 'som_clusters' - check SOM.output or rerun train.SOM/clustering.SOM")
  }
  
  # Validate replicate.mode
  allowed_rep_modes <- c("first", "representative", "average")
  if (!is.character(replicate.mode) || length(replicate.mode) != 1 || is.na(replicate.mode) || !(replicate.mode %in% allowed_rep_modes)) {
    stop("Plotting aborted: replicate.mode must be 'first', 'representative', or 'average'")
  }
  
  # Validate specified set.k
  if (!is.null(set.k) && (!is.numeric(set.k) || length(set.k) != 1 || is.na(set.k) || set.k < 1 || (set.k %% 1 != 0))) {
    stop("Plotting aborted: set.k must be NULL or single positive integer >= 1")
  }
  
  # Validate specified col.pal.neighbor.dist and col.pal.clusters
  viridis_palettes <- list(
    viridis::viridis,
    viridis::magma,
    viridis::plasma,
    viridis::inferno,
    viridis::cividis,
    viridis::rocket,
    viridis::mako,
    viridis::turbo
  )
  if (!any(vapply(viridis_palettes, identical, logical(1), col.pal.neighbor.dist))) {
    stop("Plotting aborted: col.pal.neighbor.dist must be a viridis palette - viridis, magma, plasma, inferno, cividis, rocket, mako or turbo")
  }
  if (!any(vapply(viridis_palettes, identical, logical(1), col.pal.clusters))) {
    stop("Plotting aborted: col.pal.clusters must be a viridis palette - viridis, magma, plasma, inferno, cividis, rocket, mako or turbo")
  }
  
  # Validate specified save
  if (!is.logical(save) || length(save) != 1 || is.na(save)) {
    stop("Plotting aborted: save must be TRUE or FALSE")
  }
  
  # Validate specified overwrite
  if (save && (!is.logical(overwrite) || length(overwrite) != 1 || is.na(overwrite))) {
    stop("Plotting aborted: overwrite must be TRUE or FALSE")
  }
  
  # Validate specified plot.type
  if (save) {
    allowed_types <- c("svg", "png", "jpg")
    if (!is.character(plot.type) || length(plot.type) != 1 || is.na(plot.type) || !(plot.type %in% allowed_types)) {
      stop("Plotting aborted: plot.type must be one of 'svg', 'png', or 'jpg'")
    }
  }
  
  # Validate specified file.name
  if (save && !is.null(file.name) && (!is.character(file.name) || length(file.name) != 1 || is.na(file.name))) {
    stop("Plotting aborted: file.name must be NULL or single character string")
  }
  
  # Validate specified width and height (reasonable values: 4–50 cm)
  if (save) {
    if (!is.numeric(width) || length(width) != 1 || is.na(width) || width <= 0) {
      stop("Plotting aborted: width must be a single positive number (cm)")
    }
    if (width < 4) {
      message("Warning: width is very small (", width, " cm) – plot may be hard to read")
    }
    if (width > 50) {
      message("Warning: width is very large (", width, " cm) – plot may be unwieldy")
    }
    if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
      stop("Plotting aborted: height must be a single positive number (cm)")
    }
    if (height < 4) {
      message("Warning: height is very small (", height, " cm) – plot may be hard to read")
    }
    if (height > 50) {
      message("Warning: height is very large (", height, " cm) – plot may be unwieldy")
    }
  }
  
  # Validate specified resolution (reasonable range: 72–1200 dpi)
  if (save) {
    if (!is.numeric(resolution) || length(resolution) != 1 || is.na(resolution) || resolution < 72) {
      stop("Plotting aborted: resolution must be a single number ≥ 72 (dpi)")
    }
    if (resolution > 1200) {
      message("Warning: resolution is very high (", resolution, " dpi) – file may be huge")
    }
  }
  
  # Validate specified margins (reasonable range: 0–10)
  margin.list <- c(bottom.margin, left.margin, top.margin, right.margin)
  margin.names <- c("bottom.margin", "left.margin", "top.margin", "right.margin")
  for (i in seq_along(margin.list)) {
    if (!is.numeric(margin.list[i]) || length(margin.list[i]) != 1 || is.na(margin.list[i])) {
      stop("Plotting aborted: ", margin.names[i], " must be a single numeric value")
    }
    if (margin.list[i] < 0) {
      stop("Plotting aborted: ", margin.names[i], " must be ≥ 0")
    }
    if (margin.list[i] > 10) {
      message("Warning: ", margin.names[i], " is large (", margin.list[i], ") – plot area may shrink")
    }
  }
  
  # Validate specified boundary.lwd.clusters
  if (!is.numeric(boundary.lwd.clusters) || length(boundary.lwd.clusters) != 1 || is.na(boundary.lwd.clusters) || boundary.lwd.clusters <= 0) {
    stop("Plotting aborted: boundary.lwd.clusters must be a single positive number")
  }
  
  # Validate specified point.size.clusters
  if (!is.numeric(point.size.clusters) || length(point.size.clusters) != 1 || is.na(point.size.clusters) || point.size.clusters <= 0) {
    stop("Plotting aborted: point.size.clusters must be a single positive number")
  }
  
  # Validate specified point.shape.clusters (should be integer)
  if (!is.numeric(point.shape.clusters) || length(point.shape.clusters) != 1 || is.na(point.shape.clusters) || (point.shape.clusters %% 1 != 0)) {
    stop("Plotting aborted: point.shape.clusters must be a single integer")
  }
  
  # Validate specified cluster.shape.clusters and cluster.shape.neighbor.dist
  allowed_shapes <- c("straight", "round")
  if (!is.character(cluster.shape.clusters) || length(cluster.shape.clusters) != 1 || is.na(cluster.shape.clusters) || !(cluster.shape.clusters %in% allowed_shapes)) {
    stop("Plotting aborted: cluster.shape.clusters must be 'straight' or 'round'")
  }
  if (!is.character(cluster.shape.neighbor.dist) || length(cluster.shape.neighbor.dist) != 1 || is.na(cluster.shape.neighbor.dist) || !(cluster.shape.neighbor.dist %in% allowed_shapes)) {
    stop("Plotting aborted: cluster.shape.neighbor.dist must be 'straight' or 'round'")
  }
  
  # Validate specified shift.plot.clusters (should be >= 0 and < 0.5)
  if (!is.numeric(shift.plot.clusters) || length(shift.plot.clusters) != 1 || is.na(shift.plot.clusters) || shift.plot.clusters < 0 || shift.plot.clusters >= 0.5) {
    message("Invalid shift.plot.clusters value (needs to be 0 - 0.5) - default of 0.099 is used")
    shift.plot.clusters <- 0.099
  }
  
  # Validate specified boundary.col.clusters and point.col.clusters (character)
  if (!is.character(boundary.col.clusters) || length(boundary.col.clusters) != 1 || is.na(boundary.col.clusters)) {
    stop("Plotting aborted: boundary.col.clusters must be a single character (color name or hex)")
  }
  if (!is.character(point.col.clusters) || length(point.col.clusters) != 1 || is.na(point.col.clusters)) {
    stop("Plotting aborted: point.col.clusters must be a single character (color name or hex)")
  }
  
  # Validate specified title.clusters and title.neighbor.dist (NULL or character)
  if (!is.null(title.clusters) && (!is.character(title.clusters) || length(title.clusters) != 1 || is.na(title.clusters))) {
    stop("Plotting aborted: title.clusters must be NULL or single character string")
  }
  if (!is.null(title.neighbor.dist) && (!is.character(title.neighbor.dist) || length(title.neighbor.dist) != 1 || is.na(title.neighbor.dist))) {
    stop("Plotting aborted: title.neighbor.dist must be NULL or single character string")
  }
  
  # Subset replicates by set.k (if provided)
  som_models_use <- SOM.output$som_models
  som_clusters_use <- SOM.output$som_clusters
  if (!is.null(set.k)) {
    k_vals <- vapply(som_clusters_use, function(x) suppressWarnings(max(as.integer(x), na.rm = TRUE)), integer(1))
    keep <- which(is.finite(k_vals) & !is.na(k_vals) & k_vals == as.integer(set.k))
    if (length(keep) == 0) stop("Plotting aborted: no replicates match set.k - rerun clustering or choose different set.k")
    som_models_use <- som_models_use[keep]
    som_clusters_use <- som_clusters_use[keep]
  }
  
  # Choose replicate index
  replicate.index <- if (replicate.mode == "first") 1 else choose.representative.replicate(som_clusters_use)
  som_model <- som_models_use[[replicate.index]]
  
  # Build plotting clusters + optionally averaged mapping
  if (replicate.mode == "average") {
    rep_k <- suppressWarnings(max(as.integer(som_clusters_use[[replicate.index]]), na.rm = TRUE))
    if (!is.finite(rep_k) || is.na(rep_k)) stop("Plotting aborted: representative replicate has invalid k")
    k_vals <- vapply(som_clusters_use, function(x) suppressWarnings(max(as.integer(x), na.rm = TRUE)), integer(1))
    keep_k <- which(is.finite(k_vals) & !is.na(k_vals) & k_vals == rep_k)
    if (length(keep_k) == 0) stop("Plotting aborted: no replicates share k with representative replicate")
    som_models_k <- som_models_use[keep_k]
    som_clusters_k <- som_clusters_use[keep_k]
    rep_index_k <- which(keep_k == replicate.index)
    if (length(rep_index_k) != 1) rep_index_k <- 1
    
    # Align neuron positions across replicates to representative replicate (ALL layers)
    reference_som_model <- som_models_k[[rep_index_k]]
    reference_unit_classif <- as.integer(reference_som_model$unit.classif)
    reference_unit_cluster_labels <- as.integer(som_clusters_k[[rep_index_k]])
    number_of_units <- length(reference_unit_cluster_labels)
    number_of_replicates_k <- length(som_models_k)
    
    unit_classif_aligned_matrix <- matrix(NA_integer_, nrow = length(reference_unit_classif), ncol = number_of_replicates_k)
    som_clusters_aligned <- vector("list", number_of_replicates_k)
    neighbor_distances_aligned_matrix <- matrix(NA_real_, nrow = number_of_units, ncol = number_of_replicates_k)
    
    for (replicate_index in seq_len(number_of_replicates_k)) {
      if (replicate_index == rep_index_k) {
        reference_to_aligned_map <- seq_len(number_of_units)
        aligned_to_reference_map <- seq_len(number_of_units)
      } else {
        maps <- align.unit.positions(reference_som_model, som_models_k[[replicate_index]])
        reference_to_aligned_map <- maps$reference_to_aligned_map
        aligned_to_reference_map <- maps$aligned_to_reference_map
      }
      
      unit_classif_vec <- as.integer(som_models_k[[replicate_index]]$unit.classif)
      if (anyNA(unit_classif_vec)) stop("Plotting aborted: unit.classif contains NA in a replicate")
      unit_classif_aligned_matrix[, replicate_index] <- aligned_to_reference_map[unit_classif_vec]
      
      som_clusters_aligned[[replicate_index]] <- as.integer(som_clusters_k[[replicate_index]])[reference_to_aligned_map]
      
      neighbor_distances_vec <- calc.unit.neighbor.dist(som_models_k[[replicate_index]])
      neighbor_distances_aligned_matrix[, replicate_index] <- neighbor_distances_vec[reference_to_aligned_map]
    }
    
    # Consensus sample->neuron assignment (unit.classif) in reference unit space
    preferred_unit_assignments <- unit_classif_aligned_matrix[, rep_index_k]
    som_model$unit.classif <- consensus.unit.classif(unit_classif_aligned_matrix, preferred_unit_assignments)
    
    # Consensus neuron clusters (reference replicate’s neuron cluster labels as reference, Hungarian-align labels, then majority vote)
    som_cluster <- consensus.som.clusters(som_clusters_aligned, reference_replicate_index = rep_index_k)
    
    # Average neighbor distances per neuron (aligned to reference unit space)
    nd_plot <- rowMeans(neighbor_distances_aligned_matrix, na.rm = TRUE)
    
  } else {
    som_cluster <- som_clusters_use[[replicate.index]]
    nd_plot <- calc.unit.neighbor.dist(som_model)
  }
  
  # Save file
  if (save) {
    
    # Set default file name
    if (is.null(file.name)) {
      file.name <- paste0("SOM_model_plot_", paste(SOM.output$input_data_names, collapse = "_"), ".", plot.type)
    }
    
    # Set overwrite option
    if (file.exists(file.name) && !overwrite) {
      stop(paste("File", file.name, "already exists - set overwrite = T to overwrite"))
    }
    
    # Set plot format
    if (plot.type == "svg") {
      svg(file.name, 
          width = width / 2.54, 
          height = height / 2.54)
    } else if (plot.type == "png") {
      png(file.name, 
          width = width, 
          height = height, 
          res = resolution, 
          units = "cm")
    } else if (plot.type == "jpg") {
      jpeg(file.name, 
           width = width,
           height = height, 
           res = resolution, 
           units = "cm")
    }
  }
  
  # Set plotting area
  par(mfrow = c(2, 1), 
      oma = c(0, 0, 0, 0),
      mar = c(bottom.margin, 
              left.margin, 
              top.margin, 
              right.margin))
  
  # Plot SOM neighbor distances (top plot) - ALWAYS use ALL layers via property vector
  plot(x = som_model,
       type = "property",
       property = nd_plot,
       main = title.neighbor.dist,
       shape = cluster.shape.neighbor.dist,
       palette.name = function(n) rev(col.pal.neighbor.dist(n)))
  
  # Set color palette for bottom plot
  som_cluster <- as.integer(som_cluster)
  if (anyNA(som_cluster) || !all(is.finite(som_cluster))) stop("Plotting aborted: som_cluster contains NA/Inf")
  k.cols <- col.pal.clusters(max(som_cluster))
  SOM_cluster_plot_col <- rep(NA, length(som_cluster))
  for (i in seq_len(max(som_cluster))) {
    SOM_cluster_plot_col[som_cluster == i] <- k.cols[i]
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


## Function to plot sample map with cluster assignment for each individual
plot.map.SOM <- function(SOM.output,
                         Coordinates, #coordinates as "Latitude" and "Longitude" columns in dataframe/matrix
                         save = F, #whether to save plot or not
                         overwrite = T, #whether to overwrite if file exists (only if saving plot)
                         plot.type = "svg", #plot format type options: "png", "svg", "jpg" (only if saving plot)
                         file.name = NULL, #plot file name (NULL = default file name) (only if saving plot)
                         width = 15, #plot width in cm (only if saving plot)
                         height = 20, #plot height in cm (only if saving plot)
                         resolution = 300, #plot resolution in dpi (only if saving plot)
                         lat.buffer.range = 2, #add coordinates as buffer range around latitude coordinates
                         lon.buffer.range = 2, #add coordinates as buffer range around longitude coordinates
                         pie.size = 2, #pie chart size
                         pie.col.pal = viridis::viridis, #color palette of pie charts
                         USA.add.states = T, #option to add US states (only works if range includes USA)
                         USA.add.counties = F, #option to add US counties (only works if range includes USA)
                         USA.state.lwd = 0.5, #linewidth of US state borders (only works if range includes USA)
                         USA.county.lwd = 0.5, #linewidth of US county borders (only works if range includes USA)
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
                         legend.symbol.size = 1.5 #size of legend symbols
) {
  
  # Reset plotting parameters
  old_dev <- grDevices::dev.cur()
  old_plotting_parameters <- graphics::par(no.readonly = TRUE)
  on.exit({
    if (grDevices::dev.cur() == old_dev) graphics::par(old_plotting_parameters)
  }, add = TRUE)
  
  # Validate SOM.output
  if (is.null(SOM.output$ancestry_matrix) || !is.matrix(SOM.output$ancestry_matrix)) {
    stop("Plotting aborted: ancestry_matrix of SOM.output is not valid")
  }
  
  # Validate Coordinates data
  if (!all(c("Latitude", "Longitude") %in% names(Coordinates))) {
    stop("Plotting aborted: Coordinates must contain 'Latitude' and 'Longitude' columns")
  }
  if (is.null(rownames(Coordinates))) {
    stop("Plotting aborted: Coordinates must have rownames matching rownames of ancestry_matrix")
  }
  if (!is.data.frame(Coordinates) && !is.matrix(Coordinates)) {
    stop("Plotting aborted: Coordinates must be a data frame or matrix")
  }
  
  # Validate specified color palette
  viridis_palettes <- list(viridis::viridis,
                           viridis::magma,
                           viridis::plasma,
                           viridis::inferno,
                           viridis::cividis,
                           viridis::rocket,
                           viridis::mako,
                           viridis::turbo
  )
  if (!any(vapply(viridis_palettes, identical, logical(1), pie.col.pal))) {
    stop("Plotting aborted: pie.col.pal must viridis palette - viridis, magma, plasma, inferno, cividis, rocket, mako or turbo")
  }
  
  # Validate specified save
  if (!is.logical(save) || length(save) != 1 || is.na(save)) {
    stop("Plotting aborted: save must be TRUE or FALSE")
  }
  
  # Validate specified overwrite
  if (!is.logical(overwrite) || length(overwrite) != 1 || is.na(overwrite)) {
    stop("Plotting aborted: overwrite must be TRUE or FALSE")
  }
  
  # Validate specified plot.type
  if (save) {
    allowed_plot.types <- c("svg", "png", "jpg")
    if (!is.character(plot.type) || length(plot.type) != 1 || is.na(plot.type) || !(plot.type %in% allowed_plot.types)) {
      stop("Plotting aborted: plot.type must be one of 'svg', 'png', or 'jpg'")
    }
  }
  
  # Validate specified file.name
  if (save && !is.null(file.name) && (!is.character(file.name) || length(file.name) != 1 || is.na(file.name))) {
    stop("Plotting aborted: file.name must be NULL or single character string")
  }
  
  # Validate specified width and height (reasonable values: 4–50 cm)
  if (save) {
    if (!is.numeric(width) || length(width) != 1 || is.na(width) || width <= 0) {
      stop("Plotting aborted: width must be a single positive number (cm)")
    }
    if (width < 4) {
      message("Warning: width is very small (", width, " cm) – plot may be hard to read")
    }
    if (width > 50) {
      message("Warning: width is very large (", width, " cm) – plot may be unwieldy")
    }
    if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
      stop("Plotting aborted: height must be a single positive number (cm)")
    }
    if (height < 4) {
      message("Warning: height is very small (", height, " cm) – plot may be hard to read")
    }
    if (height > 50) {
      message("Warning: height is very large (", height, " cm) – plot may be unwieldy")
    }
  }
  
  # Validate specified resolution (reasonable range: 72–1200 dpi)
  if (save) {
    if (!is.numeric(resolution) || length(resolution) != 1 || is.na(resolution) || resolution < 72) {
      stop("Plotting aborted: resolution must be a single number ≥ 72 (dpi)")
    }
    if (resolution > 1200) {
      message("Warning: resolution is very high (", resolution, " dpi) – file may be huge")
    }
  }
  
  # Validate specified lat.buffer.range
  if (!is.numeric(lat.buffer.range) || length(lat.buffer.range) != 1 || is.na(lat.buffer.range) || lat.buffer.range < 0) {
    stop("Plotting aborted: lat.buffer.range must be a single non-negative number")
  }
  
  # Validate specified lon.buffer.range
  if (!is.numeric(lon.buffer.range) || length(lon.buffer.range) != 1 || is.na(lon.buffer.range) || lon.buffer.range < 0) {
    stop("Plotting aborted: lon.buffer.range must be a single non-negative number")
  }
  
  # Validate specified pie.size
  if (!is.numeric(pie.size) || length(pie.size) != 1 || is.na(pie.size) || pie.size <= 0) {
    stop("Plotting aborted: pie.size must be a single positive number")
  }
  
  # Validate specified pie.col.pal
  viridis_palettes <- list(viridis::viridis, viridis::magma, viridis::plasma, viridis::inferno,
                           viridis::cividis, viridis::rocket, viridis::mako, viridis::turbo)
  if (!any(vapply(viridis_palettes, identical, logical(1), pie.col.pal))) {
    stop("Plotting aborted: pie.col.pal must be viridis palette - viridis, magma, plasma, inferno, cividis, rocket, mako or turbo")
  }
  
  # Validate specified USA.add.states
  if (!is.logical(USA.add.states) || length(USA.add.states) != 1 || is.na(USA.add.states)) {
    stop("Plotting aborted: USA.add.states must be TRUE or FALSE")
  }
  
  # Validate specified USA.add.counties
  if (!is.logical(USA.add.counties) || length(USA.add.counties) != 1 || is.na(USA.add.counties)) {
    stop("Plotting aborted: USA.add.counties must be TRUE or FALSE")
  }
  
  # Validate specified USA.state.lwd
  if (!is.numeric(USA.state.lwd) || length(USA.state.lwd) != 1 || is.na(USA.state.lwd) || USA.state.lwd <= 0) {
    stop("Plotting aborted: USA.state.lwd must be a single positive number")
  }
  
  # Validate specified USA.county.lwd
  if (!is.numeric(USA.county.lwd) || length(USA.county.lwd) != 1 || is.na(USA.county.lwd) || USA.county.lwd <= 0) {
    stop("Plotting aborted: USA.county.lwd must be a single positive number")
  }
  # Validate specified north.arrow.position
  if (!is.numeric(north.arrow.position) || length(north.arrow.position) != 2 || any(is.na(north.arrow.position)) ||
      any(north.arrow.position < 0) || any(north.arrow.position > 1)) {
    stop("Plotting aborted: north.arrow.position must be numeric vector of length 2 with values in between 0 and 1")
  }
  
  # Validate specified north.arrow.length
  if (!is.numeric(north.arrow.length) || length(north.arrow.length) != 1 || is.na(north.arrow.length) || north.arrow.length <= 0) {
    stop("Plotting aborted: north.arrow.length must be a single positive number")
  }
  
  # Validate specified north.arrow.lwd
  if (!is.numeric(north.arrow.lwd) || length(north.arrow.lwd) != 1 || is.na(north.arrow.lwd) || north.arrow.lwd <= 0) {
    stop("Plotting aborted: north.arrow.lwd must be a single positive number")
  }
  
  # Validate specified north.arrow.N.position
  if (!is.numeric(north.arrow.N.position) || length(north.arrow.N.position) != 1 || is.na(north.arrow.N.position) || north.arrow.N.position < 0) {
    stop("Plotting aborted: north.arrow.N.position must be a single non-negative number")
  }
  
  # Validate specified north.arrow.N.size
  if (!is.numeric(north.arrow.N.size) || length(north.arrow.N.size) != 1 || is.na(north.arrow.N.size) || north.arrow.N.size <= 0) {
    stop("Plotting aborted: north.arrow.N.size must be a single positive number")
  }
  
  # Validate specified scale.position
  if (!is.numeric(scale.position) || length(scale.position) != 2 || any(is.na(scale.position)) ||
      any(scale.position < 0) || any(scale.position > 1)) {
    stop("Plotting aborted: scale.position must be numeric vector of length 2 with values in between 0 and 1")
  }
  
  # Validate specified scale.size
  if (!is.numeric(scale.size) || length(scale.size) != 1 || is.na(scale.size) || scale.size <= 0) {
    stop("Plotting aborted: scale.size must be a single positive number")
  }
  
  # Validate specified scale.font.size
  if (!is.numeric(scale.font.size) || length(scale.font.size) != 1 || is.na(scale.font.size) || scale.font.size <= 0) {
    stop("Plotting aborted: scale.font.size must be a single positive number")
  }
  
  # Validate specified legend.position
  allowed.legend.positions <- c("topright", "topleft", "bottomright", "bottomleft", 
                                "right", "left", "top", "bottom", "center")
  if (!is.character(legend.position) || length(legend.position) != 1 || is.na(legend.position) ||
      !(legend.position %in% allowed.legend.positions)) {
    stop(paste0("Plotting aborted: legend.position must be one of ", 
                paste(allowed.legend.positions, collapse = ", ")))
  }
  
  # Validate specified legend.cluster.names
  if (!is.null(legend.cluster.names)) {
    if (!is.character(legend.cluster.names) || any(is.na(legend.cluster.names))) {
      stop("Plotting aborted: legend.cluster.names must be NULL or character vector (no NAs)")
    }
    n.clusters <- ncol(SOM.output$ancestry_matrix)
    if (length(legend.cluster.names) != n.clusters) {
      stop(paste0("Plotting aborted: length of legend.cluster.names (", length(legend.cluster.names), 
                  ") must match number of clusters (", n.clusters, ")"))
    }
  }
  
  # Validate specified legend.font.size
  if (!is.numeric(legend.font.size) || length(legend.font.size) != 1 || is.na(legend.font.size) || legend.font.size <= 0) {
    stop("Plotting aborted: legend.font.size must be a single positive number")
  }
  
  # Validate specified legend.box
  if (!is.logical(legend.box) || length(legend.box) != 1 || is.na(legend.box)) {
    stop("Plotting aborted: legend.box must be TRUE or FALSE")
  }
  
  # Validate specified legend.symbol.size
  if (!is.numeric(legend.symbol.size) || length(legend.symbol.size) != 1 || is.na(legend.symbol.size) || legend.symbol.size <= 0) {
    stop("Plotting aborted: legend.symbol.size must be a single positive number")
  }
  
  # Convert matrix to data frame if necessary
  if (is.matrix(Coordinates)) {
    Coordinates <- as.data.frame(Coordinates)
  }
  
  # Check rownames match ancestry matrix, try to reorder, remove non-matching
  coord_names <- rownames(Coordinates)
  sample_names <- rownames(SOM.output$ancestry_matrix)
  not_in_coords <- setdiff(sample_names, coord_names)
  not_in_ancestry <- setdiff(coord_names, sample_names)
  keep_names <- intersect(sample_names, coord_names)
  
  fmt_count <- function(n) if (n == 0) "No" else as.character(n)
  fmt_label <- function(n, singular, plural) if (n == 1) singular else plural
  n_not_in_coords <- length(not_in_coords)
  n_not_in_ancestry <- length(not_in_ancestry)
  n_keep <- length(keep_names)
  if (n_not_in_coords > 0 | n_not_in_ancestry > 0) {
    message(sprintf(
      "Matching samples between ancestry matrix and coordinates:\n  - %s unique %s only in ancestry_matrix\n  - %s unique %s only in Coordinates\n  - %s matching %s will be plotted",
      fmt_count(n_not_in_coords), fmt_label(n_not_in_coords, "sample", "samples"),
      fmt_count(n_not_in_ancestry), fmt_label(n_not_in_ancestry, "sample", "samples"),
      fmt_count(n_keep), fmt_label(n_keep, "sample", "samples")
    ))
  }
  if (length(keep_names) == 0) {
    stop("Plotting aborted: no matching rownames between Coordinates and ancestry_matrix")
  }
  Coordinates <- Coordinates[keep_names, , drop = FALSE]
  ancestry <- SOM.output$ancestry_matrix[keep_names, , drop = FALSE]
  
  # Remove rows with NA in Coordinates
  na_rows <- which(is.na(Coordinates$Latitude) | is.na(Coordinates$Longitude))
  if (length(na_rows) > 0) {
    message(sprintf("Dropped %d of %d rows due to NA in Coordinates", length(na_rows), nrow(Coordinates)))
    Coordinates <- Coordinates[-na_rows, , drop = FALSE]
    ancestry <- ancestry[-na_rows, , drop = FALSE]
  }
  
  # Prepare ancestry matrix
  q_matrix = as.data.frame(ancestry) #convert ancestry_matrix to dataframe
  
  # Check if number of rows (samples) in ancestry matrix matches number of Coordinates
  if (nrow(q_matrix) != nrow(Coordinates)) {
    stop("Plotting aborted: number of samples in ancestry_matrix does not match number of samples in Coordinates")
  }
  
  # Define color palette for pie charts
  k.cols = pie.col.pal(ncol(ancestry))
  
  # Set pie plot function based on make.admix.pie.plot function
  plot.admixture.pies <- function (admix.proportions, 
                                   coords, 
                                   layer.colors = NULL, 
                                   radii = 2.7, add = FALSE, 
                                   x.lim = NULL, 
                                   y.lim = NULL, 
                                   mar = c(2, 2, 2, 2)) 
  {
    K <- ncol(admix.proportions)
    N <- nrow(admix.proportions)
    layer.names <- paste0("layer_", 1:K)
    sample.names <- paste0("sample_", 1:N)
    color.tab <- caroline::nv(c(layer.colors[1:K]), layer.names)
    pie.list <- lapply(1:N, function(i) {
      caroline::nv(admix.proportions[i, ], layer.names)
    })
    names(pie.list) <- sample.names
    if (add) {
      graphics::par(new = TRUE)
    }
    else {
      graphics::par(mar = mar)
    }
    if (is.null(x.lim)) {
      x.lim <- c(min(coords[, 1]) - 1, max(coords[, 1]) + 1)
    }
    if (is.null(y.lim)) {
      y.lim <- c(min(coords[, 2]) - 1, max(coords[, 2]) + 1)
    }
    suppressWarnings(caroline::pies(pie.list, 
                                    x0 = coords[, 1], 
                                    y0 = coords[, 2], 
                                    color.table = color.tab, 
                                    border = "black", 
                                    radii = radii, 
                                    xlab = "", 
                                    ylab = "", 
                                    main = "", 
                                    lty = 1, 
                                    density = NULL, 
                                    xlim = x.lim, 
                                    ylim = y.lim))
    return(invisible(0))
  }
  
  # Define map boundaries
  lon_min = min(Coordinates$Longitude) - lon.buffer.range
  lon_max = max(Coordinates$Longitude) + lon.buffer.range
  lat_min = min(Coordinates$Latitude) - lat.buffer.range
  lat_max = max(Coordinates$Latitude) + lat.buffer.range
  
  # Create plot
  if (grDevices::dev.cur() > 1) {
    grDevices::dev.off() #close current device if open to avoid unwanted graphic distortions and other effects
  }
  
  # Set plot saving
  if (save) {
    # Set file name for saving
    if (is.null(file.name)) {
      file.name <- paste0("SOM_map_plot_", paste(SOM.output$input_data_names, collapse = "_"), ".", plot.type)
    }
    # Set overwriting 
    if (file.exists(file.name) && !overwrite) {
      stop(paste(file.name, "already exists - set overwrite = T to overwrite"))
    }
    # Set file plot.type
    if (plot.type == "svg") {
      grDevices::svg(file.name, width = width / 2.54, height = height / 2.54)
    } else if (plot.type == "png") {
      grDevices::png(file.name, width = width, height = height, res = resolution, units = "cm")
    } else if (plot.type == "jpg") {
      grDevices::jpeg(file.name, width = width, height = height, res = resolution, units = "cm")
    } else {
      stop("Plotting aborted: unsupported plot file plot.type - choose from 'svg', 'png', or 'jpg'")
    }
  }
  
  # Set layout and margins
  graphics::par(mfrow = c(1, 1), 
                oma = c(0, 0, 0, 0),
                mar = c(1, 1, 1, 1))
  
  # Create map
  maps::map("world", 
            fill = T, 
            col = "lightgrey", 
            xlim = c(lon_min, lon_max), 
            ylim = c(lat_min, lat_max))
  maps::map.axes()
  
  # Add US counties if requested
  if (USA.add.counties) {
    maps::map("county", 
              add = T, 
              col = "grey", 
              lwd = USA.county.lwd)
  }
  
  # Add US states if requested
  if (USA.add.states) {
    maps::map("state", 
              add = T, 
              col = "black", 
              lwd = USA.state.lwd)
  }
  
  # Add pie charts to map
  for (i in 1:nrow(q_matrix)) {
    coords = matrix(c(Coordinates$Longitude[i], Coordinates$Latitude[i]), ncol = 2, byrow = TRUE)
    plot.admixture.pies(
      admix.proportions = matrix(as.numeric(q_matrix[i, ]), nrow = 1),
      coords = coords,
      layer.colors = k.cols,
      radii = pie.size,
      add = TRUE
    )
  }
  
  # Define legend labels
  if (is.null(legend.cluster.names)) {
    legend.labels <- paste("Cluster", 1:length(k.cols)) # set default labels
  } else {
    legend.labels <- legend.cluster.names # use validated custom labels
  }
  
  # Set legend box
  if (legend.box) {
    legend_box <- "o"
  } else {
    legend_box <- "n"  
  }
  
  # Add legend
  graphics::legend(legend.position, 
                   legend = legend.labels,
                   pch = 21,
                   cex = legend.font.size,
                   pt.cex = legend.symbol.size,
                   pt.bg = k.cols,
                   bty = legend_box)
  
  # Add scale
  scale_position_x = scale.position[1] * (lon_max - lon_min) + lon_min
  scale_position_y = scale.position[2] * (lat_max - lat_min) + lat_min
  maps::map.scale(x = scale_position_x,
                  y = scale_position_y,
                  cex = scale.font.size,
                  relwidth = scale.size,
                  ratio = F)
  
  # Add north arrow
  north_arrow_x = north.arrow.position[1] * (lon_max - lon_min) + lon_min
  north_arrow_y = north.arrow.position[2] * (lat_max - lat_min) + lat_min
  graphics::arrows(x0 = north_arrow_x, 
                   y0 = north_arrow_y, 
                   x1 = north_arrow_x, 
                   y1 = north_arrow_y + north.arrow.length, 
                   length = 0.13, 
                   col = "black", 
                   lwd = north.arrow.lwd)
  
  # Add North "N" above north arrow
  graphics::text(x = north_arrow_x, 
                 y = north_arrow_y + north.arrow.length + north.arrow.N.position, #adjust position for "N"
                 labels = "N", 
                 cex = north.arrow.N.size, 
                 col = "black")
  
  # Close graphics device
  if (save) {
    grDevices::dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}


## Function to plot variable importance for each SOM layer (eta squared cluster separation) or map variance (codebook vectors/neuron weights)
plot.variable.importance.SOM <- function(SOM.output, 
                                         mode = "Cluster.separation", #mode (options: "Cluster.separation", "Map.variance")
                                         col.pal = viridis::turbo, #color palette
                                         save = FALSE, #option to save plot
                                         overwrite = TRUE, #option to overwrite plot if it already exists
                                         plot.type = "svg", #plot file plot.type (options: "png", "svg", "jpg")
                                         file.name = NULL, #plot file name (if NULL, default file name is used)
                                         width = 20, #plot width in cm
                                         height = 15, #plot height in cm
                                         resolution = 300, #plot resolution in dpi
                                         bottom.margin.total = 1, #bottom margin of entire plot
                                         left.margin.total = 1, #left margin of entire plot
                                         top.margin.total = 2, #top margin of entire plot
                                         right.margin.total = 0, #right margin of entire plot
                                         bottom.margin = 2.5, #bottom margin of individual plots
                                         left.margin = 4.5, #left margin of individual plots
                                         top.margin = 2, #top margin of individual plots
                                         right.margin = 2, #right margin of individual plots
                                         bars.threshold.N = 50, #threshold for leaving out bar labels
                                         title.font.size = 1.2, #font size of title
                                         matrix.label.font.size = 1, #font size of matrix label(s)
                                         bar.label.font.size = 0.65, #font size of bar labels
                                         add.boxplot.whiskers = TRUE, #whether to plot boxplot whiskers
                                         importance.threshold = 0.001, #threshold for showing variable importance
                                         set.k = NULL #if NULL, include all replicates - if integer, include only replicates where number of clusters (K) equals set.k (only if mode = "Cluster.separation")
) {
  
  # Reset plotting parameters
  old_plotting_parameters <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_plotting_parameters), add = TRUE)
  
  # Validate SOM.output
  if (is.null(SOM.output$som_models) || length(SOM.output$som_models) < 1) {
    stop("Plotting aborted: som_models not found in SOM.output - run train.SOM()")
  }
  if (is.null(SOM.output$input_data_names)) {
    stop("Plotting aborted: input_data_names not found in SOM.output - check SOM.output or rerun train.SOM()")
  }
  if (!is.character(mode) || length(mode) != 1 || is.na(mode) || !(mode %in% c("Cluster.separation", "Map.variance"))) {
    stop("Plotting aborted: mode must be 'Cluster.separation' or 'Map.variance'")
  }
  if (mode == "Cluster.separation") {
    if (is.null(SOM.output$som_clusters) || length(SOM.output$som_clusters) != length(SOM.output$som_models)) {
      stop("Plotting aborted: som_clusters not found in SOM.output or length mismatch - run clustering.SOM()")
    }
  }
  
  # Validate specified col.pal
  viridis_palettes <- list(viridis::viridis, viridis::magma, viridis::plasma,
                           viridis::inferno, viridis::cividis, viridis::rocket,
                           viridis::mako, viridis::turbo)
  if (!any(vapply(viridis_palettes, identical, logical(1), col.pal))) {
    stop("Plotting aborted: col.pal must be viridis palette - viridis, magma, plasma, inferno, cividis, rocket, mako or turbo")
  }
  
  # Validate specified save
  if (!is.logical(save) || length(save) != 1 || is.na(save)) {
    stop("Plotting aborted: save must be TRUE or FALSE")
  }
  
  # Validate specified overwrite
  if (!is.logical(overwrite) || length(overwrite) != 1 || is.na(overwrite)) {
    stop("Plotting aborted: overwrite must be TRUE or FALSE")
  }
  
  # Validate specified plot.type
  if (save) {
    allowed_plot.types <- c("svg", "png", "jpg")
    if (!is.character(plot.type) || length(plot.type) != 1 || is.na(plot.type) || !(plot.type %in% allowed_plot.types)) {
      stop("Plotting aborted: plot.type must be one of 'svg', 'png', or 'jpg'")
    }
  }
  
  # Validate specified file.name
  if (save && !is.null(file.name) && (!is.character(file.name) || length(file.name) != 1 || is.na(file.name))) {
    stop("Plotting aborted: file.name must be NULL or single character string")
  }
  
  # Validate specified width and height
  if (save) {
    if (!is.numeric(width) || length(width) != 1 || is.na(width) || width <= 0) {
      stop("Plotting aborted: width must be a single positive number (cm)")
    }
    if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
      stop("Plotting aborted: height must be a single positive number (cm)")
    }
  }
  
  # Validate specified resolution
  if (save) {
    if (!is.numeric(resolution) || length(resolution) != 1 || is.na(resolution) || resolution < 72) {
      stop("Plotting aborted: resolution must be a single number ≥ 72 (dpi)")
    }
  }
  
  # Validate specified margins
  margin.list <- c(bottom.margin.total, left.margin.total, top.margin.total, right.margin.total,
                   bottom.margin, left.margin, top.margin, right.margin)
  margin.names <- c("bottom.margin.total", "left.margin.total", "top.margin.total", "right.margin.total",
                    "bottom.margin", "left.margin", "top.margin", "right.margin")
  for (i in seq_along(margin.list)) {
    if (!is.numeric(margin.list[i]) || length(margin.list[i]) != 1 || is.na(margin.list[i])) {
      stop("Plotting aborted: ", margin.names[i], " must be a single numeric value")
    }
    if (margin.list[i] < 0) {
      stop("Plotting aborted: ", margin.names[i], " must be ≥ 0")
    }
  }
  
  # Validate specified bars.threshold.N
  if (!is.numeric(bars.threshold.N) || length(bars.threshold.N) != 1 || is.na(bars.threshold.N) ||
      bars.threshold.N < 0 || bars.threshold.N %% 1 != 0) {
    stop("Plotting aborted: bars.threshold.N must be a single non-negative integer")
  }
  
  # Validate specified title.font.size
  if (!is.numeric(title.font.size) || length(title.font.size) != 1 || is.na(title.font.size) || title.font.size <= 0) {
    stop("Plotting aborted: title.font.size must be a single positive number")
  }
  
  # Validate specified matrix.label.font.size
  if (!is.numeric(matrix.label.font.size) || length(matrix.label.font.size) != 1 || is.na(matrix.label.font.size) || matrix.label.font.size <= 0) {
    stop("Plotting aborted: matrix.label.font.size must be a single positive number")
  }
  
  # Validate specified bar.label.font.size
  if (!is.numeric(bar.label.font.size) || length(bar.label.font.size) != 1 || is.na(bar.label.font.size) || bar.label.font.size <= 0) {
    stop("Plotting aborted: bar.label.font.size must be a single positive number")
  }
  
  # Validate specified add.boxplot.whiskers
  if (!is.logical(add.boxplot.whiskers) || length(add.boxplot.whiskers) != 1 || is.na(add.boxplot.whiskers)) {
    stop("Plotting aborted: add.boxplot.whiskers must be TRUE or FALSE")
  }
  
  # Validate specified importance.threshold
  if (!is.numeric(importance.threshold) || length(importance.threshold) != 1 || is.na(importance.threshold) || importance.threshold < 0) {
    stop("Plotting aborted: importance.threshold must be a single non-negative number")
  }
  
  # Validate specified set.k
  if (!is.null(set.k) && (!is.numeric(set.k) || length(set.k) != 1 || is.na(set.k) || set.k < 1 || set.k %% 1 != 0)) {
    stop("Plotting aborted: set.k must be NULL or single positive integer")
  }
  
  # Create function to calculate weighted eta squared effect size per variable (weighted by neurons + sample counts)
  calculate.etasquared.per.variable <- function(codebook_matrix, neuron_cluster_vector, som_model, baseline_weight = 1) {
    neuron_cluster_vector <- as.integer(neuron_cluster_vector) #ensure cluster vector is integer
    valid_cluster_rows <- is.finite(neuron_cluster_vector) & !is.na(neuron_cluster_vector) #identify valid cluster rows
    n_units <- length(neuron_cluster_vector) #store original number of neurons (before filtering)
    neuron_sample_counts <- tabulate(som_model$unit.classif, nbins = n_units) #number of samples per neuron (all neurons)
    codebook_matrix <- codebook_matrix[valid_cluster_rows, , drop = FALSE] #subset codebook to valid rows
    neuron_cluster_vector <- neuron_cluster_vector[valid_cluster_rows] #subset cluster vector to valid rows
    neuron_sample_counts <- neuron_sample_counts[valid_cluster_rows] #subset to valid rows
    neuron_weights <- neuron_sample_counts + baseline_weight #neuron weight = baseline + sample support (empty neurons still count)
    etasquared_values <- apply(codebook_matrix, 2, function(variable_values) { #compute eta^2 per variable
      valid_variable_rows <- is.finite(variable_values) & !is.na(variable_values) & is.finite(neuron_weights) & !is.na(neuron_weights) #identify valid values
      variable_values <- variable_values[valid_variable_rows] #subset variable values
      cluster_labels <- neuron_cluster_vector[valid_variable_rows] #subset cluster labels
      weights <- neuron_weights[valid_variable_rows] #subset weights
      if (length(variable_values) < 2) return(NA_real_) #require at least 2 observations
      if (length(unique(cluster_labels)) < 2) return(NA_real_) #require at least 2 clusters
      if (sum(weights) <= 0) return(NA_real_) #require positive total weight
      grand_mean <- sum(weights * variable_values) / sum(weights) #weighted grand mean
      total_sum_of_squares <- sum(weights * (variable_values - grand_mean)^2) #weighted total sum of squares
      if (!is.finite(total_sum_of_squares) || total_sum_of_squares <= 0) return(0) #handle degenerate cases
      cluster_means <- tapply(weights * variable_values, cluster_labels, sum) / tapply(weights, cluster_labels, sum) #weighted cluster means
      cluster_sizes <- tapply(weights, cluster_labels, sum) #cluster "size" = sum of (baseline + hits) across neurons
      between_cluster_sum_of_squares <- sum(cluster_sizes * (cluster_means - grand_mean)^2) #weighted between-cluster sum of squares
      as.numeric(between_cluster_sum_of_squares / total_sum_of_squares) #return eta^2
    })
    etasquared_values #return vector of eta^2 values
  }
  
  # Create function to calculate weighted map variance per variable
  calculate.map.variance.per.variable <- function(codebook_matrix, som_model, baseline_weight = 1) {
    neuron_sample_counts <- tabulate(som_model$unit.classif, nbins = nrow(codebook_matrix)) #number of samples per neuron
    neuron_weights <- neuron_sample_counts + baseline_weight #empty neurons get small positive weight
    variance_values <- apply(codebook_matrix, 2, function(variable_values) {
      valid_rows <- is.finite(variable_values) & !is.na(variable_values) & is.finite(neuron_weights) & !is.na(neuron_weights) #valid rows
      variable_values <- variable_values[valid_rows] #subset values
      variable_weights <- neuron_weights[valid_rows] #subset weights
      if (length(variable_values) < 2) return(NA_real_) #require at least 2 observations
      if (sum(variable_weights) <= 0) return(NA_real_) #require positive total weight
      weighted_mean <- sum(variable_weights * variable_values) / sum(variable_weights) #weighted mean
      weighted_variance <- sum(variable_weights * (variable_values - weighted_mean)^2) / sum(variable_weights) #weighted variance
      weighted_variance
    })
    variance_values
  }
  
  # Extract matrix names from SOM.output
  matrix_names <- SOM.output$input_data_names
  
  # Determine layers from first model
  codes0 <- kohonen::getCodes(SOM.output$som_models[[1]])
  if (!is.list(codes0)) codes0 <- list(codes0)
  num_layers <- length(codes0)
  if (length(matrix_names) != num_layers) matrix_names <- paste0("layer", seq_len(num_layers))
  
  # Filter replicates by set.k (Cluster.separation) or use all replicates (Map.variance)
  if (mode == "Cluster.separation") {
    rep.k <- vapply(SOM.output$som_clusters, function(x) length(unique(x[is.finite(x) & !is.na(x)])), integer(1)) #number of clusters per replicate
    if (all(rep.k < 2L)) {
      message("Eta squared effect size (variable importance) could not be computed because all replicates produced k = 1") #no cluster separation possible
      return(invisible(NULL))
    }
    if (is.null(set.k)) {
      keep.reps <- which(is.finite(rep.k) & !is.na(rep.k) & rep.k >= 2L) #keep only k > 1 replicates
    } else {
      keep.reps <- which(rep.k == set.k) #keep only replicates matching set.k
      if (set.k < 2L) {
        message("Eta squared effect size (variable importance) could not be computed because set.k = 1") #no cluster separation possible
        return(invisible(NULL))
      }
    }
    if (length(keep.reps) == 0) {
      stop("Plotting aborted: no replicates matched set.k")
    }
  }
  if (mode == "Map.variance") {
    keep.reps <- seq_along(SOM.output$som_models)
  }
  
  # Compute per-replicate metric matrices per layer
  metric_layers <- vector("list", num_layers)
  for (i in seq_len(num_layers)) {
    metric_layers[[i]] <- list()
  }
  for (r in keep.reps) {
    som_model <- SOM.output$som_models[[r]]
    codes <- kohonen::getCodes(som_model)
    if (!is.list(codes)) codes <- list(codes)
    if (mode == "Cluster.separation") {
      som_cluster <- SOM.output$som_clusters[[r]]
    }
    for (i in seq_len(num_layers)) {
      if (is.null(colnames(codes[[i]]))) colnames(codes[[i]]) <- paste0("V", seq_len(ncol(codes[[i]]))) #ensure variable names exist
      if (mode == "Cluster.separation") {
        metric_layers[[i]][[paste0("R", r)]] <- calculate.etasquared.per.variable(codes[[i]], som_cluster, som_model)
      }
      if (mode == "Map.variance") {
        metric_layers[[i]][[paste0("R", r)]] <- calculate.map.variance.per.variable(codes[[i]], som_model)
      }
    }
  }
  
  # Build replicate x variable matrices for plotting
  all_layer_metric <- vector("list", num_layers)
  for (i in seq_len(num_layers)) {
    variable_names <- unique(unlist(lapply(metric_layers[[i]], names)))
    metric_matrix <- matrix(NA_real_, nrow = length(metric_layers[[i]]), ncol = length(variable_names))
    rownames(metric_matrix) <- names(metric_layers[[i]])
    colnames(metric_matrix) <- variable_names
    for (j in seq_along(metric_layers[[i]])) {
      metric_values <- metric_layers[[i]][[j]]
      metric_matrix[j, names(metric_values)] <- metric_values
    }
    all_layer_metric[[i]] <- metric_matrix
  }
  
  # Set plot saving
  if (save) {
    
    # Set default plotting name
    if (is.null(file.name)) {
      if (mode == "Cluster.separation") {
        file.name <- paste0("SOM_etasquared_plot_", paste(matrix_names, collapse = "_"), ".", plot.type)
      }
      if (mode == "Map.variance") {
        file.name <- paste0("SOM_map_variance_plot_", paste(matrix_names, collapse = "_"), ".", plot.type)
      }
    }
    
    # Check overwrite
    if (file.exists(file.name) && !overwrite) {
      stop(paste(file.name, "already exists - set overwrite = TRUE to overwrite"))
    }
    
    # Set plot plot.type
    if (plot.type == "svg") {
      grDevices::svg(file.name, 
                     width = width / 2.54, 
                     height = height / 2.54)
    } else if (plot.type == "png") {
      grDevices::png(file.name,
                     width = width, 
                     height = height, 
                     res = resolution, 
                     units = "cm")
    } else if (plot.type == "jpg") {
      grDevices::jpeg(file.name, 
                      width = width, 
                      height = height, 
                      res = resolution, 
                      units = "cm")
    } else {
      stop("Plotting aborted: unsupported file plot.type - choose from 'svg', 'png', or 'jpg'")
    }
  }
  
  # Set plotting area
  graphics::par(mfrow = if (num_layers <= 3) c(1, num_layers) else if (num_layers == 4) c(2, 2) else if (num_layers <= 6) c(2, 3) else if (num_layers <= 8) c(2, 4) else if (num_layers == 9) c(3, 3) else c(ceiling(num_layers / 3), 3),
                oma = c(bottom.margin.total, left.margin.total, top.margin.total, right.margin.total),
                mar = c(bottom.margin, 
                        left.margin, 
                        top.margin, 
                        right.margin))
  
  # Iterate over each matrix and generate plots
  for (i in seq_along(all_layer_metric)) {
    
    # Extract replicate x variable matrix for this layer
    var_mat <- all_layer_metric[[i]]
    if (is.null(var_mat) || nrow(var_mat) == 0 || ncol(var_mat) == 0) {
      message(paste("No values available for", matrix_names[i], "- skipping"))
      next
    }
    
    # Compute median metric per variable across replicates
    median_metric_per_variable <- apply(var_mat, 2, stats::median, na.rm = TRUE)
    median_metric_per_variable[!is.finite(median_metric_per_variable)] <- NA_real_
    
    # Filter variables by threshold
    keep_vars <- which(is.finite(median_metric_per_variable) & !is.na(median_metric_per_variable) & median_metric_per_variable > importance.threshold)
    if (length(keep_vars) == 0) {
      message(paste("No variables exceed importance.threshold of", importance.threshold, "for", matrix_names[i], " - specify lower value for importance.threshold"))
      next
    }
    var_mat <- var_mat[, keep_vars, drop = FALSE]
    
    # Sort variables by median (lowest at bottom)
    median_metric_per_variable <- median_metric_per_variable[colnames(var_mat)]
    sort_idx <- order(median_metric_per_variable, decreasing = FALSE)
    var_mat <- var_mat[, sort_idx, drop = FALSE]
    y_labels <- colnames(var_mat)
    num_bars <- ncol(var_mat)
    
    # Set layer color
    layer_colors <- setNames(col.pal(length(matrix_names)), matrix_names)
    layer.col <- layer_colors[matrix_names[i]]
    
    # Plot boxplots
    boxplot(var_mat,
            horizontal = TRUE,
            las = 1,
            notch = FALSE,
            outline = FALSE,
            col = rep(layer.col, num_bars),
            axes = FALSE,
            whisklty = if(!add.boxplot.whiskers || num_bars > bars.threshold.N) 0 else 1,
            staplelty = if(!add.boxplot.whiskers || num_bars > bars.threshold.N) 0 else 1,
            names = rep("", num_bars))
    axis(1)
    
    # Add y-axis labels as text
    axis(2,
         at = seq_len(num_bars),
         labels = if(num_bars > bars.threshold.N) rep("", num_bars) else y_labels,
         las = 1,
         tick = FALSE,
         cex.axis = bar.label.font.size)
    
    # Add layer title
    mtext(matrix_names[i], 
          side = 3,
          line = 0.3,
          cex = matrix.label.font.size,
          font = 1)
  }
  
  # Add main title in outer margin
  if (mode == "Cluster.separation") {
    mtext("Variable importance of SOM layers (eta squared cluster separation)", outer = TRUE, side = 3, line = 0, cex = title.font.size)
  }
  if (mode == "Map.variance") {
    mtext("Variable importance across SOM map (variation in neuron weights)", outer = TRUE, side = 3, line = 0, cex = title.font.size)
  }
  
  # Close graphics device
  if (save) {
    grDevices::dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}

#' Process SNP data for SOM-based species delimitation
#'
#' Filters and converts SNP data into a numeric SNP dosage matrix suitable for
#' SOM or Super-SOM analyses. The function accepts one of eight input types:
#' a VCF file, an existing `genind` object, an existing `genlight` object, an
#' already processed numeric SNP matrix/data frame, a PLINK `.raw` dosage file,
#' a NEXUS alignment, a FASTA alignment, or a PHYLIP alignment.
#'
#' @description
#' The output is a data frame with samples as rows and SNP variables as columns.
#' For diploid biallelic VCF, `genlight`, PLINK `.raw`, and numeric matrix/data
#' inputs, SNPs are expected to be coded as genotype dosages:
#' `0 = homozygous reference`, `1 = heterozygous`, and
#' `2 = homozygous alternate`. Haploid `0/1` SNP matrices are also accepted.
#' Missing genotypes are coded as `NA`.
#'
#' For biallelic sequence alignments, SNPs are coded as haploid binary
#' presence/absence values: `0 = reference allele`, `1 = alternate allele`,
#' and `NA = missing/ambiguous state`.
#'
#' If `make.biallelic = FALSE`, multiallelic loci from alignment or `genind`
#' inputs are retained by encoding each polymorphic locus as `k - 1`
#' allele-dosage columns, where `k` is the number of observed non-missing alleles
#' at that locus. For direct dosage-matrix-like inputs (`snp.matrix.input`,
#' `genlight.input`, and `plink.raw.path`), `make.biallelic` does not create
#' multiallelic encodings because those inputs are already dosage matrices.
#'
#' @details
#' Filtering is applied in the following general order:
#'
#' \enumerate{
#'   \item Remove non-biallelic loci if `make.biallelic = TRUE` and the input contains explicit allele/locus structure.
#'   \item Remove loci exceeding `missing.loci.cutoff.lenient`.
#'   \item Remove individuals exceeding `missing.individuals.cutoff`.
#'   \item Remove true singleton loci if `singleton.loci.filter = TRUE`.
#'   \item Remove loci exceeding `missing.loci.cutoff.final`.
#'   \item Remove invariant loci if `invariant.loci.filter = TRUE`.
#' }
#'
#' Removal messages are printed only when at least one locus or individual is
#' actually removed. A blank separator line is printed before the final summary
#' only if at least one filtering message has been printed.
#'
#' The fast VCF path supports standard diploid biallelic genotype strings:
#' `0/0`, `0/1`, `1/0`, `1/1`, and phased equivalents using `|`. Missing
#' genotypes can be represented as `./.`, `.`, partial-missing diploid
#' genotypes such as `0/.`, or `NA`. Other genotype codings trigger fallback
#' to the `genind` workflow.
#'
#' For NEXUS/FASTA/PHYLIP alignment input, the following symbols are treated as
#' missing or ambiguous and are excluded from allele counting:
#' `?`, `-`, `.`, `N`, `R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`, and `X`.
#' In NEXUS files, `.` is first expanded as a match character relative to the
#' first sequence; unresolved `.` symbols are then treated as missing.
#'
#' PLINK `.raw` input is expected to contain the standard metadata columns
#' `FID`, `IID`, `PAT`, `MAT`, `SEX`, and `PHENOTYPE`, followed by dosage
#' columns coded as `0`, `1`, `2`, or missing.
#'
#' @param vcf.path Character or `NULL`. Path to a VCF file.
#' @param genind.input A `genind` object or `NULL`.
#' @param genlight.input A `genlight` object or `NULL`.
#' @param snp.matrix.input Numeric matrix/data frame or `NULL`. Rows must be
#'   samples and columns must be SNP variables. Values must be `0`, `1`, `2`,
#'   or `NA`; haploid `0/1` matrices are also accepted.
#' @param plink.raw.path Character or `NULL`. Path to a PLINK `.raw` dosage file.
#' @param nexus.path Character or `NULL`. Path to a NEXUS sequence alignment.
#' @param fasta.path Character or `NULL`. Path to a FASTA sequence alignment.
#' @param phylip.path Character or `NULL`. Path to a PHYLIP sequence alignment.
#' @param phylip.format Character. PHYLIP format passed to `ape::read.dna()`.
#'   Must be `"sequential"` or `"interleaved"`.
#' @param plink.raw.metadata.columns Character vector of PLINK `.raw` metadata
#'   columns to exclude before extracting SNP dosage columns.
#' @param snp.matrix.ploidy Integer. Ploidy used for dosage-matrix-like inputs:
#'   `snp.matrix.input`, `genlight.input`, and `plink.raw.path`. Use `2` for
#'   diploid dosage matrices coded as `0/1/2`, and `1` for haploid/binary
#'   matrices coded as `0/1`. This argument is not used for VCF, `genind`,
#'   NEXUS, FASTA, or PHYLIP inputs.
#' @param make.biallelic Logical. If `TRUE`, retain only biallelic loci where
#'   explicit allele/locus structure is available. If `FALSE`, retain
#'   multiallelic loci as `k - 1` allele-dosage columns where supported.
#' @param missing.loci.cutoff.lenient Numeric between 0 and 1 or `NULL`.
#'   First-pass missing-data cutoff for loci. Loci with a missing proportion
#'   greater than this value are removed.
#' @param missing.loci.cutoff.final Numeric between 0 and 1 or `NULL`.
#'   Final stricter missing-data cutoff for loci.
#' @param missing.individuals.cutoff Numeric between 0 and 1 or `NULL`.
#'   Missing-data cutoff for individuals. Individuals with a missing proportion
#'   greater than this value are removed.
#' @param singleton.loci.filter Logical. If `TRUE`, remove true singleton loci.
#'   Invariant loci are not counted as singletons; they are handled by the
#'   invariant-locus filter.
#' @param invariant.loci.filter Logical. If `TRUE`, remove invariant loci.
#' @param verbose Logical. If `TRUE`, print filtering messages and final matrix
#'   dimensions.
#'
#' @return
#' A data frame containing the processed SNP matrix. Rows are samples and
#' columns are SNP variables. Values are numeric/integer genotype or allele
#' dosages with missing values coded as `NA`.
#'
#' @export
process.SNP.data.SOM <- function(vcf.path = NULL, #optional path to VCF file
                                 genind.input = NULL, #optional genind object
                                 genlight.input = NULL, #optional genlight object
                                 snp.matrix.input = NULL, #optional SNP dosage matrix/data frame
                                 plink.raw.path = NULL, #optional path to PLINK .raw dosage file
                                 nexus.path = NULL, #optional path to NEXUS alignment file
                                 fasta.path = NULL, #optional path to FASTA alignment file
                                 phylip.path = NULL, #optional path to PHYLIP alignment file
                                 phylip.format = "sequential", #PHYLIP format: sequential or interleaved
                                 plink.raw.metadata.columns = c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"), #PLINK .raw metadata columns
                                 snp.matrix.ploidy = 2, #ploidy for snp.matrix.input, genlight.input, and plink.raw.path
                                 make.biallelic = TRUE, #whether to restrict to biallelic loci or retain multiallelic loci as k-1 dosage columns
                                 missing.loci.cutoff.lenient = 0.7, #remove loci with > this proportion missing (1st lenient filter)
                                 missing.loci.cutoff.final = 0.5, #remove loci with > this proportion missing (2nd final stringent filter)
                                 missing.individuals.cutoff = 0.5, #remove individuals with > this proportion missing
                                 singleton.loci.filter = TRUE, #whether to remove singleton loci
                                 invariant.loci.filter = TRUE, #whether to remove invariant loci
                                 verbose = TRUE #whether to show filtering messages and final summary
) {
  
  # Track whether any filtering message was printed
  filter.messages.printed <- FALSE
  
  # Create function to print filter messages
  print.filter.message <- function(...) {
    if (isTRUE(verbose)) {
      message(...)
      filter.messages.printed <<- TRUE
    }
  }
  
  # Create function to print final summary
  print.final.summary <- function(snp.matrix) {
    if (verbose && filter.messages.printed) message("")
    if (verbose) message("Final SNP matrix: ", nrow(snp.matrix), " samples × ", ncol(snp.matrix), " loci") #summary
  }
  
  # Create function to validate arguments
  validate.arguments <- function() {
    provided.input.count <- sum(!is.null(vcf.path),
                                !is.null(genind.input),
                                !is.null(genlight.input),
                                !is.null(snp.matrix.input),
                                !is.null(plink.raw.path),
                                !is.null(nexus.path),
                                !is.null(fasta.path),
                                !is.null(phylip.path)) #count provided inputs
    if (provided.input.count != 1) stop("Provide exactly one of vcf.path, genind.input, genlight.input, snp.matrix.input, plink.raw.path, nexus.path, fasta.path, or phylip.path") #must provide only one
    if (!is.logical(make.biallelic) || length(make.biallelic) != 1 || is.na(make.biallelic)) stop("make.biallelic must be TRUE or FALSE") #validate
    if (!is.logical(singleton.loci.filter) || length(singleton.loci.filter) != 1 || is.na(singleton.loci.filter)) stop("singleton.loci.filter must be TRUE or FALSE") #validate
    if (!is.logical(invariant.loci.filter) || length(invariant.loci.filter) != 1 || is.na(invariant.loci.filter)) stop("invariant.loci.filter must be TRUE or FALSE") #validate
    if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) stop("verbose must be TRUE or FALSE") #validate
    if (!is.null(missing.loci.cutoff.lenient) && (!is.numeric(missing.loci.cutoff.lenient) || length(missing.loci.cutoff.lenient) != 1 || is.na(missing.loci.cutoff.lenient) || missing.loci.cutoff.lenient < 0 || missing.loci.cutoff.lenient > 1)) stop("missing.loci.cutoff.lenient must be NULL or a single numeric value between 0 and 1") #validate
    if (!is.null(missing.loci.cutoff.final) && (!is.numeric(missing.loci.cutoff.final) || length(missing.loci.cutoff.final) != 1 || is.na(missing.loci.cutoff.final) || missing.loci.cutoff.final < 0 || missing.loci.cutoff.final > 1)) stop("missing.loci.cutoff.final must be NULL or a single numeric value between 0 and 1") #validate
    if (!is.null(missing.individuals.cutoff) && (!is.numeric(missing.individuals.cutoff) || length(missing.individuals.cutoff) != 1 || is.na(missing.individuals.cutoff) || missing.individuals.cutoff < 0 || missing.individuals.cutoff > 1)) stop("missing.individuals.cutoff must be NULL or a single numeric value between 0 and 1") #validate
    if (!is.character(phylip.format) || length(phylip.format) != 1 || !(phylip.format %in% c("sequential", "interleaved"))) stop("phylip.format must be 'sequential' or 'interleaved'") #validate
    if (!is.character(plink.raw.metadata.columns)) stop("plink.raw.metadata.columns must be a character vector") #validate
    if (!is.null(vcf.path) && !file.exists(vcf.path)) stop("VCF file does not exist: ", vcf.path) #check vcf path
    if (!is.null(plink.raw.path) && !file.exists(plink.raw.path)) stop("PLINK .raw file does not exist: ", plink.raw.path) #check PLINK .raw path
    if (!is.numeric(snp.matrix.ploidy) || length(snp.matrix.ploidy) != 1 || is.na(snp.matrix.ploidy) || !(snp.matrix.ploidy %in% c(1, 2))) stop("snp.matrix.ploidy must be 1 or 2") #validate ploidy
    if (!is.null(nexus.path) && !file.exists(nexus.path)) stop("NEXUS file does not exist: ", nexus.path) #check nexus path
    if (!is.null(fasta.path) && !file.exists(fasta.path)) stop("FASTA file does not exist: ", fasta.path) #check fasta path
    if (!is.null(phylip.path) && !file.exists(phylip.path)) stop("PHYLIP file does not exist: ", phylip.path) #check phylip path
    if (!is.null(genlight.input) && !inherits(genlight.input, "genlight")) stop("genlight.input must be a genlight object") #check genlight input
    if (!is.null(genind.input) && !inherits(genind.input, "genind")) stop("genind.input must be a genind object") #check genind input
    if (!is.null(snp.matrix.input) && !(is.matrix(snp.matrix.input) || is.data.frame(snp.matrix.input))) stop("snp.matrix.input must be a matrix or data frame") #check SNP matrix input
    if (!is.null(vcf.path) && !requireNamespace("vcfR", quietly = TRUE)) stop("Package 'vcfR' is required for vcf.path") #check vcfR package
    if ((!is.null(genind.input) || !is.null(genlight.input) || (!is.null(vcf.path) && !isTRUE(make.biallelic))) && !requireNamespace("adegenet", quietly = TRUE)) stop("Package 'adegenet' is required for genind.input, genlight.input, or non-biallelic/fallback VCF processing") #check adegenet package
    if ((!is.null(nexus.path) || !is.null(fasta.path) || !is.null(phylip.path)) && !requireNamespace("ape", quietly = TRUE)) stop("Package 'ape' is required for NEXUS, FASTA, or PHYLIP input") #check ape package
  }
  
  # Create function to coerce and validate dosage matrix
  coerce.and.validate.dosage.matrix <- function(snp.matrix, input.name = "SNP matrix", snp.matrix.ploidy = 2) {
    if (is.data.frame(snp.matrix)) snp.matrix <- as.matrix(snp.matrix) #convert data frame
    if (!is.matrix(snp.matrix)) stop(input.name, " must be a matrix or data frame") #check matrix
    if (nrow(snp.matrix) == 0 || ncol(snp.matrix) == 0) stop(input.name, " is empty") #check empty
    
    original.dimnames <- dimnames(snp.matrix) #store dimnames
    
    if (!is.numeric(snp.matrix) && !is.integer(snp.matrix)) {
      character.matrix <- matrix(as.character(snp.matrix), nrow = nrow(snp.matrix), ncol = ncol(snp.matrix), dimnames = original.dimnames) #convert to character
      character.matrix[character.matrix %in% c("", ".", "NA", "NaN", "nan")] <- NA_character_ #standardize missing strings
      numeric.values <- suppressWarnings(as.numeric(character.matrix)) #convert to numeric
      created.missing <- is.na(numeric.values) & !is.na(as.vector(character.matrix)) #failed conversions
      if (any(created.missing)) stop(input.name, " contains non-numeric non-missing values") #stop on invalid values
      snp.matrix <- matrix(numeric.values, nrow = nrow(character.matrix), ncol = ncol(character.matrix), dimnames = original.dimnames) #restore matrix
    } else {
      storage.mode(snp.matrix) <- "numeric" #coerce numeric
    }
    
    valid.values <- if (snp.matrix.ploidy == 1) c(0, 1) else c(0, 1, 2) #valid dosage values
    invalid.values <- !is.na(snp.matrix) & !(snp.matrix %in% valid.values) #invalid dosage values
    if (any(invalid.values)) stop(input.name, " must contain only ", paste(valid.values, collapse = ", "), ", or NA values for snp.matrix.ploidy = ", snp.matrix.ploidy) #validate dosage values
    
    if (is.null(rownames(snp.matrix))) rownames(snp.matrix) <- paste0("Sample", seq_len(nrow(snp.matrix))) #fallback sample names
    if (is.null(colnames(snp.matrix))) colnames(snp.matrix) <- paste0("SNP", seq_len(ncol(snp.matrix))) #fallback SNP names
    rownames(snp.matrix) <- make.unique(as.character(rownames(snp.matrix))) #unique rownames
    colnames(snp.matrix) <- make.unique(as.character(colnames(snp.matrix))) #unique colnames
    
    storage.mode(snp.matrix) <- "integer" #store as integer dosage
    return(snp.matrix) #return dosage matrix
  }
  
  # Create function to filter SNP dosage matrix
  filter.SNP.matrix.fast <- function(snp.matrix,
                                     missing.loci.cutoff.lenient = 0.5,
                                     missing.loci.cutoff.final = 0.2,
                                     missing.individuals.cutoff = 0.4,
                                     singleton.loci.filter = TRUE,
                                     invariant.loci.filter = TRUE,
                                     snp.matrix.ploidy = 2
                                     
  ) {
    
    snp.matrix <- coerce.and.validate.dosage.matrix(snp.matrix = snp.matrix, input.name = "SNP dosage matrix", snp.matrix.ploidy = snp.matrix.ploidy) #validate matrix
    dosage.ploidy <- as.integer(snp.matrix.ploidy) #use declared ploidy
    invariant.loci.removed.total <- 0L #track total invariant loci removed
    invariant.loci.reference.count <- NA_integer_ #track denominator for invariant message
    
    # Loci missing (lenient)
    if (!is.null(missing.loci.cutoff.lenient)) {
      locus.count.before.lenient.missing.filter <- ncol(snp.matrix) #loci before filter
      keep.loci <- colMeans(is.na(snp.matrix)) <= missing.loci.cutoff.lenient #loci passing filter
      snp.matrix <- snp.matrix[, keep.loci, drop = FALSE] #filter loci
      loci.removed <- locus.count.before.lenient.missing.filter - ncol(snp.matrix) #number removed
      if (loci.removed > 0) print.filter.message(loci.removed, " of ", locus.count.before.lenient.missing.filter, " loci removed due to >", missing.loci.cutoff.lenient * 100, "% missing data (lenient filter)") #report
      if (ncol(snp.matrix) == 0) stop("All loci removed after lenient missing data filter") #stop if all gone
    }
    
    # Individuals missing
    if (!is.null(missing.individuals.cutoff)) {
      individual.count.before.missing.filter <- nrow(snp.matrix) #individuals before filter
      keep.individuals <- rowMeans(is.na(snp.matrix)) <= missing.individuals.cutoff #individuals passing filter
      snp.matrix <- snp.matrix[keep.individuals, , drop = FALSE] #filter individuals
      individuals.removed <- individual.count.before.missing.filter - nrow(snp.matrix) #number removed
      if (individuals.removed > 0) print.filter.message(individuals.removed, " of ", individual.count.before.missing.filter, " individuals removed due to >", missing.individuals.cutoff * 100, "% missing data") #report
      if (nrow(snp.matrix) == 0) stop("All individuals removed after missing data filter") #stop if all gone
    }
    
    # Singleton loci filter
    if (isTRUE(singleton.loci.filter)) {
      locus.count.before.singleton.filter <- ncol(snp.matrix) #loci before filter
      alternate.allele.counts <- colSums(snp.matrix, na.rm = TRUE) #alternate allele counts
      called.chromosome.counts <- dosage.ploidy * colSums(!is.na(snp.matrix)) #called chromosomes
      minor.allele.counts <- pmin(alternate.allele.counts, called.chromosome.counts - alternate.allele.counts) #minor allele counts
      keep.loci <- minor.allele.counts != 1L #remove only true singleton loci; keep invariant loci for final invariant filter
      snp.matrix <- snp.matrix[, keep.loci, drop = FALSE] #filter loci
      loci.removed <- locus.count.before.singleton.filter - ncol(snp.matrix) #number removed
      if (loci.removed > 0) print.filter.message(loci.removed, " of ", locus.count.before.singleton.filter, " singleton loci removed") #report
      if (ncol(snp.matrix) == 0) stop("All loci removed after singleton filter") #stop if all gone
    }
    
    # Loci missing (strict)
    if (!is.null(missing.loci.cutoff.final)) {
      locus.count.before.strict.missing.filter <- ncol(snp.matrix) #loci before filter
      keep.loci <- colMeans(is.na(snp.matrix)) <= missing.loci.cutoff.final #loci passing filter
      snp.matrix <- snp.matrix[, keep.loci, drop = FALSE] #filter loci
      loci.removed <- locus.count.before.strict.missing.filter - ncol(snp.matrix) #number removed
      if (loci.removed > 0) print.filter.message(loci.removed, " of ", locus.count.before.strict.missing.filter, " loci removed due to >", missing.loci.cutoff.final * 100, "% missing data (stricter filter)") #report
      if (ncol(snp.matrix) == 0) stop("All loci removed after final missing data filter") #stop if all gone
    }
    
    # Invariant loci filter
    if (isTRUE(invariant.loci.filter)) {
      locus.count.before.final.invariant.filter <- ncol(snp.matrix) #loci before filter
      invariant.loci.reference.count <- locus.count.before.final.invariant.filter #store denominator
      locus.minimum.values <- suppressWarnings(apply(snp.matrix, 2, min, na.rm = TRUE)) #minimum dosage per locus
      locus.maximum.values <- suppressWarnings(apply(snp.matrix, 2, max, na.rm = TRUE)) #maximum dosage per locus
      keep.loci <- is.finite(locus.minimum.values) & is.finite(locus.maximum.values) & locus.minimum.values != locus.maximum.values #variable loci
      snp.matrix <- snp.matrix[, keep.loci, drop = FALSE] #filter loci
      invariant.loci.removed.total <- invariant.loci.removed.total + locus.count.before.final.invariant.filter - ncol(snp.matrix) #update total removed
      if (invariant.loci.removed.total > 0) print.filter.message(invariant.loci.removed.total, " of ", invariant.loci.reference.count, " invariant loci removed") #report
      if (ncol(snp.matrix) == 0) stop("All loci removed after invariant filter") #stop if all gone
    }
    
    return(as.data.frame(snp.matrix, stringsAsFactors = FALSE)) #return data frame
  }
  
  # Create function to process direct dosage matrix inputs
  process.dosage.matrix.input <- function(snp.matrix, input.name) {
    snp.matrix <- coerce.and.validate.dosage.matrix(snp.matrix = snp.matrix, input.name = input.name, snp.matrix.ploidy = snp.matrix.ploidy) #validate matrix
    snp.matrix <- filter.SNP.matrix.fast(snp.matrix = snp.matrix,
                                         missing.loci.cutoff.lenient = missing.loci.cutoff.lenient,
                                         missing.loci.cutoff.final = missing.loci.cutoff.final,
                                         missing.individuals.cutoff = missing.individuals.cutoff,
                                         singleton.loci.filter = singleton.loci.filter,
                                         invariant.loci.filter = invariant.loci.filter,
                                         snp.matrix.ploidy = snp.matrix.ploidy) #filter matrix
    print.final.summary(snp.matrix) #summary
    return(snp.matrix) #return SNP matrix
  }
  
  # Create function to read PLINK .raw dosage file
  read.plink.raw.input <- function(plink.raw.path) {
    plink.raw.data <- utils::read.table(plink.raw.path,
                                        header = TRUE,
                                        stringsAsFactors = FALSE,
                                        check.names = FALSE,
                                        comment.char = "",
                                        na.strings = c("NA", "NaN", ".", "-9", "")) #read PLINK .raw file
    
    if (nrow(plink.raw.data) == 0 || ncol(plink.raw.data) == 0) stop("PLINK .raw file is empty") #check empty
    genotype.columns <- setdiff(colnames(plink.raw.data), plink.raw.metadata.columns) #SNP columns
    if (length(genotype.columns) == 0) stop("No SNP dosage columns found in PLINK .raw file") #check SNP columns
    
    snp.matrix <- as.matrix(plink.raw.data[, genotype.columns, drop = FALSE]) #extract SNP matrix
    
    if ("IID" %in% colnames(plink.raw.data)) {
      sample.names <- as.character(plink.raw.data[["IID"]]) #sample names from IID
      if (any(is.na(sample.names)) || any(sample.names == "") || any(duplicated(sample.names))) {
        if (all(c("FID", "IID") %in% colnames(plink.raw.data))) {
          sample.names <- paste(as.character(plink.raw.data[["FID"]]), as.character(plink.raw.data[["IID"]]), sep = "_") #fallback FID_IID
        } else {
          sample.names <- paste0("Sample", seq_len(nrow(plink.raw.data))) #fallback names
        }
      }
    } else {
      sample.names <- paste0("Sample", seq_len(nrow(plink.raw.data))) #fallback names
    }
    
    rownames(snp.matrix) <- make.unique(sample.names) #set rownames
    colnames(snp.matrix) <- make.unique(genotype.columns) #set colnames
    return(snp.matrix) #return matrix
  }
  
  # Create function to convert genlight object to dosage matrix
  convert.genlight.to.dosage.matrix <- function(genlight.input) {
    if (!requireNamespace("adegenet", quietly = TRUE)) stop("Package 'adegenet' is required for genlight.input") #check package
    snp.matrix <- suppressWarnings(as.matrix(genlight.input)) #convert to matrix
    if (is.null(dim(snp.matrix))) stop("genlight.input could not be converted to a matrix") #check matrix
    
    individual.names <- tryCatch(adegenet::indNames(genlight.input), error = function(error) NULL) #individual names
    locus.names <- tryCatch(adegenet::locNames(genlight.input), error = function(error) NULL) #locus names
    
    if (!is.null(individual.names) && length(individual.names) > 0 && ncol(snp.matrix) == length(individual.names) && nrow(snp.matrix) != length(individual.names)) {
      snp.matrix <- t(snp.matrix) #transpose if orientation is loci x individuals
    }
    
    if (!is.null(individual.names) && length(individual.names) == nrow(snp.matrix)) rownames(snp.matrix) <- individual.names #set rownames
    if (!is.null(locus.names) && length(locus.names) == ncol(snp.matrix)) colnames(snp.matrix) <- locus.names #set colnames
    
    return(snp.matrix) #return matrix
  }
  
  # Create function to encode multiallelic alignment as k-1 dosage columns
  encode.multiallelic.alignment.as.k.minus.one <- function(alignment.matrix,
                                                           missing.symbols = c("?", "-", "N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "X")
  ) {
    observed.alleles.per.locus <- lapply(seq_len(ncol(alignment.matrix)), function(locus.index) sort(setdiff(unique(alignment.matrix[, locus.index]), missing.symbols))) #observed alleles
    retained.alleles.per.locus <- lapply(observed.alleles.per.locus, function(observed.alleles) if (length(observed.alleles) < 2) character(0) else observed.alleles[-length(observed.alleles)]) #retained k-1 alleles
    retained.column.counts <- lengths(retained.alleles.per.locus) #number of retained columns per locus
    total.retained.columns <- sum(retained.column.counts) #total retained columns
    
    if (total.retained.columns == 0) return(data.frame(row.names = rownames(alignment.matrix))) #return empty if none
    
    multiallelic.snp.matrix <- matrix(NA_integer_,
                                      nrow = nrow(alignment.matrix),
                                      ncol = total.retained.columns,
                                      dimnames = list(rownames(alignment.matrix), NULL)) #initialize matrix
    multiallelic.snp.names <- character(total.retained.columns) #initialize names
    output.column.index <- 1L #initialize output column index
    
    for(locus.index in seq_len(ncol(alignment.matrix))) {
      retained.alleles <- retained.alleles.per.locus[[locus.index]] #retained alleles
      if (length(retained.alleles) == 0) next #skip invariant loci
      missing.values.at.locus <- alignment.matrix[, locus.index] %in% missing.symbols #missing values at locus
      
      for(retained.allele in retained.alleles) {
        multiallelic.snp.matrix[, output.column.index] <- ifelse(missing.values.at.locus, NA_integer_, ifelse(alignment.matrix[, locus.index] == retained.allele, 1L, 0L)) #haploid allele dosage
        multiallelic.snp.names[output.column.index] <- paste0("SNP", locus.index, ".", retained.allele) #set column name
        output.column.index <- output.column.index + 1L #advance column
      }
    }
    
    colnames(multiallelic.snp.matrix) <- multiallelic.snp.names #set colnames
    return(as.data.frame(multiallelic.snp.matrix, stringsAsFactors = FALSE)) #return matrix
  }
  
  # Create function to count observed alleles per locus in genind object
  count.observed.alleles.per.locus <- function(genind.object) {
    allele.count.matrix <- suppressMessages(suppressWarnings(adegenet::tab(genind.object, NA.method = "asis"))) #allele count matrix
    locus.names <- adegenet::locNames(genind.object) #locus names
    locus.membership.vector <- factor(as.character(genind.object@loc.fac), levels = locus.names) #locus factor
    locus.column.list <- split(seq_along(locus.membership.vector), locus.membership.vector, drop = TRUE) #columns per locus
    allele.present <- colSums(allele.count.matrix, na.rm = TRUE) > 0 #observed allele columns
    observed.allele.counts <- vapply(locus.column.list, function(locus.column.indices) sum(allele.present[locus.column.indices]), integer(1)) #observed alleles per locus
    return(observed.allele.counts) #return counts
  }
  
  # Create function to count minimum observed allele count per locus in genind object
  count.minor.alleles.per.locus <- function(genind.object) {
    allele.count.matrix <- suppressMessages(suppressWarnings(adegenet::tab(genind.object, NA.method = "asis"))) #allele count matrix
    locus.names <- adegenet::locNames(genind.object) #locus names
    locus.membership.vector <- factor(as.character(genind.object@loc.fac), levels = locus.names) #locus factor
    locus.column.list <- split(seq_along(locus.membership.vector), locus.membership.vector, drop = TRUE) #columns per locus
    
    minor.allele.counts <- vapply(locus.column.list, function(locus.column.indices) {
      allele.counts <- colSums(allele.count.matrix[, locus.column.indices, drop = FALSE], na.rm = TRUE) #allele counts
      observed.allele.counts <- allele.counts[allele.counts > 0] #observed allele counts
      if (length(observed.allele.counts) < 2) return(0L) #invariant locus
      return(as.integer(min(observed.allele.counts))) #minimum observed allele count
    }, integer(1)) #minor allele counts
    
    return(minor.allele.counts) #return counts
  }
  
  # Create function to convert genind object to k-1 allele dosage matrix
  convert.genind.to.k.minus.one.dosage <- function(genind.object) {
    allele.count.matrix <- suppressMessages(suppressWarnings(adegenet::tab(genind.object, NA.method = "asis"))) #allele count matrix
    locus.names <- adegenet::locNames(genind.object) #locus names
    locus.membership.vector <- factor(as.character(genind.object@loc.fac), levels = locus.names) #locus factor
    locus.column.list <- split(seq_along(locus.membership.vector), locus.membership.vector, drop = TRUE) #columns per locus
    allele.present <- colSums(allele.count.matrix, na.rm = TRUE) > 0 #observed allele columns
    
    retained.column.indices <- unlist(lapply(locus.column.list, function(locus.column.indices) {
      observed.column.indices <- locus.column.indices[allele.present[locus.column.indices]] #observed columns
      if (length(observed.column.indices) < 2) return(integer(0)) #skip invariant loci
      observed.column.indices[-length(observed.column.indices)] #retain k-1 columns
    }), use.names = FALSE) #retained columns
    
    if (length(retained.column.indices) == 0) return(data.frame(row.names = adegenet::indNames(genind.object))) #return empty if none
    multiallelic.snp.matrix <- as.data.frame(allele.count.matrix[, retained.column.indices, drop = FALSE], stringsAsFactors = FALSE) #subset once
    rownames(multiallelic.snp.matrix) <- adegenet::indNames(genind.object) #set rownames
    colnames(multiallelic.snp.matrix) <- colnames(allele.count.matrix)[retained.column.indices] #set colnames
    return(multiallelic.snp.matrix) #return matrix
  }
  
  # Create function to count singleton loci in character alignment matrix
  find.singleton.alignment.loci <- function(alignment.matrix, missing.symbols) {
    singleton.present <- vapply(seq_len(ncol(alignment.matrix)), function(locus.index) {
      allele.values <- alignment.matrix[, locus.index] #values at locus
      allele.values <- allele.values[!(allele.values %in% missing.symbols)] #remove missing and ambiguous values
      if (length(allele.values) == 0) return(FALSE) #no called alleles
      allele.counts <- table(allele.values) #allele counts
      if (length(allele.counts) < 2) return(FALSE) #do not count invariant loci as singleton loci
      any(as.integer(allele.counts) == 1L) #TRUE if any observed allele is singleton
    }, logical(1)) #singleton status per locus
    
    return(singleton.present) #return singleton status
  }
  
  # Create function to convert alignment matrix to sequence list
  convert.alignment.matrix.to.sequence.list <- function(alignment.matrix) {
    if (is.null(dim(alignment.matrix))) stop("Alignment input could not be converted to a matrix") #check matrix
    sequence.list <- lapply(seq_len(nrow(alignment.matrix)), function(row.index) alignment.matrix[row.index, ]) #convert rows to list
    if (!is.null(rownames(alignment.matrix))) {
      names(sequence.list) <- rownames(alignment.matrix) #use rownames
    } else {
      names(sequence.list) <- paste0("Sample", seq_len(nrow(alignment.matrix))) #fallback names
    }
    return(sequence.list) #return sequence list
  }
  
  # Create function to expand NEXUS match characters
  expand.nexus.match.characters <- function(alignment.matrix, match.symbol = ".") {
    if (nrow(alignment.matrix) == 0 || ncol(alignment.matrix) == 0) return(alignment.matrix) #return empty matrix unchanged
    match.position.matrix <- alignment.matrix == match.symbol #positions with match character
    if (!any(match.position.matrix, na.rm = TRUE)) return(alignment.matrix) #return if no match characters
    reference.sequence <- alignment.matrix[1, ] #first sequence is NEXUS match reference
    for(locus.index in seq_len(ncol(alignment.matrix))) {
      replacement.allele <- reference.sequence[locus.index] #reference allele at locus
      if (!is.na(replacement.allele) && replacement.allele != match.symbol) {
        alignment.matrix[match.position.matrix[, locus.index], locus.index] <- replacement.allele #replace match symbols
      }
    }
    alignment.matrix[alignment.matrix == match.symbol] <- "?" #unresolved match symbols become missing
    return(alignment.matrix) #return expanded alignment
  }
  
  # Create function to generate aligned sequence processor for NEXUS/FASTA/PHYLIP
  process.alignment.input <- function(sequence.list, file.type) {
    if (length(sequence.list) == 0) stop("No sequences found in ", file.type, " file") #check for empty
    if (is.null(names(sequence.list))) names(sequence.list) <- paste0("Sample", seq_along(sequence.list)) #fallback names
    sequence.lengths <- vapply(sequence.list, function(single.sequence) nchar(paste0(single.sequence, collapse = "")), integer(1)) #sequence lengths
    if (any(sequence.lengths == 0)) stop("Empty sequences detected in ", file.type, " file") #empty sequence
    if (length(unique(sequence.lengths)) != 1) stop(file.type, " file is not aligned: sequences have different lengths") #not aligned
    
    alignment.matrix <- t(sapply(sequence.list, function(single.sequence) strsplit(paste0(single.sequence, collapse = ""), "")[[1]])) #make samples × sites character matrix
    if (any(dim(alignment.matrix) == 0)) stop(file.type, " alignment matrix is empty or malformed") #check for empty matrix
    alignment.matrix <- toupper(alignment.matrix) #standardize case
    rownames(alignment.matrix) <- make.unique(as.character(names(sequence.list))) #set rownames
    if (identical(file.type, "NEXUS")) alignment.matrix <- expand.nexus.match.characters(alignment.matrix = alignment.matrix, match.symbol = ".") #expand NEXUS match characters
    
    missing.symbols <- c("?", "-", ".", "N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "X") #define missing and ambiguous symbols
    observed.alleles.per.site <- lapply(seq_len(ncol(alignment.matrix)), function(locus.index) sort(setdiff(unique(alignment.matrix[, locus.index]), missing.symbols))) #observed alleles per site
    total.locus.count.before.filtering <- ncol(alignment.matrix) #total loci before filtering
    
    # Biallelic alignment processing
    if (isTRUE(make.biallelic)) {
      biallelic.site.indices <- which(sapply(observed.alleles.per.site, length) == 2) #index of biallelic sites
      loci.removed <- total.locus.count.before.filtering - length(biallelic.site.indices) #number removed
      if (loci.removed > 0) print.filter.message(loci.removed, " of ", total.locus.count.before.filtering, " loci removed because they were not biallelic") #report non-biallelic filter
      if (length(biallelic.site.indices) == 0) stop("No biallelic SNPs found in alignment") #stop if none
      
      biallelic.alignment.matrix <- alignment.matrix[, biallelic.site.indices, drop = FALSE] #keep only biallelic sites
      biallelic.snp.matrix <- matrix(NA_integer_,
                                     nrow = nrow(biallelic.alignment.matrix),
                                     ncol = ncol(biallelic.alignment.matrix),
                                     dimnames = list(rownames(biallelic.alignment.matrix), paste0("SNP", seq_len(ncol(biallelic.alignment.matrix))))) #initialize matrix
      
      # Recoding to 0/1
      alleles.at.biallelic.sites <- observed.alleles.per.site[biallelic.site.indices] #alleles at biallelic sites
      ref.alleles <- vapply(alleles.at.biallelic.sites, `[`, character(1L), 1L) #reference allele per site
      alt.alleles <- vapply(alleles.at.biallelic.sites, `[`, character(1L), 2L) #alternate allele per site
      ref.broadcast.matrix <- matrix(ref.alleles, nrow = nrow(biallelic.alignment.matrix), ncol = length(ref.alleles), byrow = TRUE) #broadcast ref alleles
      alt.broadcast.matrix <- matrix(alt.alleles, nrow = nrow(biallelic.alignment.matrix), ncol = length(alt.alleles), byrow = TRUE) #broadcast alt alleles
      biallelic.snp.matrix[biallelic.alignment.matrix == ref.broadcast.matrix] <- 0L #assign reference
      biallelic.snp.matrix[biallelic.alignment.matrix == alt.broadcast.matrix] <- 1L #assign alternate
      
      invariant.loci.removed.total <- 0L #track total invariant loci removed
      invariant.loci.reference.count <- NA_integer_ #track denominator for invariant message
      
      # Loci missing (lenient)
      if (!is.null(missing.loci.cutoff.lenient)) {
        locus.count.before.lenient.missing.filter <- ncol(biallelic.snp.matrix) #loci before filter
        biallelic.snp.matrix <- biallelic.snp.matrix[, (colMeans(is.na(biallelic.snp.matrix)) <= missing.loci.cutoff.lenient), drop = FALSE] #filter loci
        loci.removed <- locus.count.before.lenient.missing.filter - ncol(biallelic.snp.matrix) #number removed
        if (loci.removed > 0) print.filter.message(loci.removed, " of ", locus.count.before.lenient.missing.filter, " loci removed due to >", missing.loci.cutoff.lenient * 100, "% missing data (lenient filter)") #report
        if (ncol(biallelic.snp.matrix) == 0) stop("All loci removed after lenient missing data filter") #stop if all gone
      }
      
      # Individuals missing
      if (!is.null(missing.individuals.cutoff)) {
        individual.count.before.missing.filter <- nrow(biallelic.snp.matrix) #individuals before filter
        biallelic.snp.matrix <- biallelic.snp.matrix[(rowMeans(is.na(biallelic.snp.matrix)) <= missing.individuals.cutoff), , drop = FALSE] #filter individuals
        individuals.removed <- individual.count.before.missing.filter - nrow(biallelic.snp.matrix) #number removed
        if (individuals.removed > 0) print.filter.message(individuals.removed, " of ", individual.count.before.missing.filter, " individuals removed due to >", missing.individuals.cutoff * 100, "% missing data") #report
        if (nrow(biallelic.snp.matrix) == 0) stop("All individuals removed after missing data filter") #stop if all gone
      }
      
      # Singleton loci filter
      if (isTRUE(singleton.loci.filter)) {
        locus.count.before.singleton.filter <- ncol(biallelic.snp.matrix) #loci before singleton filter
        called.counts <- colSums(!is.na(biallelic.snp.matrix)) #non-missing calls per locus
        alternate.allele.sums <- colSums(biallelic.snp.matrix, na.rm = TRUE) #alternate allele counts
        minor.allele.counts <- pmin(alternate.allele.sums, called.counts - alternate.allele.sums) #minor allele counts
        keep.loci <- minor.allele.counts != 1L #remove only true singleton loci
        biallelic.snp.matrix <- biallelic.snp.matrix[, keep.loci, drop = FALSE] #filter loci
        loci.removed <- locus.count.before.singleton.filter - ncol(biallelic.snp.matrix) #number removed
        if (loci.removed > 0) print.filter.message(loci.removed, " of ", locus.count.before.singleton.filter, " singleton loci removed") #report
        if (ncol(biallelic.snp.matrix) == 0) stop("All loci removed after singleton filter") #stop if all gone
      }
      
      # Loci missing (strict)
      if (!is.null(missing.loci.cutoff.final)) {
        locus.count.before.strict.missing.filter <- ncol(biallelic.snp.matrix) #loci before filter
        biallelic.snp.matrix <- biallelic.snp.matrix[, (colMeans(is.na(biallelic.snp.matrix)) <= missing.loci.cutoff.final), drop = FALSE] #filter loci
        loci.removed <- locus.count.before.strict.missing.filter - ncol(biallelic.snp.matrix) #number removed
        if (loci.removed > 0) print.filter.message(loci.removed, " of ", locus.count.before.strict.missing.filter, " loci removed due to >", missing.loci.cutoff.final * 100, "% missing data (stricter filter)") #report
        if (ncol(biallelic.snp.matrix) == 0) stop("All loci removed after final missing data filter") #stop if all gone
      }
      
      # Invariant loci filter
      if (isTRUE(invariant.loci.filter)) {
        locus.count.before.final.invariant.filter <- ncol(biallelic.snp.matrix) #loci before final invariant filter
        invariant.loci.reference.count <- locus.count.before.final.invariant.filter #store denominator
        locus.minimum.values <- suppressWarnings(apply(biallelic.snp.matrix, 2, min, na.rm = TRUE)) #minimum per locus
        locus.maximum.values <- suppressWarnings(apply(biallelic.snp.matrix, 2, max, na.rm = TRUE)) #maximum per locus
        keep.loci <- is.finite(locus.minimum.values) & is.finite(locus.maximum.values) & locus.minimum.values != locus.maximum.values #variable loci
        biallelic.snp.matrix <- biallelic.snp.matrix[, keep.loci, drop = FALSE] #remove invariant loci
        invariant.loci.removed.total <- invariant.loci.removed.total + locus.count.before.final.invariant.filter - ncol(biallelic.snp.matrix) #update total removed
        if (invariant.loci.removed.total > 0) print.filter.message(invariant.loci.removed.total, " of ", invariant.loci.reference.count, " invariant loci removed") #report
        if (ncol(biallelic.snp.matrix) == 0) stop("All loci removed after invariant filter") #stop if all gone
      }
      
      # Final summary
      biallelic.snp.matrix <- as.data.frame(biallelic.snp.matrix, stringsAsFactors = FALSE) #convert to data frame
      print.final.summary(biallelic.snp.matrix) #summary
      return(biallelic.snp.matrix) #return binary SNP matrix
    }
    
    # Multiallelic alignment processing
    polymorphic.site.indices <- which(sapply(observed.alleles.per.site, length) > 1) #keep polymorphic sites
    loci.removed <- total.locus.count.before.filtering - length(polymorphic.site.indices) #number removed
    if (loci.removed > 0) print.filter.message(loci.removed, " of ", total.locus.count.before.filtering, " loci removed because they were invariant") #report invariant filter
    if (length(polymorphic.site.indices) == 0) stop("No polymorphic SNPs found in alignment") #stop if none
    alignment.matrix <- alignment.matrix[, polymorphic.site.indices, drop = FALSE] #keep polymorphic sites
    
    # Loci missing (lenient)
    if (!is.null(missing.loci.cutoff.lenient)) {
      locus.count.before.lenient.missing.filter <- ncol(alignment.matrix) #loci before filter
      missing.value.matrix <- alignment.matrix %in% missing.symbols #missing values
      dim(missing.value.matrix) <- dim(alignment.matrix) #restore matrix dimensions
      locus.missing.proportions <- colMeans(missing.value.matrix) #missing proportion per locus
      alignment.matrix <- alignment.matrix[, locus.missing.proportions <= missing.loci.cutoff.lenient, drop = FALSE] #filter loci
      loci.removed <- locus.count.before.lenient.missing.filter - ncol(alignment.matrix) #number removed
      if (loci.removed > 0) print.filter.message(loci.removed, " of ", locus.count.before.lenient.missing.filter, " loci removed due to >", missing.loci.cutoff.lenient * 100, "% missing data (lenient filter)") #report
      if (ncol(alignment.matrix) == 0) stop("All loci removed after lenient missing data filter") #stop if all gone
    }
    
    # Individuals missing
    if (!is.null(missing.individuals.cutoff)) {
      individual.count.before.missing.filter <- nrow(alignment.matrix) #individuals before filter
      missing.value.matrix <- alignment.matrix %in% missing.symbols #missing values
      dim(missing.value.matrix) <- dim(alignment.matrix) #restore matrix dimensions
      individual.missing.proportions <- rowMeans(missing.value.matrix) #missing proportion per individual
      alignment.matrix <- alignment.matrix[individual.missing.proportions <= missing.individuals.cutoff, , drop = FALSE] #filter individuals
      individuals.removed <- individual.count.before.missing.filter - nrow(alignment.matrix) #number removed
      if (individuals.removed > 0) print.filter.message(individuals.removed, " of ", individual.count.before.missing.filter, " individuals removed due to >", missing.individuals.cutoff * 100, "% missing data") #report
      if (nrow(alignment.matrix) == 0) stop("All individuals removed after missing data filter") #stop if all gone
    }
    
    # Singleton loci filter
    if (isTRUE(singleton.loci.filter)) {
      locus.count.before.singleton.filter <- ncol(alignment.matrix) #loci before singleton filter
      singleton.present <- find.singleton.alignment.loci(alignment.matrix = alignment.matrix, missing.symbols = missing.symbols) #singleton loci
      alignment.matrix <- alignment.matrix[, !singleton.present, drop = FALSE] #filter loci
      loci.removed <- locus.count.before.singleton.filter - ncol(alignment.matrix) #number removed
      if (loci.removed > 0) print.filter.message(loci.removed, " of ", locus.count.before.singleton.filter, " singleton loci removed") #report
      if (ncol(alignment.matrix) == 0) stop("All loci removed after singleton filter") #stop if all gone
    }
    
    # Loci missing (strict)
    if (!is.null(missing.loci.cutoff.final)) {
      locus.count.before.strict.missing.filter <- ncol(alignment.matrix) #loci before filter
      missing.value.matrix <- alignment.matrix %in% missing.symbols #missing values
      dim(missing.value.matrix) <- dim(alignment.matrix) #restore matrix dimensions
      locus.missing.proportions <- colMeans(missing.value.matrix) #missing proportion per locus
      alignment.matrix <- alignment.matrix[, locus.missing.proportions <= missing.loci.cutoff.final, drop = FALSE] #filter loci
      loci.removed <- locus.count.before.strict.missing.filter - ncol(alignment.matrix) #number removed
      if (loci.removed > 0) print.filter.message(loci.removed, " of ", locus.count.before.strict.missing.filter, " loci removed due to >", missing.loci.cutoff.final * 100, "% missing data (stricter filter)") #report
      if (ncol(alignment.matrix) == 0) stop("All loci removed after final missing data filter") #stop if all gone
    }
    
    # Invariant loci filter
    if (isTRUE(invariant.loci.filter)) {
      locus.count.before.final.invariant.filter <- ncol(alignment.matrix) #loci before final invariant filter
      final.observed.alleles.per.site <- lapply(seq_len(ncol(alignment.matrix)), function(locus.index) sort(setdiff(unique(alignment.matrix[, locus.index]), missing.symbols))) #recalculate observed alleles
      polymorphic.site.indices <- which(sapply(final.observed.alleles.per.site, length) > 1) #keep polymorphic sites
      alignment.matrix <- alignment.matrix[, polymorphic.site.indices, drop = FALSE] #remove invariant loci
      loci.removed <- locus.count.before.final.invariant.filter - ncol(alignment.matrix) #number removed
      if (loci.removed > 0) print.filter.message(loci.removed, " of ", locus.count.before.final.invariant.filter, " invariant loci removed") #report
      if (ncol(alignment.matrix) == 0) stop("No polymorphic loci remain after filtering") #stop if all gone
    }
    
    # Convert to k-1 dosage matrix
    multiallelic.snp.matrix <- encode.multiallelic.alignment.as.k.minus.one(alignment.matrix = alignment.matrix, missing.symbols = missing.symbols) #encode multiallelic loci
    if (ncol(multiallelic.snp.matrix) == 0) stop("No multiallelic SNP columns remain after filtering") #stop if empty
    
    # Final summary
    print.final.summary(multiallelic.snp.matrix) #summary
    return(multiallelic.snp.matrix) #return multiallelic SNP matrix
  }
  
  # Validate arguments
  validate.arguments() #validate inputs
  
  # Process direct SNP matrix input
  if (!is.null(snp.matrix.input)) {
    return(process.dosage.matrix.input(snp.matrix = snp.matrix.input, input.name = "snp.matrix.input")) #process matrix
  }
  
  # Process genlight input
  if (!is.null(genlight.input)) {
    snp.matrix <- convert.genlight.to.dosage.matrix(genlight.input) #convert genlight
    return(process.dosage.matrix.input(snp.matrix = snp.matrix, input.name = "genlight.input")) #process matrix
  }
  
  # Process PLINK .raw input
  if (!is.null(plink.raw.path)) {
    snp.matrix <- read.plink.raw.input(plink.raw.path) #read PLINK raw
    return(process.dosage.matrix.input(snp.matrix = snp.matrix, input.name = "plink.raw.path")) #process matrix
  }
  
  # Process biallelic diploid VCF directly
  if (!is.null(vcf.path) && isTRUE(make.biallelic)) {
    fast.vcf.result <- tryCatch({
      vcf.object <- suppressWarnings(vcfR::read.vcfR(vcf.path, verbose = FALSE)) #read VCF
      fixed.matrix <- vcf.object@fix #extract fixed fields
      biallelic.variant.indices <- which(!grepl(",", fixed.matrix[, "ALT"]) & !is.na(fixed.matrix[, "ALT"]) & fixed.matrix[, "ALT"] != ".") #biallelic variants
      if (length(biallelic.variant.indices) == 0) stop("No biallelic loci found") #stop if none
      
      loci.removed <- nrow(fixed.matrix) - length(biallelic.variant.indices) #number removed
      if (loci.removed > 0) print.filter.message(loci.removed, " of ", nrow(fixed.matrix), " loci removed because they were not biallelic") #report
      vcf.object <- vcf.object[biallelic.variant.indices, ] #keep biallelic variants
      
      genotype.matrix <- vcfR::extract.gt(vcf.object, element = "GT", as.numeric = FALSE) #variants x samples
      if (is.null(dim(genotype.matrix))) stop("VCF genotype matrix could not be extracted as a matrix") #check matrix
      genotype.matrix <- t(genotype.matrix) #samples x variants
      genotype.dimnames <- dimnames(genotype.matrix) #store dimnames
      genotype.matrix <- matrix(gsub("\\|", "/", genotype.matrix),
                                nrow = nrow(genotype.matrix),
                                ncol = ncol(genotype.matrix),
                                dimnames = genotype.dimnames) #standardize phased genotypes
      
      missing.genotype.codes <- c("./.", ".", "0/.", "./0", "1/.", "./1") #missing and partial-missing genotypes
      recognized.genotype.matrix <- genotype.matrix %in% c("0/0", "0/1", "1/0", "1/1", missing.genotype.codes) | is.na(genotype.matrix) #recognized genotypes
      if (any(!recognized.genotype.matrix)) stop("UNSUPPORTED_FAST_VCF_GENOTYPES") #fallback trigger
      
      snp.matrix <- matrix(NA_integer_,
                           nrow = nrow(genotype.matrix),
                           ncol = ncol(genotype.matrix),
                           dimnames = dimnames(genotype.matrix)) #initialize SNP matrix
      
      snp.matrix[genotype.matrix == "0/0"] <- 0L #homozygous reference
      snp.matrix[genotype.matrix %in% c("0/1", "1/0")] <- 1L #heterozygous
      snp.matrix[genotype.matrix == "1/1"] <- 2L #homozygous alternate
      
      variant.names <- fixed.matrix[biallelic.variant.indices, "ID"] #variant IDs
      missing.variant.names <- is.na(variant.names) | variant.names == "." | variant.names == "" #missing IDs
      variant.names[missing.variant.names] <- paste0("SNP", seq_along(variant.names)[missing.variant.names]) #fallback names
      colnames(snp.matrix) <- make.unique(variant.names) #set unique names
      
      snp.matrix <- filter.SNP.matrix.fast(snp.matrix = snp.matrix,
                                           missing.loci.cutoff.lenient = missing.loci.cutoff.lenient,
                                           missing.loci.cutoff.final = missing.loci.cutoff.final,
                                           missing.individuals.cutoff = missing.individuals.cutoff,
                                           singleton.loci.filter = singleton.loci.filter,
                                           invariant.loci.filter = invariant.loci.filter,
                                           snp.matrix.ploidy = 2L) #filter
      
      print.final.summary(snp.matrix) #summary
      snp.matrix #return result
    }, error = function(fast.vcf.error) {
      if (identical(conditionMessage(fast.vcf.error), "UNSUPPORTED_FAST_VCF_GENOTYPES")) {
        if (!requireNamespace("adegenet", quietly = TRUE)) stop("VCF fast path was skipped because genotype coding is unsupported, and package 'adegenet' is required for fallback genind processing") #check fallback package
        if (verbose) message("VCF fast path skipped: unsupported genotype coding - falling back to genind processing") #report fallback
        return(NULL) #fallback
      }
      stop(fast.vcf.error) #do not hide real filtering/input errors
    })
    
    if (!is.null(fast.vcf.result)) return(fast.vcf.result) #return fast VCF result
  }
  
  # Process NEXUS file
  if (!is.null(nexus.path)) {
    nexus.sequence.list <- ape::read.nexus.data(nexus.path) #read aligned sequences
    return(process.alignment.input(sequence.list = nexus.sequence.list, file.type = "NEXUS")) #process alignment
  }
  
  # Process FASTA file
  if (!is.null(fasta.path)) {
    fasta.matrix <- ape::read.dna(fasta.path, format = "fasta", as.character = TRUE) #read FASTA alignment
    fasta.sequence.list <- convert.alignment.matrix.to.sequence.list(fasta.matrix) #convert to sequence list
    return(process.alignment.input(sequence.list = fasta.sequence.list, file.type = "FASTA")) #process alignment
  }
  
  # Process PHYLIP file
  if (!is.null(phylip.path)) {
    phylip.matrix <- ape::read.dna(phylip.path, format = phylip.format, as.character = TRUE) #read PHYLIP alignment
    phylip.sequence.list <- convert.alignment.matrix.to.sequence.list(phylip.matrix) #convert to sequence list
    return(process.alignment.input(sequence.list = phylip.sequence.list, file.type = "PHYLIP")) #process alignment
  }
  
  # Process genind object or VCF file
  genind.object <- NULL #initialize
  if (!is.null(genind.input)) {
    genind.object <- genind.input #use genind object if provided
  } else {
    vcf.object <- tryCatch({
      suppressWarnings(suppressMessages(vcfR::read.vcfR(vcf.path, verbose = FALSE))) #read VCF
    }, error = function(read.error) stop("VCF could not be read: ", conditionMessage(read.error))) #catch error
    genind.object <- vcfR::vcfR2genind(vcf.object) #convert to genind
  }
  
  # Biallelic loci filter
  if (isTRUE(make.biallelic)) {
    total.locus.count.before.biallelic.filter <- adegenet::nLoc(genind.object) #loci before
    biallelic.locus.indices <- which(genind.object@loc.n.all == 2) #biallelic indices
    if (length(biallelic.locus.indices) == 0) stop("No biallelic loci found") #stop if none
    genind.object <- genind.object[loc = biallelic.locus.indices, drop = TRUE] #keep biallelic
    total.locus.count.after.biallelic.filter <- adegenet::nLoc(genind.object) #loci after
    loci.removed <- total.locus.count.before.biallelic.filter - total.locus.count.after.biallelic.filter #number removed
    if (loci.removed > 0) print.filter.message(loci.removed, " of ", total.locus.count.before.biallelic.filter, " loci removed because they were not biallelic") #report
  }
  
  # Loci missing (lenient)
  if (!is.null(missing.loci.cutoff.lenient)) {
    locus.count.before.lenient.missing.filter <- adegenet::nLoc(genind.object) #loci before filter
    allele.count.matrix.lenient <- suppressMessages(suppressWarnings(adegenet::tab(genind.object, NA.method = "asis"))) #allele count matrix
    locus.factor.lenient <- factor(as.character(genind.object@loc.fac), levels = adegenet::locNames(genind.object)) #locus factor
    first.allele.column.lenient <- vapply(split(seq_along(locus.factor.lenient), locus.factor.lenient, drop = TRUE), `[`, integer(1L), 1L) #first allele column per locus
    locus.missing.proportions.lenient <- colMeans(is.na(allele.count.matrix.lenient[, first.allele.column.lenient, drop = FALSE])) #missing proportion per locus
    names(locus.missing.proportions.lenient) <- names(first.allele.column.lenient) #restore locus names
    genind.object <- genind.object[loc = names(locus.missing.proportions.lenient)[locus.missing.proportions.lenient <= missing.loci.cutoff.lenient], drop = TRUE] #filter loci
    loci.removed <- locus.count.before.lenient.missing.filter - adegenet::nLoc(genind.object) #number removed
    if (loci.removed > 0) print.filter.message(loci.removed, " of ", locus.count.before.lenient.missing.filter, " loci removed due to >", missing.loci.cutoff.lenient * 100, "% missing data (lenient filter)") #report
    if (adegenet::nLoc(genind.object) == 0) stop("All loci removed after lenient missing data filter") #stop if all gone
  }
  
  # Individuals missing
  if (!is.null(missing.individuals.cutoff)) {
    individual.count.before.missing.filter <- adegenet::nInd(genind.object) #individuals before
    allele.count.matrix.individuals <- suppressMessages(suppressWarnings(adegenet::tab(genind.object, NA.method = "asis"))) #allele count matrix
    locus.factor.individuals <- factor(as.character(genind.object@loc.fac), levels = adegenet::locNames(genind.object)) #locus factor
    first.allele.column.individuals <- vapply(split(seq_along(locus.factor.individuals), locus.factor.individuals, drop = TRUE), `[`, integer(1L), 1L) #first allele column per locus
    individual.missing.proportions <- rowMeans(is.na(allele.count.matrix.individuals[, first.allele.column.individuals, drop = FALSE])) #missing proportion per individual
    genind.object <- genind.object[adegenet::indNames(genind.object)[individual.missing.proportions <= missing.individuals.cutoff], drop = TRUE] #filter individuals
    individuals.removed <- individual.count.before.missing.filter - adegenet::nInd(genind.object) #number removed
    if (individuals.removed > 0) print.filter.message(individuals.removed, " of ", individual.count.before.missing.filter, " individuals removed due to >", missing.individuals.cutoff * 100, "% missing data") #report
    if (adegenet::nInd(genind.object) == 0) stop("All individuals removed after missing data filter") #stop if all gone
  }
  
  # Singleton loci filter
  if (isTRUE(singleton.loci.filter)) {
    locus.count.before.singleton.filter <- adegenet::nLoc(genind.object) #loci before singleton filter
    minor.allele.counts <- count.minor.alleles.per.locus(genind.object) #minor allele counts per locus
    singleton.locus.indices <- which(minor.allele.counts == 1L) #true singleton loci only
    if (length(singleton.locus.indices) > 0) {
      genind.object <- genind.object[loc = -singleton.locus.indices, drop = TRUE] #remove them
      print.filter.message(length(singleton.locus.indices), " of ", locus.count.before.singleton.filter, " singleton loci removed") #report
      if (adegenet::nLoc(genind.object) == 0) stop("All loci removed after singleton filter") #stop if all gone
    }
  }
  
  # Loci missing (strict)
  if (!is.null(missing.loci.cutoff.final)) {
    locus.count.before.strict.missing.filter <- adegenet::nLoc(genind.object) #loci before
    allele.count.matrix.strict <- suppressMessages(suppressWarnings(adegenet::tab(genind.object, NA.method = "asis"))) #allele count matrix
    locus.factor.strict <- factor(as.character(genind.object@loc.fac), levels = adegenet::locNames(genind.object)) #locus factor
    first.allele.column.strict <- vapply(split(seq_along(locus.factor.strict), locus.factor.strict, drop = TRUE), `[`, integer(1L), 1L) #first allele column per locus
    locus.missing.proportions.strict <- colMeans(is.na(allele.count.matrix.strict[, first.allele.column.strict, drop = FALSE])) #missing proportion per locus
    names(locus.missing.proportions.strict) <- names(first.allele.column.strict) #restore locus names
    genind.object <- genind.object[loc = names(locus.missing.proportions.strict)[locus.missing.proportions.strict <= missing.loci.cutoff.final], drop = TRUE] #filter loci
    loci.removed <- locus.count.before.strict.missing.filter - adegenet::nLoc(genind.object) #number removed
    if (loci.removed > 0) print.filter.message(loci.removed, " of ", locus.count.before.strict.missing.filter, " loci removed due to >", missing.loci.cutoff.final * 100, "% missing data (stricter filter)") #report
    if (adegenet::nLoc(genind.object) == 0) stop("All loci removed after final missing data filter") #stop if all gone
  }
  
  invariant.loci.removed.total <- 0L #track total invariant loci removed
  invariant.loci.reference.count <- adegenet::nLoc(genind.object) #track denominator for invariant message
  
  # Invariant loci filter
  if (isTRUE(invariant.loci.filter)) {
    observed.allele.counts <- count.observed.alleles.per.locus(genind.object) #count observed alleles per locus
    final.invariant.locus.indices <- which(observed.allele.counts <= 1) #final invariant check
    if (length(final.invariant.locus.indices) > 0) {
      genind.object <- genind.object[loc = -final.invariant.locus.indices, drop = TRUE] #remove them
      invariant.loci.removed.total <- invariant.loci.removed.total + length(final.invariant.locus.indices) #update total removed
      if (adegenet::nLoc(genind.object) == 0) stop("All loci removed after invariant filter") #stop if all gone
    }
    if (invariant.loci.removed.total > 0) print.filter.message(invariant.loci.removed.total, " of ", invariant.loci.reference.count, " invariant loci removed") #report
  }
  
  # Convert to final SNP matrix
  if (isTRUE(make.biallelic)) {
    allele.count.matrix <- suppressMessages(suppressWarnings(adegenet::tab(genind.object, NA.method = "asis"))) #allele count matrix
    locus.names <- adegenet::locNames(genind.object) #locus names
    locus.membership.vector <- factor(as.character(genind.object@loc.fac), levels = locus.names) #locus membership
    locus.column.list <- split(seq_along(locus.membership.vector), locus.membership.vector, drop = TRUE) #columns per locus
    allele.present <- colSums(allele.count.matrix, na.rm = TRUE) > 0 #observed allele columns
    
    retained.info <- lapply(names(locus.column.list), function(current.locus.name) {
      locus.column.indices <- locus.column.list[[current.locus.name]] #columns for locus
      observed.column.indices <- locus.column.indices[allele.present[locus.column.indices]] #observed allele columns
      if (length(observed.column.indices) < 2) return(NULL) #skip invariant loci
      return(list(index = observed.column.indices[length(observed.column.indices)], name = current.locus.name)) #retain one allele dosage column
    }) #retained column information
    retained.info <- Filter(Negate(is.null), retained.info) #remove skipped loci
    
    if (length(retained.info) == 0) stop("No biallelic SNP columns remain after filtering") #stop if empty
    retained.column.indices <- vapply(retained.info, function(single.info) single.info$index, integer(1)) #retained column indices
    retained.allele.dosage.names <- vapply(retained.info, function(single.info) single.info$name, character(1)) #retained locus names
    
    biallelic.snp.matrix <- as.data.frame(allele.count.matrix[, retained.column.indices, drop = FALSE], stringsAsFactors = FALSE) #subset once
    rownames(biallelic.snp.matrix) <- adegenet::indNames(genind.object) #set rownames
    colnames(biallelic.snp.matrix) <- retained.allele.dosage.names #set colnames
    
    print.final.summary(biallelic.snp.matrix) #summary
    return(biallelic.snp.matrix) #return SNP matrix
  }
  
  multiallelic.snp.matrix <- convert.genind.to.k.minus.one.dosage(genind.object) #convert to k-1 dosage matrix
  if (ncol(multiallelic.snp.matrix) == 0) stop("No multiallelic SNP columns remain after filtering") #stop if empty
  print.final.summary(multiallelic.snp.matrix) #summary
  return(multiallelic.snp.matrix) #return multiallelic SNP matrix
}

# Function to summarize and plot eta squared and/or map variance across SOM layers
plot.layer.importance.varimp.SOM <- function(SOM.output, #clustered SOM output from clustering.SOM
                                             col.pal = viridis::turbo, #color palette as in plot.layers.SOM
                                             save = FALSE, #option to save plot
                                             overwrite = TRUE, #option to overwrite plot if it already exists
                                             plot.type = "svg", #plot type options: "svg", "png", "jpg"
                                             file.name = NULL, #set file.name for saving (if NULL, default plot file.name is used)
                                             width = 20, #plot width in cm
                                             height = 12, #plot height in cm
                                             resolution = 300, #plot resolution in dpi
                                             bottom.margin = 4, #bottom margin
                                             left.margin = 5, #left margin
                                             top.margin = 3, #top margin
                                             right.margin = 2.5, #right margin
                                             etasquared.title = "Eta squared across layers", #title of eta plot
                                             mapvariance.title = "Map variance across layers", #title of map variance plot
                                             etasquared.y.axis.label = "Eta squared", #y axis label for eta plot
                                             mapvariance.y.axis.label = "Map variance", #y axis label for map variance plot
                                             title.font.size = 1.2, #font size of plot titles
                                             axis.font.size = 0.9, #font size of axis labels
                                             add.boxplot.whiskers = TRUE, #whether to show boxplot whiskers
                                             sort.by.median = TRUE, #whether to sort layers by median importance
                                             verbose = TRUE #whether to print messages
) {
  
  # Set messages
  messager <- function(...) if (isTRUE(verbose)) message(...)
  
  # Reset plotting parameters
  old_graphics_device <- dev.cur()
  old_plotting_parameters <- par(no.readonly = TRUE)
  on.exit({
    if (dev.cur() == old_graphics_device) par(old_plotting_parameters)
  }, add = TRUE)
  
  # Validate specified SOM.output
  if (is.null(SOM.output) || !is.list(SOM.output)) {
    stop("Plotting aborted: SOM.output must be a non-NULL list")
  }
  
  # Validate specified save
  if (!is.logical(save) || length(save) != 1 || is.na(save)) {
    stop("Plotting aborted: save must be TRUE or FALSE")
  }
  
  # Validate specified overwrite
  if (!is.logical(overwrite) || length(overwrite) != 1 || is.na(overwrite)) {
    stop("Plotting aborted: overwrite must be TRUE or FALSE")
  }
  
  # Validate specified plot.type
  if (!is.character(plot.type) || length(plot.type) != 1 || !(plot.type %in% c("svg", "png", "jpg"))) {
    stop("Plotting aborted: plot.type must be 'svg', 'png', or 'jpg'")
  }
  
  # Validate specified title.font.size
  if (!is.numeric(title.font.size) || length(title.font.size) != 1 || is.na(title.font.size) || title.font.size <= 0) {
    stop("Plotting aborted: title.font.size must be a single positive numeric value")
  }
  
  # Validate specified axis.font.size
  if (!is.numeric(axis.font.size) || length(axis.font.size) != 1 || is.na(axis.font.size) || axis.font.size <= 0) {
    stop("Plotting aborted: axis.font.size must be a single positive numeric value")
  }
  
  # Validate specified add.boxplot.whiskers
  if (!is.logical(add.boxplot.whiskers) || length(add.boxplot.whiskers) != 1 || is.na(add.boxplot.whiskers)) {
    stop("Plotting aborted: add.boxplot.whiskers must be TRUE or FALSE")
  }
  
  # Validate specified sort.by.median
  if (!is.logical(sort.by.median) || length(sort.by.median) != 1 || is.na(sort.by.median)) {
    stop("Plotting aborted: sort.by.median must be TRUE or FALSE")
  }
  
  # Validate specified verbose
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("Plotting aborted: verbose must be TRUE or FALSE")
  }
  
  # Extract layer names
  SOM_layer_names <- NULL
  if (!is.null(SOM.output$input_data_names)) {
    SOM_layer_names <- as.character(SOM.output$input_data_names)
  }
  if (is.null(SOM_layer_names) && !is.null(SOM.output$distance_weights_matrix)) {
    SOM_layer_names <- colnames(SOM.output$distance_weights_matrix)
  }
  if (is.null(SOM_layer_names) || length(SOM_layer_names) == 0) {
    stop("Plotting aborted: layer names could not be determined from SOM.output")
  }
  
  # Extract eta squared and map variance lists if present
  eta_squared_variable_importance_list <- NULL
  map_variance_variable_importance_list <- NULL
  if (!is.null(SOM.output$median_etasquared_variable_importance)) {
    eta_squared_variable_importance_list <- SOM.output$median_etasquared_variable_importance
    if (is.null(names(eta_squared_variable_importance_list)) && length(eta_squared_variable_importance_list) == length(SOM_layer_names)) {
      names(eta_squared_variable_importance_list) <- SOM_layer_names
    }
    if (!is.null(names(eta_squared_variable_importance_list)) && all(SOM_layer_names %in% names(eta_squared_variable_importance_list))) {
      eta_squared_variable_importance_list <- eta_squared_variable_importance_list[SOM_layer_names]
    }
    eta_squared_variable_importance_list <- lapply(eta_squared_variable_importance_list, function(variable_importance_values) {
      variable_importance_values <- as.numeric(variable_importance_values)
      variable_importance_values <- variable_importance_values[is.finite(variable_importance_values) & !is.na(variable_importance_values)]
      return(variable_importance_values)
    })
  }
  if (!is.null(SOM.output$median_map_variance_variable_importance)) {
    map_variance_variable_importance_list <- SOM.output$median_map_variance_variable_importance
    if (is.null(names(map_variance_variable_importance_list)) && length(map_variance_variable_importance_list) == length(SOM_layer_names)) {
      names(map_variance_variable_importance_list) <- SOM_layer_names
    }
    if (!is.null(names(map_variance_variable_importance_list)) && all(SOM_layer_names %in% names(map_variance_variable_importance_list))) {
      map_variance_variable_importance_list <- map_variance_variable_importance_list[SOM_layer_names]
    }
    map_variance_variable_importance_list <- lapply(map_variance_variable_importance_list, function(variable_importance_values) {
      variable_importance_values <- as.numeric(variable_importance_values)
      variable_importance_values <- variable_importance_values[is.finite(variable_importance_values) & !is.na(variable_importance_values)]
      return(variable_importance_values)
    })
  }
  
  # Check whether there are valid values
  eta_squared_available <- FALSE
  map_variance_available <- FALSE
  if (!is.null(eta_squared_variable_importance_list)) {
    eta_squared_available <- any(vapply(eta_squared_variable_importance_list, length, numeric(1)) > 0)
  }
  if (!is.null(map_variance_variable_importance_list)) {
    map_variance_available <- any(vapply(map_variance_variable_importance_list, length, numeric(1)) > 0)
  }
  
  # Report informative messages for special cases
  if (length(SOM_layer_names) == 1) {
    messager("Only one SOM layer was detected")
  }
  
  if (!eta_squared_available) {
    if (!is.null(SOM.output$optim_k_vals) && length(SOM.output$optim_k_vals) > 0 && all(SOM.output$optim_k_vals == 1, na.rm = TRUE)) {
      messager("Eta squared could not be calculated because all retained SOM replicates had k = 1")
    } else if (length(SOM_layer_names) == 1) {
      messager("Eta squared could not be calculated because only one layer was detected")
    }
  }
  
  if (!eta_squared_available && !map_variance_available) {
    stop("Plotting aborted: no valid eta squared or map variance values found across layers")
  }
  
  # Calculate summary table
  layer_importance_summary_table <- data.frame(layer = SOM_layer_names,
                                               eta.squared.mean = NA_real_,
                                               eta.squared.sd = NA_real_,
                                               map.variance.mean = NA_real_,
                                               map.variance.sd = NA_real_,
                                               stringsAsFactors = FALSE)
  
  for (layer_index in seq_along(SOM_layer_names)) {
    
    # Extract layer name
    current_layer_name <- SOM_layer_names[layer_index]
    
    # Extract eta squared values
    if (eta_squared_available) {
      current_eta_squared_values <- NULL
      if (!is.null(names(eta_squared_variable_importance_list)) && current_layer_name %in% names(eta_squared_variable_importance_list)) {
        current_eta_squared_values <- eta_squared_variable_importance_list[[current_layer_name]]
      } else if (length(eta_squared_variable_importance_list) >= layer_index) {
        current_eta_squared_values <- eta_squared_variable_importance_list[[layer_index]]
      }
      if (!is.null(current_eta_squared_values) && length(current_eta_squared_values) > 0) {
        layer_importance_summary_table$eta.squared.mean[layer_index] <- mean(current_eta_squared_values, na.rm = TRUE)
        layer_importance_summary_table$eta.squared.sd[layer_index] <- ifelse(length(current_eta_squared_values) > 1, stats::sd(current_eta_squared_values, na.rm = TRUE), 0)
      }
    }
    
    # Extract map variance values
    if (map_variance_available) {
      current_map_variance_values <- NULL
      if (!is.null(names(map_variance_variable_importance_list)) && current_layer_name %in% names(map_variance_variable_importance_list)) {
        current_map_variance_values <- map_variance_variable_importance_list[[current_layer_name]]
      } else if (length(map_variance_variable_importance_list) >= layer_index) {
        current_map_variance_values <- map_variance_variable_importance_list[[layer_index]]
      }
      if (!is.null(current_map_variance_values) && length(current_map_variance_values) > 0) {
        layer_importance_summary_table$map.variance.mean[layer_index] <- mean(current_map_variance_values, na.rm = TRUE)
        layer_importance_summary_table$map.variance.sd[layer_index] <- ifelse(length(current_map_variance_values) > 1, stats::sd(current_map_variance_values, na.rm = TRUE), 0)
      }
    }
  }
  
  # Create sorted plot lists if requested
  eta_squared_plot_list <- eta_squared_variable_importance_list
  map_variance_plot_list <- map_variance_variable_importance_list
  
  if (sort.by.median) {
    if (eta_squared_available) {
      eta_squared_layer_medians <- vapply(eta_squared_plot_list, function(x) stats::median(x, na.rm = TRUE), numeric(1))
      eta_squared_layer_medians[!is.finite(eta_squared_layer_medians) | is.na(eta_squared_layer_medians)] <- -Inf
      eta_squared_plot_list <- eta_squared_plot_list[names(sort(eta_squared_layer_medians, decreasing = TRUE))]
    }
    if (map_variance_available) {
      map_variance_layer_medians <- vapply(map_variance_plot_list, function(x) stats::median(x, na.rm = TRUE), numeric(1))
      map_variance_layer_medians[!is.finite(map_variance_layer_medians) | is.na(map_variance_layer_medians)] <- -Inf
      map_variance_plot_list <- map_variance_plot_list[names(sort(map_variance_layer_medians, decreasing = TRUE))]
    }
  }
  
  # Set file name
  if (is.null(file.name)) {
    if (eta_squared_available && map_variance_available) {
      file.name <- paste0("SOM_variable_importance_layers_both_", paste(SOM_layer_names, collapse = "_"), ".", plot.type)
    }
    if (eta_squared_available && !map_variance_available) {
      file.name <- paste0("SOM_variable_importance_layers_eta_squared_", paste(SOM_layer_names, collapse = "_"), ".", plot.type)
    }
    if (!eta_squared_available && map_variance_available) {
      file.name <- paste0("SOM_variable_importance_layers_map_variance_", paste(SOM_layer_names, collapse = "_"), ".", plot.type)
    }
  }
  
  # Check for file existence if overwrite is FALSE
  if (!overwrite && file.exists(file.name)) {
    stop(file.name, " already exists - skipping plot saving")
  }
  
  # Define color palette for layers
  layer_colors <- setNames(col.pal(length(SOM_layer_names)), SOM_layer_names)
  
  # Save plot if requested
  if (save) {
    if (eta_squared_available && map_variance_available) {
      plot_width_to_use <- width
    } else {
      plot_width_to_use <- width / 2
    }
    if (plot.type == "svg") {
      svg(file.name,
          width = plot_width_to_use / 2.54,
          height = height / 2.54)
    } else if (plot.type == "png") {
      png(file.name,
          width = plot_width_to_use,
          height = height,
          units = "cm",
          res = resolution)
    } else if (plot.type == "jpg") {
      jpeg(file.name,
           width = plot_width_to_use,
           height = height,
           units = "cm",
           res = resolution)
    }
  }
  
  # Plot both panels
  if (eta_squared_available && map_variance_available) {
    par(mfrow = c(1, 2),
        mar = c(bottom.margin,
                left.margin,
                top.margin,
                right.margin),
        mgp = c(2.5, 0.8, 0),
        xpd = FALSE)
    
    # Calculate eta squared boxplot object and y-axis limits
    eta_squared_boxplot_object <- boxplot(eta_squared_plot_list, plot = FALSE, outline = FALSE)
    eta_squared_y_axis_limits <- range(eta_squared_boxplot_object$stats, na.rm = TRUE)
    if (diff(eta_squared_y_axis_limits) == 0) {
      eta_squared_axis_padding <- max(abs(eta_squared_y_axis_limits[1]) * 0.05, 0.01)
      eta_squared_y_axis_limits <- c(eta_squared_y_axis_limits[1] - eta_squared_axis_padding,
                                     eta_squared_y_axis_limits[2] + eta_squared_axis_padding)
    } else {
      eta_squared_axis_padding <- diff(eta_squared_y_axis_limits) * 0.05
      eta_squared_y_axis_limits <- c(eta_squared_y_axis_limits[1] - eta_squared_axis_padding,
                                     eta_squared_y_axis_limits[2] + eta_squared_axis_padding)
    }
    
    # Calculate map variance boxplot object and y-axis limits
    map_variance_boxplot_object <- boxplot(map_variance_plot_list, plot = FALSE, outline = FALSE)
    map_variance_y_axis_limits <- range(map_variance_boxplot_object$stats, na.rm = TRUE)
    if (diff(map_variance_y_axis_limits) == 0) {
      map_variance_axis_padding <- max(abs(map_variance_y_axis_limits[1]) * 0.05, 0.01)
      map_variance_y_axis_limits <- c(map_variance_y_axis_limits[1] - map_variance_axis_padding,
                                      map_variance_y_axis_limits[2] + map_variance_axis_padding)
    } else {
      map_variance_axis_padding <- diff(map_variance_y_axis_limits) * 0.05
      map_variance_y_axis_limits <- c(map_variance_y_axis_limits[1] - map_variance_axis_padding,
                                      map_variance_y_axis_limits[2] + map_variance_axis_padding)
    }
    
    # Plot eta squared boxplots
    boxplot(eta_squared_plot_list,
            names = FALSE,
            col = layer_colors[names(eta_squared_plot_list)],
            las = 1,
            axes = FALSE,
            outline = FALSE,
            ylab = etasquared.y.axis.label,
            main = etasquared.title,
            cex.main = title.font.size,
            cex.lab = axis.font.size,
            ylim = eta_squared_y_axis_limits,
            whisklty = ifelse(add.boxplot.whiskers, 1, 0),
            staplelty = ifelse(add.boxplot.whiskers, 1, 0))
    axis(1,
         at = seq_along(eta_squared_plot_list),
         labels = names(eta_squared_plot_list),
         las = 2,
         tick = FALSE,
         line = -0.5,
         cex.axis = axis.font.size)
    axis(2,
         cex.axis = axis.font.size)
    box()
    
    # Plot map variance boxplots
    boxplot(map_variance_plot_list,
            names = FALSE,
            col = layer_colors[names(map_variance_plot_list)],
            las = 1,
            axes = FALSE,
            outline = FALSE,
            ylab = mapvariance.y.axis.label,
            main = mapvariance.title,
            cex.main = title.font.size,
            cex.lab = axis.font.size,
            ylim = map_variance_y_axis_limits,
            whisklty = ifelse(add.boxplot.whiskers, 1, 0),
            staplelty = ifelse(add.boxplot.whiskers, 1, 0))
    axis(1,
         at = seq_along(map_variance_plot_list),
         labels = names(map_variance_plot_list),
         las = 2,
         tick = FALSE,
         line = -0.5,
         cex.axis = axis.font.size)
    axis(2,
         cex.axis = axis.font.size)
    box()
  }
  
  # Plot eta squared only
  if (eta_squared_available && !map_variance_available) {
    par(mfrow = c(1, 1),
        mar = c(bottom.margin,
                left.margin,
                top.margin,
                right.margin),
        mgp = c(2.5, 0.8, 0),
        xpd = FALSE)
    eta_squared_boxplot_object <- boxplot(eta_squared_plot_list, plot = FALSE, outline = FALSE)
    eta_squared_y_axis_limits <- range(eta_squared_boxplot_object$stats, na.rm = TRUE)
    if (diff(eta_squared_y_axis_limits) == 0) {
      eta_squared_axis_padding <- max(abs(eta_squared_y_axis_limits[1]) * 0.05, 0.01)
      eta_squared_y_axis_limits <- c(eta_squared_y_axis_limits[1] - eta_squared_axis_padding,
                                     eta_squared_y_axis_limits[2] + eta_squared_axis_padding)
    } else {
      eta_squared_axis_padding <- diff(eta_squared_y_axis_limits) * 0.05
      eta_squared_y_axis_limits <- c(eta_squared_y_axis_limits[1] - eta_squared_axis_padding,
                                     eta_squared_y_axis_limits[2] + eta_squared_axis_padding)
    }
    boxplot(eta_squared_plot_list,
            names = FALSE,
            col = layer_colors[names(eta_squared_plot_list)],
            las = 1,
            axes = FALSE,
            outline = FALSE,
            ylab = etasquared.y.axis.label,
            main = etasquared.title,
            cex.main = title.font.size,
            cex.lab = axis.font.size,
            ylim = eta_squared_y_axis_limits,
            whisklty = ifelse(add.boxplot.whiskers, 1, 0),
            staplelty = ifelse(add.boxplot.whiskers, 1, 0))
    axis(1,
         at = seq_along(eta_squared_plot_list),
         labels = names(eta_squared_plot_list),
         las = 2,
         tick = FALSE,
         line = -0.5,
         cex.axis = axis.font.size)
    axis(2,
         cex.axis = axis.font.size)
    box()
  }
  
  # Plot map variance only
  if (!eta_squared_available && map_variance_available) {
    par(mfrow = c(1, 1),
        mar = c(bottom.margin,
                left.margin,
                top.margin,
                right.margin),
        mgp = c(2.5, 0.8, 0),
        xpd = FALSE)
    map_variance_boxplot_object <- boxplot(map_variance_plot_list, plot = FALSE, outline = FALSE)
    map_variance_y_axis_limits <- range(map_variance_boxplot_object$stats, na.rm = TRUE)
    if (diff(map_variance_y_axis_limits) == 0) {
      map_variance_axis_padding <- max(abs(map_variance_y_axis_limits[1]) * 0.05, 0.01)
      map_variance_y_axis_limits <- c(map_variance_y_axis_limits[1] - map_variance_axis_padding,
                                      map_variance_y_axis_limits[2] + map_variance_axis_padding)
    } else {
      map_variance_axis_padding <- diff(map_variance_y_axis_limits) * 0.05
      map_variance_y_axis_limits <- c(map_variance_y_axis_limits[1] - map_variance_axis_padding,
                                      map_variance_y_axis_limits[2] + map_variance_axis_padding)
    }
    boxplot(map_variance_plot_list,
            names = FALSE,
            col = layer_colors[names(map_variance_plot_list)],
            las = 1,
            axes = FALSE,
            outline = FALSE,
            ylab = mapvariance.y.axis.label,
            main = mapvariance.title,
            cex.main = title.font.size,
            cex.lab = axis.font.size,
            ylim = map_variance_y_axis_limits,
            whisklty = ifelse(add.boxplot.whiskers, 1, 0),
            staplelty = ifelse(add.boxplot.whiskers, 1, 0))
    axis(1,
         at = seq_along(map_variance_plot_list),
         labels = names(map_variance_plot_list),
         las = 2,
         tick = FALSE,
         line = -0.5,
         cex.axis = axis.font.size)
    axis(2,
         cex.axis = axis.font.size)
    box()
  }
  
  # Close graphics device
  if (save) {
    dev.off()
    messager(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
  
  return(layer_importance_summary_table)
}



# Create function to estimate and plot layer importance by replicate-matched leave-one-layer-out SOM analyses
plot.layer.importance.leaveoneout.SOM <- function(SOM_output, #clustered SOM output from full multilayer input
                                                  col.pal = viridis::turbo, #color palette
                                                  add.points = TRUE, #whether to add replicate-level jittered points
                                                  point.cex = 0.7, #point size
                                                  point.alpha = 0.5, #point transparency
                                                  bottom.margin = 8, #bottom margin
                                                  left.margin = 5.5, #left margin
                                                  top.margin = 3, #top margin
                                                  right.margin = 2, #right margin
                                                  distance.axis.label = 4, #distance of axis title from axis labels in par(mgp = c(...))
                                                  save.leave.one.layer.out.results = FALSE, #whether to save leave-one-layer-out results to file
                                                  save.leave.one.layer.out.results.name = NULL, #file name for saving leave-one-layer-out results; if NULL, default name is generated
                                                  overwrite.leave.one.layer.out.results = FALSE, #if FALSE, existing results are loaded instead of re-running
                                                  save = FALSE, #whether to save plot
                                                  overwrite = TRUE, #whether to overwrite plot if it already exists
                                                  plot.type = "svg", #plot file type (options: "svg", "png", "jpg")
                                                  file.name = NULL, #plot file name (if NULL, default file name is used)
                                                  width = 20, #plot width in cm
                                                  height = 15, #plot height in cm
                                                  resolution = 300, #plot resolution in dpi
                                                  title = NULL, #optional overall plot title
                                                  message.N.replicates = 10, #frequency of progress messages during leave-one-layer-out SOM reruns
                                                  verbose = TRUE #whether to print messages
) {
  
  # Set messages
  messager <- function(...) if (isTRUE(verbose)) message(...)
  
  # Validate specified SOM_output
  if (is.null(SOM_output) || !is.list(SOM_output)) {
    stop("Leave-one-layer-out layer importance aborted: SOM_output must be a non-NULL list")
  }
  if (is.null(SOM_output$som_models) || length(SOM_output$som_models) < 1) {
    stop("Leave-one-layer-out layer importance aborted: SOM_output must contain non-empty som_models")
  }
  if (is.null(SOM_output$cluster_assignment)) {
    stop("Leave-one-layer-out layer importance aborted: SOM_output must contain cluster_assignment")
  }
  if (is.null(SOM_output$optim_k_vals)) {
    stop("Leave-one-layer-out layer importance aborted: SOM_output must contain optim_k_vals")
  }
  
  # Recover input_data from SOM_output
  if (!is.null(SOM_output$input_data)) {
    input_data <- SOM_output$input_data
  } else {
    stop("Leave-one-layer-out layer importance aborted: SOM_output does not contain stored input_data")
  }
  
  # Recover baseline training seed from SOM_output
  if (!is.null(SOM_output$train.SOM.set.seed.N)) {
    baseline.train.SOM.set.seed.N <- SOM_output$train.SOM.set.seed.N
  } else {
    stop("Leave-one-layer-out layer importance aborted: SOM_output does not contain train.SOM.set.seed.N (required for replicate-matching)")
  }
  
  # Recover baseline clustering seed from SOM_output
  if (!is.null(SOM_output$clustering.SOM.set.seed.N)) {
    baseline.clustering.SOM.set.seed.N <- SOM_output$clustering.SOM.set.seed.N
  } else {
    stop("Leave-one-layer-out layer importance aborted: SOM_output does not contain clustering.SOM.set.seed.N (required for replicate-matching)")
  }
  
  # Recover train.SOM args from SOM_output
  if (!is.null(SOM_output$train.SOM.args)) {
    train.SOM.args <- SOM_output$train.SOM.args
  } else {
    train.SOM.args <- list()
  }
  
  # Recover clustering.SOM args from SOM_output
  if (!is.null(SOM_output$clustering.SOM.args)) {
    clustering.SOM.args <- SOM_output$clustering.SOM.args
  } else {
    clustering.SOM.args <- list()
  }
  
  # Remove arguments that must be overridden for leave-one-layer-out reruns
  train.SOM.args$set.seed.N <- NULL
  train.SOM.args$N.replicates <- NULL
  train.SOM.args$N.replicates.requested <- NULL
  train.SOM.args$N.replicates.failed <- NULL
  train.SOM.args$parallel <- NULL
  train.SOM.args$save.SOM.results <- NULL
  train.SOM.args$save.SOM.results.name <- NULL
  train.SOM.args$overwrite.SOM.results <- NULL
  clustering.SOM.args$set.seed.N <- NULL
  clustering.SOM.args$quantization.error.quantile <- NULL
  clustering.SOM.args$topographic.error.quantile <- NULL
  
  # Set or validate message frequency for leave-one-layer-out progress messages
  if (is.null(message.N.replicates)) {
    if (is.null(train.SOM.args$message.N.replicates)) {
      train.SOM.args$message.N.replicates <- 20
    }
  } else {
    if (!is.numeric(message.N.replicates) || length(message.N.replicates) != 1 || is.na(message.N.replicates) ||
        message.N.replicates < 1 || (message.N.replicates %% 1 != 0)) {
      stop("Leave-one-layer-out layer importance aborted: message.N.replicates must be NULL or a single positive integer (>= 1)")
    }
    train.SOM.args$message.N.replicates <- message.N.replicates
  }
  
  # Validate specified input_data
  if (!is.list(input_data) || length(input_data) < 2) {
    stop("Plotting aborted: function requires at least two layers")
  }
  if (is.null(names(input_data)) || any(names(input_data) == "")) {
    names(input_data) <- paste0("Layer_", seq_along(input_data))
  }
  
  # Convert input_data to matrices and validate row names
  input_data <- lapply(input_data, function(input_layer_matrix) {
    input_layer_matrix <- as.matrix(input_layer_matrix)
    if (is.null(rownames(input_layer_matrix))) {
      stop("Leave-one-layer-out layer importance aborted: all input layers must have rownames")
    }
    return(input_layer_matrix)
  })
  
  # Validate specified seeds
  if (!is.numeric(baseline.train.SOM.set.seed.N) || length(baseline.train.SOM.set.seed.N) != 1 || is.na(baseline.train.SOM.set.seed.N) || baseline.train.SOM.set.seed.N < 1 || baseline.train.SOM.set.seed.N %% 1 != 0) {
    stop("Leave-one-layer-out layer importance aborted: baseline.train.SOM.set.seed.N must be a single positive integer")
  }
  if (!is.numeric(baseline.clustering.SOM.set.seed.N) || length(baseline.clustering.SOM.set.seed.N) != 1 || is.na(baseline.clustering.SOM.set.seed.N) || baseline.clustering.SOM.set.seed.N < 1 || baseline.clustering.SOM.set.seed.N %% 1 != 0) {
    stop("Leave-one-layer-out layer importance aborted: baseline.clustering.SOM.set.seed.N must be a single positive integer")
  }
  
  # Validate specified add.points
  if (!is.logical(add.points) || length(add.points) != 1 || is.na(add.points)) {
    stop("Leave-one-layer-out layer importance aborted: add.points must be TRUE or FALSE")
  }
  
  # Validate specified distance.axis.label
  if (!is.numeric(distance.axis.label) || length(distance.axis.label) != 1 || is.na(distance.axis.label) || distance.axis.label <= 0) {
    stop("Leave-one-layer-out layer importance aborted: distance.axis.label must be a single positive numeric value")
  }
  
  # Validate specified save.leave.one.layer.out.results
  if (!is.logical(save.leave.one.layer.out.results) || length(save.leave.one.layer.out.results) != 1 || is.na(save.leave.one.layer.out.results)) {
    stop("Leave-one-layer-out layer importance aborted: save.leave.one.layer.out.results must be TRUE or FALSE")
  }
  
  # Validate specified save.leave.one.layer.out.results.name
  if (save.leave.one.layer.out.results && !is.null(save.leave.one.layer.out.results.name)) {
    if (!is.character(save.leave.one.layer.out.results.name) || length(save.leave.one.layer.out.results.name) != 1 || is.na(save.leave.one.layer.out.results.name) || trimws(save.leave.one.layer.out.results.name) == "") {
      stop("Leave-one-layer-out layer importance aborted: save.leave.one.layer.out.results.name must be a non-empty character string if provided")
    }
    if (tolower(tools::file_ext(save.leave.one.layer.out.results.name)) != "rdata") {
      stop("Leave-one-layer-out layer importance aborted: save.leave.one.layer.out.results.name must end with '.Rdata'")
    }
  }
  
  # Validate specified overwrite.leave.one.layer.out.results
  if (!is.logical(overwrite.leave.one.layer.out.results) || length(overwrite.leave.one.layer.out.results) != 1 || is.na(overwrite.leave.one.layer.out.results)) {
    stop("Leave-one-layer-out layer importance aborted: overwrite.leave.one.layer.out.results must be TRUE or FALSE")
  }
  
  # Validate specified save
  if (!is.logical(save) || length(save) != 1 || is.na(save)) {
    stop("Leave-one-layer-out layer importance aborted: save must be TRUE or FALSE")
  }
  
  # Validate specified overwrite
  if (!is.logical(overwrite) || length(overwrite) != 1 || is.na(overwrite)) {
    stop("Leave-one-layer-out layer importance aborted: overwrite must be TRUE or FALSE")
  }
  
  # Validate specified plot.type
  if (!is.character(plot.type) || length(plot.type) != 1 || is.na(plot.type) || !(tolower(plot.type) %in% c("svg", "png", "jpg"))) {
    stop("Leave-one-layer-out layer importance aborted: plot.type must be one of 'svg', 'png', or 'jpg'")
  }
  plot.type <- tolower(plot.type)
  
  # Validate specified file.name
  if (!is.null(file.name)) {
    if (!is.character(file.name) || length(file.name) != 1 || is.na(file.name) || trimws(file.name) == "") {
      stop("Leave-one-layer-out layer importance aborted: file.name must be NULL or a non-empty character string")
    }
  }
  
  # Validate specified width
  if (!is.numeric(width) || length(width) != 1 || is.na(width) || width <= 0) {
    stop("Leave-one-layer-out layer importance aborted: width must be a single positive numeric value")
  }
  
  # Validate specified height
  if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
    stop("Leave-one-layer-out layer importance aborted: height must be a single positive numeric value")
  }
  
  # Validate specified resolution
  if (!is.numeric(resolution) || length(resolution) != 1 || is.na(resolution) || resolution <= 0) {
    stop("Leave-one-layer-out layer importance aborted: resolution must be a single positive numeric value")
  }
  
  # Validate specified title
  if (!is.null(title)) {
    if (!is.character(title) || length(title) != 1 || is.na(title)) {
      stop("Leave-one-layer-out layer importance aborted: title must be NULL or a character string of length 1")
    }
  }
  
  # Validate specified verbose
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("Leave-one-layer-out layer importance aborted: verbose must be TRUE or FALSE")
  }
  
  # Create function to return mean or NA
  mean.or.NA.SOM <- function(numeric_vector) {
    numeric_vector <- numeric_vector[is.finite(numeric_vector) & !is.na(numeric_vector)]
    if (length(numeric_vector) == 0) return(NA_real_)
    return(mean(numeric_vector))
  }
  
  # Create function to return median or NA
  median.or.NA.SOM <- function(numeric_vector) {
    numeric_vector <- numeric_vector[is.finite(numeric_vector) & !is.na(numeric_vector)]
    if (length(numeric_vector) == 0) return(NA_real_)
    return(stats::median(numeric_vector))
  }
  
  # Create function to build one-hot assignment probabilities from hard labels
  create.one.hot.assignment.probabilities.SOM <- function(hard_cluster_labels) {
    if (is.null(hard_cluster_labels) || length(hard_cluster_labels) == 0) {
      return(NULL)
    }
    if (is.null(names(hard_cluster_labels))) {
      return(NULL)
    }
    hard_cluster_labels <- as.character(hard_cluster_labels)
    unique_cluster_labels <- sort(unique(hard_cluster_labels))
    one_hot_assignment_probability_matrix <- matrix(0,
                                                    nrow = length(hard_cluster_labels),
                                                    ncol = length(unique_cluster_labels),
                                                    dimnames = list(names(hard_cluster_labels),
                                                                    unique_cluster_labels))
    one_hot_assignment_probability_matrix[cbind(seq_along(hard_cluster_labels),
                                                match(hard_cluster_labels, unique_cluster_labels))] <- 1
    return(one_hot_assignment_probability_matrix)
  }
  
  # Create function to normalize assignment probability matrix
  normalize.assignment.probabilities.SOM <- function(assignment_probability_matrix) {
    if (is.null(assignment_probability_matrix)) {
      return(NULL)
    }
    assignment_probability_matrix <- as.matrix(assignment_probability_matrix) #convert to matrix
    storage.mode(assignment_probability_matrix) <- "numeric"
    if (nrow(assignment_probability_matrix) == 0 || ncol(assignment_probability_matrix) == 0) { #validate matrix dimensions
      return(NULL)
    }
    
    # Replace invalid values
    assignment_probability_matrix[!is.finite(assignment_probability_matrix)] <- NA_real_
    assignment_probability_matrix[assignment_probability_matrix < 0] <- 0
    
    # Normalize rows
    assignment_row_sums <- rowSums(assignment_probability_matrix, na.rm = TRUE)
    valid_assignment_rows <- is.finite(assignment_row_sums) & assignment_row_sums > 0
    if (!any(valid_assignment_rows)) {
      return(NULL)
    }
    assignment_probability_matrix[valid_assignment_rows, ] <- assignment_probability_matrix[valid_assignment_rows, , drop = FALSE] / assignment_row_sums[valid_assignment_rows]
    
    # Keep only valid rows
    assignment_probability_matrix <- assignment_probability_matrix[valid_assignment_rows, , drop = FALSE]
    
    # Return assignment probability matrix
    return(assignment_probability_matrix)
  }
  
  # Create function to extract hard cluster labels from cluster_assignment
  get.hard.cluster.labels.from.assignment.SOM <- function(cluster_assignment_matrix) {
    
    # Validate specified cluster_assignment_matrix
    if (is.null(cluster_assignment_matrix)) {
      stop("Hard cluster label extraction aborted: cluster_assignment_matrix is NULL")
    }
    cluster_assignment_matrix <- as.matrix(cluster_assignment_matrix)
    if (nrow(cluster_assignment_matrix) == 0 || ncol(cluster_assignment_matrix) == 0) {
      stop("Hard cluster label extraction aborted: cluster_assignment_matrix is empty")
    }
    
    # Extract hard cluster labels
    hard_cluster_labels <- cluster_assignment_matrix[, 1]
    names(hard_cluster_labels) <- rownames(cluster_assignment_matrix)
    
    return(hard_cluster_labels)
  }
  
  # Create function to extract replicate-level assignment probabilities if available
  extract.single.replicate.assignment.probabilities.SOM <- function(clustered_SOM_output,
                                                                    hard_cluster_labels = NULL,
                                                                    retained_replicate_position = NULL
  ) {
    
    # Extract sample names from hard cluster labels
    sample_names <- names(hard_cluster_labels)
    if (is.null(sample_names) || length(sample_names) == 0) {
      return(create.one.hot.assignment.probabilities.SOM(hard_cluster_labels))
    }
    
    # Create function to validate and align a candidate assignment-probability matrix
    validate.and.align.assignment.probabilities.SOM <- function(candidate_assignment_probability_matrix) {
      
      # Normalize candidate matrix
      candidate_assignment_probability_matrix <- normalize.assignment.probabilities.SOM(candidate_assignment_probability_matrix)
      if (is.null(candidate_assignment_probability_matrix)) {
        return(NULL)
      }
      
      # If row names are missing but dimensions match, assume sample order matches
      if (is.null(rownames(candidate_assignment_probability_matrix)) && nrow(candidate_assignment_probability_matrix) == length(sample_names)) {
        rownames(candidate_assignment_probability_matrix) <- sample_names
      }
      
      # Require row names after attempted repair
      if (is.null(rownames(candidate_assignment_probability_matrix))) {
        return(NULL)
      }
      
      # Require overlap with requested samples
      shared_sample_names <- intersect(sample_names, rownames(candidate_assignment_probability_matrix))
      if (length(shared_sample_names) < 2) {
        return(NULL)
      }
      
      # Align to hard-label sample order
      aligned_sample_names <- sample_names[sample_names %in% rownames(candidate_assignment_probability_matrix)]
      candidate_assignment_probability_matrix <- candidate_assignment_probability_matrix[aligned_sample_names, , drop = FALSE]
      
      # Renormalize after alignment
      candidate_assignment_probability_matrix <- normalize.assignment.probabilities.SOM(candidate_assignment_probability_matrix)
      if (is.null(candidate_assignment_probability_matrix)) {
        return(NULL)
      }
      
      return(candidate_assignment_probability_matrix)
    }
    
    # Return one-hot fallback if output is invalid
    if (is.null(clustered_SOM_output) || !is.list(clustered_SOM_output)) {
      return(create.one.hot.assignment.probabilities.SOM(hard_cluster_labels))
    }
    
    # Check common list-style fields first
    possible_list_field_names <- c("replicate_ancestry_matrices",
                                   "ancestry_matrices",
                                   "sample_ancestry_matrices",
                                   "ancestry_matrix_list",
                                   "replicate.assignment.probabilities",
                                   "sample.assignment.probabilities",
                                   "replicate.ancestry.matrices",
                                   "sample.ancestry.matrices",
                                   "replicate_assignment_probability_matrices",
                                   "replicate_assignment_probabilities")
    for (field_name in possible_list_field_names) {
      if (!is.null(clustered_SOM_output[[field_name]]) && is.list(clustered_SOM_output[[field_name]])) {
        current_assignment_probability_matrix_list <- clustered_SOM_output[[field_name]]
        if (!is.null(retained_replicate_position) && length(current_assignment_probability_matrix_list) >= retained_replicate_position) {
          current_assignment_probability_matrix <- validate.and.align.assignment.probabilities.SOM(current_assignment_probability_matrix_list[[retained_replicate_position]])
          if (!is.null(current_assignment_probability_matrix)) {
            return(current_assignment_probability_matrix)
          }
        }
        if (is.null(retained_replicate_position) && length(current_assignment_probability_matrix_list) == 1) {
          current_assignment_probability_matrix <- validate.and.align.assignment.probabilities.SOM(current_assignment_probability_matrix_list[[1]])
          if (!is.null(current_assignment_probability_matrix)) {
            return(current_assignment_probability_matrix)
          }
        }
      }
    }
    
    # For single-replicate leave-one-layer-out output, ancestry_matrix may already correspond to that single run
    if (is.null(retained_replicate_position) && !is.null(clustered_SOM_output$ancestry_matrix)) {
      current_assignment_probability_matrix <- validate.and.align.assignment.probabilities.SOM(clustered_SOM_output$ancestry_matrix)
      if (!is.null(current_assignment_probability_matrix)) {
        return(current_assignment_probability_matrix)
      }
    }
    
    # Fallback to one-hot hard assignments
    return(create.one.hot.assignment.probabilities.SOM(hard_cluster_labels))
  }
  
  # Create function to summarize k distributions
  summarize.k.distribution.SOM <- function(k_values) {
    
    # Clean k values
    k_values <- as.numeric(k_values)
    k_values <- k_values[is.finite(k_values) & !is.na(k_values)]
    if (length(k_values) == 0) return(NA_character_)
    
    # Summarize k values
    k_value_table <- table(k_values)
    k_value_proportions <- prop.table(k_value_table)
    k_distribution_summary <- paste0("k", names(k_value_proportions), "=", round(as.numeric(k_value_proportions), 3), collapse = "; ")
    
    return(k_distribution_summary)
  }
  
  # Create function to calculate total variation distance between baseline and leave-one-layer-out k distributions
  calculate.k.distribution.TVD.SOM <- function(baseline_k_values,
                                               leave_one_layer_out_k_values
  ) {
    
    # Clean k values
    baseline_k_values <- as.numeric(baseline_k_values)
    leave_one_layer_out_k_values <- as.numeric(leave_one_layer_out_k_values)
    baseline_k_values <- baseline_k_values[is.finite(baseline_k_values) & !is.na(baseline_k_values)]
    leave_one_layer_out_k_values <- leave_one_layer_out_k_values[is.finite(leave_one_layer_out_k_values) & !is.na(leave_one_layer_out_k_values)]
    if (length(baseline_k_values) == 0 || length(leave_one_layer_out_k_values) == 0) {
      return(NA_real_)
    }
    
    # Extract all observed k values
    all_observed_k_values <- sort(unique(c(baseline_k_values, leave_one_layer_out_k_values)))
    
    # Calculate proportions
    baseline_k_proportions <- rep(0, length(all_observed_k_values))
    names(baseline_k_proportions) <- all_observed_k_values
    leave_one_layer_out_k_proportions <- rep(0, length(all_observed_k_values))
    names(leave_one_layer_out_k_proportions) <- all_observed_k_values
    
    baseline_k_table <- prop.table(table(baseline_k_values))
    leave_one_layer_out_k_table <- prop.table(table(leave_one_layer_out_k_values))
    
    baseline_k_proportions[names(baseline_k_table)] <- as.numeric(baseline_k_table)
    leave_one_layer_out_k_proportions[names(leave_one_layer_out_k_table)] <- as.numeric(leave_one_layer_out_k_table)
    
    # Calculate total variation distance
    k_distribution_total_variation_distance <- 0.5 * sum(abs(baseline_k_proportions - leave_one_layer_out_k_proportions))
    
    return(k_distribution_total_variation_distance)
  }
  
  # Create function to calculate pairwise co-assignment change
  calculate.pairwise.coassignment.change.SOM <- function(baseline_cluster_labels,
                                                         leave_one_layer_out_cluster_labels
  ) {
    
    # Match samples
    shared_sample_names <- intersect(names(baseline_cluster_labels), names(leave_one_layer_out_cluster_labels))
    if (length(shared_sample_names) < 2) return(NA_real_)
    
    baseline_cluster_labels <- baseline_cluster_labels[shared_sample_names]
    leave_one_layer_out_cluster_labels <- leave_one_layer_out_cluster_labels[shared_sample_names]
    
    # Calculate pairwise co-assignment matrices
    baseline_pairwise_coassignment_matrix <- outer(baseline_cluster_labels, baseline_cluster_labels, FUN = "==")
    leave_one_layer_out_pairwise_coassignment_matrix <- outer(leave_one_layer_out_cluster_labels, leave_one_layer_out_cluster_labels, FUN = "==")
    lower_triangle_indices <- lower.tri(baseline_pairwise_coassignment_matrix)
    
    # Calculate pairwise co-assignment change
    pairwise_coassignment_change <- mean(abs(as.numeric(baseline_pairwise_coassignment_matrix[lower_triangle_indices]) - as.numeric(leave_one_layer_out_pairwise_coassignment_matrix[lower_triangle_indices])), na.rm = TRUE)
    
    return(pairwise_coassignment_change)
  }
  
  # Create function to generate permutations
  generate.permutations.SOM <- function(cluster_label_vector) {
    
    # Return permutations
    if (length(cluster_label_vector) == 1) return(list(cluster_label_vector))
    permutation_list <- list()
    permutation_index <- 1
    for (cluster_index in seq_along(cluster_label_vector)) {
      remaining_cluster_labels <- cluster_label_vector[-cluster_index]
      remaining_permutations <- generate.permutations.SOM(remaining_cluster_labels)
      for (remaining_permutation_index in seq_along(remaining_permutations)) {
        permutation_list[[permutation_index]] <- c(cluster_label_vector[cluster_index], remaining_permutations[[remaining_permutation_index]])
        permutation_index <- permutation_index + 1
      }
    }
    return(permutation_list)
  }
  
  # Create function to calculate assignment accuracy compared to original baseline replicate
  calculate.assignment.accuracy.to.original.SOM <- function(baseline_cluster_labels,
                                                            leave_one_layer_out_cluster_labels
  ) {
    
    # Match samples
    shared_sample_names <- intersect(names(baseline_cluster_labels), names(leave_one_layer_out_cluster_labels))
    if (length(shared_sample_names) < 1) return(NA_real_)
    
    baseline_cluster_labels <- baseline_cluster_labels[shared_sample_names]
    leave_one_layer_out_cluster_labels <- leave_one_layer_out_cluster_labels[shared_sample_names]
    
    # If k differs, return 0
    baseline_unique_cluster_labels <- sort(unique(as.character(baseline_cluster_labels)))
    leave_one_layer_out_unique_cluster_labels <- sort(unique(as.character(leave_one_layer_out_cluster_labels)))
    if (length(baseline_unique_cluster_labels) != length(leave_one_layer_out_unique_cluster_labels)) {
      return(0)
    }
    
    # Use exact permutation matching for small k
    if (length(baseline_unique_cluster_labels) <= 8) {
      leave_one_layer_out_cluster_label_permutations <- generate.permutations.SOM(leave_one_layer_out_unique_cluster_labels)
      best_assignment_accuracy <- 0
      
      for (permutation_index in seq_along(leave_one_layer_out_cluster_label_permutations)) {
        current_cluster_label_permutation <- leave_one_layer_out_cluster_label_permutations[[permutation_index]]
        relabel_map <- setNames(baseline_unique_cluster_labels, current_cluster_label_permutation)
        relabeled_leave_one_layer_out_cluster_labels <- relabel_map[as.character(leave_one_layer_out_cluster_labels)]
        current_assignment_accuracy <- mean(as.character(baseline_cluster_labels) == relabeled_leave_one_layer_out_cluster_labels, na.rm = TRUE)
        if (current_assignment_accuracy > best_assignment_accuracy) {
          best_assignment_accuracy <- current_assignment_accuracy
        }
      }
      
      return(best_assignment_accuracy)
    }
    
    # Fallback greedy matching for larger k
    cluster_label_contingency_table <- table(as.character(leave_one_layer_out_cluster_labels), as.character(baseline_cluster_labels))
    relabel_map <- apply(cluster_label_contingency_table, 1, function(cluster_count_vector) colnames(cluster_label_contingency_table)[which.max(cluster_count_vector)])
    relabeled_leave_one_layer_out_cluster_labels <- relabel_map[as.character(leave_one_layer_out_cluster_labels)]
    assignment_accuracy_to_original <- mean(as.character(baseline_cluster_labels) == relabeled_leave_one_layer_out_cluster_labels, na.rm = TRUE)
    
    return(assignment_accuracy_to_original)
  }
  
  # Create function to calculate mean assignment margin
  calculate.mean.assignment.margin.SOM <- function(assignment_probability_matrix) {
    
    # Return missing if matrix is unavailable
    if (is.null(assignment_probability_matrix)) return(NA_real_)
    
    # k = 1 provides no meaningful separation metric
    if (ncol(assignment_probability_matrix) <= 1) {
      return(NA_real_)
    }
    
    # Calculate mean assignment margin
    sorted_row_assignment_probabilities <- t(apply(assignment_probability_matrix, 1, sort, decreasing = TRUE))
    mean_assignment_margin <- mean(sorted_row_assignment_probabilities[, 1] - sorted_row_assignment_probabilities[, 2], na.rm = TRUE)
    
    return(mean_assignment_margin)
  }
  
  # Create function to calculate mean normalized assignment entropy
  calculate.mean.normalized.assignment.entropy.SOM <- function(assignment_probability_matrix) {
    
    # Return missing if matrix is unavailable
    if (is.null(assignment_probability_matrix)) return(NA_real_)
    
    # k = 1 provides no meaningful separation metric
    if (ncol(assignment_probability_matrix) <= 1) {
      return(NA_real_)
    }
    
    # Calculate mean normalized assignment entropy
    safe_assignment_probability_matrix <- pmax(assignment_probability_matrix, .Machine$double.eps)
    row_assignment_entropies <- -rowSums(safe_assignment_probability_matrix * log(safe_assignment_probability_matrix), na.rm = TRUE)
    normalized_row_assignment_entropies <- row_assignment_entropies / log(ncol(assignment_probability_matrix))
    mean_normalized_assignment_entropy <- mean(normalized_row_assignment_entropies, na.rm = TRUE)
    
    return(mean_normalized_assignment_entropy)
  }
  
  # Set default saving name if needed
  if (is.null(save.leave.one.layer.out.results.name)) {
    input_layer_names <- names(input_data)
    save.leave.one.layer.out.results.name <- paste0("leave_one_layer_out_results_", paste(input_layer_names, collapse = "_"), ".Rdata")
  }
  
  # Set default plot file name if needed
  if (is.null(file.name)) {
    input_layer_names <- names(input_data)
    file.name <- paste0("leave_one_layer_out_layer_importance_", paste(input_layer_names, collapse = "_"), ".", plot.type)
  }
  
  # Load existing results if requested
  if (save.leave.one.layer.out.results && !overwrite.leave.one.layer.out.results && file.exists(save.leave.one.layer.out.results.name)) {
    messager("Leave-one-layer-out results already exist - loading results from file and skipping re-run")
    load(save.leave.one.layer.out.results.name)
    if (!exists("leave.one.layer.out.results")) {
      stop("Leave-one-layer-out layer importance aborted: loaded file does not contain object 'leave.one.layer.out.results'")
    }
    if (!is.list(leave.one.layer.out.results) || is.null(leave.one.layer.out.results$layer.summary) || is.null(leave.one.layer.out.results$replicate.matched.results)) {
      stop("Leave-one-layer-out layer importance aborted: loaded leave.one.layer.out.results object is malformed")
    }
    layer_summary_table <- leave.one.layer.out.results$layer.summary
    replicate_matched_results_table <- leave.one.layer.out.results$replicate.matched.results
  } else {
    
    # Create function to fit one matched single-replicate leave-one-layer-out SOM
    fit.single.replicate.leave.one.layer.out.SOM <- function(input_data_for_SOM,
                                                             train.SOM.args = list(),
                                                             clustering.SOM.args = list(),
                                                             matched_training_seed,
                                                             matched_clustering_seed
    ) {
      
      # Override training to one replicate and non-parallel for deterministic seed matching
      train.SOM.args$N.replicates <- 1
      train.SOM.args$parallel <- FALSE
      train.SOM.args$set.seed.N <- matched_training_seed
      train.SOM.args$save.SOM.results <- FALSE
      train.SOM.args$save.SOM.results.name <- NULL
      train.SOM.args$overwrite.SOM.results <- TRUE
      train.SOM.args$verbose <- FALSE
      
      # Fit leave-one-layer-out SOM
      trained_single_replicate_SOM_output <- do.call(train.SOM,
                                                     c(list(input_data = input_data_for_SOM), train.SOM.args))
      
      # Extract clustering arguments with defaults
      max.k.current <- if (!is.null(clustering.SOM.args$max.k)) clustering.SOM.args$max.k else 10
      set.k.current <- if (!is.null(clustering.SOM.args$set.k)) clustering.SOM.args$set.k else NULL
      clustering.method.current <- if (!is.null(clustering.SOM.args$clustering.method)) clustering.SOM.args$clustering.method else stop("Leave-one-layer-out layer importance aborted: clustering.method is missing in clustering.SOM.args")
      BIC.thresh.current <- if (!is.null(clustering.SOM.args$BIC.thresh)) clustering.SOM.args$BIC.thresh else 6
      
      # Extract codebook matrix and count distinct codebook rows
      trained_codes <- kohonen::getCodes(trained_single_replicate_SOM_output$som_models[[1]])
      if (!is.list(trained_codes)) trained_codes <- list(trained_codes)
      trained_code_matrix <- do.call(cbind, lapply(trained_codes, as.matrix))
      distinct_codebook_rows <- nrow(unique(trained_code_matrix))
      
      # If only one distinct codebook row remains, return a valid k = 1 clustering result directly
      if (distinct_codebook_rows < 2) {
        single_cluster_assignment <- matrix(1,
                                            nrow = nrow(input_data_for_SOM[[1]]),
                                            ncol = 1,
                                            dimnames = list(rownames(input_data_for_SOM[[1]]),
                                                            "R1"))
        clustered_single_replicate_SOM_output <- trained_single_replicate_SOM_output
        clustered_single_replicate_SOM_output$cluster_assignment <- single_cluster_assignment
        clustered_single_replicate_SOM_output$optim_k_vals <- 1
        clustered_single_replicate_SOM_output$optim_k_summary <- data.frame(Count = 1,
                                                                            Proportion = 1,
                                                                            row.names = "k1")
      } else {
        
        # Reduce max.k / set.k if leave-one-layer-out SOM has too few distinct codebook rows
        max.k.allowed <- distinct_codebook_rows - 1
        max.k.current <- min(max.k.current, max.k.allowed)
        if (!is.null(set.k.current)) {
          set.k.current <- min(set.k.current, max.k.allowed)
        }
        
        # Fit clustering
        clustered_single_replicate_SOM_output <- tryCatch(
          withCallingHandlers(
            clustering.SOM(SOM.output = trained_single_replicate_SOM_output,
                           max.k = max.k.current,
                           set.k = set.k.current,
                           clustering.method = clustering.method.current,
                           BIC.thresh = BIC.thresh.current,
                           quantization.error.quantile = NULL,
                           topographic.error.quantile = NULL,
                           set.seed.N = matched_clustering_seed),
            warning = function(w) {
              if (grepl("^Eta squared effect size \\(variable importance\\) could not be computed because all replicates produced k = 1$",
                        conditionMessage(w))) {
                invokeRestart("muffleWarning")
              }
            }
          ),
          error = function(error_message) {
            if (grepl("mehr Clusterzentren als verschiedene Datenpunkte|more cluster centers than distinct data points",
                      conditionMessage(error_message))) {
              single_cluster_assignment <- matrix(1,
                                                 nrow = nrow(input_data_for_SOM[[1]]),
                                                  ncol = 1,
                                                  dimnames = list(rownames(input_data_for_SOM[[1]]),
                                                                  "R1"))
              clustered_single_replicate_SOM_output <- trained_single_replicate_SOM_output
              clustered_single_replicate_SOM_output$cluster_assignment <- single_cluster_assignment
              clustered_single_replicate_SOM_output$optim_k_vals <- 1
              clustered_single_replicate_SOM_output$optim_k_summary <- data.frame(Count = 1,
                                                                                  Proportion = 1,
                                                                                  row.names = "k1")
              return(clustered_single_replicate_SOM_output)
            }
            stop(error_message)
          }
        )
      }
 return(clustered_single_replicate_SOM_output)
    }
    
    # Extract baseline retained replicate indices
    if (!is.null(SOM_output$retained_replicates) && length(SOM_output$retained_replicates) == length(SOM_output$som_models)) {
      baseline_retained_replicate_indices <- SOM_output$retained_replicates
      if (is.factor(baseline_retained_replicate_indices)) {
        baseline_retained_replicate_indices <- as.character(baseline_retained_replicate_indices)
      }
      if (is.list(baseline_retained_replicate_indices)) {
        baseline_retained_replicate_indices <- unlist(baseline_retained_replicate_indices, recursive = TRUE, use.names = FALSE)
      }
      if (is.character(baseline_retained_replicate_indices)) {
        baseline_retained_replicate_indices <- sub("^R", "", baseline_retained_replicate_indices)
      }
      baseline_retained_replicate_indices <- suppressWarnings(as.integer(baseline_retained_replicate_indices))
      if (any(!is.finite(baseline_retained_replicate_indices)) || any(is.na(baseline_retained_replicate_indices)) || any(baseline_retained_replicate_indices < 1)) {
        stop("Leave-one-layer-out layer importance aborted: SOM_output$retained_replicates could not be converted to positive integer replicate indices")
      }
    } else {
      baseline_retained_replicate_indices <- seq_along(SOM_output$som_models)
    }
    
    # Extract baseline replicate-wise cluster assignments
    baseline_cluster_assignment_matrix <- as.matrix(SOM_output$cluster_assignment)
    if (ncol(baseline_cluster_assignment_matrix) != length(SOM_output$som_models)) {
      stop("Leave-one-layer-out layer importance aborted: number of cluster_assignment columns does not match number of retained som_models")
    }
    
    # Extract baseline replicate-wise optimal k values
    baseline_optimal_k_values <- as.numeric(SOM_output$optim_k_vals)
    if (length(baseline_optimal_k_values) != length(SOM_output$som_models)) {
      stop("Leave-one-layer-out layer importance aborted: number of optim_k_vals does not match number of retained som_models")
    }
    
    # Run replicate-matched leave-one-layer-out analyses
    messager("RUNNING LEAVE-ONE-LAYER-OUT ANALYSES ...")
    replicate_matched_results_list <- vector("list", length(SOM_output$som_models) * length(input_data))
    results_counter <- 1
    for (retained_replicate_position in seq_along(SOM_output$som_models)) {
      
      # Show message
      if (retained_replicate_position %% message.N.replicates == 0) messager("Running replicate: ", retained_replicate_position, " of ", length(SOM_output$som_models))
      
      # Extract baseline replicate seed indices
      baseline_training_replicate_index <- baseline_retained_replicate_indices[retained_replicate_position]
      if (is.character(baseline_training_replicate_index)) {
        baseline_training_replicate_index <- sub("^R", "", baseline_training_replicate_index)
      }
      baseline_training_replicate_index <- suppressWarnings(as.integer(baseline_training_replicate_index))
      if (!is.finite(baseline_training_replicate_index) || is.na(baseline_training_replicate_index) || baseline_training_replicate_index < 1) {
        stop("Leave-one-layer-out layer importance aborted: baseline retained replicate index could not be converted to a positive integer")
      }
      matched_training_seed <- as.integer(baseline.train.SOM.set.seed.N + baseline_training_replicate_index - 1)
      matched_clustering_seed <- as.integer(baseline.clustering.SOM.set.seed.N + retained_replicate_position - 1)
      
      # Extract baseline replicate results directly from stored output
      baseline_single_replicate_cluster_assignment <- baseline_cluster_assignment_matrix[, retained_replicate_position, drop = FALSE]
      baseline_hard_cluster_labels <- get.hard.cluster.labels.from.assignment.SOM(baseline_single_replicate_cluster_assignment)
      baseline_modal_k <- baseline_optimal_k_values[retained_replicate_position]
      baseline_sample_names <- names(baseline_hard_cluster_labels)
      
      # Extract baseline assignment probabilities if available, otherwise fallback to one-hot hard assignments
      baseline_assignment_probability_matrix <- extract.single.replicate.assignment.probabilities.SOM(clustered_SOM_output = SOM_output,
                                                                                                      hard_cluster_labels = baseline_hard_cluster_labels,
                                                                                                      retained_replicate_position = retained_replicate_position)
      baseline_mean_assignment_margin <- calculate.mean.assignment.margin.SOM(baseline_assignment_probability_matrix)
      baseline_mean_normalized_assignment_entropy <- calculate.mean.normalized.assignment.entropy.SOM(baseline_assignment_probability_matrix)
      
      # Loop through omitted layers
      for (layer_index in seq_along(input_data)) {
        
        # Extract omitted layer name
        omitted_layer_name <- names(input_data)[layer_index]
        
        # Create leave-one-layer-out input
        leave_one_layer_out_input_data <- input_data[-layer_index]
        
        # Create fresh per-run copies of argument lists
        current_train.SOM.args <- train.SOM.args
        current_clustering.SOM.args <- clustering.SOM.args
        
        # Subset layer.distance.functions to match remaining layers
        if (!is.null(current_train.SOM.args$layer.distance.functions)) {
          if (length(current_train.SOM.args$layer.distance.functions) > 1) {
            current_train.SOM.args$layer.distance.functions <- current_train.SOM.args$layer.distance.functions[-layer_index]
          }
        }
        
        # Subset manual.layer.weights to match remaining layers
        if (!is.null(current_train.SOM.args$manual.layer.weights)) {
          if (length(current_train.SOM.args$manual.layer.weights) > 1) {
            current_train.SOM.args$manual.layer.weights <- current_train.SOM.args$manual.layer.weights[-layer_index]
          }
        }
        
        # Fit leave-one-layer-out matched single-replicate SOM
        leave_one_layer_out.SOM.output <- tryCatch({
          fit.single.replicate.leave.one.layer.out.SOM(input_data_for_SOM = leave_one_layer_out_input_data,
                                                       train.SOM.args = current_train.SOM.args,
                                                       clustering.SOM.args = current_clustering.SOM.args,
                                                       matched_training_seed = matched_training_seed,
                                                       matched_clustering_seed = matched_clustering_seed)
        }, error = function(error_message) error_message)
        
        # Store failed leave-one-layer-out run
        if (inherits(leave_one_layer_out.SOM.output, "error")) {
          replicate_matched_results_list[[results_counter]] <- data.frame(retained.replicate.position = retained_replicate_position,
                                                                          retained.replicate.index = baseline_training_replicate_index,
                                                                          layer = omitted_layer_name,
                                                                          baseline.modal.k = baseline_modal_k,
                                                                          leave.one.layer.out.modal.k = NA_real_,
                                                                          signed.k.shift = NA_real_,
                                                                          absolute.k.deviation = NA_real_,
                                                                          assignment.accuracy.to.original = NA_real_,
                                                                          pairwise.coassignment.change = NA_real_,
                                                                          ARI = NA_real_,
                                                                          baseline.mean.assignment.margin = baseline_mean_assignment_margin,
                                                                          leave.one.layer.out.mean.assignment.margin = NA_real_,
                                                                          delta.mean.assignment.margin = NA_real_,
                                                                          baseline.mean.normalized.assignment.entropy = baseline_mean_normalized_assignment_entropy,
                                                                          leave.one.layer.out.mean.normalized.assignment.entropy = NA_real_,
                                                                          increase.mean.normalized.assignment.entropy = NA_real_,
                                                                          error = conditionMessage(leave_one_layer_out.SOM.output),
                                                                          stringsAsFactors = FALSE)
          results_counter <- results_counter + 1
          next
        }
        
        # Extract leave-one-layer-out replicate results
        leave_one_layer_out_cluster_assignment_matrix <- as.matrix(leave_one_layer_out.SOM.output$cluster_assignment)
        leave_one_layer_out_hard_cluster_labels <- get.hard.cluster.labels.from.assignment.SOM(leave_one_layer_out_cluster_assignment_matrix)
        leave_one_layer_out_modal_k <- as.numeric(leave_one_layer_out.SOM.output$optim_k_vals)[1]
        
        # Extract shared sample names
        shared_sample_names <- intersect(baseline_sample_names,
                                         names(leave_one_layer_out_hard_cluster_labels))
        if (length(shared_sample_names) < 2) {
          replicate_matched_results_list[[results_counter]] <- data.frame(retained.replicate.position = retained_replicate_position,
                                                                          retained.replicate.index = baseline_training_replicate_index,
                                                                          layer = omitted_layer_name,
                                                                          baseline.modal.k = baseline_modal_k,
                                                                          leave.one.layer.out.modal.k = leave_one_layer_out_modal_k,
                                                                          signed.k.shift = leave_one_layer_out_modal_k - baseline_modal_k,
                                                                          absolute.k.deviation = abs(leave_one_layer_out_modal_k - baseline_modal_k),
                                                                          assignment.accuracy.to.original = NA_real_,
                                                                          pairwise.coassignment.change = NA_real_,
                                                                          ARI = NA_real_,
                                                                          baseline.mean.assignment.margin = baseline_mean_assignment_margin,
                                                                          leave.one.layer.out.mean.assignment.margin = NA_real_,
                                                                          delta.mean.assignment.margin = NA_real_,
                                                                          baseline.mean.normalized.assignment.entropy = baseline_mean_normalized_assignment_entropy,
                                                                          leave.one.layer.out.mean.normalized.assignment.entropy = NA_real_,
                                                                          increase.mean.normalized.assignment.entropy = NA_real_,
                                                                          error = "Too few shared samples between stored baseline replicate and leave-one-layer-out output",
                                                                          stringsAsFactors = FALSE)
          results_counter <- results_counter + 1
          next
        }
        
        # Extract leave-one-layer-out assignment probabilities if available, otherwise fallback to one-hot hard assignments
        leave_one_layer_out_assignment_probability_matrix <- extract.single.replicate.assignment.probabilities.SOM(clustered_SOM_output = leave_one_layer_out.SOM.output,
                                                                                                                   hard_cluster_labels = leave_one_layer_out_hard_cluster_labels,
                                                                                                                   retained_replicate_position = NULL)
        leave_one_layer_out_mean_assignment_margin <- calculate.mean.assignment.margin.SOM(leave_one_layer_out_assignment_probability_matrix)
        leave_one_layer_out_mean_normalized_assignment_entropy <- calculate.mean.normalized.assignment.entropy.SOM(leave_one_layer_out_assignment_probability_matrix)
        
        # Calculate replicate-level metrics
        signed.k.shift <- leave_one_layer_out_modal_k - baseline_modal_k
        absolute.k.deviation <- abs(signed.k.shift)
        assignment.accuracy.to.original <- calculate.assignment.accuracy.to.original.SOM(baseline_cluster_labels = baseline_hard_cluster_labels,
                                                                                         leave_one_layer_out_cluster_labels = leave_one_layer_out_hard_cluster_labels)
        pairwise.coassignment.change <- calculate.pairwise.coassignment.change.SOM(baseline_cluster_labels = baseline_hard_cluster_labels,
                                                                                   leave_one_layer_out_cluster_labels = leave_one_layer_out_hard_cluster_labels)
        adjusted.rand.index <- tryCatch({
          mclust::adjustedRandIndex(baseline_hard_cluster_labels[shared_sample_names],
                                    leave_one_layer_out_hard_cluster_labels[shared_sample_names])
        }, error = function(error_message) NA_real_)
        
        # Calculate new continuous certainty/separation metrics
        delta.mean.assignment.margin <- baseline_mean_assignment_margin - leave_one_layer_out_mean_assignment_margin
        increase.mean.normalized.assignment.entropy <- leave_one_layer_out_mean_normalized_assignment_entropy - baseline_mean_normalized_assignment_entropy
        
        # Store results
        replicate_matched_results_list[[results_counter]] <- data.frame(retained.replicate.position = retained_replicate_position,
                                                                        retained.replicate.index = baseline_training_replicate_index,
                                                                        layer = omitted_layer_name,
                                                                        baseline.modal.k = baseline_modal_k,
                                                                        leave.one.layer.out.modal.k = leave_one_layer_out_modal_k,
                                                                        signed.k.shift = signed.k.shift,
                                                                        absolute.k.deviation = absolute.k.deviation,
                                                                        assignment.accuracy.to.original = assignment.accuracy.to.original,
                                                                        pairwise.coassignment.change = pairwise.coassignment.change,
                                                                        ARI = adjusted.rand.index,
                                                                        baseline.mean.assignment.margin = baseline_mean_assignment_margin,
                                                                        leave.one.layer.out.mean.assignment.margin = leave_one_layer_out_mean_assignment_margin,
                                                                        delta.mean.assignment.margin = delta.mean.assignment.margin,
                                                                        baseline.mean.normalized.assignment.entropy = baseline_mean_normalized_assignment_entropy,
                                                                        leave.one.layer.out.mean.normalized.assignment.entropy = leave_one_layer_out_mean_normalized_assignment_entropy,
                                                                        increase.mean.normalized.assignment.entropy = increase.mean.normalized.assignment.entropy,
                                                                        error = NA_character_,
                                                                        stringsAsFactors = FALSE)
        results_counter <- results_counter + 1
      }
    }
    
    # Combine replicate-matched results
    replicate_matched_results_table <- do.call(rbind, replicate_matched_results_list)
    
    # Create layer-level summary table
    layer_summary_list <- lapply(split(replicate_matched_results_table, replicate_matched_results_table$layer), function(current_layer_results_table) {
      
      # Extract successful runs
      successful_current_layer_results_table <- current_layer_results_table[is.na(current_layer_results_table$error), , drop = FALSE]
      
      # Calculate baseline and leave-one-layer-out k distributions
      baseline.k.distribution <- summarize.k.distribution.SOM(successful_current_layer_results_table$baseline.modal.k)
      leave.one.layer.out.k.distribution <- summarize.k.distribution.SOM(successful_current_layer_results_table$leave.one.layer.out.modal.k)
      
      # Calculate k distribution TVD
      k.distribution.TVD <- calculate.k.distribution.TVD.SOM(baseline_k_values = successful_current_layer_results_table$baseline.modal.k,
                                                             leave_one_layer_out_k_values = successful_current_layer_results_table$leave.one.layer.out.modal.k)
      
      # Return summary row
      data.frame(layer = current_layer_results_table$layer[1],
                 N.replicates = nrow(current_layer_results_table),
                 N.successful = nrow(successful_current_layer_results_table),
                 mean.signed.k.shift = mean.or.NA.SOM(successful_current_layer_results_table$signed.k.shift),
                 median.signed.k.shift = median.or.NA.SOM(successful_current_layer_results_table$signed.k.shift),
                 mean.absolute.k.deviation = mean.or.NA.SOM(successful_current_layer_results_table$absolute.k.deviation),
                 median.absolute.k.deviation = median.or.NA.SOM(successful_current_layer_results_table$absolute.k.deviation),
                 mean.assignment.accuracy.to.original = mean.or.NA.SOM(successful_current_layer_results_table$assignment.accuracy.to.original),
                 median.assignment.accuracy.to.original = median.or.NA.SOM(successful_current_layer_results_table$assignment.accuracy.to.original),
                 mean.pairwise.coassignment.change = mean.or.NA.SOM(successful_current_layer_results_table$pairwise.coassignment.change),
                 median.pairwise.coassignment.change = median.or.NA.SOM(successful_current_layer_results_table$pairwise.coassignment.change),
                 mean.ARI = mean.or.NA.SOM(successful_current_layer_results_table$ARI),
                 median.ARI = median.or.NA.SOM(successful_current_layer_results_table$ARI),
                 mean.baseline.mean.assignment.margin = mean.or.NA.SOM(successful_current_layer_results_table$baseline.mean.assignment.margin),
                 mean.leave.one.layer.out.mean.assignment.margin = mean.or.NA.SOM(successful_current_layer_results_table$leave.one.layer.out.mean.assignment.margin),
                 mean.delta.mean.assignment.margin = mean.or.NA.SOM(successful_current_layer_results_table$delta.mean.assignment.margin),
                 median.delta.mean.assignment.margin = median.or.NA.SOM(successful_current_layer_results_table$delta.mean.assignment.margin),
                 mean.baseline.mean.normalized.assignment.entropy = mean.or.NA.SOM(successful_current_layer_results_table$baseline.mean.normalized.assignment.entropy),
                 mean.leave.one.layer.out.mean.normalized.assignment.entropy = mean.or.NA.SOM(successful_current_layer_results_table$leave.one.layer.out.mean.normalized.assignment.entropy),
                 mean.increase.mean.normalized.assignment.entropy = mean.or.NA.SOM(successful_current_layer_results_table$increase.mean.normalized.assignment.entropy),
                 median.increase.mean.normalized.assignment.entropy = median.or.NA.SOM(successful_current_layer_results_table$increase.mean.normalized.assignment.entropy),
                 k.distribution.TVD = k.distribution.TVD,
                 baseline.k.distribution = baseline.k.distribution,
                 leave.one.layer.out.k.distribution = leave.one.layer.out.k.distribution,
                 stringsAsFactors = FALSE)
    })
    layer_summary_table <- do.call(rbind, layer_summary_list)
    layer_summary_table <- layer_summary_table[order(layer_summary_table$mean.absolute.k.deviation,
                                                     decreasing = TRUE), , drop = FALSE]
    rownames(layer_summary_table) <- NULL
    
    # Create results object
    leave.one.layer.out.results <- list(layer.summary = layer_summary_table,
                                        replicate.matched.results = replicate_matched_results_table)
    
    # Save results if requested
    if (save.leave.one.layer.out.results) {
      
      # Check if directory exists
      save_directory_path <- dirname(save.leave.one.layer.out.results.name)
      if (!dir.exists(save_directory_path)) {
        dir.create(save_directory_path, recursive = TRUE)
        messager(paste("Specified directory", save_directory_path, "did not exist and was created"))
      }
      
      # Save results
      save(leave.one.layer.out.results, file = save.leave.one.layer.out.results.name)
      if (!overwrite.leave.one.layer.out.results) {
        messager("Leave-one-layer-out results saved as ", save.leave.one.layer.out.results.name)
      }
      if (overwrite.leave.one.layer.out.results) {
        messager("Leave-one-layer-out results overwritten as ", save.leave.one.layer.out.results.name)
      }
    }
  }
  
  # Filter to successful replicate-level results
  successful_replicate_matched_results_table <- replicate_matched_results_table[is.na(replicate_matched_results_table$error) | replicate_matched_results_table$error == "", , drop = FALSE]
  if (nrow(successful_replicate_matched_results_table) == 0) {
    warning("Plotting skipped: no successful replicate-level results found")
    return(leave.one.layer.out.results)
  }
  
  # Calculate assignment inaccuracy to original
  successful_replicate_matched_results_table$assignment.inaccuracy.to.original <- 1 - successful_replicate_matched_results_table$assignment.accuracy.to.original
  
  # Extract layer names and colors
  SOM_layer_names <- as.character(layer_summary_table$layer)
  layer_colors <- setNames(col.pal(length(SOM_layer_names)), SOM_layer_names)
  
  # Reorder replicate-level table to match layer.summary order
  successful_replicate_matched_results_table$layer <- factor(successful_replicate_matched_results_table$layer,
                                                             levels = SOM_layer_names)
  
  # Open plot device if requested
  if (save) {
    
    # Check whether plot file already exists
    if (file.exists(file.name) && !overwrite) {
      stop("Leave-one-layer-out layer importance aborted: plot file already exists and overwrite = FALSE")
    }
    
    # Check if directory exists
    plot_directory_path <- dirname(file.name)
    if (!dir.exists(plot_directory_path)) {
      dir.create(plot_directory_path, recursive = TRUE)
      messager(paste("Specified directory", plot_directory_path, "did not exist and was created"))
    }
    
    # Open graphics device
    if (plot.type == "svg") {
      grDevices::svg(filename = file.name,
                     width = width / 2.54,
                     height = height / 2.54)
    }
    if (plot.type == "png") {
      grDevices::png(filename = file.name,
                     width = width,
                     height = height,
                     units = "cm",
                     res = resolution)
    }
    if (plot.type == "jpg") {
      grDevices::jpeg(filename = file.name,
                      width = width,
                      height = height,
                      units = "cm",
                      res = resolution)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  # Reset plotting parameters
  old_plotting_parameters <- par(no.readonly = TRUE)
  on.exit(par(old_plotting_parameters), add = TRUE)
  
  # Determine whether certainty/separation plots can be shown
  show.assignment.margin.plot <- any(is.finite(successful_replicate_matched_results_table$delta.mean.assignment.margin) &
                                       !is.na(successful_replicate_matched_results_table$delta.mean.assignment.margin))
  show.assignment.entropy.plot <- any(is.finite(successful_replicate_matched_results_table$increase.mean.normalized.assignment.entropy) &
                                        !is.na(successful_replicate_matched_results_table$increase.mean.normalized.assignment.entropy))
  
  # Set plotting layout
  if (show.assignment.margin.plot || show.assignment.entropy.plot) {
    par(mfrow = c(3, 2),
        mar = c(bottom.margin,
                left.margin,
                top.margin,
                right.margin),
        oma = if (is.null(title)) c(0, 0, 0, 0) else c(0, 0, 2, 0),
        mgp = c(distance.axis.label, 1, 0))
  } else {
    par(mfrow = c(2, 2),
        mar = c(bottom.margin,
                left.margin,
                top.margin,
                right.margin),
        oma = if (is.null(title)) c(0, 0, 0, 0) else c(0, 0, 2, 0),
        mgp = c(distance.axis.label, 1, 0))
    messager("Assignment margin and entropy plots skipped because all successful replicates had k = 1")
  }
  
  # Add overall title if requested
  if (!is.null(title)) {
    graphics::mtext(title, side = 3, outer = TRUE, line = -1.5, cex = 1.2)
  }
  
  # Create function to add jittered points
  add.jittered.points.SOM <- function(response_variable_name) {
    
    # Add points if requested
    if (isTRUE(add.points)) {
      for (layer_index in seq_along(SOM_layer_names)) {
        current_layer_name <- SOM_layer_names[layer_index]
        current_layer_values <- successful_replicate_matched_results_table[successful_replicate_matched_results_table$layer == current_layer_name, response_variable_name]
        current_layer_values <- current_layer_values[is.finite(current_layer_values) & !is.na(current_layer_values)]
        if (length(current_layer_values) == 0) next
        points(jitter(rep(layer_index, length(current_layer_values)), amount = 0.15),
               current_layer_values,
               pch = 16,
               cex = point.cex,
               col = grDevices::adjustcolor(layer_colors[current_layer_name], alpha.f = point.alpha))
      }
    }
  }
  
  # Plot absolute k deviation
  boxplot(absolute.k.deviation ~ layer,
          data = successful_replicate_matched_results_table,
          col = layer_colors[SOM_layer_names],
          outline = FALSE,
          xaxt = "n",
          las = 2,
          ylab = "Absolute k deviation",
          xlab = "",
          main = "Absolute k deviation")
  add.jittered.points.SOM("absolute.k.deviation")
  
  # Plot k distribution TVD
  barplot(height = layer_summary_table$k.distribution.TVD,
          axisnames = FALSE,
          col = layer_colors[SOM_layer_names],
          las = 2,
          ylab = "K distribution TVD",
          xlab = "",
          main = "K distribution TVD")
  
  # Plot assignment inaccuracy to original
  boxplot(assignment.inaccuracy.to.original ~ layer,
          data = successful_replicate_matched_results_table,
          col = layer_colors[SOM_layer_names],
          outline = FALSE,
          las = 2,
          ylab = "1 - assignment accuracy to original",
          xlab = "",
          main = "Assignment inaccuracy")
  add.jittered.points.SOM("assignment.inaccuracy.to.original")
  
  # Plot pairwise co-assignment change
  boxplot(pairwise.coassignment.change ~ layer,
          data = successful_replicate_matched_results_table,
          col = layer_colors[SOM_layer_names],
          outline = FALSE,
          las = 2,
          ylab = "Pairwise co-assignment change",
          xlab = "",
          main = "Pairwise co-assignment change")
  add.jittered.points.SOM("pairwise.coassignment.change")
  
  # Plot assignment margin change
  if (show.assignment.margin.plot) {
    boxplot(delta.mean.assignment.margin ~ layer,
            data = successful_replicate_matched_results_table,
            col = layer_colors[SOM_layer_names],
            outline = FALSE,
            las = 2,
            ylab = "Baseline - leave-one-layer-out margin",
            xlab = "",
            main = "Assignment margin change")
    add.jittered.points.SOM("delta.mean.assignment.margin")
    abline(h = 0, lty = 2)
  }
  
  # Plot normalized assignment entropy increase
  if (show.assignment.entropy.plot) {
    boxplot(increase.mean.normalized.assignment.entropy ~ layer,
            data = successful_replicate_matched_results_table,
            col = layer_colors[SOM_layer_names],
            outline = FALSE,
            las = 2,
            ylab = "Leave-one-layer-out - baseline entropy",
            xlab = "",
            main = "Entropy increase")
    add.jittered.points.SOM("increase.mean.normalized.assignment.entropy")
    abline(h = 0, lty = 2)
  }
  
  # Report saved plot if requested
  if (save) {
    if (overwrite) {
      messager("Leave-one-layer-out plot overwritten as ", file.name)
    }
    if (!overwrite) {
      messager("Leave-one-layer-out plot saved as ", file.name)
    }
  }
  
  # Return results
  return(leave.one.layer.out.results)
}
