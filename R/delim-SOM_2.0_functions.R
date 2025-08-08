###################################################################################
#### SOM-based species delimitation (Pyron et al. 2022, 2023) v.2.0
###################################################################################

## Install and load required R packages
# Define CRAN and Bioconductor packages
CRAN_packages <- c(
  "adegenet",    # genetic data manipulation
  "ape",         # phylogenetics
  "caroline",    # pie charts
  "dbscan",      # HDBSCAN clustering
  "doRNG",       # reproducible RNG
  "doParallel",  # parallel backend
  "dplyr",       # data manipulation
  "FactoMineR",  # multiple factor analysis
  "FNN",         # nearest neighbor methods
  "foreach",     # parallel processing
  "kernlab",     # spectral clustering
  "kohonen",     # SOM / supersom
  "maps",        # mapping
  "mclust",      # Gaussian mixture models
  "poppr",       # population genetics
  "readr",       # reading data
  "sf",          # spatial data
  "stringr",     # string manipulation
  "vcfR",        # VCF file handling
  "viridis"      # color palettes
)

Bioconductor_packages <- c("SeqArray")

# Install missing CRAN packages
for (pkg in CRAN_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install missing Bioconductor packages
for (bio_pkg in Bioconductor_packages) {
  if (!requireNamespace(bio_pkg, quietly = TRUE)) BiocManager::install(bio_pkg)
}

# Load all packages
required_packages <- c(CRAN_packages, Bioconductor_packages)
lapply(required_packages, library, character.only = TRUE)


#### Functions

## Function to train single-layer SOM (one matrix) or multi-layer Super-SOM (multiple matrices)
train.SOM <- function(input_data, #one matrix/dataframe or multiple matrices/dataframes provided as list()
                      N.steps = 100, #number of training iterations S for SOM
                      N.replicates = 30, #number of SOM runs R
                      parallel = FALSE, #run SOM training in parallel 
                      N_cores = 4, #number of cores for training SOM in parallel
                      grid.size = NULL, #grid size - specify as c(x, y) if desired
                      grid.multiplier = 5, #tuning parameter defining how “fine” or “coarse” SOM grid is relative to number of samples (recommended: 5)
                      learning.rate.initial = 0.6, #initial learning rate for SOM training
                      learning.rate.final = 0.2, #final learning rate for SOM training
                      learning.rate.tuning = TRUE, #set learning.rate tuning to choose best initial and final learning rate values
                      max.NA.row = 0.4, #maximum fraction of missing values allowed per row (sample) in input data to prevent row to be removed
                      max.NA.col = 0.7, #maximum fraction of missing values allowed per column (variable) in input data to prevent row to be removed
                      training.neighborhoods = "gaussian", #neighborhood function used for SOM training (options: "gaussian" or "bubble")
                      save.SOM.results = FALSE, #whether to save SOM results to file
                      save.SOM.results.name = NULL, #file name for saving SOM results (if NULL, default name based on input_data is generated)
                      overwrite.SOM.results = TRUE, #if FALSE, existing results are loaded instead of re-running SOM
                      verbose = TRUE, #show messages
                      message.N.replicates = 10 #frequency of progress messages during training (message is printed every message.N.replicates iterations)
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
    stop("Data processing aborted: input_data must be matrix, data.frame or list of such objects")
  }
  if (is.list(input_data) && length(input_data) == 0) {
    stop("Data processing aborted: input_data is empty list")
  }
  
  # Validate specified N.steps
  if (!is.numeric(N.steps) || length(N.steps) != 1 || is.na(N.steps) || N.steps < 1 || (N.steps %% 1 != 0)) {
    stop("Data processing aborted: N.steps must be single positive integer (>= 1)")
  }
  if (N.steps < 20) {
    messager("Warning: N.steps is low (", N.steps, ") - SOM training may be unstable (recommended: 50–200)")
  }
  if (N.steps > 500) {
    messager("Warning: N.steps is high (", N.steps, ") - computation will be slow (recommended: 50–200)")
  }
  
  # Validate specified N.replicates
  if (!is.numeric(N.replicates) || length(N.replicates) != 1 || is.na(N.replicates) || N.replicates < 1 || (N.replicates %% 1 != 0)) {
    stop("Data processing aborted: N.replicates must be single positive integer (>= 1)")
  }
  if (N.replicates < 10) {
    messager("Warning: N.replicates is low (", N.replicates, ") - results may be unreliable (recommended: 30–100)")
  }
  if (N.replicates > 150) {
    messager("Warning: N.replicates is high (", N.replicates, ") - computation will be slow (recommended: 30–100)")
  }
  
  # Validate specified parallel
  if (!is.logical(parallel) || length(parallel) != 1 || is.na(parallel)) {
    stop("Data processing aborted: parallel must be single logical value (TRUE or FALSE)")
  }
  
  # Validate specified N_cores
  if (parallel) {
    if (!is.numeric(N_cores) || length(N_cores) != 1 || is.na(N_cores) || N_cores < 1 || (N_cores %% 1 != 0)) {
      stop("Data processing aborted: N_cores must be single positive integer (>= 1)")
    }
    max_cores <- parallel::detectCores(logical = FALSE) #detect number of physical cores available
    if (N_cores > max_cores) {
      stop(sprintf("Data processing aborted: requested N_cores (%d) exceeds available cores (%d) - set N_cores <= %d", N_cores, max_cores, max_cores))
    }
  }
  
  # Validate specified grid.size
  if (!is.null(grid.size)) {
    if (!is.numeric(grid.size) || length(grid.size) != 2 || any(is.na(grid.size)) ||
        any(grid.size <= 0) || any(grid.size %% 1 != 0)) {
      stop("Input aborted: grid.size must be NULL or numeric vector of length 2 with positive integers (e.g., c(5, 5))")
    }
  }
  
  # Validate specified grid.multiplier
  if (!is.numeric(grid.multiplier) || length(grid.multiplier) != 1 || grid.multiplier < 1) {
    stop("Data processing aborted: grid.multiplier must be single numeric value (recommended: 5)")
  }
  
  # Validate specified learning rate parameters if learning rate tuning is not done
  if (!learning.rate.tuning) {
    if (!is.numeric(learning.rate.initial) || 
        length(learning.rate.initial) != 1 ||
        learning.rate.initial <= 0 || 
        learning.rate.initial > 1) { #initial learning rate
      stop("Data processing aborted: learning.rate.initial must be single numeric value between 0 and 1 (e.g., 0.1-0.6)")
    }
    if (!is.numeric(learning.rate.final) || 
        length(learning.rate.final) != 1 ||
        learning.rate.final < 0) 
    {
      stop("Data processing aborted: learning.rate.final must be single numeric value between 0 and 1 (e.g., 0.001–0.1)")
    }
    if (learning.rate.final > learning.rate.initial) {
      stop("Data processing aborted: learning.rate.final must be smaller than learning.rate.initial")
    }
  }
  
  # Validate specified learning.rate.tuning
  if (!is.logical(learning.rate.tuning) || length(learning.rate.tuning) != 1 || is.na(learning.rate.tuning)) {
    stop("Data processing aborted: learning.rate.tuning must be single logical value (TRUE or FALSE)")
  }
  
  # Validate specified NA‐max.NA.row
  if (!is.numeric(max.NA.row) || length(max.NA.row) != 1 ||
      max.NA.row < 0 || max.NA.row > 1) {
    stop("Data processing aborted: max.NA.row must be single numeric value between 0 and 1 (e.g. 0.5)")
  }
  
  # Validate specified NA‐max.NA.col
  if (!is.numeric(max.NA.col) || length(max.NA.col) != 1 ||
      max.NA.col < 0 || max.NA.col > 1) {
    stop("Data processing aborted: max.NA.col must be single numeric value between 0 and 1 (e.g. 0.5)")
  }
  
  # Validate specified training.neighborhoods
  if (!is.character(training.neighborhoods) || length(training.neighborhoods) != 1 ||
      !(training.neighborhoods %in% c("gaussian", "bubble"))) {
    stop("Data processing aborted: training.neighborhoods must be 'gaussian' or 'bubble'")
  }
  
  # Validate specified save.SOM.results
  if (!is.logical(save.SOM.results) || length(save.SOM.results) != 1 || is.na(save.SOM.results)) {
    stop("Data processing aborted: save.SOM.results must be single logical value (TRUE or FALSE)")
  }
  
  # Validate specified save.SOM.results.name
  if (save.SOM.results && !is.null(save.SOM.results.name)) {
    if (!is.character(save.SOM.results.name) ||
        length(save.SOM.results.name) != 1 ||
        is.na(save.SOM.results.name) ||
        trimws(save.SOM.results.name) == "") {
      stop("Data processing aborted: save.SOM.results.name must be non-empty character string (file path) if provided")
    }
    valid_ext <- tolower(tools::file_ext(save.SOM.results.name)) #extract extension
    if (valid_ext != "rdata") {
      stop("Data processing aborted: save.SOM.results.name must end with .Rdata") #abort if not .Rdata
    }
  }
  
  # Validate specified overwrite.SOM.results
  if (!is.logical(overwrite.SOM.results) || length(overwrite.SOM.results) != 1 || is.na(overwrite.SOM.results)) {
    stop("Data processing aborted: overwrite.SOM.results must be single logical value (TRUE or FALSE)")
  }
  
  # Validate specified verbose
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("Data processing aborted: verbose must be single logical value (TRUE or FALSE)")
  }
  
  # Validate message.N.replicates
  if (!is.numeric(message.N.replicates) || length(message.N.replicates) != 1 || message.N.replicates < 1 || (message.N.replicates %% 1 != 0)) {
    stop("Data processing aborted: message.N.replicates must be single positive integer (>= 1)")
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
    
    required_fields <- c("distance_matrix", "learning_values_list")
    if (!all(required_fields %in% names(SOM_results))) { #check if structure of SOM_results is correct
      stop("Data processing aborted: could not load SOM results (results do not contain expected objects 'distance_matrix' and 'learning_values_list') - check saved file or rerun SOM")
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
      if (n_removed <= 50) {
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
    
    # Remove non‑numeric columns
    for (i in seq_along(input_data)) {
      mat <- input_data[[i]] #extract data
      mat <- as.data.frame(mat, stringsAsFactors = FALSE) #convert to dataframe
      non_numeric_cols <- which(sapply(mat, function(x) is.factor(x) || is.character(x)))
      if (length(non_numeric_cols) > 0) { #print message if any columns were removed
        n_removed <- length(non_numeric_cols)
        n_total <- ncol(mat) #number of columns before removal
        if (n_removed <= 50) {
          messager(sprintf(
            "Removed %d of %d columns (variables) in dataset %s due to non-numeric type: %s",
            n_removed, n_total, input_data_names[[i]], paste(names(mat)[non_numeric_cols], collapse = ", ")
          ))
        } else {
          messager(sprintf(
            "Removed %d of %d columns (variables) in dataset %s due to non-numeric type",
            n_removed, n_total, input_data_names[[i]]
          ))
        }
        mat <- mat[, -non_numeric_cols, drop = FALSE] #drop non-numeric columns
      }
      if (ncol(mat) == 0) stop(sprintf( #stop with message if no columns remain
        "Data processing aborted: no columns (variables) remain in dataset %s after removing all non-numeric columns - check input data", 
        input_data_names[[i]]
      ))
      if (ncol(mat) == 1) stop(sprintf( #stop with message if only one column remains
        "Data processing aborted: only one column (variable) remains in dataset %s after removing all non-numeric columns - check input data", 
        input_data_names[[i]]
      )) 
      input_data[[i]] <- mat
    }
    
    # Filter by max.NA.row
    for (i in seq_along(input_data)) {
      mat <- input_data[[i]]
      na_frac <- rowMeans(is.na(mat))
      removed <- rownames(mat)[na_frac > max.NA.row]
      if (length(removed) > 0) { #show message if samples are removed
        if (length(removed) <= 50) {
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
      mat <- input_data[[i]]
      
      # Remove columns filled with only NA
      n_cols <- ncol(mat)
      all_na_cols <- vapply(mat, function(x) all(is.na(x)), logical(1)) #extract columns with only NA
      if (any(all_na_cols)) { #print message if any columns were removed
        n_removed <- sum(all_na_cols)
        if (n_removed <= 50) {
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
        if (is.null(dim(mat)) || ncol(mat) == 0) stop("Data processing aborted: no columns (variables) remain after removing all-NA columns - check input data")
        if (ncol(mat) == 1) stop("Data processing aborted: only one column (variable) remains after removing all-NA columns - check input data")
      }
      
      # Remove variables (columns) with > max.NA.col missing
      col_na_frac <- colMeans(is.na(mat))
      dropped_cols <- names(col_na_frac)[col_na_frac > max.NA.col]
      if (length(dropped_cols)) { #print message if any columns were removed
        n0 <- ncol(mat)
        n_dropped <- length(dropped_cols)
        if (n_dropped <= 50) {
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
          stop(sprintf(
            "Data processing aborted: only one colum (variable) remains in dataset %s after applying max.NA.col = %.2f - check input data or increase max.NA.col",
            input_data_names[i], max.NA.col
          ))
        }
      }
      
      # Remove variables (columns) with zero variance
      zero_var_cols <- vapply(mat, function(x) {
        variance <- var(x, na.rm = TRUE)
        is.finite(variance) && variance == 0
      }, logical(1))
      
      if (any(zero_var_cols)) {
        removed_cols <- names(which(zero_var_cols))
        n_removed <- length(removed_cols)
        total_cols <- ncol(mat)
        
        if (n_removed <= 50) {
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
        if (is.null(dim(mat)) || ncol(mat) == 0) stop("Data processing aborted: no columns (variables) remain after removing all columns with zero variance - check input data")
        if (ncol(mat) == 1) stop("Data processing aborted: only one column (variable) remains after removing all columns with zero variance - check input data")
      }
      
      # Remove all-NA rows again after column filtering
      all_na_rows <- rowSums(is.na(mat)) == ncol(mat)
      if (any(all_na_rows)) { #print message if any rows are removed
        removed_names <- rownames(mat)[all_na_rows]
        n_removed <- sum(all_na_rows)
        n_rows <- nrow(mat)
        if (n_removed <= 50) {
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
    
    # Normalize columns for each matrix
    for (i in seq_along(processed_data)) {
      mat <- processed_data[[i]]
      mat <- as.matrix(mat)
      mat <- apply(mat, 2, function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
      if (nrow(mat) == 0) stop(sprintf("Data processing aborted: dataset %s has no columns (variables) remaining after normalization - check input data", processed_names[i])) #check if all columns are gone after normalization
      if (is.null(dim(mat))) stop(sprintf("Data processing aborted: dataset %s has only one column (variable) remaining after normalization - check input data", processed_names[i])) #check if only one column remains after normalization
      if (nrow(mat) == 0) stop(sprintf("Data processing aborted: dataset %s has no row (sample) remaining after normalization - check input data", processed_names[i])) #check if all rows are gone after normalization
      if (nrow(mat) == 1) stop(sprintf("Data processing aborted: dataset %s has only one row (sample) remaining after normalization - check input data", processed_names[i])) #check if only one row remains after normalization
      
      # Remove any columns that are all NA after normalization
      mat <- as.data.frame(mat, stringsAsFactors = FALSE)
      all_na_cols <- vapply(mat, function(x) all(is.na(x)), logical(1))
      if (any(all_na_cols)) { #print message if any columns were removed
        n_removed <- sum(all_na_cols)
        removed_names <- colnames(mat)[all_na_cols]
        if (n_removed <= 50) {
          messager(sprintf("Removed %d columns (variables) in dataset %s after normalization because all values are NA: %s",
                           n_removed, processed_names[i], paste(removed_names, collapse = ", ")))
        } else {
          messager(sprintf("Removed %d columns (variables) in dataset %s after normalization because all values are NA",
                           n_removed, processed_names[i]))
        }
        mat <- mat[, !all_na_cols, drop = FALSE]
        if (is.null(dim(mat)) || ncol(mat) == 0) stop("Data processing aborted: no columns (variables) remain after removing all-NA columns - check input data or increase max.NA.col")
        if (ncol(mat) == 1) stop("Data processing aborted: only one column (variable) remains after removing all-NA columns - check input data or increase max.NA.col")
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
    
    # Remove non‑numeric columns
    non_numeric_cols <- which(sapply(mat, function(x) is.factor(x) || is.character(x)))
    if (length(non_numeric_cols) > 0) { #print message if any columns were removed
      n_removed <- length(non_numeric_cols)
      n_total <- ncol(mat) #number of columns before removal
      if (n_removed <= 50) {
        messager(sprintf(
          "Removed %d of %d columns (variables) in dataset %s due to non-numeric type: %s",
          n_removed,
          n_total,
          input_data_names[[1]],
          paste(names(mat)[non_numeric_cols], collapse = ", ")
        ))
      } else {
        messager(sprintf(
          "Removed %d of %d columns (variables) in dataset %s due to non-numeric type",
          n_removed,
          n_total,
          input_data_names[[1]]
        ))
      }
      mat <- mat[, -non_numeric_cols, drop = FALSE] #drop non-numeric columns
    }
    if (ncol(mat) == 0) stop("Data processing aborted: no columns (variables) remain after removing all non-numeric columns - check input data")
    if (ncol(mat) == 1) stop("Data processing aborted: only one column (variable) remains after removing all non-numeric columns - check input data")    
    
    # Remove rows with > max.NA.row missing
    na_props <- rowMeans(is.na(mat)) #fraction of NA per row
    bad_samples <- rownames(mat)[na_props > max.NA.row] #extract rows exceeding threshold
    if (length(bad_samples) > 0) { #print message if rows (samples) are removed
      if (length(bad_samples) <= 50) {
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
      if (n_removed <= 50) {
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
      variance <- var(x, na.rm = TRUE)
      is.finite(variance) && variance == 0
    }, logical(1)) #extract zero-variance columns
    if (any(zero_var_cols)) { #print message if columns are removed
      removed_cols <- names(which(zero_var_cols))
      n_removed <- length(removed_cols)
      if (n_removed <= 50) {
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
      if (n_dropped <= 50) {
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
      if (n_removed <= 50) {
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
      if (n_removed <= 50) {
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
    
    # Wrap matrix as list for consistency with multi-dataset branch
    input_data <- list(as.matrix(mat)) #final output is list of matrix
    
  }
  
  # Report number of rows used as SOM input
  messager(sprintf("Data processing completed - using %d samples for SOM", nrow(input_data[[1]])))
  
  # Create SOM output grid
  messager("")
  messager("CREATING SOM GRID ...")
  
  if (!is.null(grid.size)) {
    
    # Use user-specified grid dimensions
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
  } else {
    
    # Use default calculation based on number of samples and aspect ratio (MFA or PCA)
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
        na_row_mask <- rowSums(is.na(combined)) < ncol(combined) #identify rows not all NA
        if (sum(na_row_mask) < 5) stop("Not enough non-NA rows for MFA") #require at least 5 rows
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
      
      # Calculate SOM grid aspect ratio from input covariance matrix eigenvalues or if fails using mean-imputation on columns with less than 25% missingness
    } else {
      mat <- as.matrix(input_data[[1]]) #convert input data to matrix
      mat[is.infinite(mat)] <- NA #set infinite values to NA
      
      col_good <- which(colSums(!is.na(mat)) >= 2) #keep columns with at least 2 non-NA values
      mat1 <- mat[, col_good, drop = FALSE] #subset to those columns
      row_good <- which(rowSums(!is.na(mat1)) >= 2) #keep rows with at least 2 non-NA
      mat1 <- mat1[row_good, , drop = FALSE] #subset to those rows
      mat_cc <- mat1[complete.cases(mat1), , drop = FALSE] #subset to complete-case rows
      
      if (nrow(mat_cc) >= 2 && ncol(mat_cc) >= 2) { #enough complete rows/cols for covariance
        covariance_mat <- cov(mat_cc) #calculate covariance matrix
        eigenvalues <- tryCatch(eigen(covariance_mat)$values, error = function(e) c(1, 1)) #get eigenvalues or fallback
        if (length(eigenvalues) >= 2 && all(is.finite(eigenvalues)) && all(eigenvalues > 0)) {
          SOM_grid_aspect_ratio <- sqrt(eigenvalues[1] / eigenvalues[2]) #aspect ratio from eigenvalues
        } else {
          SOM_grid_aspect_ratio <- NA #flag to try next step
        }
      } else {
        SOM_grid_aspect_ratio <- NA #flag to try next step
      }
      
      if (is.na(SOM_grid_aspect_ratio)) { #mean-imputation step: only try if missing data not excessive (<25% missing overall)
        col_missing_prop <- colMeans(is.na(mat)) #fraction of missing data per column
        col_good2 <- which(col_missing_prop < 0.25 & colSums(!is.na(mat)) >= 2) #columns with <25% missing and at least 2 non-NA
        row_good2 <- which(rowSums(!is.na(mat[, col_good2, drop = FALSE])) >= 2) #rows with at least 2 non-NA in these columns
        mat2 <- mat[row_good2, col_good2, drop = FALSE] #subset matrix to good rows and columns
        if (nrow(mat2) >= 2 && ncol(mat2) >= 2) { #proceed if enough data
          for (j in seq_len(ncol(mat2))) { #iterate over columns for imputation
            idx <- which(is.na(mat2[, j])) #find missing values in column
            if (length(idx) > 0 && sum(!is.na(mat2[, j])) > 0) mat2[idx, j] <- mean(mat2[, j], na.rm = TRUE) #impute mean if at least one value
          }
          covariance_mat2 <- cov(mat2) #calculate covariance matrix
          eigenvalues2 <- tryCatch(eigen(covariance_mat2)$values, error = function(e) c(1, 1)) #extract eigenvalues, fallback if error
          if (length(eigenvalues2) >= 2 && all(is.finite(eigenvalues2)) && all(eigenvalues2 > 0)) {
            SOM_grid_aspect_ratio <- sqrt(eigenvalues2[1] / eigenvalues2[2]) #compute aspect ratio
          } else {
            SOM_grid_aspect_ratio <- 1 #fallback to square grid
            messager("Aspect ratio calculation failed (not enough samples or too much NA) - using square/square-like default aspect ratio for SOM grid") #warn if fails
          }
        } else {
          SOM_grid_aspect_ratio <- 1 #fallback to square grid if not enough usable data
          messager("Aspect ratio calculation failed (not enough samples or too much NA) - using square/square-like default aspect ratio for SOM grid") #warn if fails
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
  
  # Perform learning rate tuning
  if (learning.rate.tuning) {
    messager("")
    messager("TUNING LEARNING RATES ...")
    
    # Set tuning values and number of runs
    learning_rate_initial_tuning_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    learning_rate_final_tuning_values <- c(0.001, 0.01, 0.1, 0.3, 0.5)
    learning_rate_N_steps <- 80
    learning_rate_N_replicates <- 10
    
    tuning_values <- expand.grid(learning_rate_initial = learning_rate_initial_tuning_values,
                                 learning_rate_final = learning_rate_final_tuning_values,
                                 stringsAsFactors = FALSE)
    tuning_values <- tuning_values[tuning_values$learning_rate_final < tuning_values$learning_rate_initial, , drop = FALSE] #keep only valid combinations: final learning rate < initial learning rate
    tuning_values$qe_mean <- NA_real_
    tuning_values$qe_sd <- NA_real_
    
    for (i in seq_len(nrow(tuning_values))) {
      learning_rate_initial <- tuning_values$learning_rate_initial[i]
      learning_rate_final <- tuning_values$learning_rate_final[i]
      
      # Run learning_rate_N_replicates replicates and record mean QE (quantification error) each time
      set.seed(i + 500)
      qes <- replicate(learning_rate_N_replicates, {
        som_tuning <- kohonen::supersom(data = input_data,
                                        grid = SOM_output_grid,
                                        maxNA.fraction = max.NA.row,
                                        alpha = c(learning_rate_initial, learning_rate_final),
                                        rlen = learning_rate_N_steps)
        mean(som_tuning$distances)
      })
      
      tuning_values$qe_mean[i] <- mean(qes)
      tuning_values$qe_sd[i] <- sd(qes)
    }
    tuning_values_sorted <- tuning_values[order(tuning_values$qe_mean), ]
    tuning_values_best <- tuning_values_sorted[1, ]
    messager(sprintf("Best tuning parameters: initial = %.2f, final = %.2f",
                     tuning_values_best$learning_rate_initial, tuning_values_best$learning_rate_final
    ))
    learning.rate.initial <- tuning_values_best$learning_rate_initial
    learning.rate.final <- tuning_values_best$learning_rate_final
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
                                   rlen = N.steps)
    
    # Store learning values for each matrix
    learning_values_list <- lapply(seq_along(input_data), function(i) som_model$changes[, i])
    
    # Store distance weights
    d_vec <- som_model$distance.weights 
    
    # Return results
    return(list(
      d_vec = d_vec, 
      learning_values_list = learning_values_list, 
      som_model = som_model
    ))
  }
  
  # Perform SOM training and collect results for all replicates
  
  # Run in parallel
  if (parallel) { 
    messager(sprintf("Running SOM in parallel with %d cores", N_cores))
    parallel_cluster <- parallel::makeCluster(N_cores) #start PSOCK cluster
    on.exit(parallel::stopCluster(parallel_cluster), add = TRUE)
    parallel::clusterExport( #export helper and all needed globals to each worker
      parallel_cluster,
      varlist = c("replicate_som",
                  "input_data",
                  "SOM_output_grid",
                  "N.replicates",
                  "N.steps",
                  "max.NA.row",
                  "learning.rate.initial",
                  "learning.rate.final",
                  "message.N.replicates"),
      envir = environment())
    doParallel::registerDoParallel(parallel_cluster) #register cluster for foreach
    doRNG::registerDoRNG(seed = 1) #set seed
    results <- foreach::foreach( #launch in parallel (loads all your required_packages automatically)
      j = seq_len(N.replicates),
      .packages = required_packages
    ) %dopar% {
      replicate_som(j)
    }
  } else {
    
    # Run SOM normally (non-parallel)
    set.seed(1)
    results <- tryCatch(
      lapply(seq_len(N.replicates), function(j) {
        set.seed(j + 1e4) # set unique seed for each replicate
        replicate_som(j)
      }),
      error = function(e) { #print error message if SOM training fails
        stop("SOM training aborted: try reducing grid.size or increasing max.NA.col and max.NA.row or check input data")
      }
    )
    if (is.null(results)) return(invisible(NULL))
  }
  
  # Combine results from all replicates
  distance_matrix <- do.call(rbind, lapply(results, `[[`, "d_vec"))
  rownames(distance_matrix) <- paste0("R", seq_len(nrow(distance_matrix)))
  colnames(distance_matrix) <- as.character(input_data_names)
  
  learning_values_list <- lapply(1:length(input_data), function(i) {
    learning_values_combined <- do.call(cbind, lapply(results, function(res) res$learning_values_list[[i]]))
    rownames(learning_values_combined) <- paste0("S", seq_len(N.steps))
    colnames(learning_values_combined) <- paste0("R", seq_len(ncol(learning_values_combined)))
    return(learning_values_combined)
  })
  
  som_models <- lapply(results, `[[`, "som_model")
  
  # Compute variable importance for each layer
  n_layers <- length(som_models[[1]]$codes)
  all_layer_codes <- vector("list", n_layers)
  for (l in seq_len(n_layers)) {
    all_layer_codes[[l]] <- do.call(
      rbind,
      lapply(som_models, function(m) m$codes[[l]])
    )
  }
  
  # Extract median variable importance for each layer
  median_variable_importance <- vector("list", n_layers)
  for (l in seq_len(n_layers)) {
    layer_codes <- all_layer_codes[[l]] #extract all codes for this layer across all replicates
    median_variable_importance[[l]] <- apply(layer_codes, 2, median, na.rm = TRUE)
  }
  if (!is.null(colnames(all_layer_codes[[1]]))) { #add variable names if available
    for (l in seq_len(n_layers)) {
      names(median_variable_importance[[l]]) <- colnames(all_layer_codes[[l]])
    }
  }
  
  # Save results
  SOM_results <- list(
    distance_matrix = distance_matrix, 
    learning_values_list = learning_values_list, 
    input_data_names = as.character(input_data_names),
    N_steps = N.steps,
    N_replicates = N.replicates,
    learning_rate_initial = learning.rate.initial,
    learning_rate_final = learning.rate.final,
    codebook_vectors = all_layer_codes,
    median_variable_importance = median_variable_importance,
    som_models = som_models
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
  
  return(SOM_results)
}


## Function to cluster SOM codebook vectors
clustering.SOM <- function(SOM.output,
                           max.k = 7, #maximum of considered clusters K + 1
                           set.k = NULL, #set to test single value of K
                           clustering.method, #set clustering method
                           pca.codebooks = TRUE, #run PCA on codebook vectors
                           pca.codebooks.var.threshold = 0.99, #variance threshold for codebook PCA to choose how many PCs to retain
                           BIC.thresh = 6 #BIC threshold for selecting K > 1 - we suggest using Raftery (1995) ranges: 2, 6, or 10 for weak, medium or strong support
) {
  
  
  # Validate specified SOM.output
  if (!is.list(SOM.output) || is.null(SOM.output$som_models) || length(SOM.output$som_models) < 1) {
    stop("Aborted SOM clustering: SOM.output must be list from train.SOM() with non-empty $som_models")
  }
  
  # Validate specified max.k
  if (!is.numeric(max.k) || length(max.k) != 1 || is.na(max.k) || max.k < 2 || (max.k %% 1 != 0)) {
    stop("Aborted SOM clustering: max.k must be single integer >= 2")
  }
  
  # Validate specified set.k
  if (!is.null(set.k) && (!is.numeric(set.k) || length(set.k) != 1 || set.k < 1 || (set.k %% 1 != 0))) {
    stop("Aborted SOM clustering: set.k must be NULL or single positive integer >= 1")
  }
  
  # Validate specified clustering method
  valid.methods <- c(
    "GMM+BIC",
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
    "PAM+gapstat"
  )
  if (!clustering.method %in% valid.methods) {
    stop("Aborted SOM clustering: Invalid clustering.method - must be one of ", paste(valid.methods, collapse = ", "))
  }
  
  # Validate specified pca.codebooks
  if (!is.logical(pca.codebooks) || length(pca.codebooks) != 1 || is.na(pca.codebooks)) {
    stop("Aborted SOM clustering: pca.codebooks must be TRUE or FALSE")
  }
  
  # Validate specified pca.codebooks.var.threshold
  if (!is.numeric(pca.codebooks.var.threshold) || pca.codebooks.var.threshold < 0 || pca.codebooks.var.threshold > 1) {
    stop("Aborted SOM clustering: pca.codebooks.var.threshold must be numeric value between 0 and 1")
  }
  
  # Validate specified BIC.thresh
  if (!is.numeric(BIC.thresh) || length(BIC.thresh) != 1 || is.na(BIC.thresh) || BIC.thresh <= 0) {
    stop("Aborted SOM clustering: BIC.thresh must be single positive numeric value (e.g., 2, 6, or 10 for low, moderate or strong support, respectively)")
  }
  
  # Extract number of replicates
  N.replicates <- length(SOM.output$som_models)
  
  # Create function to cluster SOM models
  replicate_clust <- function(j) {
    
    # Set seed
    set.seed(j + 1000)
    
    # Extract SOM models
    som_model <- SOM.output$som_models[[j]] 
    
    # Extract code book vectors
    codes <- getCodes(som_model)
    if (!is.list(codes)) codes <- list(codes)
    som_codes <- do.call(cbind, codes)
    rownames(som_codes) <- paste0("G", seq_len(nrow(som_codes)))
    
    # Perform PCA on codebook vectors
    if (pca.codebooks) {
      nonconst <- apply(som_codes, 2, function(x) { #remove all‐constant codebook columns to allow PCA
        var(x, na.rm = TRUE) > 0
      })
      if (!all(nonconst)) {
        som_codes <- som_codes[, nonconst, drop = FALSE]
      }
      if (ncol(som_codes) == 0) {stop("Aborted SOM clustering: no variable left after removing constant codebook columns - check input data")
      }
      
      if (ncol(som_codes) == 1) {
        message("Only one non‑constant codebook dimension left - skipping PCA on codebook vectors")
      } else {
        pca_res <- prcomp(som_codes, scale. = TRUE) #perform PCA on codebook vectors
        explained_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
        if (any(explained_var >= pca.codebooks.var.threshold)) {
          keep_pcs <- which(explained_var >= pca.codebooks.var.threshold)[1] #number of PCs to keep
        } else {
          keep_pcs <- length(explained_var) #if no PC cross threshold, keep all
        }
        som_pca <- pca_res$x[, 1:keep_pcs, drop = FALSE]
        som_codes <- som_pca
        if (ncol(som_codes) < 2 || nrow(som_codes) < 2) {
          stop(sprintf("Aborted SOM clustering: codebook matrix has %d rows and %d columns after PCA on codebook vectors and constant-column removal (need at least 2 x 2)", nrow(som_codes), ncol(som_codes)))
        }
      }
    }
    
    # Ensure max.k is equal to or smaller than number of available codebook vectors (rows of som_codes)
    n_codes <- nrow(som_codes)
    if (max.k >= n_codes) {
      stop(sprintf(
        "Aborted SOM clustering: max.k = %d exceeds available codebook rows of %d after PCA - reduce max.k to ≤ %d",
        max.k, n_codes, n_codes - 1
      ))
    }
    
    # Create function to perform kmeans clustering and calculate within‐cluster sum of squares (wss) for each cluster (sum of squared Euclidean distances of SOM units to cluster center)
    calculate.wss <- function(som_codes, max.k) {
      wss <- numeric(max.k) #wss vector for k = 1 ... max.k
      wss[1] <- (nrow(som_codes) - 1) * sum(apply(som_codes, 2, var)) #wss for k = 1 (= total sum of squared distances to overall mean)
      wss[2:max.k] <- sapply(2:max.k, function(i) sum(kmeans(som_codes, centers = i, nstart = 30, iter.max = 1e5)$withinss)) #wss for k = 2 ... max.k via kmeans
      wss
    }
    
    # Create function to calculate BIC with number of parameters (≈ dimensions × clusters) penalty for each cluster using WSS (Pelleg & Moore 2000. X‐means: Extending K‐means with efficient estimation of the number of clusters. In ICML’00 (pp. 727‐734), Citeseer)
    calculate.wssBIC <- function (wss, som_codes) {
      N <- nrow(som_codes) #number of output cells
      k_vals <- seq_along(wss)
      BIC_vec <- N * log(wss / N) + log(N) * k_vals #BIC for each k
      BIC_vec
    }
    
    # Create function to determine optimal number of clusters by selecting smallest k where BIC improvement falls below BIC threshold, otherwise, choose k with minimum BIC
    select.k.BICthresh <- function(BIC_vec, BIC.thresh, set.k = NULL) {
      if (!is.null(set.k)) return(set.k) #user-specified number of clusters
      som_N_clusters <- NA #store optimal k
      for (k in 2:length(BIC_vec)) {
        if ((BIC_vec[k] - BIC_vec[k - 1]) > -BIC.thresh) { #if improvement is not at least as large as threshold, select previous k
          som_N_clusters <- k - 1
          break
        }
      }
      if (is.na(som_N_clusters)) {
        som_N_clusters <- which.min(BIC_vec) #if all differences exceed threshold, pick k with lowest BIC
      }
      som_N_clusters
    }
    
    # Create function to determine optimal number of clusters using BIC jump rule
    select.k.BICjump <- function(BIC_vec, BIC.thresh, set.k = NULL) {
      if (!is.null(set.k)) return(set.k) #user-specified k
      if (length(BIC_vec) < 2) return(1) # if there is only one or no BIC value, set default to k = 1
      if (abs(diff(BIC_vec))[1] < BIC.thresh) { #computes first ΔBIC between k = 1 and k = 2 and if that drop is smaller than BIC threshold, there is no real improvement from adding second cluster, so fix k = 1
        som_N_clusters <- 1 #k = 1
      } else {
        if (length(BIC_vec) <= 2) { #if only 2 BIC values
          som_N_clusters <- 2 #k = 2
        } else {
          deltaBIC_clusters <- cutree(hclust(dist(diff(BIC_vec)), method = "ward.D"), k = 2) #cluster ΔBIC values into two groups using hierarchical clustering
          best_group <- which.min(tapply(diff(BIC_vec), deltaBIC_clusters, mean)) #group with smallest mean ΔBIC
          som_N_clusters <- max(which(deltaBIC_clusters == best_group)) + 1 #select optimal k
        }
      }
      som_N_clusters
    }
    
    # Create function to assign clusters using spectral clustering
    assign.clusters.spectral.clustering <- function(som_codes, som_N_clusters) { 
      if (som_N_clusters < 2) { #k = 1
        return(rep(1, nrow(som_codes)))
      } else {
        spec_clust_result <- try(kernlab::specc(as.matrix(som_codes), centers = som_N_clusters), silent = TRUE) #perform default spectral clustering
        if (!inherits(spec_clust_result, "Specc")) { #if default failed, try different sigma values for spectral clustering
          sigma_grid <- c(0.001, 0.01, 0.1, 1, 10, 100, 500, 1000)
          for (sig in sigma_grid) {
            spec_clust_result <- try(
              kernlab::specc(
                as.matrix(som_codes),
                centers = som_N_clusters,
                kernel = "rbfdot",
                kpar   = list(sigma = sig)
              ),
              silent = TRUE
            )
            if (inherits(spec_clust_result, "Specc")) break
          }
        }
        if (!inherits(spec_clust_result, "Specc")) { #stop if spectral clustering fails even with all tested sigma values
          stop("Aborted SOM clustering: spectral clustering failed - check input data or try different clustering.method")
        }
        return(as.integer(spec_clust_result))
      }
    } 
    
    # Perform clustering of SOM codebook vectors based on specified method
    
    # Clustering method: kmeans + BICjump + spectral clustering
    if (clustering.method == "kmeans+BICjump+spectral") {
      wss <- calculate.wss(som_codes, max.k) #perform kmeans and calculate within‐cluster sum of squares (wss)
      BIC_vec <- calculate.wssBIC(wss, som_codes) #calculate BIC based on wss
      som_N_clusters <- select.k.BICjump(BIC_vec, BIC.thresh, set.k) #determine optimal number of clusters by selecting smallest k where BIC improvement falls below BIC threshold, otherwise, choose k with minimum BIC
      som_cluster <- assign.clusters.spectral.clustering(som_codes, som_N_clusters) #assign clusters using spectral clustering
    }
    
    # Clustering method: kmeans + BICjump
    if (clustering.method == "kmeans+BICjump") {
      wss <- calculate.wss(som_codes, max.k) #perform kmeans and calculate within‐cluster sum of squares (wss)
      BIC_vec <- calculate.wssBIC(wss, som_codes) #calculate BIC based on wss
      som_N_clusters <- select.k.BICjump(BIC_vec, BIC.thresh, set.k) #determine optimal number of clusters using BIC jump rule
      som_cluster <- kmeans(som_codes, centers = som_N_clusters, nstart = 30, iter.max = 1e5)$cluster #assign clusters using k-means with selected k
    }
    
    # Clustering method: kmeans + BICthresh
    if (clustering.method == "kmeans+BICthresh") {
      wss <- calculate.wss(som_codes, max.k) #perform kmeans and calculate within‐cluster sum of squares (wss)
      BIC_vec <- calculate.wssBIC(wss, som_codes) #calculate BIC based on wss
      som_N_clusters <- select.k.BICthresh(BIC_vec, BIC.thresh, set.k) #determine optimal number of clusters using BIC threshold rule (selecting smallest k where BIC improvement falls below BIC threshold, otherwise, choose k with minimum BIC)
      som_cluster <- kmeans(som_codes, centers = som_N_clusters, nstart = 30, iter.max = 1e5)$cluster #assign clusters using k-means with selected k
    }
    
    # Clustering method: kmeans + BICthresh + spectral clustering
    if (clustering.method == "kmeans+BICthresh+spectral") {
      wss <- calculate.wss(som_codes, max.k) #perform kmeans and calculate within‐cluster sum of squares (wss)
      BIC_vec <- calculate.wssBIC(wss, som_codes) #calculate BIC based on wss
      som_N_clusters <- select.k.BICthresh(BIC_vec, BIC.thresh, set.k) #determine optimal number of clusters using BIC threshold rule (selecting smallest k where BIC improvement falls below BIC threshold, otherwise, choose k with minimum BIC)
      som_cluster <- assign.clusters.spectral.clustering(som_codes, som_N_clusters) #assign clusters using spectral clustering
    }
    
    # Clustering method: kmeans + BICjump + hierarchical clustering
    if (clustering.method == "kmeans+BICjump+hierarch") {
      wss <- calculate.wss(som_codes, max.k) #perform kmeans and calculate within‐cluster sum of squares (wss)
      BIC_vec <- calculate.wssBIC(wss, som_codes) #calculate BIC based on wss
      som_N_clusters <- select.k.BICjump(BIC_vec, BIC.thresh, set.k) #determine optimal number of clusters using hierarchical clustering based on BIC values
      dist_som_codes <- dist(som_codes) #compute pairwise Euclidean distances between all rows of som_codes
      hierarchical_clust_som_codes <- hclust(dist_som_codes) #perform hierarchical (agglomerative) clustering on distances building dendrogram that successively merges closest clusters
      som_cluster <- cutree(hierarchical_clust_som_codes, som_N_clusters) #cut dendrogram at height which produces exactly som_N_clusters group, assigning each neuron to one of those clusters
    }
    
    # Clustering method: GMM (Gaussian mixture model) + BIC
    if (clustering.method == "GMM+BIC") { #Gaussian Mixture Model (GMM) with best K based on lowest BIC
      mclust_fit <- mclust::Mclust(som_codes, G = 1:max.k, verbose = FALSE) #fit probabilistic mixture of Gaussian (normal) distributions to codebook vectors for each k
      som_N_clusters <- mclust_fit$G #optimal k based on lowest BIC
      som_cluster <- mclust_fit$classification #cluster assignment for each codebook vector
      BIC_vec <- mclust_fit$BIC[1, 1:max.k] #BIC for each k
    }
    
    # Clustering method: GMM (Gaussian mixture model) + BIC threshold
    if (clustering.method == "GMM+BICthresh") { #GMM with best K based on selecting smallest k where BIC improvement falls below BIC threshold, otherwise, choose k with minimum BIC
      mclust_fit_all <- mclust::Mclust(som_codes, G = 1:max.k, verbose = FALSE) #fit GMMs for all k
      BIC_vec <- mclust_fit_all$BIC[1, 1:max.k] #extract BIC values
      som_N_clusters <- select.k.BICthresh(BIC_vec, BIC.thresh, set.k) #determine optimal number of clusters using BIC threshold rule (selecting smallest k where BIC improvement falls below BIC threshold, otherwise, choose k with minimum BIC)
      som_cluster <- mclust::Mclust(som_codes, G = som_N_clusters, verbose = FALSE)$classification #rerun GMM with determined k and extract cluster assignment
    }
    
    # Clustering method: kmeans + gap statistics
    if (clustering.method == "kmeans+gapstat") { #gap statistic with k-means
      set.seed(j + 10000)
      gap <- cluster::clusGap(som_codes, B = 500, verbose = FALSE, K.max = max.k, FUN = function(x, k) list(cluster = kmeans(x, centers = k, nstart = 30)$cluster)) #run kmeans and calculate gap statistic (within-cluster dispersion, i.e. how tight clusters are)
      gap_stats <- gap$Tab #extract gap statistic table
      gap_diff <- gap_stats[-nrow(gap_stats), "gap"] - gap_stats[-1, "gap"] + gap_stats[-1, "SE.sim"] #pick smallest k where gap does not improve more than its standard error if you add another cluster (Tibshirani et al. 2001 https://doi.org/10.1111/1467-9868.00293)
      gap_k <- which.max(gap_diff >= 0) #k where this condition is first met
      if (length(gap_k) == 0) gap_k <- which.max(gap_stats[, "gap"]) #if none is found, falls back to k with maximum gap statistic
      som_N_clusters <- gap_k
      som_cluster <- kmeans(som_codes, centers = som_N_clusters, nstart = 30)$cluster #cluster assignment for each codebook vector using k-means with selected k
      BIC_vec <- rep(NA_real_, max.k) #gap statistic does not use BIC so fill with NA
    }
    
    # Clustering method: PAM (Partitioning Around Medoids) + gap statistics
    if (clustering.method == "PAM+gapstat") { #gap statistic with PAM (Partitioning Around Medoids)
      set.seed(j + 10000)
      gap <- cluster::clusGap(som_codes, K.max = max.k, B = 500, verbose = FALSE, FUN = function(x, k) list(cluster = cluster::pam(x, k, cluster.only = TRUE))) #use PAM for clustering (returns cluster assignments only)
      gap_stats <- gap$Tab #extract gap statistic table
      gap_diff <- gap_stats[-nrow(gap_stats), "gap"] - gap_stats[-1, "gap"] + gap_stats[-1, "SE.sim"] #pick smallest k where gap does not improve more than its standard error if you add another cluster (Tibshirani et al. 2001 https://doi.org/10.1111/1467-9868.00293)
      gap_k <- which.max(gap_diff >= 0) #k where this condition is first met
      if (length(gap_k) == 0) gap_k <- which.max(gap_stats[, "gap"]) #if none is found, fall back to k with maximum gap statistic (global max)
      som_N_clusters <- gap_k
      som_cluster <- cluster::pam(som_codes, k = som_N_clusters)$clustering #rerun PAM with selected k and assign clusters
      BIC_vec <- rep(NA_real_, max.k) #gap statistic does not use BIC so fill with NA
    }
    
    # Clustering method: hierarchical clustering + gap statistics
    if (clustering.method == "hierarchical+gapstat") { #gap statistic with hierarchical clustering
      set.seed(j + 10000)
      gap <- cluster::clusGap(som_codes, FUN = function(x, k) list(cluster = cutree(hclust(dist(x), method = "complete"), k = k)), verbose = FALSE, K.max = max.k, B = 500) #run hierarchical clustering and obtain gap statistic
      gap_stats <- gap$Tab #extract gap statistic table
      gap_diff <- gap_stats[-nrow(gap_stats), "gap"] - gap_stats[-1, "gap"] + gap_stats[-1, "SE.sim"] #pick smallest k where gap does not improve more than its standard error if you add another cluster (Tibshirani et al. 2001 https://doi.org/10.1111/1467-9868.00293)
      gap_k <- which.max(gap_diff >= 0) #k where this condition is first met
      if (length(gap_k) == 0) gap_k <- which.max(gap_stats[, "gap"]) #if none is found, fall back to k with maximum gap statistic (global max)
      som_N_clusters <- gap_k
      som_cluster <- cutree(hclust(dist(som_codes), method = "complete"), k = som_N_clusters) #rerun hierarchical clustering using selected k and assign clusters
      BIC_vec <- rep(NA_real_, max.k) #gap statistic does not use BIC so fill with NA
    }
    
    # Clustering method: kmeans + BIC threshold + hierarchical clustering
    if (clustering.method == "kmeans+BICthresh+hierarch") {
      wss <- calculate.wss(som_codes, max.k) #perform kmeans and calculate within‐cluster sum of squares (wss)
      BIC_vec <- calculate.wssBIC(wss, som_codes) #calculate BIC based on wss
      som_N_clusters <- select.k.BICthresh(BIC_vec, BIC.thresh, set.k) #determine optimal number of clusters by selecting smallest k where BIC improvement falls below BIC threshold, otherwise, choose k with minimum BIC
      dist_som_codes <- dist(som_codes) #compute pairwise Euclidean distances between all rows of som_codes
      hierarchical_clust_som_codes <- hclust(dist_som_codes) #perform hierarchical (agglomerative) clustering on distances building dendrogram that successively merges closest clusters
      som_cluster <- cutree(hierarchical_clust_som_codes, som_N_clusters) #cut dendrogram at height which produces exactly som_N_clusters group, assigning each neuron to one of those clusters
    }
    
    # Clustering method: GMM (Gaussian mixture model) + gap statistics
    if (clustering.method == "GMM+gapstat") { #Gap statistic with GMM, include k = 1 permutation test
      set.seed(j + 10000)
      gap <- cluster::clusGap(som_codes, FUN = function(x, k) list(cluster = mclust::Mclust(x, G = k, verbose = FALSE)$classification), K.max = max.k, verbose = FALSE, B = 500) #calculate gap statistic based on GMM
      gap_stats <- gap$Tab #extract gap statistic table
      gap_diff <- gap_stats[-nrow(gap_stats), "gap"] - gap_stats[-1, "gap"] + gap_stats[-1, "SE.sim"] #pick smallest k where gap does not improve more than its SE if you add another cluster
      gap_k <- which.max(gap_diff >= 0) #k where this condition is first met
      if (length(gap_k) == 0) gap_k <- which.max(gap_stats[, "gap"]) #if none found, fallback to k with global max gap
      som_N_clusters <- gap_k
      if (som_N_clusters == 2) { #if k = 2, permutation test for k = 2 vs k = 1
        set.seed(1)
        fit_obs <- mclust::Mclust(som_codes, G = 2, verbose = FALSE)
        wss_obs <- sum(dist(som_codes[fit_obs$classification == 1, ])^2) + sum(dist(som_codes[fit_obs$classification == 2, ])^2) #within-cluster sum of squares for observed
        n_perm <- 500 #number of permutations
        wss_null <- numeric(n_perm)
        for (i in seq_len(n_perm)) { #run permutation test
          permuted <- apply(som_codes, 2, sample) #permute each column
          fit_perm <- mclust::Mclust(permuted, G = 2, verbose = FALSE)
          wss_null[i] <- sum(dist(permuted[fit_perm$classification == 1, ])^2) + sum(dist(permuted[fit_perm$classification == 2, ])^2)
        }
        pval <- mean(wss_null <= wss_obs) #lower WSS is better, so pval is permuted <= observed
        if (pval >= 0.05) {
          som_N_clusters <- 1
          som_cluster <- rep(1, nrow(som_codes))
        } else {
          som_cluster <- fit_obs$classification
        }
      } else {
        fit <- mclust::Mclust(som_codes, G = som_N_clusters, verbose = FALSE)
        som_cluster <- fit$classification
      }
      BIC_vec <- rep(NA_real_, max.k) #fill with NA
    }
    
    # Clustering method: HDBSCAN (Hierarchical Density-Based Spatial Clustering of Applications with Noise)
    if (clustering.method == "HDBSCAN") { #use hdbscan via dbscan package
      
      # Determine optimal minPts based on cluster quality (mean membership probability)
      max_minPts <- floor(nrow(som_codes) / 2) #set max_minPTs
      minPts_vals <- seq(4, max_minPts) #determine range of minPts
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
        som_cluster <- rep(1, nrow(som_codes)) #assign all points to one cluster
        som_N_clusters <- 1
        BIC_vec <- rep(NA_real_, max.k)
      } else { #run best model and evaluate reassignment strategies
        
        # Run best HDBSCAN model with selected minPts
        hdbscan_model_results <- hdbscan_model_results[valid_HDBSCAN, , drop = FALSE]
        best_minPts_row <- which.max(hdbscan_model_results$mean_mem)
        best_minPts <- hdbscan_model_results$minPts[best_minPts_row]
        best_hdbscan_model <- dbscan::hdbscan(som_codes, minPts = best_minPts)
        base_clusters <- best_hdbscan_model$cluster
        
        # Identify core and noise points
        noise_idx <- which(base_clusters == 0) #noise points (cluster = 0)
        core_idx <- which(base_clusters != 0) #core points (real clusters)
        
        # Assign noise points
        if (length(noise_idx) == 0 || length(core_idx) == 0) { #no reassignment needed
          som_cluster <- base_clusters
          som_N_clusters <- length(unique(base_clusters[base_clusters != 0]))
        } else {
          
          # Strategy 1: assign noise points to nearest cluster
          hdbscan_nn <- FNN::get.knnx(som_codes[core_idx, , drop = FALSE], som_codes[noise_idx, , drop = FALSE], k = 1)
          nearest_clusters <- base_clusters[core_idx][hdbscan_nn$nn.index]
          clusters_nearest <- base_clusters
          clusters_nearest[noise_idx] <- nearest_clusters
          k_nearest <- length(unique(clusters_nearest[clusters_nearest != 0]))
          sil_nearest <- NA_real_
          if (k_nearest >= 2) {
            sil_nearest <- tryCatch({
              mean(cluster::silhouette(clusters_nearest, dist(som_codes))[, 3])
            }, error = function(e) NA_real_)
          }
          
          # Strategy 2: assign distant noise to new singleton clusters
          dist_thresh <- quantile(hdbscan_nn$nn.dist, 0.95) #distance threshold for singleton assignment
          clusters_singleton <- base_clusters
          next_cluster <- max(base_clusters)
          for (i in seq_along(noise_idx)) {
            ni <- noise_idx[i]
            dist_to_core <- hdbscan_nn$nn.dist[i]
            if (dist_to_core > dist_thresh) {
              next_cluster <- next_cluster + 1
              clusters_singleton[ni] <- next_cluster
            } else {
              clusters_singleton[ni] <- base_clusters[core_idx][hdbscan_nn$nn.index[i]]
            }
          }
          k_singleton <- length(unique(clusters_singleton[clusters_singleton != 0]))
          sil_singleton <- NA_real_
          if (k_singleton >= 2) {
            sil_singleton <- tryCatch({
              mean(cluster::silhouette(clusters_singleton, dist(som_codes))[, 3])
            }, error = function(e) NA_real_)
          }
          
          # Compare strategies and select best based on silhouette score (or cluster count fallback)
          if (is.na(sil_singleton) && is.na(sil_nearest)) {
            if (k_singleton > k_nearest) {
              som_cluster <- clusters_singleton
            } else {
              som_cluster <- clusters_nearest
            }
          } else if (is.na(sil_singleton) || sil_nearest >= sil_singleton) {
            som_cluster <- clusters_nearest
          } else {
            som_cluster <- clusters_singleton
          }
          
          # Extract number of clusters
          if (all(som_cluster == 0)) { #all noise
            som_cluster[] <- 1
            som_N_clusters <- 1
          } else {
            som_N_clusters <- length(unique(som_cluster))
          }
        }
        
        BIC_vec <- rep(NA_real_, max.k) #fill BIC with NA (not used for HDBSCAN)
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
    
    # Return results
    return(list(cluster_assignment = cluster_assignment, 
                BIC_vec = BIC_vec, 
                som_model = som_model,
                som_cluster = som_cluster,
                som_N_clusters = som_N_clusters,
                cluster_gridcell_assignments = cluster_gridcell_assignments
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
  
  optim_k_vals <- sapply(results, `[[`, "som_N_clusters")
  optim_k_vals <- t(as.matrix(optim_k_vals))
  rownames(optim_k_vals) <- "optim_k_vals"
  colnames(optim_k_vals) <- paste0("R", seq_len(N.replicates))
  optim_k_mean <- mean(optim_k_vals, na.rm = T)
  if (all(is.na(optim_k_vals))) {
    stop("Aborted SOM clustering: all optimal K values are NA - check input data") 
  }
  k_levels <- sort(unique(as.numeric(optim_k_vals)))
  optim_k_vals_counts <- as.numeric(table(factor(as.numeric(optim_k_vals), levels = k_levels)))
  optim_k_vals_props <- optim_k_vals_counts / sum(optim_k_vals_counts)
  optim_k_summary <- cbind(Count = optim_k_vals_counts,
                           Proportion = round(optim_k_vals_props, 2)
  )
  rownames(optim_k_summary) <- paste0("k", k_levels)
  
  som_models <- lapply(results, `[[`, "som_model")
  som_clusters <- lapply(results, `[[`, "som_cluster")
  cluster_gridcell_assignments <- lapply(results, `[[`, "cluster_gridcell_assignments")
  
  # Replace NA values with 0.5 for labeling
  processed_input_data <- results[[1]]$som_model$data[[1]]
  processed_input_data[is.na(processed_input_data)] <- 0.5
  
  # Preprocess input data and generate cluster labels
  cluster_labels <- do.call(
    cbind,
    lapply(1:max.k, function(x) kmeans(processed_input_data, x)$cluster)
  )
  rownames(cluster_labels) <- rownames(processed_input_data) #set row names for cluster labels
  colnames(cluster_labels) <- paste("k", 1:max.k, sep = '') #set column names
  
  # Extract SOM cluster assignments and filter replicates based on maximum K
  all_k <- apply(cluster_assignment, 2, max, na.rm = TRUE) #calculate maximum K for each replicate (column)
  if (length(which(optim_k_vals <= max.k)) == 0) stop("Aborted SOM clustering: no replicates have k ≤ max.k - increase max.k or check input data")
  assignment_matrix <- cluster_assignment[, which(optim_k_vals <= max.k), drop = FALSE] #filter to keep replicates with K <= max_k
  
  # Relabel across replicates with Hungarian algorithm
  for (i in seq_len(ncol(assignment_matrix))) {
    run_k <- max(assignment_matrix[, i], na.rm = TRUE)
    if (is.na(run_k) || run_k < 1) next
    ref <- cluster_labels[, run_k] #reference labels for this K
    obs <- assignment_matrix[, i] #labels from this replicate
    tab <- table(factor(ref, levels = 1:run_k), factor(obs, levels = 1:run_k)) #create table (rows: ref clusters, cols: obs clusters; with levels 1:run_k)
    if (nrow(tab) != ncol(tab)) { #pad tab to square if needed for Hungarian algorithm
      new_size <- max(nrow(tab), ncol(tab)) #determine new size for square table
      tmp <- matrix(0, nrow = new_size, ncol = new_size) #create square matrix of zeros
      tmp[1:nrow(tab), 1:ncol(tab)] <- tab #copy original table into upper-left
      tab <- tmp #update tab to padded table
    }
    cost <- max(tab) - tab #build cost matrix (max overlap = minimal cost)
    mapping <- clue::solve_LSAP(cost) #solve assignment with Hungarian algorithm
    obs_levels <- sort(unique(obs)) #get sorted unique cluster labels from obs
    mapping_full <- rep(NA, max(obs_levels)) #initialize mapping vector up to highest label
    mapping_full[obs_levels] <- mapping[1:length(obs_levels)] #fill mapping for present clusters
    assignment_matrix[, i] <- as.integer(mapping_full[obs]) #assign new cluster labels to samples
  }
  
  # Build ancestry from the re-labelled matrix
  k.max <- max(assignment_matrix, na.rm = TRUE)
  prop_list <- lapply(seq_len(nrow(assignment_matrix)), function(i) {
    prop.table(table(factor(assignment_matrix[i, ], levels = seq_len(k.max))))
  })
  ancestry_matrix <- do.call(rbind, prop_list)
  colnames(ancestry_matrix) <- paste0("Cluster_", seq_len(ncol(ancestry_matrix)))
  rownames(ancestry_matrix) <- rownames(assignment_matrix)
  
  # Save results
  SOM_results <- SOM.output
  SOM_results$cluster_assignment <- cluster_assignment
  SOM_results$BIC_values <- BIC_values
  SOM_results$ancestry_matrix <- ancestry_matrix
  SOM_results$optim_k_vals <- optim_k_vals
  SOM_results$optim_k_mean <- optim_k_mean
  SOM_results$optim_k_summary <- optim_k_summary
  SOM_results$max_k <- max.k
  SOM_results$set_k <- set.k
  SOM_results$som_clusters <- som_clusters
  SOM_results$cluster_gridcell_assignments <- cluster_gridcell_assignments
  
  return(SOM_results)
}


## Function to plot Structure-like barplots
plot.Structure.SOM <- function(SOM.output,
                               col.pal = viridis::viridis, #color palette
                               save = F, #save plot
                               overwrite = T, #overwrite plot if already present (only if saving plot)
                               plot.type = "svg", #plot file type (choose: png, jpg or svg; only if saving plot)
                               file.name = NULL, #plot file name (if NULL, default name is created; only if saving plot)
                               width = 10, #plot width in cm (only if saving plot)
                               height = 15, #plot height in cm (only if saving plot)
                               resolution = 300, #plot resolution in dpi (only if saving plot)
                               margin.bottom = 5, #bottom margin
                               margin.left = 5, #left margin
                               margin.top = 1.5, #top margin
                               margin.right = 1.5, #right margin
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
      stop("Plotting aborted: width must be single positive number (cm)")
    }
    if (width < 4) {
      message("Warning: width is very small (", width, " cm) – plot may be hard to read")
    }
    if (width > 50) {
      message("Warning: width is very large (", width, " cm) – plot may be unwieldy")
    }
    if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
      stop("Plotting aborted: height must be single positive number (cm)")
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
      stop("Plotting aborted: resolution must be single number ≥ 72 (dpi)")
    }
    if (resolution > 1200) {
      message("Warning: resolution is very high (", resolution, " dpi) – file may be huge")
    }
  }
  
  # Validate specified margins (reasonable range: 0–10)
  margin.list <- c(margin.bottom, margin.left, margin.top, margin.right)
  margin.names <- c("margin.bottom", "margin.left", "margin.top", "margin.right")
  for (i in seq_along(margin.list)) {
    if (!is.numeric(margin.list[i]) || length(margin.list[i]) != 1 || is.na(margin.list[i])) {
      stop("Plotting aborted: ", margin.names[i], " must be single numeric value")
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
      stop("Plotting aborted: bar.border.lwd must be single positive number")
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
  cluster_order <- hclust(dist(ancestry_proportions), method = linkage.method)$order
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
    
    use.colors <- layer.colors[1:K][layer.order] #get colors for clusters in correct order
    
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
    
    # Optional: draw vertical lines *between* bars (samples)
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
                 mar = c(margin.bottom,
                         margin.left,
                         margin.top,
                         margin.right),
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
plot.Learning.SOM <- function(SOM.output, 
                              col.pal = viridis::turbo, #set color palette
                              save = F, #option to save plot
                              overwrite = T, #option to overwrite plot if it already exists (only if saving plot)
                              plot.type = "svg", #options: "svg", "png", "jpg" (only if saving plot)
                              file.name = NULL, #set plot file.name (if NULL, default plot file.name is used; only if saving plot)
                              width = 20, #plot width in cm (only if saving plot)
                              height = 15, #plot height in cm (only if saving plot)
                              resolution = 300, #plot resolution in dpi (only if saving plot)
                              margin.bottom = 5, #bottom margin
                              margin.left = 5, #left margin
                              margin.top = 3, #top margin
                              margin.right = 2.5, #right margin
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
      stop("Plotting aborted: width must be single positive number (cm)")
    }
    if (width < 4) {
      message("Warning: width is very small (", width, " cm) – plot may be hard to read")
    }
    if (width > 50) {
      message("Warning: width is very large (", width, " cm) – plot may be unwieldy")
    }
    if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
      stop("Plotting aborted: height must be single positive number (cm)")
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
      stop("Plotting aborted: resolution must be single number ≥ 72 (dpi)")
    }
    if (resolution > 1200) {
      message("Warning: resolution is very high (", resolution, " dpi) – file may be huge")
    }
  }
  
  # Validate specified margins (reasonable range: 0–10)
  margin.list <- c(margin.bottom, margin.left, margin.top, margin.right)
  margin.names <- c("margin.bottom", "margin.left", "margin.top", "margin.right")
  for (i in seq_along(margin.list)) {
    if (!is.numeric(margin.list[i]) || length(margin.list[i]) != 1 || is.na(margin.list[i])) {
      stop("Plotting aborted: ", margin.names[i], " must be single numeric value")
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
    stop("Plotting aborted: lines.thickness must be single positive numeric value")
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
    stop("Plotting aborted: legend.lines.thickness must be single positive numeric value")
  }
  
  # Validate specified x.axis.label and y.axis.label (must be character)
  if (!is.character(x.axis.label) || length(x.axis.label) != 1 || is.na(x.axis.label)) {
    stop("Plotting aborted: x.axis.label must be single character string")
  }
  if (!is.character(y.axis.label) || length(y.axis.label) != 1 || is.na(y.axis.label)) {
    stop("Plotting aborted: y.axis.label must be single character string")
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
    stop("Plotting aborted: matrix names not found in provided SOM.output")
  }
  
  # Check if matrix names match number of layers
  if (length(matrix_names) != length(SOM.output$learning_values_list)) {
    stop("Plotting aborted: number of matrix names does not match number of matrices")
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
            col = grDevices::adjustcolor(layer_colors[i], alpha.f = lines.alpha),
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
                            overwrite = T, #option to overwrite plot if it already exists (only if saving plot)
                            plot.type = "svg", #plot type options: "svg", "png", "jpg" (only if saving plot)
                            file.name = NULL, #set file.name for saving (if NULL, default plot file.name is used; only if saving plot)
                            width = 20, #plot width in cm (only if plot is saved)
                            height = 15, #plot width in cm (only if plot is saved)
                            resolution = 300, #plot resolution in dpi (only if plot is saved)
                            margin.bottom = 5, #bottom margin
                            margin.left = 5, #left margin
                            margin.top = 3.5, #top margin
                            margin.right = 2.5, #right margin
                            title = "Layer weights", #plot title name
                            y_axis_label = "Relative weights" #set y axis label
) {
  
  # Reset plotting parameters
  old_dev <- dev.cur()
  old_plotting_parameters <- par(no.readonly = TRUE)
  on.exit({
    if (dev.cur() == old_dev) par(old_plotting_parameters)
  }, add = TRUE)
  
  # Validate SOM.output
  if (!"distance_matrix" %in% names(SOM.output)) {
    stop("Plotting aborted: distance_matrix could not be found in SOM output")
  }
  if (length(SOM.output$input_data_names) == 1) {
    stop("Plotting aborted: single-layer SOM detected (multiple layers are needed to plot and compare layer weights)")
  }
  d.mat <- SOM.output$distance_matrix #extract distance matrix
  layer.names <- SOM.output$input_data_names #extract layer names
  if (length(layer.names) != ncol(d.mat)) {
    message("Mismatch between layer names and distance matrix columns - using generic layer names")
    layer.names <- paste0("Layer ", seq_len(ncol(d.mat)))
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
      stop("Plotting aborted: width must be single positive number (cm)")
    }
    if (width < 4) {
      message("Warning: width is very small (", width, " cm) – plot may be hard to read")
    }
    if (width > 50) {
      message("Warning: width is very large (", width, " cm) – plot may be unwieldy")
    }
    if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
      stop("Plotting aborted: height must be single positive number (cm)")
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
      stop("Plotting aborted: resolution must be single number ≥ 72 (dpi)")
    }
    if (resolution > 1200) {
      message("Warning: resolution is very high (", resolution, " dpi) – file may be huge")
    }
  }
  
  # Validate specified margins (reasonable range: 0–10)
  margin.list <- c(margin.bottom, margin.left, margin.top, margin.right)
  margin.names <- c("margin.bottom", "margin.left", "margin.top", "margin.right")
  for (i in seq_along(margin.list)) {
    if (!is.numeric(margin.list[i]) || length(margin.list[i]) != 1 || is.na(margin.list[i])) {
      stop("Plotting aborted: ", margin.names[i], " must be single numeric value")
    }
    if (margin.list[i] < 0) {
      stop("Plotting aborted: ", margin.names[i], " must be ≥ 0")
    }
    if (margin.list[i] > 10) {
      message("Warning: ", margin.names[i], " is large (", margin.list[i], ") – plot area may shrink")
    }
  }
  
  # Calculate relative layer weights
  raw.weights <- sqrt(1 / ifelse(colMeans(d.mat, na.rm = T) > 0, colMeans(d.mat, na.rm = T), NA)) #compute weights
  names(raw.weights) <- layer.names #assign names to weights
  sorted.weights <- sort(raw.weights[!is.na(raw.weights)], decreasing = T) #sort weights
  
  # Set default file name for saving based on layer.names
  if (is.null(file.name)) { 
    file.name <- paste0("SOM_layers_plot_", paste(layer.names, collapse = "_"), ".", plot.type)
  }
  
  # Check for file existence if overwrite is FALSE
  if (!overwrite && file.exists(file.name)) {
    stop(file.name, " already exists - skipping plot saving")
  }
  
  # Define color palette for layers
  layer.cols <- setNames(col.pal(length(layer.names)), layer.names) #assign colors
  
  # Save plot if requested
  if (save) {
    
    # Set plot plot.type
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
  
  # Set layout and margins for plotting
  par(mfrow = c(1, 1), 
      oma = c(0, 0, 0, 0),
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
                       margin.bottom = 4, #bottom margin
                       margin.left = 5.5, #left margin
                       margin.top = 2.5, #top margin
                       margin.right = 2.5, #right margin
                       N.axis.labels.BIC.plot = 4, #number of axis labels for BIC plot
                       N.axis.labels.deltaBIC.plot = 4, #number of axis labels for deltaBIC plot
                       round.axis.labels.BIC.plot = 0, #rounding axis labels of BIC plot to X digits
                       round.axis.labels.deltaBIC.plot = 0, #rounding axis labels of delta BIC plot to X digits
                       title = "Number of clusters (k)") #change plot title 
{
  
  # Reset plotting parameters
  old_dev <- dev.cur()
  old_plotting_parameters <- par(no.readonly = TRUE)
  on.exit({
    if (dev.cur() == old_dev) par(old_plotting_parameters)
  }, add = TRUE)
  
  # Validate SOM.output
  if (is.null(SOM.output$BIC_values)) {
    stop("Plotting aborted: BIC_values could not be found in SOM output - check if clustering.SOM was run")
  } else if (is.null(SOM.output$N_replicates)) {
    stop("Plotting aborted: N_replicates could not be found in SOM output - check if clustering.SOM was run")
  } else if (is.null(SOM.output$optim_k_vals)) {
    stop("Plotting aborted: optim_k_vals could not be found in SOM output - check if clustering.SOM was run")
  } else if (is.null(SOM.output$max_k)) {
    stop("Plotting aborted: max_k could not be found in SOM output - check if clustering.SOM was run")
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
      stop("Plotting aborted: width must be single positive number (cm)")
    }
    if (width < 4) {
      message("Warning: width is very small (", width, " cm) – plot may be hard to read")
    }
    if (width > 50) {
      message("Warning: width is very large (", width, " cm) – plot may be unwieldy")
    }
    if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
      stop("Plotting aborted: height must be single positive number (cm)")
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
      stop("Plotting aborted: resolution must be single number ≥ 72 (dpi)")
    }
    if (resolution > 1200) {
      message("Warning: resolution is very high (", resolution, " dpi) – file may be huge")
    }
  }
  
  # Validate margins (reasonable range: 0–10)
  margin.list <- c(margin.bottom, margin.left, margin.top, margin.right)
  margin.names <- c("margin.bottom", "margin.left", "margin.top", "margin.right")
  for (i in seq_along(margin.list)) {
    if (!is.numeric(margin.list[i]) || length(margin.list[i]) != 1 || is.na(margin.list[i])) {
      stop("Plotting aborted: ", margin.names[i], " must be single numeric value")
    }
    if (margin.list[i] < 0) {
      stop("Plotting aborted: ", margin.names[i], " must be ≥ 0")
    }
    if (margin.list[i] > 10) {
      message("Warning: ", margin.names[i], " is large (", margin.list[i], ") – plot area may shrink")
    }
  }
  
  # Validate specified N.axis.labels.BIC.plot (must be positive integer)
  if (!is.numeric(N.axis.labels.BIC.plot) || length(N.axis.labels.BIC.plot) != 1 ||
      is.na(N.axis.labels.BIC.plot) || N.axis.labels.BIC.plot <= 0 || N.axis.labels.BIC.plot %% 1 != 0) {
    stop("Plotting aborted: N.axis.labels.BIC.plot must be single positive integer")
  }
  
  # Validate specified N.axis.labels.deltaBIC.plot (must be positive integer)
  if (!is.numeric(N.axis.labels.deltaBIC.plot) || length(N.axis.labels.deltaBIC.plot) != 1 ||
      is.na(N.axis.labels.deltaBIC.plot) || N.axis.labels.deltaBIC.plot <= 0 || N.axis.labels.deltaBIC.plot %% 1 != 0) {
    stop("Plotting aborted: N.axis.labels.deltaBIC.plot must be single positive integer")
  }
  
  # Validate round.axis.labels.BIC.plot (must be non-negative integer)
  if (!is.numeric(round.axis.labels.BIC.plot) || length(round.axis.labels.BIC.plot) != 1 ||
      is.na(round.axis.labels.BIC.plot) || round.axis.labels.BIC.plot < 0 || round.axis.labels.BIC.plot %% 1 != 0) {
    stop("Plotting aborted: round.axis.labels.BIC.plot must be single non-negative integer")
  }
  
  # Validate specified round.axis.labels.deltaBIC.plot (must be non-negative integer)
  if (!is.numeric(round.axis.labels.deltaBIC.plot) || length(round.axis.labels.deltaBIC.plot) != 1 ||
      is.na(round.axis.labels.deltaBIC.plot) || round.axis.labels.deltaBIC.plot < 0 || round.axis.labels.deltaBIC.plot %% 1 != 0) {
    stop("Plotting aborted: round.axis.labels.deltaBIC.plot must be single non-negative integer")
  }
  
  # Validate specified title (must be NULL or character string)
  if (!is.null(title) && (!is.character(title) || length(title) != 1 || is.na(title))) {
    stop("Plotting aborted: title must be NULL or single character string")
  }
  
  # Set color palette
  if (SOM.output$max_k == 2) { #for max_k = 2
    k.cols <- col.pal(SOM.output$max_k) 
  } else {
    k.cols <- col.pal(SOM.output$max_k - 1) 
  }
  
  # Open device if needed
  if (save) {
    if (is.null(file.name)) {
      file.name <- paste0("SOM_K_evaluation_plot_", paste(SOM.output$input_data_names, collapse = "_"), ".", plot.type)
    }
    
    if (!overwrite && file.exists(file.name)) {
      stop(file.name, " already exists - skipping plot saving")
    }
    
    if (plot.type == "svg") {
      svg(file.name, width = width / 2.54, height = height / 2.54)
    } else if (plot.type == "png") {
      png(file.name, width = width, height = height, res = resolution, units = "cm")
    } else if (plot.type == "jpg") {
      jpeg(file.name, width = width, height = height, res = resolution, units = "cm")
    }
  }
  
  # Plots for max_k = 2
  if (SOM.output$max_k == 2) {
    
    # Set plotting parameters
    par(mfrow = c(2, 1), 
        bty = "n",
        mar = c(margin.bottom, margin.left, margin.top, margin.right)) #set plot margins
    
    # Plot 1: Boxplot of BIC values
    BIC_vals <- unlist(SOM.output$BIC_values)
    bic_min <- floor(min(BIC_vals, na.rm = TRUE))
    bic_max <- ceiling(max(BIC_vals, na.rm = TRUE))
    bic_breaks <- round(seq(from = bic_min, to = bic_max, length.out = N.axis.labels.BIC.plot), 
                        round.axis.labels.BIC.plot)
    
    boxplot(t(SOM.output$BIC_values)[, 1:(SOM.output$max_k)],
            outline = FALSE,
            notch = FALSE,
            axes = FALSE,
            ylab = "BIC",
            ylim = range(bic_breaks),
            col = k.cols)
    
    axis(1, at = 1:(SOM.output$max_k), labels = NA)
    axis(2, at = bic_breaks, las = 3)
    title(title, line = 0)
    
    # Plot 2: Barplot of sampling frequency (how many times was each k value selected as optimal over replicates)
    barplot_data <- table(factor(SOM.output$optim_k_vals, levels = 1:(SOM.output$max_k))) / SOM.output$N_replicates
    bar_positions <- barplot(barplot_data, 
                             ylab = "Frequency of selected K", 
                             ylim = c(0, 1), 
                             col = k.cols,
                             axes = FALSE)
    axis(2, las = 3)
    
    # Plots for max_k > 2
  } else { 
    
    # Set plotting parameters
    par(mfrow = c(3, 1), 
        bty = "n",
        mar = c(margin.bottom, margin.left, margin.top, margin.right)) #set plot margins
    
    # Plot 1: Boxplot of BIC values
    BIC_vals <- unlist(SOM.output$BIC_values)
    bic_min <- floor(min(BIC_vals, na.rm = TRUE))
    bic_max <- ceiling(max(BIC_vals, na.rm = TRUE))
    bic_breaks <- round(seq(from = bic_min, to = bic_max, length.out = N.axis.labels.BIC.plot), 
                        round.axis.labels.BIC.plot)
    
    boxplot(t(SOM.output$BIC_values)[, 1:(SOM.output$max_k - 1)],
            outline = FALSE,
            notch = FALSE,
            axes = FALSE,
            ylab = "BIC",
            ylim = range(bic_breaks),
            col = k.cols)
    
    axis(1, at = 1:(SOM.output$max_k - 1), labels = NA)
    axis(2, at = bic_breaks, las = 3)
    title(title, line = 0)
    
    # Plot 2: Boxplot of delta BIC
    d_wss_raw <- apply(SOM.output$BIC_values, 2, function(x) diff(diff(x)))
    
    if (is.null(dim(d_wss_raw))) {
      # when max_k == 3, apply() returns a vector → coerce to 1×R matrix
      d_wss <- matrix(d_wss_raw,
                      nrow = 1,
                      dimnames = list(2:(SOM.output$max_k - 1),
                                      colnames(SOM.output$BIC_values)))
    } else {
      d_wss <- d_wss_raw
      rownames(d_wss) <- 2:(SOM.output$max_k - 1)
    }
    
    # prepend an NA row so that diff levels line up
    d_wss <- rbind(NA, d_wss)
    
    delta_vals <- unlist(na.omit(d_wss))
    delta_min <- floor(min(delta_vals, na.rm = TRUE))
    delta_max <- ceiling(max(delta_vals, na.rm = TRUE))
    delta_breaks <- round(seq(from = delta_min, to = delta_max, length.out = N.axis.labels.deltaBIC.plot), 
                          round.axis.labels.deltaBIC.plot)
    
    plot(1, 
         type = "n", 
         xlim = c(0.5, SOM.output$max_k - 0.5),
         ylim = range(delta_breaks), 
         xaxt = "n", 
         yaxt = "n",
         xlab = "", 
         ylab = "delta BIC")
    abline(h = 0, lty = 2, col = "black")
    boxplot(t(d_wss),
            outline = FALSE,
            notch = FALSE,
            axes = FALSE,
            ylab = "delta BIC",
            add = T,
            col = k.cols)
    axis(1, at = 1:(SOM.output$max_k - 1), labels = NA)
    axis(2, at = delta_breaks, las = 3)
    
    # Plot 3: Barplot of sampling frequency (how many times was each k value selected as optimal over replicates)
    barplot_data <- table(factor(SOM.output$optim_k_vals, levels = 1:(SOM.output$max_k - 1))) / SOM.output$N_replicates
    bar_positions <- barplot(barplot_data, 
                             ylab = "Frequency of selected K", 
                             ylim = c(0, 1), 
                             col = k.cols,
                             axes = FALSE)
    axis(2, las = 3)
  }
  
  # Close device if saving
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}


## Function to plot model results as SOM grids (showing sample assignment to cells, cell distances and boundaries between cell clusters)
plot.Model.SOM <- function(SOM.output,
                           col.pal.neighbor.dist = viridis::cividis, #color palette of neighbor distance plot (top)
                           col.pal.clusters = viridis::viridis, #color palette of cluster plot (bottom)
                           save = F, #option to save plot
                           overwrite = T, #option to overwrite plot if already present (only if saving plot)
                           plot.type = "svg", #plot type options: "svg", "png", "jpg" (only if saving plot)
                           file.name = NULL, #set plot file name (if NULL, default name is used; only if saving plot)
                           width = 10, #plot width in cm (only if saving plot)
                           height = 15, #plot height in cm (only if saving plot)
                           resolution = 300, #plot resolution in dpi (only if saving plot)
                           margin.bottom = 0, #bottom margin
                           margin.left = 0, #left margin
                           margin.top = 1, #top margin
                           margin.right = 0, #right margin
                           boundary.col.clusters = "red", #color of cluster boundaries (in bottom plot)
                           boundary.lwd.clusters = 3, #line width of cluster boundaries (in bottom plot)
                           point.col.clusters = "white", #color of sample points (in bottom plot)
                           point.shape.clusters = 19, #shape of sample points (in bottom plot)
                           point.size.clusters = 0.8, #size of sample points (in bottom plot)
                           cluster.shape.clusters = "straight", #shape ("straight" or "round") of cluster cells (in bottom plot)
                           cluster.shape.neighbor.dist = "straight", #shape ("straight" or "round") of cluster cells (in top plot)
                           shift.plot.clusters = 0.099, #shift bottom plot slightly to the right to align with gridraster of top plot
                           title.clusters = "SOM clusters", #title of cluster plot (bottom)
                           title.neighbor.dist = "SOM neighbor distances" #title of neighbor distances plot
) {
  
  # Reset plotting parameters
  old_dev <- dev.cur()
  old_plotting_parameters <- par(no.readonly = TRUE)
  on.exit({
    if (dev.cur() == old_dev) par(old_plotting_parameters)
  }, add = TRUE)
  
  # Validate SOM.output
  if (is.null(SOM.output$som_models) || is.null(SOM.output$som_clusters)) {
    stop("Plotting aborted: SOM.output is missing 'som_models' or 'som_clusters' - check SOM.output or rerun run.SOM/clustering.SOM")
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
      stop("Plotting aborted: width must be single positive number (cm)")
    }
    if (width < 4) {
      message("Warning: width is very small (", width, " cm) – plot may be hard to read")
    }
    if (width > 50) {
      message("Warning: width is very large (", width, " cm) – plot may be unwieldy")
    }
    if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
      stop("Plotting aborted: height must be single positive number (cm)")
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
      stop("Plotting aborted: resolution must be single number ≥ 72 (dpi)")
    }
    if (resolution > 1200) {
      message("Warning: resolution is very high (", resolution, " dpi) – file may be huge")
    }
  }
  
  # Validate specified margins (reasonable range: 0–10)
  margin.list <- c(margin.bottom, margin.left, margin.top, margin.right)
  margin.names <- c("margin.bottom", "margin.left", "margin.top", "margin.right")
  for (i in seq_along(margin.list)) {
    if (!is.numeric(margin.list[i]) || length(margin.list[i]) != 1 || is.na(margin.list[i])) {
      stop("Plotting aborted: ", margin.names[i], " must be single numeric value")
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
    stop("Plotting aborted: boundary.lwd.clusters must be single positive number")
  }
  
  # Validate specified point.size.clusters
  if (!is.numeric(point.size.clusters) || length(point.size.clusters) != 1 || is.na(point.size.clusters) || point.size.clusters <= 0) {
    stop("Plotting aborted: point.size.clusters must be single positive number")
  }
  
  # Validate specified point.shape.clusters (should be integer)
  if (!is.numeric(point.shape.clusters) || length(point.shape.clusters) != 1 || is.na(point.shape.clusters) || (point.shape.clusters %% 1 != 0)) {
    stop("Plotting aborted: point.shape.clusters must be single integer")
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
    stop("Plotting aborted: boundary.col.clusters must be single character (color name or hex)")
  }
  if (!is.character(point.col.clusters) || length(point.col.clusters) != 1 || is.na(point.col.clusters)) {
    stop("Plotting aborted: point.col.clusters must be single character (color name or hex)")
  }
  
  # Validate specified title.clusters and title.neighbor.dist (NULL or character)
  if (!is.null(title.clusters) && (!is.character(title.clusters) || length(title.clusters) != 1 || is.na(title.clusters))) {
    stop("Plotting aborted: title.clusters must be NULL or single character string")
  }
  if (!is.null(title.neighbor.dist) && (!is.character(title.neighbor.dist) || length(title.neighbor.dist) != 1 || is.na(title.neighbor.dist))) {
    stop("Plotting aborted: title.neighbor.dist must be NULL or single character string")
  }
  
  # Extract cluster numbers 
  som_model <- SOM.output$som_models[[1]]
  som_cluster <- SOM.output$som_clusters[[1]]
  
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
plot.Map.SOM <- function(SOM.output,
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
  old_dev <- dev.cur()
  old_plotting_parameters <- par(no.readonly = TRUE)
  on.exit({
    if (dev.cur() == old_dev) par(old_plotting_parameters)
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
      stop("Plotting aborted: width must be single positive number (cm)")
    }
    if (width < 4) {
      message("Warning: width is very small (", width, " cm) – plot may be hard to read")
    }
    if (width > 50) {
      message("Warning: width is very large (", width, " cm) – plot may be unwieldy")
    }
    if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
      stop("Plotting aborted: height must be single positive number (cm)")
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
      stop("Plotting aborted: resolution must be single number ≥ 72 (dpi)")
    }
    if (resolution > 1200) {
      message("Warning: resolution is very high (", resolution, " dpi) – file may be huge")
    }
  }
  
  # Validate specified lat.buffer.range
  if (!is.numeric(lat.buffer.range) || length(lat.buffer.range) != 1 || is.na(lat.buffer.range) || lat.buffer.range < 0) {
    stop("Plotting aborted: lat.buffer.range must be single non-negative number")
  }
  
  # Validate specified lon.buffer.range
  if (!is.numeric(lon.buffer.range) || length(lon.buffer.range) != 1 || is.na(lon.buffer.range) || lon.buffer.range < 0) {
    stop("Plotting aborted: lon.buffer.range must be single non-negative number")
  }
  
  # Validate specified pie.size
  if (!is.numeric(pie.size) || length(pie.size) != 1 || is.na(pie.size) || pie.size <= 0) {
    stop("Plotting aborted: pie.size must be single positive number")
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
    stop("Plotting aborted: USA.state.lwd must be single positive number")
  }
  
  # Validate specified USA.county.lwd
  if (!is.numeric(USA.county.lwd) || length(USA.county.lwd) != 1 || is.na(USA.county.lwd) || USA.county.lwd <= 0) {
    stop("Plotting aborted: USA.county.lwd must be single positive number")
  }
  # Validate specified north.arrow.position
  if (!is.numeric(north.arrow.position) || length(north.arrow.position) != 2 || any(is.na(north.arrow.position)) ||
      any(north.arrow.position < 0) || any(north.arrow.position > 1)) {
    stop("Plotting aborted: north.arrow.position must be numeric vector of length 2 with values in between 0 and 1")
  }
  
  # Validate specified north.arrow.length
  if (!is.numeric(north.arrow.length) || length(north.arrow.length) != 1 || is.na(north.arrow.length) || north.arrow.length <= 0) {
    stop("Plotting aborted: north.arrow.length must be single positive number")
  }
  
  # Validate specified north.arrow.lwd
  if (!is.numeric(north.arrow.lwd) || length(north.arrow.lwd) != 1 || is.na(north.arrow.lwd) || north.arrow.lwd <= 0) {
    stop("Plotting aborted: north.arrow.lwd must be single positive number")
  }
  
  # Validate specified north.arrow.N.position
  if (!is.numeric(north.arrow.N.position) || length(north.arrow.N.position) != 1 || is.na(north.arrow.N.position) || north.arrow.N.position < 0) {
    stop("Plotting aborted: north.arrow.N.position must be single non-negative number")
  }
  
  # Validate specified north.arrow.N.size
  if (!is.numeric(north.arrow.N.size) || length(north.arrow.N.size) != 1 || is.na(north.arrow.N.size) || north.arrow.N.size <= 0) {
    stop("Plotting aborted: north.arrow.N.size must be single positive number")
  }
  
  # Validate specified scale.position
  if (!is.numeric(scale.position) || length(scale.position) != 2 || any(is.na(scale.position)) ||
      any(scale.position < 0) || any(scale.position > 1)) {
    stop("Plotting aborted: scale.position must be numeric vector of length 2 with values in between 0 and 1")
  }
  
  # Validate specified scale.size
  if (!is.numeric(scale.size) || length(scale.size) != 1 || is.na(scale.size) || scale.size <= 0) {
    stop("Plotting aborted: scale.size must be single positive number")
  }
  
  # Validate specified scale.font.size
  if (!is.numeric(scale.font.size) || length(scale.font.size) != 1 || is.na(scale.font.size) || scale.font.size <= 0) {
    stop("Plotting aborted: scale.font.size must be single positive number")
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
    stop("Plotting aborted: legend.font.size must be single positive number")
  }
  
  # Validate specified legend.box
  if (!is.logical(legend.box) || length(legend.box) != 1 || is.na(legend.box)) {
    stop("Plotting aborted: legend.box must be TRUE or FALSE")
  }
  
  # Validate specified legend.symbol.size
  if (!is.numeric(legend.symbol.size) || length(legend.symbol.size) != 1 || is.na(legend.symbol.size) || legend.symbol.size <= 0) {
    stop("Plotting aborted: legend.symbol.size must be single positive number")
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
  if (dev.cur() > 1) {
    dev.off() #close current device if open to avoid unwanted graphic distortions and other effects
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
      svg(file.name, width = width / 2.54, height = height / 2.54)
    } else if (plot.type == "png") {
      png(file.name, width = width, height = height, res = resolution, units = "cm")
    } else if (plot.type == "jpg") {
      jpeg(file.name, width = width, height = height, res = resolution, units = "cm")
    } else {
      stop("Plotting aborted: unsupported plot file plot.type - choose from 'svg', 'png', or 'jpg'")
    }
  }
  
  # Set layout and margins
  par(mfrow = c(1, 1), 
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
  maps::map.scale(x = scale_position_x,
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


## Function to plot variable importance for each SOM layer (based on codebook vectors/neuron weights)
plot.Variable.Importance.SOM <- function(SOM.output, 
                                         col.pal = viridis::turbo, #color palette
                                         save = F, #option to save plot
                                         overwrite = T, #option to overwrite plot if it already exists
                                         plot.type = "svg", #plot file plot.type (options: "png", "svg", "jpg")
                                         file.name = NULL, #plot file name (if NULL, default file name is used)
                                         width = 20, #plot width in cm
                                         height = 15, #plot height in cm
                                         resolution = 300, #plot resolution in dpi
                                         bottom.margin.total = 3, #bottom margin of entire plot
                                         left.margin.total = 1, #left margin of entire plot
                                         top.margin.total = 0, #top margin of entire plot
                                         right.margin.total = 0, #right margin of entire plot
                                         bottom.margin = 0, #bottom margin of individual plots
                                         left.margin = 3.5, #left margin of individual plots
                                         top.margin = 4, #top margin of individual plots
                                         right.margin = 2, #right margin of individual plots
                                         bars.threshold.N = 50, #threshold for leaving out bar labels
                                         title.font.size = 1.2, #font size of title
                                         matrix.label.font.size = 1, #font size of matrix label(s)
                                         bar.label.font.size = 0.65, #font size of bar labels
                                         importance.threshold = 0.001 #threshold for showing variable importance
) {
  
  # Reset plotting parameters
  old_dev <- dev.cur()
  old_plotting_parameters <- par(no.readonly = TRUE)
  on.exit({
    if (dev.cur() == old_dev) par(old_plotting_parameters)
  }, add = TRUE)
  
  # Validate SOM.output
  if (is.null(SOM.output$codebook_vectors)) {
    stop("Plotting aborted: codebook_vectors not found in SOM.output")
  }
  if (is.null(SOM.output$input_data_names)) {
    stop("Plotting aborted: input_data_names not found in SOM.output - check SOM.output or rerun run.SOM function")
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
  
  # Validate specified width and height (reasonable values: 4–50 cm)
  if (save) {
    if (!is.numeric(width) || length(width) != 1 || is.na(width) || width <= 0) {
      stop("Plotting aborted: width must be single positive number (cm)")
    }
    if (width < 4) {
      message("Warning: width is very small (", width, " cm) – plot may be hard to read")
    }
    if (width > 50) {
      message("Warning: width is very large (", width, " cm) – plot may be unwieldy")
    }
    if (!is.numeric(height) || length(height) != 1 || is.na(height) || height <= 0) {
      stop("Plotting aborted: height must be single positive number (cm)")
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
      stop("Plotting aborted: resolution must be single number ≥ 72 (dpi)")
    }
    if (resolution > 1200) {
      message("Warning: resolution is very high (", resolution, " dpi) – file may be huge")
    }
  }
  
  # Validate specified margins (reasonable range: 0–10)
  margin.list <- c(bottom.margin.total, left.margin.total, top.margin.total, right.margin.total,
                   bottom.margin, left.margin, top.margin, right.margin)
  margin.names <- c("bottom.margin.total", "left.margin.total", "top.margin.total", "right.margin.total",
                    "bottom.margin", "left.margin", "top.margin", "right.margin")
  for (i in seq_along(margin.list)) {
    if (!is.numeric(margin.list[i]) || length(margin.list[i]) != 1 || is.na(margin.list[i])) {
      stop("Plotting aborted: ", margin.names[i], " must be single numeric value")
    }
    if (margin.list[i] < 0) {
      stop("Plotting aborted: ", margin.names[i], " must be ≥ 0")
    }
    if (margin.list[i] > 10) {
      message("Warning: ", margin.names[i], " is large (", margin.list[i], ") – plot area may shrink")
    }
  }
  
  # Validate specified bars.threshold.N
  if (!is.numeric(bars.threshold.N) || length(bars.threshold.N) != 1 || is.na(bars.threshold.N) ||
      bars.threshold.N < 0 || bars.threshold.N %% 1 != 0) {
    stop("Plotting aborted: bars.threshold.N must be single non-negative integer")
  }
  
  # Validate specified title.font.size
  if (!is.numeric(title.font.size) || length(title.font.size) != 1 || is.na(title.font.size) || title.font.size <= 0) {
    stop("Plotting aborted: title.font.size must be single positive number")
  }
  
  # Validate specified matrix.label.font.size
  if (!is.numeric(matrix.label.font.size) || length(matrix.label.font.size) != 1 || is.na(matrix.label.font.size) || matrix.label.font.size <= 0) {
    stop("Plotting aborted: matrix.label.font.size must be single positive number")
  }
  
  # Validate specified bar.label.font.size
  if (!is.numeric(bar.label.font.size) || length(bar.label.font.size) != 1 || is.na(bar.label.font.size) || bar.label.font.size <= 0) {
    stop("Plotting aborted: bar.label.font.size must be single positive number")
  }
  
  # Validate specified importance.threshold
  if (!is.numeric(importance.threshold) || length(importance.threshold) != 1 || is.na(importance.threshold) || importance.threshold < 0) {
    stop("Plotting aborted: importance.threshold must be single non-negative number")
  }
  
  # Extract all codebook vectors for all replicates for each layer from SOM.output
  all_layer_codes <- SOM.output$codebook_vectors
  
  # Extract matrix names from SOM.output
  matrix_names <- SOM.output$input_data_names
  
  # Set plot saving
  if (save) {
    
    # Set default plotting name
    if (is.null(file.name)) {
      file.name <- paste0("SOM_codebook_vectors_plot_", paste(matrix_names, collapse = "_"), ".", plot.type)
    }
    
    # Check overwrite
    if (file.exists(file.name) && !overwrite) {
      stop(paste(file.name, "already exists - set overwrite = TRUE to overwrite"))
    }
    
    # Set plot plot.type
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
    } else {
      stop("Plotting aborted: unsupported file plot.type - choose from 'svg', 'png', or 'jpg'")
    }
  }
  
  # Set plot margins
  par(mar = c(bottom.margin, 
              left.margin, 
              top.margin, 
              right.margin))
  par(oma = c(bottom.margin.total, 
              left.margin.total, 
              top.margin.total, 
              right.margin.total))
  
  # Set plot layout based on number of SOM layers
  num_layers <- length(all_layer_codes)
  if (num_layers == 1) {
    par(mfrow = c(1, 1))
  } else if (num_layers == 2) {
    par(mfrow = c(1, 2))
  } else if (num_layers == 3) {
    par(mfrow = c(1, 3))
  } else if (num_layers == 4) {
    par(mfrow = c(2, 2))
  } else if (num_layers == 5) {
    par(mfrow = c(2, 3))
  } else if (num_layers == 6) {
    par(mfrow = c(2, 3))
  } else if (num_layers == 7) {
    par(mfrow = c(2, 4))
  } else if (num_layers == 8) {
    par(mfrow = c(2, 4))
  } else if (num_layers == 9) {
    par(mfrow = c(3, 3))
  } else {
    num_rows <- ceiling(num_layers / 3)
    par(mfrow = c(num_rows, 3))
  }
  
  # Iterate over each matrix and generate plots
  for (i in seq_along(all_layer_codes)) {
    var_mat <- all_layer_codes[[i]] #rows = neurons × replicates, cols = variables
    keep_vars <- which(apply(var_mat, 2, function(x) mean(abs(x)) > importance.threshold)) #drop variables with all values < threshold:
    if (length(keep_vars) == 0) {
      message(paste("No variables exceed importance.threshold of", importance.threshold, "for", matrix_names[i], " - specify lower value for importance.threshold"))
      next
    }
    var_mat <- var_mat[, keep_vars, drop = FALSE]
    y_labels <- colnames(var_mat)
    num_bars <- ncol(var_mat)
    layer.col <- col.pal(length(all_layer_codes))[i]
    
    medians <- apply(var_mat, 2, median, na.rm = TRUE)
    sort_idx <- order(medians, decreasing = FALSE) #sort by median
    var_mat <- var_mat[, sort_idx, drop = FALSE]
    y_labels <- y_labels[sort_idx]
    
    boxplot(var_mat,
            horizontal = TRUE,
            las = 1,
            notch = FALSE,
            outline = FALSE,
            col = layer.col,
            axes = FALSE,        
            whisklty = if(num_bars > bars.threshold.N) 0 else 1, #whiskers
            staplelty = if(num_bars > bars.threshold.N) 0 else 1, #“cap” lines at ends
            names = rep("", num_bars))
    axis(1)
    text(x = par("usr")[1], 
         y = seq_len(num_bars),
         labels = if(num_bars > bars.threshold.N) rep("", num_bars) else y_labels,
         pos = 2, 
         offset = 0.6, 
         cex = bar.label.font.size, 
         xpd = TRUE)
    mtext(matrix_names[i], 
          side = 3, #top margin
          line = 0, #lines out from plot
          cex = matrix.label.font.size,
          font = 1)
  }
  
  # Reset layout for title
  layout(matrix(c(1), nrow = 1, ncol = 1))
  par(mar = c(0, 0, 3, 0))
  
  # Add main title
  if (num_layers > 1) {
    title(main = "Variable importance of SOM layers", 
          cex.main = title.font.size)
  } else {
    title(main = "Variable importance across replicates", 
          cex.main = title.font.size)
  }
  
  # Close graphics device
  if (save) {
    dev.off()
    message(paste("Plot", ifelse(overwrite, "overwritten to", "saved to"), file.name))
  }
}


## Function to remove highly correlated and low-variance variables
remove.lowCV.multicollinearity.SOM <- function(df, #data.frame with numeric columns (e.g., climatic, environmental or morphological variables)
                                               CV.threshold = 0.05, #numeric, remove variables with CV ≤ this value (only for non-binary vars)
                                               cor.threshold = 0.9, #numeric, remove variables correlated above this threshold (absolute)
                                               exclude.cols = NULL, #character vector of columns to exclude from filtering (e.g. Latitude, Longitude)
                                               verbose = TRUE #logical, print messages about filtering steps
) {
  # Validate input
  if (!is.numeric(CV.threshold) || length(CV.threshold) != 1 || is.na(CV.threshold) || CV.threshold < 0) stop("CV.threshold must be single non-negative numeric")
  if (!is.numeric(cor.threshold) || length(cor.threshold) != 1 || is.na(cor.threshold) || cor.threshold < 0 || cor.threshold > 1) stop("cor.threshold must be single numeric between 0 and 1")
  if (!is.null(exclude.cols)) {if (!all(exclude.cols %in% colnames(df))) stop("All exclude.cols must be column names in df")}
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) stop("verbose must be TRUE or FALSE")
  
  # Ensure dataframe format
  if (!is.data.frame(df)) df <- as.data.frame(df)
  
  # Preserve row names
  row_names <- rownames(df)
  
  # Remove columns to exclude from filtering
  if (!is.null(exclude.cols)) {
    df_filtered <- df[, !(colnames(df) %in% exclude.cols), drop = FALSE]
  } else {
    df_filtered <- df
  }
  
  # Ensure numeric
  df_filtered <- as.data.frame(lapply(df_filtered, as.numeric))
  
  # Detect binary variables (only 2 unique values 0 and 1, ignoring NA)
  is_binary <- sapply(df_filtered, function(x) {
    ux <- unique(na.omit(x))
    length(ux) == 2 && all(sort(ux) == c(0, 1))
  })
  
  # Calculate CV for non-binary variables only
  cv_values <- sapply(df_filtered[!is_binary], function(x) {
    mu <- mean(x, na.rm = TRUE)
    sigma <- sd(x, na.rm = TRUE)
    if (is.na(mu) || abs(mu) < .Machine$double.eps) {
      return(NA_real_)
    }
    sigma / abs(mu)
  })
  
  # Variables to keep by CV: all binary vars + non-binary vars with CV > threshold and not NA
  keep_vars_cv <- c(names(is_binary)[is_binary], names(cv_values)[!is.na(cv_values) & cv_values > CV.threshold])
  
  # Variables removed by CV filtering (non-binary vars with low CV or NA CV)
  removed_by_cv <- setdiff(colnames(df_filtered), keep_vars_cv)
  
  # Subset to variables kept by CV
  df_filtered <- df_filtered[, keep_vars_cv, drop = FALSE]
  
  # Calculate absolute correlation matrix, ignoring NAs pairwise
  cor_mat <- abs(cor(df_filtered, use = "pairwise.complete.obs"))
  diag(cor_mat) <- 0
  
  # Iteratively remove variables until no correlation > threshold remains
  remove_cols <- c()
  while (any(cor_mat > cor.threshold, na.rm = TRUE)) {
    max_idx <- which(cor_mat == max(cor_mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
    col_to_remove <- colnames(cor_mat)[max_idx[2]]
    remove_cols <- c(remove_cols, col_to_remove)
    cor_mat[, col_to_remove] <- 0
    cor_mat[col_to_remove, ] <- 0
  }
  
  # Remove correlated columns
  df_filtered <- df_filtered[, !(colnames(df_filtered) %in% remove_cols)]
  
  # Bind excluded columns back in original order
  if (!is.null(exclude.cols)) {
    df_result <- cbind(df[, exclude.cols, drop = FALSE], df_filtered)
  } else {
    df_result <- as.data.frame(df_filtered) # always return a data frame
  }
  
  # Restore row names
  rownames(df_result) <- row_names
  
  # Report filtering results if verbose
  if (verbose) {
    n_removed_cv <- length(removed_by_cv)
    n_retained_cv <- ncol(df) - length(exclude.cols) - n_removed_cv
    message(n_removed_cv, ifelse(n_removed_cv == 1, " variable", " variables"), " removed due to low CV ≤ ", CV.threshold, " (only for non-binary variables)")
    message(n_retained_cv, ifelse(n_retained_cv == 1, " variable", " variables"), " retained after CV filtering")
    
    n_removed_cor <- length(remove_cols)
    n_retained_cor <- ncol(df_filtered)
    message(n_removed_cor, ifelse(n_removed_cor == 1, " variable", " variables"), " removed due to high correlation > ", cor.threshold)
    message(n_retained_cor, ifelse(n_retained_cor == 1, " variable", " variables"), " retained after correlation filtering")
    
    total_retained <- ncol(df_result)
    if (!is.null(exclude.cols)) {
      message("Number of variables left after processing (including excluded columns): ", total_retained)
    } else {
      message("Number of variables left after processing: ", total_retained)
    } 
  }
  
  return(df_result)
}


# Function to import and process SNP data from VCF, genind, or NEXUS alignment
process.SNP.data.SOM <- function(vcf.path = NULL, #optional path to VCF file
                                 genind.input = NULL, #optional genind object
                                 nexus.path = NULL, #optional path to NEXUS alignment file
                                 missing.loci.cutoff.lenient = 0.5, #remove loci with > this proportion missing (1st lenient filter)
                                 missing.loci.cutoff.final = 0.2, #remove loci with > this proportion missing (2nd final stringent filter)
                                 missing.individuals.cutoff = 0.4, #remove individuals with > this proportion missing
                                 singleton.loci.filter = TRUE, #whether to remove singleton loci
                                 invariant.loci.filter = TRUE, #whether to remove invariant loci
                                 verbose = TRUE #whether to show filtering messages and final summary
) {
  # Quiet evaluation helper
  quiet <- function(expr) {
    temp_file <- tempfile() #create temp file
    con <- file(temp_file, open = "wt") #open connection
    sink(con) #redirect stdout
    sink(con, type = "message") #redirect messages
    on.exit({
      sink(type = "message") #restore message stream
      sink() #restore output stream
      close(con) #close connection
      unlink(temp_file) #remove temp file
    }, add = TRUE)
    result <- eval(expr) #evaluate expression
    return(result) #return result
  }
  
  # Validate which input provided
  nset <- sum(!is.null(vcf.path), !is.null(genind.input), !is.null(nexus.path)) #count provided inputs
  if (nset != 1) stop("Provide exactly one of vcf.path, genind.input, or nexus.path") #must provide only one
  if (!is.null(vcf.path) && !file.exists(vcf.path)) stop("VCF file does not exist: ", vcf.path) #check vcf path
  if (!is.null(nexus.path) && !file.exists(nexus.path)) stop("NEXUS file does not exist: ", nexus.path) #check nexus path
  
  # Process NEXUS file
  if (!is.null(nexus.path)) { #if nexus input
    # Validate NEXUS content
    seq_list <- ape::read.nexus.data(nexus.path) #read aligned sequences
    if (length(seq_list) == 0) stop("No sequences found in NEXUS file") #check for empty
    seq_lengths <- sapply(seq_list, length) #sequence lengths
    if (any(seq_lengths == 0)) stop("Empty sequences detected in NEXUS file") #empty sequence
    if (length(unique(seq_lengths)) != 1) stop("NEXUS file is not aligned: sequences have different lengths") #not aligned
    seq_mat <- t(sapply(seq_list, function(x) strsplit(paste0(x, collapse = ""), "")[[1]])) #make samples × sites character matrix
    if (any(dim(seq_mat) == 0)) stop("NEXUS alignment matrix is empty or malformed") #check for empty matrix
    alleles_per_site <- apply(seq_mat, 2, function(col) setdiff(unique(col), c("?", "N", "-"))) #biallelic columns only
    n_loci_start <- ncol(seq_mat) #total loci before filtering
    keep_cols <- which(sapply(alleles_per_site, length) == 2) #index of biallelic sites
    
    # Biallelic site filter
    if (verbose) message(n_loci_start - length(keep_cols), " of ", n_loci_start, " loci removed because they were not biallelic") #report non-biallelic filter
    if (length(keep_cols) == 0) stop("No biallelic SNPs found in alignment") #stop if none
    seq_bi <- seq_mat[, keep_cols, drop = FALSE] #keep only biallelic sites
    snp_mat <- as.data.frame(seq_bi, stringsAsFactors = FALSE) #to data frame
    
    # Recoding to 0/1
    for(j in seq_len(ncol(snp_mat))) { #recode to 0/1
      a <- sort(alleles_per_site[[keep_cols[j]]]) #sort alleles
      ref <- a[1] #reference allele
      alt <- a[2] #alternate allele
      snp_mat[[j]] <- ifelse(snp_mat[[j]] == ref, 0L, ifelse(snp_mat[[j]] == alt, 1L, NA_integer_)) #assign 0/1/NA
    }
    rownames(snp_mat) <- names(seq_list) #set rownames
    colnames(snp_mat) <- paste0("SNP", seq_len(ncol(snp_mat))) #set colnames
    
    # Invariant loci filter (initial)
    n_loci_bi <- ncol(snp_mat) #biallelic loci before removing invariant
    snp_mat <- snp_mat[, sapply(snp_mat, function(x) var(x, na.rm = TRUE) > 0), drop = FALSE] #remove invariant loci
    if (verbose) message(n_loci_bi - ncol(snp_mat), " of ", n_loci_bi, " loci removed because they were invariant") #report invariant filter
    
    # Loci missing (lenient)
    if (!is.null(missing.loci.cutoff.lenient)) {
      n_loci_before <- ncol(snp_mat)
      snp_mat <- snp_mat[, (colMeans(is.na(snp_mat)) <= missing.loci.cutoff.lenient), drop = FALSE]
      if (verbose) message(n_loci_before - ncol(snp_mat), " of ", n_loci_before, " loci removed due to >", missing.loci.cutoff.lenient * 100, "% missing data") #report
    }
    
    # Individuals missing
    if (!is.null(missing.individuals.cutoff)) {
      n_ind_before <- nrow(snp_mat)
      snp_mat <- snp_mat[(rowMeans(is.na(snp_mat)) <= missing.individuals.cutoff), , drop = FALSE]
      if (verbose) message(n_ind_before - nrow(snp_mat), " of ", n_ind_before, " individuals removed due to >", missing.individuals.cutoff * 100, "% missing data") #report
    }
    
    # Loci missing (strict)
    if (!is.null(missing.loci.cutoff.final)) {
      n_loci_before <- ncol(snp_mat)
      snp_mat <- snp_mat[, (colMeans(is.na(snp_mat)) <= missing.loci.cutoff.final), drop = FALSE]
      if (verbose) message(n_loci_before - ncol(snp_mat), " of ", n_loci_before, " loci removed due to >", missing.loci.cutoff.final * 100, "% missing data (stricter filter)") #report
    }
    
    # Invariant loci filter (final pass)
    n_loci_before <- ncol(snp_mat)
    snp_mat <- snp_mat[, sapply(snp_mat, function(x) var(x, na.rm = TRUE) > 0), drop = FALSE]
    if (verbose) message(n_loci_before - ncol(snp_mat), " of ", n_loci_before, " invariant loci removed (final pass)") #report
    
    # Final summary
    if (verbose) message("Final SNP matrix: ", nrow(snp_mat), " samples × ", ncol(snp_mat), " loci") #print summary
    return(snp_mat) #return binary SNP matrix (haploid/mtDNA style)
  }
  
  # Process genind object or VCF file
  genind.object <- NULL #initialize
  if (!is.null(genind.input)) {
    genind.object <- genind.input #use genind object if provided
  } else {
    vcf.object <- tryCatch({
      suppressWarnings({
        temp_con <- file(tempfile(), open = "wt") #open temp file
        sink(temp_con, type = "output") #suppress output
        on.exit({sink(type = "output"); close(temp_con)}, add = TRUE) #restore output on exit
        read.vcfR(vcf.path) #read VCF
      })
    }, error = function(e) stop("VCF could not be read: ", conditionMessage(e))) #catch error
    genind.object <- vcfR2genind(vcf.object) #convert to genind
  }
  
  # Biallelic loci filter
  total.loci.before <- nLoc(genind.object) #loci before
  biallelic.index <- which(genind.object$loc.n.all == 2) #biallelic indices
  genind.object <- genind.object[loc = biallelic.index, drop = TRUE] #keep biallelic
  total.loci.after <- nLoc(genind.object) #loci after
  if (verbose) message(total.loci.before - total.loci.after, " of ", total.loci.before, " loci removed because they were not biallelic") #report
  if (total.loci.after == 0) stop("No biallelic loci found") #stop if none
  
  # Loci missing (lenient)
  if (!is.null(missing.loci.cutoff.lenient)) {
    number.loci.before <- nLoc(genind.object) #loci before filter
    genind.object <- suppressMessages(poppr::missingno(genind.object, type = "loci", cutoff = missing.loci.cutoff.lenient, quiet = TRUE)) #filter
    number.loci.after <- nLoc(genind.object) #after filter
    if (verbose) message(number.loci.before - number.loci.after, " of ", number.loci.before, " loci removed due to >", missing.loci.cutoff.lenient * 100, "% missing data") #report
    if (number.loci.after == 0) stop("All loci removed after lenient missing data filter") #stop if all gone
  }
  
  # Individuals missing
  if (!is.null(missing.individuals.cutoff)) {
    number.ind.before <- nInd(genind.object) #individuals before
    genind.object <- suppressMessages(poppr::missingno(genind.object, type = "geno", cutoff = missing.individuals.cutoff)) #filter
    number.ind.after <- nInd(genind.object) #after filter
    if (verbose) message(number.ind.before - number.ind.after, " of ", number.ind.before, " individuals removed due to >", missing.individuals.cutoff * 100, "% missing data") #report
    if (number.ind.after == 0) stop("All individuals removed after missing data filter") #stop if all gone
  }
  
  # Singleton loci filter
  if (isTRUE(singleton.loci.filter)) {
    singleton.loci <- which(minorAllele(genind.object) < 1 / nInd(genind.object)) #singleton loci
    if (length(singleton.loci) > 0) {
      if (verbose) message(length(singleton.loci), " of ", nLoc(genind.object), " singleton loci removed") #report
      genind.object <- genind.object[loc = -singleton.loci, drop = TRUE] #remove them
      if (nLoc(genind.object) == 0) stop("All loci removed after singleton filter") #stop if all gone
    }
  }
  
  # Invariant loci filter (first pass)
  snp.matrix <- quiet(quote(makefreq(genind.object, missing = NA))) #allele frequencies
  if (isTRUE(invariant.loci.filter)) {
    invariant.loci <- which(apply(snp.matrix, 2, var, na.rm = TRUE) == 0) #invariant loci
    if (length(invariant.loci) > 0) {
      if (verbose) message(length(invariant.loci), " of ", ncol(snp.matrix), " invariant loci removed (1st pass)") #report
      genind.object <- genind.object[loc = -invariant.loci, drop = TRUE] #remove them
      if (nLoc(genind.object) == 0) stop("All loci removed after invariant filter (1st pass)") #stop if all gone
    }
  }
  
  # Loci missing (strict)
  if (!is.null(missing.loci.cutoff.final)) {
    number.loci.before <- nLoc(genind.object) #loci before
    genind.object <- suppressMessages(poppr::missingno(genind.object, type = "loci", cutoff = missing.loci.cutoff.final, quiet = TRUE)) #filter
    number.loci.after <- nLoc(genind.object) #after
    if (verbose) message(number.loci.before - number.loci.after, " of ", number.loci.before, " loci removed due to >", missing.loci.cutoff.final * 100, "% missing data (stricter filter)") #report
    if (number.loci.after == 0) stop("All loci removed after final missing data filter") #stop if all gone
  }
  
  # Invariant loci filter (final pass)
  snp.matrix <- quiet(quote(makefreq(genind.object, missing = NA))) #allele frequencies again
  if (isTRUE(invariant.loci.filter)) {
    final.invariant.loci <- which(apply(snp.matrix, 2, var, na.rm = TRUE) == 0) #final invariant check
    if (length(final.invariant.loci) > 0) {
      if (verbose) message(length(final.invariant.loci), " of ", ncol(snp.matrix), " invariant loci removed (final pass)") #report
      genind.object <- genind.object[loc = -final.invariant.loci, drop = TRUE] #remove them
      if (nLoc(genind.object) == 0) stop("All loci removed after final invariant filter") #stop if all gone
    }
  }
  
  # Convert to dosage matrix (0, 1, 2)
  snp.dosage <- adegenet::genind2df(genind.object, sep = "/", usepop = FALSE) #genotype strings
  snp.dosage <- snp.dosage[, -1, drop = FALSE] #remove ind column
  snp.dosage.mat <- apply(snp.dosage, 2, function(x) {
    sapply(x, function(gt) {
      if (is.na(gt)) return(NA) #NA if missing
      gt_split <- unlist(strsplit(gt, "/")) #split genotype
      sum(as.numeric(gt_split)) #sum alleles (0/1/2)
    })
  })
  snp.dosage.mat <- as.data.frame(snp.dosage.mat) #to data frame
  rownames(snp.dosage.mat) <- indNames(genind.object) #set rownames
  if (verbose) message("Final SNP matrix: ", nrow(snp.dosage.mat), " samples × ", ncol(snp.dosage.mat), " loci") #summary
  return(snp.dosage.mat) #return SNP matrix
}


## Function to convert specified categorical columns into binary (0/1) indicators
make.cols.binary.SOM <- function(dataframe, #dataframe - input data frame
                                 make.binary.cols, #make.binary.cols - character vector of categorical column names to convert
                                 remove.original.cols = TRUE, #remove.original.cols - if TRUE, remove original categorical columns after processing
                                 append.to.original = FALSE #append.to.original - if TRUE, append to input; if FALSE, return only binary indicators
) {
  if (!is.data.frame(dataframe)) stop("dataframe must be a data frame") #ensure input is a data frame
  if (!is.character(make.binary.cols)) stop("make.binary.cols must be character vector of column names") #check type
  if (length(make.binary.cols) == 0) stop("make.binary.cols must contain at least one column name") #non-empty
  if (!all(make.binary.cols %in% colnames(dataframe))) { #check all exist in df
    missing_cols <- make.binary.cols[!make.binary.cols %in% colnames(dataframe)] #identify missing
    stop("The following columns are not in data frame: ", paste(missing_cols, collapse = ", ")) #error if any missing
  }
  if (!is.logical(remove.original.cols) || length(remove.original.cols) != 1) stop("remove.original.cols must be TRUE or FALSE") #check logical
  if (!is.logical(append.to.original) || length(append.to.original) != 1) stop("append.to.original must be TRUE or FALSE") #check logical
  
  dataframe_subset <- dataframe[, make.binary.cols, drop = FALSE] #extract selected columns
  noncat_cols <- sapply(dataframe_subset, function(x) !is.factor(x) && !is.character(x)) #identify non-categorical
  if (any(noncat_cols)) {
    bad_cols <- names(noncat_cols[noncat_cols]) #get bad column names
    stop("The following columns are not categorical (factor or character): ", paste(bad_cols, collapse = ", ")) #stop if any bad
  }
  
  dataframe_subset <- as.data.frame(lapply(dataframe_subset, function(x) {if (!is.factor(x)) x <- as.factor(x); return(x)})) #convert to factor
  
  high_card <- sapply(dataframe_subset, nlevels) > 30 #check number of levels
  if (any(high_card)) warning("The following columns have >30 levels: ", paste(names(dataframe_subset)[high_card], collapse = ", ")) #warn if too many levels
  
  binary_list <- list() #initialize list of binary matrices
  for (colname in colnames(dataframe_subset)) {
    col_factor <- dataframe_subset[[colname]] #get factor column
    col_factor <- factor(col_factor, exclude = NULL) #preserve NA rows
    levs <- levels(col_factor) #get levels
    if (length(levs) < 2) {message("Skipping column '", colname, "' because it has fewer than 2 levels"); next} #skip if too few levels
    model_mat <- model.matrix(~ col_factor - 1) #create dummy matrix
    colnames(model_mat) <- paste0(colname, "_", levs) #assign readable names
    binary_list[[colname]] <- model_mat #store in list
  }
  
  if (length(binary_list) == 0) stop("No binary columns could be created — all input columns had fewer than 2 levels") #stop if nothing created
  binary_dataframe <- as.data.frame(do.call(cbind, binary_list)) #combine to single data frame
  rownames(binary_dataframe) <- rownames(dataframe) #preserve rownames
  
  if (append.to.original) {
    overlap <- intersect(colnames(dataframe), colnames(binary_dataframe)) #check for collisions
    if (length(overlap) > 0) stop("Cannot append binary variables: the following column names already exist in your data frame: ", paste(overlap, collapse = ", ")) #stop if collision
    dataframe_out <- cbind(dataframe, binary_dataframe) #append to original
    if (remove.original.cols) dataframe_out <- dataframe_out[, !(colnames(dataframe_out) %in% make.binary.cols), drop = FALSE] #remove originals
  } else {
    dataframe_out <- binary_dataframe #return only binary columns
  }
  
  return(dataframe_out) #return result
}
