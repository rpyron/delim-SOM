#### Set environment and load packages #########################################

rm(list = ls()) #clear environment
gc()
setwd("C:/Users/danie/Desktop/Phd research/Manuscripts/SOM package")
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")



#### Define input, output, plot, and PCA settings #############################

## Define input and output directories
simulation.directory <- "Simulations/Simulation_set_2" #define simulation directory
vcf.directory <- file.path(simulation.directory, "fastsimcoal2_vcf_files") #define input VCF directory
results.root.directory <- file.path(simulation.directory, "fastsimcoal2_results") #define root output results directory
STRUCTURE.loglik.directory <- file.path(results.root.directory, "Structure_loglik") #define STRUCTURE log-likelihood directory


## Define analysis sets and output root
analysis.sample.n.per.population.symmetric_24 <- 24 #number of individuals per population for symmetric_24
analysis.sample.n.per.population.symmetric_16 <- 16 #number of individuals per population for symmetric_16
analysis.sample.n.per.population.symmetric_8 <- 8 #number of individuals per population for symmetric_8

analysis.sample.n.by.population.asymmetric_24 <- c(pop1 = 36, pop2 = 12) #total = 48; ratio = 3:1
analysis.sample.n.by.population.asymmetric_16 <- c(pop1 = 24, pop2 = 8) #total = 32; ratio = 3:1
analysis.sample.n.by.population.asymmetric_8 <- c(pop1 = 12, pop2 = 4) #total = 16; ratio = 3:1

analysis.sample.random.seed <- 1 #base seed for reproducible random individual sampling
dir.create(results.root.directory, recursive = TRUE, showWarnings = FALSE) #create root results directory if needed


## Set main parameters
override <- FALSE #rerun all VCFs once because the result-table structure was changed; set to FALSE after a complete successful run
max_k <- 2
expected.mig.tags <- c("0", "1e-6", "4e-6", "7e-6") #define expected migration-rate filename tags
clustering_method_SOM <- "kmeans+BICelbow"

deNovo.kmeans.BIC.thresh <- 4 #minimum BIC improvement required to choose K2
deNovo.kmeans.n.iter <- 10000 #maximum number of k-means iterations
deNovo.kmeans.n.start <- 1000 #number of random starts for de novo k-means/BIC
deNovo.kmeans.max.n.clust <- 2 #maximum number of clusters tested by de novo k-means/BIC
deNovo.kmeans.max.proportion.PCs <- 0.9 #maximum proportion of possible PCs retained for de novo k-means/BIC
deNovo.kmeans.center <- TRUE #center genotype matrix before PCA for de novo k-means/BIC
deNovo.kmeans.scale <- FALSE #scale genotype matrix before PCA for de novo k-means/BIC

sNMF.K.values <- 1:2 #candidate K values
sNMF.repetitions <- 50 #number of sNMF replicate runs per K
sNMF.ploidy <- 2 #diploid SNP data
sNMF.cross.entropy.thresh <- 0 #minimum K2 cross-entropy improvement required; 0 means any improvement
sNMF.seed <- 1 #random seed

expected.structure.replicates <- 10 #expected number of STRUCTURE replicates per VCF and K
STRUCTURE.mean.lnprob.delta.threshold <- 50 #choose K2 if mean lnprob_K2 - mean lnprob_K1 is greater than this

migration.vcf.mig.for.PCA <- 0 #migration rate to select for PCA
migration.vcf.file.number.for.PCA <- 20 #file number N among VCF files for selected migration rate

plot.k2.threshold.proportion <- 0.5 #K2 detection probability used for vertical threshold lines
plot.fitted.prediction.n <- 1000 #number of points used for fitted GLM prediction lines
plot.point.size <- 2 #point size for observed K2 proportions
plot.point.alpha <- 0.8 #point transparency for observed K2 proportions
plot.fitted.line.width <- 1.1 #line width for fitted binomial GLM curves
plot.threshold.line.width <- 1.1 #line width for 50% threshold lines
plot.threshold.line.type <- "dashed" #line type for 50% threshold lines
plot.base.size <- 9 #base ggplot font size
plot.axis.title.size <- 9 #axis title font size
plot.axis.text.size <- 7.1 #axis number/tick-label font size
plot.legend.title.size <- 9 #legend title font size
plot.legend.text.size <- 9 #legend label font size
plot.width.cm <- 16 #plot width in centimeters
plot.height.cm <- 10 #plot height in centimeters
plot.output.device <- "svg" #plot output device
plot.output.units <- "cm" #plot output units

Fst_plot_point_size <- 1.7 #Fst plot point size
Fst_plot_point_alpha <- 0.4 #Fst plot point transparency
Fst_plot_line_width <- 2 #Fst loess line width



#### Check VCF input files #####################################################

## List and check VCF files
if (!dir.exists(vcf.directory)) stop(paste("VCF directory does not exist:", vcf.directory))
vcf.file.paths <- list.files(vcf.directory, pattern = "\\.vcf$", full.names = TRUE) #list all VCF files
if (length(vcf.file.paths) == 0) stop(paste("No VCF files found in directory:", vcf.directory))
vcf.file.names <- basename(vcf.file.paths) #extract VCF file names
if (any(duplicated(vcf.file.names))) stop("Duplicate VCF file names found in VCF directory")
expected.vcf.file.pattern <- "^sim[0-9]+_tdiv[0-9.]+_mig.+\\.vcf$" #define expected VCF filename pattern
invalid.vcf.file.names <- vcf.file.names[!grepl(expected.vcf.file.pattern, vcf.file.names)] #find invalid VCF file names
if (length(invalid.vcf.file.names) > 0) {
  invalid.file.message <- paste(utils::head(invalid.vcf.file.names, 20), collapse = ", ")
  if (length(invalid.vcf.file.names) > 20) invalid.file.message <- paste0(invalid.file.message, ", ...")
  stop(paste("VCF filenames do not match expected pattern:", invalid.file.message))
}


## Parse VCF filename metadata
vcf.filename.matches <- regexec("^sim([0-9]+)_tdiv([0-9.]+)_mig(.+)\\.vcf$", vcf.file.names)
vcf.filename.parts <- regmatches(vcf.file.names, vcf.filename.matches)
vcf.metadata.table <- do.call(rbind, lapply(seq_along(vcf.filename.parts), function(vcf.index) {
  current.parts <- vcf.filename.parts[[vcf.index]]
  current.mig.tag <- current.parts[4]
  current.mig <- suppressWarnings(as.numeric(current.mig.tag))
  data.frame(file = vcf.file.names[vcf.index],
             full.path = vcf.file.paths[vcf.index],
             mig.index = as.integer(current.parts[2]),
             tdiv = as.numeric(current.parts[3]),
             mig.tag = current.mig.tag,
             mig = current.mig,
             stringsAsFactors = FALSE)
})) #create VCF metadata table
if (any(is.na(vcf.metadata.table$mig.index))) stop("At least one VCF has missing parsed migration-scenario index")
if (any(is.na(vcf.metadata.table$tdiv))) stop("At least one VCF has missing parsed divergence time")
if (any(is.na(vcf.metadata.table$mig))) stop("At least one VCF has missing parsed migration rate")


## Check duplicate divergence/migration combinations
duplicate.vcf.key <- paste(vcf.metadata.table$tdiv, vcf.metadata.table$mig.tag, sep = "_") #create unique divergence/migration key
if (any(duplicated(duplicate.vcf.key))) {
  duplicate.vcf.files <- vcf.metadata.table$file[duplicated(duplicate.vcf.key)]
  duplicate.file.message <- paste(utils::head(duplicate.vcf.files, 20), collapse = ", ")
  if (length(duplicate.vcf.files) > 20) duplicate.file.message <- paste0(duplicate.file.message, ", ...")
  stop(paste("Duplicate divergence/migration VCF combinations found:", duplicate.file.message))
}


## Check file size and VCF header
vcf.file.info <- file.info(vcf.metadata.table$full.path) #get VCF file information
empty.vcf.files <- vcf.metadata.table$file[is.na(vcf.file.info$size) | vcf.file.info$size <= 0]
if (length(empty.vcf.files) > 0) {
  empty.file.message <- paste(utils::head(empty.vcf.files, 20), collapse = ", ")
  if (length(empty.vcf.files) > 20) empty.file.message <- paste0(empty.file.message, ", ...")
  stop(paste("Empty or unreadable VCF files found:", empty.file.message))
}
vcf.first.lines <- vapply(vcf.metadata.table$full.path, function(current.vcf.path) {
  current.first.line <- tryCatch(readLines(current.vcf.path, n = 1, warn = FALSE), error = function(e) NA_character_)
  if (length(current.first.line) == 0) return(NA_character_)
  return(current.first.line[1])
}, character(1)) #read first line from each VCF
invalid.vcf.header.files <- vcf.metadata.table$file[is.na(vcf.first.lines) | !grepl("^##fileformat=VCF", vcf.first.lines)]
if (length(invalid.vcf.header.files) > 0) {
  invalid.header.message <- paste(utils::head(invalid.vcf.header.files, 20), collapse = ", ")
  if (length(invalid.vcf.header.files) > 20) invalid.header.message <- paste0(invalid.header.message, ", ...")
  stop(paste("Files do not start with a valid VCF header:", invalid.header.message))
}


## Check expected migration combinations
missing.mig.tags <- setdiff(expected.mig.tags, unique(vcf.metadata.table$mig.tag)) #find missing migration rates
if (length(missing.mig.tags) > 0) stop(paste("Missing expected migration-rate VCF combinations:", paste(missing.mig.tags, collapse = ", ")))
unexpected.mig.tags <- setdiff(unique(vcf.metadata.table$mig.tag), expected.mig.tags) #find unexpected migration-rate tags
if (length(unexpected.mig.tags) > 0) stop(paste("Unexpected migration-rate VCF combinations:", paste(unexpected.mig.tags, collapse = ", ")))


## Check that VCFs are balanced across migration combinations by divergence time
mig.tdiv.list <- split(vcf.metadata.table$tdiv, vcf.metadata.table$mig.tag) #split divergence times by migration rate
reference.tdiv.set <- sort(unique(vcf.metadata.table$tdiv)) #define reference divergence-time set
unbalanced.mig <- names(mig.tdiv.list)[vapply(mig.tdiv.list, function(current.tdiv.values) {!identical(sort(unique(current.tdiv.values)), reference.tdiv.set)}, logical(1))] #find migration groups with missing or extra divergence times
if (length(unbalanced.mig) > 0) stop(paste("VCF set is unbalanced across migration combinations:", paste(unbalanced.mig, collapse = ", ")))


## Order checked VCF metadata
vcf.metadata.table <- vcf.metadata.table[order(vcf.metadata.table$mig.index, vcf.metadata.table$tdiv, vcf.metadata.table$mig), ] #order VCF metadata table
rownames(vcf.metadata.table) <- NULL


## Summarize VCF input files and simulated divergence-time distribution
vcf.mig.file.count.table <- as.data.frame(table(vcf.metadata.table$mig.tag), stringsAsFactors = FALSE) #summarize file counts by migration rate
colnames(vcf.mig.file.count.table) <- c("mig.tag", "n.files")
simulated.tdiv.table <- data.frame(tdiv.index = seq_along(sort(unique(vcf.metadata.table$tdiv))), tdiv = sort(unique(vcf.metadata.table$tdiv)), stringsAsFactors = FALSE) #extract unique simulated divergence times
simulated.tdiv.values <- simulated.tdiv.table$tdiv #extract divergence-time values
simulated.tdiv.summary.table <- data.frame(n.simulated.tdiv = length(simulated.tdiv.values),
                                           n.unique.tdiv = length(unique(simulated.tdiv.values)),
                                           min.tdiv = min(simulated.tdiv.values),
                                           q025.tdiv = unname(stats::quantile(simulated.tdiv.values, 0.025, names = FALSE)),
                                           q250.tdiv = unname(stats::quantile(simulated.tdiv.values, 0.250, names = FALSE)),
                                           mean.tdiv = mean(simulated.tdiv.values),
                                           median.tdiv = stats::median(simulated.tdiv.values),
                                           q750.tdiv = unname(stats::quantile(simulated.tdiv.values, 0.750, names = FALSE)),
                                           q975.tdiv = unname(stats::quantile(simulated.tdiv.values, 0.975, names = FALSE)),
                                           max.tdiv = max(simulated.tdiv.values),
                                           stringsAsFactors = FALSE) #summarize simulated divergence-time distribution


## Check VCF files
cat("Total VCF files:", nrow(vcf.metadata.table), "\n")
cat("Unique simulated divergence times:", nrow(simulated.tdiv.table), "\n\n")
vcf.mig.file.count.table
simulated.tdiv.summary.table
head(simulated.tdiv.table)



#### Create additional functions ###############################################

## Create function to parse VCF filename metadata
parse.vcf.filename.metadata <- function(vcf.filename) {
  filename.matches <- regexec("^sim([0-9]+)_tdiv([0-9.]+)_mig(.+)\\.vcf$", vcf.filename)
  filename.parts <- regmatches(vcf.filename, filename.matches)[[1]]
  if (length(filename.parts) == 0) stop(paste("Filename does not match expected pattern:", vcf.filename))
  mig.tag <- filename.parts[4]
  migration.value <- suppressWarnings(as.numeric(mig.tag))
  if (is.na(migration.value)) stop(paste("Could not parse migration tag in filename:", vcf.filename))
  metadata.table <- data.frame(file = vcf.filename,
                               mig.index = as.integer(filename.parts[2]),
                               tdiv = as.numeric(filename.parts[3]),
                               mig.tag = mig.tag,
                               mig = migration.value,
                               stringsAsFactors = FALSE)
  return(metadata.table)
}


## Create function to create safe file tag for one VCF file
create.vcf.result.file.tag <- function(vcf.file.metadata) {
  file.tag <- tools::file_path_sans_ext(vcf.file.metadata$file)
  file.tag <- gsub("[^[:alnum:]_\\-]+", "_", file.tag)
  return(file.tag)
}


## Create function to create per-VCF intermediate file paths
create.vcf.intermediate.file.paths <- function(vcf.file.metadata, intermediate.results.directory) {
  file.tag <- create.vcf.result.file.tag(vcf.file.metadata = vcf.file.metadata)
  result.summary.directory <- file.path(intermediate.results.directory, "result_summary") #define result-summary subfolder
  optim.k.summary.directory <- file.path(intermediate.results.directory, "optim_k_summary") #define optim_k-summary subfolder
  fst.summary.directory <- file.path(intermediate.results.directory, "Weir_Cockerham_Fst") #define Fst-summary subfolder
  deNovo.kmeans.summary.directory <- file.path(intermediate.results.directory, "de_novo_kmeans_BIC") #define de novo k-means/BIC-summary subfolder
  sNMF.summary.directory <- file.path(intermediate.results.directory, "sNMF_summary") #define sNMF-summary subfolder
  sNMF.results.directory <- file.path(intermediate.results.directory, "sNMF_projects", file.tag) #define sNMF project subfolder
  dir.create(result.summary.directory, recursive = TRUE, showWarnings = FALSE) #create result-summary subfolder if needed
  dir.create(optim.k.summary.directory, recursive = TRUE, showWarnings = FALSE) #create optim_k-summary subfolder if needed
  dir.create(fst.summary.directory, recursive = TRUE, showWarnings = FALSE) #create Fst-summary subfolder if needed
  dir.create(deNovo.kmeans.summary.directory, recursive = TRUE, showWarnings = FALSE) #create de novo k-means/BIC-summary subfolder if needed
  dir.create(sNMF.summary.directory, recursive = TRUE, showWarnings = FALSE) #create sNMF-summary subfolder if needed
  dir.create(sNMF.results.directory, recursive = TRUE, showWarnings = FALSE) #create sNMF project subfolder if needed
  result.summary.path <- file.path(result.summary.directory, paste0(file.tag, "_result_summary.csv"))
  optim.k.summary.path <- file.path(optim.k.summary.directory, paste0(file.tag, "_optim_k_summary.csv"))
  fst.summary.path <- file.path(fst.summary.directory, paste0(file.tag, "_Weir_Cockerham_Fst.csv"))
  deNovo.kmeans.summary.path <- file.path(deNovo.kmeans.summary.directory, paste0(file.tag, "_de_novo_kmeans_BIC.csv"))
  sNMF.summary.path <- file.path(sNMF.summary.directory, paste0(file.tag, "_sNMF_summary.csv"))
  return(list(result.summary.path = result.summary.path,
              optim.k.summary.path = optim.k.summary.path,
              fst.summary.path = fst.summary.path,
              deNovo.kmeans.summary.path = deNovo.kmeans.summary.path,
              sNMF.summary.path = sNMF.summary.path,
              sNMF.results.directory = sNMF.results.directory))
}


## Create function to check whether per-VCF intermediate results already exist
check.vcf.intermediate.results.available <- function(vcf.file.metadata, intermediate.results.directory) {
  intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = vcf.file.metadata, intermediate.results.directory = intermediate.results.directory)
  files.exist <- file.exists(intermediate.file.paths$result.summary.path) & file.exists(intermediate.file.paths$optim.k.summary.path)
  return(files.exist)
}


## Create function to save per-VCF intermediate results
save.vcf.intermediate.results <- function(result.summary.row,
                                          optim.k.summary.table,
                                          vcf.file.metadata,
                                          intermediate.results.directory) {
  intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = vcf.file.metadata, intermediate.results.directory = intermediate.results.directory)
  write.csv(result.summary.row, intermediate.file.paths$result.summary.path, row.names = FALSE) #save one-row result summary
  write.csv(optim.k.summary.table, intermediate.file.paths$optim.k.summary.path, row.names = FALSE) #save full optim_k_summary table
}


## Create function to load per-VCF intermediate results
load.vcf.intermediate.results <- function(vcf.file.metadata, intermediate.results.directory) {
  intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = vcf.file.metadata, intermediate.results.directory = intermediate.results.directory)
  result.summary.row <- read.csv(intermediate.file.paths$result.summary.path, stringsAsFactors = FALSE) #load one-row result summary
  optim.k.summary.table <- read.csv(intermediate.file.paths$optim.k.summary.path, stringsAsFactors = FALSE) #load full optim_k_summary table
  return(list(result.summary.row = result.summary.row, optim.k.summary.table = optim.k.summary.table))
}


## Create function to check whether per-VCF intermediate Fst results already exist
check.vcf.intermediate.Fst.results.available <- function(vcf.file.metadata, intermediate.results.directory) {
  intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = vcf.file.metadata, intermediate.results.directory = intermediate.results.directory)
  files.exist <- file.exists(intermediate.file.paths$fst.summary.path)
  return(files.exist)
}


## Create function to save per-VCF intermediate Fst results
save.vcf.intermediate.Fst.results <- function(fst.row, vcf.file.metadata, intermediate.results.directory) {
  intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = vcf.file.metadata, intermediate.results.directory = intermediate.results.directory)
  write.csv(fst.row, intermediate.file.paths$fst.summary.path, row.names = FALSE) #save one-row Fst result
}


## Create function to load per-VCF intermediate Fst results
load.vcf.intermediate.Fst.results <- function(vcf.file.metadata, intermediate.results.directory) {
  intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = vcf.file.metadata, intermediate.results.directory = intermediate.results.directory)
  fst.row <- read.csv(intermediate.file.paths$fst.summary.path, stringsAsFactors = FALSE) #load one-row Fst result
  return(fst.row)
}



## Create function to check whether per-VCF intermediate de novo k-means/BIC results already exist
check.vcf.intermediate.deNovo.kmeans.results.available <- function(vcf.file.metadata, intermediate.results.directory) {
  intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = vcf.file.metadata, intermediate.results.directory = intermediate.results.directory)
  files.exist <- file.exists(intermediate.file.paths$deNovo.kmeans.summary.path)
  return(files.exist)
}


## Create function to save per-VCF intermediate de novo k-means/BIC results
save.vcf.intermediate.deNovo.kmeans.results <- function(deNovo.kmeans.row, vcf.file.metadata, intermediate.results.directory) {
  intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = vcf.file.metadata, intermediate.results.directory = intermediate.results.directory)
  write.csv(deNovo.kmeans.row, intermediate.file.paths$deNovo.kmeans.summary.path, row.names = FALSE) #save one-row de novo k-means/BIC result
}


## Create function to load per-VCF intermediate de novo k-means/BIC results
load.vcf.intermediate.deNovo.kmeans.results <- function(vcf.file.metadata, intermediate.results.directory) {
  intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = vcf.file.metadata, intermediate.results.directory = intermediate.results.directory)
  deNovo.kmeans.row <- read.csv(intermediate.file.paths$deNovo.kmeans.summary.path, stringsAsFactors = FALSE) #load one-row de novo k-means/BIC result
  return(deNovo.kmeans.row)
}


## Create function to check whether per-VCF intermediate sNMF results already exist
check.vcf.intermediate.sNMF.results.available <- function(vcf.file.metadata, intermediate.results.directory) {
  intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = vcf.file.metadata, intermediate.results.directory = intermediate.results.directory)
  if (!file.exists(intermediate.file.paths$sNMF.summary.path)) return(FALSE)
  sNMF.row <- tryCatch(read.csv(intermediate.file.paths$sNMF.summary.path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(sNMF.row)) return(FALSE)
  results.available <- "sNMF.status" %in% colnames(sNMF.row) && length(sNMF.row$sNMF.status) > 0 && !is.na(sNMF.row$sNMF.status[1]) && sNMF.row$sNMF.status[1] == "ok"
  return(results.available)
}


## Create function to save per-VCF intermediate sNMF results
save.vcf.intermediate.sNMF.results <- function(sNMF.row, vcf.file.metadata, intermediate.results.directory) {
  intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = vcf.file.metadata, intermediate.results.directory = intermediate.results.directory)
  write.csv(sNMF.row, intermediate.file.paths$sNMF.summary.path, row.names = FALSE) #save one-row sNMF result
}


## Create function to load per-VCF intermediate sNMF results
load.vcf.intermediate.sNMF.results <- function(vcf.file.metadata, intermediate.results.directory) {
  intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = vcf.file.metadata, intermediate.results.directory = intermediate.results.directory)
  sNMF.row <- read.csv(intermediate.file.paths$sNMF.summary.path, stringsAsFactors = FALSE) #load one-row sNMF result
  return(sNMF.row)
}


## Create function to convert vcfR object to LEA geno format
convert.vcfR.to.LEA.geno <- function(vcf.object, geno.file.path) {
  genotype.matrix <- vcfR::extract.gt(vcf.object, element = "GT", as.numeric = FALSE) #extract GT genotypes
  if (is.null(genotype.matrix)) stop("Could not extract GT genotypes from VCF")
  if (nrow(genotype.matrix) < 1) stop("No loci available for LEA geno file")
  if (ncol(genotype.matrix) < 2) stop("Fewer than two individuals available for LEA geno file")
  genotype.code.matrix <- matrix("9", nrow = nrow(genotype.matrix), ncol = ncol(genotype.matrix)) #initialize LEA genotype matrix
  rownames(genotype.code.matrix) <- rownames(genotype.matrix)
  colnames(genotype.code.matrix) <- colnames(genotype.matrix)
  for (locus.index in seq_len(nrow(genotype.matrix))) {
    current.genotypes <- genotype.matrix[locus.index, ] #extract genotypes for one locus
    current.genotypes <- sub(":.*$", "", current.genotypes) #remove FORMAT fields after GT if present
    current.genotypes[is.na(current.genotypes)] <- "./." #convert NA to missing genotype
    current.codes <- vapply(strsplit(current.genotypes, "[/|]"), function(current.alleles) {
      if (length(current.alleles) != 2) return("9")
      if (any(is.na(current.alleles)) || any(current.alleles == ".")) return("9")
      if (!all(current.alleles %in% c("0", "1"))) return("9")
      return(as.character(sum(current.alleles == "1")))
    }, character(1)) #convert genotype to 0, 1, 2, or 9
    genotype.code.matrix[locus.index, ] <- current.codes #store LEA codes
  }
  writeLines(apply(genotype.code.matrix, 1, paste0, collapse = ""), con = geno.file.path) #write one SNP per row, one individual per character
  return(geno.file.path)
}



## Create function to assign fastsimcoal2 populations from sample names
assign.fastsimcoal2.populations <- function(genind.object) {
  population.vector <- ifelse(grepl("^G_1_", adegenet::indNames(genind.object)), "pop1", ifelse(grepl("^G_2_", adegenet::indNames(genind.object)), "pop2", NA_character_)) #assign populations from sample names
  if (any(is.na(population.vector))) stop("Some individuals do not match expected G_1_ or G_2_ sample-name pattern")
  genind.object <- adegenet::`pop<-`(genind.object, value = population.vector) #store populations in genind object
  return(genind.object)
}


## Create function to create deterministic sampling seed for one simulation/population
create.fastsimcoal2.subsample.seed <- function(sampling.group, population.name, sampling.seed = 1) {
  sampling.group.character <- as.character(sampling.group) #store simulation group as character
  sampling.group.integer <- suppressWarnings(as.integer(sampling.group.character)) #try to convert simulation group to integer
  if (is.na(sampling.group.integer)) sampling.group.integer <- sum(utf8ToInt(sampling.group.character)) #fallback for non-integer group labels
  population.offset <- match(population.name, c("pop1", "pop2")) #create population-specific offset
  if (is.na(population.offset)) population.offset <- sum(utf8ToInt(as.character(population.name))) #fallback for unexpected population names
  seed.value <- as.integer((as.integer(sampling.seed) + sampling.group.integer * 1000L + population.offset * 100000L) %% .Machine$integer.max) #combine base seed, simulation group, and population
  if (seed.value < 1L) seed.value <- 1L #avoid invalid seed
  return(seed.value)
}


## Create function to select n samples from one population
select.fastsimcoal2.samples.from.population <- function(population.samples, n.per.population = NULL, sampling.group = NULL, sampling.seed = 1, population.name = NA_character_) {
  if (is.null(n.per.population) || is.na(n.per.population)) return(population.samples) #keep all samples if no target sample size is given
  if (length(population.samples) < n.per.population) stop(paste("Population", population.name, "has fewer than", n.per.population, "samples"))
  if (length(population.samples) == n.per.population) return(population.samples) #keep original order if all samples are used
  old.random.seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL #store current random seed
  on.exit({
    if (is.null(old.random.seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv) #restore missing seed state
    } else {
      assign(".Random.seed", old.random.seed, envir = .GlobalEnv) #restore previous seed state
    }
  }, add = TRUE)
  if (!is.null(sampling.group) && !is.na(sampling.group)) set.seed(create.fastsimcoal2.subsample.seed(sampling.group = sampling.group, population.name = population.name, sampling.seed = sampling.seed)) #use same random subset for the same simulation group and population
  population.samples.for.sampling <- sort(population.samples) #sort samples so matched sampling does not depend on VCF sample order
  sampled.population.samples <- sample(population.samples.for.sampling, n.per.population, replace = FALSE) #randomly sample individuals
  sampled.population.samples <- population.samples[population.samples %in% sampled.population.samples] #return selected individuals in original VCF order
  return(sampled.population.samples)
}


## Create function to read sample names from VCF header
read.vcf.sample.names.from.header <- function(vcf.file.path) {
  vcf.connection <- file(vcf.file.path, open = "r")
  on.exit(close(vcf.connection), add = TRUE)
  chrom.header.line <- NULL
  repeat {
    current.line <- readLines(vcf.connection, n = 1, warn = FALSE)
    if (length(current.line) == 0) break
    if (grepl("^#CHROM", current.line)) {
      chrom.header.line <- current.line
      break
    }
  }
  if (is.null(chrom.header.line)) stop(paste("Could not find #CHROM header line in VCF:", vcf.file.path))
  header.fields <- strsplit(chrom.header.line, "\t", fixed = TRUE)[[1]]
  if (length(header.fields) <= 9) stop(paste("VCF contains no sample columns:", vcf.file.path))
  sample.names <- header.fields[-seq_len(9)]
  return(sample.names)
}


## Create function to convert sampling target to named pop1/pop2 vector
standardize.sample.n.by.population <- function(sample.n.by.population) {
  if (length(sample.n.by.population) == 1) sample.n.by.population <- c(pop1 = sample.n.by.population, pop2 = sample.n.by.population) #convert symmetric sample size to named pop1/pop2 vector
  if (is.null(names(sample.n.by.population)) || !all(c("pop1", "pop2") %in% names(sample.n.by.population))) stop("sample.n.by.population must be either a single number or a named vector with pop1 and pop2") #check population names
  sample.n.by.population <- sample.n.by.population[c("pop1", "pop2")] #order population sample sizes
  sample.n.by.population <- stats::setNames(as.integer(sample.n.by.population), c("pop1", "pop2")) #convert to integer while preserving names
  if (any(is.na(sample.n.by.population)) || any(sample.n.by.population < 1)) stop("sample.n.by.population contains invalid values") #check sample sizes
  return(sample.n.by.population)
}


## Create function to create sample-exclusion table for one VCF and one sampling design
create.samples.to.exclude.for.one.vcf <- function(vcf.file.metadata,
                                                  sample.n.by.population,
                                                  analysis.set.name,
                                                  sampling.seed = 1) {
  sample.n.by.population <- standardize.sample.n.by.population(sample.n.by.population = sample.n.by.population)
  sample.names <- read.vcf.sample.names.from.header(vcf.file.path = vcf.file.metadata$full.path)
  population.vector <- ifelse(grepl("^G_1_", sample.names), "pop1", ifelse(grepl("^G_2_", sample.names), "pop2", NA_character_))
  if (any(is.na(population.vector))) stop(paste("Some VCF samples do not match expected G_1_ or G_2_ pattern:", vcf.file.metadata$file))
  population.names <- sort(unique(population.vector))
  if (length(population.names) != 2) stop(paste("Expected exactly two populations in VCF:", vcf.file.metadata$file))
  selected.samples <- unlist(lapply(population.names, function(population.name) {
    population.samples <- sample.names[population.vector == population.name]
    current.sample.n <- sample.n.by.population[population.name]
    select.fastsimcoal2.samples.from.population(population.samples = population.samples,
                                                n.per.population = current.sample.n,
                                                sampling.group = vcf.file.metadata$simulation.group,
                                                sampling.seed = sampling.seed,
                                                population.name = population.name)
  }), use.names = FALSE)
  excluded.samples <- sample.names[!(sample.names %in% selected.samples)]
  if (length(excluded.samples) == 0) {
    excluded.sample.table <- data.frame(file = character(0),
                                        mig.index = integer(0),
                                        tdiv = numeric(0),
                                        mig.tag = character(0),
                                        mig = numeric(0),
                                        analysis.set.name = character(0),
                                        sample.n.pop1 = integer(0),
                                        sample.n.pop2 = integer(0),
                                        sample = character(0),
                                        population = character(0),
                                        action = character(0),
                                        stringsAsFactors = FALSE)
    return(excluded.sample.table)
  }
  excluded.sample.table <- data.frame(file = vcf.file.metadata$file,
                                      mig.index = vcf.file.metadata$mig.index,
                                      tdiv = vcf.file.metadata$tdiv,
                                      mig.tag = vcf.file.metadata$mig.tag,
                                      mig = vcf.file.metadata$mig,
                                      analysis.set.name = analysis.set.name,
                                      sample.n.pop1 = unname(sample.n.by.population["pop1"]),
                                      sample.n.pop2 = unname(sample.n.by.population["pop2"]),
                                      sample = excluded.samples,
                                      population = population.vector[match(excluded.samples, sample.names)],
                                      action = "exclude",
                                      stringsAsFactors = FALSE)
  return(excluded.sample.table)
}


## Create function to create sample-exclusion table for all VCFs and one sampling design
create.samples.to.exclude.table <- function(vcf.metadata.table,
                                            sample.n.by.population,
                                            analysis.set.name,
                                            sampling.seed = 1) {
  if (!("simulation.group" %in% colnames(vcf.metadata.table))) vcf.metadata.table$simulation.group <- match(vcf.metadata.table$tdiv, sort(unique(vcf.metadata.table$tdiv)))
  exclusion.list <- lapply(seq_len(nrow(vcf.metadata.table)), function(vcf.index) {
    create.samples.to.exclude.for.one.vcf(vcf.file.metadata = vcf.metadata.table[vcf.index, ],
                                          sample.n.by.population = sample.n.by.population,
                                          analysis.set.name = analysis.set.name,
                                          sampling.seed = sampling.seed)
  })
  exclusion.table <- do.call(rbind, exclusion.list)
  rownames(exclusion.table) <- NULL
  return(exclusion.table)
}


## Create function to subset genind object to n individuals per population
subset.genind.to.n.per.population <- function(genind.object, n.per.population = NULL, sampling.group = NULL, sampling.seed = 1) {
  if (is.null(n.per.population) || length(n.per.population) == 0 || all(is.na(n.per.population))) return(genind.object) #keep all samples if no target sample size is given
  sample.n.by.population <- standardize.sample.n.by.population(sample.n.by.population = n.per.population) #convert symmetric/asymmetric target to named pop1/pop2 vector
  population.vector <- as.character(adegenet::pop(genind.object)) #extract population labels
  population.names <- sort(unique(population.vector)) #extract population names
  if (length(population.names) != 2) stop("Expected exactly two populations before subsetting")
  selected.individuals <- unlist(lapply(population.names, function(population.name) {
    population.individuals <- adegenet::indNames(genind.object)[population.vector == population.name] #extract individuals in one population
    current.sample.n <- sample.n.by.population[population.name] #extract target sample size for current population
    select.fastsimcoal2.samples.from.population(population.samples = population.individuals,
                                                n.per.population = current.sample.n,
                                                sampling.group = sampling.group,
                                                sampling.seed = sampling.seed,
                                                population.name = population.name) #randomly select individuals
  }), use.names = FALSE) #combine selected individuals
  genind.object <- genind.object[selected.individuals, ] #subset genind object
  return(genind.object)
}


## Create function to read VCF as genind and optionally subset individuals per population
read.fastsimcoal2.genind.for.analysis <- function(vcf.file.path, n.per.population = NULL, sampling.group = NULL, sampling.seed = 1) {
  vcf.object <- vcfR::read.vcfR(vcf.file.path, verbose = FALSE) #read VCF
  genind.object.raw <- vcfR::vcfR2genind(vcf.object) #convert VCF to genind
  genind.object.raw <- assign.fastsimcoal2.populations(genind.object = genind.object.raw) #assign populations
  genind.object.raw <- subset.genind.to.n.per.population(genind.object = genind.object.raw,
                                                         n.per.population = n.per.population,
                                                         sampling.group = sampling.group,
                                                         sampling.seed = sampling.seed) #optionally subset individuals per population
  return(genind.object.raw)
}


## Create function to subset vcfR object to n individuals per population
subset.vcfR.to.n.per.population <- function(vcf.object, n.per.population = NULL, sampling.group = NULL, sampling.seed = 1) {
  if (is.null(n.per.population) || length(n.per.population) == 0 || all(is.na(n.per.population))) return(vcf.object) #keep all samples if no target sample size is given
  sample.n.by.population <- standardize.sample.n.by.population(sample.n.by.population = n.per.population) #convert symmetric/asymmetric target to named pop1/pop2 vector
  sample.names <- colnames(vcf.object@gt)[-1] #extract VCF sample names
  population.vector <- ifelse(grepl("^G_1_", sample.names), "pop1", ifelse(grepl("^G_2_", sample.names), "pop2", NA_character_)) #assign populations from sample names
  if (any(is.na(population.vector))) stop("Some VCF samples do not match expected G_1_ or G_2_ sample-name pattern")
  population.names <- sort(unique(population.vector)) #extract population names
  if (length(population.names) != 2) stop("Expected exactly two populations before VCF subsetting")
  selected.samples <- unlist(lapply(population.names, function(population.name) {
    population.samples <- sample.names[population.vector == population.name] #extract samples in one population
    current.sample.n <- sample.n.by.population[population.name] #extract target sample size for current population
    select.fastsimcoal2.samples.from.population(population.samples = population.samples,
                                                n.per.population = current.sample.n,
                                                sampling.group = sampling.group,
                                                sampling.seed = sampling.seed,
                                                population.name = population.name) #randomly select samples
  }), use.names = FALSE) #combine selected samples
  vcf.object@gt <- vcf.object@gt[, c(colnames(vcf.object@gt)[1], selected.samples), drop = FALSE] #subset genotype table
  return(vcf.object)
}


## Create function to create one-row sNMF summary
create.sNMF.summary.row <- function(vcf.file.metadata,
                                    sNMF.best.k = NA_integer_,
                                    sNMF.cross.entropy.K1 = NA_real_,
                                    sNMF.cross.entropy.K2 = NA_real_,
                                    sNMF.cross.entropy.drop.K1.to.K2 = NA_real_,
                                    sNMF.best.run.K1 = NA_integer_,
                                    sNMF.best.run.K2 = NA_integer_,
                                    sNMF.n.samples.raw = NA_integer_,
                                    sNMF.n.loci.raw = NA_integer_,
                                    sNMF.n.loci.used = NA_integer_,
                                    sNMF.repetitions = NA_integer_,
                                    sNMF.status = "ok",
                                    sNMF.error = NA_character_) {
  sNMF.row <- data.frame(file = vcf.file.metadata$file,
                         sNMF.best.k = sNMF.best.k,
                         sNMF.cross.entropy.K1 = sNMF.cross.entropy.K1,
                         sNMF.cross.entropy.K2 = sNMF.cross.entropy.K2,
                         sNMF.cross.entropy.drop.K1.to.K2 = sNMF.cross.entropy.drop.K1.to.K2,
                         sNMF.best.run.K1 = sNMF.best.run.K1,
                         sNMF.best.run.K2 = sNMF.best.run.K2,
                         sNMF.n.samples.raw = sNMF.n.samples.raw,
                         sNMF.n.loci.raw = sNMF.n.loci.raw,
                         sNMF.n.loci.used = sNMF.n.loci.used,
                         sNMF.repetitions = sNMF.repetitions,
                         sNMF.status = sNMF.status,
                         sNMF.error = sNMF.error,
                         stringsAsFactors = FALSE) #create one-row sNMF result
  return(sNMF.row)
}


## Create function to run LEA/sNMF for one VCF file
run.sNMF.for.one.vcf <- function(vcf.file.path,
                                 vcf.file.metadata,
                                 sNMF.results.directory,
                                 K.values = 1:2,
                                 repetitions = 100,
                                 CPU = 1,
                                 ploidy = 2,
                                 cross.entropy.thresh = 0,
                                 seed = 1,
                                 sample.n.per.population = NULL,
                                 sampling.group = NULL,
                                 sampling.seed = 1) {
  if (!all(c(1, 2) %in% K.values)) stop("K.values must include K1 and K2")
  if (!requireNamespace("LEA", quietly = TRUE)) stop("Package LEA is required")
  dir.create(sNMF.results.directory, recursive = TRUE, showWarnings = FALSE) #create sNMF results directory if needed
  vcf.object <- vcfR::read.vcfR(vcf.file.path, verbose = FALSE) #read VCF
  vcf.object <- subset.vcfR.to.n.per.population(vcf.object = vcf.object,
                                                n.per.population = sample.n.per.population,
                                                sampling.group = sampling.group,
                                                sampling.seed = sampling.seed) #optionally subset individuals per population
  n.samples.raw <- ncol(vcf.object@gt) - 1 #store number of samples
  n.loci.raw <- nrow(vcf.object@fix) #store number of raw loci
  if (n.loci.raw < 1) stop("VCF contains no loci")
  if (n.samples.raw < 2) stop("VCF contains fewer than two samples")
  n.loci.used <- n.loci.raw #use all simulated loci
  sNMF.working.directory <- tempfile(pattern = "sNMF_") #create short temporary sNMF working directory
  dir.create(sNMF.working.directory, recursive = TRUE, showWarnings = FALSE) #create short temporary sNMF working directory
  old.working.directory <- getwd() #store current working directory
  on.exit({
    setwd(old.working.directory) #restore working directory
    unlink(sNMF.working.directory, recursive = TRUE, force = TRUE) #delete temporary sNMF working directory
  }, add = TRUE)
  sNMF.geno.file.path <- file.path(sNMF.working.directory, "input.geno") #define short geno path
  sNMF.geno.file.path <- convert.vcfR.to.LEA.geno(vcf.object = vcf.object, geno.file.path = sNMF.geno.file.path) #convert full VCF to LEA geno format
  setwd(sNMF.working.directory)
  set.seed(seed)
  suppressWarnings(suppressMessages(utils::capture.output({
    sNMF.project <- LEA::snmf(input.file = "input.geno", K = K.values, repetitions = repetitions, CPU = CPU, entropy = TRUE, ploidy = ploidy, project = "new", seed = seed)
  }, file = nullfile()))) #run unsupervised sNMF silently
  sNMF.cross.entropy.table <- do.call(rbind, lapply(K.values, function(current.K) {
    current.cross.entropy <- tryCatch({
      suppressWarnings(suppressMessages(utils::capture.output({current.cross.entropy.tmp <- LEA::cross.entropy(sNMF.project, K = current.K)}, file = nullfile())))
      current.cross.entropy.tmp
    }, error = function(e) rep(NA_real_, repetitions)) #extract cross-entropy values
    data.frame(K = current.K,
               run = seq_along(as.numeric(current.cross.entropy)),
               cross.entropy = as.numeric(current.cross.entropy),
               stringsAsFactors = FALSE)
  })) #combine cross-entropy values
  finite.cross.entropy.table <- sNMF.cross.entropy.table[is.finite(sNMF.cross.entropy.table$cross.entropy) & !is.na(sNMF.cross.entropy.table$cross.entropy), ] #keep finite values
  if (nrow(finite.cross.entropy.table) == 0) stop("No finite sNMF cross-entropy values returned")
  sNMF.best.cross.entropy.by.K <- do.call(rbind, lapply(split(finite.cross.entropy.table,
                                                              finite.cross.entropy.table$K), function(K.table) {
                                                                best.run.index <- which.min(K.table$cross.entropy)
                                                                data.frame(K = K.table$K[best.run.index], best.run = K.table$run[best.run.index], best.cross.entropy = K.table$cross.entropy[best.run.index], stringsAsFactors = FALSE)
                                                              })) #select best run per K
  rownames(sNMF.best.cross.entropy.by.K) <- NULL
  sNMF.cross.entropy.K1 <- sNMF.best.cross.entropy.by.K$best.cross.entropy[match(1, sNMF.best.cross.entropy.by.K$K)] #extract K1 cross-entropy
  sNMF.cross.entropy.K2 <- sNMF.best.cross.entropy.by.K$best.cross.entropy[match(2, sNMF.best.cross.entropy.by.K$K)] #extract K2 cross-entropy
  sNMF.best.run.K1 <- sNMF.best.cross.entropy.by.K$best.run[match(1, sNMF.best.cross.entropy.by.K$K)] #extract best K1 run
  sNMF.best.run.K2 <- sNMF.best.cross.entropy.by.K$best.run[match(2, sNMF.best.cross.entropy.by.K$K)] #extract best K2 run
  if (is.na(sNMF.cross.entropy.K1) || is.na(sNMF.cross.entropy.K2)) stop("K1 or K2 cross-entropy is missing")
  sNMF.cross.entropy.drop.K1.to.K2 <- sNMF.cross.entropy.K1 - sNMF.cross.entropy.K2 #positive means K2 improves cross-entropy
  if (sNMF.cross.entropy.drop.K1.to.K2 <= cross.entropy.thresh) {
    sNMF.best.k <- 1L #select K1 if K2 does not improve enough
  } else {
    sNMF.best.k <- 2L #select K2 if K2 improves enough
  }
  sNMF.row <- create.sNMF.summary.row(vcf.file.metadata = vcf.file.metadata,
                                      sNMF.best.k = sNMF.best.k,
                                      sNMF.cross.entropy.K1 = sNMF.cross.entropy.K1,
                                      sNMF.cross.entropy.K2 = sNMF.cross.entropy.K2,
                                      sNMF.cross.entropy.drop.K1.to.K2 = sNMF.cross.entropy.drop.K1.to.K2,
                                      sNMF.best.run.K1 = sNMF.best.run.K1,
                                      sNMF.best.run.K2 = sNMF.best.run.K2,
                                      sNMF.n.samples.raw = n.samples.raw,
                                      sNMF.n.loci.raw = n.loci.raw,
                                      sNMF.n.loci.used = n.loci.used,
                                      sNMF.repetitions = repetitions,
                                      sNMF.status = "ok",
                                      sNMF.error = NA_character_) #create one-row sNMF result
  return(list(sNMF.row = sNMF.row,
              sNMF.cross.entropy.table = sNMF.cross.entropy.table,
              sNMF.best.cross.entropy.by.K = sNMF.best.cross.entropy.by.K,
              sNMF.project = sNMF.project))
}


## Create function to calculate and append sNMF results
calculate.and.append.sNMF <- function(vcf.metadata.table,
                                      result.table,
                                      optim.k.result.table,
                                      intermediate.results.directory,
                                      override = FALSE,
                                      K.values = 1:2,
                                      repetitions = 100,
                                      CPU = 1,
                                      ploidy = 2,
                                      cross.entropy.thresh = 0,
                                      seed = 1,
                                      sample.n.per.population = NULL,
                                      sampling.seed = 1) {
  sNMF.list <- vector("list", nrow(vcf.metadata.table)) #initialize sNMF list
  for (vcf.index in seq_len(nrow(vcf.metadata.table))) {
    current.metadata <- vcf.metadata.table[vcf.index, ] #extract current metadata
    if (!override && check.vcf.intermediate.sNMF.results.available(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory)) {
      cat("\nLoading saved sNMF result for", vcf.index, "of", nrow(vcf.metadata.table), ":", current.metadata$file, "\n")
      sNMF.list[[vcf.index]] <- load.vcf.intermediate.sNMF.results(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory) #load saved per-VCF sNMF result
    } else {
      cat("\nRunning sNMF for", vcf.index, "of", nrow(vcf.metadata.table), ":", current.metadata$file, "\n")
      sNMF.list[[vcf.index]] <- tryCatch({
        intermediate.file.paths <- create.vcf.intermediate.file.paths(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory)
        unlink(intermediate.file.paths$sNMF.results.directory, recursive = TRUE, force = TRUE) #delete old sNMF project directory before rerun
        dir.create(intermediate.file.paths$sNMF.results.directory, recursive = TRUE, showWarnings = FALSE) #recreate sNMF project directory
        sNMF.run <- run.sNMF.for.one.vcf(vcf.file.path = current.metadata$full.path,
                                         vcf.file.metadata = current.metadata,
                                         sNMF.results.directory = intermediate.file.paths$sNMF.results.directory,
                                         K.values = K.values,
                                         repetitions = repetitions,
                                         CPU = CPU,
                                         ploidy = ploidy,
                                         cross.entropy.thresh = cross.entropy.thresh,
                                         seed = seed,
                                         sample.n.per.population = sample.n.per.population,
                                         sampling.group = current.metadata$simulation.group,
                                         sampling.seed = sampling.seed)
        sNMF.run$sNMF.row
      }, error = function(e) {
        create.sNMF.summary.row(vcf.file.metadata = current.metadata,
                                sNMF.repetitions = repetitions,
                                sNMF.status = "error",
                                sNMF.error = conditionMessage(e))
      })
      save.vcf.intermediate.sNMF.results(sNMF.row = sNMF.list[[vcf.index]], vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory) #save per-VCF sNMF result immediately after calculation
    }
    stop.if.status.not.ok(status.value = sNMF.list[[vcf.index]]$sNMF.status,
                          error.value = sNMF.list[[vcf.index]]$sNMF.error,
                          file.value = current.metadata$file,
                          method.label = "sNMF") #stop immediately if sNMF failed
  }
  sNMF.table <- do.call(rbind, sNMF.list) #combine sNMF rows
  rownames(sNMF.table) <- NULL
  sNMF.columns <- setdiff(colnames(sNMF.table), "file")
  result.table <- result.table[, !(colnames(result.table) %in% sNMF.columns), drop = FALSE] #remove existing sNMF columns if present
  optim.k.result.table <- optim.k.result.table[, !(colnames(optim.k.result.table) %in% sNMF.columns), drop = FALSE] #remove existing sNMF columns if present
  for (sNMF.column in sNMF.columns) {
    result.table[[sNMF.column]] <- sNMF.table[[sNMF.column]][match(result.table$file, sNMF.table$file)] #add sNMF column to main results
    optim.k.result.table[[sNMF.column]] <- sNMF.table[[sNMF.column]][match(optim.k.result.table$file, sNMF.table$file)] #add sNMF column to optim_k results
  }
  return(list(result.table = result.table, optim.k.result.table = optim.k.result.table, sNMF.table = sNMF.table))
}


## Create function to coerce optim_k_summary to data frame
coerce.optim.k.summary.to.data.frame <- function(optim.k.summary) {
  if (is.null(optim.k.summary)) return(data.frame())
  if (is.data.frame(optim.k.summary)) {
    optim.k.summary.table <- optim.k.summary
  } else if (is.matrix(optim.k.summary)) {
    optim.k.summary.table <- as.data.frame(optim.k.summary, stringsAsFactors = FALSE)
  } else {
    optim.k.summary.table <- tryCatch(as.data.frame(optim.k.summary, stringsAsFactors = FALSE), error = function(e) data.frame())
  }
  if (nrow(optim.k.summary.table) > 0) {
    optim.k.summary.table$k.label <- rownames(optim.k.summary.table) #store original row names such as k2, k3
    optim.k.summary.table$optim.k.row <- seq_len(nrow(optim.k.summary.table)) #store row index
  }
  rownames(optim.k.summary.table) <- NULL
  return(optim.k.summary.table)
}



## Create function to extract row-level k proportions from optim_k_summary
extract.k.proportions.from.summary <- function(optim.k.summary.table, max.k) {
  k.proportion.vector <- rep(0, max.k)
  names(k.proportion.vector) <- paste0("proportion.k", seq_len(max.k))
  if (is.null(optim.k.summary.table) || nrow(optim.k.summary.table) == 0) return(k.proportion.vector)
  if (!("k.label" %in% colnames(optim.k.summary.table)) || !("Proportion" %in% colnames(optim.k.summary.table))) return(k.proportion.vector)
  valid.rows <- !is.na(optim.k.summary.table$k.label) & grepl("^k[0-9]+$", optim.k.summary.table$k.label)
  if (!any(valid.rows)) return(k.proportion.vector)
  current.k.values <- as.integer(sub("^k", "", optim.k.summary.table$k.label[valid.rows]))
  current.proportions <- as.numeric(optim.k.summary.table$Proportion[valid.rows])
  for (current.index in seq_along(current.k.values)) {
    current.k <- current.k.values[current.index]
    if (!is.na(current.k) && current.k >= 1 && current.k <= max.k) k.proportion.vector[paste0("proportion.k", current.k)] <- current.proportions[current.index]
  }
  return(k.proportion.vector)
}



## Create function to stop when intermediate result contains error
stop.if.status.not.ok <- function(status.value, error.value, file.value, method.label, allowed.status = "ok") {
  if (length(status.value) == 0 || is.na(status.value[1]) || !(as.character(status.value[1]) %in% allowed.status)) {
    error.message <- if (length(error.value) == 0 || is.na(error.value[1]) || !nzchar(as.character(error.value[1]))) {
      "No error message stored"
    } else {
      as.character(error.value[1])
    }
    stop(paste(method.label, "failed for", file.value, ":", error.message))
  }
  invisible(TRUE)
}


## Create function to find STRUCTURE per-replicate log-likelihood file
find.STRUCTURE.replicate.file <- function(analysis.set.name, structure.loglik.directory) {
  expected.file.name <- paste0("structure_loglik_", analysis.set.name, ".csv") #define expected STRUCTURE per-replicate file name
  flat.path <- file.path(structure.loglik.directory, expected.file.name) #define flat STRUCTURE CSV path
  nested.path <- file.path(structure.loglik.directory, analysis.set.name, expected.file.name) #define nested STRUCTURE CSV path
  if (file.exists(flat.path)) return(flat.path)
  if (file.exists(nested.path)) return(nested.path)
  recursive.matches <- list.files(structure.loglik.directory, pattern = paste0("^", expected.file.name, "$"), recursive = TRUE, full.names = TRUE) #search recursively for exact file name
  if (length(recursive.matches) == 1) return(recursive.matches)
  if (length(recursive.matches) > 1) stop(paste("More than one STRUCTURE log-likelihood file found for", analysis.set.name, ":", paste(recursive.matches, collapse = ", ")))
  return(NA_character_)
}


## Create function to check required columns in table
check.required.columns <- function(input.table, required.columns, table.name) {
  missing.columns <- setdiff(required.columns, colnames(input.table)) #find missing required columns
  if (length(missing.columns) > 0) stop(paste(table.name, "is missing required columns:", paste(missing.columns, collapse = ", ")))
  invisible(TRUE)
}


## Create function to parse STRUCTURE file-tag metadata
parse.STRUCTURE.file.tag.metadata <- function(file.tags) {
  STRUCTURE.metadata.list <- lapply(file.tags, function(file.tag) {
    original.file.tag <- file.tag
    file.tag <- tools::file_path_sans_ext(basename(file.tag))
    file.name <- paste0(file.tag, ".vcf")
    filename.matches <- regexec("^sim([0-9]+)_tdiv([0-9.]+)_mig(.+)$", file.tag)
    filename.parts <- regmatches(file.tag, filename.matches)[[1]] #extract matched filename parts
    if (length(filename.parts) == 0) stop(paste("STRUCTURE file_tag does not match expected pattern:", file.tag))
    mig.tag <- filename.parts[4] #extract migration-rate tag
    mig.value <- suppressWarnings(as.numeric(mig.tag)) #convert migration-rate tag to numeric
    if (is.na(mig.value)) stop(paste("Could not parse migration tag from STRUCTURE file_tag:", file.tag))
    metadata.table <- data.frame(file = file.name,
                                 mig.index = as.integer(filename.parts[2]),
                                 tdiv = as.numeric(filename.parts[3]),
                                 mig.tag = mig.tag,
                                 mig = mig.value,
                                 stringsAsFactors = FALSE)
    return(metadata.table)
  })
  STRUCTURE.metadata.table <- do.call(rbind, STRUCTURE.metadata.list) #combine metadata rows
  rownames(STRUCTURE.metadata.table) <- NULL
  return(STRUCTURE.metadata.table)
}


## Create function to create STRUCTURE mean-lnP binomial table for one analysis set
create.STRUCTURE.mean.lnprob.binomial.table <- function(analysis.set.name,
                                                        structure.loglik.directory,
                                                        results.directory,
                                                        expected.structure.replicates = 20,
                                                        STRUCTURE.mean.lnprob.delta.threshold = 0) {
  structure.file <- find.STRUCTURE.replicate.file(analysis.set.name = analysis.set.name, structure.loglik.directory = structure.loglik.directory) #find STRUCTURE per-replicate CSV
  if (is.na(structure.file)) stop(paste("No STRUCTURE per-replicate CSV found for", analysis.set.name, "in", structure.loglik.directory))
  STRUCTURE.results.directory <- file.path(results.directory, "STRUCTURE_results") #define STRUCTURE result directory for this analysis set
  dir.create(STRUCTURE.results.directory, recursive = TRUE, showWarnings = FALSE) #create STRUCTURE result directory if needed
  STRUCTURE.replicate.table <- read.csv(structure.file, stringsAsFactors = FALSE, check.names = FALSE) #read STRUCTURE per-replicate table
  STRUCTURE.replicate.table$analysis.set.name <- analysis.set.name #store analysis-set name
  check.required.columns(input.table = STRUCTURE.replicate.table,
                         required.columns = c("file_tag", "k", "rep", "estimated_ln_prob_data"),
                         table.name = structure.file) #check required STRUCTURE columns
  STRUCTURE.replicate.table$file_tag <- as.character(STRUCTURE.replicate.table$file_tag) #coerce file tags
  STRUCTURE.replicate.table$k <- as.integer(STRUCTURE.replicate.table$k) #coerce K values
  STRUCTURE.replicate.table$rep <- as.integer(STRUCTURE.replicate.table$rep) #coerce replicate values
  STRUCTURE.replicate.table$estimated_ln_prob_data <- as.numeric(STRUCTURE.replicate.table$estimated_ln_prob_data) #coerce lnP values
  if (any(is.na(STRUCTURE.replicate.table$file_tag))) stop(paste("Some STRUCTURE rows have missing file_tag for", analysis.set.name))
  if (any(is.na(STRUCTURE.replicate.table$k))) stop(paste("Some STRUCTURE rows have missing k for", analysis.set.name))
  if (any(is.na(STRUCTURE.replicate.table$rep))) stop(paste("Some STRUCTURE rows have missing rep for", analysis.set.name))
  if (any(!is.finite(STRUCTURE.replicate.table$estimated_ln_prob_data))) stop(paste("Some STRUCTURE rows have missing or non-finite estimated_ln_prob_data for", analysis.set.name))
  STRUCTURE.replicate.table <- STRUCTURE.replicate.table[STRUCTURE.replicate.table$k %in% c(1L, 2L), ] #keep K1 and K2 only
  duplicate.STRUCTURE.rows <- duplicated(STRUCTURE.replicate.table[, c("analysis.set.name", "file_tag", "k", "rep")]) #check duplicate replicate rows
  if (any(duplicate.STRUCTURE.rows)) {
    write.csv(STRUCTURE.replicate.table[duplicate.STRUCTURE.rows, c("analysis.set.name", "file_tag", "k", "rep")], file.path(STRUCTURE.results.directory, "STRUCTURE_duplicate_replicate_rows.csv"), row.names = FALSE)
    stop(paste("Duplicate STRUCTURE rows found for", analysis.set.name))
  }
  STRUCTURE.k1.reps <- STRUCTURE.replicate.table[STRUCTURE.replicate.table$k == 1L, ] #subset K1 rows
  STRUCTURE.k2.reps <- STRUCTURE.replicate.table[STRUCTURE.replicate.table$k == 2L, ] #subset K2 rows
  STRUCTURE.rep.comparison.table <- merge(STRUCTURE.k1.reps[, c("analysis.set.name", "file_tag", "rep", "estimated_ln_prob_data")],
                                          STRUCTURE.k2.reps[, c("analysis.set.name", "file_tag", "rep", "estimated_ln_prob_data")],
                                          by = c("analysis.set.name", "file_tag", "rep"),
                                          suffixes = c(".K1", ".K2")) #match K1 and K2 by VCF and replicate
  all.STRUCTURE.file.table <- unique(STRUCTURE.replicate.table[, c("analysis.set.name", "file_tag")]) #list all VCFs with STRUCTURE rows
  STRUCTURE.rep.count.check.table <- aggregate(rep ~ analysis.set.name + file_tag, data = STRUCTURE.rep.comparison.table, FUN = length) #count matched K1/K2 comparisons
  colnames(STRUCTURE.rep.count.check.table)[colnames(STRUCTURE.rep.count.check.table) == "rep"] <- "n.matched.replicates"
  STRUCTURE.rep.count.check.table <- merge(all.STRUCTURE.file.table, STRUCTURE.rep.count.check.table, by = c("analysis.set.name", "file_tag"), all.x = TRUE) #include VCFs with zero matches
  STRUCTURE.rep.count.check.table$n.matched.replicates[is.na(STRUCTURE.rep.count.check.table$n.matched.replicates)] <- 0L #replace missing counts with zero
  bad.STRUCTURE.rep.count.table <- STRUCTURE.rep.count.check.table[STRUCTURE.rep.count.check.table$n.matched.replicates != expected.structure.replicates, ] #find bad replicate counts
  write.csv(STRUCTURE.rep.count.check.table, file.path(STRUCTURE.results.directory, "STRUCTURE_matched_replicate_count_check.csv"), row.names = FALSE) #save replicate count check
  if (nrow(bad.STRUCTURE.rep.count.table) > 0) {
    write.csv(bad.STRUCTURE.rep.count.table, file.path(STRUCTURE.results.directory, "STRUCTURE_bad_matched_replicate_counts.csv"), row.names = FALSE)
    stop(paste("At least one STRUCTURE VCF does not have the expected number of matched K1/K2 replicate comparisons for", analysis.set.name))
  }
  STRUCTURE.rep.comparison.table$STRUCTURE.delta.K2.minus.K1.rep <- STRUCTURE.rep.comparison.table$estimated_ln_prob_data.K2 - STRUCTURE.rep.comparison.table$estimated_ln_prob_data.K1 #calculate replicate-level lnP difference
  STRUCTURE.rep.comparison.table$STRUCTURE.rep.K2.higher.lnprob <- as.integer(STRUCTURE.rep.comparison.table$STRUCTURE.delta.K2.minus.K1.rep > 0) #diagnostic replicate-level K2 support
  STRUCTURE.mean.lnprob.by.vcf.table <- aggregate(cbind(estimated_ln_prob_data.K1,
                                                        estimated_ln_prob_data.K2,
                                                        STRUCTURE.delta.K2.minus.K1.rep,
                                                        STRUCTURE.rep.K2.higher.lnprob) ~ analysis.set.name + file_tag,
                                                  data = STRUCTURE.rep.comparison.table,
                                                  FUN = mean) #summarize replicate means per VCF
  colnames(STRUCTURE.mean.lnprob.by.vcf.table)[colnames(STRUCTURE.mean.lnprob.by.vcf.table) == "estimated_ln_prob_data.K1"] <- "STRUCTURE.mean.lnprob.K1"
  colnames(STRUCTURE.mean.lnprob.by.vcf.table)[colnames(STRUCTURE.mean.lnprob.by.vcf.table) == "estimated_ln_prob_data.K2"] <- "STRUCTURE.mean.lnprob.K2"
  colnames(STRUCTURE.mean.lnprob.by.vcf.table)[colnames(STRUCTURE.mean.lnprob.by.vcf.table) == "STRUCTURE.delta.K2.minus.K1.rep"] <- "STRUCTURE.mean.delta.K2.minus.K1"
  colnames(STRUCTURE.mean.lnprob.by.vcf.table)[colnames(STRUCTURE.mean.lnprob.by.vcf.table) == "STRUCTURE.rep.K2.higher.lnprob"] <- "STRUCTURE.proportion.replicates.K2.higher.lnprob"
  STRUCTURE.sd.delta.by.vcf.table <- aggregate(STRUCTURE.delta.K2.minus.K1.rep ~ analysis.set.name + file_tag, data = STRUCTURE.rep.comparison.table, FUN = stats::sd) #calculate SD of replicate-level delta per VCF
  colnames(STRUCTURE.sd.delta.by.vcf.table)[colnames(STRUCTURE.sd.delta.by.vcf.table) == "STRUCTURE.delta.K2.minus.K1.rep"] <- "STRUCTURE.sd.delta.K2.minus.K1"
  STRUCTURE.n.reps.by.vcf.table <- aggregate(rep ~ analysis.set.name + file_tag, data = STRUCTURE.rep.comparison.table, FUN = length) #count matched comparisons per VCF
  colnames(STRUCTURE.n.reps.by.vcf.table)[colnames(STRUCTURE.n.reps.by.vcf.table) == "rep"] <- "STRUCTURE.n.matched.replicates"
  STRUCTURE.mean.lnprob.by.vcf.table <- merge(STRUCTURE.mean.lnprob.by.vcf.table, STRUCTURE.sd.delta.by.vcf.table, by = c("analysis.set.name", "file_tag"), all.x = TRUE) #add SD delta
  STRUCTURE.mean.lnprob.by.vcf.table <- merge(STRUCTURE.mean.lnprob.by.vcf.table, STRUCTURE.n.reps.by.vcf.table, by = c("analysis.set.name", "file_tag"), all.x = TRUE) #add replicate count
  STRUCTURE.mean.lnprob.by.vcf.table$STRUCTURE.best.k.by.mean.lnprob <- ifelse(STRUCTURE.mean.lnprob.by.vcf.table$STRUCTURE.mean.delta.K2.minus.K1 > STRUCTURE.mean.lnprob.delta.threshold, 2L, 1L) #select K from mean lnP
  STRUCTURE.mean.lnprob.by.vcf.table$STRUCTURE.n.k2.vcf.call <- as.integer(STRUCTURE.mean.lnprob.by.vcf.table$STRUCTURE.best.k.by.mean.lnprob == 2L) #binary K2 call
  STRUCTURE.mean.lnprob.by.vcf.table$STRUCTURE.n.k1.vcf.call <- as.integer(STRUCTURE.mean.lnprob.by.vcf.table$STRUCTURE.best.k.by.mean.lnprob == 1L) #binary K1 call
  STRUCTURE.metadata.table <- parse.STRUCTURE.file.tag.metadata(STRUCTURE.mean.lnprob.by.vcf.table$file_tag) #parse metadata from STRUCTURE file tags
  STRUCTURE.mean.lnprob.by.vcf.table <- cbind(STRUCTURE.mean.lnprob.by.vcf.table, STRUCTURE.metadata.table[, c("file", "mig.index", "tdiv", "mig.tag", "mig"), drop = FALSE]) #add metadata columns
  STRUCTURE.mean.lnprob.by.vcf.table <- STRUCTURE.mean.lnprob.by.vcf.table[order(STRUCTURE.mean.lnprob.by.vcf.table$mig, STRUCTURE.mean.lnprob.by.vcf.table$tdiv, STRUCTURE.mean.lnprob.by.vcf.table$file_tag), ] #order VCF-level table
  rownames(STRUCTURE.mean.lnprob.by.vcf.table) <- NULL
  STRUCTURE.binomial.table <- data.frame(method = "STRUCTURE",
                                         file = STRUCTURE.mean.lnprob.by.vcf.table$file,
                                         mig.tag = STRUCTURE.mean.lnprob.by.vcf.table$mig.tag,
                                         mig = STRUCTURE.mean.lnprob.by.vcf.table$mig,
                                         tdiv = STRUCTURE.mean.lnprob.by.vcf.table$tdiv,
                                         n.k2 = STRUCTURE.mean.lnprob.by.vcf.table$STRUCTURE.n.k2.vcf.call,
                                         n.not.k2 = STRUCTURE.mean.lnprob.by.vcf.table$STRUCTURE.n.k1.vcf.call,
                                         proportion.k2 = STRUCTURE.mean.lnprob.by.vcf.table$STRUCTURE.n.k2.vcf.call,
                                         stringsAsFactors = FALSE) #create binary STRUCTURE table from mean-lnP best K
  STRUCTURE.binomial.table <- STRUCTURE.binomial.table[order(STRUCTURE.binomial.table$mig, STRUCTURE.binomial.table$tdiv, STRUCTURE.binomial.table$file), ] #order binomial table
  rownames(STRUCTURE.binomial.table) <- NULL
  write.csv(STRUCTURE.replicate.table, file.path(STRUCTURE.results.directory, "STRUCTURE_replicate_loglik_all.csv"), row.names = FALSE) #save replicate lnP table
  write.csv(STRUCTURE.rep.comparison.table, file.path(STRUCTURE.results.directory, "STRUCTURE_replicate_K1_K2_comparison.csv"), row.names = FALSE) #save replicate comparison table
  write.csv(STRUCTURE.mean.lnprob.by.vcf.table, file.path(STRUCTURE.results.directory, "STRUCTURE_mean_lnprob_best_K_by_vcf.csv"), row.names = FALSE) #save VCF-level best K table
  write.csv(STRUCTURE.binomial.table, file.path(STRUCTURE.results.directory, "STRUCTURE_binomial_table_by_vcf.csv"), row.names = FALSE) #save STRUCTURE binomial table
  return(list(STRUCTURE.binomial.table = STRUCTURE.binomial.table,
              STRUCTURE.mean.lnprob.by.vcf.table = STRUCTURE.mean.lnprob.by.vcf.table,
              STRUCTURE.rep.comparison.table = STRUCTURE.rep.comparison.table,
              STRUCTURE.rep.count.check.table = STRUCTURE.rep.count.check.table))
}


## Create function to run SOM workflow for one VCF file
run.SOM.workflow.for.one.vcf <- function(vcf.file.path, vcf.file.metadata, sample.n.per.population = NULL, sampling.group = NULL, sampling.seed = 1) {
  genind.object.raw <- read.fastsimcoal2.genind.for.analysis(vcf.file.path = vcf.file.path,
                                                             n.per.population = sample.n.per.population,
                                                             sampling.group = sampling.group,
                                                             sampling.seed = sampling.seed) #read VCF, assign populations, and optionally subset individuals per population
  n.individuals.raw <- adegenet::nInd(genind.object.raw) #store number of individuals
  n.loci.raw <- adegenet::nLoc(genind.object.raw) #store number of loci
  genind.object.processed <- process.SNP.data.SOM(genind.input = genind.object.raw, verbose = FALSE) #process SNP data for SOM
  n.samples.processed <- if (!is.null(dim(genind.object.processed))) nrow(genind.object.processed) else NA_integer_ #store number of processed samples
  n.loci.processed <- if (!is.null(dim(genind.object.processed))) ncol(genind.object.processed) else NA_integer_ #store number of processed loci
  trained.SOM.object <- train.SOM(genind.object.processed, parallel = FALSE, verbose = FALSE, grid.multiplier = SOM_grid_multiplier) #train SOM on full processed SNP matrix
  genind.SOM.object <- clustering.SOM(trained.SOM.object, verbose = FALSE, save.SOM.results = FALSE, clustering.method = clustering_method_SOM, max.k = max_k) #cluster full SOM object after QE and TE filtering
  optim.k.summary.table <- coerce.optim.k.summary.to.data.frame(genind.SOM.object$optim_k_summary) #coerce optim_k_summary to data frame
  k.proportion.vector <- extract.k.proportions.from.summary(optim.k.summary.table, max.k = max_k) #extract row-level K proportions
  result.summary.row <- data.frame(file = vcf.file.metadata$file,
                                   mig.index = vcf.file.metadata$mig.index,
                                   tdiv = vcf.file.metadata$tdiv,
                                   mig.tag = vcf.file.metadata$mig.tag,
                                   mig = vcf.file.metadata$mig,
                                   n.ind.raw = n.individuals.raw,
                                   n.loc.raw = n.loci.raw,
                                   n.samples.processed = n.samples.processed,
                                   n.loci.processed = n.loci.processed,
                                   proportion.k2 = unname(k.proportion.vector["proportion.k2"]),
                                   status = "ok",
                                   error = NA_character_,
                                   stringsAsFactors = FALSE) #create one-row SOM K2-proportion result summary
  if (nrow(optim.k.summary.table) > 0) {
    optim.k.summary.table <- cbind(vcf.file.metadata[rep(1, nrow(optim.k.summary.table)), c("file", "mig.index", "tdiv", "mig.tag", "mig")],
                                   data.frame(n.ind.raw = rep(n.individuals.raw, nrow(optim.k.summary.table)),
                                              n.loc.raw = rep(n.loci.raw, nrow(optim.k.summary.table)),
                                              n.samples.processed = rep(n.samples.processed, nrow(optim.k.summary.table)),
                                              n.loci.processed = rep(n.loci.processed, nrow(optim.k.summary.table)),
                                              status = rep("ok", nrow(optim.k.summary.table)),
                                              error = rep(NA_character_, nrow(optim.k.summary.table)),
                                              stringsAsFactors = FALSE),
                                   optim.k.summary.table) #add metadata to optim_k_summary table
  } else {
    optim.k.summary.table <- data.frame(file = vcf.file.metadata$file,
                                        mig.index = vcf.file.metadata$mig.index,
                                        tdiv = vcf.file.metadata$tdiv,
                                        mig.tag = vcf.file.metadata$mig.tag,
                                        mig = vcf.file.metadata$mig,
                                        n.ind.raw = n.individuals.raw,
                                        n.loc.raw = n.loci.raw,
                                        n.samples.processed = n.samples.processed,
                                        n.loci.processed = n.loci.processed,
                                        status = "ok",
                                        error = NA_character_,
                                        Count = NA_real_,
                                        Proportion = NA_real_,
                                        k.label = NA_character_,
                                        optim.k.row = NA_integer_,
                                        stringsAsFactors = FALSE) #create empty optim_k_summary table if missing
  }
  return(list(result.summary.row = result.summary.row, optim.k.summary.table = optim.k.summary.table)) #return results
}


## Create function to run k-means + BIC
run.kmeans.BIC.matrix <- function(genotype.matrix,
                                  n.iter,
                                  n.start,
                                  max.n.clust = 2,
                                  n.pca,
                                  center = TRUE,
                                  scale = FALSE,
                                  plot.BIC = FALSE) {
  genotype.matrix <- as.matrix(genotype.matrix)
  storage.mode(genotype.matrix) <- "numeric"
  for (genotype.column.index in seq_len(ncol(genotype.matrix))) {
    finite.values <- genotype.matrix[is.finite(genotype.matrix[, genotype.column.index]), genotype.column.index]
    if (length(finite.values) == 0) {
      genotype.matrix[, genotype.column.index] <- 0
    } else {
      column.mean <- mean(finite.values)
      genotype.matrix[!is.finite(genotype.matrix[, genotype.column.index]), genotype.column.index] <- column.mean
    }
  }
  finite.columns <- colSums(is.finite(genotype.matrix)) == nrow(genotype.matrix)
  genotype.matrix <- genotype.matrix[, finite.columns, drop = FALSE]
  variable.columns <- apply(genotype.matrix, 2, function(genotype.column) stats::var(genotype.column, na.rm = TRUE) > 0)
  variable.columns[is.na(variable.columns)] <- FALSE
  if (!any(variable.columns)) stop("No finite variable genotype columns remain after filtering")
  genotype.matrix <- genotype.matrix[, variable.columns, drop = FALSE]
  if (any(!is.finite(genotype.matrix))) stop("genotype matrix contains NA, NaN, Inf, or -Inf values")
  genotype.data.frame <- as.data.frame(genotype.matrix)
  individual.count <- nrow(genotype.data.frame)
  max.n.clust <- min(max(max.n.clust, 2), individual.count - 1)
  cluster.numbers <- 2:max.n.clust
  pca.result <- ade4::dudi.pca(genotype.data.frame, center = center, scale = scale, scannf = FALSE, nf = min(dim(genotype.data.frame)))
  retained.pca.count <- min(length(pca.result$eig), n.pca)
  if (retained.pca.count < 1) stop("No valid PCA axes are available for de novo k-means/BIC")
  pca.scores <- pca.result$li[, 1:retained.pca.count, drop = FALSE]
  kmeans.results <- vector("list", length(cluster.numbers))
  names(kmeans.results) <- paste0("K=", cluster.numbers)
  within.sum.squares <- numeric(length(cluster.numbers))
  names(within.sum.squares) <- paste0("K=", cluster.numbers)
  for (i in seq_along(cluster.numbers)) {
    kmeans.results[[i]] <- stats::kmeans(pca.scores, centers = cluster.numbers[i], iter.max = n.iter, nstart = n.start)
    within.sum.squares[i] <- kmeans.results[[i]]$tot.withinss
  }
  original.within.sum.squares <- sum(apply(pca.scores, 2, function(pca.score.column) sum((pca.score.column - mean(pca.score.column))^2)))
  BIC.values <- individual.count * log(c(original.within.sum.squares, within.sum.squares) / individual.count) + log(individual.count) * c(1, cluster.numbers)
  names(BIC.values) <- paste0("K=", c(1, cluster.numbers))
  BIC.table <- data.frame(K = c(1, cluster.numbers), BIC = as.numeric(BIC.values))
  if (plot.BIC) {
    plot(BIC.table$K, BIC.table$BIC, type = "o", xlab = "Number of clusters", ylab = "BIC", main = "De novo k-means BIC", xaxt = "n")
    axis(1, at = BIC.table$K)
  }
  result <- list(Kstat = BIC.values,
                 BIC = BIC.table,
                 kmeans.results = kmeans.results,
                 within.sum.squares = within.sum.squares,
                 original.within.sum.squares = original.within.sum.squares,
                 pca.result = pca.result,
                 pca.scores = pca.scores,
                 retained.pca.count = retained.pca.count)
  return(result)
}


## Create function to select K1 or K2 from de novo k-means/BIC
select.deNovo.kmeans.K1.K2 <- function(BIC.values, BIC.thresh = 0) {
  if (length(BIC.values) < 2) stop("K1 and K2 BIC values are required")
  if (is.na(BIC.values[1]) || is.na(BIC.values[2])) stop("K1 or K2 BIC value is NA")
  BIC.drop.K1.to.K2 <- unname(BIC.values[1] - BIC.values[2]) #positive means K2 improves BIC
  if (BIC.drop.K1.to.K2 < BIC.thresh) {
    selected.K <- 1L #select K1 if K2 does not improve BIC enough
  } else {
    selected.K <- 2L #select K2 if BIC improves enough
  }
  return(list(selected.K = selected.K, BIC.drop.K1.to.K2 = BIC.drop.K1.to.K2))
}



## Create function to create one-row de novo k-means/BIC summary
create.deNovo.kmeans.summary.row <- function(vcf.file.metadata,
                                             deNovo.kmeans.result = NULL,
                                             deNovo.kmeans.selection = NULL,
                                             n.samples.processed = NA_integer_,
                                             n.loci.processed = NA_integer_,
                                             deNovo.kmeans.status = "ok",
                                             deNovo.kmeans.error = NA_character_) {
  if (is.null(deNovo.kmeans.result) || is.null(deNovo.kmeans.selection)) {
    deNovo.kmeans.row <- data.frame(file = vcf.file.metadata$file,
                                    deNovo.kmeans.best.k = NA_integer_,
                                    deNovo.kmeans.BIC.K1 = NA_real_,
                                    deNovo.kmeans.BIC.K2 = NA_real_,
                                    deNovo.kmeans.BIC.drop.K1.to.K2 = NA_real_,
                                    deNovo.kmeans.retained.pca.count = NA_integer_,
                                    deNovo.kmeans.n.samples.processed = n.samples.processed,
                                    deNovo.kmeans.n.loci.processed = n.loci.processed,
                                    deNovo.kmeans.status = deNovo.kmeans.status,
                                    deNovo.kmeans.error = deNovo.kmeans.error,
                                    stringsAsFactors = FALSE)
  } else {
    deNovo.kmeans.row <- data.frame(file = vcf.file.metadata$file,
                                    deNovo.kmeans.best.k = deNovo.kmeans.selection$selected.K,
                                    deNovo.kmeans.BIC.K1 = unname(deNovo.kmeans.result$Kstat["K=1"]),
                                    deNovo.kmeans.BIC.K2 = unname(deNovo.kmeans.result$Kstat["K=2"]),
                                    deNovo.kmeans.BIC.drop.K1.to.K2 = deNovo.kmeans.selection$BIC.drop.K1.to.K2,
                                    deNovo.kmeans.retained.pca.count = deNovo.kmeans.result$retained.pca.count,
                                    deNovo.kmeans.n.samples.processed = n.samples.processed,
                                    deNovo.kmeans.n.loci.processed = n.loci.processed,
                                    deNovo.kmeans.status = deNovo.kmeans.status,
                                    deNovo.kmeans.error = deNovo.kmeans.error,
                                    stringsAsFactors = FALSE)
  }
  return(deNovo.kmeans.row)
}


## Create function to run de novo k-means/BIC for one VCF file
run.deNovo.kmeans.BIC.for.one.vcf <- function(vcf.file.path,
                                              vcf.file.metadata,
                                              BIC.thresh = 0,
                                              n.iter = 10,
                                              n.start = 50,
                                              max.n.clust = 2,
                                              max.proportion.PCs = 0.9,
                                              center = TRUE,
                                              scale = FALSE,
                                              sample.n.per.population = NULL,
                                              sampling.group = NULL,
                                              sampling.seed = 1) {
  genind.object.raw <- read.fastsimcoal2.genind.for.analysis(vcf.file.path = vcf.file.path,
                                                             n.per.population = sample.n.per.population,
                                                             sampling.group = sampling.group,
                                                             sampling.seed = sampling.seed) #read VCF, assign populations, and optionally subset individuals per population
  genind.object.processed <- process.SNP.data.SOM(genind.input = genind.object.raw, verbose = FALSE) #process SNP data like SOM
  n.samples.processed <- if (!is.null(dim(genind.object.processed))) nrow(genind.object.processed) else NA_integer_ #store number of processed samples
  n.loci.processed <- if (!is.null(dim(genind.object.processed))) ncol(genind.object.processed) else NA_integer_ #store number of processed loci
  n.pca.for.kmeans <- min(ncol(genind.object.processed), nrow(genind.object.processed) - 1, floor(max.proportion.PCs * (nrow(genind.object.processed) - 1))) #maximum PCs for de novo k-means/BIC
  if (n.pca.for.kmeans < 1) stop("No valid PCA axes are available for de novo k-means/BIC")
  deNovo.kmeans.result <- run.kmeans.BIC.matrix(genotype.matrix = genind.object.processed,
                                                n.iter = n.iter,
                                                n.start = n.start,
                                                max.n.clust = max.n.clust,
                                                n.pca = n.pca.for.kmeans,
                                                center = center,
                                                scale = scale,
                                                plot.BIC = FALSE) #run de novo K1 versus K2 k-means/BIC
  deNovo.kmeans.selection <- select.deNovo.kmeans.K1.K2(BIC.values = deNovo.kmeans.result$Kstat,
                                                        BIC.thresh = BIC.thresh) #select K1 or K2
  deNovo.kmeans.row <- create.deNovo.kmeans.summary.row(vcf.file.metadata = vcf.file.metadata,
                                                        deNovo.kmeans.result = deNovo.kmeans.result,
                                                        deNovo.kmeans.selection = deNovo.kmeans.selection,
                                                        n.samples.processed = n.samples.processed,
                                                        n.loci.processed = n.loci.processed,
                                                        deNovo.kmeans.status = "ok",
                                                        deNovo.kmeans.error = NA_character_) #create one-row de novo k-means/BIC summary
  return(deNovo.kmeans.row)
}


## Create function to calculate and append de novo k-means/BIC results
calculate.and.append.deNovo.kmeans <- function(vcf.metadata.table,
                                               result.table,
                                               optim.k.result.table,
                                               intermediate.results.directory,
                                               override = FALSE,
                                               BIC.thresh = 0,
                                               n.iter = 10,
                                               n.start = 50,
                                               max.n.clust = 2,
                                               max.proportion.PCs = 0.9,
                                               center = TRUE,
                                               scale = FALSE,
                                               sample.n.per.population = NULL,
                                               sampling.seed = 1) {
  deNovo.kmeans.list <- vector("list", nrow(vcf.metadata.table)) #initialize de novo k-means/BIC list
  for (vcf.index in seq_len(nrow(vcf.metadata.table))) {
    current.metadata <- vcf.metadata.table[vcf.index, ] #extract current metadata
    if (!override && check.vcf.intermediate.deNovo.kmeans.results.available(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory)) {
      cat("\nIntermediate de novo k-means/BIC result already found for", current.metadata$file, "- loading saved result\n")
      deNovo.kmeans.list[[vcf.index]] <- load.vcf.intermediate.deNovo.kmeans.results(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory) #load saved per-VCF de novo k-means/BIC result
    } else {
      cat("\nCalculating de novo k-means/BIC for", vcf.index, "of", nrow(vcf.metadata.table), ":", current.metadata$file, "\n")
      deNovo.kmeans.list[[vcf.index]] <- tryCatch({
        run.deNovo.kmeans.BIC.for.one.vcf(vcf.file.path = current.metadata$full.path,
                                          vcf.file.metadata = current.metadata,
                                          BIC.thresh = BIC.thresh,
                                          n.iter = n.iter,
                                          n.start = n.start,
                                          max.n.clust = max.n.clust,
                                          max.proportion.PCs = max.proportion.PCs,
                                          center = center,
                                          scale = scale,
                                          sample.n.per.population = sample.n.per.population,
                                          sampling.group = current.metadata$simulation.group,
                                          sampling.seed = sampling.seed)
      }, error = function(e) {
        create.deNovo.kmeans.summary.row(vcf.file.metadata = current.metadata,
                                         deNovo.kmeans.result = NULL,
                                         deNovo.kmeans.selection = NULL,
                                         deNovo.kmeans.status = "error",
                                         deNovo.kmeans.error = conditionMessage(e))
      })
      save.vcf.intermediate.deNovo.kmeans.results(deNovo.kmeans.row = deNovo.kmeans.list[[vcf.index]], vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory) #save per-VCF de novo k-means/BIC result immediately after calculation
      stop.if.status.not.ok(status.value = deNovo.kmeans.list[[vcf.index]]$deNovo.kmeans.status, error.value = deNovo.kmeans.list[[vcf.index]]$deNovo.kmeans.error, file.value = current.metadata$file,  method.label = "de novo k-means/BIC") #stop immediately if de novo k-means/BIC failed
    }
  }
  deNovo.kmeans.table <- do.call(rbind, deNovo.kmeans.list) #combine de novo k-means/BIC rows
  rownames(deNovo.kmeans.table) <- NULL
  deNovo.kmeans.columns <- setdiff(colnames(deNovo.kmeans.table), "file")
  result.table <- result.table[, !(colnames(result.table) %in% deNovo.kmeans.columns), drop = FALSE] #remove existing de novo k-means/BIC columns if present
  optim.k.result.table <- optim.k.result.table[, !(colnames(optim.k.result.table) %in% deNovo.kmeans.columns), drop = FALSE] #remove existing de novo k-means/BIC columns if present
  for (deNovo.kmeans.column in deNovo.kmeans.columns) {
    result.table[[deNovo.kmeans.column]] <- deNovo.kmeans.table[[deNovo.kmeans.column]][match(result.table$file, deNovo.kmeans.table$file)] #add de novo k-means/BIC column to main results
    optim.k.result.table[[deNovo.kmeans.column]] <- deNovo.kmeans.table[[deNovo.kmeans.column]][match(optim.k.result.table$file, deNovo.kmeans.table$file)] #add de novo k-means/BIC column to optim_k results
  }
  return(list(result.table = result.table, optim.k.result.table = optim.k.result.table, deNovo.kmeans.table = deNovo.kmeans.table))
}


## Create function to calculate overall Weir and Cockerham Fst
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


## Create function to calculate and append Weir and Cockerham Fst
calculate.and.append.Fst <- function(vcf.metadata.table, result.table, optim.k.result.table, intermediate.results.directory, override = FALSE, sample.n.per.population = NULL, sampling.seed = 1) {
  fst.list <- vector("list", nrow(vcf.metadata.table)) #initialize Fst list
  for (vcf.index in seq_len(nrow(vcf.metadata.table))) {
    current.metadata <- vcf.metadata.table[vcf.index, ] #extract current metadata
    if (!override && check.vcf.intermediate.Fst.results.available(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory)) {
      cat("\nIntermediate Fst result already found for", current.metadata$file, "- loading saved Fst result\n")
      fst.list[[vcf.index]] <- load.vcf.intermediate.Fst.results(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory) #load saved per-VCF Fst result
    } else {
      cat("\nCalculating Weir and Cockerham Fst for", vcf.index, "of", nrow(vcf.metadata.table), ":", current.metadata$file, "\n")
      fst.list[[vcf.index]] <- tryCatch({
        genind.object.raw <- read.fastsimcoal2.genind.for.analysis(vcf.file.path = current.metadata$full.path,
                                                                   n.per.population = sample.n.per.population,
                                                                   sampling.group = current.metadata$simulation.group,
                                                                   sampling.seed = sampling.seed) #read VCF, assign populations, and optionally subset individuals per population
        fst <- calculate.overall.Fst(genind.object = genind.object.raw) #calculate overall Weir and Cockerham Fst
        data.frame(file = current.metadata$file,
                   fst = fst,
                   fst.status = "ok",
                   fst.error = NA_character_,
                   stringsAsFactors = FALSE)
      }, error = function(e) {
        data.frame(file = current.metadata$file,
                   fst = NA_real_,
                   fst.status = "error",
                   fst.error = conditionMessage(e),
                   stringsAsFactors = FALSE)
      })
      save.vcf.intermediate.Fst.results(fst.row = fst.list[[vcf.index]], vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory) #save per-VCF intermediate Fst result immediately after calculation
    }
    stop.if.status.not.ok(status.value = fst.list[[vcf.index]]$fst.status, error.value = fst.list[[vcf.index]]$fst.error, file.value = current.metadata$file, method.label = "Weir and Cockerham Fst") #stop immediately if Fst failed
  }
  fst.table <- do.call(rbind, fst.list) #combine Fst rows
  rownames(fst.table) <- NULL
  result.table <- result.table[, !(colnames(result.table) %in% c("fst", "fst.status", "fst.error")), drop = FALSE] #remove existing Fst columns if present
  result.table$fst <- fst.table$fst[match(result.table$file, fst.table$file)] #add Fst values
  result.table$fst.status <- fst.table$fst.status[match(result.table$file, fst.table$file)] #add Fst status
  result.table$fst.error <- fst.table$fst.error[match(result.table$file, fst.table$file)] #add Fst error message
  optim.k.result.table <- optim.k.result.table[, !(colnames(optim.k.result.table) %in% c("fst", "fst.status", "fst.error")), drop = FALSE] #remove existing Fst columns if present
  optim.k.result.table$fst <- fst.table$fst[match(optim.k.result.table$file, fst.table$file)] #add Fst values
  optim.k.result.table$fst.status <- fst.table$fst.status[match(optim.k.result.table$file, fst.table$file)] #add Fst status
  optim.k.result.table$fst.error <- fst.table$fst.error[match(optim.k.result.table$file, fst.table$file)] #add Fst error message
  return(list(result.table = result.table, optim.k.result.table = optim.k.result.table, fst.table = fst.table))
}



#### Create STRUCTURE sample-exclusion CSVs ####################################
if (!("simulation.group" %in% colnames(vcf.metadata.table))) vcf.metadata.table$simulation.group <- match(vcf.metadata.table$tdiv, sort(unique(vcf.metadata.table$tdiv)))
samples.to.exclude.symmetric_24 <- create.samples.to.exclude.table(vcf.metadata.table = vcf.metadata.table,
                                                                   sample.n.by.population = analysis.sample.n.per.population.symmetric_24,
                                                                   analysis.set.name = "symmetric_24",
                                                                   sampling.seed = analysis.sample.random.seed)
samples.to.exclude.symmetric_16 <- create.samples.to.exclude.table(vcf.metadata.table = vcf.metadata.table,
                                                                   sample.n.by.population = analysis.sample.n.per.population.symmetric_16,
                                                                   analysis.set.name = "symmetric_16",
                                                                   sampling.seed = analysis.sample.random.seed)
samples.to.exclude.symmetric_8 <- create.samples.to.exclude.table(vcf.metadata.table = vcf.metadata.table,
                                                                  sample.n.by.population = analysis.sample.n.per.population.symmetric_8,
                                                                  analysis.set.name = "symmetric_8",
                                                                  sampling.seed = analysis.sample.random.seed)
samples.to.exclude.asymmetric_24 <- create.samples.to.exclude.table(vcf.metadata.table = vcf.metadata.table,
                                                                    sample.n.by.population = analysis.sample.n.by.population.asymmetric_24,
                                                                    analysis.set.name = "asymmetric_24",
                                                                    sampling.seed = analysis.sample.random.seed)
samples.to.exclude.asymmetric_16 <- create.samples.to.exclude.table(vcf.metadata.table = vcf.metadata.table,
                                                                    sample.n.by.population = analysis.sample.n.by.population.asymmetric_16,
                                                                    analysis.set.name = "asymmetric_16",
                                                                    sampling.seed = analysis.sample.random.seed)
samples.to.exclude.asymmetric_8 <- create.samples.to.exclude.table(vcf.metadata.table = vcf.metadata.table,
                                                                   sample.n.by.population = analysis.sample.n.by.population.asymmetric_8,
                                                                   analysis.set.name = "asymmetric_8",
                                                                   sampling.seed = analysis.sample.random.seed)
samples.to.exclude.for.STRUCTURE <- rbind(samples.to.exclude.symmetric_24,
                                          samples.to.exclude.symmetric_16,
                                          samples.to.exclude.symmetric_8,
                                          samples.to.exclude.asymmetric_24,
                                          samples.to.exclude.asymmetric_16,
                                          samples.to.exclude.asymmetric_8)
write.csv(samples.to.exclude.symmetric_24, file.path(results.root.directory, "fastsimcoal2_samples_to_exclude_symmetric_24.csv"), row.names = FALSE)
write.csv(samples.to.exclude.symmetric_16, file.path(results.root.directory, "fastsimcoal2_samples_to_exclude_symmetric_16.csv"), row.names = FALSE)
write.csv(samples.to.exclude.symmetric_8, file.path(results.root.directory, "fastsimcoal2_samples_to_exclude_symmetric_8.csv"), row.names = FALSE)
write.csv(samples.to.exclude.asymmetric_24, file.path(results.root.directory, "fastsimcoal2_samples_to_exclude_asymmetric_24.csv"), row.names = FALSE)
write.csv(samples.to.exclude.asymmetric_16, file.path(results.root.directory, "fastsimcoal2_samples_to_exclude_asymmetric_16.csv"), row.names = FALSE)
write.csv(samples.to.exclude.asymmetric_8, file.path(results.root.directory, "fastsimcoal2_samples_to_exclude_asymmetric_8.csv"), row.names = FALSE)
write.csv(samples.to.exclude.for.STRUCTURE, file.path(results.root.directory, "fastsimcoal2_samples_to_exclude_for_STRUCTURE.csv"), row.names = FALSE)
cat("STRUCTURE sample-exclusion rows for symmetric_24:", nrow(samples.to.exclude.symmetric_24), "\n")
cat("STRUCTURE sample-exclusion rows for symmetric_16:", nrow(samples.to.exclude.symmetric_16), "\n")
cat("STRUCTURE sample-exclusion rows for symmetric_8:", nrow(samples.to.exclude.symmetric_8), "\n")
cat("STRUCTURE sample-exclusion rows for asymmetric_24:", nrow(samples.to.exclude.asymmetric_24), "\n")
cat("STRUCTURE sample-exclusion rows for asymmetric_16:", nrow(samples.to.exclude.asymmetric_16), "\n")
cat("STRUCTURE sample-exclusion rows for asymmetric_8:", nrow(samples.to.exclude.asymmetric_8), "\n")
cat("Combined STRUCTURE sample-exclusion rows:", nrow(samples.to.exclude.for.STRUCTURE), "\n")




#### Run workflow for VCF files - symmetric_24 #################################
SOM_grid_multiplier <- 4
analysis.sample.n.per.population <- analysis.sample.n.per.population.symmetric_24 #set individuals per population for symmetric_24
analysis.set.name <- paste0("symmetric_", analysis.sample.n.per.population) #define analysis set name from sample size
results.directory <- file.path(results.root.directory, analysis.set.name) #define output directory for this analysis set
dir.create(results.directory, recursive = TRUE, showWarnings = FALSE) #create output directory if needed
intermediate.results.directory <- file.path(results.directory, "intermediate_results") #define intermediate results directory for this analysis set
dir.create(intermediate.results.directory, recursive = TRUE, showWarnings = FALSE) #create intermediate results directory if needed
cat("\n\nRunning analysis set:", analysis.set.name, "\n")
cat("Individuals per population:", analysis.sample.n.per.population, "\n")
cat("Results directory:", results.directory, "\n")
cat("Random sampling base seed:", analysis.sample.random.seed, "\n")


## List VCF files
vcf.file.pattern <- "^sim[0-9]+_tdiv[0-9.]+_mig.+\\.vcf$" #define VCF filename pattern
vcf.file.paths <- list.files(vcf.directory, pattern = vcf.file.pattern, full.names = TRUE) #list VCF files only
if (length(vcf.file.paths) == 0) stop("No VCF files found")
vcf.file.names <- basename(vcf.file.paths) #extract file names


## Parse filename metadata
vcf.metadata.list <- lapply(vcf.file.names, parse.vcf.filename.metadata) #parse metadata
vcf.metadata.table <- do.call(rbind, vcf.metadata.list) #combine metadata rows
vcf.metadata.table$full.path <- vcf.file.paths #add full path


## Order metadata table
vcf.metadata.table <- vcf.metadata.table[order(vcf.metadata.table$mig.index,
                                               vcf.metadata.table$tdiv,
                                               vcf.metadata.table$mig), ] #order rows
rownames(vcf.metadata.table) <- NULL
vcf.metadata.table$simulation.group <- match(vcf.metadata.table$tdiv, sort(unique(vcf.metadata.table$tdiv))) #define matched sampling group across migration rates for each divergence time
head(vcf.metadata.table) #show metadata for files


## Create empty lists to store results
result.summary.list <- vector("list", nrow(vcf.metadata.table)) #initialize result summary list
optim.k.summary.list <- vector("list", nrow(vcf.metadata.table)) #initialize optim_k_summary list


## Loop over VCF files
for (vcf.index in seq_len(nrow(vcf.metadata.table))) {
  current.metadata <- vcf.metadata.table[vcf.index, ] #extract current metadata
  cat("\nProcessing", vcf.index, "of", nrow(vcf.metadata.table), ":", current.metadata$file, "\n")
  if (!override && check.vcf.intermediate.results.available(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory)) {
    cat("Intermediate results already found for", current.metadata$file, "- loading saved results\n")
    current.output <- load.vcf.intermediate.results(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory) #load saved per-VCF results
  } else {
    current.output <- tryCatch({
      run.SOM.workflow.for.one.vcf(vcf.file.path = current.metadata$full.path,
                                   vcf.file.metadata = current.metadata,
                                   sample.n.per.population = analysis.sample.n.per.population,
                                   sampling.group = current.metadata$simulation.group,
                                   sampling.seed = analysis.sample.random.seed)
    }, error = function(e) {
      result.summary.row <- data.frame(file = current.metadata$file,
                                       mig.index = current.metadata$mig.index,
                                       tdiv = current.metadata$tdiv,
                                       mig.tag = current.metadata$mig.tag,
                                       mig = current.metadata$mig,
                                       n.ind.raw = NA_integer_,
                                       n.loc.raw = NA_integer_,
                                       n.samples.processed = NA_integer_,
                                       n.loci.processed = NA_integer_,
                                       proportion.k2 = NA_real_,
                                       status = "error",
                                       error = conditionMessage(e),
                                       stringsAsFactors = FALSE)
      optim.k.summary.table <- data.frame(file = current.metadata$file,
                                          mig.index = current.metadata$mig.index,
                                          tdiv = current.metadata$tdiv,
                                          mig.tag = current.metadata$mig.tag,
                                          mig = current.metadata$mig,
                                          n.ind.raw = NA_integer_,
                                          n.loc.raw = NA_integer_,
                                          n.samples.processed = NA_integer_,
                                          n.loci.processed = NA_integer_,
                                          status = "error",
                                          error = conditionMessage(e),
                                          Count = NA_real_,
                                          Proportion = NA_real_,
                                          k.label = NA_character_,
                                          optim.k.row = NA_integer_,
                                          stringsAsFactors = FALSE)
      list(result.summary.row = result.summary.row, optim.k.summary.table = optim.k.summary.table)
    })
    save.vcf.intermediate.results(result.summary.row = current.output$result.summary.row,
                                  optim.k.summary.table = current.output$optim.k.summary.table,
                                  vcf.file.metadata = current.metadata,
                                  intermediate.results.directory = intermediate.results.directory) #save per-VCF results immediately
  }
  result.summary.list[[vcf.index]] <- current.output$result.summary.row #store result summary row
  optim.k.summary.list[[vcf.index]] <- current.output$optim.k.summary.table #store optim_k_summary table
  stop.if.status.not.ok(status.value = current.output$result.summary.row$status,
                        error.value = current.output$result.summary.row$error,
                        file.value = current.metadata$file,
                        method.label = "SOM") #stop immediately if SOM failed
}


## Combine SOM results
result.table <- do.call(rbind, result.summary.list) #combine result summary rows
rownames(result.table) <- NULL
optim.k.result.table <- do.call(rbind, optim.k.summary.list) #combine optim_k_summary rows
rownames(optim.k.result.table) <- NULL


## Calculate Fst and append to SOM result tables
Fst.output <- calculate.and.append.Fst(vcf.metadata.table = vcf.metadata.table,
                                       result.table = result.table,
                                       optim.k.result.table = optim.k.result.table,
                                       intermediate.results.directory = intermediate.results.directory,
                                       override = override,
                                       sample.n.per.population = analysis.sample.n.per.population,
                                       sampling.seed = analysis.sample.random.seed) #calculate and append Fst
result.table <- Fst.output$result.table #update result table with Fst
optim.k.result.table <- Fst.output$optim.k.result.table #update optim_k result table with Fst
fst.table <- Fst.output$fst.table #extract Fst table


## Calculate de novo k-means/BIC and append to result tables
deNovo.kmeans.output <- calculate.and.append.deNovo.kmeans(vcf.metadata.table = vcf.metadata.table,
                                                           result.table = result.table,
                                                           optim.k.result.table = optim.k.result.table,
                                                           intermediate.results.directory = intermediate.results.directory,
                                                           override = override,
                                                           BIC.thresh = deNovo.kmeans.BIC.thresh,
                                                           n.iter = deNovo.kmeans.n.iter,
                                                           n.start = deNovo.kmeans.n.start,
                                                           max.n.clust = deNovo.kmeans.max.n.clust,
                                                           max.proportion.PCs = deNovo.kmeans.max.proportion.PCs,
                                                           center = deNovo.kmeans.center,
                                                           scale = deNovo.kmeans.scale,
                                                           sample.n.per.population = analysis.sample.n.per.population,
                                                           sampling.seed = analysis.sample.random.seed) #calculate de novo k-means/BIC
result.table <- deNovo.kmeans.output$result.table #update result table with de novo k-means/BIC
optim.k.result.table <- deNovo.kmeans.output$optim.k.result.table #update optim_k result table with de novo k-means/BIC
deNovo.kmeans.table <- deNovo.kmeans.output$deNovo.kmeans.table #extract de novo k-means/BIC table


## Calculate sNMF and append to result tables
sNMF.output <- calculate.and.append.sNMF(vcf.metadata.table = vcf.metadata.table,
                                         result.table = result.table,
                                         optim.k.result.table = optim.k.result.table,
                                         intermediate.results.directory = intermediate.results.directory,
                                         override = override,
                                         K.values = sNMF.K.values,
                                         repetitions = sNMF.repetitions,
                                         CPU = 1,
                                         ploidy = sNMF.ploidy,
                                         cross.entropy.thresh = sNMF.cross.entropy.thresh,
                                         seed = sNMF.seed,
                                         sample.n.per.population = analysis.sample.n.per.population,
                                         sampling.seed = analysis.sample.random.seed) #calculate sNMF
result.table <- sNMF.output$result.table #update result table with sNMF
optim.k.result.table <- sNMF.output$optim.k.result.table #update optim_k result table with sNMF
sNMF.table <- sNMF.output$sNMF.table #extract sNMF table


## Save result tables
write.csv(result.table, file.path(results.directory, "fastsimcoal2_SOM_results.csv"), row.names = FALSE) #save main SOM result summary
write.csv(optim.k.result.table, file.path(results.directory, "fastsimcoal2_SOM_optim_k_summary.csv"), row.names = FALSE) #save full optim_k_summary table
write.csv(fst.table, file.path(results.directory, "fastsimcoal2_Weir_Cockerham_Fst.csv"), row.names = FALSE) #save Fst table
write.csv(deNovo.kmeans.table, file.path(results.directory, "fastsimcoal2_de_novo_kmeans_BIC.csv"), row.names = FALSE) #save de novo k-means/BIC table
write.csv(sNMF.table, file.path(results.directory, "fastsimcoal2_sNMF_results.csv"), row.names = FALSE) #save sNMF table


## Inspect results
head(result.table)
optim.k.result.table
head(fst.table)
head(deNovo.kmeans.table)
head(sNMF.table)
table(result.table$status)
table(result.table$fst.status)
table(result.table$deNovo.kmeans.status)
table(result.table$deNovo.kmeans.best.k, useNA = "ifany")
table(result.table$sNMF.status)


## Plot fitted K2 proportions for SOM, de novo k-means/BIC, and sNMF
plot.result.table <- result.table[result.table$status == "ok", ] #subset successful rows
plot.result.table <- plot.result.table[order(plot.result.table$mig, plot.result.table$tdiv), ] #order rows for plotting
if (nrow(plot.result.table) == 0) stop("No successful results available for fitted plotting")
mig.level.table <- unique(plot.result.table[order(plot.result.table$mig), c("mig.tag", "mig")]) #define unique ordered migration-rate levels
mig.level.table$mig.tag <- as.character(mig.level.table$mig.tag) #force migration-rate labels to discrete character values
plot.result.table$mig.tag <- factor(as.character(plot.result.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels by numeric migration rate
mig.colors <- setNames(viridisLite::viridis(nrow(mig.level.table), option = "D"), mig.level.table$mig.tag) #define discrete ordered colors


## Create SOM binomial count table
SOM.optim.k.count.table <- optim.k.result.table[!is.na(optim.k.result.table$Count) & optim.k.result.table$file %in% plot.result.table$file, ] #keep rows with SOM replicate counts
SOM.optim.k.count.table$Count <- as.numeric(as.character(SOM.optim.k.count.table$Count)) #force counts to numeric
SOM.total.count.table <- stats::aggregate(Count ~ file + mig.tag + mig + tdiv, data = SOM.optim.k.count.table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE)) #total retained SOM replicates per VCF
colnames(SOM.total.count.table)[colnames(SOM.total.count.table) == "Count"] <- "n.total"
SOM.k2.count.table <- SOM.optim.k.count.table[SOM.optim.k.count.table$k.label == "k2", ] #keep k2 count rows
if (nrow(SOM.k2.count.table) > 0) {
  SOM.k2.count.table <- stats::aggregate(Count ~ file + mig.tag + mig + tdiv, data = SOM.k2.count.table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE)) #number of retained SOM replicates choosing k = 2
  colnames(SOM.k2.count.table)[colnames(SOM.k2.count.table) == "Count"] <- "n.k2"
} else {
  SOM.k2.count.table <- SOM.total.count.table[, c("file", "mig.tag", "mig", "tdiv")] #create empty k2 count table if no k2 rows exist
  SOM.k2.count.table$n.k2 <- 0
}
SOM.binomial.table <- merge(SOM.total.count.table, SOM.k2.count.table, by = c("file", "mig.tag", "mig", "tdiv"), all.x = TRUE) #merge total and k2 counts
SOM.binomial.table$n.k2[is.na(SOM.binomial.table$n.k2)] <- 0 #set missing k2 rows to zero
SOM.binomial.table$n.not.k2 <- SOM.binomial.table$n.total - SOM.binomial.table$n.k2 #number not choosing k = 2
SOM.binomial.table$proportion.k2 <- SOM.binomial.table$n.k2 / SOM.binomial.table$n.total #calculate K2 proportion
SOM.binomial.table$method <- "SOM" #add method label
SOM.binomial.table <- SOM.binomial.table[, c("method", "file", "mig.tag", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")] #order columns


## Create de novo k-means/BIC binomial table
deNovo.kmeans.binomial.source.table <- result.table[result.table$status == "ok" & !is.na(result.table$deNovo.kmeans.best.k), ] #subset valid de novo k-means/BIC rows
if ("deNovo.kmeans.status" %in% colnames(deNovo.kmeans.binomial.source.table)) deNovo.kmeans.binomial.source.table <- deNovo.kmeans.binomial.source.table[deNovo.kmeans.binomial.source.table$deNovo.kmeans.status == "ok", ] #keep successful de novo k-means/BIC rows
deNovo.kmeans.best.k <- as.integer(as.character(deNovo.kmeans.binomial.source.table$deNovo.kmeans.best.k)) #extract de novo k-means/BIC best K
deNovo.kmeans.binomial.table <- data.frame(method = "deNovo.kmeans.BIC",
                                           file = deNovo.kmeans.binomial.source.table$file,
                                           mig.tag = deNovo.kmeans.binomial.source.table$mig.tag,
                                           mig = deNovo.kmeans.binomial.source.table$mig,
                                           tdiv = deNovo.kmeans.binomial.source.table$tdiv,
                                           n.k2 = as.integer(deNovo.kmeans.best.k == 2L),
                                           n.not.k2 = as.integer(deNovo.kmeans.best.k != 2L),
                                           proportion.k2 = as.numeric(deNovo.kmeans.best.k == 2L),
                                           stringsAsFactors = FALSE) #create binary de novo k-means/BIC table


## Create sNMF binomial table
sNMF.binomial.source.table <- result.table[result.table$status == "ok" & !is.na(result.table$sNMF.best.k), ] #subset valid sNMF rows
if ("sNMF.status" %in% colnames(sNMF.binomial.source.table)) sNMF.binomial.source.table <- sNMF.binomial.source.table[sNMF.binomial.source.table$sNMF.status == "ok", ] #keep successful sNMF rows
sNMF.best.k <- as.integer(as.character(sNMF.binomial.source.table$sNMF.best.k)) #extract sNMF best K
sNMF.binomial.table <- data.frame(method = "sNMF",
                                  file = sNMF.binomial.source.table$file,
                                  mig.tag = sNMF.binomial.source.table$mig.tag,
                                  mig = sNMF.binomial.source.table$mig,
                                  tdiv = sNMF.binomial.source.table$tdiv,
                                  n.k2 = as.integer(sNMF.best.k == 2L),
                                  n.not.k2 = as.integer(sNMF.best.k != 2L),
                                  proportion.k2 = as.numeric(sNMF.best.k == 2L),
                                  stringsAsFactors = FALSE) #create binary sNMF table


## Create STRUCTURE binomial table
STRUCTURE.results <- create.STRUCTURE.mean.lnprob.binomial.table(analysis.set.name = analysis.set.name,
                                                                 structure.loglik.directory = STRUCTURE.loglik.directory,
                                                                 results.directory = results.directory,
                                                                 expected.structure.replicates = expected.structure.replicates,
                                                                 STRUCTURE.mean.lnprob.delta.threshold = STRUCTURE.mean.lnprob.delta.threshold) #create STRUCTURE mean-lnP best-K table
STRUCTURE.binomial.table <- STRUCTURE.results$STRUCTURE.binomial.table #extract STRUCTURE binomial table
STRUCTURE.binomial.table$mig.tag <- vapply(STRUCTURE.binomial.table$mig, function(current.mig) {
  matching.mig.index <- which(abs(mig.level.table$mig - current.mig) < sqrt(.Machine$double.eps))
  if (length(matching.mig.index) != 1) return(NA_character_)
  return(as.character(mig.level.table$mig.tag[matching.mig.index]))
}, character(1)) #standardize STRUCTURE migration labels to match VCF migration labels
if (any(is.na(STRUCTURE.binomial.table$mig.tag))) stop("Some STRUCTURE migration rates could not be matched to VCF migration labels")
STRUCTURE.mean.lnprob.by.vcf.table <- STRUCTURE.results$STRUCTURE.mean.lnprob.by.vcf.table #extract STRUCTURE VCF-level table
STRUCTURE.rep.comparison.table <- STRUCTURE.results$STRUCTURE.rep.comparison.table #extract STRUCTURE replicate comparison table


## Combine method-specific binomial tables
fitted.binomial.input.table <- rbind(SOM.binomial.table, deNovo.kmeans.binomial.table, sNMF.binomial.table, STRUCTURE.binomial.table) #combine all methods
method.levels <- c("SOM", "deNovo.kmeans.BIC", "sNMF", "STRUCTURE") #define method order
fitted.binomial.input.table$method <- factor(as.character(fitted.binomial.input.table$method), levels = method.levels) #order methods
fitted.binomial.input.table$mig.tag <- factor(as.character(fitted.binomial.input.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels


## Fit binomial GLMs separately for each method and migration rate
fitted.GLM.prediction.list <- list() #initialize prediction list
fitted.GLM.summary.list <- list() #initialize summary list
fitted.GLM.list.index <- 1 #initialize list index
for (current.method in method.levels) {
  for (current.mig.tag in mig.level.table$mig.tag) {
    current.binomial.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == current.method & fitted.binomial.input.table$mig.tag == current.mig.tag, ] #subset one method and migration rate
    current.binomial.table <- current.binomial.table[is.finite(current.binomial.table$tdiv), ] #keep rows with finite divergence time
    if (nrow(current.binomial.table) == 0) next
    current.binomial.table <- current.binomial.table[order(current.binomial.table$tdiv), ] #order rows by divergence time
    prediction.tdiv <- seq(min(current.binomial.table$tdiv), max(current.binomial.table$tdiv), length.out = plot.fitted.prediction.n) #create prediction grid
    current.warnings <- character(0) #store GLM warnings
    if (all(current.binomial.table$n.k2 == 0)) {
      current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = 0, stringsAsFactors = FALSE) #flat zero prediction when K2 is never selected
      current.summary.row <- data.frame(method = current.method,
                                        mig.tag = as.character(current.mig.tag),
                                        mig = current.binomial.table$mig[1],
                                        n.vcf = nrow(current.binomial.table),
                                        n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                        n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                        min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                        n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                        n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                        all.zero = TRUE,
                                        all.one = FALSE,
                                        glm.fitted = FALSE,
                                        glm.converged = NA,
                                        glm.intercept = NA_real_,
                                        glm.tdiv.slope = NA_real_,
                                        tdiv.at.threshold.k2 = NA_real_,
                                        tdiv.at.threshold.k2.within.range = FALSE,
                                        glm.warning = NA_character_,
                                        glm.error = NA_character_,
                                        stringsAsFactors = FALSE) #store summary row
    } else if (all(current.binomial.table$n.not.k2 == 0)) {
      current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = 1, stringsAsFactors = FALSE) #flat one prediction when K2 is always selected
      current.summary.row <- data.frame(method = current.method,
                                        mig.tag = as.character(current.mig.tag),
                                        mig = current.binomial.table$mig[1],
                                        n.vcf = nrow(current.binomial.table),
                                        n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                        n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                        min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                        n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                        n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                        all.zero = FALSE,
                                        all.one = TRUE,
                                        glm.fitted = FALSE,
                                        glm.converged = NA,
                                        glm.intercept = NA_real_,
                                        glm.tdiv.slope = NA_real_,
                                        tdiv.at.threshold.k2 = NA_real_,
                                        tdiv.at.threshold.k2.within.range = FALSE,
                                        glm.warning = NA_character_,
                                        glm.error = NA_character_,
                                        stringsAsFactors = FALSE) #store summary row
    } else {
      current.fit <- withCallingHandlers(tryCatch(stats::glm(cbind(n.k2, n.not.k2) ~ tdiv, data = current.binomial.table, family = stats::binomial(link = "logit"), control = stats::glm.control(maxit = 100)), error = function(e) e),
                                         warning = function(w) {
                                           current.warnings <<- c(current.warnings, conditionMessage(w))
                                           invokeRestart("muffleWarning")
                                         }) #fit count-based binomial GLM and store warnings
      fit.ok <- inherits(current.fit, "glm") #check whether GLM fit succeeded
      if (fit.ok) {
        current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = stats::predict(current.fit, newdata = data.frame(tdiv = prediction.tdiv), type = "response"), stringsAsFactors = FALSE) #predict fitted K2 probability
        current.coefficients <- stats::coef(current.fit) #extract coefficients
        current.glm.intercept <- unname(current.coefficients["(Intercept)"]) #extract intercept
        current.glm.tdiv.slope <- unname(current.coefficients["tdiv"]) #extract tdiv slope
        if (is.finite(current.glm.intercept) && is.finite(current.glm.tdiv.slope) && current.glm.tdiv.slope != 0) {
          current.tdiv.at.threshold.k2 <- (stats::qlogis(plot.k2.threshold.proportion) - current.glm.intercept) / current.glm.tdiv.slope #solve logit(p) for threshold
        } else {
          current.tdiv.at.threshold.k2 <- NA_real_
        }
        current.tdiv.at.threshold.k2.within.range <- is.finite(current.tdiv.at.threshold.k2) && current.tdiv.at.threshold.k2 >= min(current.binomial.table$tdiv, na.rm = TRUE) && current.tdiv.at.threshold.k2 <= max(current.binomial.table$tdiv, na.rm = TRUE) #check whether threshold is inside simulated range
        current.summary.row <- data.frame(method = current.method,
                                          mig.tag = as.character(current.mig.tag),
                                          mig = current.binomial.table$mig[1],
                                          n.vcf = nrow(current.binomial.table),
                                          n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                          n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                          min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                          n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                          n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                          all.zero = all(current.binomial.table$proportion.k2 == 0),
                                          all.one = all(current.binomial.table$proportion.k2 == 1),
                                          glm.fitted = TRUE,
                                          glm.converged = isTRUE(current.fit$converged),
                                          glm.intercept = current.glm.intercept,
                                          glm.tdiv.slope = current.glm.tdiv.slope,
                                          tdiv.at.threshold.k2 = current.tdiv.at.threshold.k2,
                                          tdiv.at.threshold.k2.within.range = current.tdiv.at.threshold.k2.within.range,
                                          glm.warning = if (length(current.warnings) == 0) NA_character_ else paste(unique(current.warnings), collapse = " | "),
                                          glm.error = NA_character_,
                                          stringsAsFactors = FALSE) #store summary row
      } else {
        current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = NA_real_, stringsAsFactors = FALSE) #return missing predictions if GLM failed
        current.summary.row <- data.frame(method = current.method,
                                          mig.tag = as.character(current.mig.tag),
                                          mig = current.binomial.table$mig[1],
                                          n.vcf = nrow(current.binomial.table),
                                          n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                          n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                          min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                          n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                          n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                          all.zero = all(current.binomial.table$proportion.k2 == 0),
                                          all.one = all(current.binomial.table$proportion.k2 == 1),
                                          glm.fitted = FALSE,
                                          glm.converged = FALSE,
                                          glm.intercept = NA_real_,
                                          glm.tdiv.slope = NA_real_,
                                          tdiv.at.threshold.k2 = NA_real_,
                                          tdiv.at.threshold.k2.within.range = FALSE,
                                          glm.warning = if (length(current.warnings) == 0) NA_character_ else paste(unique(current.warnings), collapse = " | "),
                                          glm.error = conditionMessage(current.fit),
                                          stringsAsFactors = FALSE) #store summary row
      }
    }
    fitted.GLM.prediction.list[[fitted.GLM.list.index]] <- current.prediction.table #store prediction table
    fitted.GLM.summary.list[[fitted.GLM.list.index]] <- current.summary.row #store summary row
    fitted.GLM.list.index <- fitted.GLM.list.index + 1 #advance list index
  }
}


## Combine fitted GLM outputs
fitted.GLM.prediction.table <- do.call(rbind, fitted.GLM.prediction.list) #combine prediction tables
fitted.GLM.summary.table <- do.call(rbind, fitted.GLM.summary.list) #combine summary rows
fitted.GLM.prediction.table$method <- factor(as.character(fitted.GLM.prediction.table$method), levels = method.levels) #order methods
fitted.GLM.prediction.table$mig.tag <- factor(as.character(fitted.GLM.prediction.table$mig.tag), levels = mig.level.table$mig.tag) #order migration labels
fitted.GLM.summary.table$method <- factor(as.character(fitted.GLM.summary.table$method), levels = method.levels) #order methods
fitted.GLM.summary.table$mig.tag <- factor(as.character(fitted.GLM.summary.table$mig.tag), levels = mig.level.table$mig.tag) #order migration labels
fitted.GLM.summary.table <- fitted.GLM.summary.table[order(fitted.GLM.summary.table$method, fitted.GLM.summary.table$mig), ] #order summary table
rownames(fitted.GLM.summary.table) <- NULL
fitted.GLM.threshold.table <- fitted.GLM.summary.table[fitted.GLM.summary.table$tdiv.at.threshold.k2.within.range & is.finite(fitted.GLM.summary.table$tdiv.at.threshold.k2), ] #keep valid threshold rows


## Plot SOM results
SOM.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "SOM", ] #subset SOM observations
SOM.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "SOM" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset SOM fitted predictions
SOM.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "SOM", ] #subset SOM threshold rows
SOM.k2.fitted.plot <- ggplot2::ggplot(SOM.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = SOM.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = SOM.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(SOM.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "SOM_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = SOM.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot Weir and Cockerham Fst results
fst.plot.table <- plot.result.table[is.finite(plot.result.table$fst), ] #subset rows with finite Fst values
fst.plot.table <- fst.plot.table[order(fst.plot.table$mig, fst.plot.table$tdiv), ] #order rows for plotting
fst.plot.table$mig.tag <- factor(as.character(fst.plot.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels
SOM.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "SOM" & fitted.GLM.threshold.table$tdiv.at.threshold.k2.within.range & is.finite(fitted.GLM.threshold.table$tdiv.at.threshold.k2), ] #subset valid SOM threshold rows
SOM.k2.threshold.plot.table$mig.tag <- factor(as.character(SOM.k2.threshold.plot.table$mig.tag), levels = mig.level.table$mig.tag) #order SOM threshold migration-rate labels
fst.plot <- ggplot2::ggplot(fst.plot.table, ggplot2::aes(x = tdiv, y = fst, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = Fst_plot_point_size, alpha = Fst_plot_point_alpha) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, linewidth = Fst_plot_line_width, se = FALSE) +
  ggplot2::geom_vline(data = SOM.k2.threshold.plot.table,
                      ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag),
                      linewidth = plot.threshold.line.width,
                      linetype = plot.threshold.line.type,
                      show.legend = FALSE,
                      inherit.aes = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Weir and Cockerham Fst", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(fst.plot)
ggplot2::ggsave(file.path(results.directory, "Weir_Cockerham_Fst_by_tdiv_and_mig.svg"),
                plot = fst.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot de novo k-means/BIC results
deNovo.kmeans.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "deNovo.kmeans.BIC", ] #subset de novo k-means/BIC observations
deNovo.kmeans.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "deNovo.kmeans.BIC" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset de novo k-means/BIC fitted predictions
deNovo.kmeans.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "deNovo.kmeans.BIC", ] #subset de novo k-means/BIC threshold rows
deNovo.kmeans.k2.fitted.plot <- ggplot2::ggplot(deNovo.kmeans.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = deNovo.kmeans.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = deNovo.kmeans.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(deNovo.kmeans.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "deNovo_kmeans_BIC_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = deNovo.kmeans.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot sNMF K2 results
sNMF.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "sNMF", ] #subset sNMF observations
sNMF.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "sNMF" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset sNMF fitted predictions
sNMF.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "sNMF", ] #subset sNMF threshold rows
sNMF.k2.fitted.plot <- ggplot2::ggplot(sNMF.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = sNMF.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = sNMF.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(sNMF.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "sNMF_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = sNMF.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot STRUCTURE K2 results
STRUCTURE.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "STRUCTURE", ] #subset STRUCTURE observations
STRUCTURE.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "STRUCTURE" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset STRUCTURE fitted predictions
STRUCTURE.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "STRUCTURE", ] #subset STRUCTURE threshold rows
STRUCTURE.k2.fitted.plot <- ggplot2::ggplot(STRUCTURE.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = STRUCTURE.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = STRUCTURE.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(STRUCTURE.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "STRUCTURE_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = STRUCTURE.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Inspect fitted 50% threshold results
fitted.GLM.summary.table
fitted.GLM.threshold.table



## Create PCA plot for one VCF file

## Show available migration settings
available.migration.table <- unique(vcf.metadata.table[, c("mig.tag", "mig")]) #show all available migration settings including mig = 0
available.migration.table <- available.migration.table[order(available.migration.table$mig), ]
rownames(available.migration.table) <- NULL
print(available.migration.table)
if (!(migration.vcf.mig.for.PCA %in% unique(vcf.metadata.table$mig))) stop(paste("migration.vcf.mig.for.PCA must be one of:", paste(sort(unique(vcf.metadata.table$mig)), collapse = ", ")))


## Select file number N among VCF files for selected migration rate
migration.vcf.metadata.table <- vcf.metadata.table[vcf.metadata.table$mig == migration.vcf.mig.for.PCA, ] #keep selected VCF files only
migration.vcf.metadata.table <- migration.vcf.metadata.table[order(migration.vcf.metadata.table$tdiv), ] #order selected files by divergence time
rownames(migration.vcf.metadata.table) <- NULL
if (nrow(migration.vcf.metadata.table) == 0) stop(paste("No VCF files found for migration rate", migration.vcf.mig.for.PCA))
if (migration.vcf.file.number.for.PCA < 1 || migration.vcf.file.number.for.PCA > nrow(migration.vcf.metadata.table)) stop(paste("migration.vcf.file.number.for.PCA must be between 1 and", nrow(migration.vcf.metadata.table), "for migration rate", migration.vcf.mig.for.PCA))
current.pca.metadata <- migration.vcf.metadata.table[migration.vcf.file.number.for.PCA, ] #select VCF file N for selected migration rate
cat("\nCreating PCA for VCF file", migration.vcf.file.number.for.PCA, "of", nrow(migration.vcf.metadata.table), "for migration rate", migration.vcf.mig.for.PCA, ":", current.pca.metadata$file, "\n")
cat("Migration-index tag:", current.pca.metadata$mig.index, "\n")
cat("Divergence time:", current.pca.metadata$tdiv, "generations\n")
cat("Migration tag:", current.pca.metadata$mig.tag, "\n")
cat("Migration rate:", current.pca.metadata$mig, "\n")


## Read VCF and process SNP matrix
pca.genind.object.raw <- read.fastsimcoal2.genind.for.analysis(vcf.file.path = current.pca.metadata$full.path,
                                                               n.per.population = analysis.sample.n.per.population,
                                                               sampling.group = current.pca.metadata$simulation.group,
                                                               sampling.seed = analysis.sample.random.seed) #read selected VCF, assign populations, and optionally subset individuals per population


## Calculate Weir and Cockerham Fst before PCA processing
pca.fst.value <- calculate.overall.Fst(genind.object = pca.genind.object.raw) #calculate overall Weir and Cockerham Fst
cat("Weir and Cockerham Fst between pop1 and pop2:", round(pca.fst.value, 4), "\n")


## Process SNP matrix for PCA
pca.snp.matrix <- process.SNP.data.SOM(genind.input = pca.genind.object.raw, verbose = FALSE) #process SNP data for PCA
pca.snp.matrix <- as.matrix(pca.snp.matrix) #coerce to matrix
storage.mode(pca.snp.matrix) <- "numeric" #ensure numeric matrix
if (is.null(rownames(pca.snp.matrix))) rownames(pca.snp.matrix) <- adegenet::indNames(pca.genind.object.raw)


## Remove invalid loci
valid.pca.loci <- apply(pca.snp.matrix, 2, function(locus.values) {
  non.missing.values <- locus.values[!is.na(locus.values)]
  length(non.missing.values) > 1 && stats::var(non.missing.values) > 0
}) #identify loci with nonzero variance
pca.snp.matrix <- pca.snp.matrix[, valid.pca.loci, drop = FALSE] #keep valid loci only
if (ncol(pca.snp.matrix) < 2) stop("Fewer than two valid loci available for PCA")


## Impute remaining missing values for PCA
if (any(is.na(pca.snp.matrix))) {
  pca.snp.matrix <- apply(pca.snp.matrix, 2, function(locus.values) {
    locus.values[is.na(locus.values)] <- mean(locus.values, na.rm = TRUE)
    locus.values
  }) #replace remaining missing values with locus means
  pca.snp.matrix <- as.matrix(pca.snp.matrix)
}


## Run PCA and create plot table
pca.object <- stats::prcomp(pca.snp.matrix, center = TRUE, scale. = FALSE) #run PCA
pca.variance <- pca.object$sdev^2 / sum(pca.object$sdev^2) #calculate variance explained
pca.plot.table <- data.frame(individual = rownames(pca.object$x),
                             PC1 = pca.object$x[, 1],
                             PC2 = pca.object$x[, 2],
                             population = as.character(adegenet::pop(pca.genind.object.raw))[match(rownames(pca.object$x), adegenet::indNames(pca.genind.object.raw))],
                             file = current.pca.metadata$file,
                             mig.index = current.pca.metadata$mig.index,
                             tdiv = current.pca.metadata$tdiv,
                             mig.tag = current.pca.metadata$mig.tag,
                             mig = current.pca.metadata$mig,
                             fst.pop1.pop2 = pca.fst.value,
                             stringsAsFactors = FALSE) #create PCA plotting table


## Plot PCA with population groupings
pca.plot <- ggplot2::ggplot(pca.plot.table, ggplot2::aes(x = PC1, y = PC2, color = population, group = population)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(title = paste0("PCA - ", current.pca.metadata$file, " | tdiv = ", round(current.pca.metadata$tdiv), " generations | mig = ", current.pca.metadata$mig, " | Fst = ", round(pca.fst.value, 4)),
                x = paste0("PC1 (", round(pca.variance[1] * 100, 2), "%)"),
                y = paste0("PC2 (", round(pca.variance[2] * 100, 2), "%)"),
                color = "Population") +
  ggplot2::theme_classic()
pca.plot


## Print selected PCA summary
cat("Divergence:", round(current.pca.metadata$tdiv / 1000000, 2), "million generations\n")
cat("Migration rate:", current.pca.metadata$mig, "\n")
cat("Weir and Cockerham Fst:", round(pca.fst.value, 4), "\n")




#### Run workflow for VCF files - symmetric_16 #################################
SOM_grid_multiplier <- 4
analysis.sample.n.per.population <- analysis.sample.n.per.population.symmetric_16 #set individuals per population for symmetric_16
analysis.set.name <- paste0("symmetric_", analysis.sample.n.per.population) #define analysis set name from sample size
results.directory <- file.path(results.root.directory, analysis.set.name) #define output directory for this analysis set
dir.create(results.directory, recursive = TRUE, showWarnings = FALSE) #create output directory if needed
intermediate.results.directory <- file.path(results.directory, "intermediate_results") #define intermediate results directory for this analysis set
dir.create(intermediate.results.directory, recursive = TRUE, showWarnings = FALSE) #create intermediate results directory if needed
cat("\n\nRunning analysis set:", analysis.set.name, "\n")
cat("Individuals per population:", analysis.sample.n.per.population, "\n")
cat("Results directory:", results.directory, "\n")
cat("Random sampling base seed:", analysis.sample.random.seed, "\n")


## List VCF files
vcf.file.pattern <- "^sim[0-9]+_tdiv[0-9.]+_mig.+\\.vcf$" #define VCF filename pattern
vcf.file.paths <- list.files(vcf.directory, pattern = vcf.file.pattern, full.names = TRUE) #list VCF files only
if (length(vcf.file.paths) == 0) stop("No VCF files found")
vcf.file.names <- basename(vcf.file.paths) #extract file names


## Parse filename metadata
vcf.metadata.list <- lapply(vcf.file.names, parse.vcf.filename.metadata) #parse metadata
vcf.metadata.table <- do.call(rbind, vcf.metadata.list) #combine metadata rows
vcf.metadata.table$full.path <- vcf.file.paths #add full path


## Order metadata table
vcf.metadata.table <- vcf.metadata.table[order(vcf.metadata.table$mig.index,
                                               vcf.metadata.table$tdiv,
                                               vcf.metadata.table$mig), ] #order rows
rownames(vcf.metadata.table) <- NULL
vcf.metadata.table$simulation.group <- match(vcf.metadata.table$tdiv, sort(unique(vcf.metadata.table$tdiv))) #define matched sampling group across migration rates for each divergence time
head(vcf.metadata.table) #show metadata for files


## Create empty lists to store results
result.summary.list <- vector("list", nrow(vcf.metadata.table)) #initialize result summary list
optim.k.summary.list <- vector("list", nrow(vcf.metadata.table)) #initialize optim_k_summary list


## Loop over VCF files
for (vcf.index in seq_len(nrow(vcf.metadata.table))) {
  current.metadata <- vcf.metadata.table[vcf.index, ] #extract current metadata
  cat("\nProcessing", vcf.index, "of", nrow(vcf.metadata.table), ":", current.metadata$file, "\n")
  if (!override && check.vcf.intermediate.results.available(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory)) {
    cat("Intermediate results already found for", current.metadata$file, "- loading saved results\n")
    current.output <- load.vcf.intermediate.results(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory) #load saved per-VCF results
  } else {
    current.output <- tryCatch({
      run.SOM.workflow.for.one.vcf(vcf.file.path = current.metadata$full.path,
                                   vcf.file.metadata = current.metadata,
                                   sample.n.per.population = analysis.sample.n.per.population,
                                   sampling.group = current.metadata$simulation.group,
                                   sampling.seed = analysis.sample.random.seed)
    }, error = function(e) {
      result.summary.row <- data.frame(file = current.metadata$file,
                                       mig.index = current.metadata$mig.index,
                                       tdiv = current.metadata$tdiv,
                                       mig.tag = current.metadata$mig.tag,
                                       mig = current.metadata$mig,
                                       n.ind.raw = NA_integer_,
                                       n.loc.raw = NA_integer_,
                                       n.samples.processed = NA_integer_,
                                       n.loci.processed = NA_integer_,
                                       proportion.k2 = NA_real_,
                                       status = "error",
                                       error = conditionMessage(e),
                                       stringsAsFactors = FALSE)
      optim.k.summary.table <- data.frame(file = current.metadata$file,
                                          mig.index = current.metadata$mig.index,
                                          tdiv = current.metadata$tdiv,
                                          mig.tag = current.metadata$mig.tag,
                                          mig = current.metadata$mig,
                                          n.ind.raw = NA_integer_,
                                          n.loc.raw = NA_integer_,
                                          n.samples.processed = NA_integer_,
                                          n.loci.processed = NA_integer_,
                                          status = "error",
                                          error = conditionMessage(e),
                                          Count = NA_real_,
                                          Proportion = NA_real_,
                                          k.label = NA_character_,
                                          optim.k.row = NA_integer_,
                                          stringsAsFactors = FALSE)
      list(result.summary.row = result.summary.row, optim.k.summary.table = optim.k.summary.table)
    })
    save.vcf.intermediate.results(result.summary.row = current.output$result.summary.row,
                                  optim.k.summary.table = current.output$optim.k.summary.table,
                                  vcf.file.metadata = current.metadata,
                                  intermediate.results.directory = intermediate.results.directory) #save per-VCF results immediately
  }
  result.summary.list[[vcf.index]] <- current.output$result.summary.row #store result summary row
  optim.k.summary.list[[vcf.index]] <- current.output$optim.k.summary.table #store optim_k_summary table
  stop.if.status.not.ok(status.value = current.output$result.summary.row$status,
                        error.value = current.output$result.summary.row$error,
                        file.value = current.metadata$file,
                        method.label = "SOM") #stop immediately if SOM failed
}


## Combine SOM results
result.table <- do.call(rbind, result.summary.list) #combine result summary rows
rownames(result.table) <- NULL
optim.k.result.table <- do.call(rbind, optim.k.summary.list) #combine optim_k_summary rows
rownames(optim.k.result.table) <- NULL


## Calculate Fst and append to SOM result tables
Fst.output <- calculate.and.append.Fst(vcf.metadata.table = vcf.metadata.table,
                                       result.table = result.table,
                                       optim.k.result.table = optim.k.result.table,
                                       intermediate.results.directory = intermediate.results.directory,
                                       override = override,
                                       sample.n.per.population = analysis.sample.n.per.population,
                                       sampling.seed = analysis.sample.random.seed) #calculate and append Fst
result.table <- Fst.output$result.table #update result table with Fst
optim.k.result.table <- Fst.output$optim.k.result.table #update optim_k result table with Fst
fst.table <- Fst.output$fst.table #extract Fst table


## Calculate de novo k-means/BIC and append to result tables
deNovo.kmeans.output <- calculate.and.append.deNovo.kmeans(vcf.metadata.table = vcf.metadata.table,
                                                           result.table = result.table,
                                                           optim.k.result.table = optim.k.result.table,
                                                           intermediate.results.directory = intermediate.results.directory,
                                                           override = override,
                                                           BIC.thresh = deNovo.kmeans.BIC.thresh,
                                                           n.iter = deNovo.kmeans.n.iter,
                                                           n.start = deNovo.kmeans.n.start,
                                                           max.n.clust = deNovo.kmeans.max.n.clust,
                                                           max.proportion.PCs = deNovo.kmeans.max.proportion.PCs,
                                                           center = deNovo.kmeans.center,
                                                           scale = deNovo.kmeans.scale,
                                                           sample.n.per.population = analysis.sample.n.per.population,
                                                           sampling.seed = analysis.sample.random.seed) #calculate de novo k-means/BIC
result.table <- deNovo.kmeans.output$result.table #update result table with de novo k-means/BIC
optim.k.result.table <- deNovo.kmeans.output$optim.k.result.table #update optim_k result table with de novo k-means/BIC
deNovo.kmeans.table <- deNovo.kmeans.output$deNovo.kmeans.table #extract de novo k-means/BIC table


## Calculate sNMF and append to result tables
sNMF.output <- calculate.and.append.sNMF(vcf.metadata.table = vcf.metadata.table,
                                         result.table = result.table,
                                         optim.k.result.table = optim.k.result.table,
                                         intermediate.results.directory = intermediate.results.directory,
                                         override = override,
                                         K.values = sNMF.K.values,
                                         repetitions = sNMF.repetitions,
                                         CPU = 1,
                                         ploidy = sNMF.ploidy,
                                         cross.entropy.thresh = sNMF.cross.entropy.thresh,
                                         seed = sNMF.seed,
                                         sample.n.per.population = analysis.sample.n.per.population,
                                         sampling.seed = analysis.sample.random.seed) #calculate sNMF
result.table <- sNMF.output$result.table #update result table with sNMF
optim.k.result.table <- sNMF.output$optim.k.result.table #update optim_k result table with sNMF
sNMF.table <- sNMF.output$sNMF.table #extract sNMF table


## Save result tables
write.csv(result.table, file.path(results.directory, "fastsimcoal2_SOM_results.csv"), row.names = FALSE) #save main SOM result summary
write.csv(optim.k.result.table, file.path(results.directory, "fastsimcoal2_SOM_optim_k_summary.csv"), row.names = FALSE) #save full optim_k_summary table
write.csv(fst.table, file.path(results.directory, "fastsimcoal2_Weir_Cockerham_Fst.csv"), row.names = FALSE) #save Fst table
write.csv(deNovo.kmeans.table, file.path(results.directory, "fastsimcoal2_de_novo_kmeans_BIC.csv"), row.names = FALSE) #save de novo k-means/BIC table
write.csv(sNMF.table, file.path(results.directory, "fastsimcoal2_sNMF_results.csv"), row.names = FALSE) #save sNMF table


## Inspect results
head(result.table)
optim.k.result.table
head(fst.table)
head(deNovo.kmeans.table)
head(sNMF.table)
table(result.table$status)
table(result.table$fst.status)
table(result.table$deNovo.kmeans.status)
table(result.table$deNovo.kmeans.best.k, useNA = "ifany")
table(result.table$sNMF.status)


## Plot fitted K2 proportions for SOM, de novo k-means/BIC, and sNMF
plot.result.table <- result.table[result.table$status == "ok", ] #subset successful rows
plot.result.table <- plot.result.table[order(plot.result.table$mig, plot.result.table$tdiv), ] #order rows for plotting
if (nrow(plot.result.table) == 0) stop("No successful results available for fitted plotting")
mig.level.table <- unique(plot.result.table[order(plot.result.table$mig), c("mig.tag", "mig")]) #define unique ordered migration-rate levels
mig.level.table$mig.tag <- as.character(mig.level.table$mig.tag) #force migration-rate labels to discrete character values
plot.result.table$mig.tag <- factor(as.character(plot.result.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels by numeric migration rate
mig.colors <- setNames(viridisLite::viridis(nrow(mig.level.table), option = "D"), mig.level.table$mig.tag) #define discrete ordered colors


## Create SOM binomial count table
SOM.optim.k.count.table <- optim.k.result.table[!is.na(optim.k.result.table$Count) & optim.k.result.table$file %in% plot.result.table$file, ] #keep rows with SOM replicate counts
SOM.optim.k.count.table$Count <- as.numeric(as.character(SOM.optim.k.count.table$Count)) #force counts to numeric
SOM.total.count.table <- stats::aggregate(Count ~ file + mig.tag + mig + tdiv, data = SOM.optim.k.count.table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE)) #total retained SOM replicates per VCF
colnames(SOM.total.count.table)[colnames(SOM.total.count.table) == "Count"] <- "n.total"
SOM.k2.count.table <- SOM.optim.k.count.table[SOM.optim.k.count.table$k.label == "k2", ] #keep k2 count rows
if (nrow(SOM.k2.count.table) > 0) {
  SOM.k2.count.table <- stats::aggregate(Count ~ file + mig.tag + mig + tdiv, data = SOM.k2.count.table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE)) #number of retained SOM replicates choosing k = 2
  colnames(SOM.k2.count.table)[colnames(SOM.k2.count.table) == "Count"] <- "n.k2"
} else {
  SOM.k2.count.table <- SOM.total.count.table[, c("file", "mig.tag", "mig", "tdiv")] #create empty k2 count table if no k2 rows exist
  SOM.k2.count.table$n.k2 <- 0
}
SOM.binomial.table <- merge(SOM.total.count.table, SOM.k2.count.table, by = c("file", "mig.tag", "mig", "tdiv"), all.x = TRUE) #merge total and k2 counts
SOM.binomial.table$n.k2[is.na(SOM.binomial.table$n.k2)] <- 0 #set missing k2 rows to zero
SOM.binomial.table$n.not.k2 <- SOM.binomial.table$n.total - SOM.binomial.table$n.k2 #number not choosing k = 2
SOM.binomial.table$proportion.k2 <- SOM.binomial.table$n.k2 / SOM.binomial.table$n.total #calculate K2 proportion
SOM.binomial.table$method <- "SOM" #add method label
SOM.binomial.table <- SOM.binomial.table[, c("method", "file", "mig.tag", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")] #order columns


## Create de novo k-means/BIC binomial table
deNovo.kmeans.binomial.source.table <- result.table[result.table$status == "ok" & !is.na(result.table$deNovo.kmeans.best.k), ] #subset valid de novo k-means/BIC rows
if ("deNovo.kmeans.status" %in% colnames(deNovo.kmeans.binomial.source.table)) deNovo.kmeans.binomial.source.table <- deNovo.kmeans.binomial.source.table[deNovo.kmeans.binomial.source.table$deNovo.kmeans.status == "ok", ] #keep successful de novo k-means/BIC rows
deNovo.kmeans.best.k <- as.integer(as.character(deNovo.kmeans.binomial.source.table$deNovo.kmeans.best.k)) #extract de novo k-means/BIC best K
deNovo.kmeans.binomial.table <- data.frame(method = "deNovo.kmeans.BIC",
                                           file = deNovo.kmeans.binomial.source.table$file,
                                           mig.tag = deNovo.kmeans.binomial.source.table$mig.tag,
                                           mig = deNovo.kmeans.binomial.source.table$mig,
                                           tdiv = deNovo.kmeans.binomial.source.table$tdiv,
                                           n.k2 = as.integer(deNovo.kmeans.best.k == 2L),
                                           n.not.k2 = as.integer(deNovo.kmeans.best.k != 2L),
                                           proportion.k2 = as.numeric(deNovo.kmeans.best.k == 2L),
                                           stringsAsFactors = FALSE) #create binary de novo k-means/BIC table


## Create sNMF binomial table
sNMF.binomial.source.table <- result.table[result.table$status == "ok" & !is.na(result.table$sNMF.best.k), ] #subset valid sNMF rows
if ("sNMF.status" %in% colnames(sNMF.binomial.source.table)) sNMF.binomial.source.table <- sNMF.binomial.source.table[sNMF.binomial.source.table$sNMF.status == "ok", ] #keep successful sNMF rows
sNMF.best.k <- as.integer(as.character(sNMF.binomial.source.table$sNMF.best.k)) #extract sNMF best K
sNMF.binomial.table <- data.frame(method = "sNMF",
                                  file = sNMF.binomial.source.table$file,
                                  mig.tag = sNMF.binomial.source.table$mig.tag,
                                  mig = sNMF.binomial.source.table$mig,
                                  tdiv = sNMF.binomial.source.table$tdiv,
                                  n.k2 = as.integer(sNMF.best.k == 2L),
                                  n.not.k2 = as.integer(sNMF.best.k != 2L),
                                  proportion.k2 = as.numeric(sNMF.best.k == 2L),
                                  stringsAsFactors = FALSE) #create binary sNMF table


## Create STRUCTURE binomial table
STRUCTURE.results <- create.STRUCTURE.mean.lnprob.binomial.table(analysis.set.name = analysis.set.name,
                                                                 structure.loglik.directory = STRUCTURE.loglik.directory,
                                                                 results.directory = results.directory,
                                                                 expected.structure.replicates = expected.structure.replicates,
                                                                 STRUCTURE.mean.lnprob.delta.threshold = STRUCTURE.mean.lnprob.delta.threshold) #create STRUCTURE mean-lnP best-K table
STRUCTURE.binomial.table <- STRUCTURE.results$STRUCTURE.binomial.table #extract STRUCTURE binomial table
STRUCTURE.binomial.table$mig.tag <- vapply(STRUCTURE.binomial.table$mig, function(current.mig) {
  matching.mig.index <- which(abs(mig.level.table$mig - current.mig) < sqrt(.Machine$double.eps))
  if (length(matching.mig.index) != 1) return(NA_character_)
  return(as.character(mig.level.table$mig.tag[matching.mig.index]))
}, character(1)) #standardize STRUCTURE migration labels to match VCF migration labels
if (any(is.na(STRUCTURE.binomial.table$mig.tag))) stop("Some STRUCTURE migration rates could not be matched to VCF migration labels")
STRUCTURE.mean.lnprob.by.vcf.table <- STRUCTURE.results$STRUCTURE.mean.lnprob.by.vcf.table #extract STRUCTURE VCF-level table
STRUCTURE.rep.comparison.table <- STRUCTURE.results$STRUCTURE.rep.comparison.table #extract STRUCTURE replicate comparison table


## Combine method-specific binomial tables
fitted.binomial.input.table <- rbind(SOM.binomial.table, deNovo.kmeans.binomial.table, sNMF.binomial.table, STRUCTURE.binomial.table) #combine all methods
method.levels <- c("SOM", "deNovo.kmeans.BIC", "sNMF", "STRUCTURE") #define method order
fitted.binomial.input.table$method <- factor(as.character(fitted.binomial.input.table$method), levels = method.levels) #order methods
fitted.binomial.input.table$mig.tag <- factor(as.character(fitted.binomial.input.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels


## Fit binomial GLMs separately for each method and migration rate
fitted.GLM.prediction.list <- list() #initialize prediction list
fitted.GLM.summary.list <- list() #initialize summary list
fitted.GLM.list.index <- 1 #initialize list index
for (current.method in method.levels) {
  for (current.mig.tag in mig.level.table$mig.tag) {
    current.binomial.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == current.method & fitted.binomial.input.table$mig.tag == current.mig.tag, ] #subset one method and migration rate
    current.binomial.table <- current.binomial.table[is.finite(current.binomial.table$tdiv), ] #keep rows with finite divergence time
    if (nrow(current.binomial.table) == 0) next
    current.binomial.table <- current.binomial.table[order(current.binomial.table$tdiv), ] #order rows by divergence time
    prediction.tdiv <- seq(min(current.binomial.table$tdiv), max(current.binomial.table$tdiv), length.out = plot.fitted.prediction.n) #create prediction grid
    current.warnings <- character(0) #store GLM warnings
    if (all(current.binomial.table$n.k2 == 0)) {
      current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = 0, stringsAsFactors = FALSE) #flat zero prediction when K2 is never selected
      current.summary.row <- data.frame(method = current.method,
                                        mig.tag = as.character(current.mig.tag),
                                        mig = current.binomial.table$mig[1],
                                        n.vcf = nrow(current.binomial.table),
                                        n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                        n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                        min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                        n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                        n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                        all.zero = TRUE,
                                        all.one = FALSE,
                                        glm.fitted = FALSE,
                                        glm.converged = NA,
                                        glm.intercept = NA_real_,
                                        glm.tdiv.slope = NA_real_,
                                        tdiv.at.threshold.k2 = NA_real_,
                                        tdiv.at.threshold.k2.within.range = FALSE,
                                        glm.warning = NA_character_,
                                        glm.error = NA_character_,
                                        stringsAsFactors = FALSE) #store summary row
    } else if (all(current.binomial.table$n.not.k2 == 0)) {
      current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = 1, stringsAsFactors = FALSE) #flat one prediction when K2 is always selected
      current.summary.row <- data.frame(method = current.method,
                                        mig.tag = as.character(current.mig.tag),
                                        mig = current.binomial.table$mig[1],
                                        n.vcf = nrow(current.binomial.table),
                                        n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                        n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                        min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                        n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                        n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                        all.zero = FALSE,
                                        all.one = TRUE,
                                        glm.fitted = FALSE,
                                        glm.converged = NA,
                                        glm.intercept = NA_real_,
                                        glm.tdiv.slope = NA_real_,
                                        tdiv.at.threshold.k2 = NA_real_,
                                        tdiv.at.threshold.k2.within.range = FALSE,
                                        glm.warning = NA_character_,
                                        glm.error = NA_character_,
                                        stringsAsFactors = FALSE) #store summary row
    } else {
      current.fit <- withCallingHandlers(tryCatch(stats::glm(cbind(n.k2, n.not.k2) ~ tdiv, data = current.binomial.table, family = stats::binomial(link = "logit"), control = stats::glm.control(maxit = 100)), error = function(e) e),
                                         warning = function(w) {
                                           current.warnings <<- c(current.warnings, conditionMessage(w))
                                           invokeRestart("muffleWarning")
                                         }) #fit count-based binomial GLM and store warnings
      fit.ok <- inherits(current.fit, "glm") #check whether GLM fit succeeded
      if (fit.ok) {
        current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = stats::predict(current.fit, newdata = data.frame(tdiv = prediction.tdiv), type = "response"), stringsAsFactors = FALSE) #predict fitted K2 probability
        current.coefficients <- stats::coef(current.fit) #extract coefficients
        current.glm.intercept <- unname(current.coefficients["(Intercept)"]) #extract intercept
        current.glm.tdiv.slope <- unname(current.coefficients["tdiv"]) #extract tdiv slope
        if (is.finite(current.glm.intercept) && is.finite(current.glm.tdiv.slope) && current.glm.tdiv.slope != 0) {
          current.tdiv.at.threshold.k2 <- (stats::qlogis(plot.k2.threshold.proportion) - current.glm.intercept) / current.glm.tdiv.slope #solve logit(p) for threshold
        } else {
          current.tdiv.at.threshold.k2 <- NA_real_
        }
        current.tdiv.at.threshold.k2.within.range <- is.finite(current.tdiv.at.threshold.k2) && current.tdiv.at.threshold.k2 >= min(current.binomial.table$tdiv, na.rm = TRUE) && current.tdiv.at.threshold.k2 <= max(current.binomial.table$tdiv, na.rm = TRUE) #check whether threshold is inside simulated range
        current.summary.row <- data.frame(method = current.method,
                                          mig.tag = as.character(current.mig.tag),
                                          mig = current.binomial.table$mig[1],
                                          n.vcf = nrow(current.binomial.table),
                                          n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                          n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                          min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                          n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                          n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                          all.zero = all(current.binomial.table$proportion.k2 == 0),
                                          all.one = all(current.binomial.table$proportion.k2 == 1),
                                          glm.fitted = TRUE,
                                          glm.converged = isTRUE(current.fit$converged),
                                          glm.intercept = current.glm.intercept,
                                          glm.tdiv.slope = current.glm.tdiv.slope,
                                          tdiv.at.threshold.k2 = current.tdiv.at.threshold.k2,
                                          tdiv.at.threshold.k2.within.range = current.tdiv.at.threshold.k2.within.range,
                                          glm.warning = if (length(current.warnings) == 0) NA_character_ else paste(unique(current.warnings), collapse = " | "),
                                          glm.error = NA_character_,
                                          stringsAsFactors = FALSE) #store summary row
      } else {
        current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = NA_real_, stringsAsFactors = FALSE) #return missing predictions if GLM failed
        current.summary.row <- data.frame(method = current.method,
                                          mig.tag = as.character(current.mig.tag),
                                          mig = current.binomial.table$mig[1],
                                          n.vcf = nrow(current.binomial.table),
                                          n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                          n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                          min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                          n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                          n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                          all.zero = all(current.binomial.table$proportion.k2 == 0),
                                          all.one = all(current.binomial.table$proportion.k2 == 1),
                                          glm.fitted = FALSE,
                                          glm.converged = FALSE,
                                          glm.intercept = NA_real_,
                                          glm.tdiv.slope = NA_real_,
                                          tdiv.at.threshold.k2 = NA_real_,
                                          tdiv.at.threshold.k2.within.range = FALSE,
                                          glm.warning = if (length(current.warnings) == 0) NA_character_ else paste(unique(current.warnings), collapse = " | "),
                                          glm.error = conditionMessage(current.fit),
                                          stringsAsFactors = FALSE) #store summary row
      }
    }
    fitted.GLM.prediction.list[[fitted.GLM.list.index]] <- current.prediction.table #store prediction table
    fitted.GLM.summary.list[[fitted.GLM.list.index]] <- current.summary.row #store summary row
    fitted.GLM.list.index <- fitted.GLM.list.index + 1 #advance list index
  }
}


## Combine fitted GLM outputs
fitted.GLM.prediction.table <- do.call(rbind, fitted.GLM.prediction.list) #combine prediction tables
fitted.GLM.summary.table <- do.call(rbind, fitted.GLM.summary.list) #combine summary rows
fitted.GLM.prediction.table$method <- factor(as.character(fitted.GLM.prediction.table$method), levels = method.levels) #order methods
fitted.GLM.prediction.table$mig.tag <- factor(as.character(fitted.GLM.prediction.table$mig.tag), levels = mig.level.table$mig.tag) #order migration labels
fitted.GLM.summary.table$method <- factor(as.character(fitted.GLM.summary.table$method), levels = method.levels) #order methods
fitted.GLM.summary.table$mig.tag <- factor(as.character(fitted.GLM.summary.table$mig.tag), levels = mig.level.table$mig.tag) #order migration labels
fitted.GLM.summary.table <- fitted.GLM.summary.table[order(fitted.GLM.summary.table$method, fitted.GLM.summary.table$mig), ] #order summary table
rownames(fitted.GLM.summary.table) <- NULL
fitted.GLM.threshold.table <- fitted.GLM.summary.table[fitted.GLM.summary.table$tdiv.at.threshold.k2.within.range & is.finite(fitted.GLM.summary.table$tdiv.at.threshold.k2), ] #keep valid threshold rows


## Plot SOM results
SOM.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "SOM", ] #subset SOM observations
SOM.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "SOM" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset SOM fitted predictions
SOM.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "SOM", ] #subset SOM threshold rows
SOM.k2.fitted.plot <- ggplot2::ggplot(SOM.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = SOM.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = SOM.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(SOM.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "SOM_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = SOM.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot Weir and Cockerham Fst results
fst.plot.table <- plot.result.table[is.finite(plot.result.table$fst), ] #subset rows with finite Fst values
fst.plot.table <- fst.plot.table[order(fst.plot.table$mig, fst.plot.table$tdiv), ] #order rows for plotting
fst.plot.table$mig.tag <- factor(as.character(fst.plot.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels
SOM.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "SOM" & fitted.GLM.threshold.table$tdiv.at.threshold.k2.within.range & is.finite(fitted.GLM.threshold.table$tdiv.at.threshold.k2), ] #subset valid SOM threshold rows
SOM.k2.threshold.plot.table$mig.tag <- factor(as.character(SOM.k2.threshold.plot.table$mig.tag), levels = mig.level.table$mig.tag) #order SOM threshold migration-rate labels
fst.plot <- ggplot2::ggplot(fst.plot.table, ggplot2::aes(x = tdiv, y = fst, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = Fst_plot_point_size, alpha = Fst_plot_point_alpha) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, linewidth = Fst_plot_line_width, se = FALSE) +
  ggplot2::geom_vline(data = SOM.k2.threshold.plot.table,
                      ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag),
                      linewidth = plot.threshold.line.width,
                      linetype = plot.threshold.line.type,
                      show.legend = FALSE,
                      inherit.aes = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Weir and Cockerham Fst", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(fst.plot)
ggplot2::ggsave(file.path(results.directory, "Weir_Cockerham_Fst_by_tdiv_and_mig.svg"),
                plot = fst.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot de novo k-means/BIC results
deNovo.kmeans.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "deNovo.kmeans.BIC", ] #subset de novo k-means/BIC observations
deNovo.kmeans.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "deNovo.kmeans.BIC" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset de novo k-means/BIC fitted predictions
deNovo.kmeans.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "deNovo.kmeans.BIC", ] #subset de novo k-means/BIC threshold rows
deNovo.kmeans.k2.fitted.plot <- ggplot2::ggplot(deNovo.kmeans.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = deNovo.kmeans.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = deNovo.kmeans.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(deNovo.kmeans.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "deNovo_kmeans_BIC_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = deNovo.kmeans.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot sNMF K2 results
sNMF.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "sNMF", ] #subset sNMF observations
sNMF.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "sNMF" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset sNMF fitted predictions
sNMF.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "sNMF", ] #subset sNMF threshold rows
sNMF.k2.fitted.plot <- ggplot2::ggplot(sNMF.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = sNMF.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = sNMF.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(sNMF.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "sNMF_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = sNMF.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot STRUCTURE K2 results
STRUCTURE.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "STRUCTURE", ] #subset STRUCTURE observations
STRUCTURE.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "STRUCTURE" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset STRUCTURE fitted predictions
STRUCTURE.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "STRUCTURE", ] #subset STRUCTURE threshold rows
STRUCTURE.k2.fitted.plot <- ggplot2::ggplot(STRUCTURE.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = STRUCTURE.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = STRUCTURE.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(STRUCTURE.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "STRUCTURE_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = STRUCTURE.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Inspect fitted 50% threshold results
fitted.GLM.summary.table
fitted.GLM.threshold.table



## Create PCA plot for one VCF file

## Show available migration settings
available.migration.table <- unique(vcf.metadata.table[, c("mig.tag", "mig")]) #show all available migration settings including mig = 0
available.migration.table <- available.migration.table[order(available.migration.table$mig), ]
rownames(available.migration.table) <- NULL
print(available.migration.table)
if (!(migration.vcf.mig.for.PCA %in% unique(vcf.metadata.table$mig))) stop(paste("migration.vcf.mig.for.PCA must be one of:", paste(sort(unique(vcf.metadata.table$mig)), collapse = ", ")))


## Select file number N among VCF files for selected migration rate
migration.vcf.metadata.table <- vcf.metadata.table[vcf.metadata.table$mig == migration.vcf.mig.for.PCA, ] #keep selected VCF files only
migration.vcf.metadata.table <- migration.vcf.metadata.table[order(migration.vcf.metadata.table$tdiv), ] #order selected files by divergence time
rownames(migration.vcf.metadata.table) <- NULL
if (nrow(migration.vcf.metadata.table) == 0) stop(paste("No VCF files found for migration rate", migration.vcf.mig.for.PCA))
if (migration.vcf.file.number.for.PCA < 1 || migration.vcf.file.number.for.PCA > nrow(migration.vcf.metadata.table)) stop(paste("migration.vcf.file.number.for.PCA must be between 1 and", nrow(migration.vcf.metadata.table), "for migration rate", migration.vcf.mig.for.PCA))
current.pca.metadata <- migration.vcf.metadata.table[migration.vcf.file.number.for.PCA, ] #select VCF file N for selected migration rate
cat("\nCreating PCA for VCF file", migration.vcf.file.number.for.PCA, "of", nrow(migration.vcf.metadata.table), "for migration rate", migration.vcf.mig.for.PCA, ":", current.pca.metadata$file, "\n")
cat("Migration-index tag:", current.pca.metadata$mig.index, "\n")
cat("Divergence time:", current.pca.metadata$tdiv, "generations\n")
cat("Migration tag:", current.pca.metadata$mig.tag, "\n")
cat("Migration rate:", current.pca.metadata$mig, "\n")


## Read VCF and process SNP matrix
pca.genind.object.raw <- read.fastsimcoal2.genind.for.analysis(vcf.file.path = current.pca.metadata$full.path,
                                                               n.per.population = analysis.sample.n.per.population,
                                                               sampling.group = current.pca.metadata$simulation.group,
                                                               sampling.seed = analysis.sample.random.seed) #read selected VCF, assign populations, and optionally subset individuals per population


## Calculate Weir and Cockerham Fst before PCA processing
pca.fst.value <- calculate.overall.Fst(genind.object = pca.genind.object.raw) #calculate overall Weir and Cockerham Fst
cat("Weir and Cockerham Fst between pop1 and pop2:", round(pca.fst.value, 4), "\n")


## Process SNP matrix for PCA
pca.snp.matrix <- process.SNP.data.SOM(genind.input = pca.genind.object.raw, verbose = FALSE) #process SNP data for PCA
pca.snp.matrix <- as.matrix(pca.snp.matrix) #coerce to matrix
storage.mode(pca.snp.matrix) <- "numeric" #ensure numeric matrix
if (is.null(rownames(pca.snp.matrix))) rownames(pca.snp.matrix) <- adegenet::indNames(pca.genind.object.raw)


## Remove invalid loci
valid.pca.loci <- apply(pca.snp.matrix, 2, function(locus.values) {
  non.missing.values <- locus.values[!is.na(locus.values)]
  length(non.missing.values) > 1 && stats::var(non.missing.values) > 0
}) #identify loci with nonzero variance
pca.snp.matrix <- pca.snp.matrix[, valid.pca.loci, drop = FALSE] #keep valid loci only
if (ncol(pca.snp.matrix) < 2) stop("Fewer than two valid loci available for PCA")


## Impute remaining missing values for PCA
if (any(is.na(pca.snp.matrix))) {
  pca.snp.matrix <- apply(pca.snp.matrix, 2, function(locus.values) {
    locus.values[is.na(locus.values)] <- mean(locus.values, na.rm = TRUE)
    locus.values
  }) #replace remaining missing values with locus means
  pca.snp.matrix <- as.matrix(pca.snp.matrix)
}


## Run PCA and create plot table
pca.object <- stats::prcomp(pca.snp.matrix, center = TRUE, scale. = FALSE) #run PCA
pca.variance <- pca.object$sdev^2 / sum(pca.object$sdev^2) #calculate variance explained
pca.plot.table <- data.frame(individual = rownames(pca.object$x),
                             PC1 = pca.object$x[, 1],
                             PC2 = pca.object$x[, 2],
                             population = as.character(adegenet::pop(pca.genind.object.raw))[match(rownames(pca.object$x), adegenet::indNames(pca.genind.object.raw))],
                             file = current.pca.metadata$file,
                             mig.index = current.pca.metadata$mig.index,
                             tdiv = current.pca.metadata$tdiv,
                             mig.tag = current.pca.metadata$mig.tag,
                             mig = current.pca.metadata$mig,
                             fst.pop1.pop2 = pca.fst.value,
                             stringsAsFactors = FALSE) #create PCA plotting table


## Plot PCA with population groupings
pca.plot <- ggplot2::ggplot(pca.plot.table, ggplot2::aes(x = PC1, y = PC2, color = population, group = population)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(title = paste0("PCA - ", current.pca.metadata$file, " | tdiv = ", round(current.pca.metadata$tdiv), " generations | mig = ", current.pca.metadata$mig, " | Fst = ", round(pca.fst.value, 4)),
                x = paste0("PC1 (", round(pca.variance[1] * 100, 2), "%)"),
                y = paste0("PC2 (", round(pca.variance[2] * 100, 2), "%)"),
                color = "Population") +
  ggplot2::theme_classic()
pca.plot


## Print selected PCA summary
cat("Divergence:", round(current.pca.metadata$tdiv / 1000000, 2), "million generations\n")
cat("Migration rate:", current.pca.metadata$mig, "\n")
cat("Weir and Cockerham Fst:", round(pca.fst.value, 4), "\n")




#### Run workflow for VCF files - symmetric_8 #################################
SOM_grid_multiplier <- 4
analysis.sample.n.per.population <- analysis.sample.n.per.population.symmetric_8 #set individuals per population for symmetric_8
analysis.set.name <- paste0("symmetric_", analysis.sample.n.per.population) #define analysis set name from sample size
results.directory <- file.path(results.root.directory, analysis.set.name) #define output directory for this analysis set
dir.create(results.directory, recursive = TRUE, showWarnings = FALSE) #create output directory if needed
intermediate.results.directory <- file.path(results.directory, "intermediate_results") #define intermediate results directory for this analysis set
dir.create(intermediate.results.directory, recursive = TRUE, showWarnings = FALSE) #create intermediate results directory if needed
cat("\n\nRunning analysis set:", analysis.set.name, "\n")
cat("Individuals per population:", analysis.sample.n.per.population, "\n")
cat("Results directory:", results.directory, "\n")
cat("Random sampling base seed:", analysis.sample.random.seed, "\n")


## List VCF files
vcf.file.pattern <- "^sim[0-9]+_tdiv[0-9.]+_mig.+\\.vcf$" #define VCF filename pattern
vcf.file.paths <- list.files(vcf.directory, pattern = vcf.file.pattern, full.names = TRUE) #list VCF files only
if (length(vcf.file.paths) == 0) stop("No VCF files found")
vcf.file.names <- basename(vcf.file.paths) #extract file names


## Parse filename metadata
vcf.metadata.list <- lapply(vcf.file.names, parse.vcf.filename.metadata) #parse metadata
vcf.metadata.table <- do.call(rbind, vcf.metadata.list) #combine metadata rows
vcf.metadata.table$full.path <- vcf.file.paths #add full path


## Order metadata table
vcf.metadata.table <- vcf.metadata.table[order(vcf.metadata.table$mig.index,
                                               vcf.metadata.table$tdiv,
                                               vcf.metadata.table$mig), ] #order rows
rownames(vcf.metadata.table) <- NULL
vcf.metadata.table$simulation.group <- match(vcf.metadata.table$tdiv, sort(unique(vcf.metadata.table$tdiv))) #define matched sampling group across migration rates for each divergence time
head(vcf.metadata.table) #show metadata for files


## Create empty lists to store results
result.summary.list <- vector("list", nrow(vcf.metadata.table)) #initialize result summary list
optim.k.summary.list <- vector("list", nrow(vcf.metadata.table)) #initialize optim_k_summary list


## Loop over VCF files
for (vcf.index in seq_len(nrow(vcf.metadata.table))) {
  current.metadata <- vcf.metadata.table[vcf.index, ] #extract current metadata
  cat("\nProcessing", vcf.index, "of", nrow(vcf.metadata.table), ":", current.metadata$file, "\n")
  if (!override && check.vcf.intermediate.results.available(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory)) {
    cat("Intermediate results already found for", current.metadata$file, "- loading saved results\n")
    current.output <- load.vcf.intermediate.results(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory) #load saved per-VCF results
  } else {
    current.output <- tryCatch({
      run.SOM.workflow.for.one.vcf(vcf.file.path = current.metadata$full.path,
                                   vcf.file.metadata = current.metadata,
                                   sample.n.per.population = analysis.sample.n.per.population,
                                   sampling.group = current.metadata$simulation.group,
                                   sampling.seed = analysis.sample.random.seed)
    }, error = function(e) {
      result.summary.row <- data.frame(file = current.metadata$file,
                                       mig.index = current.metadata$mig.index,
                                       tdiv = current.metadata$tdiv,
                                       mig.tag = current.metadata$mig.tag,
                                       mig = current.metadata$mig,
                                       n.ind.raw = NA_integer_,
                                       n.loc.raw = NA_integer_,
                                       n.samples.processed = NA_integer_,
                                       n.loci.processed = NA_integer_,
                                       proportion.k2 = NA_real_,
                                       status = "error",
                                       error = conditionMessage(e),
                                       stringsAsFactors = FALSE)
      optim.k.summary.table <- data.frame(file = current.metadata$file,
                                          mig.index = current.metadata$mig.index,
                                          tdiv = current.metadata$tdiv,
                                          mig.tag = current.metadata$mig.tag,
                                          mig = current.metadata$mig,
                                          n.ind.raw = NA_integer_,
                                          n.loc.raw = NA_integer_,
                                          n.samples.processed = NA_integer_,
                                          n.loci.processed = NA_integer_,
                                          status = "error",
                                          error = conditionMessage(e),
                                          Count = NA_real_,
                                          Proportion = NA_real_,
                                          k.label = NA_character_,
                                          optim.k.row = NA_integer_,
                                          stringsAsFactors = FALSE)
      list(result.summary.row = result.summary.row, optim.k.summary.table = optim.k.summary.table)
    })
    save.vcf.intermediate.results(result.summary.row = current.output$result.summary.row,
                                  optim.k.summary.table = current.output$optim.k.summary.table,
                                  vcf.file.metadata = current.metadata,
                                  intermediate.results.directory = intermediate.results.directory) #save per-VCF results immediately
  }
  result.summary.list[[vcf.index]] <- current.output$result.summary.row #store result summary row
  optim.k.summary.list[[vcf.index]] <- current.output$optim.k.summary.table #store optim_k_summary table
  stop.if.status.not.ok(status.value = current.output$result.summary.row$status,
                        error.value = current.output$result.summary.row$error,
                        file.value = current.metadata$file,
                        method.label = "SOM") #stop immediately if SOM failed
}


## Combine SOM results
result.table <- do.call(rbind, result.summary.list) #combine result summary rows
rownames(result.table) <- NULL
optim.k.result.table <- do.call(rbind, optim.k.summary.list) #combine optim_k_summary rows
rownames(optim.k.result.table) <- NULL


## Calculate Fst and append to SOM result tables
Fst.output <- calculate.and.append.Fst(vcf.metadata.table = vcf.metadata.table,
                                       result.table = result.table,
                                       optim.k.result.table = optim.k.result.table,
                                       intermediate.results.directory = intermediate.results.directory,
                                       override = override,
                                       sample.n.per.population = analysis.sample.n.per.population,
                                       sampling.seed = analysis.sample.random.seed) #calculate and append Fst
result.table <- Fst.output$result.table #update result table with Fst
optim.k.result.table <- Fst.output$optim.k.result.table #update optim_k result table with Fst
fst.table <- Fst.output$fst.table #extract Fst table


## Calculate de novo k-means/BIC and append to result tables
deNovo.kmeans.output <- calculate.and.append.deNovo.kmeans(vcf.metadata.table = vcf.metadata.table,
                                                           result.table = result.table,
                                                           optim.k.result.table = optim.k.result.table,
                                                           intermediate.results.directory = intermediate.results.directory,
                                                           override = override,
                                                           BIC.thresh = deNovo.kmeans.BIC.thresh,
                                                           n.iter = deNovo.kmeans.n.iter,
                                                           n.start = deNovo.kmeans.n.start,
                                                           max.n.clust = deNovo.kmeans.max.n.clust,
                                                           max.proportion.PCs = deNovo.kmeans.max.proportion.PCs,
                                                           center = deNovo.kmeans.center,
                                                           scale = deNovo.kmeans.scale,
                                                           sample.n.per.population = analysis.sample.n.per.population,
                                                           sampling.seed = analysis.sample.random.seed) #calculate de novo k-means/BIC
result.table <- deNovo.kmeans.output$result.table #update result table with de novo k-means/BIC
optim.k.result.table <- deNovo.kmeans.output$optim.k.result.table #update optim_k result table with de novo k-means/BIC
deNovo.kmeans.table <- deNovo.kmeans.output$deNovo.kmeans.table #extract de novo k-means/BIC table


## Calculate sNMF and append to result tables
sNMF.output <- calculate.and.append.sNMF(vcf.metadata.table = vcf.metadata.table,
                                         result.table = result.table,
                                         optim.k.result.table = optim.k.result.table,
                                         intermediate.results.directory = intermediate.results.directory,
                                         override = override,
                                         K.values = sNMF.K.values,
                                         repetitions = sNMF.repetitions,
                                         CPU = 1,
                                         ploidy = sNMF.ploidy,
                                         cross.entropy.thresh = sNMF.cross.entropy.thresh,
                                         seed = sNMF.seed,
                                         sample.n.per.population = analysis.sample.n.per.population,
                                         sampling.seed = analysis.sample.random.seed) #calculate sNMF
result.table <- sNMF.output$result.table #update result table with sNMF
optim.k.result.table <- sNMF.output$optim.k.result.table #update optim_k result table with sNMF
sNMF.table <- sNMF.output$sNMF.table #extract sNMF table


## Save result tables
write.csv(result.table, file.path(results.directory, "fastsimcoal2_SOM_results.csv"), row.names = FALSE) #save main SOM result summary
write.csv(optim.k.result.table, file.path(results.directory, "fastsimcoal2_SOM_optim_k_summary.csv"), row.names = FALSE) #save full optim_k_summary table
write.csv(fst.table, file.path(results.directory, "fastsimcoal2_Weir_Cockerham_Fst.csv"), row.names = FALSE) #save Fst table
write.csv(deNovo.kmeans.table, file.path(results.directory, "fastsimcoal2_de_novo_kmeans_BIC.csv"), row.names = FALSE) #save de novo k-means/BIC table
write.csv(sNMF.table, file.path(results.directory, "fastsimcoal2_sNMF_results.csv"), row.names = FALSE) #save sNMF table


## Inspect results
head(result.table)
optim.k.result.table
head(fst.table)
head(deNovo.kmeans.table)
head(sNMF.table)
table(result.table$status)
table(result.table$fst.status)
table(result.table$deNovo.kmeans.status)
table(result.table$deNovo.kmeans.best.k, useNA = "ifany")
table(result.table$sNMF.status)


## Plot fitted K2 proportions for SOM, de novo k-means/BIC, and sNMF
plot.result.table <- result.table[result.table$status == "ok", ] #subset successful rows
plot.result.table <- plot.result.table[order(plot.result.table$mig, plot.result.table$tdiv), ] #order rows for plotting
if (nrow(plot.result.table) == 0) stop("No successful results available for fitted plotting")
mig.level.table <- unique(plot.result.table[order(plot.result.table$mig), c("mig.tag", "mig")]) #define unique ordered migration-rate levels
mig.level.table$mig.tag <- as.character(mig.level.table$mig.tag) #force migration-rate labels to discrete character values
plot.result.table$mig.tag <- factor(as.character(plot.result.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels by numeric migration rate
mig.colors <- setNames(viridisLite::viridis(nrow(mig.level.table), option = "D"), mig.level.table$mig.tag) #define discrete ordered colors


## Create SOM binomial count table
SOM.optim.k.count.table <- optim.k.result.table[!is.na(optim.k.result.table$Count) & optim.k.result.table$file %in% plot.result.table$file, ] #keep rows with SOM replicate counts
SOM.optim.k.count.table$Count <- as.numeric(as.character(SOM.optim.k.count.table$Count)) #force counts to numeric
SOM.total.count.table <- stats::aggregate(Count ~ file + mig.tag + mig + tdiv, data = SOM.optim.k.count.table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE)) #total retained SOM replicates per VCF
colnames(SOM.total.count.table)[colnames(SOM.total.count.table) == "Count"] <- "n.total"
SOM.k2.count.table <- SOM.optim.k.count.table[SOM.optim.k.count.table$k.label == "k2", ] #keep k2 count rows
if (nrow(SOM.k2.count.table) > 0) {
  SOM.k2.count.table <- stats::aggregate(Count ~ file + mig.tag + mig + tdiv, data = SOM.k2.count.table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE)) #number of retained SOM replicates choosing k = 2
  colnames(SOM.k2.count.table)[colnames(SOM.k2.count.table) == "Count"] <- "n.k2"
} else {
  SOM.k2.count.table <- SOM.total.count.table[, c("file", "mig.tag", "mig", "tdiv")] #create empty k2 count table if no k2 rows exist
  SOM.k2.count.table$n.k2 <- 0
}
SOM.binomial.table <- merge(SOM.total.count.table, SOM.k2.count.table, by = c("file", "mig.tag", "mig", "tdiv"), all.x = TRUE) #merge total and k2 counts
SOM.binomial.table$n.k2[is.na(SOM.binomial.table$n.k2)] <- 0 #set missing k2 rows to zero
SOM.binomial.table$n.not.k2 <- SOM.binomial.table$n.total - SOM.binomial.table$n.k2 #number not choosing k = 2
SOM.binomial.table$proportion.k2 <- SOM.binomial.table$n.k2 / SOM.binomial.table$n.total #calculate K2 proportion
SOM.binomial.table$method <- "SOM" #add method label
SOM.binomial.table <- SOM.binomial.table[, c("method", "file", "mig.tag", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")] #order columns


## Create de novo k-means/BIC binomial table
deNovo.kmeans.binomial.source.table <- result.table[result.table$status == "ok" & !is.na(result.table$deNovo.kmeans.best.k), ] #subset valid de novo k-means/BIC rows
if ("deNovo.kmeans.status" %in% colnames(deNovo.kmeans.binomial.source.table)) deNovo.kmeans.binomial.source.table <- deNovo.kmeans.binomial.source.table[deNovo.kmeans.binomial.source.table$deNovo.kmeans.status == "ok", ] #keep successful de novo k-means/BIC rows
deNovo.kmeans.best.k <- as.integer(as.character(deNovo.kmeans.binomial.source.table$deNovo.kmeans.best.k)) #extract de novo k-means/BIC best K
deNovo.kmeans.binomial.table <- data.frame(method = "deNovo.kmeans.BIC",
                                           file = deNovo.kmeans.binomial.source.table$file,
                                           mig.tag = deNovo.kmeans.binomial.source.table$mig.tag,
                                           mig = deNovo.kmeans.binomial.source.table$mig,
                                           tdiv = deNovo.kmeans.binomial.source.table$tdiv,
                                           n.k2 = as.integer(deNovo.kmeans.best.k == 2L),
                                           n.not.k2 = as.integer(deNovo.kmeans.best.k != 2L),
                                           proportion.k2 = as.numeric(deNovo.kmeans.best.k == 2L),
                                           stringsAsFactors = FALSE) #create binary de novo k-means/BIC table


## Create sNMF binomial table
sNMF.binomial.source.table <- result.table[result.table$status == "ok" & !is.na(result.table$sNMF.best.k), ] #subset valid sNMF rows
if ("sNMF.status" %in% colnames(sNMF.binomial.source.table)) sNMF.binomial.source.table <- sNMF.binomial.source.table[sNMF.binomial.source.table$sNMF.status == "ok", ] #keep successful sNMF rows
sNMF.best.k <- as.integer(as.character(sNMF.binomial.source.table$sNMF.best.k)) #extract sNMF best K
sNMF.binomial.table <- data.frame(method = "sNMF",
                                  file = sNMF.binomial.source.table$file,
                                  mig.tag = sNMF.binomial.source.table$mig.tag,
                                  mig = sNMF.binomial.source.table$mig,
                                  tdiv = sNMF.binomial.source.table$tdiv,
                                  n.k2 = as.integer(sNMF.best.k == 2L),
                                  n.not.k2 = as.integer(sNMF.best.k != 2L),
                                  proportion.k2 = as.numeric(sNMF.best.k == 2L),
                                  stringsAsFactors = FALSE) #create binary sNMF table


## Create STRUCTURE binomial table
STRUCTURE.results <- create.STRUCTURE.mean.lnprob.binomial.table(analysis.set.name = analysis.set.name,
                                                                 structure.loglik.directory = STRUCTURE.loglik.directory,
                                                                 results.directory = results.directory,
                                                                 expected.structure.replicates = expected.structure.replicates,
                                                                 STRUCTURE.mean.lnprob.delta.threshold = STRUCTURE.mean.lnprob.delta.threshold) #create STRUCTURE mean-lnP best-K table
STRUCTURE.binomial.table <- STRUCTURE.results$STRUCTURE.binomial.table #extract STRUCTURE binomial table
STRUCTURE.binomial.table$mig.tag <- vapply(STRUCTURE.binomial.table$mig, function(current.mig) {
  matching.mig.index <- which(abs(mig.level.table$mig - current.mig) < sqrt(.Machine$double.eps))
  if (length(matching.mig.index) != 1) return(NA_character_)
  return(as.character(mig.level.table$mig.tag[matching.mig.index]))
}, character(1)) #standardize STRUCTURE migration labels to match VCF migration labels
if (any(is.na(STRUCTURE.binomial.table$mig.tag))) stop("Some STRUCTURE migration rates could not be matched to VCF migration labels")
STRUCTURE.mean.lnprob.by.vcf.table <- STRUCTURE.results$STRUCTURE.mean.lnprob.by.vcf.table #extract STRUCTURE VCF-level table
STRUCTURE.rep.comparison.table <- STRUCTURE.results$STRUCTURE.rep.comparison.table #extract STRUCTURE replicate comparison table


## Combine method-specific binomial tables
fitted.binomial.input.table <- rbind(SOM.binomial.table, deNovo.kmeans.binomial.table, sNMF.binomial.table, STRUCTURE.binomial.table) #combine all methods
method.levels <- c("SOM", "deNovo.kmeans.BIC", "sNMF", "STRUCTURE") #define method order
fitted.binomial.input.table$method <- factor(as.character(fitted.binomial.input.table$method), levels = method.levels) #order methods
fitted.binomial.input.table$mig.tag <- factor(as.character(fitted.binomial.input.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels


## Fit binomial GLMs separately for each method and migration rate
fitted.GLM.prediction.list <- list() #initialize prediction list
fitted.GLM.summary.list <- list() #initialize summary list
fitted.GLM.list.index <- 1 #initialize list index
for (current.method in method.levels) {
  for (current.mig.tag in mig.level.table$mig.tag) {
    current.binomial.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == current.method & fitted.binomial.input.table$mig.tag == current.mig.tag, ] #subset one method and migration rate
    current.binomial.table <- current.binomial.table[is.finite(current.binomial.table$tdiv), ] #keep rows with finite divergence time
    if (nrow(current.binomial.table) == 0) next
    current.binomial.table <- current.binomial.table[order(current.binomial.table$tdiv), ] #order rows by divergence time
    prediction.tdiv <- seq(min(current.binomial.table$tdiv), max(current.binomial.table$tdiv), length.out = plot.fitted.prediction.n) #create prediction grid
    current.warnings <- character(0) #store GLM warnings
    if (all(current.binomial.table$n.k2 == 0)) {
      current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = 0, stringsAsFactors = FALSE) #flat zero prediction when K2 is never selected
      current.summary.row <- data.frame(method = current.method,
                                        mig.tag = as.character(current.mig.tag),
                                        mig = current.binomial.table$mig[1],
                                        n.vcf = nrow(current.binomial.table),
                                        n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                        n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                        min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                        n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                        n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                        all.zero = TRUE,
                                        all.one = FALSE,
                                        glm.fitted = FALSE,
                                        glm.converged = NA,
                                        glm.intercept = NA_real_,
                                        glm.tdiv.slope = NA_real_,
                                        tdiv.at.threshold.k2 = NA_real_,
                                        tdiv.at.threshold.k2.within.range = FALSE,
                                        glm.warning = NA_character_,
                                        glm.error = NA_character_,
                                        stringsAsFactors = FALSE) #store summary row
    } else if (all(current.binomial.table$n.not.k2 == 0)) {
      current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = 1, stringsAsFactors = FALSE) #flat one prediction when K2 is always selected
      current.summary.row <- data.frame(method = current.method,
                                        mig.tag = as.character(current.mig.tag),
                                        mig = current.binomial.table$mig[1],
                                        n.vcf = nrow(current.binomial.table),
                                        n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                        n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                        min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                        n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                        n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                        all.zero = FALSE,
                                        all.one = TRUE,
                                        glm.fitted = FALSE,
                                        glm.converged = NA,
                                        glm.intercept = NA_real_,
                                        glm.tdiv.slope = NA_real_,
                                        tdiv.at.threshold.k2 = NA_real_,
                                        tdiv.at.threshold.k2.within.range = FALSE,
                                        glm.warning = NA_character_,
                                        glm.error = NA_character_,
                                        stringsAsFactors = FALSE) #store summary row
    } else {
      current.fit <- withCallingHandlers(tryCatch(stats::glm(cbind(n.k2, n.not.k2) ~ tdiv, data = current.binomial.table, family = stats::binomial(link = "logit"), control = stats::glm.control(maxit = 100)), error = function(e) e),
                                         warning = function(w) {
                                           current.warnings <<- c(current.warnings, conditionMessage(w))
                                           invokeRestart("muffleWarning")
                                         }) #fit count-based binomial GLM and store warnings
      fit.ok <- inherits(current.fit, "glm") #check whether GLM fit succeeded
      if (fit.ok) {
        current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = stats::predict(current.fit, newdata = data.frame(tdiv = prediction.tdiv), type = "response"), stringsAsFactors = FALSE) #predict fitted K2 probability
        current.coefficients <- stats::coef(current.fit) #extract coefficients
        current.glm.intercept <- unname(current.coefficients["(Intercept)"]) #extract intercept
        current.glm.tdiv.slope <- unname(current.coefficients["tdiv"]) #extract tdiv slope
        if (is.finite(current.glm.intercept) && is.finite(current.glm.tdiv.slope) && current.glm.tdiv.slope != 0) {
          current.tdiv.at.threshold.k2 <- (stats::qlogis(plot.k2.threshold.proportion) - current.glm.intercept) / current.glm.tdiv.slope #solve logit(p) for threshold
        } else {
          current.tdiv.at.threshold.k2 <- NA_real_
        }
        current.tdiv.at.threshold.k2.within.range <- is.finite(current.tdiv.at.threshold.k2) && current.tdiv.at.threshold.k2 >= min(current.binomial.table$tdiv, na.rm = TRUE) && current.tdiv.at.threshold.k2 <= max(current.binomial.table$tdiv, na.rm = TRUE) #check whether threshold is inside simulated range
        current.summary.row <- data.frame(method = current.method,
                                          mig.tag = as.character(current.mig.tag),
                                          mig = current.binomial.table$mig[1],
                                          n.vcf = nrow(current.binomial.table),
                                          n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                          n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                          min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                          n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                          n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                          all.zero = all(current.binomial.table$proportion.k2 == 0),
                                          all.one = all(current.binomial.table$proportion.k2 == 1),
                                          glm.fitted = TRUE,
                                          glm.converged = isTRUE(current.fit$converged),
                                          glm.intercept = current.glm.intercept,
                                          glm.tdiv.slope = current.glm.tdiv.slope,
                                          tdiv.at.threshold.k2 = current.tdiv.at.threshold.k2,
                                          tdiv.at.threshold.k2.within.range = current.tdiv.at.threshold.k2.within.range,
                                          glm.warning = if (length(current.warnings) == 0) NA_character_ else paste(unique(current.warnings), collapse = " | "),
                                          glm.error = NA_character_,
                                          stringsAsFactors = FALSE) #store summary row
      } else {
        current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = NA_real_, stringsAsFactors = FALSE) #return missing predictions if GLM failed
        current.summary.row <- data.frame(method = current.method,
                                          mig.tag = as.character(current.mig.tag),
                                          mig = current.binomial.table$mig[1],
                                          n.vcf = nrow(current.binomial.table),
                                          n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                          n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                          min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                          n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                          n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                          all.zero = all(current.binomial.table$proportion.k2 == 0),
                                          all.one = all(current.binomial.table$proportion.k2 == 1),
                                          glm.fitted = FALSE,
                                          glm.converged = FALSE,
                                          glm.intercept = NA_real_,
                                          glm.tdiv.slope = NA_real_,
                                          tdiv.at.threshold.k2 = NA_real_,
                                          tdiv.at.threshold.k2.within.range = FALSE,
                                          glm.warning = if (length(current.warnings) == 0) NA_character_ else paste(unique(current.warnings), collapse = " | "),
                                          glm.error = conditionMessage(current.fit),
                                          stringsAsFactors = FALSE) #store summary row
      }
    }
    fitted.GLM.prediction.list[[fitted.GLM.list.index]] <- current.prediction.table #store prediction table
    fitted.GLM.summary.list[[fitted.GLM.list.index]] <- current.summary.row #store summary row
    fitted.GLM.list.index <- fitted.GLM.list.index + 1 #advance list index
  }
}


## Combine fitted GLM outputs
fitted.GLM.prediction.table <- do.call(rbind, fitted.GLM.prediction.list) #combine prediction tables
fitted.GLM.summary.table <- do.call(rbind, fitted.GLM.summary.list) #combine summary rows
fitted.GLM.prediction.table$method <- factor(as.character(fitted.GLM.prediction.table$method), levels = method.levels) #order methods
fitted.GLM.prediction.table$mig.tag <- factor(as.character(fitted.GLM.prediction.table$mig.tag), levels = mig.level.table$mig.tag) #order migration labels
fitted.GLM.summary.table$method <- factor(as.character(fitted.GLM.summary.table$method), levels = method.levels) #order methods
fitted.GLM.summary.table$mig.tag <- factor(as.character(fitted.GLM.summary.table$mig.tag), levels = mig.level.table$mig.tag) #order migration labels
fitted.GLM.summary.table <- fitted.GLM.summary.table[order(fitted.GLM.summary.table$method, fitted.GLM.summary.table$mig), ] #order summary table
rownames(fitted.GLM.summary.table) <- NULL
fitted.GLM.threshold.table <- fitted.GLM.summary.table[fitted.GLM.summary.table$tdiv.at.threshold.k2.within.range & is.finite(fitted.GLM.summary.table$tdiv.at.threshold.k2), ] #keep valid threshold rows


## Plot SOM results
SOM.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "SOM", ] #subset SOM observations
SOM.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "SOM" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset SOM fitted predictions
SOM.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "SOM", ] #subset SOM threshold rows
SOM.k2.fitted.plot <- ggplot2::ggplot(SOM.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = SOM.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = SOM.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(SOM.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "SOM_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = SOM.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot Weir and Cockerham Fst results
fst.plot.table <- plot.result.table[is.finite(plot.result.table$fst), ] #subset rows with finite Fst values
fst.plot.table <- fst.plot.table[order(fst.plot.table$mig, fst.plot.table$tdiv), ] #order rows for plotting
fst.plot.table$mig.tag <- factor(as.character(fst.plot.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels
SOM.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "SOM" & fitted.GLM.threshold.table$tdiv.at.threshold.k2.within.range & is.finite(fitted.GLM.threshold.table$tdiv.at.threshold.k2), ] #subset valid SOM threshold rows
SOM.k2.threshold.plot.table$mig.tag <- factor(as.character(SOM.k2.threshold.plot.table$mig.tag), levels = mig.level.table$mig.tag) #order SOM threshold migration-rate labels
fst.plot <- ggplot2::ggplot(fst.plot.table, ggplot2::aes(x = tdiv, y = fst, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = Fst_plot_point_size, alpha = Fst_plot_point_alpha) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, linewidth = Fst_plot_line_width, se = FALSE) +
  ggplot2::geom_vline(data = SOM.k2.threshold.plot.table,
                      ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag),
                      linewidth = plot.threshold.line.width,
                      linetype = plot.threshold.line.type,
                      show.legend = FALSE,
                      inherit.aes = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Weir and Cockerham Fst", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(fst.plot)
ggplot2::ggsave(file.path(results.directory, "Weir_Cockerham_Fst_by_tdiv_and_mig.svg"),
                plot = fst.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot de novo k-means/BIC results
deNovo.kmeans.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "deNovo.kmeans.BIC", ] #subset de novo k-means/BIC observations
deNovo.kmeans.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "deNovo.kmeans.BIC" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset de novo k-means/BIC fitted predictions
deNovo.kmeans.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "deNovo.kmeans.BIC", ] #subset de novo k-means/BIC threshold rows
deNovo.kmeans.k2.fitted.plot <- ggplot2::ggplot(deNovo.kmeans.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = deNovo.kmeans.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = deNovo.kmeans.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(deNovo.kmeans.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "deNovo_kmeans_BIC_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = deNovo.kmeans.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot sNMF K2 results
sNMF.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "sNMF", ] #subset sNMF observations
sNMF.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "sNMF" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset sNMF fitted predictions
sNMF.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "sNMF", ] #subset sNMF threshold rows
sNMF.k2.fitted.plot <- ggplot2::ggplot(sNMF.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = sNMF.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = sNMF.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(sNMF.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "sNMF_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = sNMF.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot STRUCTURE K2 results
STRUCTURE.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "STRUCTURE", ] #subset STRUCTURE observations
STRUCTURE.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "STRUCTURE" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset STRUCTURE fitted predictions
STRUCTURE.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "STRUCTURE", ] #subset STRUCTURE threshold rows
STRUCTURE.k2.fitted.plot <- ggplot2::ggplot(STRUCTURE.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = STRUCTURE.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = STRUCTURE.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(STRUCTURE.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "STRUCTURE_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = STRUCTURE.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Inspect fitted 50% threshold results
fitted.GLM.summary.table
fitted.GLM.threshold.table



## Create PCA plot for one VCF file

## Show available migration settings
available.migration.table <- unique(vcf.metadata.table[, c("mig.tag", "mig")]) #show all available migration settings including mig = 0
available.migration.table <- available.migration.table[order(available.migration.table$mig), ]
rownames(available.migration.table) <- NULL
print(available.migration.table)
if (!(migration.vcf.mig.for.PCA %in% unique(vcf.metadata.table$mig))) stop(paste("migration.vcf.mig.for.PCA must be one of:", paste(sort(unique(vcf.metadata.table$mig)), collapse = ", ")))


## Select file number N among VCF files for selected migration rate
migration.vcf.metadata.table <- vcf.metadata.table[vcf.metadata.table$mig == migration.vcf.mig.for.PCA, ] #keep selected VCF files only
migration.vcf.metadata.table <- migration.vcf.metadata.table[order(migration.vcf.metadata.table$tdiv), ] #order selected files by divergence time
rownames(migration.vcf.metadata.table) <- NULL
if (nrow(migration.vcf.metadata.table) == 0) stop(paste("No VCF files found for migration rate", migration.vcf.mig.for.PCA))
if (migration.vcf.file.number.for.PCA < 1 || migration.vcf.file.number.for.PCA > nrow(migration.vcf.metadata.table)) stop(paste("migration.vcf.file.number.for.PCA must be between 1 and", nrow(migration.vcf.metadata.table), "for migration rate", migration.vcf.mig.for.PCA))
current.pca.metadata <- migration.vcf.metadata.table[migration.vcf.file.number.for.PCA, ] #select VCF file N for selected migration rate
cat("\nCreating PCA for VCF file", migration.vcf.file.number.for.PCA, "of", nrow(migration.vcf.metadata.table), "for migration rate", migration.vcf.mig.for.PCA, ":", current.pca.metadata$file, "\n")
cat("Migration-index tag:", current.pca.metadata$mig.index, "\n")
cat("Divergence time:", current.pca.metadata$tdiv, "generations\n")
cat("Migration tag:", current.pca.metadata$mig.tag, "\n")
cat("Migration rate:", current.pca.metadata$mig, "\n")


## Read VCF and process SNP matrix
pca.genind.object.raw <- read.fastsimcoal2.genind.for.analysis(vcf.file.path = current.pca.metadata$full.path,
                                                               n.per.population = analysis.sample.n.per.population,
                                                               sampling.group = current.pca.metadata$simulation.group,
                                                               sampling.seed = analysis.sample.random.seed) #read selected VCF, assign populations, and optionally subset individuals per population


## Calculate Weir and Cockerham Fst before PCA processing
pca.fst.value <- calculate.overall.Fst(genind.object = pca.genind.object.raw) #calculate overall Weir and Cockerham Fst
cat("Weir and Cockerham Fst between pop1 and pop2:", round(pca.fst.value, 4), "\n")


## Process SNP matrix for PCA
pca.snp.matrix <- process.SNP.data.SOM(genind.input = pca.genind.object.raw, verbose = FALSE) #process SNP data for PCA
pca.snp.matrix <- as.matrix(pca.snp.matrix) #coerce to matrix
storage.mode(pca.snp.matrix) <- "numeric" #ensure numeric matrix
if (is.null(rownames(pca.snp.matrix))) rownames(pca.snp.matrix) <- adegenet::indNames(pca.genind.object.raw)


## Remove invalid loci
valid.pca.loci <- apply(pca.snp.matrix, 2, function(locus.values) {
  non.missing.values <- locus.values[!is.na(locus.values)]
  length(non.missing.values) > 1 && stats::var(non.missing.values) > 0
}) #identify loci with nonzero variance
pca.snp.matrix <- pca.snp.matrix[, valid.pca.loci, drop = FALSE] #keep valid loci only
if (ncol(pca.snp.matrix) < 2) stop("Fewer than two valid loci available for PCA")


## Impute remaining missing values for PCA
if (any(is.na(pca.snp.matrix))) {
  pca.snp.matrix <- apply(pca.snp.matrix, 2, function(locus.values) {
    locus.values[is.na(locus.values)] <- mean(locus.values, na.rm = TRUE)
    locus.values
  }) #replace remaining missing values with locus means
  pca.snp.matrix <- as.matrix(pca.snp.matrix)
}


## Run PCA and create plot table
pca.object <- stats::prcomp(pca.snp.matrix, center = TRUE, scale. = FALSE) #run PCA
pca.variance <- pca.object$sdev^2 / sum(pca.object$sdev^2) #calculate variance explained
pca.plot.table <- data.frame(individual = rownames(pca.object$x),
                             PC1 = pca.object$x[, 1],
                             PC2 = pca.object$x[, 2],
                             population = as.character(adegenet::pop(pca.genind.object.raw))[match(rownames(pca.object$x), adegenet::indNames(pca.genind.object.raw))],
                             file = current.pca.metadata$file,
                             mig.index = current.pca.metadata$mig.index,
                             tdiv = current.pca.metadata$tdiv,
                             mig.tag = current.pca.metadata$mig.tag,
                             mig = current.pca.metadata$mig,
                             fst.pop1.pop2 = pca.fst.value,
                             stringsAsFactors = FALSE) #create PCA plotting table


## Plot PCA with population groupings
pca.plot <- ggplot2::ggplot(pca.plot.table, ggplot2::aes(x = PC1, y = PC2, color = population, group = population)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(title = paste0("PCA - ", current.pca.metadata$file, " | tdiv = ", round(current.pca.metadata$tdiv), " generations | mig = ", current.pca.metadata$mig, " | Fst = ", round(pca.fst.value, 4)),
                x = paste0("PC1 (", round(pca.variance[1] * 100, 2), "%)"),
                y = paste0("PC2 (", round(pca.variance[2] * 100, 2), "%)"),
                color = "Population") +
  ggplot2::theme_classic()
pca.plot


## Print selected PCA summary
cat("Divergence:", round(current.pca.metadata$tdiv / 1000000, 2), "million generations\n")
cat("Migration rate:", current.pca.metadata$mig, "\n")
cat("Weir and Cockerham Fst:", round(pca.fst.value, 4), "\n")




#### Run workflow for VCF files - asymmetric_24 #################################
SOM_grid_multiplier <- 4
analysis.sample.n.per.population <- analysis.sample.n.by.population.asymmetric_24 #set individuals per population for asymmetric_24
analysis.set.name <- "asymmetric_24" #define analysis set name
results.directory <- file.path(results.root.directory, analysis.set.name) #define output directory for this analysis set
dir.create(results.directory, recursive = TRUE, showWarnings = FALSE) #create output directory if needed
intermediate.results.directory <- file.path(results.directory, "intermediate_results") #define intermediate results directory for this analysis set
dir.create(intermediate.results.directory, recursive = TRUE, showWarnings = FALSE) #create intermediate results directory if needed
cat("\n\nRunning analysis set:", analysis.set.name, "\n")
cat("Individuals per population:", analysis.sample.n.per.population, "\n")
cat("Results directory:", results.directory, "\n")
cat("Random sampling base seed:", analysis.sample.random.seed, "\n")


## List VCF files
vcf.file.pattern <- "^sim[0-9]+_tdiv[0-9.]+_mig.+\\.vcf$" #define VCF filename pattern
vcf.file.paths <- list.files(vcf.directory, pattern = vcf.file.pattern, full.names = TRUE) #list VCF files only
if (length(vcf.file.paths) == 0) stop("No VCF files found")
vcf.file.names <- basename(vcf.file.paths) #extract file names


## Parse filename metadata
vcf.metadata.list <- lapply(vcf.file.names, parse.vcf.filename.metadata) #parse metadata
vcf.metadata.table <- do.call(rbind, vcf.metadata.list) #combine metadata rows
vcf.metadata.table$full.path <- vcf.file.paths #add full path


## Order metadata table
vcf.metadata.table <- vcf.metadata.table[order(vcf.metadata.table$mig.index,
                                               vcf.metadata.table$tdiv,
                                               vcf.metadata.table$mig), ] #order rows
rownames(vcf.metadata.table) <- NULL
vcf.metadata.table$simulation.group <- match(vcf.metadata.table$tdiv, sort(unique(vcf.metadata.table$tdiv))) #define matched sampling group across migration rates for each divergence time
head(vcf.metadata.table) #show metadata for files


## Create empty lists to store results
result.summary.list <- vector("list", nrow(vcf.metadata.table)) #initialize result summary list
optim.k.summary.list <- vector("list", nrow(vcf.metadata.table)) #initialize optim_k_summary list


## Loop over VCF files
for (vcf.index in seq_len(nrow(vcf.metadata.table))) {
  current.metadata <- vcf.metadata.table[vcf.index, ] #extract current metadata
  cat("\nProcessing", vcf.index, "of", nrow(vcf.metadata.table), ":", current.metadata$file, "\n")
  if (!override && check.vcf.intermediate.results.available(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory)) {
    cat("Intermediate results already found for", current.metadata$file, "- loading saved results\n")
    current.output <- load.vcf.intermediate.results(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory) #load saved per-VCF results
  } else {
    current.output <- tryCatch({
      run.SOM.workflow.for.one.vcf(vcf.file.path = current.metadata$full.path,
                                   vcf.file.metadata = current.metadata,
                                   sample.n.per.population = analysis.sample.n.per.population,
                                   sampling.group = current.metadata$simulation.group,
                                   sampling.seed = analysis.sample.random.seed)
    }, error = function(e) {
      result.summary.row <- data.frame(file = current.metadata$file,
                                       mig.index = current.metadata$mig.index,
                                       tdiv = current.metadata$tdiv,
                                       mig.tag = current.metadata$mig.tag,
                                       mig = current.metadata$mig,
                                       n.ind.raw = NA_integer_,
                                       n.loc.raw = NA_integer_,
                                       n.samples.processed = NA_integer_,
                                       n.loci.processed = NA_integer_,
                                       proportion.k2 = NA_real_,
                                       status = "error",
                                       error = conditionMessage(e),
                                       stringsAsFactors = FALSE)
      optim.k.summary.table <- data.frame(file = current.metadata$file,
                                          mig.index = current.metadata$mig.index,
                                          tdiv = current.metadata$tdiv,
                                          mig.tag = current.metadata$mig.tag,
                                          mig = current.metadata$mig,
                                          n.ind.raw = NA_integer_,
                                          n.loc.raw = NA_integer_,
                                          n.samples.processed = NA_integer_,
                                          n.loci.processed = NA_integer_,
                                          status = "error",
                                          error = conditionMessage(e),
                                          Count = NA_real_,
                                          Proportion = NA_real_,
                                          k.label = NA_character_,
                                          optim.k.row = NA_integer_,
                                          stringsAsFactors = FALSE)
      list(result.summary.row = result.summary.row, optim.k.summary.table = optim.k.summary.table)
    })
    save.vcf.intermediate.results(result.summary.row = current.output$result.summary.row,
                                  optim.k.summary.table = current.output$optim.k.summary.table,
                                  vcf.file.metadata = current.metadata,
                                  intermediate.results.directory = intermediate.results.directory) #save per-VCF results immediately
  }
  result.summary.list[[vcf.index]] <- current.output$result.summary.row #store result summary row
  optim.k.summary.list[[vcf.index]] <- current.output$optim.k.summary.table #store optim_k_summary table
  stop.if.status.not.ok(status.value = current.output$result.summary.row$status,
                        error.value = current.output$result.summary.row$error,
                        file.value = current.metadata$file,
                        method.label = "SOM") #stop immediately if SOM failed
}


## Combine SOM results
result.table <- do.call(rbind, result.summary.list) #combine result summary rows
rownames(result.table) <- NULL
optim.k.result.table <- do.call(rbind, optim.k.summary.list) #combine optim_k_summary rows
rownames(optim.k.result.table) <- NULL


## Calculate Fst and append to SOM result tables
Fst.output <- calculate.and.append.Fst(vcf.metadata.table = vcf.metadata.table,
                                       result.table = result.table,
                                       optim.k.result.table = optim.k.result.table,
                                       intermediate.results.directory = intermediate.results.directory,
                                       override = override,
                                       sample.n.per.population = analysis.sample.n.per.population,
                                       sampling.seed = analysis.sample.random.seed) #calculate and append Fst
result.table <- Fst.output$result.table #update result table with Fst
optim.k.result.table <- Fst.output$optim.k.result.table #update optim_k result table with Fst
fst.table <- Fst.output$fst.table #extract Fst table


## Calculate de novo k-means/BIC and append to result tables
deNovo.kmeans.output <- calculate.and.append.deNovo.kmeans(vcf.metadata.table = vcf.metadata.table,
                                                           result.table = result.table,
                                                           optim.k.result.table = optim.k.result.table,
                                                           intermediate.results.directory = intermediate.results.directory,
                                                           override = override,
                                                           BIC.thresh = deNovo.kmeans.BIC.thresh,
                                                           n.iter = deNovo.kmeans.n.iter,
                                                           n.start = deNovo.kmeans.n.start,
                                                           max.n.clust = deNovo.kmeans.max.n.clust,
                                                           max.proportion.PCs = deNovo.kmeans.max.proportion.PCs,
                                                           center = deNovo.kmeans.center,
                                                           scale = deNovo.kmeans.scale,
                                                           sample.n.per.population = analysis.sample.n.per.population,
                                                           sampling.seed = analysis.sample.random.seed) #calculate de novo k-means/BIC
result.table <- deNovo.kmeans.output$result.table #update result table with de novo k-means/BIC
optim.k.result.table <- deNovo.kmeans.output$optim.k.result.table #update optim_k result table with de novo k-means/BIC
deNovo.kmeans.table <- deNovo.kmeans.output$deNovo.kmeans.table #extract de novo k-means/BIC table


## Calculate sNMF and append to result tables
sNMF.output <- calculate.and.append.sNMF(vcf.metadata.table = vcf.metadata.table,
                                         result.table = result.table,
                                         optim.k.result.table = optim.k.result.table,
                                         intermediate.results.directory = intermediate.results.directory,
                                         override = override,
                                         K.values = sNMF.K.values,
                                         repetitions = sNMF.repetitions,
                                         CPU = 1,
                                         ploidy = sNMF.ploidy,
                                         cross.entropy.thresh = sNMF.cross.entropy.thresh,
                                         seed = sNMF.seed,
                                         sample.n.per.population = analysis.sample.n.per.population,
                                         sampling.seed = analysis.sample.random.seed) #calculate sNMF
result.table <- sNMF.output$result.table #update result table with sNMF
optim.k.result.table <- sNMF.output$optim.k.result.table #update optim_k result table with sNMF
sNMF.table <- sNMF.output$sNMF.table #extract sNMF table


## Save result tables
write.csv(result.table, file.path(results.directory, "fastsimcoal2_SOM_results.csv"), row.names = FALSE) #save main SOM result summary
write.csv(optim.k.result.table, file.path(results.directory, "fastsimcoal2_SOM_optim_k_summary.csv"), row.names = FALSE) #save full optim_k_summary table
write.csv(fst.table, file.path(results.directory, "fastsimcoal2_Weir_Cockerham_Fst.csv"), row.names = FALSE) #save Fst table
write.csv(deNovo.kmeans.table, file.path(results.directory, "fastsimcoal2_de_novo_kmeans_BIC.csv"), row.names = FALSE) #save de novo k-means/BIC table
write.csv(sNMF.table, file.path(results.directory, "fastsimcoal2_sNMF_results.csv"), row.names = FALSE) #save sNMF table


## Inspect results
head(result.table)
optim.k.result.table
head(fst.table)
head(deNovo.kmeans.table)
head(sNMF.table)
table(result.table$status)
table(result.table$fst.status)
table(result.table$deNovo.kmeans.status)
table(result.table$deNovo.kmeans.best.k, useNA = "ifany")
table(result.table$sNMF.status)


## Plot fitted K2 proportions for SOM, de novo k-means/BIC, and sNMF
plot.result.table <- result.table[result.table$status == "ok", ] #subset successful rows
plot.result.table <- plot.result.table[order(plot.result.table$mig, plot.result.table$tdiv), ] #order rows for plotting
if (nrow(plot.result.table) == 0) stop("No successful results available for fitted plotting")
mig.level.table <- unique(plot.result.table[order(plot.result.table$mig), c("mig.tag", "mig")]) #define unique ordered migration-rate levels
mig.level.table$mig.tag <- as.character(mig.level.table$mig.tag) #force migration-rate labels to discrete character values
plot.result.table$mig.tag <- factor(as.character(plot.result.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels by numeric migration rate
mig.colors <- setNames(viridisLite::viridis(nrow(mig.level.table), option = "D"), mig.level.table$mig.tag) #define discrete ordered colors


## Create SOM binomial count table
SOM.optim.k.count.table <- optim.k.result.table[!is.na(optim.k.result.table$Count) & optim.k.result.table$file %in% plot.result.table$file, ] #keep rows with SOM replicate counts
SOM.optim.k.count.table$Count <- as.numeric(as.character(SOM.optim.k.count.table$Count)) #force counts to numeric
SOM.total.count.table <- stats::aggregate(Count ~ file + mig.tag + mig + tdiv, data = SOM.optim.k.count.table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE)) #total retained SOM replicates per VCF
colnames(SOM.total.count.table)[colnames(SOM.total.count.table) == "Count"] <- "n.total"
SOM.k2.count.table <- SOM.optim.k.count.table[SOM.optim.k.count.table$k.label == "k2", ] #keep k2 count rows
if (nrow(SOM.k2.count.table) > 0) {
  SOM.k2.count.table <- stats::aggregate(Count ~ file + mig.tag + mig + tdiv, data = SOM.k2.count.table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE)) #number of retained SOM replicates choosing k = 2
  colnames(SOM.k2.count.table)[colnames(SOM.k2.count.table) == "Count"] <- "n.k2"
} else {
  SOM.k2.count.table <- SOM.total.count.table[, c("file", "mig.tag", "mig", "tdiv")] #create empty k2 count table if no k2 rows exist
  SOM.k2.count.table$n.k2 <- 0
}
SOM.binomial.table <- merge(SOM.total.count.table, SOM.k2.count.table, by = c("file", "mig.tag", "mig", "tdiv"), all.x = TRUE) #merge total and k2 counts
SOM.binomial.table$n.k2[is.na(SOM.binomial.table$n.k2)] <- 0 #set missing k2 rows to zero
SOM.binomial.table$n.not.k2 <- SOM.binomial.table$n.total - SOM.binomial.table$n.k2 #number not choosing k = 2
SOM.binomial.table$proportion.k2 <- SOM.binomial.table$n.k2 / SOM.binomial.table$n.total #calculate K2 proportion
SOM.binomial.table$method <- "SOM" #add method label
SOM.binomial.table <- SOM.binomial.table[, c("method", "file", "mig.tag", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")] #order columns


## Create de novo k-means/BIC binomial table
deNovo.kmeans.binomial.source.table <- result.table[result.table$status == "ok" & !is.na(result.table$deNovo.kmeans.best.k), ] #subset valid de novo k-means/BIC rows
if ("deNovo.kmeans.status" %in% colnames(deNovo.kmeans.binomial.source.table)) deNovo.kmeans.binomial.source.table <- deNovo.kmeans.binomial.source.table[deNovo.kmeans.binomial.source.table$deNovo.kmeans.status == "ok", ] #keep successful de novo k-means/BIC rows
deNovo.kmeans.best.k <- as.integer(as.character(deNovo.kmeans.binomial.source.table$deNovo.kmeans.best.k)) #extract de novo k-means/BIC best K
deNovo.kmeans.binomial.table <- data.frame(method = "deNovo.kmeans.BIC",
                                           file = deNovo.kmeans.binomial.source.table$file,
                                           mig.tag = deNovo.kmeans.binomial.source.table$mig.tag,
                                           mig = deNovo.kmeans.binomial.source.table$mig,
                                           tdiv = deNovo.kmeans.binomial.source.table$tdiv,
                                           n.k2 = as.integer(deNovo.kmeans.best.k == 2L),
                                           n.not.k2 = as.integer(deNovo.kmeans.best.k != 2L),
                                           proportion.k2 = as.numeric(deNovo.kmeans.best.k == 2L),
                                           stringsAsFactors = FALSE) #create binary de novo k-means/BIC table


## Create sNMF binomial table
sNMF.binomial.source.table <- result.table[result.table$status == "ok" & !is.na(result.table$sNMF.best.k), ] #subset valid sNMF rows
if ("sNMF.status" %in% colnames(sNMF.binomial.source.table)) sNMF.binomial.source.table <- sNMF.binomial.source.table[sNMF.binomial.source.table$sNMF.status == "ok", ] #keep successful sNMF rows
sNMF.best.k <- as.integer(as.character(sNMF.binomial.source.table$sNMF.best.k)) #extract sNMF best K
sNMF.binomial.table <- data.frame(method = "sNMF",
                                  file = sNMF.binomial.source.table$file,
                                  mig.tag = sNMF.binomial.source.table$mig.tag,
                                  mig = sNMF.binomial.source.table$mig,
                                  tdiv = sNMF.binomial.source.table$tdiv,
                                  n.k2 = as.integer(sNMF.best.k == 2L),
                                  n.not.k2 = as.integer(sNMF.best.k != 2L),
                                  proportion.k2 = as.numeric(sNMF.best.k == 2L),
                                  stringsAsFactors = FALSE) #create binary sNMF table


## Create STRUCTURE binomial table
STRUCTURE.results <- create.STRUCTURE.mean.lnprob.binomial.table(analysis.set.name = analysis.set.name,
                                                                 structure.loglik.directory = STRUCTURE.loglik.directory,
                                                                 results.directory = results.directory,
                                                                 expected.structure.replicates = expected.structure.replicates,
                                                                 STRUCTURE.mean.lnprob.delta.threshold = STRUCTURE.mean.lnprob.delta.threshold) #create STRUCTURE mean-lnP best-K table
STRUCTURE.binomial.table <- STRUCTURE.results$STRUCTURE.binomial.table #extract STRUCTURE binomial table
STRUCTURE.binomial.table$mig.tag <- vapply(STRUCTURE.binomial.table$mig, function(current.mig) {
  matching.mig.index <- which(abs(mig.level.table$mig - current.mig) < sqrt(.Machine$double.eps))
  if (length(matching.mig.index) != 1) return(NA_character_)
  return(as.character(mig.level.table$mig.tag[matching.mig.index]))
}, character(1)) #standardize STRUCTURE migration labels to match VCF migration labels
if (any(is.na(STRUCTURE.binomial.table$mig.tag))) stop("Some STRUCTURE migration rates could not be matched to VCF migration labels")
STRUCTURE.mean.lnprob.by.vcf.table <- STRUCTURE.results$STRUCTURE.mean.lnprob.by.vcf.table #extract STRUCTURE VCF-level table
STRUCTURE.rep.comparison.table <- STRUCTURE.results$STRUCTURE.rep.comparison.table #extract STRUCTURE replicate comparison table


## Combine method-specific binomial tables
fitted.binomial.input.table <- rbind(SOM.binomial.table, deNovo.kmeans.binomial.table, sNMF.binomial.table, STRUCTURE.binomial.table) #combine all methods
method.levels <- c("SOM", "deNovo.kmeans.BIC", "sNMF", "STRUCTURE") #define method order
fitted.binomial.input.table$method <- factor(as.character(fitted.binomial.input.table$method), levels = method.levels) #order methods
fitted.binomial.input.table$mig.tag <- factor(as.character(fitted.binomial.input.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels


## Fit binomial GLMs separately for each method and migration rate
fitted.GLM.prediction.list <- list() #initialize prediction list
fitted.GLM.summary.list <- list() #initialize summary list
fitted.GLM.list.index <- 1 #initialize list index
for (current.method in method.levels) {
  for (current.mig.tag in mig.level.table$mig.tag) {
    current.binomial.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == current.method & fitted.binomial.input.table$mig.tag == current.mig.tag, ] #subset one method and migration rate
    current.binomial.table <- current.binomial.table[is.finite(current.binomial.table$tdiv), ] #keep rows with finite divergence time
    if (nrow(current.binomial.table) == 0) next
    current.binomial.table <- current.binomial.table[order(current.binomial.table$tdiv), ] #order rows by divergence time
    prediction.tdiv <- seq(min(current.binomial.table$tdiv), max(current.binomial.table$tdiv), length.out = plot.fitted.prediction.n) #create prediction grid
    current.warnings <- character(0) #store GLM warnings
    if (all(current.binomial.table$n.k2 == 0)) {
      current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = 0, stringsAsFactors = FALSE) #flat zero prediction when K2 is never selected
      current.summary.row <- data.frame(method = current.method,
                                        mig.tag = as.character(current.mig.tag),
                                        mig = current.binomial.table$mig[1],
                                        n.vcf = nrow(current.binomial.table),
                                        n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                        n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                        min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                        n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                        n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                        all.zero = TRUE,
                                        all.one = FALSE,
                                        glm.fitted = FALSE,
                                        glm.converged = NA,
                                        glm.intercept = NA_real_,
                                        glm.tdiv.slope = NA_real_,
                                        tdiv.at.threshold.k2 = NA_real_,
                                        tdiv.at.threshold.k2.within.range = FALSE,
                                        glm.warning = NA_character_,
                                        glm.error = NA_character_,
                                        stringsAsFactors = FALSE) #store summary row
    } else if (all(current.binomial.table$n.not.k2 == 0)) {
      current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = 1, stringsAsFactors = FALSE) #flat one prediction when K2 is always selected
      current.summary.row <- data.frame(method = current.method,
                                        mig.tag = as.character(current.mig.tag),
                                        mig = current.binomial.table$mig[1],
                                        n.vcf = nrow(current.binomial.table),
                                        n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                        n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                        min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                        n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                        n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                        all.zero = FALSE,
                                        all.one = TRUE,
                                        glm.fitted = FALSE,
                                        glm.converged = NA,
                                        glm.intercept = NA_real_,
                                        glm.tdiv.slope = NA_real_,
                                        tdiv.at.threshold.k2 = NA_real_,
                                        tdiv.at.threshold.k2.within.range = FALSE,
                                        glm.warning = NA_character_,
                                        glm.error = NA_character_,
                                        stringsAsFactors = FALSE) #store summary row
    } else {
      current.fit <- withCallingHandlers(tryCatch(stats::glm(cbind(n.k2, n.not.k2) ~ tdiv, data = current.binomial.table, family = stats::binomial(link = "logit"), control = stats::glm.control(maxit = 100)), error = function(e) e),
                                         warning = function(w) {
                                           current.warnings <<- c(current.warnings, conditionMessage(w))
                                           invokeRestart("muffleWarning")
                                         }) #fit count-based binomial GLM and store warnings
      fit.ok <- inherits(current.fit, "glm") #check whether GLM fit succeeded
      if (fit.ok) {
        current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = stats::predict(current.fit, newdata = data.frame(tdiv = prediction.tdiv), type = "response"), stringsAsFactors = FALSE) #predict fitted K2 probability
        current.coefficients <- stats::coef(current.fit) #extract coefficients
        current.glm.intercept <- unname(current.coefficients["(Intercept)"]) #extract intercept
        current.glm.tdiv.slope <- unname(current.coefficients["tdiv"]) #extract tdiv slope
        if (is.finite(current.glm.intercept) && is.finite(current.glm.tdiv.slope) && current.glm.tdiv.slope != 0) {
          current.tdiv.at.threshold.k2 <- (stats::qlogis(plot.k2.threshold.proportion) - current.glm.intercept) / current.glm.tdiv.slope #solve logit(p) for threshold
        } else {
          current.tdiv.at.threshold.k2 <- NA_real_
        }
        current.tdiv.at.threshold.k2.within.range <- is.finite(current.tdiv.at.threshold.k2) && current.tdiv.at.threshold.k2 >= min(current.binomial.table$tdiv, na.rm = TRUE) && current.tdiv.at.threshold.k2 <= max(current.binomial.table$tdiv, na.rm = TRUE) #check whether threshold is inside simulated range
        current.summary.row <- data.frame(method = current.method,
                                          mig.tag = as.character(current.mig.tag),
                                          mig = current.binomial.table$mig[1],
                                          n.vcf = nrow(current.binomial.table),
                                          n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                          n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                          min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                          n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                          n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                          all.zero = all(current.binomial.table$proportion.k2 == 0),
                                          all.one = all(current.binomial.table$proportion.k2 == 1),
                                          glm.fitted = TRUE,
                                          glm.converged = isTRUE(current.fit$converged),
                                          glm.intercept = current.glm.intercept,
                                          glm.tdiv.slope = current.glm.tdiv.slope,
                                          tdiv.at.threshold.k2 = current.tdiv.at.threshold.k2,
                                          tdiv.at.threshold.k2.within.range = current.tdiv.at.threshold.k2.within.range,
                                          glm.warning = if (length(current.warnings) == 0) NA_character_ else paste(unique(current.warnings), collapse = " | "),
                                          glm.error = NA_character_,
                                          stringsAsFactors = FALSE) #store summary row
      } else {
        current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = NA_real_, stringsAsFactors = FALSE) #return missing predictions if GLM failed
        current.summary.row <- data.frame(method = current.method,
                                          mig.tag = as.character(current.mig.tag),
                                          mig = current.binomial.table$mig[1],
                                          n.vcf = nrow(current.binomial.table),
                                          n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                          n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                          min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                          n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                          n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                          all.zero = all(current.binomial.table$proportion.k2 == 0),
                                          all.one = all(current.binomial.table$proportion.k2 == 1),
                                          glm.fitted = FALSE,
                                          glm.converged = FALSE,
                                          glm.intercept = NA_real_,
                                          glm.tdiv.slope = NA_real_,
                                          tdiv.at.threshold.k2 = NA_real_,
                                          tdiv.at.threshold.k2.within.range = FALSE,
                                          glm.warning = if (length(current.warnings) == 0) NA_character_ else paste(unique(current.warnings), collapse = " | "),
                                          glm.error = conditionMessage(current.fit),
                                          stringsAsFactors = FALSE) #store summary row
      }
    }
    fitted.GLM.prediction.list[[fitted.GLM.list.index]] <- current.prediction.table #store prediction table
    fitted.GLM.summary.list[[fitted.GLM.list.index]] <- current.summary.row #store summary row
    fitted.GLM.list.index <- fitted.GLM.list.index + 1 #advance list index
  }
}


## Combine fitted GLM outputs
fitted.GLM.prediction.table <- do.call(rbind, fitted.GLM.prediction.list) #combine prediction tables
fitted.GLM.summary.table <- do.call(rbind, fitted.GLM.summary.list) #combine summary rows
fitted.GLM.prediction.table$method <- factor(as.character(fitted.GLM.prediction.table$method), levels = method.levels) #order methods
fitted.GLM.prediction.table$mig.tag <- factor(as.character(fitted.GLM.prediction.table$mig.tag), levels = mig.level.table$mig.tag) #order migration labels
fitted.GLM.summary.table$method <- factor(as.character(fitted.GLM.summary.table$method), levels = method.levels) #order methods
fitted.GLM.summary.table$mig.tag <- factor(as.character(fitted.GLM.summary.table$mig.tag), levels = mig.level.table$mig.tag) #order migration labels
fitted.GLM.summary.table <- fitted.GLM.summary.table[order(fitted.GLM.summary.table$method, fitted.GLM.summary.table$mig), ] #order summary table
rownames(fitted.GLM.summary.table) <- NULL
fitted.GLM.threshold.table <- fitted.GLM.summary.table[fitted.GLM.summary.table$tdiv.at.threshold.k2.within.range & is.finite(fitted.GLM.summary.table$tdiv.at.threshold.k2), ] #keep valid threshold rows


## Plot SOM results
SOM.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "SOM", ] #subset SOM observations
SOM.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "SOM" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset SOM fitted predictions
SOM.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "SOM", ] #subset SOM threshold rows
SOM.k2.fitted.plot <- ggplot2::ggplot(SOM.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = SOM.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = SOM.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(SOM.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "SOM_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = SOM.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot Weir and Cockerham Fst results
fst.plot.table <- plot.result.table[is.finite(plot.result.table$fst), ] #subset rows with finite Fst values
fst.plot.table <- fst.plot.table[order(fst.plot.table$mig, fst.plot.table$tdiv), ] #order rows for plotting
fst.plot.table$mig.tag <- factor(as.character(fst.plot.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels
SOM.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "SOM" & fitted.GLM.threshold.table$tdiv.at.threshold.k2.within.range & is.finite(fitted.GLM.threshold.table$tdiv.at.threshold.k2), ] #subset valid SOM threshold rows
SOM.k2.threshold.plot.table$mig.tag <- factor(as.character(SOM.k2.threshold.plot.table$mig.tag), levels = mig.level.table$mig.tag) #order SOM threshold migration-rate labels
fst.plot <- ggplot2::ggplot(fst.plot.table, ggplot2::aes(x = tdiv, y = fst, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = Fst_plot_point_size, alpha = Fst_plot_point_alpha) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, linewidth = Fst_plot_line_width, se = FALSE) +
  ggplot2::geom_vline(data = SOM.k2.threshold.plot.table,
                      ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag),
                      linewidth = plot.threshold.line.width,
                      linetype = plot.threshold.line.type,
                      show.legend = FALSE,
                      inherit.aes = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Weir and Cockerham Fst", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(fst.plot)
ggplot2::ggsave(file.path(results.directory, "Weir_Cockerham_Fst_by_tdiv_and_mig.svg"),
                plot = fst.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot de novo k-means/BIC results
deNovo.kmeans.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "deNovo.kmeans.BIC", ] #subset de novo k-means/BIC observations
deNovo.kmeans.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "deNovo.kmeans.BIC" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset de novo k-means/BIC fitted predictions
deNovo.kmeans.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "deNovo.kmeans.BIC", ] #subset de novo k-means/BIC threshold rows
deNovo.kmeans.k2.fitted.plot <- ggplot2::ggplot(deNovo.kmeans.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = deNovo.kmeans.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = deNovo.kmeans.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(deNovo.kmeans.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "deNovo_kmeans_BIC_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = deNovo.kmeans.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot sNMF K2 results
sNMF.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "sNMF", ] #subset sNMF observations
sNMF.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "sNMF" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset sNMF fitted predictions
sNMF.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "sNMF", ] #subset sNMF threshold rows
sNMF.k2.fitted.plot <- ggplot2::ggplot(sNMF.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = sNMF.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = sNMF.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(sNMF.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "sNMF_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = sNMF.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot STRUCTURE K2 results
STRUCTURE.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "STRUCTURE", ] #subset STRUCTURE observations
STRUCTURE.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "STRUCTURE" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset STRUCTURE fitted predictions
STRUCTURE.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "STRUCTURE", ] #subset STRUCTURE threshold rows
STRUCTURE.k2.fitted.plot <- ggplot2::ggplot(STRUCTURE.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = STRUCTURE.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = STRUCTURE.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(STRUCTURE.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "STRUCTURE_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = STRUCTURE.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Inspect fitted 50% threshold results
fitted.GLM.summary.table
fitted.GLM.threshold.table



## Create PCA plot for one VCF file

## Show available migration settings
available.migration.table <- unique(vcf.metadata.table[, c("mig.tag", "mig")]) #show all available migration settings including mig = 0
available.migration.table <- available.migration.table[order(available.migration.table$mig), ]
rownames(available.migration.table) <- NULL
print(available.migration.table)
if (!(migration.vcf.mig.for.PCA %in% unique(vcf.metadata.table$mig))) stop(paste("migration.vcf.mig.for.PCA must be one of:", paste(sort(unique(vcf.metadata.table$mig)), collapse = ", ")))


## Select file number N among VCF files for selected migration rate
migration.vcf.metadata.table <- vcf.metadata.table[vcf.metadata.table$mig == migration.vcf.mig.for.PCA, ] #keep selected VCF files only
migration.vcf.metadata.table <- migration.vcf.metadata.table[order(migration.vcf.metadata.table$tdiv), ] #order selected files by divergence time
rownames(migration.vcf.metadata.table) <- NULL
if (nrow(migration.vcf.metadata.table) == 0) stop(paste("No VCF files found for migration rate", migration.vcf.mig.for.PCA))
if (migration.vcf.file.number.for.PCA < 1 || migration.vcf.file.number.for.PCA > nrow(migration.vcf.metadata.table)) stop(paste("migration.vcf.file.number.for.PCA must be between 1 and", nrow(migration.vcf.metadata.table), "for migration rate", migration.vcf.mig.for.PCA))
current.pca.metadata <- migration.vcf.metadata.table[migration.vcf.file.number.for.PCA, ] #select VCF file N for selected migration rate
cat("\nCreating PCA for VCF file", migration.vcf.file.number.for.PCA, "of", nrow(migration.vcf.metadata.table), "for migration rate", migration.vcf.mig.for.PCA, ":", current.pca.metadata$file, "\n")
cat("Migration-index tag:", current.pca.metadata$mig.index, "\n")
cat("Divergence time:", current.pca.metadata$tdiv, "generations\n")
cat("Migration tag:", current.pca.metadata$mig.tag, "\n")
cat("Migration rate:", current.pca.metadata$mig, "\n")


## Read VCF and process SNP matrix
pca.genind.object.raw <- read.fastsimcoal2.genind.for.analysis(vcf.file.path = current.pca.metadata$full.path,
                                                               n.per.population = analysis.sample.n.per.population,
                                                               sampling.group = current.pca.metadata$simulation.group,
                                                               sampling.seed = analysis.sample.random.seed) #read selected VCF, assign populations, and optionally subset individuals per population


## Calculate Weir and Cockerham Fst before PCA processing
pca.fst.value <- calculate.overall.Fst(genind.object = pca.genind.object.raw) #calculate overall Weir and Cockerham Fst
cat("Weir and Cockerham Fst between pop1 and pop2:", round(pca.fst.value, 4), "\n")


## Process SNP matrix for PCA
pca.snp.matrix <- process.SNP.data.SOM(genind.input = pca.genind.object.raw, verbose = FALSE) #process SNP data for PCA
pca.snp.matrix <- as.matrix(pca.snp.matrix) #coerce to matrix
storage.mode(pca.snp.matrix) <- "numeric" #ensure numeric matrix
if (is.null(rownames(pca.snp.matrix))) rownames(pca.snp.matrix) <- adegenet::indNames(pca.genind.object.raw)


## Remove invalid loci
valid.pca.loci <- apply(pca.snp.matrix, 2, function(locus.values) {
  non.missing.values <- locus.values[!is.na(locus.values)]
  length(non.missing.values) > 1 && stats::var(non.missing.values) > 0
}) #identify loci with nonzero variance
pca.snp.matrix <- pca.snp.matrix[, valid.pca.loci, drop = FALSE] #keep valid loci only
if (ncol(pca.snp.matrix) < 2) stop("Fewer than two valid loci available for PCA")


## Impute remaining missing values for PCA
if (any(is.na(pca.snp.matrix))) {
  pca.snp.matrix <- apply(pca.snp.matrix, 2, function(locus.values) {
    locus.values[is.na(locus.values)] <- mean(locus.values, na.rm = TRUE)
    locus.values
  }) #replace remaining missing values with locus means
  pca.snp.matrix <- as.matrix(pca.snp.matrix)
}


## Run PCA and create plot table
pca.object <- stats::prcomp(pca.snp.matrix, center = TRUE, scale. = FALSE) #run PCA
pca.variance <- pca.object$sdev^2 / sum(pca.object$sdev^2) #calculate variance explained
pca.plot.table <- data.frame(individual = rownames(pca.object$x),
                             PC1 = pca.object$x[, 1],
                             PC2 = pca.object$x[, 2],
                             population = as.character(adegenet::pop(pca.genind.object.raw))[match(rownames(pca.object$x), adegenet::indNames(pca.genind.object.raw))],
                             file = current.pca.metadata$file,
                             mig.index = current.pca.metadata$mig.index,
                             tdiv = current.pca.metadata$tdiv,
                             mig.tag = current.pca.metadata$mig.tag,
                             mig = current.pca.metadata$mig,
                             fst.pop1.pop2 = pca.fst.value,
                             stringsAsFactors = FALSE) #create PCA plotting table


## Plot PCA with population groupings
pca.plot <- ggplot2::ggplot(pca.plot.table, ggplot2::aes(x = PC1, y = PC2, color = population, group = population)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(title = paste0("PCA - ", current.pca.metadata$file, " | tdiv = ", round(current.pca.metadata$tdiv), " generations | mig = ", current.pca.metadata$mig, " | Fst = ", round(pca.fst.value, 4)),
                x = paste0("PC1 (", round(pca.variance[1] * 100, 2), "%)"),
                y = paste0("PC2 (", round(pca.variance[2] * 100, 2), "%)"),
                color = "Population") +
  ggplot2::theme_classic()
pca.plot


## Print selected PCA summary
cat("Divergence:", round(current.pca.metadata$tdiv / 1000000, 2), "million generations\n")
cat("Migration rate:", current.pca.metadata$mig, "\n")
cat("Weir and Cockerham Fst:", round(pca.fst.value, 4), "\n")




#### Run workflow for VCF files - asymmetric_16 #################################
SOM_grid_multiplier <- 4
analysis.sample.n.per.population <- analysis.sample.n.by.population.asymmetric_16 #set individuals per population for asymmetric_16
analysis.set.name <- "asymmetric_16" #define analysis set name
results.directory <- file.path(results.root.directory, analysis.set.name) #define output directory for this analysis set
dir.create(results.directory, recursive = TRUE, showWarnings = FALSE) #create output directory if needed
intermediate.results.directory <- file.path(results.directory, "intermediate_results") #define intermediate results directory for this analysis set
dir.create(intermediate.results.directory, recursive = TRUE, showWarnings = FALSE) #create intermediate results directory if needed
cat("\n\nRunning analysis set:", analysis.set.name, "\n")
cat("Individuals per population:", analysis.sample.n.per.population, "\n")
cat("Results directory:", results.directory, "\n")
cat("Random sampling base seed:", analysis.sample.random.seed, "\n")


## List VCF files
vcf.file.pattern <- "^sim[0-9]+_tdiv[0-9.]+_mig.+\\.vcf$" #define VCF filename pattern
vcf.file.paths <- list.files(vcf.directory, pattern = vcf.file.pattern, full.names = TRUE) #list VCF files only
if (length(vcf.file.paths) == 0) stop("No VCF files found")
vcf.file.names <- basename(vcf.file.paths) #extract file names


## Parse filename metadata
vcf.metadata.list <- lapply(vcf.file.names, parse.vcf.filename.metadata) #parse metadata
vcf.metadata.table <- do.call(rbind, vcf.metadata.list) #combine metadata rows
vcf.metadata.table$full.path <- vcf.file.paths #add full path


## Order metadata table
vcf.metadata.table <- vcf.metadata.table[order(vcf.metadata.table$mig.index,
                                               vcf.metadata.table$tdiv,
                                               vcf.metadata.table$mig), ] #order rows
rownames(vcf.metadata.table) <- NULL
vcf.metadata.table$simulation.group <- match(vcf.metadata.table$tdiv, sort(unique(vcf.metadata.table$tdiv))) #define matched sampling group across migration rates for each divergence time
head(vcf.metadata.table) #show metadata for files


## Create empty lists to store results
result.summary.list <- vector("list", nrow(vcf.metadata.table)) #initialize result summary list
optim.k.summary.list <- vector("list", nrow(vcf.metadata.table)) #initialize optim_k_summary list


## Loop over VCF files
for (vcf.index in seq_len(nrow(vcf.metadata.table))) {
  current.metadata <- vcf.metadata.table[vcf.index, ] #extract current metadata
  cat("\nProcessing", vcf.index, "of", nrow(vcf.metadata.table), ":", current.metadata$file, "\n")
  if (!override && check.vcf.intermediate.results.available(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory)) {
    cat("Intermediate results already found for", current.metadata$file, "- loading saved results\n")
    current.output <- load.vcf.intermediate.results(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory) #load saved per-VCF results
  } else {
    current.output <- tryCatch({
      run.SOM.workflow.for.one.vcf(vcf.file.path = current.metadata$full.path,
                                   vcf.file.metadata = current.metadata,
                                   sample.n.per.population = analysis.sample.n.per.population,
                                   sampling.group = current.metadata$simulation.group,
                                   sampling.seed = analysis.sample.random.seed)
    }, error = function(e) {
      result.summary.row <- data.frame(file = current.metadata$file,
                                       mig.index = current.metadata$mig.index,
                                       tdiv = current.metadata$tdiv,
                                       mig.tag = current.metadata$mig.tag,
                                       mig = current.metadata$mig,
                                       n.ind.raw = NA_integer_,
                                       n.loc.raw = NA_integer_,
                                       n.samples.processed = NA_integer_,
                                       n.loci.processed = NA_integer_,
                                       proportion.k2 = NA_real_,
                                       status = "error",
                                       error = conditionMessage(e),
                                       stringsAsFactors = FALSE)
      optim.k.summary.table <- data.frame(file = current.metadata$file,
                                          mig.index = current.metadata$mig.index,
                                          tdiv = current.metadata$tdiv,
                                          mig.tag = current.metadata$mig.tag,
                                          mig = current.metadata$mig,
                                          n.ind.raw = NA_integer_,
                                          n.loc.raw = NA_integer_,
                                          n.samples.processed = NA_integer_,
                                          n.loci.processed = NA_integer_,
                                          status = "error",
                                          error = conditionMessage(e),
                                          Count = NA_real_,
                                          Proportion = NA_real_,
                                          k.label = NA_character_,
                                          optim.k.row = NA_integer_,
                                          stringsAsFactors = FALSE)
      list(result.summary.row = result.summary.row, optim.k.summary.table = optim.k.summary.table)
    })
    save.vcf.intermediate.results(result.summary.row = current.output$result.summary.row,
                                  optim.k.summary.table = current.output$optim.k.summary.table,
                                  vcf.file.metadata = current.metadata,
                                  intermediate.results.directory = intermediate.results.directory) #save per-VCF results immediately
  }
  result.summary.list[[vcf.index]] <- current.output$result.summary.row #store result summary row
  optim.k.summary.list[[vcf.index]] <- current.output$optim.k.summary.table #store optim_k_summary table
  stop.if.status.not.ok(status.value = current.output$result.summary.row$status,
                        error.value = current.output$result.summary.row$error,
                        file.value = current.metadata$file,
                        method.label = "SOM") #stop immediately if SOM failed
}


## Combine SOM results
result.table <- do.call(rbind, result.summary.list) #combine result summary rows
rownames(result.table) <- NULL
optim.k.result.table <- do.call(rbind, optim.k.summary.list) #combine optim_k_summary rows
rownames(optim.k.result.table) <- NULL


## Calculate Fst and append to SOM result tables
Fst.output <- calculate.and.append.Fst(vcf.metadata.table = vcf.metadata.table,
                                       result.table = result.table,
                                       optim.k.result.table = optim.k.result.table,
                                       intermediate.results.directory = intermediate.results.directory,
                                       override = override,
                                       sample.n.per.population = analysis.sample.n.per.population,
                                       sampling.seed = analysis.sample.random.seed) #calculate and append Fst
result.table <- Fst.output$result.table #update result table with Fst
optim.k.result.table <- Fst.output$optim.k.result.table #update optim_k result table with Fst
fst.table <- Fst.output$fst.table #extract Fst table


## Calculate de novo k-means/BIC and append to result tables
deNovo.kmeans.output <- calculate.and.append.deNovo.kmeans(vcf.metadata.table = vcf.metadata.table,
                                                           result.table = result.table,
                                                           optim.k.result.table = optim.k.result.table,
                                                           intermediate.results.directory = intermediate.results.directory,
                                                           override = override,
                                                           BIC.thresh = deNovo.kmeans.BIC.thresh,
                                                           n.iter = deNovo.kmeans.n.iter,
                                                           n.start = deNovo.kmeans.n.start,
                                                           max.n.clust = deNovo.kmeans.max.n.clust,
                                                           max.proportion.PCs = deNovo.kmeans.max.proportion.PCs,
                                                           center = deNovo.kmeans.center,
                                                           scale = deNovo.kmeans.scale,
                                                           sample.n.per.population = analysis.sample.n.per.population,
                                                           sampling.seed = analysis.sample.random.seed) #calculate de novo k-means/BIC
result.table <- deNovo.kmeans.output$result.table #update result table with de novo k-means/BIC
optim.k.result.table <- deNovo.kmeans.output$optim.k.result.table #update optim_k result table with de novo k-means/BIC
deNovo.kmeans.table <- deNovo.kmeans.output$deNovo.kmeans.table #extract de novo k-means/BIC table


## Calculate sNMF and append to result tables
sNMF.output <- calculate.and.append.sNMF(vcf.metadata.table = vcf.metadata.table,
                                         result.table = result.table,
                                         optim.k.result.table = optim.k.result.table,
                                         intermediate.results.directory = intermediate.results.directory,
                                         override = override,
                                         K.values = sNMF.K.values,
                                         repetitions = sNMF.repetitions,
                                         CPU = 1,
                                         ploidy = sNMF.ploidy,
                                         cross.entropy.thresh = sNMF.cross.entropy.thresh,
                                         seed = sNMF.seed,
                                         sample.n.per.population = analysis.sample.n.per.population,
                                         sampling.seed = analysis.sample.random.seed) #calculate sNMF
result.table <- sNMF.output$result.table #update result table with sNMF
optim.k.result.table <- sNMF.output$optim.k.result.table #update optim_k result table with sNMF
sNMF.table <- sNMF.output$sNMF.table #extract sNMF table


## Save result tables
write.csv(result.table, file.path(results.directory, "fastsimcoal2_SOM_results.csv"), row.names = FALSE) #save main SOM result summary
write.csv(optim.k.result.table, file.path(results.directory, "fastsimcoal2_SOM_optim_k_summary.csv"), row.names = FALSE) #save full optim_k_summary table
write.csv(fst.table, file.path(results.directory, "fastsimcoal2_Weir_Cockerham_Fst.csv"), row.names = FALSE) #save Fst table
write.csv(deNovo.kmeans.table, file.path(results.directory, "fastsimcoal2_de_novo_kmeans_BIC.csv"), row.names = FALSE) #save de novo k-means/BIC table
write.csv(sNMF.table, file.path(results.directory, "fastsimcoal2_sNMF_results.csv"), row.names = FALSE) #save sNMF table


## Inspect results
head(result.table)
optim.k.result.table
head(fst.table)
head(deNovo.kmeans.table)
head(sNMF.table)
table(result.table$status)
table(result.table$fst.status)
table(result.table$deNovo.kmeans.status)
table(result.table$deNovo.kmeans.best.k, useNA = "ifany")
table(result.table$sNMF.status)


## Plot fitted K2 proportions for SOM, de novo k-means/BIC, and sNMF
plot.result.table <- result.table[result.table$status == "ok", ] #subset successful rows
plot.result.table <- plot.result.table[order(plot.result.table$mig, plot.result.table$tdiv), ] #order rows for plotting
if (nrow(plot.result.table) == 0) stop("No successful results available for fitted plotting")
mig.level.table <- unique(plot.result.table[order(plot.result.table$mig), c("mig.tag", "mig")]) #define unique ordered migration-rate levels
mig.level.table$mig.tag <- as.character(mig.level.table$mig.tag) #force migration-rate labels to discrete character values
plot.result.table$mig.tag <- factor(as.character(plot.result.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels by numeric migration rate
mig.colors <- setNames(viridisLite::viridis(nrow(mig.level.table), option = "D"), mig.level.table$mig.tag) #define discrete ordered colors


## Create SOM binomial count table
SOM.optim.k.count.table <- optim.k.result.table[!is.na(optim.k.result.table$Count) & optim.k.result.table$file %in% plot.result.table$file, ] #keep rows with SOM replicate counts
SOM.optim.k.count.table$Count <- as.numeric(as.character(SOM.optim.k.count.table$Count)) #force counts to numeric
SOM.total.count.table <- stats::aggregate(Count ~ file + mig.tag + mig + tdiv, data = SOM.optim.k.count.table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE)) #total retained SOM replicates per VCF
colnames(SOM.total.count.table)[colnames(SOM.total.count.table) == "Count"] <- "n.total"
SOM.k2.count.table <- SOM.optim.k.count.table[SOM.optim.k.count.table$k.label == "k2", ] #keep k2 count rows
if (nrow(SOM.k2.count.table) > 0) {
  SOM.k2.count.table <- stats::aggregate(Count ~ file + mig.tag + mig + tdiv, data = SOM.k2.count.table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE)) #number of retained SOM replicates choosing k = 2
  colnames(SOM.k2.count.table)[colnames(SOM.k2.count.table) == "Count"] <- "n.k2"
} else {
  SOM.k2.count.table <- SOM.total.count.table[, c("file", "mig.tag", "mig", "tdiv")] #create empty k2 count table if no k2 rows exist
  SOM.k2.count.table$n.k2 <- 0
}
SOM.binomial.table <- merge(SOM.total.count.table, SOM.k2.count.table, by = c("file", "mig.tag", "mig", "tdiv"), all.x = TRUE) #merge total and k2 counts
SOM.binomial.table$n.k2[is.na(SOM.binomial.table$n.k2)] <- 0 #set missing k2 rows to zero
SOM.binomial.table$n.not.k2 <- SOM.binomial.table$n.total - SOM.binomial.table$n.k2 #number not choosing k = 2
SOM.binomial.table$proportion.k2 <- SOM.binomial.table$n.k2 / SOM.binomial.table$n.total #calculate K2 proportion
SOM.binomial.table$method <- "SOM" #add method label
SOM.binomial.table <- SOM.binomial.table[, c("method", "file", "mig.tag", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")] #order columns


## Create de novo k-means/BIC binomial table
deNovo.kmeans.binomial.source.table <- result.table[result.table$status == "ok" & !is.na(result.table$deNovo.kmeans.best.k), ] #subset valid de novo k-means/BIC rows
if ("deNovo.kmeans.status" %in% colnames(deNovo.kmeans.binomial.source.table)) deNovo.kmeans.binomial.source.table <- deNovo.kmeans.binomial.source.table[deNovo.kmeans.binomial.source.table$deNovo.kmeans.status == "ok", ] #keep successful de novo k-means/BIC rows
deNovo.kmeans.best.k <- as.integer(as.character(deNovo.kmeans.binomial.source.table$deNovo.kmeans.best.k)) #extract de novo k-means/BIC best K
deNovo.kmeans.binomial.table <- data.frame(method = "deNovo.kmeans.BIC",
                                           file = deNovo.kmeans.binomial.source.table$file,
                                           mig.tag = deNovo.kmeans.binomial.source.table$mig.tag,
                                           mig = deNovo.kmeans.binomial.source.table$mig,
                                           tdiv = deNovo.kmeans.binomial.source.table$tdiv,
                                           n.k2 = as.integer(deNovo.kmeans.best.k == 2L),
                                           n.not.k2 = as.integer(deNovo.kmeans.best.k != 2L),
                                           proportion.k2 = as.numeric(deNovo.kmeans.best.k == 2L),
                                           stringsAsFactors = FALSE) #create binary de novo k-means/BIC table


## Create sNMF binomial table
sNMF.binomial.source.table <- result.table[result.table$status == "ok" & !is.na(result.table$sNMF.best.k), ] #subset valid sNMF rows
if ("sNMF.status" %in% colnames(sNMF.binomial.source.table)) sNMF.binomial.source.table <- sNMF.binomial.source.table[sNMF.binomial.source.table$sNMF.status == "ok", ] #keep successful sNMF rows
sNMF.best.k <- as.integer(as.character(sNMF.binomial.source.table$sNMF.best.k)) #extract sNMF best K
sNMF.binomial.table <- data.frame(method = "sNMF",
                                  file = sNMF.binomial.source.table$file,
                                  mig.tag = sNMF.binomial.source.table$mig.tag,
                                  mig = sNMF.binomial.source.table$mig,
                                  tdiv = sNMF.binomial.source.table$tdiv,
                                  n.k2 = as.integer(sNMF.best.k == 2L),
                                  n.not.k2 = as.integer(sNMF.best.k != 2L),
                                  proportion.k2 = as.numeric(sNMF.best.k == 2L),
                                  stringsAsFactors = FALSE) #create binary sNMF table


## Create STRUCTURE binomial table
STRUCTURE.results <- create.STRUCTURE.mean.lnprob.binomial.table(analysis.set.name = analysis.set.name,
                                                                 structure.loglik.directory = STRUCTURE.loglik.directory,
                                                                 results.directory = results.directory,
                                                                 expected.structure.replicates = expected.structure.replicates,
                                                                 STRUCTURE.mean.lnprob.delta.threshold = STRUCTURE.mean.lnprob.delta.threshold) #create STRUCTURE mean-lnP best-K table
STRUCTURE.binomial.table <- STRUCTURE.results$STRUCTURE.binomial.table #extract STRUCTURE binomial table
STRUCTURE.binomial.table$mig.tag <- vapply(STRUCTURE.binomial.table$mig, function(current.mig) {
  matching.mig.index <- which(abs(mig.level.table$mig - current.mig) < sqrt(.Machine$double.eps))
  if (length(matching.mig.index) != 1) return(NA_character_)
  return(as.character(mig.level.table$mig.tag[matching.mig.index]))
}, character(1)) #standardize STRUCTURE migration labels to match VCF migration labels
if (any(is.na(STRUCTURE.binomial.table$mig.tag))) stop("Some STRUCTURE migration rates could not be matched to VCF migration labels")
STRUCTURE.mean.lnprob.by.vcf.table <- STRUCTURE.results$STRUCTURE.mean.lnprob.by.vcf.table #extract STRUCTURE VCF-level table
STRUCTURE.rep.comparison.table <- STRUCTURE.results$STRUCTURE.rep.comparison.table #extract STRUCTURE replicate comparison table


## Combine method-specific binomial tables
fitted.binomial.input.table <- rbind(SOM.binomial.table, deNovo.kmeans.binomial.table, sNMF.binomial.table, STRUCTURE.binomial.table) #combine all methods
method.levels <- c("SOM", "deNovo.kmeans.BIC", "sNMF", "STRUCTURE") #define method order
fitted.binomial.input.table$method <- factor(as.character(fitted.binomial.input.table$method), levels = method.levels) #order methods
fitted.binomial.input.table$mig.tag <- factor(as.character(fitted.binomial.input.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels


## Fit binomial GLMs separately for each method and migration rate
fitted.GLM.prediction.list <- list() #initialize prediction list
fitted.GLM.summary.list <- list() #initialize summary list
fitted.GLM.list.index <- 1 #initialize list index
for (current.method in method.levels) {
  for (current.mig.tag in mig.level.table$mig.tag) {
    current.binomial.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == current.method & fitted.binomial.input.table$mig.tag == current.mig.tag, ] #subset one method and migration rate
    current.binomial.table <- current.binomial.table[is.finite(current.binomial.table$tdiv), ] #keep rows with finite divergence time
    if (nrow(current.binomial.table) == 0) next
    current.binomial.table <- current.binomial.table[order(current.binomial.table$tdiv), ] #order rows by divergence time
    prediction.tdiv <- seq(min(current.binomial.table$tdiv), max(current.binomial.table$tdiv), length.out = plot.fitted.prediction.n) #create prediction grid
    current.warnings <- character(0) #store GLM warnings
    if (all(current.binomial.table$n.k2 == 0)) {
      current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = 0, stringsAsFactors = FALSE) #flat zero prediction when K2 is never selected
      current.summary.row <- data.frame(method = current.method,
                                        mig.tag = as.character(current.mig.tag),
                                        mig = current.binomial.table$mig[1],
                                        n.vcf = nrow(current.binomial.table),
                                        n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                        n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                        min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                        n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                        n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                        all.zero = TRUE,
                                        all.one = FALSE,
                                        glm.fitted = FALSE,
                                        glm.converged = NA,
                                        glm.intercept = NA_real_,
                                        glm.tdiv.slope = NA_real_,
                                        tdiv.at.threshold.k2 = NA_real_,
                                        tdiv.at.threshold.k2.within.range = FALSE,
                                        glm.warning = NA_character_,
                                        glm.error = NA_character_,
                                        stringsAsFactors = FALSE) #store summary row
    } else if (all(current.binomial.table$n.not.k2 == 0)) {
      current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = 1, stringsAsFactors = FALSE) #flat one prediction when K2 is always selected
      current.summary.row <- data.frame(method = current.method,
                                        mig.tag = as.character(current.mig.tag),
                                        mig = current.binomial.table$mig[1],
                                        n.vcf = nrow(current.binomial.table),
                                        n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                        n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                        min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                        n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                        n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                        all.zero = FALSE,
                                        all.one = TRUE,
                                        glm.fitted = FALSE,
                                        glm.converged = NA,
                                        glm.intercept = NA_real_,
                                        glm.tdiv.slope = NA_real_,
                                        tdiv.at.threshold.k2 = NA_real_,
                                        tdiv.at.threshold.k2.within.range = FALSE,
                                        glm.warning = NA_character_,
                                        glm.error = NA_character_,
                                        stringsAsFactors = FALSE) #store summary row
    } else {
      current.fit <- withCallingHandlers(tryCatch(stats::glm(cbind(n.k2, n.not.k2) ~ tdiv, data = current.binomial.table, family = stats::binomial(link = "logit"), control = stats::glm.control(maxit = 100)), error = function(e) e),
                                         warning = function(w) {
                                           current.warnings <<- c(current.warnings, conditionMessage(w))
                                           invokeRestart("muffleWarning")
                                         }) #fit count-based binomial GLM and store warnings
      fit.ok <- inherits(current.fit, "glm") #check whether GLM fit succeeded
      if (fit.ok) {
        current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = stats::predict(current.fit, newdata = data.frame(tdiv = prediction.tdiv), type = "response"), stringsAsFactors = FALSE) #predict fitted K2 probability
        current.coefficients <- stats::coef(current.fit) #extract coefficients
        current.glm.intercept <- unname(current.coefficients["(Intercept)"]) #extract intercept
        current.glm.tdiv.slope <- unname(current.coefficients["tdiv"]) #extract tdiv slope
        if (is.finite(current.glm.intercept) && is.finite(current.glm.tdiv.slope) && current.glm.tdiv.slope != 0) {
          current.tdiv.at.threshold.k2 <- (stats::qlogis(plot.k2.threshold.proportion) - current.glm.intercept) / current.glm.tdiv.slope #solve logit(p) for threshold
        } else {
          current.tdiv.at.threshold.k2 <- NA_real_
        }
        current.tdiv.at.threshold.k2.within.range <- is.finite(current.tdiv.at.threshold.k2) && current.tdiv.at.threshold.k2 >= min(current.binomial.table$tdiv, na.rm = TRUE) && current.tdiv.at.threshold.k2 <= max(current.binomial.table$tdiv, na.rm = TRUE) #check whether threshold is inside simulated range
        current.summary.row <- data.frame(method = current.method,
                                          mig.tag = as.character(current.mig.tag),
                                          mig = current.binomial.table$mig[1],
                                          n.vcf = nrow(current.binomial.table),
                                          n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                          n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                          min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                          n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                          n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                          all.zero = all(current.binomial.table$proportion.k2 == 0),
                                          all.one = all(current.binomial.table$proportion.k2 == 1),
                                          glm.fitted = TRUE,
                                          glm.converged = isTRUE(current.fit$converged),
                                          glm.intercept = current.glm.intercept,
                                          glm.tdiv.slope = current.glm.tdiv.slope,
                                          tdiv.at.threshold.k2 = current.tdiv.at.threshold.k2,
                                          tdiv.at.threshold.k2.within.range = current.tdiv.at.threshold.k2.within.range,
                                          glm.warning = if (length(current.warnings) == 0) NA_character_ else paste(unique(current.warnings), collapse = " | "),
                                          glm.error = NA_character_,
                                          stringsAsFactors = FALSE) #store summary row
      } else {
        current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = NA_real_, stringsAsFactors = FALSE) #return missing predictions if GLM failed
        current.summary.row <- data.frame(method = current.method,
                                          mig.tag = as.character(current.mig.tag),
                                          mig = current.binomial.table$mig[1],
                                          n.vcf = nrow(current.binomial.table),
                                          n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                          n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                          min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                          n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                          n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                          all.zero = all(current.binomial.table$proportion.k2 == 0),
                                          all.one = all(current.binomial.table$proportion.k2 == 1),
                                          glm.fitted = FALSE,
                                          glm.converged = FALSE,
                                          glm.intercept = NA_real_,
                                          glm.tdiv.slope = NA_real_,
                                          tdiv.at.threshold.k2 = NA_real_,
                                          tdiv.at.threshold.k2.within.range = FALSE,
                                          glm.warning = if (length(current.warnings) == 0) NA_character_ else paste(unique(current.warnings), collapse = " | "),
                                          glm.error = conditionMessage(current.fit),
                                          stringsAsFactors = FALSE) #store summary row
      }
    }
    fitted.GLM.prediction.list[[fitted.GLM.list.index]] <- current.prediction.table #store prediction table
    fitted.GLM.summary.list[[fitted.GLM.list.index]] <- current.summary.row #store summary row
    fitted.GLM.list.index <- fitted.GLM.list.index + 1 #advance list index
  }
}


## Combine fitted GLM outputs
fitted.GLM.prediction.table <- do.call(rbind, fitted.GLM.prediction.list) #combine prediction tables
fitted.GLM.summary.table <- do.call(rbind, fitted.GLM.summary.list) #combine summary rows
fitted.GLM.prediction.table$method <- factor(as.character(fitted.GLM.prediction.table$method), levels = method.levels) #order methods
fitted.GLM.prediction.table$mig.tag <- factor(as.character(fitted.GLM.prediction.table$mig.tag), levels = mig.level.table$mig.tag) #order migration labels
fitted.GLM.summary.table$method <- factor(as.character(fitted.GLM.summary.table$method), levels = method.levels) #order methods
fitted.GLM.summary.table$mig.tag <- factor(as.character(fitted.GLM.summary.table$mig.tag), levels = mig.level.table$mig.tag) #order migration labels
fitted.GLM.summary.table <- fitted.GLM.summary.table[order(fitted.GLM.summary.table$method, fitted.GLM.summary.table$mig), ] #order summary table
rownames(fitted.GLM.summary.table) <- NULL
fitted.GLM.threshold.table <- fitted.GLM.summary.table[fitted.GLM.summary.table$tdiv.at.threshold.k2.within.range & is.finite(fitted.GLM.summary.table$tdiv.at.threshold.k2), ] #keep valid threshold rows


## Plot SOM results
SOM.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "SOM", ] #subset SOM observations
SOM.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "SOM" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset SOM fitted predictions
SOM.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "SOM", ] #subset SOM threshold rows
SOM.k2.fitted.plot <- ggplot2::ggplot(SOM.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = SOM.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = SOM.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(SOM.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "SOM_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = SOM.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot Weir and Cockerham Fst results
fst.plot.table <- plot.result.table[is.finite(plot.result.table$fst), ] #subset rows with finite Fst values
fst.plot.table <- fst.plot.table[order(fst.plot.table$mig, fst.plot.table$tdiv), ] #order rows for plotting
fst.plot.table$mig.tag <- factor(as.character(fst.plot.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels
SOM.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "SOM" & fitted.GLM.threshold.table$tdiv.at.threshold.k2.within.range & is.finite(fitted.GLM.threshold.table$tdiv.at.threshold.k2), ] #subset valid SOM threshold rows
SOM.k2.threshold.plot.table$mig.tag <- factor(as.character(SOM.k2.threshold.plot.table$mig.tag), levels = mig.level.table$mig.tag) #order SOM threshold migration-rate labels
fst.plot <- ggplot2::ggplot(fst.plot.table, ggplot2::aes(x = tdiv, y = fst, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = Fst_plot_point_size, alpha = Fst_plot_point_alpha) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, linewidth = Fst_plot_line_width, se = FALSE) +
  ggplot2::geom_vline(data = SOM.k2.threshold.plot.table,
                      ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag),
                      linewidth = plot.threshold.line.width,
                      linetype = plot.threshold.line.type,
                      show.legend = FALSE,
                      inherit.aes = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Weir and Cockerham Fst", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(fst.plot)
ggplot2::ggsave(file.path(results.directory, "Weir_Cockerham_Fst_by_tdiv_and_mig.svg"),
                plot = fst.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot de novo k-means/BIC results
deNovo.kmeans.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "deNovo.kmeans.BIC", ] #subset de novo k-means/BIC observations
deNovo.kmeans.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "deNovo.kmeans.BIC" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset de novo k-means/BIC fitted predictions
deNovo.kmeans.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "deNovo.kmeans.BIC", ] #subset de novo k-means/BIC threshold rows
deNovo.kmeans.k2.fitted.plot <- ggplot2::ggplot(deNovo.kmeans.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = deNovo.kmeans.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = deNovo.kmeans.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(deNovo.kmeans.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "deNovo_kmeans_BIC_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = deNovo.kmeans.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot sNMF K2 results
sNMF.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "sNMF", ] #subset sNMF observations
sNMF.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "sNMF" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset sNMF fitted predictions
sNMF.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "sNMF", ] #subset sNMF threshold rows
sNMF.k2.fitted.plot <- ggplot2::ggplot(sNMF.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = sNMF.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = sNMF.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(sNMF.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "sNMF_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = sNMF.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot STRUCTURE K2 results
STRUCTURE.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "STRUCTURE", ] #subset STRUCTURE observations
STRUCTURE.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "STRUCTURE" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset STRUCTURE fitted predictions
STRUCTURE.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "STRUCTURE", ] #subset STRUCTURE threshold rows
STRUCTURE.k2.fitted.plot <- ggplot2::ggplot(STRUCTURE.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = STRUCTURE.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = STRUCTURE.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(STRUCTURE.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "STRUCTURE_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = STRUCTURE.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Inspect fitted 50% threshold results
fitted.GLM.summary.table
fitted.GLM.threshold.table



## Create PCA plot for one VCF file

## Show available migration settings
available.migration.table <- unique(vcf.metadata.table[, c("mig.tag", "mig")]) #show all available migration settings including mig = 0
available.migration.table <- available.migration.table[order(available.migration.table$mig), ]
rownames(available.migration.table) <- NULL
print(available.migration.table)
if (!(migration.vcf.mig.for.PCA %in% unique(vcf.metadata.table$mig))) stop(paste("migration.vcf.mig.for.PCA must be one of:", paste(sort(unique(vcf.metadata.table$mig)), collapse = ", ")))


## Select file number N among VCF files for selected migration rate
migration.vcf.metadata.table <- vcf.metadata.table[vcf.metadata.table$mig == migration.vcf.mig.for.PCA, ] #keep selected VCF files only
migration.vcf.metadata.table <- migration.vcf.metadata.table[order(migration.vcf.metadata.table$tdiv), ] #order selected files by divergence time
rownames(migration.vcf.metadata.table) <- NULL
if (nrow(migration.vcf.metadata.table) == 0) stop(paste("No VCF files found for migration rate", migration.vcf.mig.for.PCA))
if (migration.vcf.file.number.for.PCA < 1 || migration.vcf.file.number.for.PCA > nrow(migration.vcf.metadata.table)) stop(paste("migration.vcf.file.number.for.PCA must be between 1 and", nrow(migration.vcf.metadata.table), "for migration rate", migration.vcf.mig.for.PCA))
current.pca.metadata <- migration.vcf.metadata.table[migration.vcf.file.number.for.PCA, ] #select VCF file N for selected migration rate
cat("\nCreating PCA for VCF file", migration.vcf.file.number.for.PCA, "of", nrow(migration.vcf.metadata.table), "for migration rate", migration.vcf.mig.for.PCA, ":", current.pca.metadata$file, "\n")
cat("Migration-index tag:", current.pca.metadata$mig.index, "\n")
cat("Divergence time:", current.pca.metadata$tdiv, "generations\n")
cat("Migration tag:", current.pca.metadata$mig.tag, "\n")
cat("Migration rate:", current.pca.metadata$mig, "\n")


## Read VCF and process SNP matrix
pca.genind.object.raw <- read.fastsimcoal2.genind.for.analysis(vcf.file.path = current.pca.metadata$full.path,
                                                               n.per.population = analysis.sample.n.per.population,
                                                               sampling.group = current.pca.metadata$simulation.group,
                                                               sampling.seed = analysis.sample.random.seed) #read selected VCF, assign populations, and optionally subset individuals per population


## Calculate Weir and Cockerham Fst before PCA processing
pca.fst.value <- calculate.overall.Fst(genind.object = pca.genind.object.raw) #calculate overall Weir and Cockerham Fst
cat("Weir and Cockerham Fst between pop1 and pop2:", round(pca.fst.value, 4), "\n")


## Process SNP matrix for PCA
pca.snp.matrix <- process.SNP.data.SOM(genind.input = pca.genind.object.raw, verbose = FALSE) #process SNP data for PCA
pca.snp.matrix <- as.matrix(pca.snp.matrix) #coerce to matrix
storage.mode(pca.snp.matrix) <- "numeric" #ensure numeric matrix
if (is.null(rownames(pca.snp.matrix))) rownames(pca.snp.matrix) <- adegenet::indNames(pca.genind.object.raw)


## Remove invalid loci
valid.pca.loci <- apply(pca.snp.matrix, 2, function(locus.values) {
  non.missing.values <- locus.values[!is.na(locus.values)]
  length(non.missing.values) > 1 && stats::var(non.missing.values) > 0
}) #identify loci with nonzero variance
pca.snp.matrix <- pca.snp.matrix[, valid.pca.loci, drop = FALSE] #keep valid loci only
if (ncol(pca.snp.matrix) < 2) stop("Fewer than two valid loci available for PCA")


## Impute remaining missing values for PCA
if (any(is.na(pca.snp.matrix))) {
  pca.snp.matrix <- apply(pca.snp.matrix, 2, function(locus.values) {
    locus.values[is.na(locus.values)] <- mean(locus.values, na.rm = TRUE)
    locus.values
  }) #replace remaining missing values with locus means
  pca.snp.matrix <- as.matrix(pca.snp.matrix)
}


## Run PCA and create plot table
pca.object <- stats::prcomp(pca.snp.matrix, center = TRUE, scale. = FALSE) #run PCA
pca.variance <- pca.object$sdev^2 / sum(pca.object$sdev^2) #calculate variance explained
pca.plot.table <- data.frame(individual = rownames(pca.object$x),
                             PC1 = pca.object$x[, 1],
                             PC2 = pca.object$x[, 2],
                             population = as.character(adegenet::pop(pca.genind.object.raw))[match(rownames(pca.object$x), adegenet::indNames(pca.genind.object.raw))],
                             file = current.pca.metadata$file,
                             mig.index = current.pca.metadata$mig.index,
                             tdiv = current.pca.metadata$tdiv,
                             mig.tag = current.pca.metadata$mig.tag,
                             mig = current.pca.metadata$mig,
                             fst.pop1.pop2 = pca.fst.value,
                             stringsAsFactors = FALSE) #create PCA plotting table


## Plot PCA with population groupings
pca.plot <- ggplot2::ggplot(pca.plot.table, ggplot2::aes(x = PC1, y = PC2, color = population, group = population)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(title = paste0("PCA - ", current.pca.metadata$file, " | tdiv = ", round(current.pca.metadata$tdiv), " generations | mig = ", current.pca.metadata$mig, " | Fst = ", round(pca.fst.value, 4)),
                x = paste0("PC1 (", round(pca.variance[1] * 100, 2), "%)"),
                y = paste0("PC2 (", round(pca.variance[2] * 100, 2), "%)"),
                color = "Population") +
  ggplot2::theme_classic()
pca.plot


## Show PCA plot and table
pca.file.tag <- tools::file_path_sans_ext(current.pca.metadata$file) #create PCA file tag


## Print selected PCA summary
cat("Divergence:", round(current.pca.metadata$tdiv / 1000000, 2), "million generations\n")
cat("Migration rate:", current.pca.metadata$mig, "\n")
cat("Weir and Cockerham Fst:", round(pca.fst.value, 4), "\n")




#### Run workflow for VCF files - asymmetric_8 #################################
SOM_grid_multiplier <- 4
analysis.sample.n.per.population <- analysis.sample.n.by.population.asymmetric_8 #set individuals per population for asymmetric_8
analysis.set.name <- "asymmetric_8" #define analysis set name
results.directory <- file.path(results.root.directory, analysis.set.name) #define output directory for this analysis set
dir.create(results.directory, recursive = TRUE, showWarnings = FALSE) #create output directory if needed
intermediate.results.directory <- file.path(results.directory, "intermediate_results") #define intermediate results directory for this analysis set
dir.create(intermediate.results.directory, recursive = TRUE, showWarnings = FALSE) #create intermediate results directory if needed
cat("\n\nRunning analysis set:", analysis.set.name, "\n")
cat("Individuals per population:", analysis.sample.n.per.population, "\n")
cat("Results directory:", results.directory, "\n")
cat("Random sampling base seed:", analysis.sample.random.seed, "\n")


## List VCF files
vcf.file.pattern <- "^sim[0-9]+_tdiv[0-9.]+_mig.+\\.vcf$" #define VCF filename pattern
vcf.file.paths <- list.files(vcf.directory, pattern = vcf.file.pattern, full.names = TRUE) #list VCF files only
if (length(vcf.file.paths) == 0) stop("No VCF files found")
vcf.file.names <- basename(vcf.file.paths) #extract file names


## Parse filename metadata
vcf.metadata.list <- lapply(vcf.file.names, parse.vcf.filename.metadata) #parse metadata
vcf.metadata.table <- do.call(rbind, vcf.metadata.list) #combine metadata rows
vcf.metadata.table$full.path <- vcf.file.paths #add full path


## Order metadata table
vcf.metadata.table <- vcf.metadata.table[order(vcf.metadata.table$mig.index,
                                               vcf.metadata.table$tdiv,
                                               vcf.metadata.table$mig), ] #order rows
rownames(vcf.metadata.table) <- NULL
vcf.metadata.table$simulation.group <- match(vcf.metadata.table$tdiv, sort(unique(vcf.metadata.table$tdiv))) #define matched sampling group across migration rates for each divergence time
head(vcf.metadata.table) #show metadata for files


## Create empty lists to store results
result.summary.list <- vector("list", nrow(vcf.metadata.table)) #initialize result summary list
optim.k.summary.list <- vector("list", nrow(vcf.metadata.table)) #initialize optim_k_summary list


## Loop over VCF files
for (vcf.index in seq_len(nrow(vcf.metadata.table))) {
  current.metadata <- vcf.metadata.table[vcf.index, ] #extract current metadata
  cat("\nProcessing", vcf.index, "of", nrow(vcf.metadata.table), ":", current.metadata$file, "\n")
  if (!override && check.vcf.intermediate.results.available(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory)) {
    cat("Intermediate results already found for", current.metadata$file, "- loading saved results\n")
    current.output <- load.vcf.intermediate.results(vcf.file.metadata = current.metadata, intermediate.results.directory = intermediate.results.directory) #load saved per-VCF results
  } else {
    current.output <- tryCatch({
      run.SOM.workflow.for.one.vcf(vcf.file.path = current.metadata$full.path,
                                   vcf.file.metadata = current.metadata,
                                   sample.n.per.population = analysis.sample.n.per.population,
                                   sampling.group = current.metadata$simulation.group,
                                   sampling.seed = analysis.sample.random.seed)
    }, error = function(e) {
      result.summary.row <- data.frame(file = current.metadata$file,
                                       mig.index = current.metadata$mig.index,
                                       tdiv = current.metadata$tdiv,
                                       mig.tag = current.metadata$mig.tag,
                                       mig = current.metadata$mig,
                                       n.ind.raw = NA_integer_,
                                       n.loc.raw = NA_integer_,
                                       n.samples.processed = NA_integer_,
                                       n.loci.processed = NA_integer_,
                                       proportion.k2 = NA_real_,
                                       status = "error",
                                       error = conditionMessage(e),
                                       stringsAsFactors = FALSE)
      optim.k.summary.table <- data.frame(file = current.metadata$file,
                                          mig.index = current.metadata$mig.index,
                                          tdiv = current.metadata$tdiv,
                                          mig.tag = current.metadata$mig.tag,
                                          mig = current.metadata$mig,
                                          n.ind.raw = NA_integer_,
                                          n.loc.raw = NA_integer_,
                                          n.samples.processed = NA_integer_,
                                          n.loci.processed = NA_integer_,
                                          status = "error",
                                          error = conditionMessage(e),
                                          Count = NA_real_,
                                          Proportion = NA_real_,
                                          k.label = NA_character_,
                                          optim.k.row = NA_integer_,
                                          stringsAsFactors = FALSE)
      list(result.summary.row = result.summary.row, optim.k.summary.table = optim.k.summary.table)
    })
    save.vcf.intermediate.results(result.summary.row = current.output$result.summary.row,
                                  optim.k.summary.table = current.output$optim.k.summary.table,
                                  vcf.file.metadata = current.metadata,
                                  intermediate.results.directory = intermediate.results.directory) #save per-VCF results immediately
  }
  result.summary.list[[vcf.index]] <- current.output$result.summary.row #store result summary row
  optim.k.summary.list[[vcf.index]] <- current.output$optim.k.summary.table #store optim_k_summary table
  stop.if.status.not.ok(status.value = current.output$result.summary.row$status,
                        error.value = current.output$result.summary.row$error,
                        file.value = current.metadata$file,
                        method.label = "SOM") #stop immediately if SOM failed
}


## Combine SOM results
result.table <- do.call(rbind, result.summary.list) #combine result summary rows
rownames(result.table) <- NULL
optim.k.result.table <- do.call(rbind, optim.k.summary.list) #combine optim_k_summary rows
rownames(optim.k.result.table) <- NULL


## Calculate Fst and append to SOM result tables
Fst.output <- calculate.and.append.Fst(vcf.metadata.table = vcf.metadata.table,
                                       result.table = result.table,
                                       optim.k.result.table = optim.k.result.table,
                                       intermediate.results.directory = intermediate.results.directory,
                                       override = override,
                                       sample.n.per.population = analysis.sample.n.per.population,
                                       sampling.seed = analysis.sample.random.seed) #calculate and append Fst
result.table <- Fst.output$result.table #update result table with Fst
optim.k.result.table <- Fst.output$optim.k.result.table #update optim_k result table with Fst
fst.table <- Fst.output$fst.table #extract Fst table


## Calculate de novo k-means/BIC and append to result tables
deNovo.kmeans.output <- calculate.and.append.deNovo.kmeans(vcf.metadata.table = vcf.metadata.table,
                                                           result.table = result.table,
                                                           optim.k.result.table = optim.k.result.table,
                                                           intermediate.results.directory = intermediate.results.directory,
                                                           override = override,
                                                           BIC.thresh = deNovo.kmeans.BIC.thresh,
                                                           n.iter = deNovo.kmeans.n.iter,
                                                           n.start = deNovo.kmeans.n.start,
                                                           max.n.clust = deNovo.kmeans.max.n.clust,
                                                           max.proportion.PCs = deNovo.kmeans.max.proportion.PCs,
                                                           center = deNovo.kmeans.center,
                                                           scale = deNovo.kmeans.scale,
                                                           sample.n.per.population = analysis.sample.n.per.population,
                                                           sampling.seed = analysis.sample.random.seed) #calculate de novo k-means/BIC
result.table <- deNovo.kmeans.output$result.table #update result table with de novo k-means/BIC
optim.k.result.table <- deNovo.kmeans.output$optim.k.result.table #update optim_k result table with de novo k-means/BIC
deNovo.kmeans.table <- deNovo.kmeans.output$deNovo.kmeans.table #extract de novo k-means/BIC table


## Calculate sNMF and append to result tables
sNMF.output <- calculate.and.append.sNMF(vcf.metadata.table = vcf.metadata.table,
                                         result.table = result.table,
                                         optim.k.result.table = optim.k.result.table,
                                         intermediate.results.directory = intermediate.results.directory,
                                         override = override,
                                         K.values = sNMF.K.values,
                                         repetitions = sNMF.repetitions,
                                         CPU = 1,
                                         ploidy = sNMF.ploidy,
                                         cross.entropy.thresh = sNMF.cross.entropy.thresh,
                                         seed = sNMF.seed,
                                         sample.n.per.population = analysis.sample.n.per.population,
                                         sampling.seed = analysis.sample.random.seed) #calculate sNMF
result.table <- sNMF.output$result.table #update result table with sNMF
optim.k.result.table <- sNMF.output$optim.k.result.table #update optim_k result table with sNMF
sNMF.table <- sNMF.output$sNMF.table #extract sNMF table


## Save result tables
write.csv(result.table, file.path(results.directory, "fastsimcoal2_SOM_results.csv"), row.names = FALSE) #save main SOM result summary
write.csv(optim.k.result.table, file.path(results.directory, "fastsimcoal2_SOM_optim_k_summary.csv"), row.names = FALSE) #save full optim_k_summary table
write.csv(fst.table, file.path(results.directory, "fastsimcoal2_Weir_Cockerham_Fst.csv"), row.names = FALSE) #save Fst table
write.csv(deNovo.kmeans.table, file.path(results.directory, "fastsimcoal2_de_novo_kmeans_BIC.csv"), row.names = FALSE) #save de novo k-means/BIC table
write.csv(sNMF.table, file.path(results.directory, "fastsimcoal2_sNMF_results.csv"), row.names = FALSE) #save sNMF table


## Inspect results
head(result.table)
optim.k.result.table
head(fst.table)
head(deNovo.kmeans.table)
head(sNMF.table)
table(result.table$status)
table(result.table$fst.status)
table(result.table$deNovo.kmeans.status)
table(result.table$deNovo.kmeans.best.k, useNA = "ifany")
table(result.table$sNMF.status)


## Plot fitted K2 proportions for SOM, de novo k-means/BIC, and sNMF
plot.result.table <- result.table[result.table$status == "ok", ] #subset successful rows
plot.result.table <- plot.result.table[order(plot.result.table$mig, plot.result.table$tdiv), ] #order rows for plotting
if (nrow(plot.result.table) == 0) stop("No successful results available for fitted plotting")
mig.level.table <- unique(plot.result.table[order(plot.result.table$mig), c("mig.tag", "mig")]) #define unique ordered migration-rate levels
mig.level.table$mig.tag <- as.character(mig.level.table$mig.tag) #force migration-rate labels to discrete character values
plot.result.table$mig.tag <- factor(as.character(plot.result.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels by numeric migration rate
mig.colors <- setNames(viridisLite::viridis(nrow(mig.level.table), option = "D"), mig.level.table$mig.tag) #define discrete ordered colors


## Create SOM binomial count table
SOM.optim.k.count.table <- optim.k.result.table[!is.na(optim.k.result.table$Count) & optim.k.result.table$file %in% plot.result.table$file, ] #keep rows with SOM replicate counts
SOM.optim.k.count.table$Count <- as.numeric(as.character(SOM.optim.k.count.table$Count)) #force counts to numeric
SOM.total.count.table <- stats::aggregate(Count ~ file + mig.tag + mig + tdiv, data = SOM.optim.k.count.table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE)) #total retained SOM replicates per VCF
colnames(SOM.total.count.table)[colnames(SOM.total.count.table) == "Count"] <- "n.total"
SOM.k2.count.table <- SOM.optim.k.count.table[SOM.optim.k.count.table$k.label == "k2", ] #keep k2 count rows
if (nrow(SOM.k2.count.table) > 0) {
  SOM.k2.count.table <- stats::aggregate(Count ~ file + mig.tag + mig + tdiv, data = SOM.k2.count.table, FUN = function(x) sum(as.numeric(x), na.rm = TRUE)) #number of retained SOM replicates choosing k = 2
  colnames(SOM.k2.count.table)[colnames(SOM.k2.count.table) == "Count"] <- "n.k2"
} else {
  SOM.k2.count.table <- SOM.total.count.table[, c("file", "mig.tag", "mig", "tdiv")] #create empty k2 count table if no k2 rows exist
  SOM.k2.count.table$n.k2 <- 0
}
SOM.binomial.table <- merge(SOM.total.count.table, SOM.k2.count.table, by = c("file", "mig.tag", "mig", "tdiv"), all.x = TRUE) #merge total and k2 counts
SOM.binomial.table$n.k2[is.na(SOM.binomial.table$n.k2)] <- 0 #set missing k2 rows to zero
SOM.binomial.table$n.not.k2 <- SOM.binomial.table$n.total - SOM.binomial.table$n.k2 #number not choosing k = 2
SOM.binomial.table$proportion.k2 <- SOM.binomial.table$n.k2 / SOM.binomial.table$n.total #calculate K2 proportion
SOM.binomial.table$method <- "SOM" #add method label
SOM.binomial.table <- SOM.binomial.table[, c("method", "file", "mig.tag", "mig", "tdiv", "n.k2", "n.not.k2", "proportion.k2")] #order columns


## Create de novo k-means/BIC binomial table
deNovo.kmeans.binomial.source.table <- result.table[result.table$status == "ok" & !is.na(result.table$deNovo.kmeans.best.k), ] #subset valid de novo k-means/BIC rows
if ("deNovo.kmeans.status" %in% colnames(deNovo.kmeans.binomial.source.table)) deNovo.kmeans.binomial.source.table <- deNovo.kmeans.binomial.source.table[deNovo.kmeans.binomial.source.table$deNovo.kmeans.status == "ok", ] #keep successful de novo k-means/BIC rows
deNovo.kmeans.best.k <- as.integer(as.character(deNovo.kmeans.binomial.source.table$deNovo.kmeans.best.k)) #extract de novo k-means/BIC best K
deNovo.kmeans.binomial.table <- data.frame(method = "deNovo.kmeans.BIC",
                                           file = deNovo.kmeans.binomial.source.table$file,
                                           mig.tag = deNovo.kmeans.binomial.source.table$mig.tag,
                                           mig = deNovo.kmeans.binomial.source.table$mig,
                                           tdiv = deNovo.kmeans.binomial.source.table$tdiv,
                                           n.k2 = as.integer(deNovo.kmeans.best.k == 2L),
                                           n.not.k2 = as.integer(deNovo.kmeans.best.k != 2L),
                                           proportion.k2 = as.numeric(deNovo.kmeans.best.k == 2L),
                                           stringsAsFactors = FALSE) #create binary de novo k-means/BIC table


## Create sNMF binomial table
sNMF.binomial.source.table <- result.table[result.table$status == "ok" & !is.na(result.table$sNMF.best.k), ] #subset valid sNMF rows
if ("sNMF.status" %in% colnames(sNMF.binomial.source.table)) sNMF.binomial.source.table <- sNMF.binomial.source.table[sNMF.binomial.source.table$sNMF.status == "ok", ] #keep successful sNMF rows
sNMF.best.k <- as.integer(as.character(sNMF.binomial.source.table$sNMF.best.k)) #extract sNMF best K
sNMF.binomial.table <- data.frame(method = "sNMF",
                                  file = sNMF.binomial.source.table$file,
                                  mig.tag = sNMF.binomial.source.table$mig.tag,
                                  mig = sNMF.binomial.source.table$mig,
                                  tdiv = sNMF.binomial.source.table$tdiv,
                                  n.k2 = as.integer(sNMF.best.k == 2L),
                                  n.not.k2 = as.integer(sNMF.best.k != 2L),
                                  proportion.k2 = as.numeric(sNMF.best.k == 2L),
                                  stringsAsFactors = FALSE) #create binary sNMF table


## Create STRUCTURE binomial table
STRUCTURE.results <- create.STRUCTURE.mean.lnprob.binomial.table(analysis.set.name = analysis.set.name,
                                                                 structure.loglik.directory = STRUCTURE.loglik.directory,
                                                                 results.directory = results.directory,
                                                                 expected.structure.replicates = expected.structure.replicates,
                                                                 STRUCTURE.mean.lnprob.delta.threshold = STRUCTURE.mean.lnprob.delta.threshold) #create STRUCTURE mean-lnP best-K table
STRUCTURE.binomial.table <- STRUCTURE.results$STRUCTURE.binomial.table #extract STRUCTURE binomial table
STRUCTURE.binomial.table$mig.tag <- vapply(STRUCTURE.binomial.table$mig, function(current.mig) {
  matching.mig.index <- which(abs(mig.level.table$mig - current.mig) < sqrt(.Machine$double.eps))
  if (length(matching.mig.index) != 1) return(NA_character_)
  return(as.character(mig.level.table$mig.tag[matching.mig.index]))
}, character(1)) #standardize STRUCTURE migration labels to match VCF migration labels
if (any(is.na(STRUCTURE.binomial.table$mig.tag))) stop("Some STRUCTURE migration rates could not be matched to VCF migration labels")
STRUCTURE.mean.lnprob.by.vcf.table <- STRUCTURE.results$STRUCTURE.mean.lnprob.by.vcf.table #extract STRUCTURE VCF-level table
STRUCTURE.rep.comparison.table <- STRUCTURE.results$STRUCTURE.rep.comparison.table #extract STRUCTURE replicate comparison table


## Combine method-specific binomial tables
fitted.binomial.input.table <- rbind(SOM.binomial.table, deNovo.kmeans.binomial.table, sNMF.binomial.table, STRUCTURE.binomial.table) #combine all methods
method.levels <- c("SOM", "deNovo.kmeans.BIC", "sNMF", "STRUCTURE") #define method order
fitted.binomial.input.table$method <- factor(as.character(fitted.binomial.input.table$method), levels = method.levels) #order methods
fitted.binomial.input.table$mig.tag <- factor(as.character(fitted.binomial.input.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels


## Fit binomial GLMs separately for each method and migration rate
fitted.GLM.prediction.list <- list() #initialize prediction list
fitted.GLM.summary.list <- list() #initialize summary list
fitted.GLM.list.index <- 1 #initialize list index
for (current.method in method.levels) {
  for (current.mig.tag in mig.level.table$mig.tag) {
    current.binomial.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == current.method & fitted.binomial.input.table$mig.tag == current.mig.tag, ] #subset one method and migration rate
    current.binomial.table <- current.binomial.table[is.finite(current.binomial.table$tdiv), ] #keep rows with finite divergence time
    if (nrow(current.binomial.table) == 0) next
    current.binomial.table <- current.binomial.table[order(current.binomial.table$tdiv), ] #order rows by divergence time
    prediction.tdiv <- seq(min(current.binomial.table$tdiv), max(current.binomial.table$tdiv), length.out = plot.fitted.prediction.n) #create prediction grid
    current.warnings <- character(0) #store GLM warnings
    if (all(current.binomial.table$n.k2 == 0)) {
      current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = 0, stringsAsFactors = FALSE) #flat zero prediction when K2 is never selected
      current.summary.row <- data.frame(method = current.method,
                                        mig.tag = as.character(current.mig.tag),
                                        mig = current.binomial.table$mig[1],
                                        n.vcf = nrow(current.binomial.table),
                                        n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                        n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                        min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                        n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                        n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                        all.zero = TRUE,
                                        all.one = FALSE,
                                        glm.fitted = FALSE,
                                        glm.converged = NA,
                                        glm.intercept = NA_real_,
                                        glm.tdiv.slope = NA_real_,
                                        tdiv.at.threshold.k2 = NA_real_,
                                        tdiv.at.threshold.k2.within.range = FALSE,
                                        glm.warning = NA_character_,
                                        glm.error = NA_character_,
                                        stringsAsFactors = FALSE) #store summary row
    } else if (all(current.binomial.table$n.not.k2 == 0)) {
      current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = 1, stringsAsFactors = FALSE) #flat one prediction when K2 is always selected
      current.summary.row <- data.frame(method = current.method,
                                        mig.tag = as.character(current.mig.tag),
                                        mig = current.binomial.table$mig[1],
                                        n.vcf = nrow(current.binomial.table),
                                        n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                        n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                        min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                        n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                        n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                        n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                        all.zero = FALSE,
                                        all.one = TRUE,
                                        glm.fitted = FALSE,
                                        glm.converged = NA,
                                        glm.intercept = NA_real_,
                                        glm.tdiv.slope = NA_real_,
                                        tdiv.at.threshold.k2 = NA_real_,
                                        tdiv.at.threshold.k2.within.range = FALSE,
                                        glm.warning = NA_character_,
                                        glm.error = NA_character_,
                                        stringsAsFactors = FALSE) #store summary row
    } else {
      current.fit <- withCallingHandlers(tryCatch(stats::glm(cbind(n.k2, n.not.k2) ~ tdiv, data = current.binomial.table, family = stats::binomial(link = "logit"), control = stats::glm.control(maxit = 100)), error = function(e) e),
                                         warning = function(w) {
                                           current.warnings <<- c(current.warnings, conditionMessage(w))
                                           invokeRestart("muffleWarning")
                                         }) #fit count-based binomial GLM and store warnings
      fit.ok <- inherits(current.fit, "glm") #check whether GLM fit succeeded
      if (fit.ok) {
        current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = stats::predict(current.fit, newdata = data.frame(tdiv = prediction.tdiv), type = "response"), stringsAsFactors = FALSE) #predict fitted K2 probability
        current.coefficients <- stats::coef(current.fit) #extract coefficients
        current.glm.intercept <- unname(current.coefficients["(Intercept)"]) #extract intercept
        current.glm.tdiv.slope <- unname(current.coefficients["tdiv"]) #extract tdiv slope
        if (is.finite(current.glm.intercept) && is.finite(current.glm.tdiv.slope) && current.glm.tdiv.slope != 0) {
          current.tdiv.at.threshold.k2 <- (stats::qlogis(plot.k2.threshold.proportion) - current.glm.intercept) / current.glm.tdiv.slope #solve logit(p) for threshold
        } else {
          current.tdiv.at.threshold.k2 <- NA_real_
        }
        current.tdiv.at.threshold.k2.within.range <- is.finite(current.tdiv.at.threshold.k2) && current.tdiv.at.threshold.k2 >= min(current.binomial.table$tdiv, na.rm = TRUE) && current.tdiv.at.threshold.k2 <= max(current.binomial.table$tdiv, na.rm = TRUE) #check whether threshold is inside simulated range
        current.summary.row <- data.frame(method = current.method,
                                          mig.tag = as.character(current.mig.tag),
                                          mig = current.binomial.table$mig[1],
                                          n.vcf = nrow(current.binomial.table),
                                          n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                          n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                          min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                          n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                          n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                          all.zero = all(current.binomial.table$proportion.k2 == 0),
                                          all.one = all(current.binomial.table$proportion.k2 == 1),
                                          glm.fitted = TRUE,
                                          glm.converged = isTRUE(current.fit$converged),
                                          glm.intercept = current.glm.intercept,
                                          glm.tdiv.slope = current.glm.tdiv.slope,
                                          tdiv.at.threshold.k2 = current.tdiv.at.threshold.k2,
                                          tdiv.at.threshold.k2.within.range = current.tdiv.at.threshold.k2.within.range,
                                          glm.warning = if (length(current.warnings) == 0) NA_character_ else paste(unique(current.warnings), collapse = " | "),
                                          glm.error = NA_character_,
                                          stringsAsFactors = FALSE) #store summary row
      } else {
        current.prediction.table <- data.frame(method = current.method, mig.tag = as.character(current.mig.tag), mig = current.binomial.table$mig[1], tdiv = prediction.tdiv, fitted.proportion.k2 = NA_real_, stringsAsFactors = FALSE) #return missing predictions if GLM failed
        current.summary.row <- data.frame(method = current.method,
                                          mig.tag = as.character(current.mig.tag),
                                          mig = current.binomial.table$mig[1],
                                          n.vcf = nrow(current.binomial.table),
                                          n.total.observations = sum(current.binomial.table$n.k2 + current.binomial.table$n.not.k2, na.rm = TRUE),
                                          n.k2 = sum(current.binomial.table$n.k2, na.rm = TRUE),
                                          min.proportion.k2 = min(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          max.proportion.k2 = max(current.binomial.table$proportion.k2, na.rm = TRUE),
                                          n.zero.proportion = sum(current.binomial.table$proportion.k2 == 0, na.rm = TRUE),
                                          n.one.proportion = sum(current.binomial.table$proportion.k2 == 1, na.rm = TRUE),
                                          n.mixed.proportion = sum(current.binomial.table$proportion.k2 > 0 & current.binomial.table$proportion.k2 < 1, na.rm = TRUE),
                                          all.zero = all(current.binomial.table$proportion.k2 == 0),
                                          all.one = all(current.binomial.table$proportion.k2 == 1),
                                          glm.fitted = FALSE,
                                          glm.converged = FALSE,
                                          glm.intercept = NA_real_,
                                          glm.tdiv.slope = NA_real_,
                                          tdiv.at.threshold.k2 = NA_real_,
                                          tdiv.at.threshold.k2.within.range = FALSE,
                                          glm.warning = if (length(current.warnings) == 0) NA_character_ else paste(unique(current.warnings), collapse = " | "),
                                          glm.error = conditionMessage(current.fit),
                                          stringsAsFactors = FALSE) #store summary row
      }
    }
    fitted.GLM.prediction.list[[fitted.GLM.list.index]] <- current.prediction.table #store prediction table
    fitted.GLM.summary.list[[fitted.GLM.list.index]] <- current.summary.row #store summary row
    fitted.GLM.list.index <- fitted.GLM.list.index + 1 #advance list index
  }
}


## Combine fitted GLM outputs
fitted.GLM.prediction.table <- do.call(rbind, fitted.GLM.prediction.list) #combine prediction tables
fitted.GLM.summary.table <- do.call(rbind, fitted.GLM.summary.list) #combine summary rows
fitted.GLM.prediction.table$method <- factor(as.character(fitted.GLM.prediction.table$method), levels = method.levels) #order methods
fitted.GLM.prediction.table$mig.tag <- factor(as.character(fitted.GLM.prediction.table$mig.tag), levels = mig.level.table$mig.tag) #order migration labels
fitted.GLM.summary.table$method <- factor(as.character(fitted.GLM.summary.table$method), levels = method.levels) #order methods
fitted.GLM.summary.table$mig.tag <- factor(as.character(fitted.GLM.summary.table$mig.tag), levels = mig.level.table$mig.tag) #order migration labels
fitted.GLM.summary.table <- fitted.GLM.summary.table[order(fitted.GLM.summary.table$method, fitted.GLM.summary.table$mig), ] #order summary table
rownames(fitted.GLM.summary.table) <- NULL
fitted.GLM.threshold.table <- fitted.GLM.summary.table[fitted.GLM.summary.table$tdiv.at.threshold.k2.within.range & is.finite(fitted.GLM.summary.table$tdiv.at.threshold.k2), ] #keep valid threshold rows


## Plot SOM results
SOM.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "SOM", ] #subset SOM observations
SOM.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "SOM" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset SOM fitted predictions
SOM.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "SOM", ] #subset SOM threshold rows
SOM.k2.fitted.plot <- ggplot2::ggplot(SOM.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = SOM.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = SOM.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(SOM.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "SOM_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = SOM.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot Weir and Cockerham Fst results
fst.plot.table <- plot.result.table[is.finite(plot.result.table$fst), ] #subset rows with finite Fst values
fst.plot.table <- fst.plot.table[order(fst.plot.table$mig, fst.plot.table$tdiv), ] #order rows for plotting
fst.plot.table$mig.tag <- factor(as.character(fst.plot.table$mig.tag), levels = mig.level.table$mig.tag) #order migration-rate labels
SOM.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "SOM" & fitted.GLM.threshold.table$tdiv.at.threshold.k2.within.range & is.finite(fitted.GLM.threshold.table$tdiv.at.threshold.k2), ] #subset valid SOM threshold rows
SOM.k2.threshold.plot.table$mig.tag <- factor(as.character(SOM.k2.threshold.plot.table$mig.tag), levels = mig.level.table$mig.tag) #order SOM threshold migration-rate labels
fst.plot <- ggplot2::ggplot(fst.plot.table, ggplot2::aes(x = tdiv, y = fst, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = Fst_plot_point_size, alpha = Fst_plot_point_alpha) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, linewidth = Fst_plot_line_width, se = FALSE) +
  ggplot2::geom_vline(data = SOM.k2.threshold.plot.table,
                      ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag),
                      linewidth = plot.threshold.line.width,
                      linetype = plot.threshold.line.type,
                      show.legend = FALSE,
                      inherit.aes = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Weir and Cockerham Fst", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(fst.plot)
ggplot2::ggsave(file.path(results.directory, "Weir_Cockerham_Fst_by_tdiv_and_mig.svg"),
                plot = fst.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot de novo k-means/BIC results
deNovo.kmeans.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "deNovo.kmeans.BIC", ] #subset de novo k-means/BIC observations
deNovo.kmeans.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "deNovo.kmeans.BIC" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset de novo k-means/BIC fitted predictions
deNovo.kmeans.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "deNovo.kmeans.BIC", ] #subset de novo k-means/BIC threshold rows
deNovo.kmeans.k2.fitted.plot <- ggplot2::ggplot(deNovo.kmeans.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = deNovo.kmeans.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = deNovo.kmeans.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(deNovo.kmeans.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "deNovo_kmeans_BIC_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = deNovo.kmeans.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot sNMF K2 results
sNMF.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "sNMF", ] #subset sNMF observations
sNMF.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "sNMF" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset sNMF fitted predictions
sNMF.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "sNMF", ] #subset sNMF threshold rows
sNMF.k2.fitted.plot <- ggplot2::ggplot(sNMF.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = sNMF.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = sNMF.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(sNMF.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "sNMF_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = sNMF.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Plot STRUCTURE K2 results
STRUCTURE.k2.observed.plot.table <- fitted.binomial.input.table[fitted.binomial.input.table$method == "STRUCTURE", ] #subset STRUCTURE observations
STRUCTURE.k2.prediction.plot.table <- fitted.GLM.prediction.table[fitted.GLM.prediction.table$method == "STRUCTURE" & is.finite(fitted.GLM.prediction.table$fitted.proportion.k2), ] #subset STRUCTURE fitted predictions
STRUCTURE.k2.threshold.plot.table <- fitted.GLM.threshold.table[fitted.GLM.threshold.table$method == "STRUCTURE", ] #subset STRUCTURE threshold rows
STRUCTURE.k2.fitted.plot <- ggplot2::ggplot(STRUCTURE.k2.observed.plot.table, ggplot2::aes(x = tdiv, y = proportion.k2, color = mig.tag, group = mig.tag)) +
  ggplot2::geom_point(size = plot.point.size, alpha = plot.point.alpha) +
  ggplot2::geom_line(data = STRUCTURE.k2.prediction.plot.table, ggplot2::aes(x = tdiv, y = fitted.proportion.k2, color = mig.tag, group = mig.tag), linewidth = plot.fitted.line.width, inherit.aes = FALSE) +
  ggplot2::geom_vline(data = STRUCTURE.k2.threshold.plot.table, ggplot2::aes(xintercept = tdiv.at.threshold.k2, color = mig.tag), linewidth = plot.threshold.line.width, linetype = plot.threshold.line.type, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mig.colors, breaks = mig.level.table$mig.tag, labels = mig.level.table$mig.tag) +
  ggplot2::labs(x = "Divergence time (generations)", y = "Proportion choosing k = 2", color = "Migration rate") +
  ggplot2::theme_classic(base_size = plot.base.size) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = plot.axis.title.size, face = "bold"),
                 axis.text = ggplot2::element_text(size = plot.axis.text.size),
                 legend.title = ggplot2::element_text(size = plot.legend.title.size, face = "bold"),
                 legend.text = ggplot2::element_text(size = plot.legend.text.size))
print(STRUCTURE.k2.fitted.plot)
ggplot2::ggsave(file.path(results.directory, "STRUCTURE_proportion_k2_fitted_by_tdiv_and_mig.svg"),
                plot = STRUCTURE.k2.fitted.plot,
                device = plot.output.device,
                width = plot.width.cm,
                height = plot.height.cm,
                units = plot.output.units)


## Inspect fitted 50% threshold results
fitted.GLM.summary.table
fitted.GLM.threshold.table



## Create PCA plot for one VCF file

## Show available migration settings
available.migration.table <- unique(vcf.metadata.table[, c("mig.tag", "mig")]) #show all available migration settings including mig = 0
available.migration.table <- available.migration.table[order(available.migration.table$mig), ]
rownames(available.migration.table) <- NULL
print(available.migration.table)
if (!(migration.vcf.mig.for.PCA %in% unique(vcf.metadata.table$mig))) stop(paste("migration.vcf.mig.for.PCA must be one of:", paste(sort(unique(vcf.metadata.table$mig)), collapse = ", ")))


## Select file number N among VCF files for selected migration rate
migration.vcf.metadata.table <- vcf.metadata.table[vcf.metadata.table$mig == migration.vcf.mig.for.PCA, ] #keep selected VCF files only
migration.vcf.metadata.table <- migration.vcf.metadata.table[order(migration.vcf.metadata.table$tdiv), ] #order selected files by divergence time
rownames(migration.vcf.metadata.table) <- NULL
if (nrow(migration.vcf.metadata.table) == 0) stop(paste("No VCF files found for migration rate", migration.vcf.mig.for.PCA))
if (migration.vcf.file.number.for.PCA < 1 || migration.vcf.file.number.for.PCA > nrow(migration.vcf.metadata.table)) stop(paste("migration.vcf.file.number.for.PCA must be between 1 and", nrow(migration.vcf.metadata.table), "for migration rate", migration.vcf.mig.for.PCA))
current.pca.metadata <- migration.vcf.metadata.table[migration.vcf.file.number.for.PCA, ] #select VCF file N for selected migration rate
cat("\nCreating PCA for VCF file", migration.vcf.file.number.for.PCA, "of", nrow(migration.vcf.metadata.table), "for migration rate", migration.vcf.mig.for.PCA, ":", current.pca.metadata$file, "\n")
cat("Migration-index tag:", current.pca.metadata$mig.index, "\n")
cat("Divergence time:", current.pca.metadata$tdiv, "generations\n")
cat("Migration tag:", current.pca.metadata$mig.tag, "\n")
cat("Migration rate:", current.pca.metadata$mig, "\n")


## Read VCF and process SNP matrix
pca.genind.object.raw <- read.fastsimcoal2.genind.for.analysis(vcf.file.path = current.pca.metadata$full.path,
                                                               n.per.population = analysis.sample.n.per.population,
                                                               sampling.group = current.pca.metadata$simulation.group,
                                                               sampling.seed = analysis.sample.random.seed) #read selected VCF, assign populations, and optionally subset individuals per population


## Calculate Weir and Cockerham Fst before PCA processing
pca.fst.value <- calculate.overall.Fst(genind.object = pca.genind.object.raw) #calculate overall Weir and Cockerham Fst
cat("Weir and Cockerham Fst between pop1 and pop2:", round(pca.fst.value, 4), "\n")


## Process SNP matrix for PCA
pca.snp.matrix <- process.SNP.data.SOM(genind.input = pca.genind.object.raw, verbose = FALSE) #process SNP data for PCA
pca.snp.matrix <- as.matrix(pca.snp.matrix) #coerce to matrix
storage.mode(pca.snp.matrix) <- "numeric" #ensure numeric matrix
if (is.null(rownames(pca.snp.matrix))) rownames(pca.snp.matrix) <- adegenet::indNames(pca.genind.object.raw)


## Remove invalid loci
valid.pca.loci <- apply(pca.snp.matrix, 2, function(locus.values) {
  non.missing.values <- locus.values[!is.na(locus.values)]
  length(non.missing.values) > 1 && stats::var(non.missing.values) > 0
}) #identify loci with nonzero variance
pca.snp.matrix <- pca.snp.matrix[, valid.pca.loci, drop = FALSE] #keep valid loci only
if (ncol(pca.snp.matrix) < 2) stop("Fewer than two valid loci available for PCA")


## Impute remaining missing values for PCA
if (any(is.na(pca.snp.matrix))) {
  pca.snp.matrix <- apply(pca.snp.matrix, 2, function(locus.values) {
    locus.values[is.na(locus.values)] <- mean(locus.values, na.rm = TRUE)
    locus.values
  }) #replace remaining missing values with locus means
  pca.snp.matrix <- as.matrix(pca.snp.matrix)
}


## Run PCA and create plot table
pca.object <- stats::prcomp(pca.snp.matrix, center = TRUE, scale. = FALSE) #run PCA
pca.variance <- pca.object$sdev^2 / sum(pca.object$sdev^2) #calculate variance explained
pca.plot.table <- data.frame(individual = rownames(pca.object$x),
                             PC1 = pca.object$x[, 1],
                             PC2 = pca.object$x[, 2],
                             population = as.character(adegenet::pop(pca.genind.object.raw))[match(rownames(pca.object$x), adegenet::indNames(pca.genind.object.raw))],
                             file = current.pca.metadata$file,
                             mig.index = current.pca.metadata$mig.index,
                             tdiv = current.pca.metadata$tdiv,
                             mig.tag = current.pca.metadata$mig.tag,
                             mig = current.pca.metadata$mig,
                             fst.pop1.pop2 = pca.fst.value,
                             stringsAsFactors = FALSE) #create PCA plotting table


## Plot PCA with population groupings
pca.plot <- ggplot2::ggplot(pca.plot.table, ggplot2::aes(x = PC1, y = PC2, color = population, group = population)) +
  ggplot2::geom_point(size = 3) + 
  ggplot2::labs(title = paste0("PCA - ", current.pca.metadata$file, " | tdiv = ", round(current.pca.metadata$tdiv), " generations | mig = ", current.pca.metadata$mig, " | Fst = ", round(pca.fst.value, 4)),
                x = paste0("PC1 (", round(pca.variance[1] * 100, 2), "%)"),
                y = paste0("PC2 (", round(pca.variance[2] * 100, 2), "%)"),
                color = "Population") +
  ggplot2::theme_classic()
pca.plot


## Print selected PCA summary
cat("Divergence:", round(current.pca.metadata$tdiv / 1000000, 2), "million generations\n")
cat("Migration rate:", current.pca.metadata$mig, "\n")
cat("Weir and Cockerham Fst:", round(pca.fst.value, 4), "\n")
