# Fully integrative species delimitation with the *delimSOM* 2.0 *R* package

`delimSOM` is an *R* package for fully integrative species delimitation using single- and multilayer self-organizing maps (SOMs). 

SOMs are unsupervised machine-learning models that organize complex, high-dimensional observations onto a two-dimensional grid while preserving local relationships among similar observations (Kohonen 1998, 2014). SOMs were first applied to individual-based species delimitation using genetic data by Pyron et al. (2023) and were subsequently extended to multilayer analyses capable of integrating different sources of biological information (Pyron 2023). In `delimSOM`, each data type is retained as a separate layer, but all layers jointly contribute to a shared SOM through layer-specific distance functions and weighting. During training, individuals with similar combinations of measurements are represented by neighboring map units, while more dissimilar individuals become separated across the map. The resulting map units are represented by codebook vectors that summarize the multivariate information of the individuals assigned to them. These codebook vectors are then clustered to infer candidate lineages, and replicate SOM analyses are used to evaluate support for alternative numbers of clusters (`K`) and the stability of individual assignments. The framework does not require predefined species assignments and explicitly allows `K = 1`. Any data that can be represented quantitatively for a common set of individuals can be incorporated, including genomic, morphological, behavioral, ecological, spatial, and other continuous, binary, categorical, or count data.

The motivation behind `delimSOM` is that speciation is a multidimensional evolutionary process in which genetic differentiation, phenotypic evolution, ecological specialization, reproductive isolation, and geographic isolation can arise through different mechanisms and accumulate at different rates (de Queiroz 2007). Robust species hypotheses are therefore more likely to emerge when multiple complementary sources of evidence are considered together rather than relying on a single property of lineage divergence (Dayrat 2005; Padial et al. 2010; Schlick-Steiner et al. 2010). Although integrative species delimitation has been advocated for decades, most existing workflows still analyze different sources of evidence separately and compare their resulting species hypotheses afterward, or use one data type to define candidate lineages that are subsequently evaluated with additional evidence. `delimSOM` addresses this limitation by allowing heterogeneous biological data to jointly inform species delimitation within a single unsupervised analysis. Variables are normalized before training, differences in the distance scales of layers are automatically balanced, and users can additionally specify data-type-specific distance functions or layer weights when appropriate. This makes it possible to integrate information that differs substantially in scale, dimensionality, missingness, and biological interpretation while retaining each source of evidence as a distinct component of the analysis.

## Main advantages of the approach

- Does not require predefined species assignments.
- Robust to missing data because it retains partially incomplete individuals and variables through missing-data-aware SOM training.
- Automatically normalizes variables and balances contributions among SOM layers through layer-specific distance normalization.
- Supports user-defined layer weights .
- Trains replicate SOMs to quantify stability and support for inferred lineage structure.
- Explicitly allows `K = 1`, so a dataset is not forced to contain multiple candidate species.
- Includes preprocessing functions for SNP/sequence data, categorical variables, and low-variation/high-correlation filtering.
- Provides six alternative clustering and `K`-selection approaches.
- Includes diagnostics and visualizations for SOM learning, map quality, `K` support, replicate-consensus assignments, variable importance, and layer importance.


## Development status

We have now released version 2.0 of this method!

The framework is described in a preprint (xxx) and the manuscript is currently in review.

Current R package version: `2.0.0.9000`

For bug reports, feedback, or questions, please contact Daniel Schönberger: daniel.schoenberger@uky.edu.

# Tutorial

## Installation

Install and load the *delimSOM* *R* package:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("rpyron/delim-SOM", ref = "dev2.0")

library(delimSOM)
packageVersion("delimSOM")
```


## Input data

The basic SOM input is either:

1. one numeric matrix/data frame for a single-layer analysis, or
2. a named list of numeric matrices/data frames for a multilayer analysis.

Rows represent individuals and columns represent variables.

For multilayer analyses, use individual identifiers as row names in every layer. `train.SOM()` identifies the shared individuals across layers and retains the common set in matching order.

Missing values should be represented as `NA`.

A simple multilayer input therefore looks like:

```r
SOM_data <- list(
  Genomic = genomic_matrix,
  Morphology = morphology_matrix,
  Environmental = environmental_matrix,
  Spatial = spatial_matrix
)
```

### Genetic and sequence data

`process.SNP.data.SOM()` converts several common genetic inputs into matrices suitable for SOM analysis. Supported inputs include VCF, `genind`, `genlight`, numeric SNP-dosage matrices, PLINK `.raw`, NEXUS, FASTA, and PHYLIP.

For example:

```r
SNP_data <- process.SNP.data.SOM(
  vcf.path = "data.vcf",
  missing.loci.cutoff.lenient = 0.7,
  missing.loci.cutoff.final = 0.5,
  missing.individuals.cutoff = 0.5
)
```

### Categorical data

Categorical variables can be converted to binary indicator variables with `make.cols.binary.SOM()`:

```r
binary_data <- make.cols.binary.SOM(
  dataframe = categorical_data,
  make.binary.cols = c("Host", "Habitat")
)
```

### Continuous numeric data

Continuous environmental, morphological, physiological, spatial, or other numeric matrices can be filtered with `remove.lowCV.multicollinearity.SOM()`:

```r
continuous_data <- remove.lowCV.multicollinearity.SOM(
  continuous_data,
  CV.threshold = 0.05,
  cor.threshold = 0.9
)
```

`train.SOM()` additionally performs its own row matching, missing-data filtering, zero-variance filtering, and variable-wise min-max normalization before SOM training.

## Distance functions

When `layer.distance.functions = NULL`, `train.SOM()` automatically infers a distance type from each layer:

- binary `0/1` data: `"tanimoto"`
- dosage/count-like `0/1/2` data: `"manhattan"`
- other numeric continuous data: `"sumofsquares"`

Users should verify that the automatically inferred distance is biologically appropriate. Distances can be supplied manually when needed:

```r
SOM_tr <- train.SOM(
  input_data = SOM_data,
  layer.distance.functions = c(
    "manhattan",
    "sumofsquares",
    "sumofsquares",
    "sumofsquares"
  )
)
```

## Recommended workflow

The core analysis consists of training replicate SOMs and then clustering their codebook vectors.

```r
SOM_tr <- train.SOM(
  input_data = SOM_data
)

SOM_results <- clustering.SOM(
  SOM.output = SOM_tr,
  clustering.method = "kmeans+BICelbow"
)
```

We recommend starting with `"kmeans+BICelbow"` because the simulation analyses associated with delim-SOM 2.0 found that the k-means/BIC approaches provided a strong balance of lineage recovery, assignment accuracy, and computational efficiency.

The automatically selected `K` should not be interpreted in isolation. Inspect the full `K`-support profile, replicate-consensus assignments, SOM topology, and—when biologically appropriate—hierarchical analyses within recovered candidate lineages.

# Empirical tutorial: Polygonia anglewing butterflies

The following example reproduces the main empirical focus study using western Canadian *Polygonia* anglewing butterflies from Dupuis et al. (2018):

> Dupuis, J. R. et al. (2018). [Genomics confirms surprising ecological divergence and isolation in an enigmatic butterfly species complex](https://doi.org/10.1093/zoolinnean/zlx081).

The original study recognized four species:

- *Polygonia faunus*
- *Polygonia gracilis*
- *Polygonia progne*
- *Polygonia satyrus*

For the delim-SOM 2.0 reanalysis, 200 shared individuals were analyzed across six complementary layers:

1. genome-wide GBS SNPs
2. mitochondrial COI
3. continuous morphology (RGB wing-color measurements)
4. categorical morphology (visually scored wing characters and morphotype)
5. environmental variables
6. spatial variables

The example data are stored in:

```text
Empirical_examples/Dupuis_et_al_2018/
```

Because `Empirical_examples/` is retained in the GitHub repository but excluded from the installed R package, clone or download the GitHub repository if you want to run this empirical example locally.

For example:

```text
git clone --branch dev2.0 https://github.com/rpyron/delim-SOM.git
```

Then run the tutorial from the repository root.

The environmental table included in the repository has already been extracted from geographic coordinates. The full manuscript analysis used `NicheDiv` for environmental preprocessing, so install it if you want to reproduce that step exactly:

```r
if (!requireNamespace("NicheDiv", quietly = TRUE)) {
  remotes::install_github("Daniel-1232/NicheDiv")
}
```

## 1. Set paths

```r
#### Set environment ###########################################################

example_dir <- file.path(
  "Empirical_examples",
  "Dupuis_et_al_2018"
)
```

## 2. Process genome-wide SNP data

The SNP layer contains GBS-derived biallelic SNPs.

```r
#### Process SNP data ##########################################################

Polygonia_SNP <- process.SNP.data.SOM(
  vcf.path = file.path(
    example_dir,
    "Polygonia_961SNPs.vcf"
  ),
  missing.loci.cutoff.lenient = 0.7,
  missing.loci.cutoff.final = 0.5,
  missing.individuals.cutoff = 0.5
)

rownames(Polygonia_SNP) <- sub(
  ".*?(\\d+)$",
  "\\1",
  rownames(Polygonia_SNP)
)

dim(Polygonia_SNP)
```

The manuscript analysis retained 961 SNP variables.

## 3. Process mitochondrial COI data

`process.SNP.data.SOM()` can also process aligned nucleotide sequences directly.

```r
#### Process COI data ##########################################################

Polygonia_COI <- process.SNP.data.SOM(
  nexus.path = file.path(
    example_dir,
    "Polygonia_COI.nex"
  ),
  missing.loci.cutoff.lenient = 0.7,
  missing.loci.cutoff.final = 0.5,
  missing.individuals.cutoff = 0.5
)

Polygonia_COI_numeric_rownames <- sub(
  ".*?(\\d+)$",
  "\\1",
  rownames(Polygonia_COI)
)

Polygonia_COI <- Polygonia_COI[
  !duplicated(Polygonia_COI_numeric_rownames),
  ,
  drop = FALSE
]

rownames(Polygonia_COI) <- Polygonia_COI_numeric_rownames[
  !duplicated(Polygonia_COI_numeric_rownames)
]

dim(Polygonia_COI)
```

The manuscript analysis retained 213 COI variables after processing.

## 4. Prepare continuous wing-color morphology

The continuous morphology layer contains RGB measurements from dorsal and ventral wing regions.

```r
#### Prepare continuous morphology ############################################

Polygonia_RGB <- read.delim(
  file.path(
    example_dir,
    "Polygonia_RGB_characters.txt"
  ),
  stringsAsFactors = FALSE
)

rownames(Polygonia_RGB) <- Polygonia_RGB$Species

Polygonia_RGB <- Polygonia_RGB[
  ,
  !names(Polygonia_RGB) %in% c(
    "Name",
    "Species"
  ),
  drop = FALSE
]
```

## 5. Prepare categorical wing morphology

The second morphology layer contains visually scored wing characters. Wing character 8 is treated as nominal rather than ordinal and is therefore expanded into binary indicator columns.

```r
#### Prepare categorical morphology ###########################################

Polygonia_wing_scores <- read.delim(
  file.path(
    example_dir,
    "Polygonia_visually_scored.txt"
  ),
  stringsAsFactors = FALSE
)

rownames(Polygonia_wing_scores) <- Polygonia_wing_scores$Name

Polygonia_wing_scores <- Polygonia_wing_scores[
  ,
  !names(Polygonia_wing_scores) %in% c(
    "Name",
    "Species"
  ),
  drop = FALSE
]

colnames(Polygonia_wing_scores)[
  match(
    paste0("Ch", 1:10),
    colnames(Polygonia_wing_scores)
  )
] <- paste0(
  "Wing_character_",
  1:10
)

Polygonia_wing_scores$Wing_character_8 <- factor(
  Polygonia_wing_scores$Wing_character_8,
  levels = 1:4
)

Wing_character_8_states <- stats::model.matrix(
  ~ Wing_character_8 - 1,
  data = Polygonia_wing_scores
)

colnames(Wing_character_8_states) <- paste0(
  "Wing_character_8_state_",
  1:4
)

Polygonia_wing_scores <- cbind(
  Polygonia_wing_scores,
  Wing_character_8_states
)

Polygonia_wing_scores$Wing_character_8 <- NULL
```

## 6. Import metadata and morphotype

```r
#### Import metadata ###########################################################

Polygonia_metadata <- read.csv(
  file.path(
    example_dir,
    "Polygonia_metadata.csv"
  ),
  header = TRUE,
  sep = ";"
)

rownames(Polygonia_metadata) <- Polygonia_metadata$ID

Polygonia_morphotype <- Polygonia_metadata[
  ,
  "Morphotype",
  drop = FALSE
]
```

## 7. Combine and separate morphology layers

The manuscript analysis separates continuous RGB morphology from categorical/ordinal wing characters because these variables have different data structures.

```r
#### Combine morphology datasets ##############################################

rownames(Polygonia_RGB) <- as.character(
  rownames(Polygonia_RGB)
)

rownames(Polygonia_wing_scores) <- as.character(
  rownames(Polygonia_wing_scores)
)

rownames(Polygonia_morphotype) <- as.character(
  rownames(Polygonia_morphotype)
)

Polygonia_RGB_wing_scores <- merge(
  Polygonia_RGB,
  Polygonia_wing_scores,
  by = "row.names",
  all = FALSE
)

rownames(Polygonia_RGB_wing_scores) <- Polygonia_RGB_wing_scores$Row.names
Polygonia_RGB_wing_scores$Row.names <- NULL

Polygonia_morphology_all <- merge(
  Polygonia_RGB_wing_scores,
  Polygonia_morphotype,
  by = "row.names",
  all = FALSE
)

rownames(Polygonia_morphology_all) <- Polygonia_morphology_all$Row.names
Polygonia_morphology_all$Row.names <- NULL

Polygonia_morphology_all <- make.cols.binary.SOM(
  Polygonia_morphology_all,
  make.binary.cols = "Morphotype",
  append.to.original = TRUE
)

Polygonia_morphology_all$Morphotype <- NULL

non.continuous.cols <- grepl(
  "^Wing_character_|^Morphotype_",
  colnames(Polygonia_morphology_all)
)

Polygonia_morphology_categorical <- Polygonia_morphology_all[
  ,
  non.continuous.cols,
  drop = FALSE
]

Polygonia_morphology <- Polygonia_morphology_all[
  ,
  !non.continuous.cols,
  drop = FALSE
]

Polygonia_morphology <- remove.lowCV.multicollinearity.SOM(
  Polygonia_morphology,
  CV.threshold = 0.05,
  cor.threshold = 0.9
)

dim(Polygonia_morphology)
dim(Polygonia_morphology_categorical)
```

The final manuscript analysis contained 15 continuous and 15 categorical morphology variables.

## 8. Prepare environmental and spatial layers

The repository includes the environmental variables previously extracted from the specimen coordinates. To keep this tutorial reproducible without another live GIS download, we use the coordinate/elevation columns in that prepared file as the spatial layer.

```r
#### Prepare environmental and spatial data ###################################

Polygonia_environmental <- read.csv(
  file.path(
    example_dir,
    "Polygonia_environmental.csv"
  ),
  row.names = 1,
  header = TRUE
)

Polygonia_spatial <- Polygonia_environmental[
  ,
  c(
    "Latitude",
    "Longitude",
    "Elevation"
  ),
  drop = FALSE
]

Polygonia_environmental <- Polygonia_environmental[
  ,
  !names(Polygonia_environmental) %in% c(
    "Latitude",
    "Longitude",
    "Elevation"
  ),
  drop = FALSE
]

Polygonia_environmental[] <- lapply(
  Polygonia_environmental,
  as.numeric
)

Polygonia_environmental <- (
  NicheDiv::transform.skewed.variables(
    Polygonia_environmental
  )
)$transformed

Polygonia_environmental <- remove.lowCV.multicollinearity.SOM(
  Polygonia_environmental,
  CV.threshold = 0.05,
  cor.threshold = 0.9
)

dim(Polygonia_environmental)
dim(Polygonia_spatial)
```

The manuscript analysis retained 125 environmental variables and three spatial variables.

## 9. Match individuals across all six layers

All SOM layers must refer to the same final set of individuals.

```r
#### Match individuals across layers ##########################################

Polygonia_metadata <- Polygonia_metadata[
  ,
  c(
    "Species",
    "ID"
  ),
  drop = FALSE
]

Polygonia_common_IDs <- Reduce(
  intersect,
  list(
    rownames(Polygonia_morphology),
    rownames(Polygonia_morphology_categorical),
    rownames(Polygonia_SNP),
    rownames(Polygonia_COI),
    rownames(Polygonia_spatial),
    rownames(Polygonia_environmental),
    rownames(Polygonia_metadata)
  )
)

Polygonia_morphology <- Polygonia_morphology[
  Polygonia_common_IDs,
  ,
  drop = FALSE
]

Polygonia_morphology_categorical <- Polygonia_morphology_categorical[
  Polygonia_common_IDs,
  ,
  drop = FALSE
]

Polygonia_SNP <- Polygonia_SNP[
  Polygonia_common_IDs,
  ,
  drop = FALSE
]

Polygonia_COI <- Polygonia_COI[
  Polygonia_common_IDs,
  ,
  drop = FALSE
]

Polygonia_spatial <- Polygonia_spatial[
  Polygonia_common_IDs,
  ,
  drop = FALSE
]

Polygonia_environmental <- Polygonia_environmental[
  Polygonia_common_IDs,
  ,
  drop = FALSE
]

Polygonia_metadata <- Polygonia_metadata[
  Polygonia_common_IDs,
  ,
  drop = FALSE
]

length(Polygonia_common_IDs)
```

This should return 200 individuals for the manuscript dataset.

For easier interpretation of downstream plots, append the recognized species name to each individual identifier:

```r
Polygonia_species_vec <- Polygonia_metadata[
  Polygonia_common_IDs,
  "Species"
]

Polygonia_new_rownames <- paste(
  Polygonia_common_IDs,
  Polygonia_species_vec,
  sep = "_"
)

rownames(Polygonia_morphology) <- Polygonia_new_rownames
rownames(Polygonia_morphology_categorical) <- Polygonia_new_rownames
rownames(Polygonia_SNP) <- Polygonia_new_rownames
rownames(Polygonia_COI) <- Polygonia_new_rownames
rownames(Polygonia_spatial) <- Polygonia_new_rownames
rownames(Polygonia_environmental) <- Polygonia_new_rownames
rownames(Polygonia_metadata) <- Polygonia_new_rownames
```

## 10. Combine the six SOM layers

```r
#### Combine SOM layers ########################################################

Polygonia_all_data <- list(
  Morphology = Polygonia_morphology,
  Morphology_2 = Polygonia_morphology_categorical,
  SNP = Polygonia_SNP,
  COI = Polygonia_COI,
  Environmental = Polygonia_environmental,
  Spatial = Polygonia_spatial
)

sapply(
  Polygonia_all_data,
  dim
)
```

Each layer remains separate during training but jointly contributes to the shared SOM.

## 11. Train replicate SOMs

The default workflow trains 110 replicate SOMs with 100 training steps per replicate.

```r
#### Train SOM #################################################################

Polygonia_SOM_tr <- train.SOM(
  input_data = Polygonia_all_data,
  max.NA.row = 0.5,
  max.NA.col = 0.5,
  save.SOM.results = FALSE
)
```

During training, `train.SOM()` reports the inferred layer types and assigned distance functions. Check these messages to confirm that the chosen distances are appropriate for your data.

For the manuscript Polygonia analysis, SOM training took approximately four minutes on the benchmark system. Runtime depends on hardware and parallel settings.

## 12. Cluster SOM codebook vectors

We use the recommended k-means + BIC-elbow approach:

```r
#### Cluster SOM ###############################################################

Polygonia_SOM <- clustering.SOM(
  SOM.output = Polygonia_SOM_tr,
  clustering.method = "kmeans+BICelbow",
  save.SOM.results = FALSE
)

Polygonia_SOM$optim_k_summary
```

In the manuscript analysis, the replicate summary was:

```text
K = 3: 90%
K = 4: 10%
```

The dominant `K = 3` solution recovered *P. faunus* and *P. satyrus* as separate candidate lineages and combined the closely related *P. gracilis* and *P. progne* into a third candidate lineage. The latter two species were nevertheless separated in 10% of SOM replicates and were resolved as separate lineages in a subsequent hierarchical reanalysis.

## 13. Evaluate SOM training and `K` support

### Learning trajectories

```r
plot.learning.SOM(
  Polygonia_SOM
)
```

A rapid decline followed by a stable plateau indicates that SOM learning has converged. Erratic trajectories or continued large changes late in training can indicate that more training steps are needed.

### Layer distance scales

```r
plot.layer.distance.scale.SOM(
  Polygonia_SOM
)
```

This is a diagnostic of raw whole-layer distance scale before internal normalization. It should not be interpreted as biological layer importance.

### Support for alternative `K`

```r
plot.K.SOM(
  Polygonia_SOM
)
```

For BIC-based methods, this plot shows the support profile across candidate values of `K`, successive changes in BIC, and the frequency with which each `K` was selected across retained SOM replicates.

Always inspect the full support profile rather than relying only on the modal `K`.

## 14. Visualize SOM topology and candidate lineages

```r
plot.model.SOM(
  Polygonia_SOM,
  replicate.mode = "representative"
)
```

The neighbor-distance panel shows distances among adjacent SOM units. Darker/high-distance regions indicate stronger discontinuities or "ridges" in map space. The cluster panel shows the inferred candidate-lineage boundaries.

A specific value of `K` can also be inspected:

```r
plot.model.SOM(
  Polygonia_SOM,
  replicate.mode = "representative",
  set.k = 3
)

plot.model.SOM(
  Polygonia_SOM,
  replicate.mode = "representative",
  set.k = 4
)
```

## 15. Plot replicate-consensus assignments

`plot.structure.SOM()` provides a STRUCTURE-like visualization of replicate-consensus assignment coefficients:

```r
plot.structure.SOM(
  Polygonia_SOM,
  bottom.margin = 8
)
```

Each bar represents an individual. Assignment distributed across more than one cluster indicates instability in the recovered grouping across SOM replicates and may reflect weak differentiation, admixture, recent divergence, or conflicting signal among data layers.

## 16. Plot geographic assignments

```r
plot.map.SOM(
  SOM.output = Polygonia_SOM,
  Coordinates = Polygonia_spatial[
    ,
    c(
      "Latitude",
      "Longitude"
    )
  ],
  lat.buffer.range = 4,
  lon.buffer.range = 5,
  north.arrow.position = c(
    0.04,
    0.87
  ),
  north.arrow.length = 0.7,
  north.arrow.N.position = 0.3,
  north.arrow.N.size = 1,
  scale.position = c(
    0.75,
    0.05
  )
)
```

The map displays replicate-consensus candidate-lineage assignment coefficients at the geographic coordinates of each individual.

## 17. Evaluate variable importance

`clustering.SOM()` calculates two complementary variable-importance summaries.

### Cluster-separation importance

```r
plot.variable.importance.SOM(
  Polygonia_SOM,
  mode = "Cluster.separation",
  bottom.margin = 2,
  left.margin = 5.8
)
```

This uses an ANOVA-like eta-squared effect size to quantify how strongly each variable is associated with separation among inferred clusters.

### Map-variance importance

```r
plot.variable.importance.SOM(
  Polygonia_SOM,
  mode = "Map.variance",
  left.margin = 5
)
```

Map variance measures how strongly a variable changes across the trained SOM map, irrespective of whether that variation aligns exactly with the final cluster boundaries.

The stored values can also be inspected directly:

```r
head(
  sort(
    Polygonia_SOM$median_etasquared_variable_importance$Morphology,
    decreasing = TRUE
  ),
  15
)
```

## 18. Evaluate layer importance

### Summarize variable importance by layer

```r
plot.layer.importance.varimp.SOM(
  Polygonia_SOM,
  bottom.margin = 4
)
```

### Leave-one-layer-out analysis

For a stronger test of layer importance, rerun the analysis while removing one layer at a time:

```r
plot.layer.importance.leaveoneout.SOM(
  Polygonia_SOM,
  bottom.margin = 7
)
```

This compares replicate-matched full and reduced analyses and evaluates changes in:

- the inferred number of clusters
- cluster composition
- individual assignment confidence

This analysis is substantially more computationally expensive because SOM training and clustering are repeated for each omitted layer.

In the Polygonia manuscript analysis, the SNP layer had the strongest effect on the inferred clustering, followed by categorical wing-pattern morphology and mitochondrial COI. Environmental and spatial layers had weaker effects.

## 19. Optional: compare clustering methods

`clustering.SOM()` currently implements six clustering and `K`-selection methods:

```r
clustering_methods <- c(
  "kmeans+BICelbow",
  "kmeans+BICthreshold",
  "GMM+BICthreshold",
  "hierarchical+DB",
  "HDBSCAN",
  "OPTICS+Silhouette"
)
```

For example:

```r
Polygonia_SOM_HDBSCAN <- clustering.SOM(
  SOM.output = Polygonia_SOM_tr,
  clustering.method = "HDBSCAN",
  save.SOM.results = FALSE
)

Polygonia_SOM_HDBSCAN$optim_k_summary
```

No single clustering method is expected to perform best for every possible data structure. The k-means/BIC methods are recommended as the primary starting point, while alternative methods can provide useful complementary analyses.

## 20. Optional: hierarchical reanalysis

A conservative primary analysis may combine weakly differentiated lineages. Hierarchical reanalysis can then test for additional structure within each recovered candidate lineage.

For the Polygonia example:

```r
#### Assign individuals to dominant candidate lineage #########################

Polygonia_clusters <- apply(
  Polygonia_SOM$ancestry_matrix,
  1,
  which.max
)

Polygonia_clusters <- paste0(
  "cluster",
  Polygonia_clusters
)

table(
  Polygonia_clusters
)

Polygonia_cluster_samples <- split(
  rownames(
    Polygonia_SOM$ancestry_matrix
  ),
  Polygonia_clusters
)
```

The main Polygonia analysis recovered three clusters containing approximately 79, 75, and 46 individuals. The middle cluster contained *P. gracilis* + *P. progne*.

Subset that candidate lineage and rerun the workflow:

```r
Polygonia_cluster2_data <- lapply(
  Polygonia_SOM$input_data,
  function(x) {
    x[
      Polygonia_cluster_samples$cluster2,
      ,
      drop = FALSE
    ]
  }
)

Polygonia_SOM_tr_cluster2 <- train.SOM(
  input_data = Polygonia_cluster2_data,
  max.NA.row = 0.5,
  max.NA.col = 0.5,
  save.SOM.results = FALSE
)

Polygonia_SOM_cluster2 <- clustering.SOM(
  SOM.output = Polygonia_SOM_tr_cluster2,
  clustering.method = "kmeans+BICelbow",
  max.k = 5,
  save.SOM.results = FALSE
)

Polygonia_SOM_cluster2$optim_k_summary
```

The manuscript analysis recovered `K = 2` with 100% support in this hierarchical analysis, separating *P. gracilis* and *P. progne*. No further subdivision was supported within the primary *P. faunus* or *P. satyrus* lineages.

# Applying delimSOM to your own data

For most projects, the workflow can be reduced to:

```r
#### Prepare layers ############################################################

SOM_data <- list(
  Layer_1 = matrix_1,
  Layer_2 = matrix_2,
  Layer_3 = matrix_3
)



#### Train SOM #################################################################

SOM_tr <- train.SOM(
  input_data = SOM_data
)



#### Cluster SOM ###############################################################

SOM_results <- clustering.SOM(SOM.output = SOM_tr,
                              clustering.method = "kmeans+BICelbow")



#### Evaluate results ##########################################################

SOM_results$optim_k_summary

plot.learning.SOM(SOM_results)

plot.K.SOM(SOM_results)

plot.model.SOM(
  SOM_results,
  replicate.mode = "representative"
)

plot.structure.SOM(SOM_results)
```

Important considerations:

- Use the same biological individuals across data layers whenever possible.
- Ensure that row names uniquely identify individuals.
- Encode missing values as `NA`.
- Inspect training convergence and model-quality diagnostics.
- Inspect support across alternative values of `K`, not only the modal solution.
- Treat inferred clusters as candidate-lineage hypotheses rather than definitive taxonomic conclusions.
- Consider hierarchical reanalysis when weaker structure may occur within a strongly supported candidate lineage.
- Species-rich systems and datasets with small or strongly uneven within-lineage sample sizes require more cautious interpretation.

# Main functions

## Data preparation

```r
process.SNP.data.SOM()
make.cols.binary.SOM()
remove.lowCV.multicollinearity.SOM()
```

## SOM training and clustering

```r
train.SOM()
clustering.SOM()
```

## Visualization and model evaluation

```r
plot.learning.SOM()
plot.layer.distance.scale.SOM()
plot.K.SOM()
plot.model.SOM()
plot.structure.SOM()
plot.map.SOM()
plot.variable.importance.SOM()
plot.layer.importance.varimp.SOM()
plot.layer.importance.leaveoneout.SOM()
```

# Citation

Please cite the *delimSOM* framework as follows

Schönberger, D., Pyron, R. A., & Dupuis, J. R. *delim-SOM 2.0*: Fully integrative species delimitation with machine learning and flexible diverse data types of biological and other data. *bioRxiv*


# References

Dupuis, J. R. et al. (2018). Genomics confirms surprising ecological divergence and isolation in an enigmatic butterfly species complex. *Zoological Journal of the Linnean Society*. https://doi.org/10.1093/zoolinnean/zlx081

Kohonen, T. (1998). The self-organizing map. *Neurocomputing*.

Kohonen, T. (2014). *MATLAB Implementations and Applications of the Self-Organizing Map*.

Pyron, R. A. (2023). Unsupervised machine learning for species delimitation, integrative taxonomy, and biodiversity conservation. *Molecular Phylogenetics and Evolution*, 189, 107939. https://doi.org/10.1016/j.ympev.2023.107939

Pyron, R. A., O’Connell, K. A., Duncan, S. C., Burbrink, F. T., & Beamer, D. A. (2023). Speciation hypotheses from phylogeographic delimitation yield an integrative taxonomy for seal salamanders (*Desmognathus monticola*). *Systematic Biology*, 72(1), 179–197. https://doi.org/10.1093/sysbio/syac065