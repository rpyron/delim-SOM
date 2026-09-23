# delimSOM

`delimSOM` implements Self-Organizing Maps (SOMs) for integrative species delimitation using single-layer and multi-layer biological datasets.

The package can jointly analyze heterogeneous data such as genomic, morphological, environmental, spatial, binary, categorical, count, and other multivariate data while retaining individual datasets as separate SOM layers.

## Installation

The development version can be installed from GitHub with:

```r
install.packages("remotes")
remotes::install_github("rpyron/delim-SOM", ref = "dev2.0")
```

Then load the package with:

```r
library(delimSOM)
```

## Main workflow

The two primary analysis functions are:

```r
train.SOM()
clustering.SOM()
```

`train.SOM()` trains single-layer or multi-layer SOM models. `clustering.SOM()` subsequently clusters the SOM codebook vectors to infer delimitation structure.

## Data preparation

Functions for preparing input data include:

```r
process.SNP.data.SOM()
make.cols.binary.SOM()
remove.lowCV.multicollinearity.SOM()
```

## Visualization and model evaluation

Available plotting functions include:

```r
plot.structure.SOM()
plot.learning.SOM()
plot.layer.distance.scale.SOM()
plot.K.SOM()
plot.model.SOM()
plot.map.SOM()
plot.variable.importance.SOM()
plot.layer.importance.varimp.SOM()
plot.layer.importance.leaveoneout.SOM()
```

Detailed arguments and examples are available through the R help pages, for example:

```r
?train.SOM
?clustering.SOM
?plot.structure.SOM
```

## Development status

Version 2.0 is currently under development on the `dev2.0` branch.
