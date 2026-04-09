# Fully integrative species delimitation with unsupervised machine learning using SOMs 

## Introduction
Our <code>delim-SOM</code> R package implements a Self-Organizing Map (SOM) framework for fully integrative species delimitation using genomic, phenotypic, environmental, spatial, and other multivariate data. The package implementation and approach are desribed in detail in Schönberger et al. (2026), building up on the single- and multilayer SOM workflow introduced by Pyron (2023) and Pyron et al. (2023), and relying on the <code>kohonen</code> R package (Wehrens and Buydens 2007, 2018)
The approach consists of two main steps. First, the original high-dimensional data are regularized into a two-dimensional topology-preserving grid of SOM codebook vectors (i.e., prototypes). Second, those codebook vectors are clustered to infer the delimitation structure. This two-stage strategy follows the classical use of SOMs as competitive, unsupervised neural networks that represent high-dimensional observations by prototype vectors arranged on a low-dimensional grid, thereby emphasizing local neighborhood preservation rather than global distance preservation (Kohonen, 1990, 1998, 2001; Wehrens & Buydens, 2007). Because the method is unsupervised, species labels or other a priori group assignments do not need to be provided during map training, but delimitation structure is inferred directly from multivariate similarity among observations. During training, each observation is matched to the map unit whose prototype is most similar to it (best-matching unit), and nearby units are updated as well so that similar observations become organized into neighboring regions of the map. In this framework, neighboring map units represent relatively similar observations, whereas more distant regions of the map summarize broader differences in multivariate structure (Kohonen, 1998; Wehrens & Buydens, 2007). Because SOM codebook vectors are local prototype averages rather than individual observations, the SOM acts as a regularizing representation that can reduce sensitivity to noise and individual-level irregularities while preserving major topological relationships in the data (Kohonen, 2001, 2014; Vesanto & Alhoniemi, 2000). This makes SOMs especially useful for complex biological datasets in which structure may be nonlinear, weak, partially conflicting across variables, or unevenly distributed across different data types (Lek et al., 1996; Lek & Guégan, 2000; Zamani et al., 2013). SOMs are also well-suited to integrative delimitation because multilayer SOMs can accommodate heterogeneous data types through layer-specific distance functions and weighting, allowing complementary signals from genetics, phenotype, environment, and space to contribute jointly to the learned map without requiring a single explicit generative model for all data types (Pyron, 2023; Wehrens & Buydens, 2007; Wehrens & Kruisselbrink, 2018). More generally, SOMs are exploratory, structure-discovering methods that can reveal implicit clustering tendencies and hierarchical relationships in multivariate data, while subsequent clustering of codebook vectors provides a second abstraction level for delimitation inference rather than forcing discrete groups during map training (Akman et al., 2019; Kohonen, 2013; Vesanto & Alhoniemi, 2000).

The contribution of each layer to the final model output is recorded, along with the clustering assignment of each individual over multiple learning replicates. The results therefore resemble a 'STRUCTURE'-type analysis including admixture estimates, but represent a unified delimitation model that incorporates various dimensions of ecological and evolutionary divergence for fully integative taxonomy. 
Overall, the method is extremely flexible and can take almost any data type or format (e.g., continuous, binary, categorical, count, SNP), as long as it is provided to function as a matrix/matrices or dataframe(s).

## How does SOM work - lay summary
A Self-Organizing Map (SOM) is a way of turning a complicated table of measurements into a two-dimensional map that is easier to understand. You can think of the map as a grid of little cells, where each cell stores a representative “average pattern” (called codebook vector or prototype). Each specimen is then compared to all these cells, and it is placed closest to the cell whose pattern is most similar to its own overall combination of measurements (called best-matching unit). During training, not only the winning cell, but also nearby cells on the map, are adjusted to become a little more similar to the specimen. After many repeated updates across all specimens, the grid becomes organized so that similar specimens tend to occupy the same or neighboring parts of the map, whereas more different specimens tend to fall farther apart. In this sense, the SOM does not try to preserve every distance exactly. Instead, it mainly preserves local similarity relationships (topology preservation), so nearby map positions are usually the most informative. Another useful feature is that the map summarizes many individual observations into a smaller set of representative local averages, which can reduce sensitivity to noise and make broad structure easier to see. It also tends to devote more map space to dense, information-rich parts of the dataset and less space to sparse regions. This makes SOMs especially helpful for complex biological datasets in which relationships may be nonlinear, partly conflicting across variables, or distributed across different data types. For example, if several individuals are genetically similar, occur in similar environments, and share similar phenotypic traits, the SOM will tend to place them in the same region of the map; individuals with different combinations of those features will tend to be placed in other regions. After that map has been learned, the representative map cells themselves can be grouped into larger clusters, which provides a second and more stable level of inference for delimitation rather than forcing final groups directly during the initial mapping step.

## How to load the package
Load the R package by running:
```
source("https://raw.githubusercontent.com/rpyron/delim-SOM/refs/heads/dev2.0/R/2026_04_07_delim-SOM_2.0_functions.R")
```

## Package functions
There are two important functions for training SOMs and clustering the results (i.e., the SOM codebook vectors)
```
train.SOM()
clustering.SOM()
```

We provide two functions to facilitate preparing the data for input 

```
process.SNP.data.SOM()
make.cols.binary.SOM()

```

Multiple functions are provided to evaluate and visualize the results

A STRUCTURE-type barplot:
```
plot.structure.SOM()
```

And the variable importance across layers:

```
plot.layer.importance.varimp.SOM()
```

## Empirical example: seal salamander (_Desmognathus monticola_: Plethodontidae) 

![Pyron_et_al_Figure_3](https://github.com/rpyron/delim-SOM/assets/583099/9f28c7f0-0790-4a47-a7a4-25c98024a087)

_Original SOM results from Pyron et al. (2023)_

We provide a sample dataset and analysis for Seal Salamanders (_Desmognathus monticola_), which now represents two species in the Appalachian Mountains and Piedmont of the eastern United States, based on four datasets comprising a SNP matrix from Genotype-By-Sequencing (GBS) analysis, long/lat/elevation (xyz), environment (climate, hydrology, and ecoregion), and phenotype (17 linear morphometric measurements and larval spot count).

All of the code here is given in _'monticola_models.R'_. The genetic, spatial, and environmental data come from 71 individuals from 71 sites, while the phenotypic data are expanded to include up to 163 specimens from those sites, with the mean of each measurement taken after size correction. The allele frequencies come from a GBS matrix of 5,174 SNPs and 10,526 alleles after trimming to 80% completeness.

The climate variables are Level IV Ecoregion (https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states), HUC4 watershed (https://www.usgs.gov/national-hydrography/watershed-boundary-dataset), ENVIREM - monthCountByTemp10 (https://envirem.github.io/), and WorldClim - BIO15 (https://www.worldclim.org/data/bioclim.html). These are one-hot encoded (dummy variables by factor level) in the supplied dataset; the original values are in './data/seal_data.csv.'

The phenotype variables are 17 linear morphometric measurements to 0.01mm precision: SVL (snout-vent length), TL (tail length), AG (axilla-groin length), CW (chest width), FL (femur length [rather than hindlimb length]), HL (humerus length [rather than forelimb length]), SG (snout-gular length), TW(tail width at rear of vent), TO (length of third toe), FI (length of third finger), HW (head width), ED (eye diameter), IN (internarial distance), ES (eye-snout distance), ON (orbito-narial distance), IO (inter-orbital distance), and IC (inter-canthal distance). Here, I size-correct these by SVL using pooled groups ("population2") in 'GroupStruct' (Chan and Grismer 2021, 2022): https://github.com/chankinonn/GroupStruct, then take the mean by site. I also include the mean left and right larvat-spot counts by site for 66 individuals from 40 sites; the remaining rows are _'NA'_.

# Running the Code

With all of the code and data files in the package directories ('R' and 'data'), we can simply execute monticola_example.R from the base directory. The main pieces are as follows.

```
a <- read.structure("./data/seal_in_c90.str",
                    n.ind = 71,
                    n.loc = 7809,
                    onerowperind = FALSE,
                    col.lab = 1,
                    col.pop = 0,
                    col.others = 0,
                    row.marknames = 0,
                    NA.char = -9)

#Trim missingness
a = missingno(a, type = "loci", cutoff = 0.20)
a#trimmed to 20% missing data

#Convert allele frequences to matrix
alleles <- makefreq(a)
```

The alleles matrix can be in nearly any format, with individuals in rows and allele frequencies or counts in columns. Here, I am simply loading in the STRUCTURE-formatted file from ipyrad as a genind object in 'adegenet' (Jombart 2008: https://adegenet.r-forge.r-project.org/), trimming it to 20% missing data, converting the counts to frequencies, and converting it to a matrix.

```
dat <- read.csv("./data/seal_data.csv",header=T,row.names=1)
xyz <- dat[,2:4]

space <- as.matrix(read.csv("./data/seal_space.csv",header=T,row.names=1))
climate <- as.matrix(read.csv("./data/seal_climate.csv",header=T,row.names=1))
traits <- as.matrix(read.csv("./data/seal_traits.csv",header=T,row.names=1))
```

Similarly, these are just the long, lat, elevation, climate (including ecoregions and river drainages), and trait data from Pyron et al. (2023). The morphological is trimmed to the 163 specimens from the 71 genetic localities including larval spot counts, with linear traits corrected for allometry (by SVL) using 'GroupStruct' (Chan and Grismer 2022), taking the mean by site, and min-max normalizing the resulting matrix.

```
#Size of Grid
g <- round(sqrt(5*sqrt(length(rownames(alleles)))))#common rule of thumb

#Create an output grid of size sqrt(n)
som_grid <- somgrid(xdim = g,
                    ydim = g,
                    topo="hexagonal",
                    neighbourhood.fct = "gaussian")

#Number of Replicates - can increase if you like
n <- 100

#Number of steps - doesn't usually matter beyond ~100
m <- 100
```

Next, we set the parameters for the SOM estimates. We can evaluate _m_ using the hyperparameter exploration below, but these are good default values.

```
res <- Trait.SOM()
```

We use Trait.SOM() to produce our estimates. We could also use DNA.SOM(), Space.SOM(), or Climate.SOM() if we only had those layers. The output is then stored to local variables for later analysis.

```
plotLearning.Traits(res)
```

![image](https://github.com/rpyron/delim-SOM/assets/583099/2d4bfddc-0048-4844-af81-54a51be17ab0)

We can then plot our learning estimates across the runs. The shape of this learning curve (slow decline, then sudden plateau) is inherent in the way the algorithm learns; longer runs produce the same shape, rather than a longer plateau. The scale of each variable determines the location of its plateau; we expect each matrix to stabilize, but not necessarily to converge to the same relative distance to closest unit.

```
plotLayers(res)
```

![image](https://github.com/rpyron/delim-SOM/assets/583099/b109bd9f-8e6a-468f-8b1f-70853419fd5c)

Next, we can see the layer weights. Unsurprisingly, alleles dwarf everything else, but climate is more important than traits, and both are greater than space alone. Note that this does not match Pyron (2023) exactly; the factor levels for the categorical variables were not being handled properly before, but are now incorporated correctly with one-hot encoding. Check this in your data!

```
plotK(res)
```

![image](https://github.com/rpyron/delim-SOM/assets/583099/0f0e97bc-5053-4572-a001-9020f5b5ab57)

Then, we can see the optimal values of _K_. In this case, only _K_=2 was sampled across the 100 replicates.

```
set.seed(1)
labels <- match.labels(alleles)#get DAPC labels
q_mat <- match.k(res,labels)#get admixture coefficients

par(mfrow=c(1,1),
    mar=c(0,0,0,0))
xy <- xyz[,1:2]
maps::map(database = 'county', xlim = range(xy[,1]) + c(-0.5,0.5), ylim = range(xy[,2]) + c(-0.5,0.5), col="white")
map.axes()
maps::map(database = 'county', xlim = range(xy[,1]) + c(-0.5,0.5), ylim = range(xy[,2]) + c(-0.5,0.5), col="gray",add=T)
maps::map(database = 'state', xlim = range(xy[,1]) + c(-0.5,0.5), ylim = range(xy[,2]) + c(-0.5,0.5), add = T)
make.admix.pie.plot(q_mat,xy,layer.colors = k.cols,radii=2.5,add = T)
legend(-88,38,legend=c(expression(italic("D. cheaha")),
         expression(italic("D. monticola"))),cex=2,pt.bg=k.cols[2:1],pch=21)
map.scale(-81.2,31.1)
```

![image](https://github.com/rpyron/delim-SOM/assets/583099/fbcbfaad-ee9b-4ba4-b66e-e6037d4ace16)

We can also produce a basic sample map. The match.k() function uses a CLUMPP-like algorithm to synchronize the cluster labels to the DAPC results from earlier.

```
x <- q_mat[order(q_mat[,1]),]
z <- hclust(dist(x),"single")$order
make.structure.plot(admix.proportions = x[z,], 
                    sample.names = rownames(x[z,]), 
                    mar = c(8,4,2,2), 
                    layer.colors = k.cols, 
                    sort.by = 1)
```

![image](https://github.com/rpyron/delim-SOM/assets/583099/9fca61b7-412e-49ac-818b-f778136b9cd1)

A STRUCTURE-type plot, organized hierarchically by dominant cluster membership. Given the extensive differences between these two species in terms of genetics, geography, ecology, and phenotype, the species coefficients are sharply bimodal comapred to the individual ancestry coefficients estimated from the SNP matrix alone (see Pyron et al. 2023).

```
#Example outputs from one model#
plotModel(res)
```

![image](https://github.com/rpyron/delim-SOM/assets/583099/98378184-0158-4238-a58d-5676766d1448)

An example SOM plot looking at the results from one model in terms of sample assignment to cells, cell distances, and boundaries between cell clusters.

```
#Variable Importance
Trait.SOM.varImp(res)
```

![image](https://github.com/rpyron/delim-SOM/assets/583099/d2570754-e8c6-4f1b-b73e-596e8f432db0)

**New feature!** Variable importance estimates (codebook vectors/neuron weights on the interval [0,1]) for each layer, allowing you to identify the features with the greatest impact on assigning individuals to output grid cells. This takes the median estimate across cells for each input variable from the final model, and returns a named list thereof. The plot shows the features with varImp > 0.001 and is typically dominated by alleles, which are therefore not named individually but only counted instead. Individual patterns in allele importance could be broken down by extracting the names of important loci from the list. Similar functions are available as 'DNA.SOM.varImp(),' 'Space.SOM.varImp(),' and 'Climate.SOM.varImp().'

![Pyron_UML_Graphical_Abstract](https://github.com/rpyron/delim-SOM/assets/583099/f1a64348-832e-49f9-bf28-b6c81bf7a30f)

A nice summary figure for publication (from Pyron 2023)!

```
save.image(file="Trait_SuperSOM.RData")
```

Finally, we can save all of our results.

# Hyperparameters

Many default implementations of SOMs have paid little attention to the hyperparameters of neighborhood type, learning rate, and run length (number of steps), finding them to have relatively small impacts on learning outcomes (Wehrens and Kruisselbrink 2018; Pyron et al. 2023). Given the importance of accurately quantifying layer contributions, we expand on this concern here. Previous studies have found that Gaussian neighborhoods (rather than bubble or heuristic) and linear learning-rates typically yield optimal results under a wide range of conditions (Stefanovič and Kurasova 2011; Natita et al. 2016). Consequently, we employ these as our default conditions. 

![image](https://github.com/rpyron/delim-SOM/assets/583099/c1a352c3-b14f-4ffe-9e8e-326f16372eda)

The kohonen_hyper.R file contains a brief exploration of the hyperparameters, including the learning rates and run length. Generally speaking, these don't have much of an impact for rlen > 100 and alpha > 0.1 for the alleles matrix in the _D. monticola_ dataset. Longer runs (rlen/m) may have a small impact on final learning estimates and precision. These curves are estimated from 100 estimates and may normalize after additional replicates.

# Comparison to other methods

We can compare DNA-only SOM and trait-based SuperSOM estimates of individual ancestry (DNA-only) or species coefficients (SuperSOM) to the admixture values from sNMF presented by Pyron et al. (2023).

```
#Trait-based SuperSOM
res <- Trait.SOM()
q_mat <- match.k(res)#get admixture coefficients

#DNA-only SOM results
res1 <- DNA.SOM()
q_mat.DNA <- match.k(res1)#get admixture coefficients

#sNMF values from previous analysis
sNMF_q_mat <- read.csv("./seal_q.mat.csv",row.names=1)

par(mfrow=c(2,1),mar=c(1,4.5,0.5,2),mgp=c(2,0.5,0))
plot(sNMF_q_mat[,1],q_mat.DNA[,1],pch=21,bg=layer.cols[1],xaxt='n',xlab=NA,
     ylab="DNA SOM",xlim=c(0,1),ylim=c(0,1),cex=2)
axis(side = 1, at = c(0,0.2,0.4,0.6,0.8,1), labels = FALSE)
b.sp <- cor(sNMF_q_mat[,1],q_mat.DNA[,1],method=c("spearman"))
b.pe <- cor(sNMF_q_mat[,1],q_mat.DNA[,1],method=c("pearson"))
text(0.2,0.9,paste0("Spearman's = ",round(b.sp,2)))
text(0.2,0.85,paste0("Pearson's = ",round(b.pe,2)))
text(0.675,0.05,"a) Individual Ancestries",font=2,cex=1.25)

par(mar=c(4,4.5,0.5,2))
plot(sNMF_q_mat[,1],q_mat[,1],pch=21,bg=layer.cols[4],
     xlab="sNMF Admixture Estimates",ylab="Trait SuperSOM",xlim=c(0,1),ylim=c(0,1),cex=2)
c.sp <- cor(sNMF_q_mat[,1],q_mat[,1],method=c("spearman"))
c.pe <- cor(sNMF_q_mat[,1],q_mat[,1],method=c("pearson"))
text(0.2,0.9,paste0("Spearman's = ",round(c.sp,2)))
text(0.2,0.84,paste0("Pearson's = ",round(c.pe,2)))
text(0.675,0.05,"b) Species Coefficients",font=2,cex=1.25)
```

![image](https://github.com/rpyron/delim-SOM/assets/583099/ef19a3f1-553f-43f2-87ff-352c8102081d)

Estimates from the DNA-only SOM are linear between ~30-70% ancestry but essentially binary outside of that range, as parental cluster assignment is less variable in the tails. Species coefficients from the trait-based SuperSOM are essentially binary for all individual. This is unsurprising given that even admixed populations near the hybrid zone tend to be either montane or Piedmont and have the strongly diagnostic character of 4–5 versus 6–7 larval spots in D. monticola compared to D. cheaha (Pyron et al., 2023). I do note that the sample size here is relatively small, as individual ancestry coefficients were already sharply bimodal with relatively few hybrid or admixed individuals. One specimen with relatively high genomic ancestry is more mixed in the integrative model, due to its intermediate morphology and geographic position.

# Simulations

Some basic simulations to demonstrate desirable performance under a wide variety of conditions. Requires the 'simulMGF' package of Sikorska et al. (2013): https://cran.r-project.org/web/packages/simulMGF/index.html.
```
###########################
#Simulate a K=1 SNP matrix#
###########################
simGeno(71, 5174)
alleles <- simG/2
rownames(alleles) <- paste("a",1:71,sep="")
colnames(alleles) <- paste("snp",1:5174,sep="")
```
First, simulate a _K_=1 SNP matrix.

```
################################
#Get space, climate, and traits#
#from Desmognathus dataset     #
################################

#Sample data
dat <- read.csv("../data/seal_data.csv",header=T,row.names=1)
xyz <- dat[,2:4]

###SPATIAL, CLIMATIC, AND TRAIT DATA
space <- as.matrix(read.csv("../data/seal_space.csv",header=T,row.names=1))
climate <- as.matrix(read.csv("../data/seal_climate.csv",header=T,row.names=1))
traits <- as.matrix(read.csv("../data/seal_traits.csv",header=T,row.names=1))
```
Second, load the empirical space, climate, and trait layers.

```
##########
##########
#DNA ONLY#
##########
##########
res <- DNA.SOM()

#Plot Learning#
plotLearning.DNA(res)

#Optimize K#
plotK(res)
```

Run a DNA-only SOM to show strong support for _K_=1 using the standard model-selection criteria under minimum BIC.

```
##########
##########
#TRAITS  #
##########
##########
res1 <- Trait.SOM()

#Plot Learning#
plotLearning.Traits(res1)

#Layer Weights#
plotLayers(res1)

#Optimize K#
plotK(res1)
```

Repeat including the other empirical data layers.

```
############
############
#TRAITS/dim#
############
############
alleles <- alleles/1000#reduce signal from alleles

#get labels for different K values
labels <- data.frame(V1=rep(NA,dim(alleles)[1]),row.names = rownames(alleles))
for (i in 1:10){labels[,i] <- find.clusters(alleles,n.clust=i,n.pca = dim(alleles)[1])$grp}

###Kohonen maps###
res2 <- Trait.SOM()

#Plot Learning#
plotLearning.Traits(res2)

#Layer Weights#
plotLayers(res2)

#Optimize K#
plotK(res2)

#Sample Map#
labels <- match.labels(alleles)#get DAPC labels
q_mat <- match.k(res2)#get admixture coefficients

par(mfrow=c(1,1),
    mar=c(0,0,0,0))
xy <- xyz[,1:2]
maps::map(database = 'county', xlim = range(xy[,1]) + c(-0.5,0.5), ylim = range(xy[,2]) + c(-0.5,0.5), col="white")
map.axes()
maps::map(database = 'county', xlim = range(xy[,1]) + c(-0.5,0.5), ylim = range(xy[,2]) + c(-0.5,0.5), col="gray",add=T)
maps::map(database = 'state', xlim = range(xy[,1]) + c(-0.5,0.5), ylim = range(xy[,2]) + c(-0.5,0.5), add = T)
make.admix.pie.plot(q_mat,xy,layer.colors = k.cols,radii=2.5,add = T)
map.scale(-81.2,31.1)

#Structure Plot#
x <- q_mat[order(q_mat[,1]),]
z <- hclust(dist(x),"single")$order
make.structure.plot(admix.proportions = x[z,], 
                    sample.names = rownames(x[z,]), 
                    mar = c(8,4,2,2), 
                    layer.colors = k.cols, 
                    sort.by = 1)
```

Finally, reduce impact of null alleles to ~0 to demonstrate inclusion of signal from the three other empirical layers.

![image](https://github.com/rpyron/delim-SOM/assets/583099/047234fc-250c-4e5e-8b8e-96262703bd1f)

The method strongly supports _K_=1 when little genetic structure is present, even when space, climate, and traits vary (first and second columns). In contrast, when the signal of the alleles layer is reduced to ~0, the impact of the other three layers is reflected in the output, still estimating a roughly binary division corresponding to D. cheaha in the south and D. monticola in the north while sampling multiple possible ancestries (third and fourth columns). These simulations reveal that the SOM/SuperSOM approach can detect _K_=1, does not over-split, reflects contributions from all layers with signal, and does not allow layer size (e.g., large-scale genetic matrices) to overwhelm other datasets. Note that these results are slightly different from Pyron (2023) as the factor levels are now handled more appropriately with one-hot encoding.

# Possible Improvements

I have not explored different data types for molecular datasets (e.g., allele counts, structural variants). As currently implemented, it can support _K_=1 (see Janes et al. 2017), but this should be explored in depth for this and other methods (e.g., Derkarabetian et al. 2019). Another possibility is conditioning this step directly on the delta-BIC metric, similar to the delta-K method of Evanno et al. (2005) but on a per-step basis rather than averaged across runs. This may end up being essentially identical to the 'diffNgroup' method of large vs. small changes from Jombart et al. (2010) that is currently used. Finally, adding support for the fifth conservation layer is simply an extension of the *.SOM() functions to add containers for a fifth set of variables, along with the attendant functions for plotting learning and variable importance.

**If you have a dataset that's amenable to a fifth conservation layer, please email me (rpyron@gwu.edu)! I would love to work with you to get this implemented if it's useful.**


# References
K. O. Chan and L. L. Grismer (2021). A standardized and statistically defensible framework for quantitative morphological analyses in taxonomic studies. _Zootaxa_, 5023: 293-300.

S. Derkarabetian, S. Castillo, P. K. Koo, S. Ovchinnikov and M. Hedin (2019). A demonstration of unsupervised machine learning in species delimitation. _Molecular Phylogenetics and Evolution_, 139: 106562. https://doi.org/10.1016/j.ympev.2019.106562

G. Evanno, S. Regnaut and J. Goudet (2005). Detecting the number of clusters of individuals using the software STRUCTURE: a simulation study. _Molecular Ecology_, 14: 2611-2620. https://doi.org/10.1111/j.1365-294X.2005.02553.x

M. Jakobsson and N. A. Rosenberg (2007). CLUMPP: a cluster matching and permutation program for dealing with label switching and multimodality in analysis of population structure. _Bioinformatics_, 23: 1801-1806. https://doi.org/10.1093/bioinformatics/btm233

J. K. Janes, J. M. Miller, J. R. Dupuis, R. M. Malenfant, J. C. Gorrell, C. I. Cullingham and R. L. Andrew (2017). The K=2 conundrum. _Molecular Ecology_, 26: 3594-3602. https://doi.org/10.1111/mec.14187

T. Jombart (2008). adegenet: a R package for the multivariate analysis of genetic markers. _Bioinformatics_, 24: 1403-1405. https://doi.org/10.1093/bioinformatics/btn129

W. Natita, W. Wiboonsak and S. Dusadee (2016). Appropriate learning rate and neighborhood function of Self-organizing Map (SOM) for specific humidity pattern classification over Southern Thailand. _International Journal of Modeling and Optimization_, 6: 61-65. https://doi.org/10.7763/IJMO.2016.V6.504

R. A. Pyron, K. A. O’Connell, S. C. Duncan, F. T. Burbrink and D. A. Beamer (2023). Speciation hypotheses from phylogeographic delimitation yield an integrative taxonomy for Seal Salamanders (_Desmognathus monticola_). _Systematic Biology_, 72: 179-197. https://doi.org/10.1093/sysbio/syac065

K. Sikorska, E. Lesaffre, P. F. Groenen et al. (2013). GWAS on your notebook: fast semi-parallel linear and logistic regression for genome-wide association studies. _BMC Bioinformatics_, 14: 166. https://doi.org/10.1186/1471-2105-14-166

P. Stefanovič and O. Kurasova (2011). Influence of learning rates and neighboring functions on Self-Organizing Maps. _Advances in Self-Organizing Maps_, 6731: 141-150. https://doi.org/10.1007/978-3-642-21566-7_14

R. Wehrens and L. M. Buydens (2007). Self-and super-organizing maps in R: the Kohonen package. _Journal of Statistical Software_, 21: 1-19. https://doi.org/10.18637/jss.v021.i05

R. Wehrens and J. Kruisselbrink (2018). Flexible Self-Organizing Maps in kohonen 3.0. _Journal of Statistical Software_, 87: 1-18. https://doi.org/10.18637/jss.v087.i07
