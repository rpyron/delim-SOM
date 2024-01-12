# delim-SOM
This package uses multi-layer Kohonen Self-Organizing Maps ("SuperSOMs") to delimit species and produce integrative taxonomies using Unsupervised Machine Learning (UML) as decribed in Pyron (2023). This repository expands the use of single-layer SOMs as described in Pyron et al. (2023). It relies on the R package 'kohonen' (Wehrens and Buydens 2007) to delimit species based on allelic, spatial, climatic, and phenotypic data.

The contribution of each layer to the final model output is recorded, along with the clustering assignment of each individual over multiple learning replicates. The results therefore mirror a 'STRUCTURE'-type analysis including admixture estimates, but represent an unified delimitation model that incorporates various dimensions of ecological and evolutionary divergence for integative taxonomy. If only allelic data are used with a 'DNA.SOM()' model, then the assignment probabilities approximate individual ancestry coefficients. If multiple layers are used, we treat them as "species coefficients," which might be useful for testing a variety of ecological and evolutionary hypotheses.

The requisite functions are in the './R/kohonen_code.R' file, which loads the various dependencies:

```
source("./R/kohonen_code.R")
set.seed(1)

library(adegenet);library(maps);library(scales)
library(conStruct);library(poppr);library(kohonen)
library(lsr);library(combinat);library(viridis)
set.seed(1)
```

Some of these may have to be installed manually or from various non-CRAN sources.

Overall, the method is extremely flexible and can take almost any data type or format, as long as it is introduced as a matrix in R. The matrices must be added in order, named alleles, space, climate, and traits. Adding each of those matrices in sequence allows one to run SOMs based on DNA, DNA + xyz, DNA + xyz + environment, and DNA + xyz + environment + phenotypes. 

**I also have it set to delimit a maximum of 10 species;** this can be changed by altering the code in various places (email me if needed: rpyron@gwu.edu), but it's unknown how the method will perform at larger scales. The primary requirement is to have individuals in rows in the same order in each matrix, and variables in columns, with <90% missing data and the same set of individuals in each matrix. I also min-max normalize the space, climate, and traits matrices to be on the same scale as the allele frequencies. You could modify the code to allow different missing data percentages (maxNA.frac) if necessary, but the effects are unknown.

# Run this on your data

SuperSOMs require (up to) four data layers as input matrices, called 'alleles,' 'space,' 'climate,' and 'traits.' These should each have the same number of rows (individuals, specimens, or populations), and any number of columns (however many SNPs or other variables you have). I suggest using allele frequencies for 'alleles,' and normalizing the other matrices to [0,1] to the same scale. This would include one-hot encoding factors. The available functions are 'DNA.SOM(),' 'Space.SOM(),' 'Climate.SOM(),' and 'Trait.SOM().'

```
alleles <- matrix()#Molecular data as allele frequencies per locus
space <- matrix()#Typically long, lat, and elevation
climate <- matrix()#Relevant environmental data
traits <- matrix()#Phenotypic data (e.g., morphometrics, behavior)
```

You will also want baseline clustering estimates from the molecular data to guide cluster label synchronizing later on using the CLUMPP-like algorithm (Jakobson and Rosenberg 2007). This is achieved using the match.k() function (see below) on the 'alleles' object you've created. A simple way to do this is just to create a genind object called 'a' and use that for 'alleles' such as:

```
a <- df2genind(X)#Where X is a SNP matrix imported into R as a data.frame
alleles <- makefreq(a)#Calculates allele frequencies from the genind object 'a'
labels <- match.labels(alleles)#get DAPC labels
```


Then, construct a SOM grid:

```
###Parameters for runs
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

Run a SuperSOM:

```
##############
###Run SOMs###
##############
res <- Trait.SOM()
```

Visualize the output:

 ```
 #Plot Learning#
plotLearning.Traits(res)

#Layer Weights#
plotLayers(res)

#Optimize K#
plotK(res)

#Look at a representative SOM grid#
plotModels(res)
```

Review a map, where _'xyz'_ is your long/lat/elevation matrix:

```
#Get species coefficients#
q_mat <- match.k(res,labels)

#Sample Map#
par(mar=c(0,0,0,0))
xy <- xyz[,1:2]
maps::map(database = 'world', xlim = range(xy[,1]) + c(-1,1), ylim = range(xy[,2]) + c(-1,1), col="white")
map.axes()
maps::map(database = 'world', xlim = range(xy[,1]) + c(-1,1), ylim = range(xy[,2]) + c(-1,1), add = T)
make.admix.pie.plot(q_mat,xy,layer.colors = k.cols,radii=2.5,add = T)
map.scale()
```

And a STRUCTURE-type barplot:
```
make.structure.plot(admix.proportions = q_mat, 
                    sample.names = rownames(q_mat), 
                    mar = c(8,4,2,2), 
                    layer.colors = k.cols, 
                    sort.by = 1)
```

# Example: _Desmognathus monticola_, the Seal Salamander

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

**New feature!** Variable importance estimates (codebook vectors/neuron weights on the interval [0,1]) for each layer, allowing you to identify the features with the greatest impact on cell assignment - i.e., delimitation clustering. This takes the median estimate across cells for each input variable from the final model, and returns a named list thereof. The plot shows the features with varImp > 0.001 and is typically dominated by alleles, which are therefore not named individually but only counted instead. Individual patterns in allele importance could be broken down by extracting the names of important loci from the list. Similar functions are available as 'DNA.SOM.varImp(),' 'Space.SOM.varImp(),' and 'Climate.SOM.varImp().'


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

I have not explored different data types for molecular datasets (e.g., allele counts, structural variants). I also haven't explored the impact of normalization, although this should be minimal. Selection of _K_ for each run might deserve a bit more scrutiny. Some possibilities are discussed here: https://bgstieber.github.io/post/an-introduction-to-the-kmeans-algorithm/. As currently implemented, it can support _K_=1 (see Janes et al. 2017), but this should be explored in depth for this and other methods (e.g., Derkarabetian et al. 2019). Another possibility is conditioning this step directly on the delta-BIC metric, similar to the delta-K method of Evanno et al. (2005) but on a per-step basis rather than averaged across runs. This may end up being essentially identical to the 'diffNgroup' method of large vs. small changes from Jombart et al. (2010) that is currently used. Finally, adding support for the fifth conservation layer is simply an extension of the *.SOM() functions to add containers for a fifth set of variables, along with the attendant functions for plotting learning.

**If you have a dataset that's amenable to a fifth conservation layer, please email me (rpyron@gwu.edu)! I would love to work with you to get this implemented if it's useful.**


# References

Chan, K.O. & Grismer, L. L. (2021). A standardized and statistically defensible framework for quantitative morphological analyses in taxonomic studies. Zootaxa, 5023: 293-300.

Chan, K.O. and Grismer, L.L., 2022. GroupStruct: an R package for allometric size correction. Zootaxa, 5124(4), pp.471-482.

Derkarabetian, S., Castillo, S., Koo, P.K., Ovchinnikov, S. and Hedin, M., 2019. A demonstration of unsupervised machine learning in species delimitation. Molecular phylogenetics and evolution, 139, p.106562.

Evanno, G., Regnaut, S. and Goudet, J., 2005. Detecting the number of clusters of individuals using the software STRUCTURE: a simulation study. Molecular ecology, 14(8), pp.2611-2620.

Jakobsson, M. and Rosenberg, N.A., 2007. CLUMPP: a cluster matching and permutation program for dealing with label switching and multimodality in analysis of population structure. Bioinformatics, 23(14), pp.1801-1806.

Janes, J.K., Miller, J.M., Dupuis, J.R., Malenfant, R.M., Gorrell, J.C., Cullingham, C.I. and Andrew, R.L., 2017. The K=2 conundrum. Molecular Ecology, 26(14), pp.3594-3602.

Jombart, T., 2008. adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24(11), pp.1403-1405.

Natita W., Wiboonsak W., Dusadee S. 2016. Appropriate Learning Rate and Neighborhood Function of Self-organizing Map (SOM) for Specific Humidity Pattern Classification over Southern Thailand. IJMO. 6:61–65.

Pyron, R.A., O’Connell, K.A., Duncan, S.C., Burbrink, F.T. and Beamer, D.A., 2023. Speciation hypotheses from phylogeographic delimitation yield an integrative taxonomy for Seal Salamanders (Desmognathus monticola). Systematic Biology, 72(1), pp.179-197.

Sikorska, K., Lesaffre, E., Groenen, P.F. et al. GWAS on your notebook: fast semi-parallel linear and logistic regression for genome-wide association studies. BMC Bioinformatics 14, 166 (2013). https://doi.org/10.1186/1471-2105-14-166

Stefanovič P., Kurasova O. 2011. Influence of Learning Rates and Neighboring Functions on Self-Organizing Maps. Advances in Self-Organizing Maps.:141–150.

Wehrens, R. and Buydens, L.M., 2007. Self-and super-organizing maps in R: the Kohonen package. Journal of Statistical Software, 21, pp.1-19.

Wehrens R., Kruisselbrink J. 2018. Flexible Self-Organizing Maps in kohonen 3.0. J. Stat. Soft. 87.

