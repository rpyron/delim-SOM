# delim-SOM
Using multi-layer Kohonen Self-Organizing Maps ("SuperSOMs") to delimit species and produce integrative taxonomies using Unsupervised Machine Learning (UML). This is accomplished primarily using the R package 'kohonen' (Wehrens and Buydens 2007): https://cran.r-project.org/web/packages/kohonen/index.html.

This repository expands the use of Kohonen maps as described in Pyron et al. (2023). It uses multi-layer Self-Organizing Maps ("SuperSOMs") in the R package 'kohonen' to delimit species based on allelic, spatial, climatic, and phenotypic data.

The contribution of each layer to the final model output is recorded, along with the clustering assignment of each individual over multiple learning replicates. The results therefore mirror a 'STRUCTURE'-type analysis including admixture estimates, but represent an unified delimitation model that incorporates various dimensions of ecological and evolutionary divergence for integative taxonomy. If only allelic data are used with a 'DNA.SOM()' model, then the assignment probabilities approximate individual ancestry coefficients. If multiple layers are used, we treat them as "species coefficients," which might be useful for testing a variety of ecological and evolutionary hypotheses.

The requisite functions are in the 'kohonen_code.R' file, which loads the various dependencies:

```
library(adegenet); library(maps); library(viridis);library(scales)
library(LEA); library(conStruct);library(poppr);library(kohonen)
library(lsr);library(combinat);library(caret);library(GroupStruct)
library(vcfR);library(dartR)
```

Some of these may have to be installed manually or from various non-CRAN sources.

Overall, the method is extremely flexible and can take almost any data type or format, as long as it is introduced as a matrix in R. The matrices must be added in order, named alleles, space, climate, and traits. Adding each of those matrices in sequence allows one to run SOMs based on DNA, DNA + xyz, DNA + xyz + environment, and DNA + xyz + environment + phenotypes. The primary requirement is to have individuals in rows in the same order in each matrix, and variables in columns, with <10% missing data and the same set of individuals in each matrix. I also min-max normalize the space, climate, and traits matrices to be on the same scale as the allele frequencies. You could modify the code to allow >10% missing data (maxNA.frac) if necessary, but the effects are unknown.

# Run this on your data

SuperSOMs require (up to) four data layers as input matrices, called 'alleles,' 'space,' 'climate,' and 'traits.' These should each have the same number of rows (individuals, specimens, or populations), and any number of columns (however many SNPs or other variables you have). I suggest using allele frequencies for 'alleles,' and normalizing the other matrices to [0,1] for the same scale.

```
alleles <- matrix()#Molecular data as allele frequencies per locus
space <- matrix()#Typically long, lat, and elevation
climate <- matrix()#Relevant environmental data
traits <- matrix()#Phenotypic data (e.g., morphometrics, behavior)
```

You will also want baseline clustering estimates from the molecular data to guide label synchronization:

```
#############################
###Baseline DAPC assignments#
#for synchronizing clusters #
#############################

#get labels for different K values
labels <- data.frame(V1=rep(NA,dim(alleles)[1]),row.names = rownames(alleles))
for (i in 1:10){labels[,i] <- find.clusters(a,n.clust=i,n.pca = dim(alleles)[1])$grp}
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
```

Look at a representative SOM grid:

```
#Example outputs from one model#
par(mfrow=c(2,1))

#Cell distances #
plot(res$som_model, type="dist.neighbours", main = "", palette.name = viridis)
title("SOM neighbour distances",line=-0.1)

#SOM Clusters#
som.cols <- setNames(k.cols[labels[,max(res$c_mat)]],res$cluster_assignment)#Get colors to match original SOM clusters
som.cols <- unique(som.cols[sort(names(som.cols))])#Set to refactored labels

#plot cluster
plot(res$som_model, shape="straight", type="mapping", bgcol = som.cols[res$som_cluster], main = "", pch=19, col="red")
add.cluster.boundaries(res$som_model, res$som_cluster,col="red");title("SOM clusters",line=-0.1)
```

Review a map, where _'xyz'_ is your long/lat/elevation matrix:

```
#Sample Map#
q_mat <- match.k(res)#get admixture coefficients

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

The climate variables are Level IV Ecoregion (https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states), HUC4 watershed (https://www.usgs.gov/national-hydrography/watershed-boundary-dataset), ENVIREM - monthCountByTemp10 (https://envirem.github.io/), and WorldClim - BIO15 (https://www.worldclim.org/data/bioclim.html).

The phenotype variables are 17 linear morphometric measurements to 0.01mm precision: SVL (snout-vent length), TL (tail length), AG (axilla-groin length), CW (chest width), FL (femur length [rather than hindlimb length]), HL (humerus length [rather than forelimb length]), SG (snout-gular length), TW(tail width at rear of vent), TO (length of third toe), FI (length of third finger), HW (head width), ED (eye diameter), IN (internarial distance), ES (eye-snout distance), ON (orbito-narial distance), IO (inter-orbital distance), and IC (inter-canthal distance). Here, I size-correct these by SVL using pooled groups ("population2") in 'GroupStruct' (Chan and Grismer 2021, 2022): https://github.com/chankinonn/GroupStruct, then take the mean by site. I also include the mean left and right larvat-spot counts by site for 66 individuals from 40 sites; the remaining rows are _'NA'_.

# Running the Code

With all of the files (kohonen_code.R, monticola_models.R, seal_in.str, seal_clim.csv, and seal_morph.csv) in the same directory, we can simply execute monticola_models.R. The main pieces are as follows.

```
a <- read.structure("./seal_in.str",
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

#Convert genind object from adegenet to allele frequencies
struc <- makefreq(a)

#Convert allele frequences to matrix
alleles <- matrix(unlist(as.numeric(struc)), nrow=nrow(struc))
rownames(alleles) <- rownames(struc);colnames(alleles) <- colnames(struc)
```

The alleles matrix can be in nearly any format, with individuals in rows and allele frequencies or counts in columns. Here, I am simply loading in the STRUCTURE-formatted file from ipyrad as a genind object in 'adegenet' (Jombart 2008: https://adegenet.r-forge.r-project.org/), trimming it to 20% missing data, converting the counts to frequencies, and converting it to a matrix.

```
dat <- read.csv("./seal_clim.csv",row.names = "Specimen")
xyz <- dat[,c("LONG","LAT","elevation")]
space <- data.frame(lon=xyz$LONG,lat=xyz$LAT,elev=xyz$elev)
space <- apply(space,2,minmax)
climate <- apply(dat[,c(5:9)],2,minmax)#These variables were identified as most important
```

Similarly, these are just the long, lat, elevation, and climate data from Pyron et al. (2023).

```
morph <- read.csv("./seal_morph.csv",row.names=1)#Read in trait data, 163 specimens from 71 sites with 17 measurements 
morph_log <- data.frame(pop="seal",exp(morph[,2:18]))#un-log the measurements for GroupStruct
morph_allom <- allom(morph_log,"population2")#correct for allometry using GroupStruct
morph_mean <- aggregate(morph_allom[,-1],list(morph$pop),mean)#take mean by locality - this is a simplistic approach
morph_norm <- apply(morph_mean[,-1],2,minmax)#could also use PC1-3 or similar transformation

#larval spot count
spots <- read.csv("seal_spots.csv",row.names=1)
spot_mean <- aggregate(sqrt(spots[,1:2]),list(spots$pop),mean)
spot_norm <- spot_mean[,-1]
spot_norm[-which(is.na(spot_mean[,2:3])),] <- apply(na.omit(spot_mean)[,-1],2,minmax)

#merge into traits
traits <- as.matrix(cbind(morph_norm,spot_norm))
```

Finally, the morphological dataset from Pyron et al. (2023) trimmed to the 163 specimens from the 71 genetic localities. The log-transformed data are read in, exponentiated for analysis, assigned the same "species" of "seal," corrected for allometry using 'GroupStruct,' taking the mean by site, and min-max normalizing the resulting matrix.

```
labels <- data.frame(V1=rep(NA,dim(alleles)[1]),row.names = rownames(alleles))
for (i in 1:10){labels[,i] <- find.clusters(a,n.clust=i,n.pca = dim(alleles)[1])$grp}
```

Next, we will produce a baseline clustering estimate under various _K_ values to more easily synchronize cluster labels later on using the CLUMPP-like algorithm (Jakobson and Rosenberg 2007).

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

![image](https://github.com/rpyron/delim-SOM/assets/583099/50b2dd49-b89c-4881-a710-3dd8ef5465fb)

We can then plot our learning estimates across the runs. The shape of this learning curve (slow decline, then sudden plateau) is inherent in the way the algorithm learns; longer runs produce the same shape, rather than a longer plateau. The scale of each variable determines the location of its plateau; we expect each matrix to stabilize, but not necessarily to converge to the same relative distance to closest unit.

```
plotLayers(res)
```

![image](https://github.com/rpyron/delim-SOM/assets/583099/1ed8cea5-cc42-4902-80c2-cf639e83bd75)

Next, we can see the layer weights. Unsurprisingly, alleles dwarf everything else, but traits are more important than climate, and both are greater than space alone.

```
plotK(res)
```

![image](https://github.com/rpyron/delim-SOM/assets/583099/8d00f902-2436-443f-b397-6bed47e21e9a)

Then, we can see the optimal values of _K_. In this case, only _K_=2 was sampled across the 100 replicates.

```
q_mat <- match.k(res)#get admixture coefficients

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

![image](https://github.com/rpyron/delim-SOM/assets/583099/8f020274-4b86-4043-8d6e-e1f07ee966c1)

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

![image](https://github.com/rpyron/delim-SOM/assets/583099/d9d472b1-2f9d-4cbe-914b-bb897cda5faf)

A STRUCTURE-type plot, organized hierarchically by dominant cluster membership. Given the extensive differences between these two species in terms of genetics, geography, ecology, and phenotype, the species coefficients are sharply bimodal comapred to the individual ancestry coefficients estimated from the SNP matrix alone (see Pyron et al. 2023).

```
png("./Pyron_UML_Graphical_Abstract.png",1328,531)
#SEE CODE IN FILE!
dev.off()
```

![Pyron_UML_Graphical_Abstract](https://github.com/rpyron/delim-SOM/assets/583099/f1a64348-832e-49f9-bf28-b6c81bf7a30f)

A nice summary figure for publication!

```
save.image(file="Trait_SuperSOM.RData")
```

Finally, we can save all of our results.

# Hyperparameters

Many default implementations of SOMs have paid little attention to the hyperparameters of neighborhood type, learning rate, and run length (number of steps), finding them to have relatively small impacts on learning outcomes (Wehrens and Kruisselbrink 2018; Pyron et al. 2023). Given the importance of accurately quantifying layer contributions, we expand on this concern here. Previous studies have found that Gaussian neighborhoods (rather than bubble or heuristic) and linear learning-rates typically yield optimal results under a wide range of conditions (Stefanovič and Kurasova 2011; Natita et al. 2016). Consequently, we employ these as our default conditions. 

![image](https://github.com/rpyron/delim-SOM/assets/583099/c1a352c3-b14f-4ffe-9e8e-326f16372eda)

The kohonen_hyper.R file contains a brief exploration of the hyperparameters, including the learning rates and run length. Generally speaking, these don't have much of an impact for rlen > 100 and alpha > 0.1 for the alleles matrix in the _D. monticola_ dataset. Longer runs (rlen/m) may have a small impact on final learning estimates and precision. These curves are estimated from 100 estimates and may normalize after additional replicates.

# Simulations

Some basic simulations to demonstrate desirable performance under a wide variety of conditions.
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

###SPATIAL & CLIMATIC DATA
#read in data
dat <- read.csv("./seal_clim.csv",row.names = "Specimen")
xyz <- dat[,c("LONG","LAT","elevation")]
space <- data.frame(lon=xyz$LONG,lat=xyz$LAT,elev=xyz$elev)
space <- apply(space,2,minmax)
climate <- apply(dat[,c(5:9)],2,minmax)#These variables were identified as most important by SDM

###PHENOTYPIC DATA
#linear morphometrics
morph <- read.csv("./seal_morph.csv",row.names=1)#Read in trait data, 163 specimens from 71 sites with 17 measurements 
morph_log <- data.frame(pop="seal",exp(morph[,2:18]))#un-log the measurements for GroupStruct
morph_allom <- allom(morph_log,"population2")#correct for allometry using GroupStruct
morph_mean <- aggregate(morph_allom[,-1],list(morph$pop),mean)#take mean by locality
morph_norm <- apply(morph_mean[,-1],2,minmax)#could also use PC1-3 or similar transformation

#larval spot count
spots <- read.csv("seal_spots.csv",row.names=1)
spot_mean <- aggregate(sqrt(spots[,1:2]),list(spots$pop),mean)
spot_norm <- spot_mean[,-1]
spot_norm[-which(is.na(spot_mean[,2:3])),] <- apply(na.omit(spot_mean)[,-1],2,minmax)

#merge into traits
traits <- as.matrix(cbind(morph_norm,spot_norm))
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

![Figure_5](https://github.com/rpyron/delim-SOM/assets/583099/420212e0-185b-4e0e-a5bf-b1cab55802e0)

# Possible Improvements

I have not explored different data types for molecular datasets (e.g., allele counts, structural variants). I also haven't explored the impact of normalization, although this should be minimal. Selection of _K_ for each run might deserve a bit more scrutiny. Some possibilities are discussed here: https://bgstieber.github.io/post/an-introduction-to-the-kmeans-algorithm/. As currently implemented, it can support _K_=1 (see Jaynes et al. 2017), but this should be explored in depth for this and other methods (e.g., Derkarabetian et al. 2019).


# References

Chan, K.O. & Grismer, L. L. (2021). A standardized and statistically defensible framework for quantitative morphological analyses in taxonomic studies. Zootaxa, 5023: 293-300.

Chan, K.O. and Grismer, L.L., 2022. GroupStruct: an R package for allometric size correction. Zootaxa, 5124(4), pp.471-482.

Derkarabetian, S., Castillo, S., Koo, P.K., Ovchinnikov, S. and Hedin, M., 2019. A demonstration of unsupervised machine learning in species delimitation. Molecular phylogenetics and evolution, 139, p.106562.

Jakobsson, M. and Rosenberg, N.A., 2007. CLUMPP: a cluster matching and permutation program for dealing with label switching and multimodality in analysis of population structure. Bioinformatics, 23(14), pp.1801-1806.

Janes, J.K., Miller, J.M., Dupuis, J.R., Malenfant, R.M., Gorrell, J.C., Cullingham, C.I. and Andrew, R.L., 2017. The K= 2 conundrum. Molecular Ecology, 26(14), pp.3594-3602.

Jombart, T., 2008. adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24(11), pp.1403-1405.

Natita W., Wiboonsak W., Dusadee S. 2016. Appropriate Learning Rate and Neighborhood Function of Self-organizing Map (SOM) for Specific Humidity Pattern Classification over Southern Thailand. IJMO. 6:61–65.

Pyron, R.A., O’Connell, K.A., Duncan, S.C., Burbrink, F.T. and Beamer, D.A., 2023. Speciation hypotheses from phylogeographic delimitation yield an integrative taxonomy for Seal Salamanders (Desmognathus monticola). Systematic Biology, 72(1), pp.179-197.

Stefanovič P., Kurasova O. 2011. Influence of Learning Rates and Neighboring Functions on Self-Organizing Maps. Advances in Self-Organizing Maps.:141–150.

Wehrens, R. and Buydens, L.M., 2007. Self-and super-organizing maps in R: the Kohonen package. Journal of Statistical Software, 21, pp.1-19.

Wehrens R., Kruisselbrink J. 2018. Flexible Self-Organizing Maps in kohonen 3.0. J. Stat. Soft. 87.

