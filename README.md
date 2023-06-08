# delim-SOM
Using multi-layer Kohonen Self-Organizing Maps ("SuperSOMs") to delimit species and produce integrative taxonomies using Unsupervised Machine Learning (UML). This is accomplished primarily using the R package 'kohonen' (Wehrens and Buydens 2007): https://cran.r-project.org/web/packages/kohonen/index.html.

This repository expands the use of Kohonen maps as described in Pyron et al. (2023). It uses multi-layer Self-Organizing Maps ("SuperSOMs") in the R package 'kohonen' to delimit species based on allelic, spatial, climatic, and phenotypic data.

The contribution of each layer to the final model output is recorded, along with the clustering assignment of each individual over multiple learning replicates. The results therefore mirror a 'STRUCTURE'-type analysis including admixture estimates, but represent an unified delimitation model that incorporates various dimensions of ecological and evolutionary divergence for integative taxonomy. If only allelic data are used with a 'DNA.SOM()' model, then the assignment probabilities approximate individual ancestry coefficients. If multiple layers are used, we treat them as "species coefficients," which might be useful for testing a variety of ecological and evolutionary hypotheses.

The requisite functions are in the 'kohonen_code.R' file, which loads the various dependencies:

```
library(adegenet); library(maps); library(viridis);library(scales)
library(LEA); library(conStruct);library(poppr); library(kohonen)
library(lsr);library(combinat);library(caret);library(elevatr)
library(caret);library(GroupStruct);library(vcfR);library(dartR)
```

Some of these may have to be installed manually or from various non-CRAN sources.

Overall, the method is extremely flexible and can take almost any data type or format, as long as it is introduced as a matrix in R. The matrices must be added in order, named alleles, space, climate, and traits. Adding each of those matrices in sequence allows one to run SOMs based on DNA, DNA + xyz, DNA + xyz + environment, and DNA + xyz + environment + phenotypes. The primary requirement is to have individuals in rows in the same order in each matrix, and variables in columns, with <10% missing data and the same set of individuals in each matrix. I also min-max normalize the space, climate, and traits matrices to be on the same scale as the allele frequencies. You could modify the code to allow >10% missing data (maxNA.frac) if necessary, but the effects are unknown.

# Example: _Desmognathus monticola_, the Seal Salamander

![Pyron_et_al_Figure_3](https://github.com/rpyron/delim-SOM/assets/583099/9f28c7f0-0790-4a47-a7a4-25c98024a087)

_Original SOM results from Pyron et al. (2023)_

We provide a sample dataset and analysis for Seal Salamanders (_Desmognathus monticola_), which now represents two species in the Appalachian Mountains and Piedmont of the eastern United States, based on four datasets comprising a SNP matrix from Genotype-By-Sequencing (GBS) analysis, lat/long/elevation (xyz), environment (climate, hydrology, and ecoregion), and phenotype (17 linear morphometric measurements and larval spot count).

The genetic, spatial, and environmental data come from 71 individuals from 71 sites, while the phenotypic data are expanded to include up to 163 specimens from those sites, with the mean of each measurement taken after size correction. The allele frequencies come from a GBS matrix of 5,174 SNPs and 10,526 alleles after trimming to 80% completeness.

The climate variables are Level IV Ecoregion (https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states), HUC4 watershed (https://www.usgs.gov/national-hydrography/watershed-boundary-dataset), ENVIREM - monthCountByTemp10 (https://envirem.github.io/), and WorldClim - BIO15 (https://www.worldclim.org/data/bioclim.html).

The phenotype variables are 17 linear morphometric measurements to 0.01mm precision: SVL (snout-vent length), TL (tail length), AG (axilla-groin length), CW (chest width), FL (femur length [rather than hindlimb length]), HL (humerus length [rather than forelimb length]), SG (snout-gular length), TW(tail width at rear of vent), TO (length of third toe), FI (length of third finger), HW (head width), ED (eye diameter), IN (internarial distance), ES (eye-snout distance), ON (orbito-narial distance), IO (inter-orbital distance), and IC (inter-canthal distance). Here, I size-correct these by SVL using pooled groups ("population2") in 'GroupStruct' (Chan and Grismer 2021, 2022): https://github.com/chankinonn/GroupStruct, then take the mean by site.

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

##Dimensions
locmiss = propTyped(a, by = "loc")
barplot(sort(1-locmiss), ylim = c(0,1), ylab = "Complete genotypes (proportion)", xlab = "Locus", las = 2, cex.names = 0.7)
a = missingno(a, type = "loci", cutoff = 0.20);abline(h=0.2,col="red",lty=2)
a#trimmed to 20% missing data

#Convert genind object from adegenet to allele frequencies
struc <- makefreq(a)

#Convert allele frequences to matrix
alleles <- matrix(unlist(as.numeric(struc)), nrow=nrow(struc))
```

The alleles matrix can be in nearly any format, with individuals in rows and allele frequencies or counts in columns. Here, I am simply loading in the STRUCTURE-formatted file from ipyrad as a genind object in 'adegenet' (Jombart 2008: https://adegenet.r-forge.r-project.org/), trimming it to 20% missing data, converting the counts to frequencies, and converting it to a matrix.

```
#read in data
dat <- read.csv("./seal_clim.csv",row.names = "Specimen")
xyz <- dat[,c("LONG","LAT","elevation")]
space <- data.frame(lon=xyz$LONG,lat=xyz$LAT,elev=xyz$elev)
space <- apply(space,2,minmax)
climate <- apply(dat[,c(5:9)],2,minmax)#These variables were identified as most important
```

Similarly, these are just the lat, long, elevation, and climate data from Pyron et al. (2023).

```
#linear morphometrics
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
#choose a max K
k.max <- 10

#get labels for different K values
labels <- data.frame(V1=rep(NA,dim(a$tab)[1]),row.names = rownames(a$tab))
for (i in 1:k.max){labels[,i] <- find.clusters(a,n.clust=i,n.pca = dim(a$tab)[1])$grp}
```

Next, we will produce a baseline clustering estimate under various _K_ values to more easily synchronize cluster labels later on using the CLUMPP-like algorithm (Jakobson and Rosenberg 2007).

```
#Size of Grid
g <- round(sqrt(5*sqrt(length(rownames(a$tab)))))#common rule of thumb

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

#Output
c_mat              <- res$c_mat#Get classifiers
d_mat              <- res$d_mat#Get weights
l_mat1             <- res$l_mat1#Get learning
l_mat2             <- res$l_mat2#Get learning
l_mat3             <- res$l_mat3#Get learning
l_mat4             <- res$l_mat4#Get learning
w_mat              <- res$w_mat#Get WSS
som_model          <- res$som_model#Get last SOM
som_cluster        <- res$som_cluster#Get last clusters
cluster_assignment <- res$cluster_assignment#Get last assignments
```

We use Trait.SOM() to produce our estimates. We could also use DNA.SOM(), Space.SOM(), or Climate.SOM() if we only had those layers. The output is then stored to local variables for later analysis.

```
par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(som_model, type="changes", axes=F, 
     ylim=range(unlist(l_mat1[,1]))*c(0.9,1.1), main=NA, col="white")
for(i in 1:n){lines(l_mat1[,i],col=alpha("black",0.1))}

par(new = TRUE)
plot(l_mat2[,1], type="l", col="white",axes=F, 
     ylim=range(unlist(l_mat2[,1]))*c(0.9,1.1), main=NA,xlab=NA,ylab=NA,bty="n")
for(i in 1:n){lines(l_mat2[,i],col=alpha("red",0.1))}

par(new = TRUE)
plot(l_mat3[,1], type="l", col="white",axes=F, 
     ylim=range(unlist(l_mat3[,1]))*c(0.9,1.1), main=NA,xlab=NA,ylab=NA,bty="n")
for(i in 1:n){lines(l_mat3[,i],col=alpha("green",0.1))}

par(new = TRUE)
plot(l_mat4[,1], type="l", col="white",axes=F,
     ylim=range(unlist(l_mat4[,1]))*c(0.9,1.1), main=NA,xlab=NA,ylab=NA,bty="n")
for(i in 1:n){lines(l_mat4[,i],col=alpha("blue",0.1))}

axis(1);title("Training progress (steps)",line=-1)
axis(2,at=round(range(unlist(l_mat1[,1])),3)*c(0.9,1.75),las=3)
```

![learning](https://github.com/rpyron/delim-SOM/assets/583099/9799612c-c364-4c23-9c52-7e513bfd30e7)

We can then plot our learning estimates across the runs. The shape of this learning curve (slow decline, then sudden plateau) is inherent in the way the algorithm learns; longer runs produce the same shape, rather than a longer plateau. The scale of each variable determines the location of its plateau; we expect each matrix to stabilize, but not necessarily to converge to the same relative distance to closest unit.

```
colMeans(d_mat)
layers <- rev(sort(sqrt(1/colMeans(d_mat))))
layer.cols <- setNames(c("black","red","green","blue"),
                       c("alleles","space","climate","traits"))
barplot(layers,main="Layer Weights",col=layer.cols[names(layers)])
```

![weights](https://github.com/rpyron/delim-SOM/assets/583099/a0297d45-72e0-436b-824f-b233d15570ca)

Next, we can see the layer weights. Unsurprisingly, alleles dwarf everything else, but traits are slightly more important than climate, and both are greater than space alone.

```
par(mfrow=c(3,1),mar=c(0.5,4,1,0.5))

#by absolute BIC
boxplot(t(w_mat),outline=F,notch=T,axes=F, ylab="WSS",ylim=range(unlist(w_mat)))
axis(1,at=1:(k.max+1),labels=NA);axis(2,at=round(range(unlist(w_mat))),las=3);title("Number of clusters (k)", line=0)

#by delta BIC
d_wss <- apply(w_mat,2,function(x){diff(diff(x))});rownames(d_wss)<-2:k.max;plot_dwss <- rbind(NA,d_wss,NA)
boxplot(t(plot_dwss),outline=F,notch=T,axes=F, ylab="dWSS")
axis(1,at=1:(k.max+1),labels=NA);axis(2,at=sort(c(0,round(range(unlist(d_wss))))),las=3)

#by frequency
par(mar=c(4,4,1,0.5))
all_k <- apply(c_mat,2,max) # Get the K for each run
table(all_k)#How many different Ks were learned?
max.K <- max(c_mat)#What is the highest K?
cols <- viridis(max.K)#Set colors
barplot(table(factor(all_k,levels=1:k.max)),ylab="Posterior Samples")
```

![clusters](https://github.com/rpyron/delim-SOM/assets/583099/1cb3bda5-4150-47ab-ae82-1b6abdf51f27)

Then, we can see the optimal values of _K_.

```
q_mat <- match.k()#get admixture coefficients

par(mfrow=c(1,1),
    mar=c(0,0,0,0))
xy <- xyz[,1:2]
maps::map(database = 'county', xlim = range(xy[,1]) + c(-0.5,0.5), ylim = range(xy[,2]) + c(-0.5,0.5), col="white")
map.axes()
maps::map(database = 'county', xlim = range(xy[,1]) + c(-0.5,0.5), ylim = range(xy[,2]) + c(-0.5,0.5), col="gray",add=T)
maps::map(database = 'state', xlim = range(xy[,1]) + c(-0.5,0.5), ylim = range(xy[,2]) + c(-0.5,0.5), add = T)
make.admix.pie.plot(q_mat,xy,layer.colors = cols,radii=2.5,add = T)
legend(-88,38,legend=c(expression(italic("D. cheaha")),
         expression(italic("D. monticola"))),cex=2,pt.bg=rev(viridis(2)),pch=21)
map.scale(-81.2,31.1)
```

![map](https://github.com/rpyron/delim-SOM/assets/583099/a7386dcd-6435-4442-96ab-88f0c7abfd3d)

We can also produce a basic sample map. The match.k() function uses a CLUMPP-like algorithm to synchronize the cluster labels to the DAPC results from earlier.

```
x <- q_mat[order(q_mat[,1]),]
n <- hclust(dist(x),"single")$order
make.structure.plot(admix.proportions = x[n,], sample.names = rownames(x[n,]), mar = c(8,4,2,2), layer.colors = cols, sort.by = 1)
```

![struc](https://github.com/rpyron/delim-SOM/assets/583099/32a53ff0-5525-4ba2-a2f7-0977ff41e272)

A STRUCTURE-type plot, organized hierarchically by dominant cluster membership. Given the extensive differences between these two species in terms of genetics, geography, ecology, and phenotype, the species coefficients are sharply bimodal comapred to the individual ancestry coefficients estimated from the SNP matrix alone (see Pyron et al. 2023).

```
par(mfrow=c(2,1))

#Cell distances #
plot(som_model, type="dist.neighbours", main = "", palette.name = magma)
title("SOM neighbour distances",line=-0.1)

#SOM Clusters#
som.cols <- setNames(cols[labels[,max.K]],cluster_assignment)#Get colors to match original SOM clusters
som.cols <- unique(som.cols[sort(names(som.cols))])#Set to refactored labels

#plot cluster
plot(som_model, shape="straight", type="mapping", bgcol = som.cols[som_cluster], main = "", pch=19, col="red")
add.cluster.boundaries(som_model, som_cluster,col="red");title("SOM clusters",line=-0.1)
```

![model](https://github.com/rpyron/delim-SOM/assets/583099/58d600a3-d167-4533-874f-4eba055f4d35)

An example SOM grid from one learned output model.

```
#map
par(mar=c(0,0,0,0),mgp=c(3,0.5,0))
layout(mat = matrix(c(1, 2, 1, 3), 
                    nrow = 2, 
                    ncol = 2))
maps::map(database = 'world', xlim = range(xy[,1]) + c(-7,7), ylim = range(xy[,2]) + c(-0.5,0.5), col="white")
map.axes()
maps::map(database = 'county', xlim = range(xy[,1]) + c(-7,7), ylim = range(xy[,2]) + c(-0.5,0.5), col="gray",add=T)
maps::map(database = 'state', xlim = range(xy[,1]) + c(-7,7), ylim = range(xy[,2]) + c(-0.5,0.5), add = T)
maps::map(database = 'world', xlim = range(xy[,1]) + c(-7,7), ylim = range(xy[,2]) + c(-0.5,0.5), add = T)
make.admix.pie.plot(q_mat,xy,layer.colors = cols,radii=2,add = T)
legend(-92.5,38,legend=c(expression(italic("D. cheaha")),
                       expression(italic("D. monticola"))),
       cex=2,pt.bg=rev(viridis(2)),pch=21);map.scale(-80,32)

#K by frequency
par(mar=c(4,4,1,0.5))
all_k <- apply(c_mat,2,max) # Get the K for each run
table(all_k)#How many different Ks were learned?
max.K <- max(c_mat)#What is the highest K?
cols <- viridis(max.K)#Set colors
barplot(table(factor(all_k,levels=1:k.max)),ylab="Posterior Samples (%)",col="black",cex.names=1,cex.axis=1,main=expression(paste(bold("Clusters ("),bolditalic("K"),")")))

#weights
colMeans(d_mat)
layers <- rev(sort(sqrt(1/colMeans(d_mat))))
layer.cols <- setNames(c("black","red","green","blue"),
                       c("alleles","space","climate","traits"))
barplot(layers,main="Layer Weights",col=layer.cols[names(layers)],ylab="Relative Weight - sqrt(1/w)")
```

![pub](https://github.com/rpyron/delim-SOM/assets/583099/3c6b70f0-a055-41e5-9eb7-0de652f3261a)

A nice summary figure for publication!

```
save.image(file="Trait_SuperSOM.RData")
```

Finally, we can save all of our results.

# Hyperparameters

Many default implementations of SOMs have paid little attention to the hyperparameters of neighborhood type, learning rate, and run length (number of steps), finding them to have relatively small impacts on learning outcomes (Wehrens and Kruisselbrink 2018; Pyron et al. 2023). Given the importance of accurately quantifying layer contributions, we expand on this concern here. Previous studies have found that Gaussian neighborhoods (rather than bubble or heuristic) and linear learning-rates typically yield optimal results under a wide range of conditions (Stefanovič and Kurasova 2011; Natita et al. 2016). Consequently, we employ these as our default conditions. 

![hypers](https://github.com/rpyron/delim-SOM/assets/583099/0b286a66-07d1-4517-94d0-5182f14ec190)

The kohonen_hyper.R file contains a brief exploration of the hyperparameters, including the learning rates and run length. Generally speaking, these don't have much of an impact for rlen > 100 and alpha > 0.1 for the alleles matrix in the _D. monticola_ dataset. Longer runs (rlen/m) may have a small impact on final learning estimates and precision. These curves are estimated from only 25 estimates and may normalize after additional replicates.

# Possible Improvements

I have not explored different data types for molecular datasets (e.g., allele counts, structural variants). I also haven't explored the impact of normalization, although this should be minimal. Selection of _K_ for each run might deserve a bit more scrutiny, including the use of AIC/BIC instead of the raw Weighted Sum of Squares (WSS) values. Some possibilities are discussed here: https://bgstieber.github.io/post/an-introduction-to-the-kmeans-algorithm/. As currently implemented, it can theoretically support _K_=1 (see Jaynes et al. 2017), but this should be explored in depth for this and other methods (e.g., Derkarabetian et al. 2019).


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

