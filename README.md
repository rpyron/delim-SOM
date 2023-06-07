# delim-SOM
Using multi-layer Kohonen Self-Organizing Maps ("SuperSOMs") to delimit species and produce integrative taxonomies using Unsupervised Machine Learning (UML). This is accomplished primarily using the R package 'kohonen' (Wehrens and Buydens 2007): https://cran.r-project.org/web/packages/kohonen/index.html.

This repository expands the use of Kohonen maps as described in Pyron et al. (2023). It uses multi-layer Self-Organizing Maps ("SuperSOMs") in the R package 'kohonen' to delimit species based on allelic, spatial, climatic, and phenotypic data.

The contribution of each layer to the final model output is recorded, along with the clustering assignment of each individual over multiple learning replicates. The results therefore mirror a 'STRUCTURE'-type analysis including individual ancestry coefficients, but represent an unified delimitation model that incorporates various dimensions of ecological and evolutionary divergence for integative taxonomy. If only allelic data are used ('DNA.SOM()'), then the assignment probabilities approximate individual ancestry coefficients. If multiple layers are used, we treat them as "species coefficients," which might be useful for testing a variety of ecological and evolutionary hypotheses.

The requisite functions are in the 'kohonen_code.R' file, which loads the various dependencies:

```
library(adegenet); library(maps); library(viridis);library(scales)
library(LEA); library(conStruct);library(poppr); library(kohonen)
library(lsr);library(combinat);library(caret);library(elevatr)
library(caret);library(GroupStruct);library(vcfR);library(dartR)
```

Some of these may have to be installed manually or from various non-CRAN sources.

Overall, the method is extremely flexible and can take almost any data type or format, as long as it is introduced as a matrix in R. The matrices must be added in order, named struc_matrix, space, climate, and traits. Adding each of those matrices in sequence allows one to run SOMs based on DNA, DNA + xyz, DNA + xyz + environment, and DNA + xyz + environment + phenotypes. The primary requirement is to have individuals in rows in the same order in each matrix, and variables in columns, with no missing data and the same set of individuals in each matrix. I also min-max normalize the space, climate, and traits matrices to be on the same scale as the allele frequencies.

# Example: _Desmognathus monticola_, the Seal Salamander

![Pyron_et_al_Figure_3](https://github.com/rpyron/delim-SOM/assets/583099/9f28c7f0-0790-4a47-a7a4-25c98024a087)

Original SOM results from Pyron et al. (2023)

We provide a sample dataset and analysis for Seal Salamanders (_Desmognathus monticola_), which represent two species in the Appalachian Mountains and Piedmont of the eastern United States, based on four datasets comprising a SNP matrix from Genotype-By-Sequencing (GBS) analysis, lat/long/elevation (xyz), environment (climate, hydrology, and ecoregion), and phenotype (17 linear morphometric measurements).

The genetic, spatial, and environmental data come from 71 individuals from 71 sites, while the phenotypic data are expanded to include 156 specimens from those sites, with the mean of each measurement taken after size correction. The allele frequencies come from a GBS matrix of 5,174 SNPs and 10,526 alleles after trimming to 80% completeness.

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
dat <- read.csv("./seal_clim.csv",row.names = "Specimen")
xyz <- dat[,c("LONG","LAT","elevation")]
space <- data.frame(lon=xyz$LONG,lat=xyz$LAT,elev=xyz$elev)
space <- apply(space,2,minmax)
climate <- apply(dat[,c(5:9)],2,minmax)#These variables were identified as most important
```

Similarly, these are just the lat, long, elevation, and climate data from Pyron et al. (2023).

```
morph <- read.csv("./seal_morph.csv",row.names="Specimen")#Read in trait data, 156 specimens from 71 sites with 17 measurements 
morph_log <- data.frame(pop="seal",exp(morph[,2:18]))#un-log the measurements for GroupStruct
morph_allom <- allom(morph_log,"population2")#correct for allometry using GroupStruct
morph_mean <- aggregate(morph_allom[,-1],list(morph$pop),mean)#take mean by locality - this is a simplistic approach
traits <- apply(morph_mean[,-1],2,minmax)#could also use PC1-3 or similar transformation
```

Finally, the morphological dataset from Pyron et al. (2023) trimmed to the 156 specimens from the 71 genetic localities. The log-transformed data are read in, exponentiated for analysis, assigned the same "species" of "seal," corrected for allometry using 'GroupStruct,' taking the mean by site, and min-max normalizing the resulting matrix.

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

# Hyperparameters



# References

Chan, K.O. & Grismer, L. L. (2021). A standardized and statistically defensible framework for quantitative morphological analyses in taxonomic studies. Zootaxa, 5023: 293-300.

Chan, K.O. and Grismer, L.L., 2022. GroupStruct: an R package for allometric size correction. Zootaxa, 5124(4), pp.471-482.

Jakobsson, M. and Rosenberg, N.A., 2007. CLUMPP: a cluster matching and permutation program for dealing with label switching and multimodality in analysis of population structure. Bioinformatics, 23(14), pp.1801-1806.

Jombart, T., 2008. adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24(11), pp.1403-1405.

Pyron, R.A., O’Connell, K.A., Duncan, S.C., Burbrink, F.T. and Beamer, D.A., 2023. Speciation hypotheses from phylogeographic delimitation yield an integrative taxonomy for Seal Salamanders (Desmognathus monticola). Systematic Biology, 72(1), pp.179-197.

Wehrens, R. and Buydens, L.M., 2007. Self-and super-organizing maps in R: the Kohonen package. Journal of Statistical Software, 21, pp.1-19.
