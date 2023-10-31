##########################################
# SOM/UML SPECIES DELIMITATION - kohonen #
##########################################
source("./kohonen_code.R")
set.seed(1)

#################################################
#Data from Pyron et al. 2023, Systematic Biology# 
#71 genetic  and 163 morphological specimens    #
#from 71 sites across the range of D. monticola #
#################################################

###MOLECULAR DATA
#Load the *.str file from PEA23
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
alleles <- makefreq(a)

###SPATIAL & CLIMATIC DATA
#read in data
dat <- read.csv("./seal_clim.csv",row.names = "Specimen")
xyz <- dat[,c("LONG","LAT","elevation")]
space <- data.frame(lon=xyz$LONG,lat=xyz$LAT,elev=xyz$elev)
space <- apply(space,2,minmax)
climate <- apply(dat[,c(5:9)],2,minmax)#These variables were identified as most important

###PHENOTYPIC DATA
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


##################
###Kohonen maps###
##################

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


##############
###Run SOMs###
##############
res <- Trait.SOM()

#Plot Learning#
plotLearning.Traits(res)

#Layer Weights#
plotLayers(res)

#Optimize K#
plotK(res)

#Sample Map#
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

#Structure Plot#
x <- q_mat[order(q_mat[,1]),]
z <- hclust(dist(x),"single")$order
make.structure.plot(admix.proportions = x[z,], 
                    sample.names = rownames(x[z,]), 
                    mar = c(8,4,2,2), 
                    layer.colors = k.cols, 
                    sort.by = 1)

#Example outputs from one model#
par(mfrow=c(2,1))

#Cell distances #
plot(res$som_model, type="dist.neighbours", main = "", palette.name = viridis)
title("SOM neighbour distances",line=-0.1)

#SOM Clusters#
som.cols <- setNames(k.cols[max.col(q_mat)],res$cluster_assignment)#Get colors to match original SOM clusters
som.cols <- unique(som.cols[sort(names(som.cols))])#Set to refactored labels

#plot cluster
plot(res$som_model, shape="straight", type="mapping", bgcol = som.cols[res$som_cluster], main = "", pch=19, col="red")
add.cluster.boundaries(res$som_model, res$som_cluster,col="red");title("SOM clusters",line=-0.1)


#################
#Compare to sNMF#
#################
res1 <- DNA.SOM()#Make a DNA SOM for direct comparisons of admixture
q_mat.DNA <- match.k(res1)#get admixture coefficients

sNMF_q_mat <- read.csv("seal_q.mat.csv",row.names=1)#Admixture estimates from sNMF

#Compare sNMF and DNA SOM admixture estimates
par(mfrow=c(2,1),mar=c(1,4.5,0.5,2),mgp=c(2,0.5,0))
plot(sNMF_q_mat[,1],q_mat.DNA[,1],pch=21,bg=layer.cols[1],xaxt='n',xlab=NA,
     ylab="DNA SOM",xlim=c(0,1),ylim=c(0,1),cex=2)
axis(side = 1, at = c(0,0.2,0.4,0.6,0.8,1), labels = FALSE)
b.sp <- cor(sNMF_q_mat[,1],q_mat.DNA[,1],method=c("spearman"))
b.pe <- cor(sNMF_q_mat[,1],q_mat.DNA[,1],method=c("pearson"))
text(0.2,0.9,paste0("Spearman's = ",round(b.sp,2)))
text(0.2,0.85,paste0("Pearson's = ",round(b.pe,2)))
text(0.675,0.05,"a) Individual Ancestries",font=2,cex=1.25)

#Compare SuperSOM species coefficients to sNMF admixture
par(mar=c(4,4.5,0.5,2))
plot(sNMF_q_mat[,1],q_mat[,1],pch=21,bg=layer.cols[4],
     xlab="sNMF Admixture Estimates",ylab="Trait SuperSOM",xlim=c(0,1),ylim=c(0,1),cex=2)
c.sp <- cor(sNMF_q_mat[,1],q_mat[,1],method=c("spearman"))
c.pe <- cor(sNMF_q_mat[,1],q_mat[,1],method=c("pearson"))
text(0.2,0.9,paste0("Spearman's = ",round(c.sp,2)))
text(0.2,0.84,paste0("Pearson's = ",round(c.pe,2)))
text(0.675,0.05,"b) Species Coefficients",font=2,cex=1.25)


#############
#Save output#
#############

save.image(file="Trait_SuperSOM.RData")

