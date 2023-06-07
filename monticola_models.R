##########################################
# SOM/UML SPECIES DELIMITATION - kohonen #
##########################################
source("./kohonen_code.R")
set.seed(1)

#################################################
#Data from Pyron et al. 2023, Systematic Biology# 
#71 genetic  and 156 morphological specimens    #
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

##Dimensions
locmiss = propTyped(a, by = "loc")
barplot(sort(1-locmiss), ylim = c(0,1), ylab = "Complete genotypes (proportion)", xlab = "Locus", las = 2, cex.names = 0.7)
a = missingno(a, type = "loci", cutoff = 0.20);abline(h=0.2,col="red",lty=2)
a#trimmed to 20% missing data

#Convert genind object from adegenet to allele frequencies
struc <- makefreq(a)

#Convert allele frequences to matrix
alleles <- matrix(unlist(as.numeric(struc)), nrow=nrow(struc))

###SPATIAL & CLIMATIC DATA
#read in data
dat <- read.csv("./seal_clim.csv",row.names = "Specimen")
xyz <- dat[,c("LONG","LAT","elevation")]
space <- data.frame(lon=xyz$LONG,lat=xyz$LAT,elev=xyz$elev)
space <- apply(space,2,minmax)
climate <- apply(dat[,c(5:9)],2,minmax)#These variables were identified as most important

###PHENOTYPIC DATA
morph <- read.csv("./seal_morph.csv",row.names="Specimen")#Read in trait data, 156 specimens from 71 sites with 17 measurements 
morph_log <- data.frame(pop="seal",exp(morph[,2:18]))#un-log the measurements for GroupStruct
morph_allom <- allom(morph_log,"population2")#correct for allometry using GroupStruct
morph_mean <- aggregate(morph_allom[,-1],list(morph$pop),mean)#take mean by locality - this is a simplistic approach
traits <- apply(morph_mean[,-1],2,minmax)#could also use PC1-3 or similar transformation


#############################
###Baseline DAPC assignments#
#for synchronizing clusters #
#############################

#choose a max K
k.max <- 10

#get labels for different K values
labels <- data.frame(V1=rep(NA,dim(a$tab)[1]),row.names = rownames(a$tab))
for (i in 1:k.max){labels[,i] <- find.clusters(a,n.clust=i,n.pca = dim(a$tab)[1])$grp}


##################
###Kohonen maps###
##################

###Parameters for runs
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


##################
###Kohonen maps###
##################
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


###############
#Plot Learning#
###############
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


###############
#Layer Weights#
###############
colMeans(d_mat)
layers <- rev(sort(sqrt(1/colMeans(d_mat))))
layer.cols <- setNames(c("black","red","green","blue"),
                       c("alleles","space","climate","traits"))
barplot(layers,main="Layer Weights",col=layer.cols[names(layers)])


############
#Optimize K#
############
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


############
#Sample Map#
############
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
############


################
#Structure Plot#
################
x <- q_mat[order(q_mat[,1]),]
n <- hclust(dist(x),"single")$order
make.structure.plot(admix.proportions = x[n,], sample.names = rownames(x[n,]), mar = c(8,4,2,2), layer.colors = cols, sort.by = 1)
################


################################
#Example outputs from one model#
################################
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


################################
#Summary figure for publication#
################################
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

save.image(file="Trait_SuperSOM.RData")






