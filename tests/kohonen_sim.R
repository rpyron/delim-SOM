###Basic simulations for kohonen species delimitation
source("../R/kohonen_code.R")
library(simulMGF)
set.seed(1)

###########################
#Simulate a K=1 SNP matrix#
###########################
simGeno(71, 5174)
alleles <- simG/2
rownames(alleles) <- paste("a",1:71,sep="")
colnames(alleles) <- paste("snp",1:5174,sep="")

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

#########################
###Parameters for runs###
#########################
#Size of Grid
g <- round(sqrt(5*sqrt(dim(alleles)[1])))#common rule of thumb

#Create an output grid of size sqrt(n)
som_grid <- somgrid(xdim = g,
                    ydim = g,
                    topo="hexagonal",
                    neighbourhood.fct = "gaussian")

n <- 100#Number of Replicates - can increase if you like
m <- 100#Number of steps - doesn't usually matter beyond ~100


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
q_mat <- match.k(res2,labels)#get admixture coefficients

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


###SIM FIGURE
# res, res1, res2 K & weights/map for res2
layout(mat = matrix(c(1,1,4,4,7,7,10,10,10,
                      1,1,4,4,7,7,10,10,10,
                      2,2,5,5,8,8,11,11,11,
                      2,2,5,5,8,8,11,11,11,
                      3,3,6,6,9,9,12,12,12,
                      3,3,6,6,9,9,12,12,12),ncol=9,byrow=T))
par(mar=c(2.5,4.5,2,1),mgp=c(2,0.5,0))

#####RESK
#by absolute WSS
w_mat <- res$w_mat#Get WSS
boxplot(t(w_mat)[,-11],outline=F,notch=T,axes=F, ylab="BIC",ylim=range(unlist(w_mat[,-11])),col=k.cols)
axis(1,at=1:(10),labels=NA);axis(2,at=round(range(unlist(w_mat[,-11]))),las=3)
text(7,195,"a) Null alleles",font=2,cex=1.5)

#by delta WSS
d_wss <- apply(w_mat,2,function(x){diff(diff(x))});rownames(d_wss)<-2:10;plot_dwss <- rbind(NA,d_wss)
boxplot(t(plot_dwss),outline=F,notch=T,axes=F, ylab="delta BIC",col=k.cols);abline(h=0,lty=2,col="red")
axis(1,at=1:(10),labels=NA);axis(2,at=sort(c(0,round(range(unlist(d_wss))))),las=3)
text(1,-0.45,"b)",font=2,cex=1.5,xpd=T)

#by frequency
all_k <- apply(res$c_mat,2,max) # Get the K for each run
barplot(table(factor(all_k,levels=1:10))/n,ylab="Sampling Frequency",ylim=c(0,1),col=k.cols)
text(3,0.9,"c)",font=2,cex=1.5)

#####RES1K
par(mar=c(2.5,4.5,2,1))
#by absolute WSS
w_mat <- res1$w_mat#Get WSS
boxplot(t(w_mat)[,-11],outline=F,notch=T,axes=F, ylab=NA,ylim=range(unlist(w_mat[,-11])),col=k.cols)
axis(1,at=1:(10),labels=NA);axis(2,at=round(range(unlist(w_mat[,-11]))),las=3);title("Number of clusters (k)", line=0)
text(7,189,"d) Null + 3 layers",font=2,cex=1.5)

#by delta WSS
d_wss <- apply(w_mat,2,function(x){diff(diff(x))});rownames(d_wss)<-2:10;plot_dwss <- rbind(NA,d_wss)
boxplot(t(plot_dwss),outline=F,notch=T,axes=F, ylab=NA,col=k.cols);abline(h=0,lty=2,col="red")
axis(1,at=1:(10),labels=NA);axis(2,at=sort(c(0,round(range(unlist(d_wss))))),las=3)
text(1,-0.4,"e)",font=2,cex=1.5,xpd=T)

#by frequency
all_k <- apply(res1$c_mat,2,max) # Get the K for each run
barplot(table(factor(all_k,levels=1:10))/n,ylab=NA,ylim=c(0,1),col=k.cols)
text(3,0.9,"f)",font=2,cex=1.5)

#####RES2K
par(mar=c(2.5,4.5,2,1))
#by absolute WSS
w_mat <- res2$w_mat#Get WSS
boxplot(t(w_mat)[,-11],outline=F,notch=T,axes=F, ylab=NA,ylim=range(unlist(w_mat[,-11])),col=k.cols)
axis(1,at=1:(10),labels=NA);axis(2,at=round(range(unlist(w_mat[,-11]))),las=3)
text(7,-6,"g) ~0 alleles\n + 3 layers",font=2,cex=1.5)

#by delta WSS
d_wss <- apply(w_mat,2,function(x){diff(diff(x))});rownames(d_wss)<-2:10;plot_dwss <- rbind(NA,d_wss)
boxplot(t(plot_dwss),outline=F,notch=T,axes=F, ylab=NA,col=k.cols);abline(h=0,lty=2,col="red")
axis(1,at=1:(10),labels=NA);axis(2,at=sort(c(0,round(range(unlist(d_wss))))),las=3)
text(1,-1,"h)",font=2,cex=1.5,xpd=T)

#by frequency
all_k <- apply(res2$c_mat,2,max) # Get the K for each run
barplot(table(factor(all_k,levels=1:10))/n,ylab=NA,ylim=c(0,1),col=k.cols)
text(3,0.9,"i)",font=2,cex=1.5)

par(mar=c(2,4.5,2,2))

#LAYERS
plotLayers(res2);text(4,1.1,"j)",font=2,cex=1.5)

par(mar=c(0,0,0,0))
#MAP
xy <- xyz[,1:2]
maps::map(database = 'county', xlim = range(xy[,1]) + c(-1.75,1.75), ylim = range(xy[,2]) + c(-0.5,0.5), col="white")
map.axes()
maps::map(database = 'county', xlim = range(xy[,1]) + c(-1.75,1.75), ylim = range(xy[,2]) + c(-0.5,0.5), col="gray",add=T)
maps::map(database = 'state', xlim = range(xy[,1]) + c(-1.75,1.75), ylim = range(xy[,2]) + c(-0.5,0.5), add = T)
make.admix.pie.plot(q_mat,xy,layer.colors = k.cols,radii=2,add = T);text(-80,31.5,"k)",font=2,cex=1.5)

#STRUCTURE
x <- q_mat[order(q_mat[,1]),]
z <- hclust(dist(x),"single")$order
make.structure.plot(admix.proportions = x[z,], 
                    sample.names = rownames(x[z,]), 
                    mar = c(2.5,4.5,2,2), 
                    layer.colors = k.cols, 
                    sort.by = 1)
text(-4,-0.125,"l)",font=2,cex=1.5,xpd=T)
