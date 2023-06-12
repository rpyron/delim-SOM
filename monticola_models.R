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
alleles <- matrix(unlist(as.numeric(struc)), nrow=nrow(struc))
rownames(alleles) <- rownames(struc);colnames(alleles) <- colnames(struc)

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


#############################
###Baseline DAPC assignments#
#for synchronizing clusters #
#############################

#get labels for different K values
labels <- data.frame(V1=rep(NA,dim(alleles)[1]),row.names = rownames(alleles))
for (i in 1:10){labels[,i] <- find.clusters(a,n.clust=i,n.pca = dim(alleles)[1])$grp}


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
som.cols <- setNames(k.cols[labels[,max(res$c_mat)]],res$cluster_assignment)#Get colors to match original SOM clusters
som.cols <- unique(som.cols[sort(names(som.cols))])#Set to refactored labels

#plot cluster
plot(res$som_model, shape="straight", type="mapping", bgcol = som.cols[res$som_cluster], main = "", pch=19, col="red")
add.cluster.boundaries(res$som_model, res$som_cluster,col="red");title("SOM clusters",line=-0.1)


################################
#Summary figure for publication#
################################
png("../../Pyron_UML_Graphical_Abstract.png",1328,531)
#map
par(mar=c(0,0,0,0),mgp=c(3,0.75,0))
layout(mat = matrix(c(1,1,1,1,1,2,2,3,3,
                      1,1,1,1,1,2,2,3,3,
                      1,1,1,1,1,4,4,5,5,
                      1,1,1,1,1,4,4,5,5), 
                    ncol = 9,byrow=T))
maps::map(database = 'world', xlim = range(xy[,1]) + c(-4,4), ylim = range(xy[,2]) + c(-1,1), col="white")
map.axes(cex.axis=1.5);title(main="Self-Organizing Maps for Integrative Species Delimitation",line=1,cex.main=2.5)
maps::map(database = 'county', xlim = range(xy[,1]) + c(-4,4), ylim = range(xy[,2]) + c(-1,1), col="gray",add=T)
maps::map(database = 'state', xlim = range(xy[,1]) + c(-4,4), ylim = range(xy[,2]) + c(-1,1), add = T)
maps::map(database = 'world', xlim = range(xy[,1]) + c(-4,4), ylim = range(xy[,2]) + c(-1,1), add = T)
make.admix.pie.plot(q_mat,xy,layer.colors = k.cols,radii=2.5,add = T)
legend(-91.5,38,legend=c(expression(italic("D. cheaha")),
                       expression(italic("D. monticola"))),
       cex=3,pt.bg=k.cols[2:1],pch=21);map.scale(-80,31.5,cex=2)

#K by frequency
par(mar=c(4,4,1,0.5))
all_k <- apply(res$c_mat,2,max) # Get the K for each run
table(all_k)#How many different Ks were learned?
max.K <- max(res$c_mat)#What is the highest K?
barplot(table(factor(all_k,levels=1:10))/n,ylab="Sampling Proportion",ylim=c(0,1),
        col="black",cex.lab=1.5,cex.axis=1.5,cex.names=1.5)
title(main=expression(paste(bold("Clusters ("),bolditalic("K"),")")),
      cex.main=2.5,line=-0.25)

#weights
par(mar=c(4,4.5,1,0.5),mgp=c(2.5,0.75,0))
print(colMeans(res$d_mat))
layers <- rev(sort(sqrt(1/colMeans(res$d_mat))))
barplot(layers,col=layer.cols[names(layers)],ylab="Relative Weights - sqrt(1/w)",cex.lab=1.5,cex.names=1.5,cex.axis=1.5)
title(main="Layer Weights",cex.main=2.5,adj=0.65,line=-0.5)

#learning
par(mar = c(4, 5, 2, 2) + 0.3,mgp=c(2.5,0.75,0))  # Leave space for z axis
plot(res$som_model, type="changes", axes=F, cex.lab=1.5,
     ylim=range(unlist(res$l_mat1[,1]))*c(0.9,1.1), main=NA, col="white")
for(i in 1:n){lines(res$l_mat1[,i],col=alpha(layer.cols[1],0.1))}

par(new = TRUE)
plot(res$l_mat2[,1], type="l", col="white",axes=F, 
     ylim=range(unlist(res$l_mat2[,1]))*c(0.9,1.1), main=NA,xlab=NA,ylab=NA,bty="n")
for(i in 1:n){lines(res$l_mat2[,i],col=alpha(layer.cols[2],0.1))}

par(new = TRUE)
plot(res$l_mat3[,1], type="l", col="white",axes=F, 
     ylim=range(unlist(res$l_mat3[,1]))*c(0.9,1.1), main=NA,xlab=NA,ylab=NA,bty="n")
for(i in 1:n){lines(res$l_mat3[,i],col=alpha(layer.cols[3],0.1))}

par(new = TRUE)
plot(res$l_mat4[,1], type="l", col="white",axes=F,
     ylim=range(unlist(res$l_mat4[,1]))*c(0.9,1.1), main=NA,xlab=NA,ylab=NA,bty="n")
for(i in 1:n){lines(res$l_mat4[,i],col=alpha(layer.cols[4],0.1))}

axis(1,cex.axis=1.5,cex=1.5);title("Training progress",line=0,cex.main=2.5)
axis(2,cex.axis=1.5,cex=1.5)

#SOM
som.cols <- setNames(k.cols[labels[,max(res$c_mat)]],res$cluster_assignment)#Get colors to match original SOM clusters
som.cols <- unique(som.cols[sort(names(som.cols))])#Set to refactored labels

#plot cluster
plot(res$som_model, shape="straight", type="mapping", bgcol = som.cols[res$som_cluster], main = "", pch=19, col="red")
add.cluster.boundaries(res$som_model, res$som_cluster, col="red");title("SOM clusters",line=-0.1,cex.main=2.5)
dev.off()


#################
#Compare to sNMF#
#################
res1 <- DNA.SOM()

q_mat.DNA <- match.k(res1)#get admixture coefficients
sNMF_q_mat <- read.csv("seal_q.mat.csv",row.names=1)

pdf("../../Figure_admix.pdf",5,8)
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
dev.off()

#############
#Save output#
#############

save.image(file="Trait_SuperSOM.RData")


#############
#MPE Figures#
#############

#som grid, som neighbors, learning, weights 
pdf("../../Figure_monticola_SOM.pdf",10,10)
par(mfrow=c(2,2))

#Cell Distances#
plot(res$som_model, type="dist.neighbours", main = "", palette.name = viridis)
title("SOM neighbour distances",line=-0.1);text(0,0,"a)",font=2,cex=1.25)

#Learning
plotLearning.Traits(res);text(0,0.013,"b)",font=2,cex=1.25,xpd=T)

#SOM Clusters#
som.cols <- setNames(k.cols[labels[,max(res$c_mat)]],res$cluster_assignment)#Get colors to match original SOM clusters
som.cols <- unique(som.cols[sort(names(som.cols))])#Set to refactored labels
plot(res$som_model, shape="straight", type="mapping", bgcol = som.cols[res$som_cluster], main = "", pch=19, col="red")
add.cluster.boundaries(res$som_model, res$som_cluster,col="red");title("SOM clusters",line=1)
text(0.8,0.1,"c)",font=2,cex=1.25,xpd=T)

#Weights
plotLayers(res);text(0,-1.5,"d)",font=2,cex=1.25,xpd=T)
dev.off()


#k, map, struc
pdf("../../Figure_monticola_k_map.pdf",8,8)
layout(mat = matrix(c(1,1,4,4,4,
                      1,1,4,4,4,
                      2,2,4,4,4,
                      2,2,5,5,5,
                      3,3,5,5,5,
                      3,3,5,5,5),ncol=5,byrow=T))
par(mar=c(2.5,4.5,2,2))
#by absolute WSS
w_mat <- res$w_mat#Get WSS
boxplot(t(w_mat)[,-11],outline=F,notch=T,axes=F, ylab="BIC",ylim=range(unlist(w_mat[,-11])),col=k.cols)
axis(1,at=1:(10),labels=NA);axis(2,at=round(range(unlist(w_mat[,-11]))),las=3);title("Number of clusters (k)", line=0)
text(9,150,"a)",font=2,cex=1.5)

#by delta WSS
d_wss <- apply(w_mat,2,function(x){diff(diff(x))});rownames(d_wss)<-2:10;plot_dwss <- rbind(NA,d_wss)
boxplot(t(plot_dwss),outline=F,notch=T,axes=F, ylab="delta BIC",col=k.cols);abline(h=0,lty=2,col="red")
axis(1,at=1:(10),labels=NA);axis(2,at=sort(c(0,round(range(unlist(d_wss))))),las=3)
text(9,25,"b)",font=2,cex=1.5,xpd=T)

#by frequency
all_k <- apply(res$c_mat,2,max) # Get the K for each run
barplot(table(factor(all_k,levels=1:10))/n,ylab="Sampling Frequency",ylim=c(0,1),col=k.cols)
text(10,0.9,"c)",font=2,cex=1.5)

#map
xy <- xyz[,1:2]
maps::map(database = 'county', xlim = range(xy[,1]) + c(-1.5,1.5), ylim = range(xy[,2]) + c(-0.5,0.5), col="white")
map.axes()
maps::map(database = 'county', xlim = range(xy[,1]) + c(-1.5,1.5), ylim = range(xy[,2]) + c(-0.5,0.5), col="gray",add=T)
maps::map(database = 'state', xlim = range(xy[,1]) + c(-1.5,1.5), ylim = range(xy[,2]) + c(-0.5,0.5), add = T)
make.admix.pie.plot(q_mat,xy,layer.colors = k.cols,radii=2,add = T)
legend(-89,38,legend=c(expression(italic("D. cheaha")),
                       expression(italic("D. monticola"))),
       cex=1.5,pt.bg=k.cols[2:1],pch=21,bg="white")
map.scale(-81,31.6,relwidth = 0.15,ratio=FALSE);text(-80,32,"d)",font=2,cex=1.5)

#Structure Plot#
x <- q_mat[order(q_mat[,1]),]
z <- hclust(dist(x),"single")$order
make.structure.plot(admix.proportions = x[z,], 
                    sample.names = rownames(x[z,]), 
                    mar = c(8,4,2,2), 
                    layer.colors = k.cols, 
                    sort.by = 1)
text(-2.5,-0.15,"e)",font=2,cex=1.5,xpd=T)
dev.off()
