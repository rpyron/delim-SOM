##########################################
# SOM/UML SPECIES DELIMITATION - kohonen #
##########################################

#Load libraries
library(adegenet);library(maps);library(scales)
library(conStruct);library(poppr);library(kohonen)
library(lsr);library(combinat);library(GroupStruct)
library(viridis);set.seed(1)

#auxiliary function to normalize input features
#other normalizations are possible, but we will
#stick with minmax for the time being
minmax <- function(x){(x-min(x))/(max(x)-min(x))}

#colors
k.cols <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
layer.cols <- setNames(c("#603e95","#009da1","#fac22b","#d7255d"),c("alleles","space","climate","traits"))
#############################################################
#DNA only: #Single-layer SOM based on the allele frequencies#
#############################################################

DNA.SOM <- function()
  {
  ##Output matrices
  c_mat <- data.frame(row.names = row.names(alleles))#to hold classifiers for DNAsom
  l_mat <- data.frame(row.names = 1:m)#to hold learning values
  w_mat <- data.frame(row.names = 1:(11))#to hold dWSS for K
  
#Train the SOMs
for (j in 1:n)
{som_model <- som(alleles,
                  grid=som_grid,
                  maxNA.fraction=0.9,
                  alpha=c(0.5,0.1),
                  rlen=m)
l_mat[,j] <- som_model$changes

#calculate wss to select k
mydata <- getCodes(som_model)#pull codebook vectors for wss
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))#Get weighted sum of squares
for (i in 2:(11)) {wss[i] <- sum(kmeans(mydata,centers=i,nstart=25,iter.max=1e5)$withinss)}#Repeat kmeans to get best classification
N <- dim(mydata)[1]#sample size for bic
w_mat[,j] <- N*log(wss/N)+log(N)*(1:(11))#store bic

if(which(w_mat[,j]==min(w_mat[,j]))==1){num_clusters <- 1}else
{#Choose K by diffNgroup from DAPC
temp <- cutree(hclust(dist(diff(wss)), method="ward.D"), k=2)
goodgrp <- which.min(tapply(diff(wss), temp, mean))
num_clusters <- max(which(temp==goodgrp))+1}

# use hierarchical clustering to cluster the codebook vectors
som_cluster <- cutree(hclust(dist(mydata)), num_clusters)

# get vector with cluster value for each original data sample
cluster_assignment <- som_cluster[som_model$unit.classif]

# add assignment to classifier matrix
c_mat[,j] <- cluster_assignment

print(j)}
  
  #collect results
  return(list(c_mat=c_mat,
              l_mat=l_mat,
              w_mat=w_mat,
              som_model=som_model,
              som_cluster=som_cluster,
              cluster_assignment=cluster_assignment))
  }


#################################################
#Space: #Two-layer SuperSOM based on DNA + x,y,z#
#################################################

Space.SOM <- function()
{
  ##Output matrices
  c_mat <- data.frame(row.names = row.names(alleles))#to hold classifiers for DNAsom
  d_mat <- data.frame(alleles=NA, space=NA)#to hold distance weights
  l_mat1 <- data.frame(row.names = 1:m)#to hold learning values
  l_mat2 <- data.frame(row.names = 1:m)#to hold learning values
  w_mat <- data.frame(row.names = 1:(11))#to hold dWSS for K
  
  #Train the SOMs
  for (j in 1:n)
  {som_model <- supersom(data=list(alleles,
                                   space),
                         grid=som_grid,
                         maxNA.fraction=0.9,
                         alpha=c(0.5,0.1),
                         rlen=m)
  l_mat1[,j] <- som_model$changes[,1]; l_mat2[,j] <- som_model$changes[,2]#store learning
  d_mat[j,] <- som_model$distance.weights#store weights
  
  #calculate wss/bic to select k
  mydata <- cbind(getCodes(som_model)[[1]],
                  getCodes(som_model)[[2]])#pull codebook vectors for wss
  wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))#Get weighted sum of squares
  for (i in 2:(11)) {wss[i] <- sum(kmeans(mydata,centers=i,nstart=25,iter.max=1e5)$withinss)}#Repeat kmeans to get best classification
  N <- dim(mydata)[1]#sample size for bic
  w_mat[,j] <- N*log(wss/N)+log(N)*(1:(11))#store bic
  
  if(which(w_mat[,j]==min(w_mat[,j]))==1){num_clusters <- 1}else
  {#Choose K by diffNgroup from DAPC
    temp <- cutree(hclust(dist(diff(wss)), method="ward.D"), k=2)
    goodgrp <- which.min(tapply(diff(wss), temp, mean))
    num_clusters <- max(which(temp==goodgrp))+1}
  
  # use hierarchical clustering to cluster the codebook vectors
  som_cluster <- cutree(hclust(dist(mydata)), num_clusters)
  
  # get vector with cluster value for each original data sample
  cluster_assignment <- som_cluster[som_model$unit.classif]
  
  # add assignment to classifier matrix
  c_mat[,j] <- cluster_assignment
  
  print(j)}
  
  #collect results
  return(list(c_mat=c_mat,
              d_mat=d_mat,
              l_mat1=l_mat1,
              l_mat2=l_mat2,
              w_mat=w_mat,
              som_model=som_model,
              som_cluster=som_cluster,
              cluster_assignment=cluster_assignment))
}


###########################################################
#Climate: #Three-layer SuperSOM based on DNA + x,y,z + env#
###########################################################

Climate.SOM <- function()
{
  ##Output matrices
  c_mat <- data.frame(row.names = row.names(alleles))#to hold classifiers for DNAsom
  d_mat <- data.frame(alleles=NA, space=NA, climate=NA)#to hold distance weights
  l_mat1 <- data.frame(row.names = 1:m)#to hold learning values
  l_mat2 <- data.frame(row.names = 1:m)#to hold learning values
  l_mat3 <- data.frame(row.names = 1:m)#to hold learning values
  w_mat <- data.frame(row.names = 1:(11))#to hold dWSS for K
  
  #Train the SOMs
  for (j in 1:n)
  {som_model <- supersom(data=list(alleles,
                                   space,
                                   climate),
                         grid=som_grid,
                         maxNA.fraction=0.9,
                         alpha=c(0.5,0.1),
                         rlen=m)
  l_mat1[,j] <- som_model$changes[,1]; l_mat2[,j] <- som_model$changes[,2]; l_mat3[,j] <- som_model$changes[,3]#store learning
  d_mat[j,] <- som_model$distance.weights#store weights
  
  #calculate wss to select k
  mydata <- cbind(getCodes(som_model)[[1]],
                  getCodes(som_model)[[2]],
                  getCodes(som_model)[[3]])#pull codebook vectors for wss
  wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))#Get weighted sum of squares
  for (i in 2:(11)) {wss[i] <- sum(kmeans(mydata,centers=i,nstart=25,iter.max=1e5)$withinss)}#Repeat kmeans to get best classification
  N <- dim(mydata)[1]#sample size for bic
  w_mat[,j] <- N*log(wss/N)+log(N)*(1:(11))#store bic
  
  if(which(w_mat[,j]==min(w_mat[,j]))==1){num_clusters <- 1}else
  {#Choose K by diffNgroup from DAPC
    temp <- cutree(hclust(dist(diff(wss)), method="ward.D"), k=2)
    goodgrp <- which.min(tapply(diff(wss), temp, mean))
    num_clusters <- max(which(temp==goodgrp))+1}
  
  # use hierarchical clustering to cluster the codebook vectors
  som_cluster <- cutree(hclust(dist(mydata)), num_clusters)
  
  # get vector with cluster value for each original data sample
  cluster_assignment <- som_cluster[som_model$unit.classif]
  
  # add assignment to classifier matrix
  c_mat[,j] <- cluster_assignment
  
  print(j)}
  
  #collect results
  return(list(c_mat=c_mat,
              d_mat=d_mat,
              l_mat1=l_mat1,
              l_mat2=l_mat2,
              l_mat3=l_mat3,
              w_mat=w_mat,
              som_model=som_model,
              som_cluster=som_cluster,
              cluster_assignment=cluster_assignment))
}


################################################################
#Trait: #Four-layer SuperSOM based on DNA + x,y,z + env + pheno#
################################################################

Trait.SOM <- function()
{
  ##Output matrices
  c_mat <- data.frame(row.names = row.names(alleles))#to hold classifiers for DNAsom
  d_mat <- data.frame(alleles=NA, space=NA, climate=NA, traits=NA)#to hold distance weights
  l_mat1 <- data.frame(row.names = 1:m)#to hold learning values
  l_mat2 <- data.frame(row.names = 1:m)#to hold learning values
  l_mat3 <- data.frame(row.names = 1:m)#to hold learning values
  l_mat4 <- data.frame(row.names = 1:m)#to hold learning values
  w_mat <- data.frame(row.names = 1:(11))#to hold dWSS for K
  
  #Train the SOMs
  for (j in 1:n)
  {som_model <- supersom(data=list(alleles,
                                   space,
                                   climate,
                                   traits),
                         grid=som_grid,
                         maxNA.fraction=0.9,
                         alpha=c(0.5,0.1),
                         rlen=m)
  l_mat1[,j] <- som_model$changes[,1]; l_mat2[,j] <- som_model$changes[,2]; l_mat3[,j] <- som_model$changes[,3]; l_mat4[,j] <- som_model$changes[,4]#store learning
  d_mat[j,] <- som_model$distance.weights#store weights
  
  #calculate wss to select k
  mydata <- cbind(getCodes(som_model)[[1]],
                  getCodes(som_model)[[2]],
                  getCodes(som_model)[[3]],
                  getCodes(som_model)[[4]])#pull codebook vectors for wss
  wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))#Get weighted sum of squares
  for (i in 2:(11)) {wss[i] <- sum(kmeans(mydata,centers=i,nstart=25,iter.max=1e5)$withinss)}#Repeat kmeans to get best classification
  N <- dim(mydata)[1]#sample size for bic
  w_mat[,j] <- N*log(wss/N)+log(N)*(1:(11))#store bic
  
  if(which(w_mat[,j]==min(w_mat[,j]))==1){num_clusters <- 1}else
  {#Choose K by diffNgroup from DAPC
    temp <- cutree(hclust(dist(diff(wss)), method="ward.D"), k=2)
    goodgrp <- which.min(tapply(diff(wss), temp, mean))
    num_clusters <- max(which(temp==goodgrp))+1}
  
  # use hierarchical clustering to cluster the codebook vectors
  som_cluster <- cutree(hclust(dist(mydata)), num_clusters)
  
  # get vector with cluster value for each original data sample
  cluster_assignment <- som_cluster[som_model$unit.classif]
  
  # add assignment to classifier matrix
  c_mat[,j] <- cluster_assignment
  
  print(j)}
  
  #collect results
  return(list(c_mat=c_mat,
              d_mat=d_mat,
              l_mat1=l_mat1,
              l_mat2=l_mat2,
              l_mat3=l_mat3,
              l_mat4=l_mat4,
              w_mat=w_mat,
              som_model=som_model,
              som_cluster=som_cluster,
              cluster_assignment=cluster_assignment))
}

#############################
#Summarize across K for Qmat#
#Synchronize cluster labels #
#############################
#This implements something like the
#non-shortcut version of the CLUMPP
#algorithm from Jakobsson & Rosenberg (2007)
#but scaled to the DAPC labels as a reference
#thx to Gideon Bradburd for the tip

match.k <- function(res)
  {
#get labels for different K values
labels <- rbind.data.frame(lapply(1:10,function(x){find.clusters(a,n.clust=x,n.pca = dim(alleles)[1])$grp}))
rownames(labels) <- rownames(alleles);colnames(labels) <- paste("K",1:10,sep='')

#get clusters
c_mat <- res$c_mat
all_k <- apply(c_mat,2,max) # Get the K for each run
  
# Only save runs with knowable K (i.e., <10)
cc <- c_mat[,which(all_k<=10)]

#Create duplicate c_mat for translation
cca <- cc 

#CLUMPP-like labels
for(i in 1:n)
{
  run.k <- max(cca[,i])#What is the K of each run
  run.labels <- as.numeric(labels[,run.k])#pull the DAPC labels that match K
  refactor <- data.frame(row.names=rownames(labels))#Make a data frame to hold the possible relabelings
  label.perm <- permn(1:run.k)#How many label permutations are needed?
  for (j in 1:length(label.perm)){refactor[,j] <- cca[,i]}#Fill a data frame with K repeats for relabeling
  for(k in 1:length(label.perm)){refactor[,k] <- as.numeric(permuteLevels(factor(refactor[,k]),
                                                                          perm=label.perm[[k]]))}#Relabel each repeat with all possible permutations of cluster labels
  cc[,i] <- refactor[,which.min(apply(refactor,2,function(x){sum(abs(x-run.labels))}))]#Save closest relabeling to DAPC
}

#Summary stats for admixture
q_mat <- t(apply(cc,1,FUN=function(x){table(factor(unlist(x),levels=1:max(all_k)))})/dim(cca)[2])

#give matrix
q_mat
}

########################
###Plotting Functions###
########################
plotLearning.DNA <- function(res)
{
plot(res$som_model, type="changes", axes=T, 
     ylim=range(unlist(res$l_mat[,1]))*c(0.9,1.1), col="white")
for(i in 1:n){lines(res$l_mat[,i],col=alpha(layer.cols[1],0.1))}
}

plotLearning.Traits <- function(res)
{
  par(mar = c(5, 6, 4, 4) + 0.3)  # Leave space for z axis
  plot(res$som_model, type="changes", axes=F, 
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
  
  axis(1);title("Training progress (steps)",line=0)
  axis(2)
  #axis(2,at=round(range(unlist(res$l_mat1[,1])),3)*c(0.9,1.75),las=3)
}


plotK <- function(res)
{
  w_mat <- res$w_mat#Get BIC
  par(mfrow=c(3,1),mar=c(0.5,4,1,0.5))
  
  #by absolute BIC
  boxplot(t(w_mat)[,-11],outline=F,notch=T,axes=F, ylab="BIC",ylim=range(unlist(w_mat[,-11])),col=k.cols)
  axis(1,at=1:(10),labels=NA);axis(2,at=round(range(unlist(w_mat[,-11]))),las=3);title("Number of clusters (k)", line=0)
  
  #by delta BIC
  d_wss <- apply(w_mat,2,function(x){diff(diff(x))});rownames(d_wss)<-2:10;plot_dwss <- rbind(NA,d_wss)
  boxplot(t(plot_dwss),outline=F,notch=T,axes=F, ylab="delta BIC",col=k.cols);abline(h=0,lty=2,col="red")
  axis(1,at=1:(10),labels=NA);axis(2,at=sort(c(0,round(range(unlist(d_wss))))),las=3)
  
  #by frequency
  par(mar=c(4,4,1,0.5))
  all_k <- apply(res$c_mat,2,max) # Get the K for each run
  barplot(table(factor(all_k,levels=1:10))/n,ylab="Sampling Frequency",ylim=c(0,1),col=k.cols)
}


plotLayers <- function(res)
{
  print(colMeans(res$d_mat))
  layers <- rev(sort(sqrt(1/colMeans(res$d_mat))))
  barplot(layers,main="Layer Weights",col=layer.cols[names(layers)],ylab="Relative Weights - sqrt(1/w)")
}

