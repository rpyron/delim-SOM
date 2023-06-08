##########################################
# SOM/UML SPECIES DELIMITATION - kohonen #
##########################################

#Load libraries
library(adegenet); library(maps); library(viridis);library(scales)
library(LEA); library(conStruct);library(poppr); library(kohonen)
library(lsr);library(combinat);library(caret);library(elevatr)
library(caret);library(GroupStruct);library(vcfR);library(dartR)
set.seed(1)

#auxiliary function to normalize input features
#other normalizations are possible, but we will
#stick with minmax for the time being
minmax <- function(x){(x-min(x))/(max(x)-min(x))}

#############################################################
#DNA only: #Single-layer SOM based on the allele frequencies#
#############################################################

DNA.SOM <- function()
  {
  ##Output matrices
  c_mat <- data.frame(row.names = row.names(a$tab))#to hold classifiers for DNAsom
  l_mat <- data.frame(row.names = 1:m)#to hold learning values
  w_mat <- data.frame(row.names = 1:(k.max+1))#to hold dWSS for K
  
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
for (i in 2:(k.max+1)) {wss[i] <- sum(kmeans(mydata,centers=i,nstart=25,iter.max=1e5)$withinss)}#Repeat kmeans to get best classification
w_mat[,j] <- wss#store wss

#Choose K by diffNgroup from DAPC
temp <- cutree(hclust(dist(diff(wss)), method="ward.D"), k=2)
goodgrp <- which.min(tapply(diff(wss), temp, mean))
num_clusters <- max(which(temp==goodgrp))+1

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
  c_mat <- data.frame(row.names = row.names(a$tab))#to hold classifiers for DNAsom
  d_mat <- data.frame(alleles=NA, space=NA)#to hold distance weights
  l_mat1 <- data.frame(row.names = 1:m)#to hold learning values
  l_mat2 <- data.frame(row.names = 1:m)#to hold learning values
  w_mat <- data.frame(row.names = 1:(k.max+1))#to hold dWSS for K
  
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
  for (i in 2:(k.max+1)) {wss[i] <- sum(kmeans(mydata,centers=i,nstart=25,iter.max=1e5)$withinss)}#Repeat kmeans to get best classification
  w_mat[,j] <- wss#store wss
  
  #Choose K by diffNgroup from DAPC
  temp <- cutree(hclust(dist(diff(wss)), method="ward.D"), k=2)
  goodgrp <- which.min(tapply(diff(wss), temp, mean))
  num_clusters <- max(which(temp==goodgrp))+1
  
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
  c_mat <- data.frame(row.names = row.names(a$tab))#to hold classifiers for DNAsom
  d_mat <- data.frame(alleles=NA, space=NA, climate=NA)#to hold distance weights
  l_mat1 <- data.frame(row.names = 1:m)#to hold learning values
  l_mat2 <- data.frame(row.names = 1:m)#to hold learning values
  l_mat3 <- data.frame(row.names = 1:m)#to hold learning values
  w_mat <- data.frame(row.names = 1:(k.max+1))#to hold dWSS for K
  
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
  for (i in 2:(k.max+1)) {wss[i] <- sum(kmeans(mydata,centers=i,nstart=25,iter.max=1e5)$withinss)}#Repeat kmeans to get best classification
  w_mat[,j] <- wss#store wss
  
  #Choose K by diffNgroup from DAPC
  temp <- cutree(hclust(dist(diff(wss)), method="ward.D"), k=2)
  goodgrp <- which.min(tapply(diff(wss), temp, mean))
  num_clusters <- max(which(temp==goodgrp))+1
  
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
  c_mat <- data.frame(row.names = row.names(a$tab))#to hold classifiers for DNAsom
  d_mat <- data.frame(alleles=NA, space=NA, climate=NA, traits=NA)#to hold distance weights
  l_mat1 <- data.frame(row.names = 1:m)#to hold learning values
  l_mat2 <- data.frame(row.names = 1:m)#to hold learning values
  l_mat3 <- data.frame(row.names = 1:m)#to hold learning values
  l_mat4 <- data.frame(row.names = 1:m)#to hold learning values
  w_mat <- data.frame(row.names = 1:(k.max+1))#to hold dWSS for K
  
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
  for (i in 2:(k.max+1)) {wss[i] <- sum(kmeans(mydata,centers=i,nstart=25,iter.max=1e5)$withinss)}#Repeat kmeans to get best classification
  w_mat[,j] <- wss#store wss
  
  #Choose K by diffNgroup from DAPC
  temp <- cutree(hclust(dist(diff(wss)), method="ward.D"), k=2)
  goodgrp <- which.min(tapply(diff(wss), temp, mean))
  num_clusters <- max(which(temp==goodgrp))+1
  
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

match.k <- function()
  {
# Only save runs with knowable K
cc <- c_mat[,which(all_k<=k.max)]

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
