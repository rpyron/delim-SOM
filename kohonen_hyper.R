##########################################
# SOM/UML SPECIES DELIMITATION - kohonen #
##########################################

#Load libraries and set seed
library(adegenet); library(maps); library(viridis);library(scales)
library(LEA); library(conStruct);library(poppr); library(kohonen)
library(lsr);library(combinat)
set.seed(1)

#####################################################################
#Here, I test a DNA-only SOM (alleles) under various hyperparameters#
#####################################################################
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
a = missingno(a, type = "loci", cutoff = 0.20)
a#trimmed to 20% missing data

#Convert genind object from adegenet to allele frequencies
struc <- makefreq(a)

#Convert allele frequences to matrix
alleles <- matrix(unlist(as.numeric(struc)), nrow=nrow(struc))


#################
#Hyperparameters#
#################

#Create an output grid of size sqrt(n)
g <- round(sqrt(5*sqrt(length(rownames(a$tab)))))
som_grid <- somgrid(xdim = g,
                    ydim = g,
                    topo="hexagonal",
                    neighbourhood.fct = "gaussian")

#test rlen and alpha across a range of values
#this is a quick and dirty approach to show
#that the effects are small beyond extremes
#this sampling scheme could be expanded as needed

######################
#m=100,alpha(0.1,0.1)#
######################
m <- 100#number of steps
n <- 100#number of replicates
l_mat0.1_0.1 <- data.frame(row.names = 1:m)#to hold learning values

for (j in 1:n)
  {som_model <- som(alleles,grid=som_grid,maxNA.fraction=0.85,alpha=c(0.1,0.1),rlen=m)
  l_mat0.1_0.1[,j] <- som_model$changes
  print(j)}

######################
#m=100,alpha(0.5,0.1)#
######################
m <- 100#number of steps
n <- 100#number of replicates
l_mat0.5_0.1 <- data.frame(row.names = 1:m)#to hold learning values

for (j in 1:n)
  {som_model <- som(alleles,grid=som_grid,maxNA.fraction=0.85,alpha=c(0.5,0.1),rlen=m)
  l_mat0.5_0.1[,j] <- som_model$changes
  print(j)}

######################
#m=100,alpha(0.9,0.1)#
######################
m <- 100#number of steps
n <- 100#number of replicates
l_mat0.9_0.1 <- data.frame(row.names = 1:m)#to hold learning values

for (j in 1:n)
  {som_model <- som(alleles,grid=som_grid,maxNA.fraction=0.85,alpha=c(0.9,0.1),rlen=m)
  l_mat0.9_0.1[,j] <- som_model$changes
  print(j)}

#######################
#m=100,alpha(0.9,0.01)#
#######################
m <- 100#number of steps
n <- 100#number of replicates
l_mat0.9_0.01 <- data.frame(row.names = 1:m)#to hold learning values

for (j in 1:n)
  {som_model <- som(alleles,grid=som_grid,maxNA.fraction=0.85,alpha=c(0.9,0.01),rlen=m)
  l_mat0.9_0.01[,j] <- som_model$changes
  print(j)}

########################
#m=100,alpha(0.05,0.01)#
########################
m <- 100#number of steps
n <- 100#number of replicates
l_mat0.05_0.01 <- data.frame(row.names = 1:m)#to hold learning values

for (j in 1:n)
  {som_model <- som(alleles,grid=som_grid,maxNA.fraction=0.85,alpha=c(0.05,0.01),rlen=m)
  l_mat0.05_0.01[,j] <- som_model$changes
  print(j)}

######################
#m=100,alpha(0.5,0.5)#
######################
m <- 100#number of steps
n <- 100#number of replicates
l_mat0.5_0.5 <- data.frame(row.names = 1:m)#to hold learning values

for (j in 1:n)
{som_model <- som(alleles,grid=som_grid,maxNA.fraction=0.85,alpha=c(0.5,0.5),rlen=m)
  l_mat0.5_0.5[,j] <- som_model$changes
  print(j)}

######################
#m=200,alpha(0.5,0.1)#
######################
m <- 200#number of steps
n <- 100#number of replicates
l_mat0.5_0.1m200 <- data.frame(row.names = 1:m)#to hold learning values

for (j in 1:n)
  {som_model <- som(alleles,grid=som_grid,maxNA.fraction=0.85,alpha=c(0.5,0.1),rlen=m)
  l_mat0.5_0.1m200[,j] <- som_model$changes
  print(j)}

#m=1000,alpha(0.5,0.1)
##################
m <- 1000
n <- 100#number of replicates
#Output matrices
l_mat0.5_0.1m1000 <- data.frame(row.names = 1:m)#to hold learning values

for (j in 1:n)
{#train
  som_model <- som(alleles,grid=som_grid,maxNA.fraction=0.9,alpha=c(0.5,0.1),rlen=m)
  l_mat0.5_0.1m1000[,j] <- som_model$changes
  print(j)}

#################################################
#Plot various learning and length combinations  #
#to evaluate impacts on final learning estimates#
#################################################
pdf("../../Figure_3.pdf",8,8)
hyper.colors <- c("#ff6961","#ffb480","#f8f38d","#42d6a4","#08cad1","#59adf6","#9d94ff","#c780e8")
plot(density(unlist(l_mat0.1_0.1[100,])),xlim=range(unlist(l_mat0.1_0.1[100,]))*c(0.9,1.1),ylim=c(0,2500),
     col=hyper.colors[1],main="Final Learning Distances",lwd=2,xlab="Relative Distance to Closest Unit")
lines(density(unlist(l_mat0.5_0.1[100,])),col=hyper.colors[2],lwd=2)
lines(density(unlist(l_mat0.9_0.1[100,])),col=hyper.colors[3],lwd=2)
lines(density(unlist(l_mat0.9_0.01[100,])),col=hyper.colors[4],lwd=2)
lines(density(unlist(l_mat0.05_0.01[100,])),col=hyper.colors[5],lwd=2)
lines(density(unlist(l_mat0.5_0.5[100,])),col=hyper.colors[6],lwd=2)
lines(density(unlist(l_mat0.5_0.1m200[200,])),col=hyper.colors[7],lwd=2)
lines(density(unlist(l_mat0.5_0.1m1000[1000,])),col=hyper.colors[8],lwd=2)
legend(0.00675,2500,
  legend=c("alpha0.1/0.1m100",
           "alpha0.5/0.1m100",
           "alpha0.9/0.1m100",
           "alpha0.9/0.01m100",
           "alpha0.05/0.01m100",
           "alpha0.5/0.5m100",
           "alpha0.5/0.1m200",
           "alpha0.5/0.1m1000"),fill=hyper.colors)
#Seemingly very little impact for rlen>100 and alpha>0.1
dev.off()

save.image(file="SOM_Hyperparameters.RData")
