#Code for standard metabolomics data analysis

#Packages upload
# install.packages("BiocManager")
# install.packages("readxl")
# install.packages("readr")
# BiocManager::install('clusterProfiler')
# BiocManager::install('DOSE')
# BiocManager::install('org.Hs.eg.db')
#install.packages('factoextra')
library(readxl)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(tidyr)
library(edgeR)
library(tibble)
library(mixOmics)
library(tidyverse)  
library(cluster)    
library(factoextra) 
library(dendextend)

#File opening
data<- read_excel("Table_1_Integrated Metabolomics.xlsx")
data2 <- data[, c(1, 5:76)]

#Obtain a vector with the two groups 
classes <- data2[c(1), c(2:73)]
class(classes)
classes <- as.vector(classes)
classes <- unlist(classes)


#Obtain a factor with the two classes
classes_names <- classes
classes_names <- as.factor(classes_names)
Groups = as.factor(classes_names)

#Replace zeros with half of the minimum value of the intensities of all the compounds of each sample
data[data==0]= NA
for (x in 5:length(data)) {
  min_col = min(data[,x], na.rm= TRUE)/2
  data[,x][is.na(data[,x])] = min_col
}
#Data frame transposal
tran_data <- as_tibble(t(data[5:76]))
tran_data <- add_column(tran_data, "Groups" = classes)

#Define the two groups
matrix_tdata <- as.matrix(tran_data)
classes_names <- classes
loc <- 1
for (i in classes_names){
  ifelse(i==0, classes_names[loc] <- c("Control"), NA)
  ifelse(i==1, classes_names[loc] <- c("Parkinsons"), NA)
  loc <- loc + 1
}

#Statistical analysis

#Principal component analysis (PCA) analysis
my_pca <- spca(matrix_tdata, scale = TRUE)

#PCA plot
plotIndiv(my_pca, group = classes_names, 
          ellipse = TRUE, centroid = TRUE, 
          abline = FALSE, legend = TRUE, title = 'PCA on Metabolome Data')



### Partial Least Square Discirminate Analysis (PLS-DA)
matrix_tdata <- as.matrix(tran_data) #removing the "Group" column

MePLSDA.splsda <- splsda(matrix_tdata, classes_names) # 1 Run the method

#PLS-DA plot
plotIndiv(MePLSDA.splsda, ind.names = TRUE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'sPLS-DA on Metabolome Data',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')
plotVar(MePLSDA.splsda, var.names=T) 

#VIP scores from PLS-DA
vip_plsda=vip(MePLSDA.splsda )
selectVar(MePLSDA.splsda, comp=1)$value #selects the most relevant variables
plotLoadings(MePLSDA.splsda, ndisplay = 20, contrib = 'max', method = 'mean')

#Hierarchical clustering 

#Data scaling
scaled_tdata=scale(tran_data)

#Distance calculation
distance <- dist(tran_data, method = "euclidean")
hc1 <- hclust(distance, method = "ward.D" )
title_h= 'Hierarchical clustering'

#Hierarchical clustering plot
hclust=plot(hc1, cex = 0.6, hang = -1, labels= names(data[5:76]), main = title_h, col = "#0000FF")


#Compound annotation

#m/z to mass
for (x in 1:nrow(data)){
  ifelse(data[x,4]=='positive', data[x,2]<- (data[x,2]-1.008),NA)
  ifelse(data[x,4]=='negative', data[x,2]<-(data[x,2]+1.008 ),NA)
}

id <- read.csv("hmdb_complete.csv")
id <- drop_na(id,2)
data_b= id[,1:2]
experimental = data[2:nrow(data),2]
experimental["empty"] = ""
x= 0.01

#Filtering database rows according to experimental mass range
max_exp = max(experimental[,1])
min_exp = min(experimental[,1])
data_b=subset(data_b, data_b[,2]<= max_exp & data_b[,2]>= min_exp )


for (i in 1:nrow(experimental)){
  exp_mass = as.numeric(experimental[i,1])
  a= (data_b[,2]>= exp_mass-x &  data_b[,2]<= exp_mass+x)
  idx_true= which(a)
  hmdb_list = as.list(data_b[idx_true,1])
  if (length(hmdb_list) >= 1){
    experimental$empty[i] = list(hmdb_list)}
}



