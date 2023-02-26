## Metabolite Annotation
install.packages("devtools")
devtools::install_github("jabiru/tictoc")
library("tictoc")

mdata <- read.csv("hmdb_complete.xlxs.csv")
mdata <- mdata[, c(1:2)]
data <- read_excel("Table_1_Integrated Metabolomics and Proteomics Analysis Reveals Plasma Lipid Metabolic Disturbance in Patients With Parkinsonâ__s Disease.xlsx")
data <- data[, c(2:4)]


#m/z to mass
for (x in 1:nrow(data)){
  ifelse(data[x,3]=='positive', data[x,1]<- (data[x,1]-1.008),NA)
  ifelse(data[x,3]=='negative', data[x,1]<-(data[x,1]+1.008 ),NA)
}

data["Metabolite"] <- ""

mdata <- mdata %>% na.omit()


##Checking the run-time of about 100 samples

tic("Run-time of 100 samples:")
for (m in 2:100){
  upper <- as.vector(data[m, 1] + 0.1)
  lower <- as.vector(data[m, 1] - 0.1)
  for (t in 1:nrow(mdata)){
    #print(as.vector(mdata[t, 2]))
    if (between(as.vector(mdata[t, 2]), lower, upper)){
      data[m,4] <- as.character(mdata[t,1])
    }
  }
}
toc()

mz_data <- as.vector(data[2:nrow(data), 1])
for (i in mz_data){
  print(class(i))
}



for (t in 1:nrow(mdata)){
  if (between(as.vector(mdata[t, 2]), lower, upper)){
    data[m,4] <- as.character(mdata[t,1])
  }
}





##code to search for all metabolites 

for (m in 2:nrow(data)){
  upper <- as.vector(data[m, 1] + 0.1)
  lower <- as.vector(data[m, 1] - 0.1)
  for (t in 1:nrow(mdata)){
    if (between(as.vector(mdata[t, 2]), lower, upper)){
      data[m,4] <- as.character(mdata[t,1])
    }
  }
}


