###Cluster Immune Subtypes
##Need to find number of representative clusters and assign each cell to a cluster
library(clue)
library(GGally)
library(ggplot2)
library(mclust)
library(factoextra)
library(parallel)

#load in functions from Immune LandscapePaper-SubtypeCluster.R

#load in WGCNA data
library(WGCNA)
setwd("~/Desktop/PatelLab/NewAnalysis")
load("~/Desktop/PatelLab/NewAnalysis/WGCNA/wgcna_results.rda")
#visualizisation
moduleLabels = net$colors


moduleColors = labels2colors(net$colors)
names(moduleColors)<-names(moduleLabels)
table(moduleLabels)
table(moduleColors)

#create data frame with eigensignature values for every cell
MEdf<-net$MEs
MEdf<-scale(MEdf)
MEdf<-MEdf[,-6]

#load functions from ImmuneLandscapepaper
#use pred1 and predStrenght fxn to determine number of klus
determineK(MEdf)
#[1] 9 9 9 9 9 9 9 9 9 9
pred1(MEdf, 9, 0.25)
#[1] 0.9284362 0.6054422 0.9434997 0.9351536 0.7661932 0.9436620 0.9904679 0.4987605 0.5664467
predStrength(MEdf, 9, reps = 10)

#[,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9]
#[1,] 0.4485166 0.4780269 0.4892368 0.5161863 0.5443925 0.9481627 0.9485615 0.9540128 0.9644588
#[2,] 0.5720040 0.5955538 0.8360996 0.8475751 0.8619644 0.8645304 0.8730159 0.9006648 0.9292697
#[3,] 0.4736412 0.6027923 0.6061947 0.7626872 0.7729358 0.8200717 0.9060606 0.9113747 0.9129944
#[4,] 0.4245122 0.4911550 0.5389134 0.7669548 0.8962264 0.9010989 0.9150162 0.9582851 0.9777382
#[5,] 0.4322840 0.6071069 0.6241799 0.6503597 0.7976085 0.8892332 0.9372197 0.9724335 0.9821718
#[6,] 0.4582221 0.5576826 0.6663851 0.9130435 0.9248844 0.9356298 0.9523425 0.9661458 0.9806608
#[7,] 0.5368852 0.8273543 0.8595840 0.8655793 0.8740554 0.9090242 0.9350962 0.9502605 0.9838230
#[8,] 0.6912929 0.8467909 0.8714319 0.9370306 0.9414173 0.9434918 0.9663770 0.9837209 0.9900990
#[9,] 0.5815324 0.6270197 0.6534091 0.6661355 0.7780679 0.8467213 0.9016455 0.9328679 0.9582463
#[10,] 0.6534831 0.6672755 0.6935933 0.8861502 0.9376083 0.9521127 0.9638672 0.9705320 0.9809853


#run model ensemble
model<-modelEnsemble(MEdf, klus = 9)

#use models to predict classification
modelclass<-ensemblePredict(model,MEdf, cores = 2)

#use cluster ensemble
consensus<-consensusEnsemble(model, MEdf)

preds <- ensemblePredict(model,MEdf,"list", 2)
partitions <- lapply(preds, function(a) as.cl_partition(a))
clpart <- cl_ensemble(list=partitions)

consensus <- cl_consensus(clpart, method="GV1") #Do I need to change the method? It's what the immune paper used. 
consensus$.Data
classification<-data.frame(consensus$.Data[,1:9])
colnames(classification)<-1:9

#assign a classification to each cell using >0.5
majclass<-""
for(x in 1:nrow(classification)){
  a<-classification[x,]
  maxclass<-rownames(classification)[which(a == max(a))]
  maxclass<-maxclass[1]
  if(a[maxclass] > 0.5){
    majclass<-c(majclass, maxclass)
  } else {
    majclass<-c(majclass, NA)
  }
}
majclass<-majclass[-1]


##can load in data here

subtypecluster <- read.csv("~/Desktop/PatelLab/NewAnalysis/subtypecluster.csv", row.names=1)
majclass<-subtypecluster$x
majclass<-as.factor(majclass)

#append majclass vector to MEdf
MEdf<-as.data.frame(MEdf)
MEdf$majclass<-majclass

#names(majclass)<-row.names(MEdf)
#majclass
#write.csv(majclass, file = "~/Desktop/PatelLab/NewAnalysis/subtypecluster.csv")



#change module name
###classifications: MEgreen:C5, MEbrown:C3, MEturquoise:C1, C2:MEblue, C4:MEyellow, C0:MEgrey
newnames<-c("C5", "C3", "C1", "C2", "C4", "majclass")
colnames(MEdf)<-newnames


#create bar plot with average expression scores for eigen modules by cluster
#need to create a new dataframe format in order to input into ggploot
for(x in 1:5){
  if(x == 1){
  bardf<-data.frame("normalized score" = MEdf[,x], "EigenModule" = rep(colnames(MEdf)[x], nrow(MEdf)), "subtypecluster"= MEdf[,6])
  } else {
    bardf<-rbind(bardf, data.frame("normalized score" = MEdf[,x], 
                                   "EigenModule" = rep(colnames(MEdf)[x], nrow(MEdf)), "subtypecluster"= MEdf[,6])
    )
    
  }
}

ggplot(data = bardf, aes(x = EigenModule, y = normalized.score, fill = subtypecluster))+
  geom_bar(stat = "summary",
           fun = "mean" , 
           position= position_dodge())


#create bar plot with standard deviation
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

sumdf<-data_summary(bardf, varname = "normalized.score", groupnames = c("EigenModule",
                                                                        "subtypecluster"))

ggplot(data = sumdf, aes(x = EigenModule, y = normalized.score, fill = subtypecluster))+
  geom_bar(stat = "identity",position= position_dodge())+
  geom_errorbar(aes(ymin=normalized.score-sd, ymax=normalized.score+sd), width=.2,
                position=position_dodge(.9))


#create bar plot with eigen expression data by tumor class
#load in maj class
subtypecluster <- read.csv("~/Desktop/PatelLab/NewAnalysis/subtypecluster.csv", row.names=1)

#load in seurat data
load("072523_allLocations_immunecellsonly_modulescores.rda")


#create vector with tumor class
tumorclass<-as.factor(seuratObj$TumorType)
bardf<-data.frame("cluster" = subtypecluster, "tumor class" = tumorclass)
colnames(bardf)<- c("cluster", 
                  "tumor.class")
bardf$cluster<-as.factor(bardf$cluster)

#create df with frequency


# Create a data frame with tumor class and cluster
tumorclass <- as.factor(seuratObj$TumorType)
bardf <- data.frame("cluster" = majclass, "tumor class" = tumorclass)
colnames(bardf) <- c("cluster", "tumor.class")
bardf$cluster <- as.factor(bardf$cluster)


relative_freq <- prop.table(table(bardf$tumor.class, bardf$cluster), margin = 1)

# Convert the relative frequencies into a data frame
relative_freq_df <- as.data.frame.matrix(relative_freq)
relative_freq_df$tumor.class <- rownames(relative_freq_df)

# Reshape the data for plotting
relative_freq_long <- tidyr::gather(relative_freq_df, key = "cluster", value = "relative_freq", -tumor.class)

ggplot(data = relative_freq_long, aes(x = tumor.class, y = relative_freq, fill = cluster)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Tumor Class", y = "Relative Frequency", fill = "Cluster") +
  theme_minimal()

#create seurat object with clusters
subtypecluster <- read.csv("~/Desktop/PatelLab/NewAnalysis/subtypecluster.csv", row.names=1)
load("072523_allLocations_immunecellsonly_modulescores.rda")
seuratObj$subtypecluster<-subtypecluster

sum(is.na(seuratObj$subtypecluster))

#remove NA cells
subclustcells<-subset(seuratObj, cells=colnames(seuratObj)[!is.na(seuratObj$subtypecluster)])
save(subclustcells,file = "~/Desktop/PatelLab/NewAnalysis/subtypecluster_cells.rda")
