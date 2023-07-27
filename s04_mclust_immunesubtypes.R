###Cluster Immune Subtypes
##Need to find number of representative clusters and assign each cell to a cluster
library(clue)
library(ggplot2)
library(mclust)
library(factoextra)

#load in functions from Immune LandscapePaper-SubtypeCluster.R

#load in WGCNA data
library(WGCNA)
load("~/Desktop/PatelLab/Analysis_Results/WGCNA/wgcna_results.rda")
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

#use pred1 and predStrenght fxn to determine number of klus
determineK(MEdf)
#[1] 9 9 9 9 9 9 9 9 9 9
pred1(MEdf, 9, 0.25)
#[1] 0.6458924 0.6403326 1.0000000 0.7592593 0.4930966 0.6983471 0.9665552 0.7852761 0.7783133 with 9 clusters
predStrength(MEdf, 9, reps = 10)

#           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9}
#[1,] 0.5033408 0.5183486 0.5570470 0.6403101 0.7924528 0.9439252 0.9690722 0.9861111 0.9947090
#[2,] 0.3910615 0.4974359 0.5101010 0.5401146 0.6569579 0.8307292 0.9506934 0.9632000 0.9921875
#[3,] 0.6867470 0.7410526 0.7480000 0.8400000 0.8677918 0.9628040 0.9771574 0.9920949 0.9960784
#[4,] 0.3794466 0.7187500 0.8807157 0.8812589 0.9065850 0.9414894 0.9612188 0.9801136 1.0000000
#[5,] 0.3560976 0.4104938 0.5218391 0.7494033 0.8099352 0.8978102 0.9085366 0.9644809 0.9718876
#[6,] 0.4075472 0.4814815 0.5858586 0.6430868 0.8105263 0.8641618 0.8783383 0.9808307 0.9848024
#[7,] 0.4422111 0.5436242 0.7360595 0.7673410 0.8090278 0.9401042 0.9525862 0.9661836 1.0000000
#[8,] 0.3791209 0.5177305 0.5409836 0.6012658 0.6304348 0.7537313 0.8744650 0.9834437 0.9875000
#[9,] 0.5151515 0.5212121 0.7984645 0.8016878 0.8251534 0.8815029 0.9539683 0.9792453 0.9957082
#[10,] 0.5300699 0.6911765 0.7709163 0.7914692 0.8459215 0.8875562 0.8931818 0.9672727 0.9920000


#run model ensemble
model<-modelEnsemble(MEdf, klus = 9)

#use models to predict classification
modelclass<-ensemblePredict(model,MEdf, cores = 2)

#use cluster ensemble
consensus<-consensusEnsemble(model, MEdf)

preds <- ensemblePredict(model,MEdf,"list", 2)
partitions <- lapply(preds, function(a) as.cl_partition(a))
clpart <- cl_ensemble(list=partitions)

consensus <- cl_consensus(clpart, method="GV1")

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

#append majclass vector to MEdf
MEdf<-as.data.frame(MEdf)
MEdf$majclass<-majclass


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
