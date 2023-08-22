##Develop Module Deconvolution Ratio Classifer
#use randomForest fxn
library(randomForest)
library(caret)


##load in deconvolution data
Moduledeconv<-read.csv("~/Desktop/PatelLab/NewAnalysis/Deconvolution/Cibersortx_BCMsc_Results.csv")

Moduledeconv<-Moduledeconv[,1:9]

load("~/Desktop/PatelLab/031023_Inhouse_Raleigh_Hongkong_GeneSymbol.rda")

#add new prediction annotations
sampleclass<-as.data.frame(samplesheet$class)
row.names(sampleclass)<-row.names(samplesheet)
sampleclass$ID<-row.names(sampleclass)
colnames(sampleclass)<-c("TumorClass", "ID")

Moduledeconv<-merge(Moduledeconv,sampleclass, by.x = "Mixture", by.y = "ID")
row.names(Moduledeconv)<-Moduledeconv[,1]
Moduledeconv<-Moduledeconv[,-1]
Moduledeconv$TumorClass<-as.factor(Moduledeconv$TumorClass)
#create train and data matrices
#train matrix
set.seed(113)
trainIndex <- createDataPartition(Moduledeconv$TumorClass, p = 0.8, list = FALSE)
trainData <- Moduledeconv[trainIndex, ]
testData <- Moduledeconv[-trainIndex, ]
data<-Moduledeconv

#randomforest function
model <- randomForest(TumorClass~.,  data=trainData,method='class', ntree = 5000)
 ## all the samples will be used in prediction
class<-predict(model, newdata=testData, type='response')
prob<-predict(model, newdata=testData, type='prob')
#train is a matrix with rows module ids, columns training samples, values are ratios.
#data is a matrix with rows module id, columns all samples, values are ratios.

TumorClass<-data.frame(testData$TumorClass, class)
colnames(TumorClass)<-c("Meta Predictions", "randomForest Class")
table(TumorClass$`Meta Predictions`,TumorClass$`randomForest Class`)

write.table(TumorClass, file = "~/Desktop/PatelLab/NewAnalysis/randomtreeTumorClasspredictions.txt", quote = F, sep = "\t")
write.table(prob, file = "~/Desktop/PatelLab/NewAnalysis/randomtreeprobability.txt", quote = F, sep = "\t")
