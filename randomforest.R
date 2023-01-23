##Develop Module Deconvolution Ratio Classifer
#use randomForest fxn
install.packages("randomForest")
library(randomForest)

##load in deconvolution data. Cibersort_Module_deconvolution comes from Cibersort_Module_deconvolution.csv
Moduledeconv<-Cibersort_Module_deconvolution
ID<-row.names(newpredsamples)  ##newpredsamples can be found in deconvolution.R
index<-""
for(x in 1:length(ID)){
  index<-c(index, which(Moduledeconv[,1]==ID[x]))
}
index<-index[-1]
index<-as.numeric(index)
Moduledeconv<-Moduledeconv[index,]

row.names(Moduledeconv)<-Moduledeconv[,1]
Moduledeconv<-Moduledeconv[,-1]
Moduledeconv$rowname<-row.names((Moduledeconv))

#add new prediction annotations
sampleclass<-as.data.frame(newpredsamples$predictions_meta_class)
row.names(sampleclass)<-row.names(newpredsamples)
sampleclass$ID<-row.names(sampleclass)
colnames(sampleclass)<-c("TumorClass", "ID")

Moduledeconv<-merge(Moduledeconv,sampleclass, by.x = "rowname", by.y = "ID")
row.names(Moduledeconv)<-Moduledeconv[,1]
Moduledeconv<-Moduledeconv[,-1]
Moduledeconv$TumorClass<-as.factor(Moduledeconv$TumorClass)
#create train and data matrices
#train matrix
trainddf<-newpredsamples[newpredsamples$isTrain=="yes",]
trainindex<-""
ID<-row.names(trainddf)
for(x in 1:length(ID)){
  trainindex<-c(trainindex, which(row.names(Moduledeconv)==ID[x]))
}
trainindex<-trainindex[-1]
trainindex<-as.numeric(trainindex)
train<-Moduledeconv[trainindex,]
train<-as.data.frame(train)

train$TumorClass<-as.factor(train$TumorClass)

#example
library(dplyr)
data_train <- read.csv("https://raw.githubusercontent.com/guru99-edu/R-Programming/master/train.csv")


#data matrix
data<-Moduledeconv

#randomforest function
model <- randomForest(TumorClass~.,  data=train,method='class', ntree = 5000)
 ## all the samples will be used in prediction
class<-predict(model, newdata=data.frame(data), type='response')
prob<-predict(model, newdata=data.frame(data), type='prob')
#train is a matrix with rows module ids, columns training samples, values are ratios.
#data is a matrix with rows module id, columns all samples, values are ratios.

TumorClass<-data.frame(newpredsamples$predictions_meta_class, class)
colnames(TumorClass)<-c("Meta Predictions", "randomForest Class")

write.table(TumorClass, file = "~/Desktop/PatelLab/FinalFigureAnalysis/randomtreeTumorClasspredictions.txt", quote = F, sep = "\t")
write.table(prob, file = "~/Desktop/PatelLab/FinalFigureAnalysis/randomtreeprobability.txt", quote = F, sep = "\t")
