##Develop Module Deconvolution Ratio Classifer
#use randomForest fxn
install.packages("randomForest")
library(randomForest)

##load in deconvolution data
Moduledeconv<-read.csv("~/Desktop/PatelLab/Analysis_Results/Deconvolution/Cibersort_Module_deconvolution.csv")

#load InhouseMeningioma_Clippr_results and bulk sample data
InhouseMeningioma_CLIPPR_results <- read.csv("InhouseMeningioma_CLIPPR_results.csv")
load("~/Desktop/PatelLab/09192022_GE_RAW_Normalized_N330PrimarySamples_GeneSymbol.rda")


#load in newpredsamples
newpred<-InhouseMeningioma_CLIPPR_results[,c(1,8)]
head(newpred)
intersecting<-intersect(samples$Sample.Name, newpred$ID)
index<-""
for(x in 1:length(intersecting)){
  index<-c(index,which(newpred$ID==intersecting[x]))
}
index<-index[-1]
index<-as.numeric(index)

newpred<-newpred[-index,]


index<-""
for(x in 1:length(intersecting)){
  index<-c(index,which(samples$Sample.Name==intersecting[x]))
}
index<-index[-1]
index<-as.numeric(index)

samples<-samples[-index,]

newpred<-newpred[c(-3,-4),]
newpred$samplename<-samples[1:2,1]
excludedsamples<-newpred
##two samples whose names needed manually fixing. Two samples did not overlap with the samples from the bulk data
InhouseMeningioma_CLIPPR_results[1,1]<-excludedsamples[1,3]
InhouseMeningioma_CLIPPR_results[2,1]<-excludedsamples[2,3]

##After above steps, InhouseMeningioma... should have all overlapping genes
##reload sample data
load("~/Desktop/PatelLab/09192022_GE_RAW_Normalized_N330PrimarySamples_GeneSymbol.rda")
#continue analysis
newpred<-InhouseMeningioma_CLIPPR_results[,c(1,2,8)]
samples$rowname<-row.names(samples)
newpredsamples<-merge(samples,newpred, by.x = "Sample.Name", by.y = "ID")
row.names(newpredsamples)<-newpredsamples$rowname

#index Moduledeconv with overlapping samples in newpredsamples
ID<-row.names(newpredsamples) 
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
