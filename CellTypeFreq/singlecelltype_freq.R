#Cell Type Frequency
#generate celltypist and singleR cell type frequency among different classes
#ignore the cells in our tumor cluster
##must switch cluster annoations first (Cluster_annotations.R)
#seurat obj be found in UT box

##load in seurat obj
load("~/Desktop/PatelLab/102022_SCELL_UCSF_BCM_WITH_newannotations_addmodulescore_seuratObj.rda")

sub<- subset(seuratObj, cells=colnames(seuratObj)[seuratObj$orig.ident %in% c("BCM_Front", "BCM_Post", "MSC3", "MSC5", "MSC1", "MSC4")])


#My predicted annotations of each cluster is:


new.cluster.ids <- c("other", #0
                     
                     "tumor", #1
                     
                     "microglia/monocyte", #2
                     
                     "tumor", #3
                     
                     "T-cell", #4
                     
                     "tumor",#5
                     
                     "EC", #6
                     
                     "other",#7
                     
                     "other",#8
                     
                     "microglia/monocyte",#9
                     
                     "tumor",#10
                     
                     "other",#11
                     
                     "tumor",#12
                     
                     "T-cell", #13
                     
                     "EC",#14
                     
                     "tumor",#15
                     
                     "other" ,#16
                     
                     "other",#17
                     
                     "other",#18
                     
                     "microglia/monocyte",#19
                     
                     "microglia/monocyte", #20
                     
                     "other", #21
                     
                     "T-cell", #22
                     
                     "microglia/monocyte", #23
                     
                     "other") #24

names(new.cluster.ids) <- levels(seuratObj)

seuratObj <- Seurat::RenameIdents(seuratObj, new.cluster.ids)#seurat obj be found in UT box
#subset immune cells
nontumor<- subset(seuratObj, cells=colnames(seuratObj)[seuratObj@active.ident %in% c("other","microglia/monocyte","T-cell","EC")])
immunecellscounts <- as.matrix(Seurat::GetAssayData(nontumor, slot = "counts"))

nontumor$TumorType <- rep("", length(nontumor$seurat_clusters))

nontumor$TumorType[nontumor$orig.ident %in% c("BCM_Front", "BCM_Post", "MSC3")] <- "ClassC"

nontumor$TumorType[nontumor$orig.ident %in% c("MSC5")] <- "ClassB"

nontumor$TumorType[nontumor$orig.ident %in% c("MSC1", "MSC4")] <- "ClassA"


#####Cell Typist
###Majority Voting
majorityvoting<-predicted_labels[,4]
majorityvoting<-factor(majorityvoting)
Seurat::Idents(nontumor)<-majorityvoting


#Class A
classA<-subset(nontumor, subset = TumorType == "ClassA")

types<-levels(classA)
freq<-1:length(levels(classA))
for (ident in 1:length((levels(classA)))) {
  freq[ident]<-sum(classA@active.ident==types[ident])/length(classA@active.ident)
}


ClassAMajVot<-data.frame(types,freq)

pie(ClassAMajVot$freq,labels = ClassAMajVot$types)


#Class B
classB<-subset(nontumor, subset = TumorType == "ClassB")

types<-levels(classB)
freq<-1:length(levels(classB))
for (ident in 1:length((levels(classB)))) {
  freq[ident]<-sum(classB@active.ident==types[ident])/length(classB@active.ident)
}


ClassBMajVot<-data.frame(types,freq)
pie(ClassBMajVot$freq,labels = ClassBMajVot$types)


#Class C
classC<-subset(nontumor, subset = TumorType == "ClassC")

types<-levels(classC)
freq<-1:length(levels(classC))
for (ident in 1:length((levels(classC)))) {
  freq[ident]<-sum(classC@active.ident==types[ident])/length(classC@active.ident)
}


ClassCMajVot<-data.frame(types,freq)
pie(ClassCMajVot$freq,labels = ClassCMajVot$types)

#COMBINE TABLES
combined<-merge(ClassAMajVot,ClassBMajVot, by = "types")
majvotingfreqtable<-merge(combined,ClassCMajVot, by = "types")
colnames(majvotingfreqtable)[2:4]<-c("Class A","Class B","Class C")
write.csv(majvotingfreqtable,"~/Desktop/PatelLab/freqtable_majvoting.csv")


########Single R
#make predictions using single R
#Load in reference data set
hpca.se <- BlueprintEncodeData()

#Create object with predicted labels
pred.annotations <- SingleR(test = immunecellscounts, ref = hpca.se, assay.type.test=1,
                            labels = hpca.se$label.main)

singleR<-pred.annotations$pruned.labels
singleR<-factor(singleR)
Seurat::Idents(nontumor)<-singleR

pred.annotations$pruned.labels

#Class A
classA<-subset(nontumor, subset = TumorType == "ClassA")

types<-levels(classA)
freq<-1:length(levels(classA))
for (ident in 1:length((levels(classA)))) {
  freq[ident]<-sum(classA@active.ident==types[ident])/length(classA@active.ident)
}


ClassAMajVot<-data.frame(types,freq)

pie(ClassAMajVot$freq,labels = ClassAMajVot$types)


#Class B
classB<-subset(nontumor, subset = TumorType == "ClassB")

types<-levels(classB)
freq<-1:length(levels(classB))
for (ident in 1:length((levels(classB)))) {
  freq[ident]<-sum(classB@active.ident==types[ident])/length(classB@active.ident)
}


ClassBMajVot<-data.frame(types,freq)
pie(ClassBMajVot$freq,labels = ClassBMajVot$types)

table(nontumor@meta.data$TumorType)

#Class C
ClassC<-subset(nontumor, subset = TumorType == "ClassC")

types<-levels(ClassC)
freq<-1:length(levels(ClassC))
for (ident in 1:length((levels(ClassC)))) {
  freq[ident]<-sum(ClassC@active.ident==types[ident], na.rm = T)/length(ClassC@active.ident)
}

ClassCMajVot<-data.frame(types,freq)
pie(ClassCMajVot$freq,labels = ClassCMajVot$types)

levels(classC@active.ident)

#COMBINE TABLES
combined<-merge(ClassAMajVot,ClassBMajVot, by = "types")
singleRfreqtable<-merge(combined,ClassCMajVot, by = "types")
colnames(majvotingfreqtable)[2:4]<-c("Class A","Class B","Class C")
write.csv(majvotingfreqtable,"~/Desktop/PatelLab/freqtable_singleR.csv")


