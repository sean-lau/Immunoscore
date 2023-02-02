####Cluster Annotations
#seurat obj be found in UT box
load("~/Desktop/PatelLab/102022_SCELL_UCSF_BCM_WITH_newannotations_addmodulescore_seuratObj.rda")

sub<- subset(seuratObj, cells=colnames(seuratObj)[seuratObj$orig.ident %in% c("BCM_Front", "BCM_Post", "MSC3", "MSC5", "MSC1", "MSC4")])

sub$TumorType <- rep("", length(sub$seurat_clusters))

sub$TumorType[sub$orig.ident %in% c("BCM_Front", "BCM_Post", "MSC3")] <- "ClassC"

sub$TumorType[sub$orig.ident %in% c("MSC5")] <- "ClassB"

sub$TumorType[sub$orig.ident %in% c("MSC1", "MSC4")] <- "ClassA"



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

seuratObj <- Seurat::RenameIdents(seuratObj, new.cluster.ids)

##########

#Cell Typist Predicted annotations 
predictedlabels<-celltypist[,1]
predictedlabels<-factor(predictedlabels)
indexnum<-factor(1:83)
levels(predictedlabels) <- indexnum
print(predictedlabels)
Seurat::Idents(seuratObj)<-predictedlabels

#Condensed Annotations
condensedannotation<-CellTypistlegend[,4]
names(condensedannotation)<-levels(seuratObj)
seuratObj <- Seurat::RenameIdents(seuratObj, condensedannotation)

#Legend for Cell Typist prediccted annotations
predictedlabels<-celltypist[,1]
levels(predictedlabels)
CellTypistmain_legend<- data.frame(levels(predictedlabels),1:83)
colnames(CellTypistmain_legend)<- c("Cell Type", "Index Number")
write.csv(CellTypistmain_legend, file = "~/Desktop/PatelLab/CellTypistlegend.csv")

#######

#Cell Typist majority_voting
predictedlabels<-celltypist[,3]
predictedlabels<-factor(predictedlabels)
indexnum<-factor(1:83)
levels(predictedlabels) <- indexnum
print(predictedlabels)
Seurat::Idents(seuratObj)<-predictedlabels


#######

#SingleR
predictedlabels<-SingleRAnalysis[,1]
predictedlabels<-factor(predictedlabels)
indexnum<-factor(1:83)
levels(predictedlabels) <- indexnum
print(predictedlabels)
Seurat::Idents(seuratObj)<-predictedlabels

