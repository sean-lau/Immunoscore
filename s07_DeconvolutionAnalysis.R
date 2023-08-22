#Deconvolution Analysis
#Cibersort analysis################################################################################
########################################################################################################
#create a pheatmap with both LM22 and Module ratios
Cibserort_deconvolution_LM22_bulk <- read.csv("~/Desktop/PatelLab/NewAnalysis/Deconvolution//CIBERSORTx_LM22_Results.csv", row.names=1)
Cibersort_Module_deconvolution <- read.csv("~/Desktop/PatelLab/NewAnalysis/Deconvolution/CIBERSORTx_BCMsc_Results.csv", row.names=1)

#combined bulk data
load("~/Desktop/PatelLab/031023_Inhouse_Raleigh_Hongkong_GeneSymbol.rda")
data_GS<-data
GS<-row.names(data_GS)



#combine the two dataframes
Cibersort_Module_deconvolution<-Cibersort_Module_deconvolution[,c(-9,-10,-11)]
Cibserort_deconvolution_LM22_bulk<-Cibserort_deconvolution_LM22_bulk[,c(-24,-25,-23)]
deconv<-data.frame(Cibersort_Module_deconvolution,Cibserort_deconvolution_LM22_bulk)
deconv<-t(deconv)
deconv<-as.data.frame(deconv)


#Extract class data for samples"
#hongkong class data
totclass<-samplesheet$class
totclass<-as.data.frame(totclass)
names(totclass)<-"TumorClass"
row.names(totclass)<-row.names(samplesheet)

#combine class data from both data sets
sampleclass<-totclass

#create pheatmap
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
breaks <- seq(-2, 2, length = 16)
color <- colorRampPalette(rev(brewer.pal(11,'RdBu')))(length(breaks))
annoCol<-list(TumorClass=c(A='#28A860', B ='#0042ED', C='#FE4733'))
r<-pheatmap(deconv, annotation_col = sampleclass, show_colnames = F, cluster_rows = F, scale = "row",
            color= color,breaks = breaks, annotation_colors = annoCol)

ggsave('Heatmap_LM22andModule_bulkdatadeconvolution.pdf', plot = r, width=10, height=5)

##correlational matrix between LM22 and Module
Moduledeconv<-Cibersort_Module_deconvolution
LM22deconv<-Cibserort_deconvolution_LM22_bulk
cormatrix<-cor(Moduledeconv,LM22deconv)

#pheatmap cormatrix
p<-pheatmap(cormatrix)
ggsave("heatmap_cormatrix_Lm22vsModule.pdf", plot = p, width = 10, height = 5)
write.table(cormatrix, file = '~/Desktop/PatelLab/NewAnalysis/Deconvolution/cormatrix.txt', quote = F, sep = "\t")



#########################################Combined pheatmap analysis
##Cibersort Deconvolution
# load in Cibersort Module data
Moduledeconv<- read.csv("~/Desktop/PatelLab/NewAnalysis/Deconvolution/CIBERSORTx_BCMsc_Results.csv", row.names=1)
Moduledeconv<-Moduledeconv[,1:8]

Moduledeconv<-t(Moduledeconv)

row.names(Moduledeconv)<-gsub("^cluster", "C", row.names(Moduledeconv))
row.names(Moduledeconv)<-gsub("$", "_Cibersort", row.names(Moduledeconv))

#load in bayesprism deconvolution
bayesprism_deconvolution <- read.csv("~/Desktop/PatelLab/NewAnalysis/Deconvolution/bayesprism_combined_deconvolution.csv", row.names=1)
colnames(bayesprism_deconvolution)<-gsub("$", "_BP", colnames(bayesprism_deconvolution))
bayesprism_deconvolution<-t(bayesprism_deconvolution)

#merge bayes prism and cibersort deconvolution
BPandCB<-rbind(Moduledeconv, bayesprism_deconvolution)

#new prediction annotations
#load in sample sheet with predition annotations
load("~/Desktop/PatelLab/031023_Inhouse_Raleigh_Hongkong_GeneSymbol.rda")
samplesheet<-samplesheet[,2]
samplesheet<-as.data.frame(samplesheet)
colnames(samplesheet)<-"TumorClass"
row.names(samplesheet)<-colnames(Moduledeconv)

#recover free time annotations
#load in recover time data
HK_UCSF_BCM <- read.csv("~/Desktop/PatelLab/HK_UCSF_BCM.csv")
recover<-HK_UCSF_BCM[,c(24,4) ]
row.names(recover)<-recover[,2]
recover<-as.data.frame(recover[,1])
colnames(recover)<-"Recurrence Free Time"
row.names(recover)<-HK_UCSF_BCM$Original_name


########add tumor class to recover###################
recover$names<-row.names(recover)
samplesheet$ID<-row.names(samplesheet)
recovertumorclass<-merge(recover, samplesheet, by.x = "names", by.y = "ID")
row.names(recovertumorclass)<-recovertumorclass$names
recovertumorclass<-recovertumorclass[,-1]
recovertumorclass$names<-row.names(recovertumorclass)
recovertumorclass<-recovertumorclass[,-3]


##subset BPand CB to match reccovertumorclass
index<-""
for(x in 1:nrow(recovertumorclass)){
  index<-c(index, which(colnames(BPandCB)==row.names(recovertumorclass)[x]))
}
index<- index[-1]
index<-as.numeric(index)
BPandCB<-BPandCB[,index]

#create pheatmap
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
breaks <- seq(-2, 2, length = 16)
color <- colorRampPalette(rev(brewer.pal(11,'RdBu')))(length(breaks))
annoCol<-list(TumorClass=c(A='#28A860', B ='#0042ED', C='#FE4733'))

r<-pheatmap(BPandCB, annotation_col = recovertumorclass, scale = "row", color = color, breaks = breaks, fontsize = 8,
            annotation_colors = annoCol, cluster_rows = F, labels_col = F)
setwd("~/Desktop/PatelLab/NewAnalysis/")
ggsave('heatmap_modulebulkdata_deconvolution_withrecover.pdf', plot = r, width=30, height=5)

###Create geom boxplot (ggplot: facet_wrap with C1_Ciber, C1_BP, ... ###################################################
#ALSO x-axis with my Tumor classes, and y axis the ratios?#############################################################
#add tumorclass 
Moduledeconv<-BPandCB
Moduledeconv<-t(Moduledeconv)
neworder<-(sort(row.names(Moduledeconv)))
Moduledeconv<-Moduledeconv[neworder,]

neworder<-(sort(row.names(recovertumorclass)))
recovertumorclass<-recovertumorclass[neworder,]
Moduledeconv<-as.data.frame(Moduledeconv)

Moduledeconv$tumorclass<-as.factor(recovertumorclass$TumorClass)

#create new dataframe with three columns only: ratios, tumorclass, and module/software
for(x in 1:16){
  if(x ==1){
  boxplotdf<-data.frame("Ratios" = Moduledeconv[,x], "TumorClass" = Moduledeconv[,17], 
                        "Module_Program" = rep(colnames(Moduledeconv)[x], nrow(Moduledeconv)))
  }else{ 
    new<-data.frame("Ratios" = Moduledeconv[,x], "TumorClass" = Moduledeconv[,17], 
                    "Module_Program" = rep(colnames(Moduledeconv)[x], nrow(Moduledeconv)))
    boxplotdf<-rbind(boxplotdf, new)
  }
}

#create boxplot
library(ggplot2)
# Basic box plot
b<-ggplot(boxplotdf, aes(x=TumorClass, y=Ratios, fill = TumorClass)) + 
  geom_boxplot()+ 
  facet_wrap(~ Module_Program) + 
  labs(x = "Tumor Class", y = "Deconvolution Ratios")+
  scale_fill_manual(values=c('#28A860','#0042ED','#FE4733'))+
  theme(legend.position = "none")

b
ggsave('boxplot_deconvolutionmodulesbyclass.pdf', plot = b, width=10, height=8)
