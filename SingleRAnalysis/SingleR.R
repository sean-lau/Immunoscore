#SingleR WorkSpace

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("celldex")
library(celldex)
BiocManager::install("SingleR")
library(SingleR)

#Load in reference data set
hpca.se <- BlueprintEncodeData()

#Create object with predicted labels
pred.annotations <- SingleR(test = immunecellscounts, ref = hpca.se, assay.type.test=1,
                            labels = hpca.se$label.main)

#Subset object to consolidate analysis
SingleRAnalysis<-data.frame(pred.annotations@listData$pruned.labels, pred.annotations@listData$labels)
rownames(SingleRAnalysis)<-pred.annotations@rownames
write.csv(SingleRAnalysis,"~/Desktop/PatelLab/SingleRAnalysis.csv")
