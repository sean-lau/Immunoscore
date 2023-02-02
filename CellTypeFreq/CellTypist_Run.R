#load in data
install.packages("Seurat")
library(Seurat)
load("~/Desktop/PatelLab/050922_SCELL_UCSF_BCM_WITH_CNV_seuratObj_ann.rda")

#export and write matrix to a csv file
immunecellscounts <- as.matrix(Seurat::GetAssayData(seuratObj, slot = "counts"))
write.csv(immunecellscounts,"~/Desktop/PatelLab/immunecellcounts.csv")

#LAST STEP:
#navigate to this website https://www.celltypist.org/ and have the option to upload csv file to their online server in order to run
#cell typist. Also, you have the option to download CellTypist via Github for larger files; however, the online server worked 
#well for the current matrix above (~500 MB took ~30 min). Received the 3 result files via email.