##HEAT Map
#make a heatmap where rows are the gene set ids, columns are 
#clusters (seuratObj$seurat_clusters) together with tumor class, for 
#example cluster1_A, cluster1_B, cluster_1_C, cluster2_A, …. and so on, 
#the matrix value is the mean module score of those cells, you can use 
#pheatmap. It would be a better way of visualizing the results. You can 
#also do the same for celltypist and singleR analysis, but now instead of 
#scores you can plot the number of cells annotated with that cell type.
#But here you have to normalize the numbers with the total number of cells 
#within each class.

install.packages('pheatmap')
library(pheatmap)


#Create a matrix
metadata<-sub@meta.data
clusters<-levels(sub$seurat_clusters)

table(sub$seurat_clusters)
#remove cluster that have 0 or 1 cells
clusters<-clusters[c(-18,-21,-22,-23,-24)]
#removed clusters 17, 20, 21, 22 (22 had 7 cells only), 23 (only class C)
print(clusters)
sub$TumorType<-factor(sub$TumorType)
TumorTypes<-levels(sub$TumorType)


ModuleScoreIndex<-which(colnames(metadata)=="GeneSet1"):122
for(x in 1:length(clusters)){
  if(x == 1){
    cluster<-subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters[x]])
  for(class in 1:length(TumorTypes)) {
    tumor<-subset(cluster, cells=colnames(cluster)[cluster$TumorType %in% TumorTypes[class]])
    if(class == 1){
  ##create vector with mean module index scores
     for(y in 1:length(ModuleScoreIndex)){
    if(y == 1){
      meansetscores<-mean(tumor@meta.data[,ModuleScoreIndex[y]])
    } else{
      meansetscores<-c(meansetscores,mean(tumor@meta.data[,ModuleScoreIndex[y]]))
    }
     }
      df<-data.frame(meansetscores)
    } else {
      for(y in 1:length(ModuleScoreIndex)){
        if(y == 1){
          meansetscores<-mean(tumor@meta.data[,ModuleScoreIndex[y]])
        } else{
          meansetscores<-c(meansetscores,mean(tumor@meta.data[,ModuleScoreIndex[y]]))
        }
        
    }
      df<-cbind(df,meansetscores)
    }
  }
  } else {
    cluster<-subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters[x]])
    for(class in 1:length(TumorTypes)){
      tumor<-subset(cluster, cells=colnames(cluster)[cluster$TumorType %in% TumorTypes[class]])
      for(y in 1:length(ModuleScoreIndex)){
        if(y == 1){
          meansetscores<-mean(tumor@meta.data[,ModuleScoreIndex[y]])
          print(ModuleScoreIndex[y])
        } else{
          meansetscores<-c(meansetscores,mean(tumor@meta.data[,ModuleScoreIndex[y]]))
        }
      }
      df<-cbind(df,meansetscores)
    }
   
  }
  print(clusters[x])
}

##assignn column and row names
#rownames variable "sets" is in the AddModuleScore.R script
rownames(df)<-sets

column<-"start"
Type<-c("A","B","C")
for(x in 1:length(clusters)){
  for(a in 1:3){
    column<-c(column,paste0("cluster","_",clusters[x],"_",Type[a]))
  }  
}
column<-column[-1]
colnames(df)<-column

## add clusters that have at least one class with cells
##cluster 23 C only (76 cells)
clusters<-levels(sub$seurat_clusters)
cluster<-subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters[24]])
cluster<-subset(cluster, cells=colnames(cluster)[cluster$TumorType %in% "ClassC"])
for(y in 1:length(ModuleScoreIndex)){
  if(y == 1){
    meansetscores<-mean(cluster@meta.data[,ModuleScoreIndex[y]])
    print(ModuleScoreIndex[y])
  } else{
    meansetscores<-c(meansetscores,mean(cluster@meta.data[,ModuleScoreIndex[y]]))
  }
}
cluster_23_C<-meansetscores
df<-cbind(df,cluster_23_C)


#clusters + types with 1 cell only: 24A, 22A, 22B, 20C, 17C, 15A

#rearrange columns by Tumor Class
num<-1:20
ClassAIndex<-num*3-2
ClassBIndex<-num*3-1
ClassCIndex<-num*3

colA<-df[,ClassAIndex]
colB<-df[,ClassBIndex]
colC<-df[,ClassCIndex]

dfTumorClass<-cbind(colA,colB,colC,"cluster_23_C" = df[,61])
##15A + 24A have only one cell
pheatmap(df, fontsize = 6, cluster_cols = FALSE)
pheatmap(dfTumorClass, fontsize = 6, cluster_cols = FALSE)


####Cell Typist and Single R##########
###Celltypist

#Define columns
clusters<-levels(sub$seurat_clusters)
clusters<-clusters[c(-18,-21,-22,-23,-24)]
sub$TumorType<-factor(sub$TumorType)
TumorTypes<-levels(sub$TumorType)
annotations<-levels(factor(seuratObj$CellTypistAnnotation))

#for loop to produce dataframe
for(x in 1:length(clusters)){
  if(x == 1){
    cluster<-subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters[x]])
    for(class in 1:length(TumorTypes)) {
      tumor<-subset(cluster, cells=colnames(cluster)[cluster$TumorType %in% TumorTypes[class]])
      if(class == 1){
        ##create vector with mean module index scores
        for(y in 1:length(annotations)){
          if(y == 1){
            freq<-sum(tumor$CellTypistAnnotation==annotations[y])/ncol(tumor)
          } else{
            freq<-c(freq,sum(tumor$CellTypistAnnotation==annotations[y])/ncol(tumor))
          }
        }
        df<-data.frame(freq)
      } else {
        for(y in 1:length(annotations)){
          if(y == 1){
            freq<-sum(tumor$CellTypistAnnotation==annotations[y])/ncol(tumor)
          } else{
            freq<-c(freq,sum(tumor$CellTypistAnnotation==annotations[y])/ncol(tumor))
          }
          
        }
        df<-cbind(df,freq)
      }
    }
  } else {
    cluster<-subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters[x]])
    for(class in 1:length(TumorTypes)){
      tumor<-subset(cluster, cells=colnames(cluster)[cluster$TumorType %in% TumorTypes[class]])
      for(y in 1:length(annotations)){
        if(y == 1){
          freq<-sum(tumor$CellTypistAnnotation==annotations[y])/ncol(tumor)
        } else{
          freq<-c(freq,sum(tumor$CellTypistAnnotation==annotations[y])/ncol(tumor))
        }
      }
      df<-cbind(df,freq)
    }
    
  }
  print(clusters[x])
}

##assignn column and row names

column<-"start"
Type<-c("A","B","C")
for(x in 1:length(clusters)){
  for(a in 1:3){
    column<-c(column,paste0("cluster","_",clusters[x],"_",Type[a]))
  }  
}
column<-column[-1]
colnames(df)<-column
rownames(df)<-annotations

##cluster 23 C only (76 cells)
clusters<-levels(sub$seurat_clusters)
cluster<-subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters[24]])
cluster<-subset(cluster, cells=colnames(cluster)[cluster$TumorType %in% "ClassC"])

for(y in 1:length(annotations)){
  if(y == 1){
    freq<-sum(cluster$CellTypistAnnotation==annotations[y])/ncol(cluster)
    
  } else{
    freq<-freq<-c(freq,sum(cluster$CellTypistAnnotation==annotations[y])/ncol(cluster))
  }
}
cluster_23_C<-freq
df<-cbind(df,cluster_23_C)
##cluster 22 C only (5 cells)
cluster<-subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters[23]])
cluster<-subset(cluster, cells=colnames(cluster)[cluster$TumorType %in% "ClassC"])

for(y in 1:length(annotations)){
  if(y == 1){
    freq<-sum(cluster$CellTypistAnnotation==annotations[y])/ncol(cluster)
    
  } else{
    freq<-freq<-c(freq,sum(cluster$CellTypistAnnotation==annotations[y])/ncol(cluster))
  }
}
cluster_22_C<-freq
df<-cbind(df,cluster_22_C)

#mark clusters with only 1 cell 5A and 24A?
colnames(df)<-gsub("cluster_5_A", "cluster_5_A**", colnames(df))
colnames(df)<-gsub("cluster_24_A", "cluster_24_A**", colnames(df))
#create pheatmap
pheatmap(df, fontsize = 10)

###Single R

#Define columns
clusters<-levels(sub$seurat_clusters)
clusters<-clusters[c(-18,-21,-22,-23,-24)]
sub$TumorType<-factor(sub$TumorType)
TumorTypes<-levels(sub$TumorType)
annotations<-levels(factor(seuratObj$SingleR))

#for loop to produce dataframe
for(x in 1:length(clusters)){
  if(x == 1){
    cluster<-subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters[x]])
    for(class in 1:length(TumorTypes)) {
      tumor<-subset(cluster, cells=colnames(cluster)[cluster$TumorType %in% TumorTypes[class]])
      if(class == 1){
        ##create vector with mean module index scores
        for(y in 1:length(annotations)){
          if(y == 1){
            freq<-sum(tumor$SingleR==annotations[y], na.rm = T)/ncol(tumor)
          } else{
            freq<-c(freq,sum(tumor$SingleR==annotations[y], na.rm = T)/ncol(tumor))
          }
        }
        df<-data.frame(freq)
      } else {
        for(y in 1:length(annotations)){
          if(y == 1){
            freq<-sum(tumor$SingleR==annotations[y], na.rm = T)/ncol(tumor)
          } else{
            freq<-c(freq,sum(tumor$SingleR==annotations[y], na.rm = T)/ncol(tumor))
          }
        }
        df<-cbind(df,freq)
      }
    }
  } else {
    cluster<-subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters[x]])
    for(class in 1:length(TumorTypes)){
      tumor<-subset(cluster, cells=colnames(cluster)[cluster$TumorType %in% TumorTypes[class]])
      for(y in 1:length(annotations)){
        if(y == 1){
          freq<-sum(tumor$SingleR==annotations[y], na.rm = T)/ncol(tumor)
        } else{
          freq<-c(freq,sum(tumor$SingleR==annotations[y], na.rm = T)/ncol(tumor))
          
        }
      }
      df<-cbind(df,freq)
    }
    
  }
  print(clusters[x])
}

#assign row and column names
rownames(df)<-annotations
column<-"start"
Type<-c("A","B","C")
for(x in 1:length(clusters)){
  for(a in 1:3){
    column<-c(column,paste0("cluster","_",clusters[x],"_",Type[a]))
  }  
}
column<-column[-1]
colnames(df)<-column

##cluster 23 C only (76 cells)
clusters<-levels(sub$seurat_clusters)
cluster<-subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters[24]])
cluster<-subset(cluster, cells=colnames(cluster)[cluster$TumorType %in% "ClassC"])

for(y in 1:length(annotations)){
  if(y == 1){
    freq<-sum(cluster$SingleR==annotations[y])/ncol(cluster)
    
  } else{
    freq<-c(freq,sum(cluster$SingleR==annotations[y])/ncol(cluster))
  }
}
cluster_23_C<-freq
df<-cbind(df,cluster_23_C)
##cluster 22 C only (5 cells)
clusters<-levels(sub$seurat_clusters)
cluster<-subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters[23]])
cluster<-subset(cluster, cells=colnames(cluster)[cluster$TumorType %in% "ClassC"])

for(y in 1:length(annotations)){
  if(y == 1){
    freq<-sum(cluster$SingleR==annotations[y])/ncol(cluster)
    
  } else{
    freq<-c(freq,sum(cluster$SingleR==annotations[y])/ncol(cluster))
  }
}
cluster_22_C<-freq
df<-cbind(df,cluster_22_C)

#mark clusters with only 1 cell 5A and 24A?
colnames(df)<-gsub("cluster_5_A", "cluster_5_A**", colnames(df))
colnames(df)<-gsub("cluster_24_A", "cluster_24_A**", colnames(df))


#create pheatmap
pheatmap(df, fontsize = 10)
