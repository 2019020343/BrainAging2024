rm(list = ls())
#devtools::install_github('satijalab/seurat-data')
#BiocManager::install("DescTools")
#BiocManager::install("glmGamPoi", type = "binary")

library(DescTools)

library(glmGamPoi)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(dplyr)
library(hdf5r)
library(patchwork)
dir = c('D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF241_C_ST\\outs/outs1/',#1
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF242_C_ST\\outs/outs1/',#2
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF242_T_ST\\outs/outs1/',#3
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF243_T_ST\\outs/outs1/',#4
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF248_C_ST\\outs/outs1/',#5
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF248_T_ST\\outs/outs1/',#6
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF251_T_ST\\outs/outs1/',#7
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF255_T_ST\\outs/outs1/',#8
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF256_C_ST\\outs/outs1/',#9
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF256_TC_ST\\outs/outs1/',#10
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF256_TI_ST\\outs/outs1/',#11
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF259_C_ST\\outs/outs1/',#12
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF259_T_ST\\outs/outs1/',#13
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF260_T_ST\\outs/outs1/',#14
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF262_T_ST\\outs/outs1/',#15
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF265_C_ST\\outs/outs1/',#16
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF265_T_ST\\outs/outs1/',#17
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF266_T_ST\\outs/outs1/',#18
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF268_IDHMutant_T_ST\\outs/outs1/',#19
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF269_T_ST\\outs/outs1/',#20
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF270_IDHMutant_T_ST\\outs/outs1/',#21
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF275_T_ST\\outs/outs1/',#22
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF296_T_ST\\outs/outs1/',#23
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF304_T_ST\\outs/outs1/',#24
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF313_C_ST\\outs/outs1/',#25
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF313_T_ST\\outs/outs1/',#26
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF334_C_ST\\outs/outs1/',#27
        'D:\\aging\\GBM\\cancer cell\\10XVisium_2\\10XVisium 2\\#UKF334_T_ST\\outs/outs1/'#28
)
names(dir) = c('#UKF241_C_ST', '#UKF242_C_ST', '#UKF242_T_ST', '#UKF243_T_ST',
               '#UKF248_C_ST', '#UKF248_T_ST', '#UKF251_T_ST', '#UKF255_T_ST',
               '#UKF256_C_ST',
               '#UKF256_TC_ST',
               '#UKF256_TI_ST','#UKF259_C_ST','#UKF259_T_ST', '#UKF260_T_ST', 
               '#UKF262_T_ST', 
               '#UKF265_C_ST','#UKF265_T_ST', 
               '#UKF266_T_ST','#UKF268_IDHMutant_T_ST','#UKF269_T_ST',
               '#UKF270_IDHMutant_T_ST','#UKF275_T_ST','#UKF296_T_ST','#UKF304_T_ST',
               '#UKF313_C_ST', '#UKF313_T_ST', '#UKF334_C_ST', '#UKF334_T_ST')




brain <- list()
for(i in 1:length(dir)){
  brain[[i]] <-Seurat::Load10X_Spatial(data.dir = dir[i])
  brain[[i]]@meta.data$orig.ident <-names(dir)[i]
}
for (i in 11:length(brain)) {
  brain[[i]][["percent.mt"]] <- PercentageFeatureSet(brain[[i]], pattern = "^MT[-]")
  brain[[i]] <- subset(brain[[i]], subset = nFeature_Spatial > 400 & nCount_Spatial > 1000 & percent.mt <50)
  brain[[i]] <- SCTransform(brain[[i]], assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)} ##SCT标准化

setwd("D:\\aging\\RESULT\\7mergedata") #设置输出路径


brain.merge <- merge(brain[[1]], y=c(brain[[2]],brain[[3]],brain[[4]],brain[[5]],
                                     brain[[6]],brain[[7]],brain[[8]],brain[[9]],
                                     #brain[[10]],
                                     brain[[11]],brain[[12]],brain[[13]],
                                     brain[[14]],brain[[15]],brain[[16]],brain[[17]],
                                     brain[[18]],brain[[19]],brain[[20]],brain[[21]],
                                     brain[[22]],brain[[23]],brain[[24]],brain[[25]],
                                     brain[[26]],brain[[27]],brain[[28]]
                                     
))




dim(brain.merge)
table(brain.merge@meta.data$orig.ident)
DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain[[1]]), VariableFeatures(brain[[2]]),
                                   VariableFeatures(brain[[3]]), VariableFeatures(brain[[4]]),
                                   VariableFeatures(brain[[5]]), VariableFeatures(brain[[6]]),
                                   VariableFeatures(brain[[7]]), VariableFeatures(brain[[8]]),
                                   VariableFeatures(brain[[9]]), 
                                   #VariableFeatures(brain[[10]]),
                                   VariableFeatures(brain[[11]]), VariableFeatures(brain[[12]]),
                                   VariableFeatures(brain[[13]]), VariableFeatures(brain[[14]]),
                                   VariableFeatures(brain[[15]]), VariableFeatures(brain[[16]]),
                                   VariableFeatures(brain[[17]]), VariableFeatures(brain[[18]]),
                                   VariableFeatures(brain[[19]]), VariableFeatures(brain[[20]]),
                                   VariableFeatures(brain[[21]]), VariableFeatures(brain[[22]]),
                                   VariableFeatures(brain[[23]]), VariableFeatures(brain[[24]]),
                                   VariableFeatures(brain[[25]]), VariableFeatures(brain[[26]]), 
                                   VariableFeatures(brain[[27]]),VariableFeatures(brain[[28]])
                                   )


brain.merge <- RunPCA(brain.merge,  assay = "SCT", verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, reduction = "pca", dims = 1:30)
brain.merge <- FindClusters(brain.merge,verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, reduction = "pca", dims = 1:30)


#save(brain,brain.merge,file = "Brain_ST.Rdata")



#brain.merge <- RunTSNE(brain.merge, reduction = "pca", dims = 1:30)
pdf("merge_umap_samples.pdf",width=10,height=5)
DimPlot(brain.merge, reduction = "umap", group.by = "orig.ident",cols = c("#CB6D5D","#D56B19","#FAEACE","#F8B9AB",
                                                                          "#85AD9B","#747D90","#D93B45","#E32B2D",
                                                                          "#AFBED2","#827964","#577465",
                                                                          "#ED6767","#EA4840","#F9A27E","#FB9C42",
                                                                          "#C89B60","#504646","#b288ab","#d8d8dc",
                                                                          "#C7BA8E","#BE8B89","#C2B5CD","#C3CAAC",
                                                                          "#ac536a","#735d96","#f1a7a4","#cf1223"))



dev.off()
pdf("merge_umap_group.pdf",width=10,height=5)
DimPlot(brain.merge, reduction = "umap", group.by = "orig.ident",cols = c("#a3ced6","#1f2a6b","#bb322b","#e8a361",
                                                                          "#446799","#e8a361","#bb322b","#bb322b",
                                                                          "#1f2a6b",
                                                                          "#bb322b",
                                                                          "#1f2a6b",
                                                                          "#bb322b","#bb322b","#e8a361","#446799",
                                                                          "#e8a361","#bb322b","#e8a361","#e8a361",
                                                                          "#e8a361","#bb322b","#e8a361","#bb322b",
                                                                          "#446799","#e8a361","#1f2a6b","#bb322b"))

dev.off()

##################  cortex  samples #################
load("Brain_ST.Rdata")
brain.merge_c <- merge(brain[[1]], y=c(brain[[5]],
                                      brain[[15]],
                                      brain[[24]],
                                      brain[[2]],
                                      brain[[9]],
                                      brain[[11]],
                                      brain[[26]]
))
DefaultAssay(brain.merge_c) <- "SCT"
VariableFeatures(brain.merge_c) <- c(VariableFeatures(brain[[1]]),VariableFeatures(brain[[5]]),VariableFeatures(brain[[15]]), VariableFeatures(brain[[24]]), # Y M
                                    VariableFeatures(brain[[2]]),VariableFeatures(brain[[9]]), VariableFeatures(brain[[11]]), VariableFeatures(brain[[26]]) # L
)



brain.merge_c <- RunPCA(brain.merge_c,  assay = "SCT", verbose = FALSE)
brain.merge_c <- FindNeighbors(brain.merge_c, reduction = "pca", dims = 1:30)
brain.merge_c <- FindClusters(brain.merge_c,verbose = FALSE)
brain.merge_c <- RunUMAP(brain.merge_c, reduction = "pca", dims = 1:30)


DimPlot(brain.merge_c, reduction = "umap", group.by = "orig.ident",cols = c("#FCF7F7","#CCCACA","#CCCACA","#CCCACA",
                                                                            "#000000","#000000","#000000","#000000"))





library(VISION)
hk_gene1<-read.table("D:\\aging\\RESULT\\2-3-CAS\\B_247genes.txt",header=T,sep = "\t", quote = "")

hk_gene1_sum<-apply(hk_gene1[,-1], 1, sum)
up<-cbind(hk_gene1[hk_gene1_sum>0,1],-1 ) 
down<-cbind(hk_gene1[hk_gene1_sum<0,1],1 ) 
sigData1<-t(rbind(up,down ))

sigData2<-as.numeric(sigData1[2,])
names(sigData2)<-sigData1[1,]

sig <- createGeneSignature(name = "Common Aging Score", sigData = sigData2)
mySignatures <- c(sig)


# Read in expression counts (Genes X Cells)
counts <-as.data.frame(brain.merge@assays[["SCT"]]@counts)
n.umi <- colSums(counts)

scaled_counts <- t(t(counts) / n.umi) * median(n.umi)

vis <- Vision(data = scaled_counts, signatures = mySignatures)


# Set the number of threads when running parallel computations
# On Windows, this must either be omitted or set to 1
options(mc.cores = 1)

vis <- analyze(vis)


tsne <- getProjections(vis)[["tSNE30"]]
sigScores <- getSignatureScores(vis)[, "Common Aging Score"]
sigScores<-as.data.frame(sigScores)
predictions.res <- t(sigScores ) 
rownames(predictions.res)<-"CAS"

cell_types <- rownames(predictions.res)
brain.merge <- AddMetaData(object= brain.merge, metadata = predictions.res[cell_types,], col.name = cell_types)



p1<-DimPlot(brain.merge, reduction = "umap", group.by = "orig.ident")

p2<-FeaturePlot(object = brain.merge, reduction = "umap", features = "CAS",cols = c("#45436E","#D8C7D0", "#9F4E4C")) 

p1+p2


SpatialFeaturePlot(brain.merge_c, features = "CAS")
             




##################  7 samples #################
setwd("D:\\aging\\RESULT\\7mergedata") #设置输出路径
load("Brain_cas.Rdata")
brain.merge1 <- merge(brain[[2]], y=c(brain[[3]],
                                      brain[[5]],brain[[6]],
                                      brain[[9]],brain[[10]],
                                      brain[[11]],brain[[12]],
                                      brain[[15]],brain[[16]],
                                      brain[[24]],brain[[25]],
                                      brain[[26]],brain[[27]]
))
DefaultAssay(brain.merge1) <- "SCT"
VariableFeatures(brain.merge1) <- c(VariableFeatures(brain[[2]]),VariableFeatures(brain[[3]]),
                                    VariableFeatures(brain[[5]]), VariableFeatures(brain[[6]]),
                                    VariableFeatures(brain[[9]]), VariableFeatures(brain[[10]]),
                                    VariableFeatures(brain[[11]]), VariableFeatures(brain[[12]]),
                                    VariableFeatures(brain[[15]]),VariableFeatures(brain[[16]]),
                                    VariableFeatures(brain[[24]]), VariableFeatures(brain[[25]]),
                                    VariableFeatures(brain[[26]]), VariableFeatures(brain[[27]])
)



brain.merge1 <- RunPCA(brain.merge1,  assay = "SCT", verbose = FALSE)
brain.merge1 <- FindNeighbors(brain.merge1, reduction = "pca", dims = 1:30)
brain.merge1 <- FindClusters(brain.merge1,verbose = FALSE)
brain.merge1 <- RunUMAP(brain.merge1, reduction = "pca", dims = 1:30)
brain.merge1 <- RunTSNE(brain.merge1, reduction = "pca", dims = 1:30)


library(VISION)
hk_gene1<-read.table("D:\\aging\\RESULT\\2-3-CAS\\B_247genes.txt",header=T,sep = "\t", quote = "")

hk_gene1_sum<-apply(hk_gene1[,-1], 1, sum)
up<-cbind(hk_gene1[hk_gene1_sum>0,1],-1 ) 
down<-cbind(hk_gene1[hk_gene1_sum<0,1],1 ) 
sigData1<-t(rbind(up,down ))

sigData2<-as.numeric(sigData1[2,])
names(sigData2)<-sigData1[1,]

sig <- createGeneSignature(name = "Common Aging Score", sigData = sigData2)
mySignatures <- c(sig)


# Read in expression counts (Genes X Cells)
counts <-as.data.frame(brain.merge1@assays[["SCT"]]@counts)
n.umi <- colSums(counts)

scaled_counts <- t(t(counts) / n.umi) * median(n.umi)

vis <- Vision(data = scaled_counts, signatures = mySignatures)


# Set the number of threads when running parallel computations
# On Windows, this must either be omitted or set to 1
options(mc.cores = 1)

vis <- analyze(vis)


tsne <- getProjections(vis)[["tSNE30"]]
sigScores <- getSignatureScores(vis)[, "Common Aging Score"]
sigScores<-as.data.frame(sigScores)
predictions.res <- t(sigScores ) 
rownames(predictions.res)<-"CAS"

cell_types <- rownames(predictions.res)
brain.merge1 <- AddMetaData(object= brain.merge1, metadata = predictions.res[cell_types,], col.name = cell_types)

umap2<- brain.merge1@reductions[["umap"]]@cell.embeddings

sample_group<-as.data.frame(cbind(brain.merge1@assays[["SCT"]]@data@Dimnames[[2]],brain.merge1@meta.data[["orig.ident"]],brain.merge1@meta.data[["CAS"]])) 
rownames(sample_group)<-sample_group$V1
umap3<-merge(umap2,sample_group,by="row.names")
write.table(umap3,"D:\\aging\\fig_data\\A-C.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

write.table(umap3,"D:\\aging\\fig_data\\E.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


p2<-DimPlot(brain.merge1, reduction = "umap", label = T, repel=T, alpha = 0.5,group.by = "orig.ident",cols = c("#ec95b3","#f18e25","#f5c51e",
                                                                                                      "#57a8d7","#ac536a","#735d96",
                                                                                                      "#5A7EB3","#cf1223","#f4da9a",
                                                                                                      "#f8bf89","#AFBED2","#827964",
                                                                                                      "#577465","#b288ab"))


p3<-DimPlot(brain.merge1, reduction = "umap", label = F, alpha = 0.5,group.by = "orig.ident",cols = c("#1f2a6b","#bb322b",#2-3
                                                                               "#446799","#e8a361",#5-6
                                                                               "#1f2a6b","#bb322b",#9 10 11
                                                                               "#1f2a6b","#bb322b",#9 10 11
                                                                               "#446799","#e8a361",#16  17
                                                                               "#446799","#e8a361",#25  26
                                                                               "#1f2a6b","#bb322b" #27  28
))




p4<-FeaturePlot(object = brain.merge1, reduction = "umap", alpha = 0.5,features = "CAS",cols = c("#5D779B", "#FAF3BE", "#C9301D")) 

p5<-DimPlot(brain.merge1, reduction = "umap", label = F, alpha = 0.5,group.by = "orig.ident",cols = c("#EAA3A0","#EAA3A0",#2-3
                                                                                                      "#EAA3A0","#EAA3A0",#5-6
                                                                                                      "#C3B3CD","#C3B3CD",#9 10
                                                                                                      "#C3B3CD","#C3B3CD",#11 12
                                                                                                      "#C3B3CD","#C3B3CD",#16  17
                                                                                                      "#C3B3CD","#C3B3CD",#24  25
                                                                                                      "#EAA3A0","#EAA3A0" #26  27
))
pdf("umap_7samples.pdf",width=7,height=5)  
p2
dev.off()

pdf("umap_7samples_group.pdf",width=7,height=5)  
p3
dev.off()

pdf("D:\\aging\\RESULT\\7mergedata\\umap_7samplesCAS.pdf",width=7,height=5)  
p4
dev.off()


pdf("D:\\aging\\RESULT\\7mergedata\\umap_7samples_region.pdf",width=7,height=5)  
p5
dev.off()

#save(brain.merge1,file = "Brain_7.Rdata")

################  cas #####################
library(circlize)

setwd("D:\\aging\\RESULT\\7mergedata") #设置输出路径
load("Brain_ST.Rdata")

library(VISION)
hk_gene1<-read.table("D:\\aging\\RESULT\\2-3-CAS\\B_247genes.txt",header=T,sep = "\t", quote = "")

hk_gene1_sum<-apply(hk_gene1[,-1], 1, sum)
up<-cbind(hk_gene1[hk_gene1_sum>0,1],-1 ) 
down<-cbind(hk_gene1[hk_gene1_sum<0,1],1 ) 
sigData1<-t(rbind(up,down ))
sigData2<-as.numeric(sigData1[2,])
names(sigData2)<-sigData1[1,]
sig <- createGeneSignature(name = "Common Aging Score", sigData = sigData2)
mySignatures <- c(sig)

for (i in 1:27) {
  # Read in expression counts (Genes X Cells)
  counts <-as.data.frame(brain[[i]]@assays[["SCT"]]@counts)
  
  # Scale counts within a sample
  n.umi <- colSums(counts)
  scaled_counts <- t(t(counts) / n.umi) * median(n.umi)
  vis <- Vision(data = scaled_counts, signatures = mySignatures)
  
  # Set the number of threads when running parallel computations
  # On Windows, this must either be omitted or set to 1
  options(mc.cores = 1)
  vis <- analyze(vis)
  tsne <- getProjections(vis)[["tSNE30"]]
  sigScores <- getSignatureScores(vis)[, "Common Aging Score"]
  sigScores1<-cbind(names(sigScores),sigScores)
  sigScores1<-as.data.frame( sigScores1[,2])
  colnames(sigScores1)<-"CAS"
  predictions.res <-t(sigScores1)
  rownames(predictions.res) <- make.names(rownames(predictions.res))
  cell_types <- rownames(predictions.res)
  
  
  brain[[i]] <- AddMetaData(object= brain[[i]], metadata =as.numeric( predictions.res[cell_types,] ) , col.name = cell_types)
  
  brain[[i]] <- RunPCA(brain[[i]],  assay = "SCT", verbose = FALSE)
  brain[[i]] <- FindNeighbors(brain[[i]], reduction = "pca", dims = 1:30)
  brain[[i]] <- FindClusters(brain[[i]],verbose = FALSE)
  brain[[i]] <- RunUMAP(brain[[i]], reduction = "pca", dims = 1:30)
  #brain[[i]] <- RunTSNE(brain[[i]], reduction = "pca", dims = 1:30)

}
#save(brain,file = "Brain_cas.Rdata")

min_num<-0
max_num<-0
for (i in 1:27) {
  min_num_i<- min(brain[[1]]@meta.data[["CAS"]])
  max_num_i<- max(brain[[1]]@meta.data[["CAS"]])
  min1<-c(min_num,min_num_i)
  min_num<-min1[which.min(min1)] #-0.1815214
  max1<-c(max_num,max_num_i)
  max_num<-max1[which.max(max1)] #0.7625968
  
}


for (i in 1:27) {

#DimPlot(brain[[i]], reduction = "tsne", label = T)
cols<-colorRampPalette(c("#5D779B", "#FAF3BE", "#C9301D"))(200)
pdf(paste(substr(brain[[i]]@meta.data[["orig.ident"]][1],2,9) ,".pdf",sep = "") ,width=10,height=5)  
p1<-SpatialPlot(brain[[i]], features  = "CAS",image.alpha = 0.3)+
  scale_fill_gradientn(limits=c(min_num, max_num), colours = cols)+
  labs(title = brain[[i]]@meta.data[["orig.ident"]][1])　

print(p1)
dev.off()
}

setwd("D:\\aging\\RESULT\\7mergedata") #设置输出路径
load("Brain_cas.Rdata")

CAS_number<-rbind(cbind(brain[[5]]@meta.data[["CAS"]],substr(brain[[5]]@meta.data[["orig.ident"]][1],9,9),substr(brain[[5]]@meta.data[["orig.ident"]][1],2,7) ),
                  cbind(brain[[6]]@meta.data[["CAS"]],substr(brain[[6]]@meta.data[["orig.ident"]][1],9,9),substr(brain[[6]]@meta.data[["orig.ident"]][1],2,7) ) ,
                  cbind(brain[[15]]@meta.data[["CAS"]],substr(brain[[15]]@meta.data[["orig.ident"]][15],9,9),substr(brain[[15]]@meta.data[["orig.ident"]][1],2,7) ),
                  cbind(brain[[16]]@meta.data[["CAS"]],substr(brain[[16]]@meta.data[["orig.ident"]][16],9,9),substr(brain[[16]]@meta.data[["orig.ident"]][1],2,7) ),
                  cbind(brain[[24]]@meta.data[["CAS"]],substr(brain[[24]]@meta.data[["orig.ident"]][24],9,9),substr(brain[[24]]@meta.data[["orig.ident"]][1],2,7) ),
                  cbind(brain[[25]]@meta.data[["CAS"]],substr(brain[[25]]@meta.data[["orig.ident"]][25],9,9),substr(brain[[25]]@meta.data[["orig.ident"]][1],2,7) ),
                  cbind(brain[[9]]@meta.data[["CAS"]],substr(brain[[9]]@meta.data[["orig.ident"]][9],9,9),substr(brain[[9]]@meta.data[["orig.ident"]][1],2,7) ),
                  cbind(brain[[10]]@meta.data[["CAS"]],substr(brain[[10]]@meta.data[["orig.ident"]][10],9,9),substr(brain[[10]]@meta.data[["orig.ident"]][1],2,7) ),
                  cbind(brain[[26]]@meta.data[["CAS"]],substr(brain[[26]]@meta.data[["orig.ident"]][26],9,9),substr(brain[[26]]@meta.data[["orig.ident"]][1],2,7) ),
                  cbind(brain[[27]]@meta.data[["CAS"]],substr(brain[[27]]@meta.data[["orig.ident"]][27],9,9),substr(brain[[27]]@meta.data[["orig.ident"]][1],2,7) ),
                  cbind(brain[[11]]@meta.data[["CAS"]],substr(brain[[11]]@meta.data[["orig.ident"]][11],9,9),substr(brain[[11]]@meta.data[["orig.ident"]][1],2,7) ),
                  cbind(brain[[12]]@meta.data[["CAS"]],substr(brain[[12]]@meta.data[["orig.ident"]][12],9,9),substr(brain[[12]]@meta.data[["orig.ident"]][1],2,7) ),
                  cbind(brain[[2]]@meta.data[["CAS"]],substr(brain[[2]]@meta.data[["orig.ident"]][2],9,9),substr(brain[[2]]@meta.data[["orig.ident"]][1],2,7) ),
                  cbind(brain[[3]]@meta.data[["CAS"]],substr(brain[[3]]@meta.data[["orig.ident"]][3],9,9),substr(brain[[3]]@meta.data[["orig.ident"]][1],2,7) )
                  
                  ) 
colnames(CAS_number)<-c("CAS","postion","sample")
write.table(CAS_number,"D:\\aging\\fig_data\\D.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

library(ggpubr) 
library(patchwork) 
library(ggsci)
library(tidyverse)
pdf("cas_boxplot.pdf" ,width=7,height=5)  
ggplot(CAS_number,aes(x = sample,y =as.numeric(CAS)  ,color = postion))+
  geom_boxplot(aes(fill=postion),alpha=0.5)+
  #geom_jitter(position = position_jitterdodge(jitter.height=0.02, # 散点抖动高度
  #                                            jitter.width = 0.15, # 散点抖动宽度
  #                                            dodge.width = 0.75),alpha=0.2)+ # x轴方向上的闪避量  #
  scale_fill_manual(values = c("#26496E","#CF521C"))+
  scale_color_manual(values = c("#26496E","#CF521C"))+
  theme_bw()+ 
  scale_x_discrete(limits=c("UKF248","UKF265","UKF313","UKF256","UKF334","UKF259","UKF242"))+
  theme(panel.grid = element_blank())+
  stat_compare_means(
    aes(group = postion),
    label = "p.signif",
    show.legend = F,
    label.y = 0.85)

dev.off()

i=19

umap3<-as.data.frame(brain[[19]]@meta.data) 
write.table(umap3,"D:\\aging\\fig_data\\F-H.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

p1 <- DimPlot(brain[[i]], reduction = "umap", label = TRUE,cols=c("#ec95b3","#f18e25","#f5c51e",
                                                                  "#57a8d7","#ac536a","#735d96",
                                                                  "#5A7EB3","#cf1223","#f4da9a",
                                                                  "#f8bf89"))
pdf("269_clusters.pdf",width=7,height=5)  
p1
dev.off()


p2 <- SpatialDimPlot(brain[[i]], label = TRUE, label.size = 3)+
  scale_fill_manual(values =c("#ec95b3","#f18e25","#f5c51e",
                              "#57a8d7","#ac536a","#735d96",
                              "#5A7EB3","#cf1223","#f4da9a",
                              "#f8bf89"))
pdf("269_clusters_HE.pdf",width=7,height=5)  
p2
dev.off()

CAS_269 <-rbind( cbind(subset(brain[[i]], idents = 6)@meta.data[["CAS"]],"C6") ,
cbind(subset(brain[[i]], idents = 2)@meta.data[["CAS"]],"C2") ,
cbind(subset(brain[[i]], idents = 8)@meta.data[["CAS"]],"C8") ,
cbind(subset(brain[[i]], idents = 3)@meta.data[["CAS"]],"C3") ,
cbind(subset(brain[[i]], idents = 7)@meta.data[["CAS"]],"C7") ,
cbind(subset(brain[[i]], idents = 5)@meta.data[["CAS"]],"C5") ,
 cbind(subset(brain[[i]], idents = 0)@meta.data[["CAS"]],"C0") ,
 cbind(subset(brain[[i]], idents = 4)@meta.data[["CAS"]],"C4") ,
 cbind(subset(brain[[i]], idents = 1)@meta.data[["CAS"]],"C1") ,
 cbind(subset(brain[[i]], idents = 9)@meta.data[["CAS"]],"C9") )

CAS_269<-as.data.frame(CAS_269)
colnames(CAS_269)<-c("CAS","Cluster")
CAS_269$CAS<-as.numeric(CAS_269$CAS)
write.table(CAS_269,"D:\\aging\\fig_data\\I.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

pdf("269_clusters_box.pdf",width=7,height=5)  
ggboxplot(CAS_269, x="Cluster", y="CAS", color="Cluster", add="jitter", legend="none") +
  rotate_x_text(angle = 45) +
  scale_fill_manual(values = c("#5A7EB3","#f5c51e","#f4da9a","#57a8d7",
                               "#cf1223","#735d96","#ec95b3","#ac536a",
                               "#f18e25","#f8bf89"))+
  scale_color_manual(values = c("#5A7EB3","#f5c51e","#f4da9a","#57a8d7",
                                "#cf1223","#735d96","#ec95b3","#ac536a",
                                "#f18e25","#f8bf89"))+
  theme_classic()+ 
  scale_x_discrete(limits=c("C6","C2","C8","C3",
                            "C7","C5","C0","C4",
                            "C1","C9"))+
  geom_hline(yintercept = mean(CAS_269$CAS),linetype=2) + # 添加base mean的水平线
  stat_compare_means(method = "anova", label.y = 0.62)+ # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") # Pairwise comparison against all
dev.off()

###############   CAS genes ############# 
setwd("D:\\aging\\RESULT\\7mergedata") #设置输出路径
load("Brain_cas.Rdata")







hk_gene1<-read.table("D:\\aging\\RESULT\\2-3-CAS\\B_247genes.txt",header=T,sep = "\t", quote = "")
MHCII_gene<-c("CD74","CTSD","FCER1G","HLA-DMA","HLA-DMB","HLA-DPA1","HLA-DRA","TREM2")
genelist<-intersect(as.character(hk_gene1$Var1), MHCII_gene)
j=5
for (j in 1:length(genelist) ) {
  genedata<-c()
  for (i in c(5,6,15,16,24,25,9,10,26,27,11,12,2,3)) {
    expdata<-as.data.frame(brain[[i]]@assays[["SCT"]]@data ) 
    meannum<-apply(expdata, 1, mean)
    meannum_order<-as.data.frame(meannum[ order(-meannum)] ) 
    datai<-cbind(brain[[i]]@meta.data[["orig.ident"]][1],as.numeric(expdata[which(rownames(expdata)==genelist[j]),] ) )
    genedata<-rbind(genedata,datai)
    print(i)
  }
  genedata<-as.data.frame(genedata)
  genedata$postion<-substr(genedata$V1,9,9)
  genedata$group<-substr(genedata$V1,2,7)
  maxnum<- 0
  for (k in c("T","C")) {
    data<-genedata[genedata$postion==k, ]
    data$V2<-as.numeric(data$V2)
    df1 <- data%>% group_by(group)%>%
      summarise(mean= mean(V2), sd= sd(V2))
    maxnum<-max(maxnum,max(df1$mean) ) 
    if(k=="T"){
      col1<-"#CF521C"
      
    }else{
      col1<-"#26496E"
      
    }
    df1$sample[df1$group == "UKF248"]="UKF248(44)"
    df1$sample[df1$group == "UKF265"]="UKF265(55)"
    df1$sample[df1$group == "UKF313"]="UKF313(57)"
    df1$sample[df1$group == "UKF256"]="UKF256(64)"
    df1$sample[df1$group == "UKF334"]="UKF334(73)"
    df1$sample[df1$group == "UKF259"]="UKF259(75)"
    df1$sample[df1$group == "UKF242"]="UKF242(81)"

    
      #c("#B3332B","#1C2C63","#EDB68A","#B2522D","#436493","#A0C9D1", "#E1A061")
    p1<- ggplot()+ 
      geom_bar(df1,mapping=aes(x=sample,y=mean), fill = col1,
               size = 1.5,color =col1,position="identity",
               stat="identity",width = 0.5)+
      scale_color_manual(values=col1)+ 
      theme_classic(base_line_size = 1)+
      labs(x="",y=genelist[j])+
      ylim(c(0,round(maxnum,2) ))+
      theme(axis.text.x = element_text(angle = 45))+#字体大,
      scale_x_discrete(limits=c("UKF248(44)","UKF265(55)","UKF313(57)","UKF256(64)","UKF334(73)","UKF259(75)","UKF242(81)"))
    write.table(df1,paste("D:\\aging\\fig_data\\M_",k,".txt",sep = "") ,col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
    
    #ggsave(paste("D:\\aging\\RESULT\\7mergedata\\",genelist[j],"_",k,".pdf", sep= ""), p1,width = 10, height = 10, units = "cm" )
    
  }
  
}

#rm(list = ls())
setwd("D:\\aging\\RESULT\\7mergedata") #设置输出路径
load("Brain_cas.Rdata")
hk_gene1<-read.table("D:\\aging\\RESULT\\2-3-CAS\\B_247genes.txt",header=T,sep = "\t", quote = "")

sample_cor_data<-c()
for (i in c(5,15,24,9,26,11,2)) {
  expdata1<-as.data.frame(brain[[i]]@assays[["SCT"]]@data ) 
  expdata2<-as.data.frame(brain[[i+1]]@assays[["SCT"]]@data ) 
  genelist1<-intersect(as.character(hk_gene1$Var1),rownames(expdata1)  )
  genelist<-intersect(genelist1,rownames(expdata2)  )
  set.seed(1)
  if(ncol(expdata1)>ncol(expdata2)){
  index<-sample(1:ncol(expdata1),ncol(expdata2), replace = F)
  expdata1<-expdata1[,index]
  }else{
    index<-sample(1:ncol(expdata2),ncol(expdata1), replace = F)
    expdata2<-expdata2[,index]}  
  gene_cor_data<-c()
  for (j in 1:length(genelist) ) {
  expdata1_genej<-expdata1[rownames(expdata1)==genelist[j],]
  expdata2_genej<-expdata2[rownames(expdata2)==genelist[j],]
  re<-cor.test(as.numeric(expdata1_genej) ,as.numeric(expdata2_genej))
  pvalue<- re$p.value
  rvalue<- re$estimate
  gene_cor_data1<-c(brain[[i]]@meta.data$orig.ident[1],brain[[i+1]]@meta.data$orig.ident[1],genelist[j],pvalue,rvalue)
  gene_cor_data<-rbind(gene_cor_data,gene_cor_data1)
  print(paste(j,length(genelist),sep = "/") )
  }
  sample_cor_data<-rbind(sample_cor_data,gene_cor_data)
  print(paste(i,j,length(genelist),sep = "/") )
}
sample_cor_data<-as.data.frame(sample_cor_data)
sample_cor_data_p<-drop_na(sample_cor_data)

colnames(sample_cor_data_p)<-c("cortex","tumor","gene","P","R")
sample_cor_data_p$sample<-substring(sample_cor_data_p$cortex,2,7 )
sample_cor_data_p$logp<- -log10(as.numeric(sample_cor_data_p$P ) )
write.table(sample_cor_data_p,"D:\\aging\\cor\\sample_cor_data.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

sample_cor_data_p$sample <- factor(sample_cor_data_p$sample, level=c("UKF248", "UKF265", "UKF313","UKF256", "UKF334", "UKF259", "UKF242"))


p1<- ggplot(sample_cor_data_p,aes(x=as.numeric(logp) ,fill=sample))+
  geom_density(position="stack")+
  scale_fill_manual(values=c("UKF248" = "#FCE3D7","UKF265"= "#F6B89F", "UKF313" = "#F08D70",
                             "UKF256" = "#EC674A","UKF334"= "#E13D2E", "UKF259" = "#BF1E20", "UKF242" = "#911D22"))+
  #guides(fill=guide_legend(title=NULL)) +
  theme_classic()+
  geom_vline(xintercept = 1.30103)

p1 
  
sample_cor_data_R<-sample_cor_data_p[sample_cor_data_p$P<0.05 ,]


p2<- ggplot(sample_cor_data_R,aes(x=as.numeric(R) ,fill=sample))+geom_density(position="stack")+
  scale_fill_manual(values=c("UKF248" = "#FCE3D7","UKF265"= "#F6B89F", "UKF313" = "#F08D70",
                             "UKF256" = "#EC674A","UKF334"= "#E13D2E", "UKF259" = "#BF1E20", "UKF242" = "#911D22"))+
  #scale_fill_brewer(palette = "Reds")+
  #guides(fill=guide_legend(title=NULL)) +
  theme_classic()

p2 


  ggsave(p1, file="D:\\aging\\cor\\P.pdf", width=12, height=10)
  ggsave(p2, file="D:\\aging\\cor\\R.pdf", width=12, height=10)
  
