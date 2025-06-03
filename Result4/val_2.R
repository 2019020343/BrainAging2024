rm(list=ls())
library(tidyverse)
library(SingleCellExperiment)
library(DropletUtils)
library(dplyr)
library(Seurat)
library(patchwork) 
#################################Fig.5-A-C  #########################

setwd("D:\\aging\\RESULT\\4val\\snrna\\GENECOUNT")
a<-list.files()
expmean_gene<-read.table(a[1],header=T,sep = "\t", quote = "",row.names = 1)
for(i in 2:length(a)){
  expmean_gene_i<-read.table(a[i],header=T,sep = "\t", quote = "",row.names = 1)
  expmean_gene<-cbind(expmean_gene,expmean_gene_i[,2])
 print(i)
}
colnames(expmean_gene)<-c("Symbols",as.character(sapply(a, function(x){strsplit(x,"_" )[[1]][1]})) ) 

genelist<-unique(expmean_gene[,1])
expmean_gene_1<-c()
for (k in 1:length(genelist)) {
  geneindex<-expmean_gene[which(expmean_gene[,1]==genelist[k]),-1]
  geneindex1<-matrix(as.numeric(unlist(geneindex) ),ncol = ncol(geneindex)) 
  genei_mean<-cbind(genelist[k] ,t(apply(geneindex1 , 2, mean)) )
  expmean_gene_1<-rbind(expmean_gene_1,genei_mean)
  print(paste(k,length(genelist),sep = "---------"))
}
expmean_gene_1<-as.data.frame(expmean_gene_1)
colnames(expmean_gene_1)<-colnames(expmean_gene)
expmean_gene_1<-expmean_gene_1[-which(is.na(expmean_gene_1[,1])==1) ,]
write.table(expmean_gene_1,"D:\\aging\\RESULT\\4val\\snrna\\expmean_gene.txt",col.names = T, row.names = T,sep = "\t" ,append = FALSE, quote = F)

expmean_gene_1 <-read.table("D:\\aging\\RESULT\\4val\\snrna\\expmean_gene.txt",header=T,sep = "\t", quote = "",fill = T)
rownames(expmean_gene_1)<-expmean_gene_1[,1]
expmean_gene_1<-expmean_gene_1[,-1]
sample_group<-read.table("D:\\aging\\RESULT\\4val\\snrna\\group.txt",header=T,sep = "\t", quote = "")
rownames(sample_group)<-sample_group$Sample_geo_accession
scRNA = CreateSeuratObject(counts=expmean_gene_1,
                           meta.data = sample_group,
                           min.cells = 3, 
                           min.features = 100)
head(scRNA@meta.data, 5)
#nCount_RNA与nFeature_RNA的相关性
pbmc <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 40)

i=40
pbmc <- FindNeighbors(pbmc, dims = 1:i)
pbmc <- FindClusters(pbmc, resolution = 0.8)
# clustree(pbmc@meta.data,prefix = "RNA_snn_res.")
pbmc1 <- RunUMAP(pbmc, dims = 1:i, verbose = T)
umap2<- pbmc1@reductions[["umap"]]@cell.embeddings
umap3<-merge(umap2,sample_group,by="row.names")
write.table(umap3,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\G.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


#DimPlot(pbmc1, reduction = "umap")
p3=DimPlot(pbmc1, reduction = "umap", pt.size = 2, group.by = 'celltype',label = T,cols = c("#D43738","#6B4892"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))
p3


############################### CAS ########################
library(limma)
library(gplots)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
library(dittoSeq)
library(ggpubr)
setwd("D:\\aging\\RESULT\\4val\\snrna\\")
expmean_gene_1<-read.table("expmean_gene.txt",header=T,sep = "\t", quote = "")
sample_group<-read.table("D:\\aging\\RESULT\\4val\\snrna\\group.txt",header=T,sep = "\t", quote = "")
CAS_genes<-read.table("D:\\aging\\RESULT\\3CAS\\B_247genes.txt",header=T,sep = "\t", quote = "")
CAS_genes_EXP<-merge(expmean_gene_1,CAS_genes,by.x="ID",by.y="Var1")[,1:ncol(gene)]
rownames(CAS_genes_EXP)<-CAS_genes_EXP[,1]
CAS_genes_EXP<-CAS_genes_EXP[,-1]
#CAS_genes_EXP<-cbind(CAS_genes_EXP[,1:2],log2(CAS_genes_EXP[,-(1:2)]+0.000001))
library(VISION)
CAS_genes_sum<-apply(CAS_genes[,-1], 1, sum)
up<-cbind(CAS_genes[CAS_genes_sum>0,1],-1 ) 
down<-cbind(CAS_genes[CAS_genes_sum<0,1],1 ) 
sigData1<-t(rbind(up,down ))
sigData2<-as.numeric(sigData1[2,])
names(sigData2)<-sigData1[1,]
sig <- createGeneSignature(name = "Common Aging Score", sigData = sigData2)
mySignatures <- c(sig)
vis <- Vision(data = CAS_genes_EXP, signatures = mySignatures,projection_methods="tSNE10")
options(mc.cores = 1)
vis <- analyze(vis)
tsne <- getProjections(vis)[["tSNE10"]]
sigScores <- getSignatureScores(vis)[, "Common Aging Score"]
sigScores1<-cbind(names(sigScores),sigScores)
colnames(sigScores1)<-c("Sample_geo_accession","sigScores")
sigScores_group<-merge(sigScores1,sample_group,by= "Sample_geo_accession" )
write.table(sigScores_group,"sigScores_group.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


sigScores_group<-read.table("sigScores_group.txt",header=T,sep = "\t", quote = "")
sigScores_group<-sigScores_group[sigScores_group$TEST_region == c("M1C","V1C","DFC"),]
write.table(sigScores_group,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\H.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
ggplot(sigScores_group,aes(x = celltype,y = sigScores)) +
  geom_violin(aes(color = celltype),cex=0.8,width=0.3)+
  geom_boxplot(aes(color = celltype),outlier.colour="black",width=0.1,cex=0.8)+
  geom_jitter(aes(color = celltype),width = 0.3,size=1.5)+
  #scale_fill_manual(values =  c("#D43738","#6B4892"))+
  scale_color_manual(values =  c("#D43738","#6B4892"))+
  #scale_x_discrete(limits=c("Young","Middle","Late"))+
  theme_bw()+
  facet_wrap(~TEST_region)+
  labs(x = "celltype", y = "CAS" ) +
  theme(panel.grid = element_blank())+
  stat_compare_means(method = "anova" )

species_lst <- unique(sigScores_group$celltype)
p_lst <- list()
data_lst <- list()
formula <- y ~ x
for (i in species_lst) {
  extract_data <- sigScores_group[which(sigScores_group$celltype==i),]
  p <- ggplot(extract_data,aes(age,sigScores)) + 
    geom_point() + 
    geom_smooth( method = "lm") + 
    stat_cor(label.y = max(extract_data$sigScores)*0.90) + 
    stat_poly_eq(aes(label = ..eq.label..),label.y = max(extract_data$sigScores)*0.95,
                 formula = formula,  parse = TRUE) + 
    theme_bw() + 
    ggtitle(i) + 
    theme(plot.title = element_text(hjust = 0.5,size=15))
  p_lst[[i]] <- p
  raw_furmula <- ggplot_build(p)$data[[4]]$label
  data_lst[[i]] <- c(i,gsub('[italic()~`* ]','',raw_furmula),
                     ggplot_build(p)$data[[3]]$p.value,ggplot_build(p)$data[[3]]$estimate)
}

extract_ggplotdata <- as.data.frame(t(as.data.frame(data_lst)))

write.table(extract_ggplotdata,"D:\\aging\\RESULT\\4val\\less40_slope.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
cols<-colorRampPalette(c("#5D779B", "#FAF3BE", "#C9301D"))(200)
extract_ggplotdata_less40_test<-read.table("D:\\aging\\RESULT\\4val\\less40_slope.txt",header=T,sep = "\t", quote = "")



