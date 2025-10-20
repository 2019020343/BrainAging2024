rm(list=ls())
library("Seurat")
library(gplots)
library(ggplot2)
#install.packages("tidydr")
#devtools::install_github('junjunlab/scRNAtoolVis')
library(clustree)
library(scRNAtoolVis)
library(tidydr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
setwd("D:\\aging\\RESULT\\4val\\25region_rnaseq")
#################################Fig.5-A-C  #########################
sample_group <-read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\group.txt",header=T,sep = "\t", quote = "",fill = T,row.names = 1)
region_map<-read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\region_map.txt",header=T,sep = "\t", quote = "")
sample_group<-merge(region_map,sample_group,by.x="VAL_region",by.y="Sample_region")


sample_group$region1<-"NCX"
sample_group$region1[sample_group$TEST_region=="AMY"]="AMY"
sample_group$region1[sample_group$TEST_region=="MD"]="MD"
sample_group$region1[sample_group$TEST_region=="HIP"]="HIP"
sample_group$region1[sample_group$TEST_region=="STR"]="STR"
sample_group$region1[sample_group$TEST_region=="CBC"]="CBC"


expmean_gene<-read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\GSE211792_TPM_matrix.txt",header=T,sep = "\t", quote = "")
index<-sapply(sample_group$Sample_title,function(x){which(colnames(expmean_gene)==x)})
expmean_gene<-expmean_gene[,c(1,index)]
gene_symbol <- bitr(expmean_gene[,1], fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
genes_EXP<-merge(gene_symbol,expmean_gene,by.x="ENSEMBL",by.y="X")
gene_symbol_list<-unique(genes_EXP$SYMBOL)
expmean_gene_1<-c()
for (k in 1:length(gene_symbol_list)) {
  geneindex<-genes_EXP[which(genes_EXP$SYMBOL==gene_symbol_list[k]),-c(1,2)]
  geneindex1<-matrix(as.numeric(unlist(geneindex) ),ncol = ncol(geneindex)) 
  genei_mean<-cbind(gene_symbol_list[k] ,t(apply(geneindex1 , 2, mean)) )
  expmean_gene_1<-rbind(expmean_gene_1,genei_mean)
  print(paste(k,length(gene_symbol_list),sep = "---------"))
}
expmean_gene_1<-as.data.frame(expmean_gene_1)
colnames(expmean_gene_1)<-colnames(genes_EXP)[-1]

write.table(expmean_gene_1,"D:\\aging\\RESULT\\4val\\25region_rnaseq\\expmean_gene.txt",col.names = T, row.names = T,sep = "\t" ,append = FALSE, quote = F)
expmean_gene_1<-read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\expmean_gene.txt",header=T,sep = "\t", quote = "")

expmean_gene_2<-expmean_gene_1[,-1]
rownames(expmean_gene_2)<-expmean_gene_1[,1]

scRNA = CreateSeuratObject(counts=expmean_gene_2,
                           meta.data = sample_group,
                           min.cells = 3, 
                           min.features = 100)
head(scRNA@meta.data, 5)
#nCount_RNA与nFeature_RNA的相关性
pbmc <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


i=40
pbmc <- FindNeighbors(pbmc, dims = 1:i)
pbmc <- FindClusters(pbmc, resolution = 0.8)
# clustree(pbmc@meta.data,prefix = "RNA_snn_res.")
pbmc1 <- RunUMAP(pbmc, dims = 1:i, verbose = T)
umap2<- pbmc1@reductions[["umap"]]@cell.embeddings
rownames(sample_group)<-sample_group$Sample_title
umap3<-merge(umap2,sample_group,by="row.names")
write.table(umap3,"D:\\aging\\fig_data\\C.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


#dims = 1:10 即选取前10个主成分来分类细胞。
#查看前5个细胞的分类ID

#DimPlot(pbmc1, reduction = "umap")
p2=DimPlot(pbmc1, reduction = "umap",  group.by = 'TEST_region',label = F,cols = c("#12803b","#ec95b3","#f18e25","#f5c51e",
                                                                                   "#211d1e","#57a8d7","#ac536a","#735d96",
                                                                                   "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                                                                                   "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))

p3=DimPlot(pbmc1, reduction = "umap",  group.by = 'Sample_neuron',label = T,cols = c("#D43738","#6B4892"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))

p4=DimPlot(pbmc1, reduction = "umap",  group.by = 'region1',label = T,cols = c("#ec95b3","#f18e25","#211d1e","#cacbd0","#cf1223","#136cb6"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))

p5=DimPlot(pbmc1, reduction = "umap",  group.by = 'Sample_sex',label = F,cols = c("#cf1223","#136cb6"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))

p6=DimPlot(pbmc1, reduction = "umap",  group.by = 'Sample_age',label = T,label.color = "black",cols = c("#CB6D5D","#D56B19","#FAEACE","#F8B9AB",
                                                                                                 "#85AD9B","#747D90","#D93B45","#E32B2D",
                                                                                                 "#AFBED2","#8DBCCA","#827964","#577465",
                                                                                                 "#ED6767","#EA4840","#F9A27E","#FB9C42",
                                                                                                 "#C89B60","#504646","#b288ab","#d8d8dc",
                                                                                                 "#C7BA8E","#BE8B89","#C2B5CD","#C3CAAC"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))
p=p2+p3+p4+p5+p6

#p3
ggsave("D:\\aging\\RESULT\\4val\\25region_rnaseq\\all_umap.pdf",width = 20,height = 10)

############################### CAS ########################
library(limma)
library(gplots)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
# devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(ggplot2)
#install.packages("cols4all")
library(cols4all)
#BiocManager::install("dittoSeq")
library(dittoSeq)
setwd("D:\\aging\\RESULT\\4val\\25region_rnaseq\\")
gene<-read.table("GSE211792_TPM_matrix.txt",header=T,sep = "\t", quote = "")
sample_group<-read.table("group.txt",header=T,sep = "\t", quote = "")
region_map<-read.table("region_map.txt",header=T,sep = "\t", quote = "")
sample_group<-merge(region_map,sample_group,by.x="VAL_region",by.y="Sample_region")


gene1<-gene[,c(1,sapply(sample_group$Sample_title, function(x){which(colnames(gene)==x)  }))]

CAS_genes<-read.table("D:\\aging\\RESULT\\3CAS\\B_247genes.txt",header=T,sep = "\t", quote = "")
gene_symbol <- bitr(CAS_genes[,1], fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
CAS_genes_EXP<-merge(gene_symbol,gene1,by.x="ENSEMBL",by.y="X")
rownames(CAS_genes_EXP)<-CAS_genes_EXP[,2]
CAS_genes_EXP<-CAS_genes_EXP[,-(1:2)]
#CAS_genes_EXP1<-cbind(CAS_genes_EXP[,1:2],log2(CAS_genes_EXP[,-(1:2)]+0.000001))
library(VISION)
CAS_genes_sum<-apply(CAS_genes[,-1], 1, sum)
up<-cbind(CAS_genes[CAS_genes_sum>0,1],-1 ) 
down<-cbind(CAS_genes[CAS_genes_sum<0,1],1 ) 
sigData1<-t(rbind(up,down ))
sigData2<-as.numeric(sigData1[2,])
names(sigData2)<-sigData1[1,]
sig <- createGeneSignature(name = "Common Aging Score", sigData = sigData2)
mySignatures <- c(sig)
vis <- Vision(data = CAS_genes_EXP, signatures = mySignatures)
options(mc.cores = 1)
vis <- analyze(vis)
tsne <- getProjections(vis)[["tSNE30"]]
sigScores <- getSignatureScores(vis)[, "Common Aging Score"]
sigScores1<-cbind(names(sigScores),sigScores)
colnames(sigScores1)<-c("Sample_title","sigScores")
sigScores_group<-merge(sigScores1,sample_group,by= "Sample_title" )
write.table(sigScores_group,"D:\\aging\\RESULT\\4val\\sigScores_group.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

P3<-ggplot(sigScores_group,aes(x=as.character(Sample_age) ,y=as.numeric(sigScores) ,fill=as.character(Sample_age)))+ 
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(color="azure4",outlier.colour="red",
               outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
               )+  
  facet_wrap(~TEST_region)+
  stat_compare_means(method = "anova" )+
  #scale_fill_material_d()+
  theme_classic()+  
  #ylim(c(0.5,1))+
  scale_fill_manual(values = c("#B3332B","#B2522D","#E1A061","#A0C9D1","#436493"))+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Group", y="CAS")+ 
  theme(axis.title.x =element_text(size=12,face = "bold"), axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size = 12,face = "bold"))+
  #stat_summary(fun.y = "mean", geom = "point", size = 0.5) +
  stat_summary(fun.y = "mean", geom = "line", 
               aes(group = 1), 
               size = 0.5)

write.table(sigScores_group,"D:\\aging\\fig_data\\SF5\\C.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

################ regression  -test #########################

sigScores_group1<-read.table("D:\\aging\\RESULT\\2-3-CAS\\sigScores_group.txt",header=T,sep = "\t", quote = "")
sigScores_group_test<-sigScores_group1[sigScores_group1$age<40, ]

# 逻辑回归
library(PerformanceAnalytics)#加载包
lm1 <- lm(age~sigScores, data = sigScores_group_test)
su1 = summary(lm1)

ggplot(sigScores_group_test, aes(x=age, y=sigScores,color = region)) + 
  theme_classic()+
  ggtitle('linear fit')+
  geom_smooth(method = "lm", formula = y ~ x,size = 2,se = F)+ ##二项式拟合
  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                               "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))
write.table(sigScores_group_test,"D:\\aging\\fig_data\\SF5\\A.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

species_lst <- unique(sigScores_group_test$region)
p_lst <- list()
data_lst <- list()
formula <- y ~ x
for (i in species_lst) {
  extract_data <- sigScores_group_test[which(sigScores_group_test$region==i),]
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


################ regression  -val #########################
sigScores_group<-read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\sigScores_group.txt",header=T,sep = "\t", quote = "")
lm1 <- lm(Sample_age~sigScores, data = sigScores_group)
su1 = summary(lm1)
f1 <- y ~ x #定义回归方程
ggplot(sigScores_group, aes(x=Sample_age, y=as.numeric( sigScores),color = TEST_region)) + 
  theme_classic()+
  ggtitle('linear fit')+
  geom_smooth(method = "lm", formula = y ~ x,size = 2,se = F)+ ##二项式拟合
  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                             "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))
write.table(sigScores_group,"D:\\aging\\fig_data\\SF5\\A_V.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

library(ggpubr)
library(ggpmisc)
species_lst <- unique(sigScores_group$TEST_region)
p_lst <- list()
data_lst <- list()
formula <- y ~ x
for (i in species_lst) {
  extract_data <- sigScores_group[which(sigScores_group$TEST_region==i),]
  p <- ggplot(extract_data,aes(Sample_age,as.numeric(sigScores))) + 
    geom_point() + 
    geom_smooth( method = "lm") + 
    stat_cor(label.y = max(as.numeric(extract_data$sigScores))*0.90) + 
    stat_poly_eq(aes(label = ..eq.label..),label.y = max(as.numeric(extract_data$sigScores))*0.95,
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
#write.table(extract_ggplotdata,"D:\\aging\\RESULT\\4val\\slope.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
cols<-colorRampPalette(c("#5D779B", "#FAF3BE", "#C9301D"))(200)
extract_ggplotdata_test<-cbind(read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\slope.txt",header=T,sep = "\t", quote = ""),"test")
extract_ggplotdata_val<-cbind(read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\less40_slope.txt",header=T,sep = "\t", quote = ""),"val")
colnames(extract_ggplotdata_test)<-c("region","V2","p.value","estimate","slope","dataset")
colnames(extract_ggplotdata_val)<-c("region","V2","p.value","estimate","slope","dataset")
extract_ggplotdata<-rbind(extract_ggplotdata_test,extract_ggplotdata_val)
ggplot(extract_ggplotdata, aes(x =reorder(region,-slope) , y=slope*1000, fill=dataset)) +
  geom_bar(stat="identity", position="stack")+
  # geom_bar(position=position_dodge(), stat="identity",width = 0.8, color = NA, size = 0.5)+
  scale_fill_manual(values = c("#7491C3","#DE5F30"))+ 
  #scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  ylim(-6,15)+
  labs(x="",y="CAS slope (x 10 -3)")+　
  theme(axis.text.x = element_text(angle = 45))
write.table(extract_ggplotdata,"D:\\aging\\fig_data\\B.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
t.test(extract_ggplotdata$slope[extract_ggplotdata$dataset=="test"],extract_ggplotdata$slope[extract_ggplotdata$dataset=="val"],alternative = "two.sided")

##############################  neuron VS NON-neuron #################

sigScores_group<-read.table("D:\\aging\\RESULT\\4val\\sigScores_group.txt",header=T,sep = "\t", quote = "")
sigScores_group_neuron<-sigScores_group[sigScores_group$Sample_neuron=="group: non-neuron",]
library(ggpubr)
library(ggpmisc)
species_lst <- unique(sigScores_group_neuron$TEST_region)
p_lst <- list()
data_lst <- list()
formula <- y ~ x
for (i in species_lst) {
  extract_data <- sigScores_group_neuron[which(sigScores_group_neuron$TEST_region==i),]
  p <- ggplot(extract_data,aes(Sample_age,as.numeric(sigScores))) + 
    geom_point() + 
    geom_smooth( method = "lm") + 
    stat_cor(label.y = max(as.numeric(extract_data$sigScores))*0.90) + 
    stat_poly_eq(aes(label = ..eq.label..),label.y = max(as.numeric(extract_data$sigScores))*0.95,
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
#write.table(extract_ggplotdata,"D:\\aging\\RESULT\\4val\\slope_non-neuron.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
extract_ggplotdata_neuron<-cbind(read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\slope_neuron.txt",header=T,sep = "\t", quote = ""),"neuron")
extract_ggplotdata_non_neuron<-cbind(read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\slope_non-neuron.txt",header=T,sep = "\t", quote = ""),"non_neuron")
colnames(extract_ggplotdata_neuron)<-c("region","V2","p.value","estimate","slope","cell")
colnames(extract_ggplotdata_non_neuron)<-c("region","V2","p.value","estimate","slope","cell")
extract_ggplotdata<-rbind(extract_ggplotdata_neuron,extract_ggplotdata_non_neuron)
ggplot(extract_ggplotdata, aes(x =reorder(region,-slope) , y=slope*1000, fill=cell)) +
  geom_bar(stat="identity", position="stack")+
  # geom_bar(position=position_dodge(), stat="identity",width = 0.8, color = NA, size = 0.5)+
  scale_fill_manual(values = c("#C83937","#6B4892"))+ 
  #scale_y_continuous(expand = c(0,0)) +
  ylim(-6,15)+
  theme_classic()+
  labs(x="",y="CAS slope (x 10 -3)")+　
  theme(axis.text.x = element_text(angle = 45))
write.table(extract_ggplotdata,"D:\\aging\\fig_data\\D.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
t.test(extract_ggplotdata$slope[extract_ggplotdata$cell=="neuron"],extract_ggplotdata$slope[extract_ggplotdata$cell=="non_neuron"],alternative = "two.sided")

##############################  Sex #################

sigScores_group<-read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\sigScores_group.txt",header=T,sep = "\t", quote = "")
sigScores_group_Female<-sigScores_group[sigScores_group$Sample_sex=="Female",]
lm1 <- lm(Sample_age~sigScores, data = sigScores_group_Female)
su1 = summary(lm1)
f1 <- y ~ x #定义回归方程
ggplot(sigScores_group_Female, aes(x=Sample_age, y=as.numeric( sigScores),color = TEST_region)) + 
  theme_classic()+
  ggtitle('linear fit')+
  geom_smooth(method = "lm", formula = y ~ x,size = 2,se = F)+ ##二项式拟合
  theme(legend.position = 'bottom')+
  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                               "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))

library(ggpubr)
library(ggpmisc)
species_lst <- unique(sigScores_group_Female$TEST_region)[-12]
p_lst <- list()
data_lst <- list()
formula <- y ~ x
for (i in species_lst) {
  extract_data <- sigScores_group_Female[which(sigScores_group_Female$TEST_region==i),]
  p <- ggplot(extract_data,aes(Sample_age,as.numeric(sigScores))) + 
    geom_point() + 
    geom_smooth( method = "lm") + 
    stat_cor(label.y = max(as.numeric(extract_data$sigScores))*0.90) + 
    stat_poly_eq(aes(label = ..eq.label..),label.y = max(as.numeric(extract_data$sigScores))*0.95,
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
write.table(extract_ggplotdata,"D:\\aging\\RESULT\\4val\\25region_rnaseq\\slope_Male.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
cols<-colorRampPalette(c("#5D779B", "#FAF3BE", "#C9301D"))(200)


extract_ggplotdata_Female<-cbind(read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\slope_Female.txt",header=T,sep = "\t", quote = ""),"Female")
extract_ggplotdata_Male<-cbind(read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\slope_Male.txt",header=T,sep = "\t", quote = ""),"Male")
colnames(extract_ggplotdata_Male)<-c("region","V2","p.value","estimate","slope","sex")
colnames(extract_ggplotdata_Female)<-c("region","V2","p.value","estimate","slope","sex")
extract_ggplotdata<-rbind(extract_ggplotdata_Female,extract_ggplotdata_Male)
ggplot(extract_ggplotdata, aes(x =reorder(region,-slope) , y=slope*1000, fill=sex)) +
  geom_bar(stat="identity", position="stack")+
  # geom_bar(position=position_dodge(), stat="identity",width = 0.8, color = NA, size = 0.5)+
  scale_fill_manual(values = c("#D8996D","#77977A"))+ 
  #scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  ylim(-6,15)+
  labs(x="",y="CAS slope (x 10 -3)")+　
  theme(axis.text.x = element_text(angle = 45))
t.test(extract_ggplotdata$slope[extract_ggplotdata$sex=="Female"],extract_ggplotdata$slope[extract_ggplotdata$sex=="Male"])


########### CAS GENE ###################
CAS_genes_EXP1<-t(CAS_genes_EXP)
CAS_genes_EXP1<-cbind(rownames(CAS_genes_EXP1),CAS_genes_EXP1)
CAS_sample_group<-merge(CAS_genes_EXP1,sample_group,by.x="V1",by.y="Sample_title")
for (i in 2:221) {
  ggplot(CAS_sample_group, aes(fill=as.character(Sample_age) , y=as.numeric(CAS_sample_group[,i]), x=reorder(as.character(Sample_age),-as.numeric(CAS_sample_group[,i]) )))+
    geom_bar(position=position_dodge(),
             stat="summary",
             width=0.9,
             size=1)+ 
    stat_summary(fun.data = 'mean_se', 
                 geom = "errorbar", 
                 colour = "black",
                 width = 0.2,
                 position=position_dodge(0.7))+
    theme(legend.position = 'none')+
    labs(title = colnames(CAS_sample_group)[i], y="TPM expression", x = "")+
    scale_x_discrete(limits=c("20","30","34","35","39"))+
    scale_fill_manual(values = c("#B3332B","#B2522D","#E1A061","#A0C9D1","#436493"))+ 
    theme_classic()+
    theme(legend.position = 'none')+
    geom_signif(comparisons=list(c("20","30"),
                                 c("30","34"),
                                 c("34","35"),
                                 c("35","39")
                                 ),
                map_signif_level=T,##是否使用*显示显著性
                tip_length=c(0,0,0,0),
                y_position=c(50,60,70,80),
                size=1,textsize=2,
                test="t.test")+
    facet_wrap(~TEST_region)+
    geom_jitter(data = CAS_sample_group, aes(y = as.numeric(CAS_sample_group[,i]),x=reorder(as.character(Sample_age),-as.numeric(CAS_sample_group[,i]))),
                size = 3, shape = 16,
                color="grey60",
                stroke = 0.15, show.legend = FALSE, 
                position = position_jitterdodge(jitter.height=0.5,
                                                jitter.width = 0.1,
                                                dodge.width = 0.8))
  ggsave(paste("D:\\aging\\RESULT\\4val\\casgene\\",colnames(CAS_sample_group)[i],".pdf",sep="") )
}

###############xcell ###############
library(xCell)
expmean_gene<-read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\expmean_gene.txt",header=T,row.names = 2, sep = "\t", quote = "")
expmean_gene<-expmean_gene[,-1]
sigScores_group<-read.table("D:\\aging\\RESULT\\4val\\25region_rnaseq\\sigScores_group.txt",header=T,sep = "\t", quote = "")
scores <-  xCellAnalysis(expmean_gene,rnaseq = T)
xCell_scores<-as.data.frame(t(scores))
xCell_scores1<-cbind(rownames(xCell_scores),xCell_scores)
xCell_group<-merge(xCell_scores1,sigScores_group,by.x = "rownames(xCell_scores)",by.y ="Sample_title" )
library(tidyverse)       
library(ggsci)           
library(ggExtra)         
library(ggpmisc)        
library(palmerpenguins) 
library(ggpubr)
outdata<-c()
for (i in 2:68) {
  re<-cor.test( xCell_group[,i],xCell_group$sigScores,paired = T)
  R<-re$estimate
  P<-re$p.value
  out<-cbind(colnames(xCell_group)[i],R,P)
  outdata<-rbind(outdata,out)
}
outdata<-as.data.frame(outdata)
outdata_sig<-outdata[as.numeric(outdata$P) <0.01 & abs(as.numeric(outdata$R))>0.6, ]
outdata_sig<-outdata_sig[outdata_sig$V1!="MicroenvironmentScore" &  outdata_sig$V1!="StromaScore",]
ggplot(outdata_sig, aes(x =reorder(V1,-as.numeric(R)) , y=as.numeric(R), fill=V1)) +
  scale_fill_manual( values = c("#AFBED2","#D93B45","#CB6D5D","#F8B9AB","#85AD9B",
                                "#747D90","#8DBCCA","#577465","#ED6767","#827964",
                                "#FB9C42","#F9A27E","#C89B60","#504646","#b288ab",
                                "#d8d8dc","#C7BA8E","#BE8B89","#C2B5CD","#C3CAAC"))+
  #scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  #ylim(-5,10)+
  labs(x="",y="cor")+
  theme(axis.text.x = element_text(angle = 45))+#字体大小
  geom_col(width = 0.6, color = NA, size = 0.5)
write.table(outdata_sig,"D:\\aging\\fig_data\\E.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)










