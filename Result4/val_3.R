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
#################################Fig.5-A-C  #########################
sample_group <-read.table("D:\\aging\\RESULT\\4val\\DFC\\group.txt",header=T,sep = "\t", quote = "",fill = T)
expmean_gene<-read.table("D:\\aging\\RESULT\\4val\\DFC\\gene_all.txt",header=T,sep = "\t", quote = "")
index<-sapply(sample_group$Sample_geo_accession,function(x){which(colnames(expmean_gene)==x)})
expmean_gene<-expmean_gene[,c(1,index)]
gene_symbol_list<-unique( expmean_gene$ID_REF)
gene_symbol_list<-gene_symbol_list[is.na(gene_symbol_list)!=1 ]
expmean_gene_1<-c()
for (k in 1:length(gene_symbol_list)) {
  geneindex<-expmean_gene[which(expmean_gene$ID_REF==gene_symbol_list[k]),-c(1)]
  geneindex1<-matrix(as.numeric(unlist(geneindex) ),ncol = ncol(geneindex)) 
  genei_mean<-cbind(gene_symbol_list[k] ,t(apply(geneindex1 , 2, mean)) )
  expmean_gene_1<-rbind(expmean_gene_1,genei_mean)
  print(paste(k,length(gene_symbol_list),sep = "---------"))
}
expmean_gene_1<-as.data.frame(expmean_gene_1)
colnames(expmean_gene_1)<-colnames(expmean_gene)
write.table(expmean_gene_1,"D:\\aging\\RESULT\\4val\\DFC\\expmean_gene.txt",col.names = T, row.names = T,sep = "\t" ,append = FALSE, quote = F)


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
setwd("D:\\aging\\RESULT\\4val\\DFC\\")
gene<-read.table("expmean_gene.txt",header=T,sep = "\t", quote = "")
sample_group<-read.table("group.txt",header=T,sep = "\t", quote = "")
CAS_genes<-read.table("D:\\aging\\RESULT\\3CAS\\B_247genes.txt",header=T,sep = "\t", quote = "")
CAS_genes_EXP<-merge(gene,CAS_genes,by.x="ID_REF",by.y="Var1")[,1:ncol(gene)]
rownames(CAS_genes_EXP)<-CAS_genes_EXP[,1]
CAS_genes_EXP<-CAS_genes_EXP[,-1]
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
vis <- Vision(data = CAS_genes_EXP, signatures = mySignatures,projection_methods="tSNE10")
options(mc.cores = 1)
vis <- analyze(vis)
tsne <- getProjections(vis)[["tSNE10"]]
sigScores <- getSignatureScores(vis)[, "Common Aging Score"]
sigScores1<-cbind(names(sigScores),sigScores)
colnames(sigScores1)<-c("Sample_geo_accession","sigScores")
sigScores_group<-merge(sigScores1,sample_group,by= "Sample_geo_accession" )

write.table(sigScores_group,"sigScores_group.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
sigScores_group <-read.table("D:\\aging\\RESULT\\4val\\DFC\\sigScores_group.txt",header=T,sep = "\t", quote = "",fill = T)
sigScores_group$group1=70
sigScores_group$group1[sigScores_group$Sample_characteristics_ch1.1 < 70]="60"
sigScores_group$group1[sigScores_group$Sample_characteristics_ch1.1 < 60]="50"
sigScores_group$group1[sigScores_group$Sample_characteristics_ch1.1 < 50]="40"
sigScores_group$group1[sigScores_group$Sample_characteristics_ch1.1 < 40]="30"
sigScores_group$group1[sigScores_group$Sample_characteristics_ch1.1 < 30]="20"

write.table(sigScores_group,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\J.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

P3<-ggplot(sigScores_group,aes(x=as.character(group1) ,y=as.numeric(sigScores) ,fill=as.character(group1)))+ 
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(color="azure4",outlier.colour="red",
               outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
  )+  
  
  #facet_wrap(~Sample_characteristics_ch1.2)+
  stat_compare_means(method = "anova" )+
  #scale_fill_material_d()+
  theme_classic()+  
  #ylim(c(0.5,1))+
  scale_fill_manual(values = c("#B3332B",
                               "#B2522D",
                               "#E1A061",
                               "#A0C9D1",
                               "#436493",
                               "#1C2C63"))+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Group", y="CAS")+ 
  theme(axis.title.x =element_text(size=12,face = "bold"), axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size = 12,face = "bold"))+
  #stat_summary(fun.y = "mean", geom = "point", size = 0.5) +
  stat_summary(fun.y = "mean", geom = "line", 
               aes(group = 1), 
               size = 0.5)
P3
P4<-ggplot(sigScores_group,aes(x=group1,y=as.numeric(sigScores) ,fill=Sample_characteristics_ch1.2))+ 
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(color="azure4",outlier.colour="red",
               outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
  )+  
  
  #facet_wrap(~region)+
  #stat_compare_means( method="t.test")+
  #scale_fill_material_d()+
  theme_classic()+  
  #ylim(c(0.5,1))+
  scale_fill_manual(values = c("#D8996D",
                               "#77977A"))+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Group", y="CAS")+ 
  theme(axis.title.x =element_text(size=12,face = "bold"), axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size = 12,face = "bold"))+
  #stat_summary(fun.y = "mean", geom = "point", size = 0.5) +
  stat_summary(fun.y = "mean", geom = "line", 
               aes(group = Sample_characteristics_ch1.2), 
               size = 0.5)

P4
write.table(sigScores_group,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF5\\E.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

###############xcell ###############
library(xCell)
library(reshape2)

expmean_gene<-read.table("D:\\aging\\RESULT\\4val\\DFC\\expmean_gene.txt",header=T,row.names = 2, sep = "\t", quote = "")
expmean_gene<-expmean_gene[,-1]
sigScores_group<-read.table("D:\\aging\\RESULT\\4val\\DFC\\sigScores_group.txt",header=T,sep = "\t", quote = "")
scores <-  xCellAnalysis(expmean_gene,rnaseq = T)
xCell_scores<-as.data.frame(t(scores))
xCell_scores1<-cbind(rownames(xCell_scores),xCell_scores)
xCell_scores2 <- melt(xCell_scores1,id.vars = c("rownames(xCell_scores)"))
xCell_group<-merge(xCell_scores2,sigScores_group,by.x = "rownames(xCell_scores)",by.y ="Sample_geo_accession" )
xCell_group$group[xCell_group$Sample_characteristics_ch1.1<12]="childhood"
xCell_group$group[xCell_group$Sample_characteristics_ch1.1>=12 & xCell_group$Sample_characteristics_ch1.1<20]="Adolescence"
xCell_group$group[xCell_group$Sample_characteristics_ch1.1>=20 & xCell_group$Sample_characteristics_ch1.1<40]="Young"
xCell_group$group[xCell_group$Sample_characteristics_ch1.1>=40 & xCell_group$Sample_characteristics_ch1.1<60]="Middle"
xCell_group$group[xCell_group$Sample_characteristics_ch1.1>60]="Late"
write.table(xCell_group,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF5\\F.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


library(ggplot2)
library(tidyverse)
library(gghalves) 
celltype<-unique(xCell_group$variable)

for (i in 1:67) {
  xCell_group1<-xCell_group[which(as.character(xCell_group$variable) ==celltype[i]),]

   ggplot(xCell_group1,aes(x = group,y = value)) +
     geom_violin(aes(fill = group),cex=0.8)+
    geom_boxplot(aes(fill = group),outlier.colour="black",width=0.1,cex=0.8)+
    geom_jitter(aes(color = group),width = 0.3,size=1.5)+
    scale_fill_manual(values =  c("#1f2a6b","#446799","#a3ced6"))+
    scale_color_manual(values =  c("#1f2a6b","#446799","#a3ced6"))+
    #scale_x_discrete(limits=c("Young","Middle","Late"))+
    theme_bw()+
    labs(x = "", y = celltype[i] ) +
    theme(panel.grid = element_blank())+
    stat_compare_means(method = "kruskal.test",
                        label = "p.format",
                        label.x = 2,
                        #label.y = 0.01,
                        show.legend = F)

  
  ggsave(paste("D:\\aging\\RESULT\\4val\\DFC\\xcell\\",celltype[i],".pdf",sep=""), width = 20, height = 20, units = "cm" )
  
}

write.table(xCell_group,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\K.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


#install.packages("ggExtra")
library(tidyverse)       
library(ggsci)           
library(ggExtra)         
library(ggpmisc)        
library(palmerpenguins) 
library(ggpubr)
xCell_group_point<-merge(xCell_scores1,sigScores_group,by.x = "rownames(xCell_scores)",by.y ="Sample_geo_accession" )
xCell_sigScores_cor<-c()
for (i in 2:68) {
  cor_re<-cor.test(xCell_group_point[,i],xCell_group_point$sigScores)
  p<-cor_re[[3]]
  r<-cor_re[[4]]
  out<-cbind(colnames(xCell_group_point)[i],p,r)
  xCell_sigScores_cor<-rbind(xCell_sigScores_cor,out)
}

xCell_sigScores_cor<-as.data.frame(xCell_sigScores_cor)
xCell_sigScores_cor_sig<-xCell_sigScores_cor[xCell_sigScores_cor$p<0.05,]
cols<-colorRampPalette(c("#5D779B", "#FAF3BE", "#C9301D"))(200)

write.table(xCell_sigScores_cor_sig,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\L.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
ggplot(xCell_sigScores_cor_sig, aes(x =reorder(V1,-as.numeric(r)) , y=as.numeric(r), fill=V1)) +
  scale_fill_manual( values = c("#CB6D5D","#F8B9AB","#85AD9B","#747D90",
                                "#D93B45","#AFBED2","#8DBCCA","#827964",
                                "#577465","#ED6767","#F9A27E","#FB9C42",
                                "#C89B60","#504646","#b288ab","#d8d8dc",
                                "#C7BA8E","#BE8B89","#C2B5CD","#C3CAAC"))+
  #scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  #ylim(-5,10)+
  labs(x="",y="cor")+
  theme(axis.text.x = element_text(angle = 45))+#字体大小
  geom_col(width = 0.6, color = NA, size = 0.5)
i=7
   ggplot(xCell_group_point, aes(sigScores, xCell_group_point[,i],color=Sample_characteristics_ch1)) +
      geom_point() +  # 添加散点图层，点的大小表示体重
       scale_colour_manual(values=c("#85AD9B")) +
       labs(x = "CAS", y = colnames(xCell_group_point)[i]) +  # 设置坐标轴标签
       theme_classic()  +
       stat_cor( geom = "text",label.x = -0.5)+
       theme(axis.title = element_text(family = "sans", face = "bold"))

i=10
   ggplot(xCell_group_point, aes(sigScores, xCell_group_point[,i],color=Sample_characteristics_ch1)) +
     geom_point() +  # 添加散点图层，点的大小表示体重
     scale_colour_manual(values=c("#D93B45")) +
     labs(x = "CAS", y = colnames(xCell_group_point)[i]) +  # 设置坐标轴标签
     theme_classic()  +
     stat_cor( geom = "text",label.x = -0.5)+
     theme(axis.title = element_text(family = "sans", face = "bold"))
write.table(xCell_group_point,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\L_POINT.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
   
