rm(list=ls())
library(ggplot2)
library(ggalluvial)
library(tidyverse)
library(ggbreak)
library(dplyr)
################ SFig3-1 A #############
expmean_gene<-read.table("D:\\aging\\data\\geo\\expmean_gene.txt",header=T,sep = "\t", quote = "")
sample_group<-read.table("D:\\aging\\data\\geo\\sample_group.txt",header=T,sep = "\t", quote = "")

region <- unique(sample_group$region)

region1_sample_16<-expmean_gene[,1]
region1_sample1<-expmean_gene[,1]
#setwd("D:\\aging\\data\\geo\\region_DEG\\allgene\\F")
setwd("D:\\aging\\data\\geo\\region_DEG\\allgene\\M")


a<-list.files()

region_DEG_sig_name_out<-c()
region1_mean_out<-c()
ann_col1<-c()
for (i in 1:16) {
  region_DEG<-read.table(a[i],header=T,sep = "\t", quote = "")
  p<-region_DEG[,c(3,6,9)]
  region_DEG_sig<-region_DEG[which(apply(p,1,function(x){sum(as.numeric(x) <0.05) }  )>=2),]
  region_DEG_sig_name<-cbind(region_DEG_sig[,1],strsplit(as.character(a[i]),"_")[[1]][1] ) 
  region_DEG_sig_name_out<-rbind(region_DEG_sig_name_out,region_DEG_sig_name)
  
  
  region1<-sample_group[sample_group$region==region[i],1]
  region1_sample<-expmean_gene[,sapply(region1,function(x){which(colnames(expmean_gene)==x)} )]
  region1_sample1<-cbind(region1_sample1,region1_sample)
  
  ann_col<-rep(region[i],ncol(region1_sample) ) 
  ann_col1<-c(ann_col1,ann_col)
}
region_fre<-as.data.frame(table(region_DEG_sig_name_out[,1]))

region_fre1<-as.data.frame(table(region_fre[,2]))
#write.table(region_fre1,"D:\\aging\\RESULT\\2-3-CAS\\region_fre1_M.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



p1<-ggplot(region_fre1, aes(x = Var1, y=as.numeric(Freq))) +
  #scale_fill_manual(values = c()) +
  #scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  labs(x="the number of regions in which a DEG was detected",y="# DEGs")+
  #geom_col(width = 0.6, fill= "#D8996D", color = "#D8996D", size = 0.5)
  
  geom_col(width = 0.6, fill= "#77977A", color = "#77977A", size = 0.5)

p1

region_fre2<-region_fre1[as.numeric( region_fre1$Var1)>=10,]

p2<-ggplot(region_fre2, aes(x = Var1, y=as.numeric(Freq))) +
  #scale_fill_manual(values = c()) +
  #scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  labs(x="the number of regions in which a DEG was detected",y="# DEGs")+
  #geom_col(width = 0.6, fill= "#D8996D", color = "#D8996D", size = 0.5)

 geom_col(width = 0.6, fill= "#77977A", color = "#77977A", size = 0.5)
p2




#################  SFig3-1  B ##################v
setwd("D:\\aging\\data\\geo\\region_DEG\\allgene\\F")
a<-list.files()

F_OUT<-c()
p_OUT<-c()
for (i in 1:16) {
  region_DEG_F<-read.table(paste("D:\\aging\\data\\geo\\region_DEG\\allgene\\F",a[i],sep = "\\"),header=T,sep = "\t", quote = "")
  region_DEG_M<-read.table(paste("D:\\aging\\data\\geo\\region_DEG\\allgene\\M",a[i],sep = "\\"),header=T,sep = "\t", quote = "")
  YM_F_sig_up<-region_DEG_F[region_DEG_F$YM_P.Value<0.05 & region_DEG_F$YM_logFC >0 ,1]
  YM_F_sig_down<-region_DEG_F[region_DEG_F$YM_P.Value<0.05 & region_DEG_F$YM_logFC <0 ,1]
  YM_M_sig_up<-region_DEG_M[region_DEG_M$YM_P.Value<0.05 & region_DEG_M$YM_logFC >0 ,1]
  YM_M_sig_down<-region_DEG_M[region_DEG_M$YM_P.Value<0.05 & region_DEG_M$YM_logFC <0 ,1]
  UU<-length(intersect(YM_F_sig_up,YM_M_sig_up)) 
  UD<-length(intersect(YM_F_sig_up,YM_M_sig_down)) 
  DU<-length(intersect(YM_F_sig_down,YM_M_sig_up)) 
  DD<-length(intersect(YM_F_sig_down,YM_M_sig_down)) 
  
  data <- matrix(c(UU, DU, UD, DD), nrow = 2)
  
  # 进行R Fisher精确检验
  YM_result <- fisher.test(data)
  YM_p.value<-YM_result$p.value
  YM_Fraction<-(UU+DD)/(UU+DU+UD+DD)
  
  
  
  YL_F_sig_up<-region_DEG_F[region_DEG_F$YL_P.Value<0.05 & region_DEG_F$YL_logFC >0 ,1]
  YL_F_sig_down<-region_DEG_F[region_DEG_F$YL_P.Value<0.05 & region_DEG_F$YL_logFC <0 ,1]
  YL_M_sig_up<-region_DEG_M[region_DEG_M$YL_P.Value<0.05 & region_DEG_M$YL_logFC >0 ,1]
  YL_M_sig_down<-region_DEG_M[region_DEG_M$YL_P.Value<0.05 & region_DEG_M$YL_logFC <0 ,1]
  UU<-length(intersect(YL_F_sig_up,YL_M_sig_up)) 
  UD<-length(intersect(YL_F_sig_up,YL_M_sig_down)) 
  DU<-length(intersect(YL_F_sig_down,YL_M_sig_up)) 
  DD<-length(intersect(YL_F_sig_down,YL_M_sig_down)) 
  
  data <- matrix(c(UU, DU, UD, DD), nrow = 2)
  
  # 进行R Fisher精确检验
  YL_result <- fisher.test(data)
  YL_p.value<-YL_result$p.value
  YL_Fraction<-(UU+DD)/(UU+DU+UD+DD)

  
  
  ML_F_sig_up<-region_DEG_F[region_DEG_F$ML_P.Value<0.05 & region_DEG_F$ML_logFC >0 ,1]
  ML_F_sig_down<-region_DEG_F[region_DEG_F$ML_P.Value<0.05 & region_DEG_F$ML_logFC <0 ,1]
  ML_M_sig_up<-region_DEG_M[region_DEG_M$ML_P.Value<0.05 & region_DEG_M$ML_logFC >0 ,1]
  ML_M_sig_down<-region_DEG_M[region_DEG_M$ML_P.Value<0.05 & region_DEG_M$ML_logFC <0 ,1]
  UU<-length(intersect(ML_F_sig_up,ML_M_sig_up)) 
  UD<-length(intersect(ML_F_sig_up,ML_M_sig_down)) 
  DU<-length(intersect(ML_F_sig_down,ML_M_sig_up)) 
  DD<-length(intersect(ML_F_sig_down,ML_M_sig_down)) 
  
  data <- matrix(c(UU, DU, UD, DD), nrow = 2)
  
  # 进行R Fisher精确检验
  ML_result <- fisher.test(data)
  ML_p.value<-ML_result$p.value
  ML_Fraction<-(UU+DD)/(UU+DU+UD+DD)

  p_OUT1<-rbind(YM_p.value,YL_p.value,ML_p.value)
  colnames(p_OUT1)<-strsplit(as.character(a[i]),"_")[[1]][1]
  p_OUT<-cbind(p_OUT,p_OUT1)
  
  
  
  F_OUT1<-rbind(YM_Fraction,YL_Fraction,ML_Fraction)
  colnames(F_OUT1)<-strsplit(as.character(a[i]),"_")[[1]][1]
  F_OUT<-cbind(F_OUT,F_OUT1)

}
write.table(p_OUT,"D:\\aging\\RESULT\\3CAS\\p_OUT.txt",col.names = T, row.names = T,sep = "\t" ,append = FALSE, quote = F)

write.table(F_OUT,"D:\\aging\\RESULT\\3CAS\\F_OUT.txt",col.names = T, row.names = T,sep = "\t" ,append = FALSE, quote = F)


library(pheatmap)
myColor <- colorRampPalette(c("#FAF3BE", "#C9301D"))(200)
pheatmap(F_OUT, 
         #annotation_col = annotation_col, 
         #annotation_row = annotation_row,
         cluster_rows = F,
         cluster_cols = F,
         cellwidth = 20, cellheight = 25,
         color = myColor,
         #breaks=myBreaks, 
         annotation_legend = FALSE)  

##############  SFig3-1  C   ###############

region_DEG_F<-read.table("D:\\aging\\data\\geo\\region_DEG\\allgene\\F\\A1C_re_out5.txt",header=T,sep = "\t", quote = "")
region_DEG_M<-read.table("D:\\aging\\data\\geo\\region_DEG\\allgene\\M\\A1C_re_out5.txt",header=T,sep = "\t", quote = "")

YM_F_sig<-region_DEG_F[region_DEG_F$YM_P.Value<0.05 ,c(1,2)] 
YM_M_sig<-region_DEG_M[region_DEG_M$YM_P.Value<0.05 ,c(1,2)]

ss<-as.data.frame(intersect(YM_F_sig[,1],YM_M_sig[,1]) ) 
region_DEG_s<-merge(ss,region_DEG_F,by.x="intersect(YM_F_sig[, 1], YM_M_sig[, 1])",by.y="gene")[,1:2]
region_DEG_ss<-merge(region_DEG_s,region_DEG_M,by.x="intersect(YM_F_sig[, 1], YM_M_sig[, 1])",by.y="gene")[,1:3]
colnames(region_DEG_ss)<-c("gene","logFC_F","logFC_M")
region_DEG_ss$group<-region_DEG_ss$logFC_F*region_DEG_ss$logFC_M
region_DEG_ss$group1[region_DEG_ss$group>0]="c"
region_DEG_ss$group1[region_DEG_ss$group<0]="d"


p1<-ggplot(region_DEG_ss, aes(logFC_F, logFC_M,colour = group1)) + 
  geom_point(size = 2)+
  theme_classic()+  
  ylim(c(-2,2))+
  labs(title = "A1C_YM")+ 
  theme(legend.position = "none")+
  scale_color_manual(values = c("#B2522D",
                               "#436493"))

p1
write.table(region_DEG_ss,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\B_A1C_YM.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



YL_F_sig<-region_DEG_F[region_DEG_F$YL_P.Value<0.05 ,c(1,2)] 
YL_M_sig<-region_DEG_M[region_DEG_M$YL_P.Value<0.05 ,c(1,2)]

YL_ss<-as.data.frame(intersect(YL_F_sig[,1],YL_M_sig[,1]) ) 
YL_region_DEG_s<-merge(YL_ss,region_DEG_F,by.x="intersect(YL_F_sig[, 1], YL_M_sig[, 1])",by.y="gene")[,c(1,5)]
YL_region_DEG_ss<-merge(YL_region_DEG_s,region_DEG_M,by.x="intersect(YL_F_sig[, 1], YL_M_sig[, 1])",by.y="gene")[,c(1,2,6)]
colnames(YL_region_DEG_ss)<-c("gene","logFC_F","logFC_M")
YL_region_DEG_ss$group<-YL_region_DEG_ss$logFC_F*YL_region_DEG_ss$logFC_M
YL_region_DEG_ss$group1[YL_region_DEG_ss$group>0]="c"
YL_region_DEG_ss$group1[YL_region_DEG_ss$group<0]="d"


p2<-ggplot(YL_region_DEG_ss, aes(logFC_F, logFC_M,colour = group1)) + 
  geom_point(size = 2)+
  theme_classic()+  
  ylim(c(-2,2))+
  labs(title = "A1C_YL")+ 
  theme(legend.position = "none")+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))
write.table(YL_region_DEG_ss,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\B_A1C_YL.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

p2
ML_F_sig<-region_DEG_F[region_DEG_F$ML_P.Value<0.05 ,c(1,2)] 
ML_M_sig<-region_DEG_M[region_DEG_M$ML_P.Value<0.05 ,c(1,2)]

ss<-as.data.frame(intersect(ML_F_sig[,1],ML_M_sig[,1]) ) 
region_DEG_s<-merge(ss,region_DEG_F,by.x="intersect(ML_F_sig[, 1], ML_M_sig[, 1])",by.y="gene")[,c(1,8)]
region_DEG_ss<-merge(region_DEG_s,region_DEG_M,by.x="intersect(ML_F_sig[, 1], ML_M_sig[, 1])",by.y="gene")[,c(1,2,9)]
colnames(region_DEG_ss)<-c("gene","logFC_F","logFC_M")
region_DEG_ss$group<-region_DEG_ss$logFC_F*region_DEG_ss$logFC_M
region_DEG_ss$group1[region_DEG_ss$group>0]="c"
region_DEG_ss$group1[region_DEG_ss$group<0]="d"


p3<-ggplot(region_DEG_ss, aes(logFC_F, logFC_M,colour = group1)) + 
  geom_point(size = 2)+
  theme_classic()+  
  ylim(c(-2,2))+
  labs(title = "A1C_ML")+ 
  theme(legend.position = "none")+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))


p3

write.table(region_DEG_ss,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\B_A1C_ML.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

#############  clusterProfiler  #########################
rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(topGO)
setwd("D:\\aging\\data\\geo\\region_DEG\\allgene\\F")

setwd("D:\\aging\\data\\geo\\region_DEG\\allgene\\M")

a<-list.files()

expmean_gene<-read.table("D:\\aging\\data\\geo\\expmean_gene.txt",header=T,sep = "\t", quote = "")
sample_group<-read.table("D:\\aging\\data\\geo\\sample_group.txt",header=T,sep = "\t", quote = "")

expmean_gene1<-expmean_gene[,c(1,sapply(sample_group$samples,function(x){which(colnames(expmean_gene)==x)}))] 


region <- unique(sample_group$region)
region_DEG_sig_name_out<-c()
region1_mean_out<-c()
region1_sample_16<-expmean_gene1[,1]
region1_sample1<-expmean_gene1[,1]
ann_col1<-c()
for (i in 1:16) {
  region_DEG<-read.table(a[i],header=T,sep = "\t", quote = "")
  p<-region_DEG[,c(3,6,9)]
  region_DEG_sig<-region_DEG[which(apply(p,1,function(x){sum(as.numeric(x) <0.05) }  )>=2),]
  region_DEG_sig_name<-cbind(region_DEG_sig[,1],strsplit(as.character(a[i]),"_")[[1]][1] ) 
  region_DEG_sig_name_out<-rbind(region_DEG_sig_name_out,region_DEG_sig_name)
  
  
  region1<-sample_group[sample_group$region==region[i],1]
  region1_sample<-expmean_gene1[,sapply(region1,function(x){which(colnames(expmean_gene1)==x)} )]
  region1_sample1<-cbind(region1_sample1,region1_sample)
  
  ann_col<-rep(region[i],ncol(region1_sample) ) 
  ann_col1<-c(ann_col1,ann_col)
}
region_fre<-as.data.frame(table(region_DEG_sig_name_out[,1]))
region_fre_10<-region_fre[region_fre$Freq>=10,]

region_fre_10_F<-region_fre_10
region_fre_10_M<-region_fre_10
region_fre_10_all<-read.table("D:\\aging\\RESULT\\2-3-CAS\\B_247genes.txt" , sep = "\t", header = T,stringsAsFactors = F)


  
write.table(region_fre_10_F,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\A_VENN_F.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
write.table(region_fre_10_M,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\A_VENN_M.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
write.table(region_fre_10_all,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\A_VENN_ALL.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

library(ggvenn)
dat <- list( Females =as.character(region_fre_10_F$Var1)  ,  Males =as.character(region_fre_10_M$Var1) , CAS_gene = region_fre_10_all$Var1)
ggvenn(dat,show_percentage = F,
       stroke_color = "BLACK",
       stroke_size = 0.5,
       fill_color = c("#D2966C","#759479","#42618D"),
       set_name_color =c("#D2966C","#759479","#42618D"), 
       set_name_size = 15,text_size=6)





region_DEG1_out<-region_fre_10
col_names<-c()
for (i in 1:16) {
  region_DEG<-read.table(a[i],header=T,sep = "\t", quote = "")
  region_DEG1<-as.data.frame(merge(region_fre_10,region_DEG,by.x = "Var1",by.y = "gene")[,c(3,6,9)]) 
  col_names<-c(col_names,strsplit(as.character(a[i]),"_")[[1]][1])
  region_DEG1_out<-cbind(region_DEG1_out,region_DEG1)
}

hk_gene1<-region_DEG1_out[,-2]
hk_gene1_sum<-apply(hk_gene1[,-1], 1, sum)
up_gene<-hk_gene1[which(hk_gene1_sum>=0),]

down_gene<-hk_gene1[which(hk_gene1_sum<0),]

up_gene_symbol <- bitr(geneID = up_gene[,1],  #感兴趣的基因集
                       fromType="SYMBOL",   #输入ID的类型
                       toType= "ENTREZID",   #输出ID的类型，可为多个
                       OrgDb="org.Hs.eg.db")  #物种注释数据库



up_gene <- up_gene_symbol[,2]
up_BP <- enrichGO(gene = up_gene,  #基因列表(转换的ID)
                  keyType = "ENTREZID",  #指定的基因ID类型，默认为ENTREZID
                  OrgDb=org.Hs.eg.db,  #物种对应的org包
                  ont = "BP",   #CC细胞组件，MF分子功能，BP生物学过程
                  #pvalueCutoff = 0.01,  #p值阈值
                  pAdjustMethod = "fdr",  #多重假设检验校正方式
                  minGSSize = 1,   #注释的最小基因集，默认为10
                  maxGSSize = 500,  #注释的最大基因集，默认为500
                  #qvalueCutoff = 0.005,  #q值阈值
                  readable = TRUE)  #基因ID转换为基因名




down_gene_symbol <- bitr(geneID = down_gene[,1],  #感兴趣的基因集
                         fromType="SYMBOL",   #输入ID的类型
                         toType= "ENTREZID",   #输出ID的类型，可为多个
                         OrgDb="org.Hs.eg.db")  #物种注释数据库



down_gene <- down_gene_symbol[,2]
down_BP <- enrichGO(gene = down_gene,  #基因列表(转换的ID)
                    keyType = "ENTREZID",  #指定的基因ID类型，默认为ENTREZID
                    OrgDb=org.Hs.eg.db,  #物种对应的org包
                    ont = "BP",   #CC细胞组件，MF分子功能，BP生物学过程
                    #pvalueCutoff = 0.01,  #p值阈值
                    pAdjustMethod = "fdr",  #多重假设检验校正方式
                    minGSSize = 1,   #注释的最小基因集，默认为10
                    maxGSSize = 500,  #注释的最大基因集，默认为500
                    #qvalueCutoff = 0.005,  #q值阈值
                    readable = TRUE)  #基因ID转换为基因名
down_BP1<-down_BP@result
down_BP2<-down_BP1[down_BP1$qvalue <0.0001,]

up_BP1<-up_BP@result
#up_BP2<-up_BP1[up_BP1$pvalue<0.0001,]
up_BP2<-up_BP1[up_BP1$pvalue<0.0009,]


BP1<-rbind(as.data.frame(down_BP2),as.data.frame(up_BP2))  
BP<-BP1
write.table(BP,"D:\\aging\\RESULT\\2-3-CAS\\C_BP_F.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 #axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))



p4 <- ggplot() +
  geom_point(data = BP,
             aes(x = -log10(pvalue),
                 y =reorder(Description ,-log10(pvalue)) ,
                 size = Count,
                 color = -log10(pvalue))) +
  scale_size_continuous(range=c(2,5)) + #调整气泡大小范围(默认尺寸部分过小)
  #scale_y_continuous(expand = c(0,0.1),limits = c(0,52)) +
  #scale_x_continuous(limits = c(0.57,2.5)) +
  scale_colour_distiller(palette = "Reds", direction = 1) + #更改配色
  labs(x = "-log10(pvalue)",
       y = "") +
  theme_bw() +
  mytheme

p4

