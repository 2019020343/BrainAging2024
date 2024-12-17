rm(list=ls())
library(ggplot2)
library(ggalluvial)
library(tidyverse)
library(ggbreak)
library(dplyr)
################ Fig3 A #############
setwd("D:\\aging\\data\\geo\\region_DEG\\allgene")
a<-list.files()[-c(1,6,10)]

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

region_fre1<-as.data.frame(table(region_fre[,2]))

ggplot(region_fre1, aes(x = Var1, y=as.numeric(Freq))) +
  #scale_fill_manual(values = c()) +
  #scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  labs(x="the number of regions in which a DEG was detected",y="# DEGs")+
  geom_col(width = 0.6, color = NA, size = 0.5)





#write.table(region_fre1,"D:\\aging\\RESULT\\3CAS\\A_region_fre1.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
################### Fig3 B ##############

library(pheatmap)
region_fre_10<-region_fre[region_fre$Freq>=10,]
write.table(region_fre_10,"D:\\aging\\RESULT\\9online tool\\4\\region_fre_10.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


region_DEG1_out<-region_fre_10
col_names<-c()

region_DEG2_out<-c()
for (i in 1:16) {
  region_DEG<-read.table(a[i],header=T,sep = "\t", quote = "")
  region_DEG1<-as.data.frame(merge(region_fre_10,region_DEG,by.x = "Var1",by.y = "gene")[,c(3,6,9)]) 
  col_names<-c(col_names,strsplit(as.character(a[i]),"_")[[1]][1])
  region_DEG1_out<-cbind(region_DEG1_out,region_DEG1)
  
  
  
#  region_DEG2<-cbind(as.data.frame(merge(region_fre_10,region_DEG,by.x = "Var1",by.y = "gene")[,c(1,3,6,9)]) , strsplit(as.character(a[i]),"_")[[1]][1]) 
#  region_DEG2_out<-rbind(region_DEG2_out, region_DEG2)
  
  
}

#write.table(region_DEG2_out,"D:\\aging\\RESULT\\9online tool\\2\\247gene_heatmap_table.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)




hk_gene1<-region_DEG1_out[,-2]
hk_gene1_sum<-apply(hk_gene1[,-1], 1, sum)

quantile_value <- quantile(hk_gene1_sum, probs = 0.95)

quantile_value1 <- quantile(hk_gene1_sum, probs = 0.05)

hk_gene<-hk_gene1[c(which(hk_gene1_sum>=quantile_value),which(hk_gene1_sum<=quantile_value1)),]



rownames(hk_gene)<-hk_gene[,1]
hk_gene<-hk_gene[,-1]

#hk_gene<-scale(hk_gene)
# 根据指定的样本顺序对数据进行重新排列
paletteLength <- 20
myColor <- colorRampPalette(c("#154182", "#FFFFFF", "#99131D"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(hk_gene), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(hk_gene)/paletteLength, max(hk_gene), length.out=floor(paletteLength/2)))
annotation_col<-as.data.frame( rep(col_names,each=3))
colnames(annotation_col)<-"rengion"
row.names(annotation_col) <- colnames(hk_gene)
ann_color<-list(rengion=c("A1C"="#12803b","AMY"="#ec95b3","CBC"="#f18e25","DFC"="#f5c51e",
                          "HIP"="#211d1e","IPC"="#57a8d7","ITC"="#ac536a","M1C"="#735d96",
                          "MD" ="#cacbd0","MFC"="#f1a7a4","OFC"="#cf1223","S1C"="#f4da9a",
                          "STC"="#f8bf89","STR"="#136cb6","V1C"="#c5b5d1","VFC"="#5A7EB3"))

pheatmap(hk_gene,
         #scale = "row",
         annotation_col = annotation_col, 
         #annotation_row = annotation_row,
         cluster_rows = T,
         cluster_cols = F,
         #cellwidth = 20, cellheight = 15,
         annotation_colors = ann_color,
         color = myColor,
         breaks=myBreaks, 
         #show_colnames = F,
         labels_col = c("","A1C","","","AMY","","","CBC","","","DFC","","",
                        "HIP","","","IPC","","","ITC","","","M1C","","",
                        "MD" ,"","","MFC","","","OFC","","","S1C","","",
                        "STC","","","STR","","","V1C","","","VFC",""),
         angle_col = "45",
         border = F,
         annotation_legend = FALSE,gaps_col = c(3,6,9,12,15,18,21,24,27,30,33,36,39,42,45))

colnames(hk_gene)<-c("A1C_1","A1C_2","A1C_3",  "AMY_1","AMY_2","AMY_3","CBC_1","CBC_2","CBC_3",  "DFC_1","DFC_2","DFC_3",
                     "HIP_1","HIP_2","HIP_3",  "IPC_1","IPC_2","IPC_3","ITC_1","ITC_2","ITC_3",  "M1C_1","M1C_2","M1C_3",
                     "MD_1","MD_2","MD_3",  "MFC_1","MFC_2","MFC_3","OFC_1","OFC_2","OFC_3",  "S1C_1","S1C_2","S1C_3",
                     "STC_1","STC_2","STC_3",  "STR_1","STR_2","STR_3",
                     "V1C_1","V1C_2","V1C_3",  "VFC_1","VFC_2","VFC_3")






write.table(hk_gene,"D:\\aging\\RESULT\\2-3-CAS\\B_HEATMAP.txt",col.names = T, row.names = T,sep = "\t" ,append = FALSE, quote = F)
#write.table(hk_gene1,"D:\\aging\\RESULT\\3CAS\\B_247genes.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

################### Fig3 C ##############


library(clusterProfiler)
library(org.Hs.eg.db)
library(topGO)

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

BP1<-rbind(as.data.frame(down_BP2),as.data.frame(up_BP))  
BP<-BP1


mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))



p3 <- ggplot() +
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

p3
#write.table(BP,"D:\\aging\\RESULT\\3CAS\\C_BP.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


library(tidyverse)
# devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(ggplot2)
#install.packages("cols4all")
library(cols4all)
#BiocManager::install("dittoSeq")
library(dittoSeq)

sankey_data<-c()
for (j in 1:nrow(BP)) {
  out<-cbind(BP$Description[j],unlist( strsplit(BP$geneID[j],"/")) )
  sankey_data<-rbind(sankey_data,out)
}

sankey<-as.data.frame(sankey_data)
colnames(sankey)<-c("pathNames","metamolites")
df <- sankey %>%
  make_long(metamolites, pathNames)



#指定绘图顺序（转换为因子）：
df$node <- factor(df$node,levels = c(sankey$pathNames %>% unique()%>% rev(),
                                     sankey$metamolites %>% unique() %>% rev()))
#自定义配色：
c4a_gui()
mycol <- c(rep("#C07A7F",3),rep("#98AAC5",33) ,  colorRampPalette(c("#5D779B", "#FAF3BE", "#C9301D"))(nrow(BP))  )

#绘图：
ggplot(df,aes(y = car, axis1 = inccat, axis2 =ed,axis3 = marital)) +
  geom_alluvium(aes(fill = gender)) +
  geom_stratum(width = 1/6, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  ggtitle("收入和购买汽车关系")


p4 <- ggplot(df, aes(x = x,
                     next_x = next_x,
                     node = node,
                     next_node = next_node,
                     fill = node,
                     label = node)) +
  geom_sankey(flow.alpha = 0.5,
              flow.fill = 'grey',
              flow.color = 'grey80', #条带描边色
              node.fill = mycol, #节点填充色
              smooth = 8,
              width = 0.08) +
  geom_sankey_text(size = 2,
                   color = "black")+
  theme_void() +
  theme(legend.position = 'none')
p4

#write.table(df,"D:\\aging\\RESULT\\3CAS\\C_sankey.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

p5 <- p4 + theme(plot.margin = unit(c(0,5,0,0),units="cm"))
p6 <-ggdraw() + draw_plot(p5) + draw_plot(p3, scale = 0.5, x = 0.62, y=-0.21, width=0.48, height=1.37)





######################### CAS ###################
rm(list=ls())

require(devtools)
#install_github("YosefLab/VISION",force = TRUE)

# Load VISION
library(VISION)
hk_gene1<-read.table("D:\\aging\\RESULT\\2-3-CAS\\B_247genes.txt",header=T,sep = "\t", quote = "")
sample_group<-read.table("D:\\aging\\data\\geo\\sample_group.txt",header=T,sep = "\t", quote = "")

hk_gene1_sum<-apply(hk_gene1[,-1], 1, sum)
up<-cbind(hk_gene1[hk_gene1_sum>0,1],-1 ) 
down<-cbind(hk_gene1[hk_gene1_sum<0,1],1 ) 
sigData1<-t(rbind(up,down ))

sigData2<-as.numeric(sigData1[2,])
names(sigData2)<-sigData1[1,]

sig <- createGeneSignature(name = "Common Aging Score", sigData = sigData2)
mySignatures <- c(sig)


# Read in expression counts (Genes X Cells)
counts <- read.table("D:\\aging\\data\\geo\\expmean_gene_NOlog.txt",
                     header = TRUE,
                     sep = '\t',
                     row.names = 1)



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

colnames(sigScores1)<-c("Samples","CAS")

write.table(sigScores1,"D:\\aging\\RESULT\\9online tool\\4\\CAS_out.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



sigScores1<-cbind(names(sigScores),sigScores)
sigScores_group<-merge(sigScores1,sample_group,by.x="V1",by.y= "samples" )
sigScores_group<-sigScores_group[sigScores_group$age>=20,]


sigScores_group$group1=70
sigScores_group$group1[sigScores_group$age < 70]="60"
sigScores_group$group1[sigScores_group$age < 60]="50"
sigScores_group$group1[sigScores_group$age < 50]="40"
sigScores_group$group1[sigScores_group$age < 40]="30"
sigScores_group$group1[sigScores_group$age < 30]="20"



sigScores_group$group <- factor(sigScores_group$group,levels = c('childhood','Adolescence','Young','Middle','Late'))

write.table(sigScores_group,"D:\\aging\\RESULT\\2-3-CAS\\sigScores_group.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


P3<-ggplot(sigScores_group,aes(x=group1,y=as.numeric(sigScores) ,fill=group1))+ 
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(color="azure4",outlier.colour="red",
               outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
  )+  
  
  facet_wrap(~region)+
  #stat_compare_means( method="anova")+
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
sigScores_group<-read.table("D:\\aging\\RESULT\\2-3-CAS\\sigScores_group.txt",header=T,sep = "\t", quote = "")

sigScores_group2<-sigScores_group[sigScores_group$Sex=="Sex: F",]
write.table(sigScores_group2,"D:\\aging\\RESULT\\2-3-CAS\\sigScores_group_F.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

sigScores_group3<-sigScores_group[sigScores_group$Sex=="Sex: M",]
write.table(sigScores_group3,"D:\\aging\\RESULT\\2-3-CAS\\sigScores_group_M.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



P3<-ggplot(sigScores_group3,aes(x=group1,y=as.numeric(sigScores) ,fill=Sex))+ 
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(color="azure4",outlier.colour="red",
               outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
  )+  
  
  facet_wrap(~region)+
  #stat_compare_means( method="anova")+
  #scale_fill_material_d()+
  theme_classic()+  
  #ylim(c(0.5,1))+
  scale_fill_manual(values = c("#77977A"))+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Group", y="CAS")+ 
  theme(axis.title.x =element_text(size=12,face = "bold"), axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size = 12,face = "bold"))+
  #stat_summary(fun.y = "mean", geom = "point", size = 0.5) +
  stat_summary(fun.y = "mean", geom = "line", 
               aes(group = 1), 
               size = 0.5)

P3



write.table(sigScores_group2,"D:\\aging\\RESULT\\2-3-CAS\\sigScores_group2.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



P3<-ggplot(sigScores_group2,aes(x=group1,y=as.numeric(sigScores) ,fill=Sex))+ 
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(color="azure4",outlier.colour="red",
               outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
  )+  
  
  facet_wrap(~region)+
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
               aes(group = Sex), 
               size = 0.5)

P3

###############xcell ###############
library(xCell)
expmean_gene<-read.table("D:\\aging\\data\\geo\\expmean_gene.txt",header=T,row.names = 1, sep = "\t", quote = "")
sigScores_group<-read.table("D:\\aging\\RESULT\\3CAS\\sigScores_group.txt",header=T,sep = "\t", quote = "")
expmean_gene1<-expmean_gene[,sapply(sigScores_group$V1,function(x){which(colnames(expmean_gene)==x)})]
scores <-  xCellAnalysis(expmean_gene1,rnaseq = T)
xCell_scores<-as.data.frame(t(scores))
xCell_scores1<-cbind(rownames(xCell_scores),xCell_scores)
xCell_group<-merge(xCell_scores1,sigScores_group,by.x = "rownames(xCell_scores)",by.y ="V1" )
install.packages("ggExtra")
library(tidyverse)       
library(ggsci)           
library(ggExtra)         
library(ggpmisc)        
library(palmerpenguins) 
library(ggpubr)
p <- ggplot(xCell_group, aes(sigScores, xCell_group[,11],color=region)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                               "#f8bf89","#136cb6","#c5b5d1","#5A7EB3")) +
  labs(x = "CAS", y = "cd4_tem") +  # 设置坐标轴标签
  theme_classic()  +
  
  stat_cor(aes(color = region), label.x = 0)
p

p1 <- ggplot(xCell_group, aes(sigScores, xCell_group[,15],color=region)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                               "#f8bf89","#136cb6","#c5b5d1","#5A7EB3")) +
  labs(x = "CAS", y = "cd8_tem") +  # 设置坐标轴标签
  theme_classic()  +
  
  stat_cor(aes(color = region), label.x = 0)
p1

for (i in 2:68) {
  ggplot(xCell_group, aes(sigScores, xCell_group[,i],color=region)) +
    geom_point() +  # 添加散点图层，点的大小表示体重
    scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                                 "#211d1e","#57a8d7","#ac536a","#735d96",
                                 "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                                 "#f8bf89","#136cb6","#c5b5d1","#5A7EB3")) +
    labs(x = "CAS", y = colnames(xCell_group)[i]) +  # 设置坐标轴标签
    theme_classic()  +
    stat_cor(aes(color = region), geom = "text",label.x = 0)+
    theme(axis.title = element_text(family = "sans", 
                                  face = "bold"))
        
         
  ggsave(paste("D:\\aging\\RESULT\\3CAS\\xcell\\",colnames(xCell_group)[i],".pdf",sep="") )
  
}


region_list<-unique(xCell_group$region)
values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
         "#211d1e","#57a8d7","#ac536a","#735d96",
         "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
         "#f8bf89","#136cb6","#c5b5d1","#5A7EB3")
for (j in 1:16) {
  
  cd4_tem_V1C<-xCell_group[xCell_group$region==region_list[j],]
  ggplot(cd4_tem_V1C, aes(sigScores, cd4_tem_V1C[,11],color=region)) +
    geom_point() +  # 添加散点图层，点的大小表示体重
    scale_colour_manual(values=c("#997DAD")) +
    labs(x = "CAS", y = "cd4_tem") +  # 设置坐标轴标签
    theme_classic()  +
    stat_cor(aes(color = region))+ 
    geom_smooth(method="lm", se=FALSE)
    
  ggsave(paste("D:\\aging\\RESULT\\3CAS\\xcell\\cd4_tem_",region_list[j],".pdf",sep="") )
  
}






################ regression   #########################
#install.packages("caret", dependencies = c("Depends", "Suggests"))


sigScores_group<-read.table("D:\\aging\\RESULT\\2-3-CAS\\sigScores_group.txt",header=T,sep = "\t", quote = "")
# 逻辑回归
library(PerformanceAnalytics)#加载包
head(mtcars)#查看前5行，前10列
chart.Correlation(sigScores_group[,c(2,4)], histogram=TRUE, pch=19)


lm1 <- lm(age~sigScores, data = sigScores_group)
su1 = summary(lm1)

library('ggpmisc')
f1 <- y ~ x #定义回归方程


library('ggplot2')

ggplot(sigScores_group, aes(x=age, y=as.numeric(sigScores) ,color = region)) + 
  theme_classic()+
  
  ggtitle('loess fit')+
  
  geom_smooth(method = c('loess'), se=FALSE,  span=0.8) + ##不规则拟合

  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                               "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))




ggplot(sigScores_group, aes(x=age, y=as.numeric(sigScores) ,color = region)) + 
  theme_classic()+
  
  ggtitle('linear fit')+
  
  geom_smooth(method = "lm", formula = y ~ x,size = 2,se = F)+ ##二项式拟合
  
  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                               "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))



sigScores_group2<-sigScores_group[sigScores_group$Sex=="Sex: F", ]


ggplot(sigScores_group2, aes(x=age, y=sigScores,color = region)) + 
  theme_classic()+
  
  ggtitle('loess fit')+
  
  geom_smooth(method = c('loess'), se=FALSE,  span=0.8) + ##不规则拟合
  
  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                               "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))




ggplot(sigScores_group2, aes(x=age, y=sigScores,color = region)) + 
  theme_classic()+
  
  ggtitle('linear fit')+
  
  geom_smooth(method = "lm", formula = y ~ x,size = 2,se = F)+ ##二项式拟合
  
  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                               "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))

############# Sex fit  ######################

sigScores_group3<-sigScores_group[sigScores_group$Sex=="Sex: M", ]


ggplot(sigScores_group3, aes(x=age, y=sigScores,color = region)) + 
  theme_classic()+
  
  ggtitle('loess fit')+
  
  geom_smooth(method = c('loess'), se=FALSE,  span=0.8) + ##不规则拟合
  
  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                               "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))




ggplot(sigScores_group3, aes(x=age, y=sigScores,color = region)) + 
  theme_classic()+
  
  ggtitle('linear fit')+
  
  geom_smooth(method = "lm", formula = y ~ x,size = 2,se = F)+ ##二项式拟合
  
  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                               "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))
######################### slope #########################
library(ggplot2)
library(ggpubr)
library(ggpmisc)

sigScores_group$sigScores<-as.numeric(sigScores_group$sigScores)
species_lst <- unique(sigScores_group$region)
p_lst <- list()
data_lst <- list()
formula <- y ~ x
for (i in species_lst) {
  extract_data <- sigScores_group[which(sigScores_group$region==i),]
  p <- ggplot(extract_data,aes(age,sigScores)) + 
    geom_point() + 
    geom_smooth( method = "loess") + 
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

write.table(extract_ggplotdata,"D:\\aging\\RESULT\\3CAS\\slope_loess.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


cols<-colorRampPalette(c("#5D779B", "#FAF3BE", "#C9301D"))(200)
extract_ggplotdata<-read.table("D:\\aging\\RESULT\\3CAS\\slope.txt",header=T,sep = "\t", quote = "")
ggplot(extract_ggplotdata, aes(x =reorder(region,-slope) , y=slope*1000, fill=slope)) +
  scale_fill_gradientn(limits=c(0, 0.01), colours = cols)+
  #scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  labs(x="",y="CAS slope (x 10 -3)")+
  theme(axis.text.x = element_text(angle = 45))+#字体大小
  geom_col(width = 0.6, color = NA, size = 0.5)

############## SEX ###############
sigScores_group<-sigScores_group2
#sigScores_group<-sigScores_group3
species_lst <- unique(sigScores_group$region)
p_lst <- list()
data_lst <- list()
formula <- y ~ x
for (i in species_lst) {
  extract_data <- sigScores_group[which(sigScores_group$region==i),]
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

#write.table(extract_ggplotdata,"D:\\aging\\RESULT\\3CAS\\slope_M.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

cols<-colorRampPalette(c("#5D779B", "#FAF3BE", "#C9301D"))(200)
extract_ggplotdata1<-cbind(read.table("D:\\aging\\RESULT\\2-3-CAS\\slope_F.txt",header=T,sep = "\t", quote = "")[,c(1,5)],"F")
extract_ggplotdata2<-cbind(read.table("D:\\aging\\RESULT\\2-3-CAS\\slope_M.txt",header=T,sep = "\t", quote = "")[,c(1,5)],"M")
colnames(extract_ggplotdata1)<-c("region","slope","sex")
colnames(extract_ggplotdata2)<-c("region","slope","sex")
  
  
extract_ggplotdata_sex<-rbind(extract_ggplotdata1,extract_ggplotdata2)

ggplot(extract_ggplotdata_sex, aes(x =reorder(region,-slope) , y=slope*1000, fill=sex)) +
  geom_bar(position=position_dodge(), stat="identity")+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+

  scale_fill_manual(values = c("#D8996D","#77977A"))+ 
  #scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  labs(x="",y="CAS slope (x 10 -3)")+　
  theme(axis.text.x = element_text(angle = 45))

  t.test(extract_ggplotdata_sex$slope[extract_ggplotdata_sex$sex=="F"],extract_ggplotdata_sex$slope[extract_ggplotdata_sex$sex=="M"])
  0.00480/0.00444
################# loess slope ################
  install.packages("segmented")
 library(segmented)
  library(ggplot2)
  sigScores_group<-read.table("D:\\aging\\RESULT\\3CAS\\sigScores_group.txt",header=T,sep = "\t", quote = "")
  region_list<-unique( sigScores_group$region) 
  slope_list<-c()
  slope_list_F<-c()
  slope_list_M<-c()
  
  
  for (i in 1:16) {
    sigScores_group_regioni<-sigScores_group[which( sigScores_group$region==region_list[i] ), ]
    x=sigScores_group_regioni$age
    y=sigScores_group_regioni$sigScores
    
    #Obtain points for the smooth curve
    temp = loess.smooth(x, y, evaluation = 50) #Use higher evaluation for more points
    #Obtain slope of the smooth curve
    slopes = diff(temp$y)/diff(temp$x)
    loess_slope<-cbind(temp$x[-50],slopes,region_list[i])
    colnames(loess_slope)<-c("x","slope","region")
    slope_list<-rbind(slope_list,loess_slope)
    
    sigScores_group_regioni_F<-sigScores_group_regioni[which( sigScores_group_regioni$Sex=="Sex: F"), ]
    x_F=sigScores_group_regioni_F$age
    y_F=sigScores_group_regioni_F$sigScores
    #Obtain points for the smooth curve
    temp_F = loess.smooth(x_F, y_F, evaluation = 50) #Use higher evaluation for more points
    #Obtain slope of the smooth curve
    slopes_F = diff(temp_F$y)/diff(temp_F$x)
    loess_slope_F<-cbind(temp_F$x[-50],slopes_F,region_list[i],"F")
    colnames(loess_slope_F)<-c("x","slope","region","Sex")
    slope_list_F<-rbind(slope_list_F,loess_slope_F)
    
    sigScores_group_regioni_M<-sigScores_group_regioni[which( sigScores_group_regioni$Sex=="Sex: M"), ]
    x_M=sigScores_group_regioni_M$age
    y_M=sigScores_group_regioni_M$sigScores
    #Obtain points Mor the smooth curve
    temp_M = loess.smooth(x_M, y_M, evaluation = 50) #Use higher evaluation for more points
    #Obtain slope of the smooth curve
    slopes_M = diff(temp_M$y)/diff(temp_M$x)
    loess_slope_M<-cbind(temp_M$x[-50],slopes_M,region_list[i],"M")
    colnames(loess_slope_M)<-c("x","slope","region","Sex")
    slope_list_M<-rbind(slope_list_M,loess_slope_M)

    
    
  }

  ggplot(slope_list, aes(x =as.numeric(x), y=as.numeric(slope)*1000 ,color=region)) +
    geom_line(position=position_dodge(), stat="identity")+
    
    theme_classic()+
    labs(x="",y="CAS slope (x 10 -3)")+　
    geom_text(aes(label = region), hjust = -0.2, vjust = 0.5)+
    
    scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                                 "#211d1e","#57a8d7","#ac536a","#735d96",
                                 "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                                 "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))
  
  ggplot(slope_list_F, aes(x =as.numeric(x), y=as.numeric(slope)*1000 ,color=region)) +
    geom_line(position=position_dodge(), stat="identity")+
    geom_text(aes(label = region), hjust = -0.2, vjust = 0.5)+
    
    theme_classic()+
    labs(x="",y="CAS slope (x 10 -3)")+　

    scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                                 "#211d1e","#57a8d7","#ac536a","#735d96",
                                 "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                                 "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))

  ggplot(slope_list_M, aes(x =as.numeric(x), y=as.numeric(slope)*1000 ,color=region)) +
    geom_line(position=position_dodge(), stat="identity")+
    geom_text(aes(label = region), hjust = -0.2, vjust = 0.5)+
    
    theme_classic()+
    labs(x="",y="CAS slope (x 10 -3)")+　
    
    scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                                 "#211d1e","#57a8d7","#ac536a","#735d96",
                                 "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                                 "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))
  
################ brain view ###############################
  options(repos = c(
    ggseg = 'https://ggseg.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))
  
  # Install some packages
  #install.packages('ggsegBrodmann')
  library(ggseg)
  #> Warning: package 'ggseg' was built under R version 4.1.1
  #> Loading required package: ggplot2
  library(ggseg3d)
  library(ggsegBrodmann)
  slope_BA<-read.table("D:\\aging\\RESULT\\3CAS\\slope_BA.txt",header=T,sep = "\t", quote = "",fill = T)
  
  brodmann1<-as.data.frame(brodmann)
  slope_BA1<-merge(slope_BA,brodmann1,by.x = "BA",by.y = "region")
  cols<-colorRampPalette(c("#5D779B", "#FAF3BE", "#C9301D"))(200)
  slope_BA1 %>%
    ggplot() +
    geom_brain(atlas = brodmann, 
               position = position_brain(hemi ~ side),
               aes(fill = slope)) +
    theme_classic()+
    scale_fill_gradientn(limits=c(0, 0.01), colours = cols)
########### sex ############
  slope_BA_F<-read.table("D:\\aging\\RESULT\\3CAS\\slope_F_BA.txt",header=T,sep = "\t", quote = "",fill = T)
  
  brodmann1<-as.data.frame(brodmann)
  slope_BA1_F<-merge(slope_BA_F,brodmann1,by.x = "BA",by.y = "region")
  cols<-colorRampPalette(c("#5D779B", "#FAF3BE", "#C9301D"))(200)
  slope_BA1_F %>%
    ggplot() +
    geom_brain(atlas = brodmann, 
               position = position_brain(hemi ~ side),
               aes(fill = slope)) +
    theme_classic()+
    scale_fill_gradientn(limits=c(0, 0.01), colours = cols)
  
  
  slope_BA_M<-read.table("D:\\aging\\RESULT\\3CAS\\slope_M_BA.txt",header=T,sep = "\t", quote = "",fill = T)
  
  brodmann1<-as.data.frame(brodmann)
  slope_BA1_M<-merge(slope_BA_M,brodmann1,by.x = "BA",by.y = "region")
  cols<-colorRampPalette(c("#5D779B", "#FAF3BE", "#C9301D"))(200)
  slope_BA1_M %>%
    ggplot() +
    geom_brain(atlas = brodmann, 
               position = position_brain(hemi ~ side),
               aes(fill = slope)) +
    theme_classic()+
    scale_fill_gradientn(limits=c(0, 0.01), colours = cols)
  
  slope_BA1_M %>%
    ggplot() +
    geom_brain(atlas = brodmann, 
               position = position_brain(hemi ~ side),
               aes(fill = region)) +
    #geom_text(aes(label = region))+
    theme_classic()
  #  scale_fill_gradientn(limits=c(0, 0.01), colours = cols)
  
  
  slope_BA1_AIC<-slope_BA1[slope_BA1$region=="V1C",]
  
  slope_BA1_AIC %>%
    ggplot() +
    geom_brain(atlas = brodmann, 
               position = position_brain(hemi ~ side),
               aes(fill = slope)) +
    #geom_text(aes(label = region))+
    theme_classic()
  #  scale_fill_gradientn(limits=c(0, 0.01), colours = cols)