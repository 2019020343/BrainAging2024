rm(list=ls())
library("Seurat")
library(gplots)
library(ggplot2)
#install.packages("tidydr")
#devtools::install_github('junjunlab/scRNAtoolVis')
library(clustree)
library(scRNAtoolVis)
library(tidydr)
#################################Fig.2-A-C  #########################
sample_group <-read.table("D:\\aging\\data\\geo\\sample_group.txt",header=T,sep = "\t", quote = "",fill = T,row.names = 1)
sample_group$group[sample_group$age<12]="Childhood"
sample_group$group[sample_group$age>=12 & sample_group$age<20]="Adolescence"
sample_group$group[sample_group$age>=20 & sample_group$age<40]="Young"
sample_group$group[sample_group$age>=40 & sample_group$age<60]="Middle"
sample_group$group[sample_group$age>60]="Late"

sample_group$region1<-"NCX"
sample_group$region1[sample_group$region=="region: AMY"]="AMY"
sample_group$region1[sample_group$region=="region: MD"]="MD"
sample_group$region1[sample_group$region=="region: HIP"]="HIP"
sample_group$region1[sample_group$region=="region: STR"]="STR"
sample_group$region1[sample_group$region=="region: CBC"]="CBC"

sample_group$region2[grepl("L$", sample_group$Sample_description) ]<-"Left"
sample_group$region2[grepl("R$", sample_group$Sample_description) ]<-"Right"


expmean_gene<-read.table("D:\\aging\\data\\geo\\expmean_gene.txt",header=T,sep = "\t", quote = "",row.names = 1)
index<-sapply(rownames(sample_group),function(x){which(colnames(expmean_gene)==x)})
expmean_gene<-expmean_gene[,index]

expmean_gene1<-2^expmean_gene

#write.table(expmean_gene1,"D:\\aging\\data\\geo\\expmean_gene_NOlog.txt",col.names = T, row.names = T,sep = "\t" ,append = FALSE, quote = F)

scRNA = CreateSeuratObject(counts=expmean_gene1,
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

#dims = 1:10 即选取前10个主成分来分类细胞。
#查看前5个细胞的分类ID

#DimPlot(pbmc1, reduction = "umap")
p1=DimPlot(pbmc1, reduction = "umap",  group.by = 'group',label = T,cols = c("#e8a361","#bb322b","#1f2a6b","#446799","#a3ced6"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))

p3=DimPlot(pbmc1, reduction = "umap",  group.by = 'region2',label = T,cols = c("#f18e25","#12803b"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))
p4=DimPlot(pbmc1, reduction = "umap",  group.by = 'region1',label = T,cols = c("#ec95b3","#f18e25","#211d1e","#cacbd0","#cf1223","#136cb6"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))
p2=DimPlot(pbmc1, reduction = "umap",  group.by = 'region',label = F,cols = c("#12803b","#ec95b3","#f18e25","#f5c51e",
                                                                              "#211d1e","#57a8d7","#ac536a","#735d96",
                                                                              "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                                                                              "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))


p5=DimPlot(pbmc1, reduction = "umap",  group.by = 'Sex',label = F,cols = c("#cf1223","#136cb6"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))
p6=DimPlot(pbmc1, reduction = "umap",  group.by = 'age',label = T,label.color = "black",cols = c("#CB6D5D","#D56B19","#FAEACE","#F8B9AB",
                                                                                                 "#85AD9B","#747D90","#D93B45","#E32B2D",
                                                                                                 "#AFBED2","#8DBCCA","#827964","#577465",
                                                                                                 "#ED6767","#EA4840","#F9A27E","#FB9C42",
                                                                                                 "#C89B60","#504646","#b288ab","#d8d8dc",
                                                                                                 "#C7BA8E","#BE8B89","#C2B5CD","#C3CAAC"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))
p=p1+p2+p3+p4+p5+p6
#  ggsave(paste("D:\\aging\\RESULT\\Spatiotemporal quantification\\UMAP\\umap\\resolution_",j,".png",sep = ""),width = 16,height = 8)



################### single region #########################
NCX_group<-sample_group[sample_group$region1=="NCX",]

NCX_index<-sapply(rownames(NCX_group),function(x){which(colnames(expmean_gene)==x)})
NCX_expmean_gene<-expmean_gene[,NCX_index]

expmean_gene1<-2^NCX_expmean_gene


scRNA = CreateSeuratObject(counts=expmean_gene1,
                           meta.data = sample_group,
                           min.cells = 3, 
                           min.features = 50)
head(scRNA@meta.data, 5)
#nCount_RNA与nFeature_RNA的相关性
pbmc <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.8)
#dims = 1:10 即选取前10个主成分来分类细胞。
#查看前5个细胞的分类ID
head(Idents(pbmc), 5)
pbmc1 <- RunUMAP(pbmc, dims = 1:40, verbose = T)
DimPlot(pbmc1, reduction = "umap")
DimPlot(pbmc1, reduction = "umap",  group.by = 'group',label = T,cols = c("#e8a361","#bb322b","#1f2a6b","#446799","#a3ced6"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))
DimPlot(pbmc1, reduction = "umap",  group.by = 'region2',label = T,cols = c("#f18e25","#12803b"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))



DimPlot(pbmc1, reduction = "umap",  group.by = 'region1',label = T,cols = c("#cf1223"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))
DimPlot(pbmc1, reduction = "umap",  group.by = 'region',label = F,cols = c("#12803b","#ec95b3","#f18e25","#f5c51e",
                                                                           "#211d1e","#57a8d7","#ac536a","#735d96",
                                                                           "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                                                                           "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))
DimPlot(pbmc1, reduction = "umap",  group.by = 'Sex',label = F,cols = c("#cf1223","#136cb6"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))
DimPlot(pbmc1, reduction = "umap",  group.by = 'age',label = T,cols = c("#CB6D5D","#D56B19","#FAEACE","#F8B9AB",
                                                                        "#85AD9B","#747D90","#D93B45","#E32B2D",
                                                                        "#AFBED2","#8DBCCA","#827964","#577465",
                                                                        "#ED6767","#EA4840","#F9A27E","#FB9C42",
                                                                        "#C89B60","#504646","#b288ab","#d8d8dc",
                                                                        "#C7BA8E","#BE8B89","#C2B5CD","#C3CAAC"))+
  theme_dr(xlength = 0.2,ylength=0.2,arrow=arrow(length = unit(0.2,"inches"),type = "closed")  )+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 3,hjust = 0.03))

#################################Fig.2-D  #########################
library(ggplot2)
library(ggalt)
library(reshape2)
setwd("D:\\aging\\data\\geo\\region_DEG\\allgene\\5")
a<-list.files()

region_number<-c()
r9_data1_out<-c()
for (i in 1:length(a)) {
  region_name<-strsplit(a[i],"_")[[1]][1]
  region <-read.table(a[i],header=T,sep = "\t", quote = "",fill = T)
  Adolescence_data<-region[,1:4]
  Young_data<-region[,c(1,5:7) ] 
  Middle_data<-region[, c(1,8:10)] 
  Late_data<-region[, c(1,11:13)] 
  colnames(Adolescence_data)<-c("gene", "logFC","P.Value","adj.P.Val")
  colnames(Young_data)<-c("gene", "logFC","P.Value","adj.P.Val")
  colnames(Middle_data)<-c("gene", "logFC","P.Value","adj.P.Val")
  colnames(Late_data)<-c("gene", "logFC","P.Value","adj.P.Val")
  r9_data1<-rbind( Adolescence_data,Young_data,Middle_data,Late_data )
  r9_data1$regulation<-"up"
  r9_data1$regulation[r9_data1$logFC <= 0]<-"down"
  r9_data1$Comparison<-rep(c("Adolescence","Young","Middle","Late"),each= nrow(Late_data) )
  r9_data1<-r9_data1[r9_data1$adj.P.Val<0.05,]
  
  
  r9_data1$regoin<-region_name
  r9_data1_out<-rbind(r9_data1_out,r9_data1)

  
 write.table(r9_data1_out,paste("D:\\aging\\RESULT\\9online tool\\1\\1-1\\ALL_gene_table.txt",sep="" )  ,col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

  
  data<-rbind(cbind(sum(region$Adolescence_adj.P.Val <0.05 & region$Adolescence_logFC>0),"up","Adolescence",strsplit(a[i],"_")[[1]][1]),
              cbind(-sum(region$Adolescence_adj.P.Val<0.05 & region$Adolescence_logFC<0),"down","Adolescence",strsplit(a[i],"_")[[1]][1]),
              cbind(sum(region$Young_adj.P.Val<0.05 & region$Young_logFC>0),"up","Young",strsplit(a[i],"_")[[1]][1]),
              cbind(-sum(region$Young_adj.P.Val<0.05 & region$Young_logFC<0),"down","Young",strsplit(a[i],"_")[[1]][1]),
              cbind(sum(region$Middle_adj.P.Val<0.05 & region$Middle_logFC>0),"up","Middle",strsplit(a[i],"_")[[1]][1]),
              cbind(-sum(region$Middle_adj.P.Val<0.05 & region$Middle_logFC<0),"down","Middle",strsplit(a[i],"_")[[1]][1]),
              cbind(sum(region$Late_adj.P.Val<0.05 & region$Late_logFC>0),"up","Late",strsplit(a[i],"_")[[1]][1]),
              cbind(-sum(region$Late_adj.P.Val<0.05 & region$Late_logFC<0),"down","Late",strsplit(a[i],"_")[[1]][1]))
  
  region_number<-rbind(data,region_number)

}
region_number<-as.data.frame(region_number )
colnames(region_number)<-c("numbers","regulation","Comparison","regoin")
region_number2<-rbind(cbind(sum(as.numeric( region_number[region_number$Comparison=="Adolescence" & region_number$regulation=="up",1])),"up","Adolescence") ,
                      cbind(sum(as.numeric( region_number[region_number$Comparison=="Adolescence" & region_number$regulation=="down",1])),"down","Adolescence") ,
                      cbind(sum(as.numeric( region_number[region_number$Comparison=="Young" & region_number$regulation=="up",1])),"up","Young") ,
                      cbind(sum(as.numeric( region_number[region_number$Comparison=="Young" & region_number$regulation=="down",1])),"down","Young") ,
                      cbind(sum(as.numeric( region_number[region_number$Comparison=="Middle" & region_number$regulation=="up",1])),"up","Middle") ,
                      cbind(sum(as.numeric( region_number[region_number$Comparison=="Middle" & region_number$regulation=="down",1])),"down","Middle") ,
                      cbind(sum(as.numeric( region_number[region_number$Comparison=="Late" & region_number$regulation=="up",1])),"up","Late") ,
                      cbind(sum(as.numeric( region_number[region_number$Comparison=="Late" & region_number$regulation=="down",1])),"down","Late")
) 
region_number2<-as.data.frame(region_number2)

colnames(region_number2)<-c("numbers","regulation","Comparison")

#########  up lines #########
region_number_up<-region_number[region_number$regulation=="up",]
region_number_up$v1<-rownames(region_number_up)  
region_number_up1<-region_number_up[,c(5,4,3,1)]
region_number_up2<-dcast(region_number_up1,   regoin ~ Comparison ) ## 长矩阵转短矩阵
rownames(region_number_up2)<-region_number_up2[,1]
region_number_up3<-region_number_up2[,-1]
region_number_up3<-matrix(as.numeric(unlist(region_number_up3) ),ncol = ncol(region_number_up3)) ## 字符型转数值型
rownames(region_number_up3)<-rownames(region_number_up2)
colnames(region_number_up3)<-colnames(region_number_up2)[-1]

p=ggplot(data=region_number_up,
         aes(x=Comparison,y=as.numeric(numbers),
             group=regoin,colour=regoin))+
  geom_point(size=3)+
  labs(x="Comparison", y=" # of DEGs")+
  ylim(c(-100,max(as.numeric(region_number_up$numbers))+50))+
  geom_xspline(spline_shape = -0.5)+
  geom_text(aes(label = regoin), hjust = -0.2, vjust = 0.5)+
  scale_x_discrete(limits=c("Adolescence","Young","Middle","Late"))+
  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                               "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))+
  theme_bw()

p

#########  down lines #########
region_number_down<-region_number[region_number$regulation=="down",]
region_number_down$v1<-rownames(region_number_down)  
region_number_down1<-region_number_down[,c(5,4,3,1)]
region_number_down2<-dcast(region_number_down1,   regoin ~ Comparison )
rownames(region_number_down2)<-region_number_down2[,1]
region_number_down3<-region_number_down2[,-1]
region_number_down3<-matrix(as.numeric(unlist(region_number_down3) ),ncol = ncol(region_number_down3))
rownames(region_number_down3)<-rownames(region_number_down2)
colnames(region_number_down3)<-colnames(region_number_down2)[-1]

p1=ggplot(data=region_number_down,
         aes(x=Comparison,y=as.numeric(numbers),
             group=regoin,colour=regoin))+
  geom_point(size=3)+
  ylim(c(min(as.numeric(region_number_down$numbers))-10,100))+
  labs(x="Comparison", y=" # of DEGs")+
  geom_xspline(spline_shape = -0.5)+
  geom_text(aes(label = regoin), hjust = -0.2, vjust = 0.5)+
  scale_x_discrete(limits=c("Adolescence","Young","Middle","Late"))+
  scale_colour_manual(values=c("#12803b","#ec95b3","#f18e25","#f5c51e",
                               "#211d1e","#57a8d7","#ac536a","#735d96",
                               "#cacbd0","#f1a7a4","#cf1223","#f4da9a",
                               "#f8bf89","#136cb6","#c5b5d1","#5A7EB3"))+
  theme_bw()

p1


F4data<-rbind( cbind(rownames(region_number_up3), region_number_up3,"Up"),cbind(rownames(region_number_down3), region_number_down3,"Down")) 

write.table(F4data,"D:\\aging\\RESULT\\1-Spatiotemporal quantification\\F4data.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



#################################Fig.2-E  #########################
region_number_down$v1<-rownames(region_number_down)  
region_number_down1<-region_number_down[,c(5,4,3,1)]
region_number_down2<-dcast(region_number_down1,   regoin ~ Comparison )


region_number_updown3<-cbind(region_number_down3,region_number_up3)
library(pheatmap)
colnames(region_number_updown3)<-c("Adolescence","Young","Middle","Late","Adolescence1","Young1","Middle1","Late1")
# 定义样本顺序
sample_order <-c("Adolescence","Young","Middle","Late","Adolescence1","Young1","Middle1","Late1")

# 根据指定的样本顺序对数据进行重新排列
region_number_updown3 <- region_number_updown3[,match(sample_order, colnames(region_number_updown3)) ]
paletteLength <- 50
myColor <- colorRampPalette(c("#154182", "#FFFFFF", "#99131D"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(region_number_updown3), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(region_number_updown3)/paletteLength, max(region_number_updown3), length.out=floor(paletteLength/2)))
pheatmap(region_number_updown3, 
         #annotation_col = annotation_col, 
         #annotation_row = annotation_row,
         cluster_rows = T,
         cluster_cols = F,
         cellwidth = 20, cellheight = 15,
         color = myColor,
         breaks=myBreaks, 
         annotation_legend = FALSE, gaps_col = c(4))
#################################Fig.2-F  #########################
setwd("D:\\aging\\data\\geo\\region_DEG\\group\\5")
b<-list.files()
sample_inf <-read.table("D:\\aging\\data\\geo\\sample_group.txt",header=T,sep = "\t", quote = "",fill = T)
expmean_gene<-read.table("D:\\aging\\data\\geo\\expmean_gene.txt",header=T,sep = "\t", quote = "")
expmean_gene1<-as.data.frame(t(expmean_gene))
expmean_gene1<-cbind(rownames(expmean_gene1),expmean_gene1)
colnames(expmean_gene1)<-expmean_gene1[1,]
expmean_gene1<-expmean_gene1[-1,]



re_out_p<-c()
dataout<-c()
for (i in 1:length(b)) {
  group <-read.table(b[i],header=T,sep = "\t", quote = "",fill = T)
  group1 <- merge(group,sample_inf,by.x = "childhood",by.y = "Sample_geo_accession")
  group2 <- merge(group1,expmean_gene1,by.x = "childhood",by.y = "ID_REF")
  pos_count<-0
  neg_count<-0
  for (j in 7:ncol(group2)) {
    re<-cor.test(group2$age,as.numeric(group2[,j]), method =  "spearman",exact=FALSE)
    pvalue<-re[[3]]
    rvalue<-re[[4]]
    if(pvalue<0.05 & rvalue>0){
      pos_count<-pos_count+1
      
    }else if(pvalue<0.05 & rvalue<0){
      neg_count<-neg_count+1
    }
    
    dataouti<-cbind(colnames(group2)[j],pvalue,rvalue,strsplit(b[i],"_")[[1]][1])
    dataout<-rbind(dataout,dataouti)
    print(paste(i,j,ncol(group2),sep = "____"))
    
  }
  out1<-cbind(strsplit(b[i],"_")[[1]][1], pos_count,neg_count)
  
  re_out_p<-rbind(re_out_p,out1) 

  
}

re_out_p<-as.data.frame(re_out_p)
dataout<-as.data.frame(dataout)
write.table(re_out_p,"D:\\aging\\data\\geo\\region_DEG\\Spearman5_re_out_p.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

write.table(dataout,"D:\\aging\\data\\geo\\region_DEG\\Spearman5.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


dataout<-read.table("D:\\aging\\data\\geo\\region_DEG\\Spearman5.txt",header=T,sep = "\t", quote = "",fill = T)
dataout$regulation<-"Positive"
dataout$regulation[dataout$rvalue<=0]<-"Negative"

dataout<-dataout[dataout$pvalue<0.05,]


#################### BARPLOT #####################
re_out_p<-read.table("D:\\aging\\data\\geo\\region_DEG\\Spearman5_re_out_p.txt",header=T,sep = "\t", quote = "",fill = T)

re_out_p1<-as.data.frame(rbind(cbind(re_out_p$V1,re_out_p$pos_count,"pos"),
                               cbind(re_out_p$V1,re_out_p$neg_count,"neg"))) 

library(ggplot2)
library(ggalluvial)
library(tidyverse)
library(cols4all)


ggplot(re_out_p1, aes(x = reorder(V1,-as.numeric(V2)), y=as.numeric(V2), fill = V3,
                     stratum = V3, alluvium = V3)) +
  scale_fill_manual(values = c("#154182", "#99131D")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  labs(x="",y="Age-correlated genes")+theme(axis.text.x = element_text(angle = 45))+#字体大小
  geom_col(width = 0.6,
           color = NA, size = 0.5) + #去掉柱子的描边
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0.35,
            color = 'white', size = 0.5) +
  geom_alluvium(width = 0.6, alpha = 1, knot.pos = 0.35,
                fill = NA, color = 'white', size = 0.5) #再叠加一层白色描边加强一下效果
#################################Fig.2-G  #########################
rm(list=ls())
library(WGCNA)


## 转换为样品在行，基因在列的矩阵
expmean_gene <- read.table("D:\\aging\\data\\geo\\expmean_gene.txt",header=T,sep = "\t", quote = "")
expmean_gene1<-as.data.frame(t(expmean_gene))
expmean_gene1<-cbind(rownames(expmean_gene1),expmean_gene1)
colnames(expmean_gene1)<-expmean_gene1[1,]
expmean_gene2<-expmean_gene1[-1,]
setwd("D:\\aging\\data\\geo\\region_DEG\\group\\5")
b<-list.files()


poweri<-c()
for (i in 15:length(b)) {
  group <-read.table(b[i],header=T,sep = "\t", quote = "",fill = T)
  group2 <- merge(group,expmean_gene2,by.x = "childhood",by.y = "ID_REF")
  rownames(group2)<-group2[,1]
  LUAD_Expr<-group2[,-c(1,2)]
  datExpr<-matrix(as.numeric(unlist(LUAD_Expr) ),ncol = ncol(LUAD_Expr))
  colnames(datExpr)<-colnames(LUAD_Expr)
  rownames(datExpr)<-rownames(LUAD_Expr)
  traitData  <-read.table("D:\\aging\\data\\geo\\sample_group.txt",header=T,sep = "\t", quote = "")
  dim(traitData)
  # 形成一个类似于表达数据的数据框架，以保存临床特征
  # 提取行名
  femaleSamples = rownames(datExpr)
  # 数据匹配 返回匹配行
  traitRows = match(femaleSamples, traitData$Sample_geo_accession);
  # 提取指定要求行
  datTraits =as.data.frame( traitData[traitRows, 2]) 
  colnames(datTraits)<-c("age")
  # 提取行名
  rownames(datTraits) = traitData[traitRows, 5];
  # 垃圾回收
  collectGarbage();
  # Re-cluster samples
  # 画聚类图
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  # 画表型的热图
  # 将特征转换为颜色表示：白色表示低，红色表示高，灰色表示缺少条目
  # 如果signed为true 以绿色开头代表最大负值，以白色开头代表零附近的值，然后变为红色代表正值
  traitColors = numbers2colors(datTraits, signed =FALSE);
  pdf(file = paste("D:\\aging\\RESULT\\Spatiotemporal quantification\\wgcna\\",strsplit(b[i],"_")[[1]][1],"//",strsplit(b[i],"_")[[1]][1],"-sampleClustering.pdf", sep= "") , width = 12, height = 9);
  # Plot the sample dendrogram and the colors underneath.
  # 绘制出树状图和下面的颜色 
  plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits),main = "Sample dendrogram and trait heatmap")
  dev.off()
  
  # Choose a set of soft-thresholding powers
  # 给出候选的β值，c(1:10)表示1到10；seq(from = 12, to=20, by=2)表示从12开始间隔两个数到20
  powers = c(c(1:30))
  powers
  # Call the network topology analysis function 调用网络拓扑分析函数
  # verbose表示输出结果详细程度
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 0);
  if(is.na(sft$powerEstimate)==1 ){
    sft$powerEstimate=9
  }
  poweri_1<-cbind(strsplit(b[i],"_")[[1]][1],sft$powerEstimate) 
  poweri<-rbind(poweri,poweri_1)
  
  
  
  # datExpr表达数据，TOMType拓扑重叠矩阵计算方式，minModuleSize用于模块检测的最小模块尺寸,
  # reassignThreshold 是否在模块之间重新分配基因的p值比率阈值，mergeCutHeight 树状图切割高度
  # numericLabels 返回的模块应该用颜色（FALSE）还是数字（TRUE）标记,pamRespectsDendro树状图相关参数
  # saveTOMs 字符串的向量，saveTOMFileBase 包含包含共识拓扑重叠文件的文件名库的字符串
  net = blockwiseModules(datExpr, power = sft$powerEstimate,TOMType = "unsigned", minModuleSize = 30,reassignThreshold = 0, 
                         mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,
                         saveTOMFileBase = "TOM",verbose = 3)
  table(net$colors)
  # open a graphics window
  # sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  # 将标签转化为绘图颜色
  mergedColors = labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath
  # 绘制树状图和下面的模块颜色
  # dendroLabels树状图标签。设置为FALSE完全禁用树状图标签；设置为NULL使用的行标签datExpr
  # addGuide是否应在树状图中添加垂直的“指导线”？线条使识别单个样本的颜色代码更加容易。
  pdf(file = paste("D:\\aging\\RESULT\\Spatiotemporal quantification\\wgcna\\",strsplit(b[i],"_")[[1]][1],"//",strsplit(b[i],"_")[[1]][1],"-ClusterDendrogram.pdf", sep= "") , width = 12, height = 9);
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                      dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
  dev.off()
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  
  
  # Define numbers of genes and samples
  # 获得基因数和样本数
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  
  # Recalculate MEs with color labels
  # 用彩色标签重新计算MEs
  # 在给定的单个数据集中计算模块的模块本征基因
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  # 对给定的（特征）向量进行重新排序，以使相似的向量（通过相关性度量）彼此相邻
  MEs = orderMEs(MEs0)
  
  # 计算module的ME值与表型的相关系数
  moduleTraitCor = cor(MEs, datTraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  # sizeGrWindow(10,6)
  # 显示相关性及其p值
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  pdf(file = paste("D:\\aging\\RESULT\\Spatiotemporal quantification\\wgcna\\",strsplit(b[i],"_")[[1]][1],"//",strsplit(b[i],"_")[[1]][1],"-ageRalationship.pdf", sep= "") , width = 12, height = 9);
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = colorRampPalette(c("#879f68", "#FFFFFF", "#e4c24b"))(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  # Define variable weight containing the weight column of datTrait
  # 定义包含数据特征权重列的变量权重
  age = as.data.frame(datTraits$age);
  names(age) = "age"
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  # 模块的名称(颜色) substring提取文本从第3个字母开始
  modNames = substring(names(MEs), 3)
  # 基因和模块的相关系数
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  #gene和性状的关系
  geneTraitSignificance = as.data.frame(cor(datExpr, age, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(age), sep="");
  names(GSPvalue) = paste("p.GS.", names(age), sep="");
  
  
  
  geneInfo0 = data.frame(moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)
  # Order modules by their significance for weight
  modOrder = order(-abs(cor(MEs, age, use = "p")));
  # Add module membership information in the chosen order
  for (mod in 1:ncol(geneModuleMembership))
  {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
  # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.age));
  geneInfo = geneInfo0[geneOrder, ]
  write.table(geneInfo, file = paste("D:\\aging\\RESULT\\Spatiotemporal quantification\\wgcna\\",strsplit(b[i],"_")[[1]][1],"//",strsplit(b[i],"_")[[1]][1],"geneInfo.txt", sep= ""),col.names = T, row.names = T,sep = "\t" ,append = FALSE, quote = F)
  
  
  # Calculate topological overlap anew: this could be done more efficiently by saving the TOM
  # calculated during module detection, but let us do it again here.
  # 重新计算拓扑重叠：通过保存TOM可以更有效地完成此操作
  # 是在模块检测期间计算的，但让我们在这里再次进行。
  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate);
  # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
  # 变换dissTOM
  plotTOM = dissTOM^7;
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA;
  # Call the plot function
  # sizeGrWindow(9,9)
  # 基因的聚类树聚类时的距离为1-TOM值结合基因间的距离，即1-TOM值，用热图展示
  # TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  
  nSelect = 400
  # For reproducibility, we set the random seed
  set.seed(10);
  select = sample(nGenes, size = nSelect);
  selectTOM = dissTOM[select, select];
  # There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
  # 重新画聚类图
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select];
  # Open a graphical window
  # sizeGrWindow(9,9)
  # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
  # the color palette; setting the diagonal to NA also improves the clarity of the plot
  plotDiss = selectTOM^7;
  diag(plotDiss) = NA;
  pdf(file = paste("D:\\aging\\RESULT\\Spatiotemporal quantification\\wgcna\\",strsplit(b[i],"_")[[1]][1],"//",strsplit(b[i],"_")[[1]][1],"-Networkheatmap.pdf", sep= "") , width = 12, height = 9);
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
  dev.off()
  
  # Recalculate module eigengenes
  # 重新计算基因特征值
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  # Isolate weight from the clinical traits
  age = as.data.frame(datTraits$age);
  names(age) = "age"
  # Add the weight to existing module eigengenes
  MET = orderMEs(cbind(MEs, age))
  # Plot the relationships among the eigengenes and the trait
  #sizeGrWindow(5,7.5);
  par(cex = 0.9)
  # 画树形图
  # marDendro给出树状图的边距设置，marHeatmap热图边距设置
  pdf(file = paste("D:\\aging\\RESULT\\Spatiotemporal quantification\\wgcna\\",strsplit(b[i],"_")[[1]][1],"//",strsplit(b[i],"_")[[1]][1],"-marHeatmap.pdf", sep= "") , width = 12, height = 9);
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
  dev.off()
  
  moduleColorsnumber<-unique(moduleColors)
  for (k in 1:length(moduleColorsnumber)) {
    module = moduleColorsnumber[k]
    # 匹配列
    column = match(module, modNames);
    moduleGenes = moduleColors==module;
    #Recalculate topological overlap if needed
    TOM =1-dissTOM
    #################
    dimnames(TOM) <- list(colnames(datExpr), colnames(datExpr))
    modTOM <- TOM[moduleGenes, moduleGenes]
    
    cyt <- exportNetworkToCytoscape(
      modTOM,
      edgeFile = paste("D:\\aging\\RESULT\\Spatiotemporal quantification\\wgcna\\",strsplit(b[i],"_")[[1]][1],"//",strsplit(b[i],"_")[[1]][1],"_",module,"_","edges.txt", sep= ""),
      nodeFile = paste("D:\\aging\\RESULT\\Spatiotemporal quantification\\wgcna\\",strsplit(b[i],"_")[[1]][1],"//",strsplit(b[i],"_")[[1]][1],"_",module,"_","nodes.txt", sep= ""),
      weighted = TRUE,
      threshold = 0
    )
    print(paste("region",i,"module",length(moduleColorsnumber),k,sep = "-----"))
  }
  
  
  
  
}

write.table(poweri,"D:\\aging\\RESULT\\Spatiotemporal quantification\\wgcna\\power.txt",col.names = F, row.names = F,sep = "\t" ,append = FALSE, quote = F)


################ volcano  ###############
library(ggplot2)
library(tidyverse)
library(ggrepel)


setwd("D:\\aging\\RESULT\\1-Spatiotemporal quantification\\wgcna\\region\\5")
c<-list.files()
df<-c()
for (i in 1:length(c)) {
  geneInfo <- read.table(paste(c[i],"\\",c[i],"geneInfo.txt",sep=""),header=T,sep = "\t", quote = "",row.names=NULL )[,1:4]
  geneInfo1<-cbind(geneInfo,c[i])
  df<-rbind(df,geneInfo1)
}
colnames(df)<-c("gene","module","r","p","cluster")
df <- as.data.frame(df)

#添加显著性标签：
df$label <- ifelse(df$r<0,"negative","positive")
dt<-df[df$p<0.01,]
#write.table(dt,"D:\\aging\\RESULT\\R9online tool\\1-1\\WGCNA_point.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


#根据图p中log2FC区间确定背景柱长度：
dfbar<-read.table("D:\\aging\\RESULT\\1-Spatiotemporal quantification\\wgcna\\module-age-re5.txt",header=T,sep = "\t", quote = "")
dfbarsig<-dfbar[dfbar$P<0.05,]

dfbarsigup<-as.data.frame(table(dfbarsig[dfbarsig$R>=0,]$region))
dfbarsigup$Freq<- dfbarsigup$Freq/max(dfbarsigup$Freq)
dfbarsigdown<-as.data.frame(table(dfbarsig[dfbarsig$R<=0,]$region))
dfbarsigdown$Freq<- -dfbarsigdown$Freq/max(dfbarsigdown$Freq)


dfbarsigup1<-as.data.frame(table(dfbarsig[dfbarsig$R>=0,]$region))
dfbarsigdown1<-as.data.frame(table(dfbarsig[dfbarsig$R<=0,]$region))
dfbarsigdown1$Freq<- -dfbarsigdown1$Freq

dfbarsigup$total<-dfbarsigup$Freq-dfbarsigdown$Freq

#绘制背景柱：
p1 <- ggplot()+
  geom_col(data = dfbarsigup1,mapping = aes(x = Var1,y = Freq),fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbarsigdown1, mapping = aes(x = Var1,y = Freq), fill = "#dcdcdc",alpha = 0.6)+
  theme_minimal()+
  #guides(fill=FALSE)+
  theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
        axis.line.y = element_line(color = "black",size = 1.2),
        #axis.line.x = element_blank(),
        #axis.text.x = element_blank(),
        panel.grid = element_blank() )+
  scale_x_discrete(limits=c("IPC","V1C","ITC","S1C",
                            "DFC","VFC","MD","M1C",
                            "AMY","STR","STC","CBC",
                            "A1C","OFC","HIP","MFC"))
p1


#把散点火山图叠加到背景柱上：
p2 <- ggplot()+
  geom_col(data = dfbarsigup, mapping = aes(x = reorder(Var1,-total) ,y = Freq), fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbarsigdown,mapping = aes(x = Var1,y = Freq),fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,aes(x = cluster, y = r, color = label), size = 0.85,width =0.4)
#scale_x_discrete(limits=c("IPC","V1C","ITC","S1C",
#                         "DFC","VFC","MD","M1C",
#                        "AMY","STR","STC","CBC",
#                       "A1C","OFC","HIP","MFC"))
p2
#添加X轴的cluster色块标签：
dfcol<-data.frame(x=c(1:16),
                  y=0,
                  label=c("A1C","STR","CBC","HIP","M1C",
                          "V1C","VFC","DFC","S1C","ITC",
                          "AMY","STC","IPC","OFC","MD","MFC"))



mycol <- c("#12803b","#136cb6","#f18e25","#211d1e",
           "#735d96","#c5b5d1","#5A7EB3","#f5c51e",
           "#f4da9a","#ac536a","#ec95b3","#f8bf89",
           "#57a8d7","#cf1223","#cacbd0","#f1a7a4")


p3 <- p2 + geom_tile(data = dfcol,aes(x=x,y=y),height=0.3, color = "black", fill = mycol,alpha = 0.6,show.legend = F)


p5 <- p3 +
  scale_color_manual(name=NULL,values = c("#55739e","#af545b"))

#修改X/Y轴标题和添加cluster数字：
p6 <- p5+
  labs(x="Regoins",y="Correlation of Module Genes with Age")+
  geom_text(data=dfcol,aes(x=x,y=y,label=label),size =4, color ="white")



#自定义主题美化：
p7 <- p6+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,color = "black", face = "bold"),
    axis.line.y = element_line(color = "black", size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position="none" )
p7
################# Fig2.HI #####################
rm(list = ls())
library(pheatmap)
sample_group <-read.table("D:\\aging\\data\\geo\\sample_group.txt",header=T,sep = "\t", quote = "",fill = T)
expmean_gene<-read.table("D:\\aging\\data\\geo\\expmean_gene.txt",header=T,sep = "\t", quote = "")
sample_group <-cbind(rownames(sample_group),sample_group)
group_list<-c("childhood","Adolescence","Young","Middle","Late")
for (k in 1:5) {
  group<-sample_group[sample_group$group==group_list[k],]
  region1<-c("CBC","MD","STR","NCX","HIP","AMY")
  #region1<-c("V1C","ITC","IPC","A1C","STC","MFC",
  #           "OFC","DFC","VFC","M1C","S1C")
  group_mat<-expmean_gene[,1]
  data<-matrix(0,ncol = length(region1),nrow = length(region1))
  for (i in 1:length(region1)) {
    region_gene<-expmean_gene[,sapply(group$samples[group$region1==region1[i]],function(x){which(colnames(expmean_gene)==x)})]
    #region_gene<-expmean_gene[,sapply(group$samples[group$region==region1[i]],function(x){which(colnames(expmean_gene)==x)})]
    
    region_mean<-as.data.frame(apply(region_gene,1,mean)) 
    group_mat<-cbind(group_mat,region_mean)

    
    for (j in 1:i) {
      re<-cor.test(group_mat[,i+1],group_mat[,j+1],method =  "spearman")
      data[i,j]=re$estimate
      data[j,i]=re$estimate
      

    }
  }
  colnames(data)<-region1
  rownames(data)<-region1
  myColor <- colorRampPalette(c("#154182", "#FFFFFF", "#99131D"))(length(seq(0.98,1,by=0.001)))
  myBreaks <- seq(0.98,1,by=0.001)
  #write.table(data,paste( "D:\\aging\\RESULT\\9online tool\\1-1\\",group_list[k],"_ Heatmap matrix_6.txt",sep="" ),col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
  
  pheatmap(data, 
           #annotation_col = annotation_col, 
           #annotation_row = annotation_row,
          
          cluster_rows = F,
          cluster_cols = F,
           #cellwidth = 20, cellheight = 15,
          color = myColor,
          breaks=myBreaks, 
          annotation_legend = FALSE,
          filename=paste("D:\\aging\\RESULT\\Spatiotemporal quantification\\",group_list[k],"_ Heatmap matrix.pdf"),
          main = group_list[k])
  

}







