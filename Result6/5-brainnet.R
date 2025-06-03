###############HOA point ###############

rm(list=ls())
library(tidyverse)       
library(ggsci)           
library(ggExtra)         
library(ggpmisc)        
library(palmerpenguins) 
library(ggpubr)
library(reshape2)
Node<-read.table("D:\\aging\\sample_stats\\DLBS\\MSN\\HOAatlas\\0028640.nii.txt",header=F,sep = "\t", quote = "")[1,]
Node<-as.data.frame(t(Node))
DLBS_HOAatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\DLBS\\MSN\\HOAatlas\\DLBS_HOAatlas_Period_13_cor.txt",header=F,sep = "\t", quote = "")
DLBS_HOAatlas_Period_13_cor1<-cbind(Node,DLBS_HOAatlas_Period_13_cor)

colnames(DLBS_HOAatlas_Period_13_cor1)<-c("regions",Node)
DLBS_HOAatlas_Period_13_cor2 <- melt(DLBS_HOAatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_HOAatlas_Period_13_cor2$value<-(DLBS_HOAatlas_Period_13_cor2$value+1)/2

DLBS_HOAatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\DLBS\\MIND\\HOAatlas\\DLBS_HOAatlas_Period_13_cor.txt",header=F,sep = "\t", quote = "")
DLBS_HOAatlas_Period_13_cor1_mind<-cbind(Node,DLBS_HOAatlas_Period_13_cor_mind)
colnames(DLBS_HOAatlas_Period_13_cor1_mind)<-c("regions",Node)
DLBS_HOAatlas_Period_13_cor2_mind <- melt(DLBS_HOAatlas_Period_13_cor1_mind,id.vars = c("regions"))

data<-cbind(DLBS_HOAatlas_Period_13_cor2,DLBS_HOAatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"

data1<-cbind(data$r_MSN,data$r_MIND,data$group) 
colnames(data1)<-c("r_MSN","r_MIND","group")
write.table(data1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF6_HOA\\B_Y.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


ggplot(data, aes( r_MSN, r_MIND,color=group)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm",color = "#A0C9D1") + 
  stat_cor( geom = "text",label.x = 0)+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(family = "sans", face = "bold"))

ggplot(data, aes(r_MSN, r_MIND)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm") + 
  stat_cor( geom = "text",label.x = 0)+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))


DLBS_HOAatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\DLBS\\MSN\\HOAatlas\\DLBS_HOAatlas_Period_14_cor.txt",header=F,sep = "\t", quote = "")
DLBS_HOAatlas_Period_13_cor1<-cbind(Node,DLBS_HOAatlas_Period_13_cor)
colnames(DLBS_HOAatlas_Period_13_cor1)<-c("regions",Node)
DLBS_HOAatlas_Period_13_cor2 <- melt(DLBS_HOAatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_HOAatlas_Period_13_cor2$value<-(DLBS_HOAatlas_Period_13_cor2$value+1)/2

DLBS_HOAatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\DLBS\\MIND\\HOAatlas\\DLBS_HOAatlas_Period_14_cor.txt",header=F,sep = "\t", quote = "")
DLBS_HOAatlas_Period_13_cor1_mind<-cbind(Node,DLBS_HOAatlas_Period_13_cor_mind)
colnames(DLBS_HOAatlas_Period_13_cor1_mind)<-c("regions",Node)
DLBS_HOAatlas_Period_13_cor2_mind <- melt(DLBS_HOAatlas_Period_13_cor1_mind,id.vars = c("regions"))

data<-cbind(DLBS_HOAatlas_Period_13_cor2,DLBS_HOAatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"
data1<-cbind(data$r_MSN,data$r_MIND,data$group) 
colnames(data1)<-c("r_MSN","r_MIND","group")
write.table(data1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF6_HOA\\B_M.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


ggplot(data, aes(r_MSN, r_MIND,color=group)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm",color = "#436493") + 
  stat_cor( geom = "text",label.x = 0)+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))

ggplot(data, aes(r_MSN, r_MIND)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm") + 
  stat_cor( geom = "text",label.x = 0)+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))

DLBS_HOAatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\DLBS\\MSN\\HOAatlas\\DLBS_HOAatlas_Period_15_cor.txt",header=F,sep = "\t", quote = "")
DLBS_HOAatlas_Period_13_cor1<-cbind(Node,DLBS_HOAatlas_Period_13_cor)
colnames(DLBS_HOAatlas_Period_13_cor1)<-c("regions",Node)
DLBS_HOAatlas_Period_13_cor2 <- melt(DLBS_HOAatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_HOAatlas_Period_13_cor2$value<-(DLBS_HOAatlas_Period_13_cor2$value+1)/2

DLBS_HOAatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\DLBS\\MIND\\HOAatlas\\DLBS_HOAatlas_Period_15_cor.txt",header=F,sep = "\t", quote = "")
DLBS_HOAatlas_Period_13_cor1_mind<-cbind(Node,DLBS_HOAatlas_Period_13_cor_mind)
colnames(DLBS_HOAatlas_Period_13_cor1_mind)<-c("regions",Node)
DLBS_HOAatlas_Period_13_cor2_mind <- melt(DLBS_HOAatlas_Period_13_cor1_mind,id.vars = c("regions"))

data<-cbind(DLBS_HOAatlas_Period_13_cor2,DLBS_HOAatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"
data1<-cbind(data$r_MSN,data$r_MIND,data$group) 
colnames(data1)<-c("r_MSN","r_MIND","group")
write.table(data1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF6_HOA\\B_L.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



ggplot(data, aes(r_MSN, r_MIND,color=group)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm",color = "#1C2C63") + 
  stat_cor( geom = "text",label.x = 0)+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))
ggplot(data, aes(r_MSN, r_MIND)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm") + 
  stat_cor( geom = "text",label.x = 0)+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))


################ gbm #################
DLBS_HOAatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\GBM\\MSN\\HOAatlas\\GBM_HOAatlas_Period_13_cor.txt",header=F,sep = "\t", quote = "")
DLBS_HOAatlas_Period_13_cor1<-cbind(Node,DLBS_HOAatlas_Period_13_cor)
colnames(DLBS_HOAatlas_Period_13_cor1)<-c("regions",Node)
DLBS_HOAatlas_Period_13_cor2 <- melt(DLBS_HOAatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_HOAatlas_Period_13_cor2$value<-(DLBS_HOAatlas_Period_13_cor2$value+1)/2


DLBS_HOAatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\GBM\\MIND\\HOAatlas\\GBM_HOAatlas_Period_13_cor.txt",header=F,sep = "\t", quote = "")
DLBS_HOAatlas_Period_13_cor1_mind<-cbind(Node,DLBS_HOAatlas_Period_13_cor_mind)
colnames(DLBS_HOAatlas_Period_13_cor1_mind)<-c("regions",Node)
DLBS_HOAatlas_Period_13_cor2_mind <- melt(DLBS_HOAatlas_Period_13_cor1_mind,id.vars = c("regions"))




data<-cbind(DLBS_HOAatlas_Period_13_cor2,DLBS_HOAatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"
data1<-cbind(data$r_MSN,data$r_MIND,data$group) 
colnames(data1)<-c("r_MSN","r_MIND","group")
write.table(data1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF6_HOA\\D_Y.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


ggplot(data, aes(r_MSN, r_MIND,color=group)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm",color = "#A0C9D1") + 
  stat_cor( geom = "text",label.x = 0)+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))
ggplot(data, aes(r_MSN, r_MIND)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm") + 
  stat_cor( geom = "text",label.x = 0)+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))

DLBS_HOAatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\GBM\\MSN\\HOAatlas\\GBM_HOAatlas_Period_14_cor.txt",header=F,sep = "\t", quote = "")
DLBS_HOAatlas_Period_13_cor1<-cbind(Node,DLBS_HOAatlas_Period_13_cor)
colnames(DLBS_HOAatlas_Period_13_cor1)<-c("regions",Node)
DLBS_HOAatlas_Period_13_cor2 <- melt(DLBS_HOAatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_HOAatlas_Period_13_cor2$value<-(DLBS_HOAatlas_Period_13_cor2$value+1)/2


DLBS_HOAatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\GBM\\MIND\\HOAatlas\\GBM_HOAatlas_Period_14_cor.txt",header=F,sep = "\t", quote = "")
DLBS_HOAatlas_Period_13_cor1_mind<-cbind(Node,DLBS_HOAatlas_Period_13_cor_mind)
colnames(DLBS_HOAatlas_Period_13_cor1_mind)<-c("regions",Node)
DLBS_HOAatlas_Period_13_cor2_mind <- melt(DLBS_HOAatlas_Period_13_cor1_mind,id.vars = c("regions"))

data<-cbind(DLBS_HOAatlas_Period_13_cor2,DLBS_HOAatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"
data1<-cbind(data$r_MSN,data$r_MIND,data$group) 
colnames(data1)<-c("r_MSN","r_MIND","group")
write.table(data1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF6_HOA\\D_M.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


ggplot(data, aes(r_MSN, r_MIND,color=group)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm",color = "#436493") + 
  stat_cor( geom = "text",label.x = 0)+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))
ggplot(data, aes(r_MSN, r_MIND)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm") + 
  stat_cor( geom = "text",label.x = 0)+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))

DLBS_HOAatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\GBM\\MSN\\HOAatlas\\GBM_HOAatlas_Period_15_cor.txt",header=F,sep = "\t", quote = "")
DLBS_HOAatlas_Period_13_cor1<-cbind(Node,DLBS_HOAatlas_Period_13_cor)
colnames(DLBS_HOAatlas_Period_13_cor1)<-c("regions",Node)
DLBS_HOAatlas_Period_13_cor2 <- melt(DLBS_HOAatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_HOAatlas_Period_13_cor2$value<-(DLBS_HOAatlas_Period_13_cor2$value+1)/2


DLBS_HOAatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\GBM\\MIND\\HOAatlas\\GBM_HOAatlas_Period_15_cor.txt",header=F,sep = "\t", quote = "")
DLBS_HOAatlas_Period_13_cor1_mind<-cbind(Node,DLBS_HOAatlas_Period_13_cor_mind)
colnames(DLBS_HOAatlas_Period_13_cor1_mind)<-c("regions",Node)
DLBS_HOAatlas_Period_13_cor2_mind <- melt(DLBS_HOAatlas_Period_13_cor1_mind,id.vars = c("regions"))
data<-cbind(DLBS_HOAatlas_Period_13_cor2,DLBS_HOAatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"


data1<-cbind(data$r_MSN,data$r_MIND,data$group) 
colnames(data1)<-c("r_MSN","r_MIND","group")
write.table(data1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF6_HOA\\D_L.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

ggplot(data, aes(r_MSN, r_MIND,color=group)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm",color = "#1C2C63") + 
  stat_cor( geom = "text",label.x = 0)+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(family = "sans", face = "bold"))
ggplot(data, aes(r_MSN, r_MIND)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm") + 
  stat_cor( geom = "text",label.x = 0)+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))

###############DK point ###############


rm(list=ls())
library(tidyverse)       
library(ggsci)           
library(ggExtra)         
library(ggpmisc)        
library(palmerpenguins) 
library(ggpubr)
library(reshape2)
Node<-read.table("D:\\aging\\sample_stats\\DLBS\\MSN\\DKTatlas\\0028640.nii.txt",header=F,sep = "\t", quote = "")[1,]
Node<-as.data.frame(t(Node))
DLBS_DKTatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\DLBS\\MSN\\DKTatlas\\DLBS_DKTatlas_Period_13_cor.txt",header=F,sep = "\t", quote = "")
DLBS_DKTatlas_Period_13_cor1<-cbind(Node,DLBS_DKTatlas_Period_13_cor)
colnames(DLBS_DKTatlas_Period_13_cor1)<-c("regions",Node)
DLBS_DKTatlas_Period_13_cor2 <- melt(DLBS_DKTatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_DKTatlas_Period_13_cor2$value<-(DLBS_DKTatlas_Period_13_cor2$value+1)/2
DLBS_DKTatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\DLBS\\MIND\\DKTatlas\\DLBS_DKTatlas_Period_13_cor.txt",header=F,sep = "\t", quote = "")
DLBS_DKTatlas_Period_13_cor1_mind<-cbind(Node,DLBS_DKTatlas_Period_13_cor_mind)
colnames(DLBS_DKTatlas_Period_13_cor1_mind)<-c("regions",Node)
DLBS_DKTatlas_Period_13_cor2_mind <- melt(DLBS_DKTatlas_Period_13_cor1_mind,id.vars = c("regions"))
data<-cbind(DLBS_DKTatlas_Period_13_cor2,DLBS_DKTatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"
data1<-cbind(data$r_MSN,data$r_MIND,data$group) 
colnames(data1)<-c("r_MSN","r_MIND","group")
write.table(data1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF7_DK\\B_Y.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

ggplot(data, aes(r_MSN, r_MIND,color=group)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm",color = "#A0C9D1") + 
  stat_cor( geom = "text",label.x = 0)+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(family = "sans", face = "bold"))
ggplot(data, aes(r_MSN, r_MIND)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm") + 
  stat_cor( geom = "text",label.x = 0)+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(family = "sans", face = "bold"))

DLBS_DKTatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\DLBS\\MSN\\DKTatlas\\DLBS_DKTatlas_Period_14_cor.txt",header=F,sep = "\t", quote = "")
DLBS_DKTatlas_Period_13_cor1<-cbind(Node,DLBS_DKTatlas_Period_13_cor)
colnames(DLBS_DKTatlas_Period_13_cor1)<-c("regions",Node)
DLBS_DKTatlas_Period_13_cor2 <- melt(DLBS_DKTatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_DKTatlas_Period_13_cor2$value<-(DLBS_DKTatlas_Period_13_cor2$value+1)/2
DLBS_DKTatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\DLBS\\MIND\\DKTatlas\\DLBS_DKTatlas_Period_14_cor.txt",header=F,sep = "\t", quote = "")
DLBS_DKTatlas_Period_13_cor1_mind<-cbind(Node,DLBS_DKTatlas_Period_13_cor_mind)
colnames(DLBS_DKTatlas_Period_13_cor1_mind)<-c("regions",Node)
DLBS_DKTatlas_Period_13_cor2_mind <- melt(DLBS_DKTatlas_Period_13_cor1_mind,id.vars = c("regions"))
data<-cbind(DLBS_DKTatlas_Period_13_cor2,DLBS_DKTatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"
data1<-cbind(data$r_MSN,data$r_MIND,data$group) 
colnames(data1)<-c("r_MSN","r_MIND","group")
write.table(data1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF7_DK\\B_M.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



ggplot(data, aes(r_MSN, r_MIND,color=group)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm",color = "#436493") + 
  stat_cor( geom = "text",label.x = 0)+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(family = "sans", face = "bold"))
ggplot(data, aes(r_MSN, r_MIND)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm") + 
  stat_cor( geom = "text",label.x = 0)+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(family = "sans", face = "bold"))


DLBS_DKTatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\DLBS\\MSN\\DKTatlas\\DLBS_DKTatlas_Period_15_cor.txt",header=F,sep = "\t", quote = "")
DLBS_DKTatlas_Period_13_cor1<-cbind(Node,DLBS_DKTatlas_Period_13_cor)
colnames(DLBS_DKTatlas_Period_13_cor1)<-c("regions",Node)
DLBS_DKTatlas_Period_13_cor2 <- melt(DLBS_DKTatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_DKTatlas_Period_13_cor2$value<-(DLBS_DKTatlas_Period_13_cor2$value+1)/2
DLBS_DKTatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\DLBS\\MIND\\DKTatlas\\DLBS_DKTatlas_Period_15_cor.txt",header=F,sep = "\t", quote = "")
DLBS_DKTatlas_Period_13_cor1_mind<-cbind(Node,DLBS_DKTatlas_Period_13_cor_mind)
colnames(DLBS_DKTatlas_Period_13_cor1_mind)<-c("regions",Node)
DLBS_DKTatlas_Period_13_cor2_mind <- melt(DLBS_DKTatlas_Period_13_cor1_mind,id.vars = c("regions"))
data<-cbind(DLBS_DKTatlas_Period_13_cor2,DLBS_DKTatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"
data1<-cbind(data$r_MSN,data$r_MIND,data$group) 
colnames(data1)<-c("r_MSN","r_MIND","group")
write.table(data1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF7_DK\\B_L.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

ggplot(data, aes(r_MSN, r_MIND,color=group)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm",color = "#1C2C63") + 
  stat_cor( geom = "text",label.x = 0)+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(family = "sans", face = "bold"))
ggplot(data, aes(r_MSN, r_MIND)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm") + 
  stat_cor( geom = "text",label.x = 0)+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(family = "sans", face = "bold"))


################ gbm #################
DLBS_DKTatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\GBM\\MSN\\DKTatlas\\GBM_DKTatlas_Period_13_cor.txt",header=F,sep = "\t", quote = "")
DLBS_DKTatlas_Period_13_cor1<-cbind(Node,DLBS_DKTatlas_Period_13_cor)
colnames(DLBS_DKTatlas_Period_13_cor1)<-c("regions",Node)
DLBS_DKTatlas_Period_13_cor2 <- melt(DLBS_DKTatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_DKTatlas_Period_13_cor2$value<-(DLBS_DKTatlas_Period_13_cor2$value+1)/2
DLBS_DKTatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\GBM\\MIND\\DKTatlas\\GBM_DKTatlas_Period_13_cor.txt",header=F,sep = "\t", quote = "")
DLBS_DKTatlas_Period_13_cor1_mind<-cbind(Node,DLBS_DKTatlas_Period_13_cor_mind)
colnames(DLBS_DKTatlas_Period_13_cor1_mind)<-c("regions",Node)
DLBS_DKTatlas_Period_13_cor2_mind <- melt(DLBS_DKTatlas_Period_13_cor1_mind,id.vars = c("regions"))
data<-cbind(DLBS_DKTatlas_Period_13_cor2,DLBS_DKTatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"
data1<-cbind(data$r_MSN,data$r_MIND,data$group) 
colnames(data1)<-c("r_MSN","r_MIND","group")
write.table(data1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF7_DK\\D_Y.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
ggplot(data, aes(r_MSN, r_MIND,color=group)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm",color = "#A0C9D1") + 
  stat_cor( geom = "text",label.x = 0)+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))
ggplot(data, aes(r_MSN, r_MIND)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm") + 
  stat_cor( geom = "text",label.x = 0)+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))

DLBS_DKTatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\GBM\\MSN\\DKTatlas\\GBM_DKTatlas_Period_14_cor.txt",header=F,sep = "\t", quote = "")
DLBS_DKTatlas_Period_13_cor1<-cbind(Node,DLBS_DKTatlas_Period_13_cor)
colnames(DLBS_DKTatlas_Period_13_cor1)<-c("regions",Node)
DLBS_DKTatlas_Period_13_cor2 <- melt(DLBS_DKTatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_DKTatlas_Period_13_cor2$value<-(DLBS_DKTatlas_Period_13_cor2$value+1)/2


DLBS_DKTatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\GBM\\MIND\\DKTatlas\\GBM_DKTatlas_Period_14_cor.txt",header=F,sep = "\t", quote = "")
DLBS_DKTatlas_Period_13_cor1_mind<-cbind(Node,DLBS_DKTatlas_Period_13_cor_mind)
colnames(DLBS_DKTatlas_Period_13_cor1_mind)<-c("regions",Node)
DLBS_DKTatlas_Period_13_cor2_mind <- melt(DLBS_DKTatlas_Period_13_cor1_mind,id.vars = c("regions"))

data<-cbind(DLBS_DKTatlas_Period_13_cor2,DLBS_DKTatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"
data1<-cbind(data$r_MSN,data$r_MIND,data$group) 
colnames(data1)<-c("r_MSN","r_MIND","group")
write.table(data1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF7_DK\\D_M.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



ggplot(data, aes(r_MSN, r_MIND,color=group)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm",color = "#436493") + 
  stat_cor( geom = "text",label.x = 0)+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))
ggplot(data, aes(r_MSN, r_MIND)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm") + 
  stat_cor( geom = "text",label.x = 0)+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))

DLBS_DKTatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\GBM\\MSN\\DKTatlas\\GBM_DKTatlas_Period_15_cor.txt",header=F,sep = "\t", quote = "")
DLBS_DKTatlas_Period_13_cor1<-cbind(Node,DLBS_DKTatlas_Period_13_cor)
colnames(DLBS_DKTatlas_Period_13_cor1)<-c("regions",Node)
DLBS_DKTatlas_Period_13_cor2 <- melt(DLBS_DKTatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_DKTatlas_Period_13_cor2$value<-(DLBS_DKTatlas_Period_13_cor2$value+1)/2


DLBS_DKTatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\GBM\\MIND\\DKTatlas\\GBM_DKTatlas_Period_15_cor.txt",header=F,sep = "\t", quote = "")
DLBS_DKTatlas_Period_13_cor1_mind<-cbind(Node,DLBS_DKTatlas_Period_13_cor_mind)
colnames(DLBS_DKTatlas_Period_13_cor1_mind)<-c("regions",Node)
DLBS_DKTatlas_Period_13_cor2_mind <- melt(DLBS_DKTatlas_Period_13_cor1_mind,id.vars = c("regions"))
data<-cbind(DLBS_DKTatlas_Period_13_cor2,DLBS_DKTatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"

data1<-cbind(data$r_MSN,data$r_MIND,data$group) 
colnames(data1)<-c("r_MSN","r_MIND","group")
write.table(data1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF7_DK\\D_L.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


ggplot(data, aes(r_MSN, r_MIND,color=group)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm",color = "#1C2C63") + 
  stat_cor( geom = "text",label.x = 0)+
  scale_color_manual(values = c("#B2522D",
                                "#436493"))+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(family = "sans", face = "bold"))
ggplot(data, aes(r_MSN, r_MIND)) +
  geom_point() +  # 添加散点图层，点的大小表示体重
  geom_smooth( method = "lm") + 
  stat_cor( geom = "text",label.x = 0)+
  labs(x = "MSN", y = "MIND") +  # 设置坐标轴标签
  theme_classic()  +
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme(legend.position = "none")+
  
  theme(axis.title = element_text(family = "sans", face = "bold"))

 ########## bar plot  ################
rm(list = ls())
setwd("D:\\aging\\RESULT\\5MS_Tumor\\")
Node_Brodmann64<-read.table("D:\\aging\\sample_stats\\Node_Brodmann64.txt",header=F,sep = "\t", quote = "")
list0<-c("DLBS","GBM")
for (j in 1:2) {
  list1<-c("13","14","15")
  DLBS_BAatlas_Period_cor4<-c()
  for (i in 1:3) {
    DLBS_BAatlas_Period_cor<-read.table(paste("D:\\aging\\sample_stats\\",list0[j],"\\MSN\\BAatlas\\",list0[j],"_BAatlas_Period_",list1[i],"_cor.txt",sep = ""),header=F,sep = "\t", quote = "")
    DLBS_BAatlas_Period_cor1<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_cor)
    colnames(DLBS_BAatlas_Period_cor1)<-c("regions",Node_Brodmann64$V6)
    DLBS_BAatlas_Period_cor2 <- melt(DLBS_BAatlas_Period_cor1,id.vars = c("regions"))
    DLBS_BAatlas_Period_cor2$value<-(DLBS_BAatlas_Period_cor2$value+1)/2
    DLBS_BAatlas_Period_cor2_MSN<-cbind(DLBS_BAatlas_Period_cor2,paste("Period_",list1[i],sep="") ,"MSN")
    colnames(DLBS_BAatlas_Period_cor2_MSN)<-c("region1","region2","R","Age","Net")
    
    DLBS_BAatlas_Period_cor<-read.table(paste("D:\\aging\\sample_stats\\",list0[j],"\\MIND\\BAatlas\\",list0[j],"_BAatlas_Period_",list1[i],"_cor.txt",sep = ""),header=F,sep = "\t", quote = "")
    DLBS_BAatlas_Period_cor1<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_cor)
    colnames(DLBS_BAatlas_Period_cor1)<-c("regions",Node_Brodmann64$V6)
    DLBS_BAatlas_Period_cor2 <- melt(DLBS_BAatlas_Period_cor1,id.vars = c("regions"))
    DLBS_BAatlas_Period_cor2$value<-(DLBS_BAatlas_Period_cor2$value+1)/2
    DLBS_BAatlas_Period_cor2_MIND<-cbind(DLBS_BAatlas_Period_cor2,paste("Period_",list1[i],sep=""),"MIND")
    colnames(DLBS_BAatlas_Period_cor2_MIND)<-c("region1","region2","R","Age","Net")
    
    DLBS_BAatlas_Period_cor3<-rbind(DLBS_BAatlas_Period_cor2_MSN,DLBS_BAatlas_Period_cor2_MIND)
    DLBS_BAatlas_Period_cor4<-rbind(DLBS_BAatlas_Period_cor4,DLBS_BAatlas_Period_cor3)
    
  }
  color_DLBS_MSN<-c("#9BBCD5","#5D7C95","#3E5368")
  color_DLBS_MIND<-c("#9798CA","#5D5EA8","#2F326A")
  
  color_GBM_MSN<-c("#E65E6A","#BC2531","#76181C")
  color_GBM_MIND<-c("#E7994C","#BD7124","#885224")
  
  if(j==1){
    color_MSN<-color_DLBS_MSN
    color_MIND<-color_DLBS_MIND
  }else{
    color_MSN<-color_GBM_MSN
    color_MIND<-color_GBM_MIND
    
  }
  DLBS_BAatlas_Period_cor4_1<-DLBS_BAatlas_Period_cor4[DLBS_BAatlas_Period_cor4$Net=="MSN", ]
  pdf(paste(list0[j], "_allregion_MSN_box.pdf",sep = "") ,width=15,height=5)  
  p1<-ggplot(DLBS_BAatlas_Period_cor4_1,aes(x=reorder(region2,-R),y=as.numeric(R) ,fill=Age))+ 
    stat_boxplot(geom = "errorbar")+
    geom_boxplot(color="azure4",outlier.colour="red",
                 outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
    )+  
    theme_classic()+  
    scale_fill_manual(values = color_MSN)+ 
    ylim(c(0,1))+
    
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="region", y="R")+ 
    theme(axis.title.x =element_text(size=12,face = "bold"), axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size = 8,face = "bold",angle = 45))
  print(p1)
  dev.off()
  
  
  DLBS_BAatlas_Period_cor4_2<-DLBS_BAatlas_Period_cor4[DLBS_BAatlas_Period_cor4$Net=="MIND", ]
  pdf(paste(list0[j], "_allregion_MIND_box.pdf",sep = "") ,width=15,height=5)  
  p2<-ggplot(DLBS_BAatlas_Period_cor4_2,aes(x=reorder(region2,-R),y=as.numeric(R) ,fill=Age))+ 
    stat_boxplot(geom = "errorbar")+
    geom_boxplot(color="azure4",outlier.colour="red",
                 outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
    )+  
    theme_classic()+  
    scale_fill_manual(values = color_MIND)+ 
    ylim(c(0,1))+
    
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="region", y="R")+ 
    theme(axis.title.x =element_text(size=12,face = "bold"), axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size = 8,face = "bold",angle = 45))
  print(p2)
  
  dev.off()
  
  
  
  
  regionlist<-c("1L","2L","3L","4L","9L","11L",
                "17L","20L","22L","24L","32L","33L",
                "40L","41L","44L","45L","46L",
                "1R","2R","3R","4R","9R","11R",
                "17R","20R","22R","24R","32R","33R",
                "40R","41R","44R","45R","46R")
  
  DLBS_BAatlas_Period_cor5<-DLBS_BAatlas_Period_cor4[unlist(sapply(regionlist, function(x){which(DLBS_BAatlas_Period_cor4$region2==x)}) ), ]
  DLBS_BAatlas_Period_cor5$region_16<-rep(c("S1C(1L)","S1C(2L)","S1C(3L)","M1C(4L)","DFC(9L)","OFC(11L)",
                                            "V1C(17L)","ITC(20L)","STC(22L)","MFC(24L)","MFC(32L)",
                                            "IPC(40L)","A1C(41L)","VFC(44L)","VFC(45L)",
                                            "S1C(1R)","S1C(2R)","S1C(3R)","M1C(4R)","DFC(9R)","OFC(11R)",
                                            "V1C(17R)","ITC(20R)","STC(22R)","MFC(24R)","MFC(32R)",
                                            "IPC(40R)","A1C(41R)","VFC(44R)","VFC(45R)"),each = 384)
  DLBS_BAatlas_Period_cor5$region_16_1<-substr(DLBS_BAatlas_Period_cor5$region_16,1,3)
  DLBS_BAatlas_Period_cor5_1<-DLBS_BAatlas_Period_cor5[DLBS_BAatlas_Period_cor5$Net=="MSN", ]
  DLBS_BAatlas_Period_cor5_2<-DLBS_BAatlas_Period_cor5[DLBS_BAatlas_Period_cor5$Net=="MIND", ]
  
  
  pdf(paste(list0[j], "_16region_MSN_box.pdf",sep = "") ,width=15,height=5)  
  p3<-ggplot(DLBS_BAatlas_Period_cor5_1,aes(x=reorder(region_16_1,-R),y=as.numeric(R) ,fill=Age))+ 
    stat_boxplot(geom = "errorbar")+
    geom_boxplot(color="azure4",outlier.colour="red",
                 outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
    )+  
    theme_classic()+  
    ylim(c(0,1))+
    stat_compare_means( method="anova")+
    
    scale_fill_manual(values = color_MSN)+ 
    ylim(c(0,1))+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="region", y="R")+ 
    theme(axis.title.x =element_text(size=12,face = "bold"), axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size = 8,face = "bold",angle = 45))
  print(p3)
  dev.off()
  
  pdf(paste(list0[j], "_16region_MIND_box.pdf",sep = "") ,width=15,height=5)  
  p4<-ggplot(DLBS_BAatlas_Period_cor5_2,aes(x=reorder(region_16_1,-R),y=as.numeric(R) ,fill=Age))+ 
    stat_boxplot(geom = "errorbar")+
    geom_boxplot(color="azure4",outlier.colour="red",
                 outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
    )+  
    theme_classic()+  
    scale_fill_manual(values = color_MIND)+ 
    ylim(c(0,1))+
    stat_compare_means( method="anova")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="region", y="R")+ 
    theme(axis.title.x =element_text(size=12,face = "bold"), axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size = 8,face = "bold",angle = 45))
  print(p4)
  dev.off()
  
  write.table(DLBS_BAatlas_Period_cor5_1,paste("D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\6\\",list0[j],"_16region_MSN_box.txt",sep="") ,col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
  write.table(DLBS_BAatlas_Period_cor5_2,paste("D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\6\\",list0[j],"_16region_MIND_box.txt",sep="") ,col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
  
  
}




rm(list=ls())
#install.packages("tidyverse")
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(readxl)
library(tidyverse)
library(ggsignif)
list0<-c("DKTatlas","HOAatlas","BAatlas")#Destrieux--

for (m in 1:3) {
  Periodlist<-c("Period_13","Period_14","Period_15","Group")
  Periodname<-c("Young","Middle","Late","All")
  Period<-c()
  brain2<-c()
  for (k in 1:4) {
    case_control<-c("GBM","DLBS")
    brain1<-c()
    for (n in 1:2) {
      netlist<-c("MSN","MIND")
      net1<-c()
      for (h in 1:2) {
        if(k!=4){
          Period_k<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\",netlist[h],"\\",list0[m],"\\",case_control[n],"_",list0[m],"_",Periodlist[k],"_cor.txt",sep = "") , sep = "\t", header = F,stringsAsFactors = F,check.names=F)
        }else{
          Period_k<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\",netlist[h],"\\",list0[m],"\\","Group_cor.txt",sep = "") , sep = "\t", header = F,stringsAsFactors = F,check.names=F)
        }
        Period_k1<-cbind(colnames(Period_k),Period_k)
        Period_k12<-cbind(melt(as.data.frame(Period_k1), id.vars = c("colnames(Period_k)") ),paste(netlist[h],Periodname[k],sep = "_"), case_control[n])
        colnames(Period_k12)<-c("Area1","Area2","R","Group","case_control")
        Period<-rbind(Period_k12,Period)
        if(m==1){
          leftbrain<-Period_k[1:31,1:31]
          rightbrain<-Period_k[32:62,32:62]
          interbrain<-Period_k[1:62,32:62]
        }else if(m==2){
          
          leftbrain<-Period_k[1:39,1:39]
          rightbrain<-Period_k[40:73,40:73]
          interbrain<-Period_k[1:39,40:73]
        }else{
          leftbrain<-Period_k[1:32,1:32]
          rightbrain<-Period_k[33:64,33:64]
          interbrain<-Period_k[1:32,33:64]
        }
        
        leftbrain1<-cbind(rownames(leftbrain),leftbrain)
        leftbrain12<-cbind(melt(as.data.frame(leftbrain1), id.vars = c("rownames(leftbrain)") ),"leftbrain")
        colnames(leftbrain12)<-c("Area1","Area2","R","Group")
        
        rightbrain1<-cbind(rownames(rightbrain),rightbrain)
        rightbrain12<-cbind(melt(as.data.frame(rightbrain1), id.vars = c("rownames(rightbrain)") ),"rightbrain")
        colnames(rightbrain12)<-c("Area1","Area2","R","Group")
        
        interbrain1<-cbind(rownames(interbrain),interbrain)
        interbrain12<-cbind(melt(as.data.frame(interbrain1), id.vars = c("rownames(interbrain)") ),"interbrain")
        colnames(interbrain12)<-c("Area1","Area2","R","Group")
        
        brain<-rbind(leftbrain12,rightbrain12,interbrain12)
        brain<-cbind(brain,paste(netlist[h],case_control[n],sep = "_") ,Periodname[k])
        colnames(brain)<-c("Area1","Area2","R","Group","case_control","Period")
        brain1<-rbind(brain,brain1)  
        
      }  
      
      colnames(brain1)<-c("Area1","Area2","R","Group","case_control","Period")
      brain2<-rbind(brain2,brain1)
    }   
    
  }
  #write.table(Period, paste("D:\\aging\\sample_stats\\",list0[m],"_Period.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
  
  Period_GBM<-Period[Period$case_control=="GBM" & Period$Group!="MIND_All" & Period$Group!="MSN_All",]
  p1<- ggplot(Period_GBM,aes(x=R,fill=Group))+geom_density(position="stack")+
    scale_fill_manual(values=c("MSN_Young" = "#f4606e","MSN_Middle"= "#c62232", "MSN_Late" = "#7a171c",
                               "MIND_Young" = "#ef9d4b","MIND_Middle"= "#c47322", "MIND_Late" = "#8c5218"))+
    #scale_fill_brewer(palette = "Reds")
    guides(fill=guide_legend(title=NULL)) +theme_bw()+theme(legend.position = c(0.2,0.8))+
    labs( title=paste(case_control[n],"_",list0[m],"_cor",sep = ""))
  ggsave(p1, file=paste("D:\\aging\\RESULT\\","GBM_",list0[m],"_cor.pdf",sep = ""), width=12, height=10)
  write.table(Period_GBM,paste("D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\6\\",list0[m],"_F.txt",sep = ""),col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
  
  
  Period_DLBS<-Period[Period$case_control=="DLBS" & Period$Group!="MIND_All" & Period$Group!="MSN_All",]
  p2<- ggplot(Period_DLBS,aes(x=R,fill=Group))+geom_density(position="stack")+
    scale_fill_manual(values=c("MSN_Young" = "#9dc1db","MSN_Middle"= "#5e7f9a", "MSN_Late" = "#3f5469",
                               "MIND_Young" = "#9b9bd3","MIND_Middle"= "#6262b5", "MIND_Late" = "#33336d"))+
    #scale_fill_brewer(palette = "blue")
    guides(fill=guide_legend(title=NULL)) +theme_bw()+theme(legend.position = c(0.2,0.8))+
    labs( title=paste(case_control[n],"_",list0[m],"_cor",sep = ""))
  ggsave(p2, file=paste("D:\\aging\\RESULT\\","DLBS_",list0[m],"_cor.pdf",sep = ""), width=12, height=10) 
  
  write.table(Period_DLBS,paste("D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\6\\",list0[m],"_E.txt",sep = "") ,col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
  
  
}

