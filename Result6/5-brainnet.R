rm(list=ls())
library(tidyverse)       
library(ggsci)           
library(ggExtra)         
library(ggpmisc)        
library(palmerpenguins) 
library(ggpubr)
library(reshape2)
Node_Brodmann64<-read.table("D:\\aging\\sample_stats\\Node_Brodmann64.txt",header=F,sep = "\t", quote = "")

DLBS_BAatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\DLBS\\MSN\\BAatlas\\DLBS_BAatlas_Period_13_cor.txt",header=F,sep = "\t", quote = "")
DLBS_BAatlas_Period_13_cor1<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_13_cor)
colnames(DLBS_BAatlas_Period_13_cor1)<-c("regions",Node_Brodmann64$V6)
DLBS_BAatlas_Period_13_cor2 <- melt(DLBS_BAatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_BAatlas_Period_13_cor2$value<-(DLBS_BAatlas_Period_13_cor2$value+1)/2

DLBS_BAatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\DLBS\\MIND\\BAatlas\\DLBS_BAatlas_Period_13_cor.txt",header=F,sep = "\t", quote = "")
DLBS_BAatlas_Period_13_cor1_mind<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_13_cor_mind)
colnames(DLBS_BAatlas_Period_13_cor1_mind)<-c("regions",Node_Brodmann64$V6)
DLBS_BAatlas_Period_13_cor2_mind <- melt(DLBS_BAatlas_Period_13_cor1_mind,id.vars = c("regions"))

data<-cbind(DLBS_BAatlas_Period_13_cor2,DLBS_BAatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"



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


DLBS_BAatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\DLBS\\MSN\\BAatlas\\DLBS_BAatlas_Period_14_cor.txt",header=F,sep = "\t", quote = "")
DLBS_BAatlas_Period_13_cor1<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_13_cor)
colnames(DLBS_BAatlas_Period_13_cor1)<-c("regions",Node_Brodmann64$V6)
DLBS_BAatlas_Period_13_cor2 <- melt(DLBS_BAatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_BAatlas_Period_13_cor2$value<-(DLBS_BAatlas_Period_13_cor2$value+1)/2

DLBS_BAatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\DLBS\\MIND\\BAatlas\\DLBS_BAatlas_Period_14_cor.txt",header=F,sep = "\t", quote = "")
DLBS_BAatlas_Period_13_cor1_mind<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_13_cor_mind)
colnames(DLBS_BAatlas_Period_13_cor1_mind)<-c("regions",Node_Brodmann64$V6)
DLBS_BAatlas_Period_13_cor2_mind <- melt(DLBS_BAatlas_Period_13_cor1_mind,id.vars = c("regions"))

data<-cbind(DLBS_BAatlas_Period_13_cor2,DLBS_BAatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"



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



DLBS_BAatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\DLBS\\MSN\\BAatlas\\DLBS_BAatlas_Period_15_cor.txt",header=F,sep = "\t", quote = "")
DLBS_BAatlas_Period_13_cor1<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_13_cor)
colnames(DLBS_BAatlas_Period_13_cor1)<-c("regions",Node_Brodmann64$V6)
DLBS_BAatlas_Period_13_cor2 <- melt(DLBS_BAatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_BAatlas_Period_13_cor2$value<-(DLBS_BAatlas_Period_13_cor2$value+1)/2

DLBS_BAatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\DLBS\\MIND\\BAatlas\\DLBS_BAatlas_Period_15_cor.txt",header=F,sep = "\t", quote = "")
DLBS_BAatlas_Period_13_cor1_mind<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_13_cor_mind)
colnames(DLBS_BAatlas_Period_13_cor1_mind)<-c("regions",Node_Brodmann64$V6)
DLBS_BAatlas_Period_13_cor2_mind <- melt(DLBS_BAatlas_Period_13_cor1_mind,id.vars = c("regions"))

data<-cbind(DLBS_BAatlas_Period_13_cor2,DLBS_BAatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"



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



################ gbm #################
DLBS_BAatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\GBM\\MSN\\BAatlas\\GBM_BAatlas_Period_13_cor.txt",header=F,sep = "\t", quote = "")
DLBS_BAatlas_Period_13_cor1<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_13_cor)
colnames(DLBS_BAatlas_Period_13_cor1)<-c("regions",Node_Brodmann64$V6)
DLBS_BAatlas_Period_13_cor2 <- melt(DLBS_BAatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_BAatlas_Period_13_cor2$value<-(DLBS_BAatlas_Period_13_cor2$value+1)/2


DLBS_BAatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\GBM\\MIND\\BAatlas\\GBM_BAatlas_Period_13_cor.txt",header=F,sep = "\t", quote = "")
DLBS_BAatlas_Period_13_cor1_mind<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_13_cor_mind)
colnames(DLBS_BAatlas_Period_13_cor1_mind)<-c("regions",Node_Brodmann64$V6)
DLBS_BAatlas_Period_13_cor2_mind <- melt(DLBS_BAatlas_Period_13_cor1_mind,id.vars = c("regions"))




data<-cbind(DLBS_BAatlas_Period_13_cor2,DLBS_BAatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"



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


DLBS_BAatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\GBM\\MSN\\BAatlas\\GBM_BAatlas_Period_14_cor.txt",header=F,sep = "\t", quote = "")
DLBS_BAatlas_Period_13_cor1<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_13_cor)
colnames(DLBS_BAatlas_Period_13_cor1)<-c("regions",Node_Brodmann64$V6)
DLBS_BAatlas_Period_13_cor2 <- melt(DLBS_BAatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_BAatlas_Period_13_cor2$value<-(DLBS_BAatlas_Period_13_cor2$value+1)/2


DLBS_BAatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\GBM\\MIND\\BAatlas\\GBM_BAatlas_Period_14_cor.txt",header=F,sep = "\t", quote = "")
DLBS_BAatlas_Period_13_cor1_mind<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_13_cor_mind)
colnames(DLBS_BAatlas_Period_13_cor1_mind)<-c("regions",Node_Brodmann64$V6)
DLBS_BAatlas_Period_13_cor2_mind <- melt(DLBS_BAatlas_Period_13_cor1_mind,id.vars = c("regions"))

data<-cbind(DLBS_BAatlas_Period_13_cor2,DLBS_BAatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"



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


DLBS_BAatlas_Period_13_cor<-read.table("D:\\aging\\sample_stats\\GBM\\MSN\\BAatlas\\GBM_BAatlas_Period_15_cor.txt",header=F,sep = "\t", quote = "")
DLBS_BAatlas_Period_13_cor1<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_13_cor)
colnames(DLBS_BAatlas_Period_13_cor1)<-c("regions",Node_Brodmann64$V6)
DLBS_BAatlas_Period_13_cor2 <- melt(DLBS_BAatlas_Period_13_cor1,id.vars = c("regions"))
DLBS_BAatlas_Period_13_cor2$value<-(DLBS_BAatlas_Period_13_cor2$value+1)/2


DLBS_BAatlas_Period_13_cor_mind<-read.table("D:\\aging\\sample_stats\\GBM\\MIND\\BAatlas\\GBM_BAatlas_Period_15_cor.txt",header=F,sep = "\t", quote = "")
DLBS_BAatlas_Period_13_cor1_mind<-cbind(Node_Brodmann64$V6,DLBS_BAatlas_Period_13_cor_mind)
colnames(DLBS_BAatlas_Period_13_cor1_mind)<-c("regions",Node_Brodmann64$V6)
DLBS_BAatlas_Period_13_cor2_mind <- melt(DLBS_BAatlas_Period_13_cor1_mind,id.vars = c("regions"))
data<-cbind(DLBS_BAatlas_Period_13_cor2,DLBS_BAatlas_Period_13_cor2_mind)
colnames(data)<-c("from_MSN","to_MSN","r_MSN","from_MIND","to_MIND","r_MIND")
data$group="Inconsistent"
data$group[which((data$r_MSN>=0.5 & data$r_MIND>=0.5) | (data$r_MSN<=0.5 & data$r_MIND<=0.5))]<-"Consistent"



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
  
  
}





