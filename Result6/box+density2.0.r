
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
  
  
  Period_DLBS<-Period[Period$case_control=="DLBS" & Period$Group!="MIND_All" & Period$Group!="MSN_All",]
  p2<- ggplot(Period_DLBS,aes(x=R,fill=Group))+geom_density(position="stack")+
    scale_fill_manual(values=c("MSN_Young" = "#9dc1db","MSN_Middle"= "#5e7f9a", "MSN_Late" = "#3f5469",
                               "MIND_Young" = "#9b9bd3","MIND_Middle"= "#6262b5", "MIND_Late" = "#33336d"))+
    #scale_fill_brewer(palette = "blue")
    guides(fill=guide_legend(title=NULL)) +theme_bw()+theme(legend.position = c(0.2,0.8))+
    labs( title=paste(case_control[n],"_",list0[m],"_cor",sep = ""))
  ggsave(p2, file=paste("D:\\aging\\RESULT\\","DLBS_",list0[m],"_cor.pdf",sep = ""), width=12, height=10) 
  
  
  colnames(brain2)<-c("Area1","Area2","R","Group","case_control","Period")
  write.table(brain2, paste("D:\\aging\\sample_stats\\",list0[m],"_MIND_dif.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
  
  
  
  brain3<-brain2[brain2$case_control=="MIND_DLBS" | brain2$case_control=="MSN_DLBS",]
  #brain3<-brain2[brain2$case_control=="MIND_DLBS",]
  brain3<-brain3[brain3$Period!="All",]
  p3<-
    ggplot(brain3,aes(x=Period,y=R,fill=case_control))+ 
    stat_boxplot(geom = "errorbar")+
    geom_boxplot(color="azure4",outlier.colour="red",
                 outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
    )+  
    
    facet_wrap(~case_control)+
    #stat_compare_means(comparisons = my_comparisons, method="t.test")+
    #stat_compare_means() +
    #scale_fill_material_d()+
    theme_bw()+  
    ylim(c(0,1.2))+
    scale_x_discrete(limits=c("Young","Middle","Late"))+
    scale_fill_manual(values = c("#6262b5",
                                 "#264568"))+ 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Group", y="Similarity value")+ 
    theme(axis.text.x = element_text(angle = 45))+
    theme(axis.title.x =element_text(size=12,face = "bold"), axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size = 12,face = "bold"))
  ggplot_build(p3)$data[[1]] %>% select(x,ymin,ymax) -> errorbar.df
  ggsave(p3, file=paste("D:\\aging\\RESULT\\",list0[m],"_DLBS_F1.pdf",sep=""), width=12, height=10)
  
  
  brain4<-brain2[brain2$case_control=="MIND_GBM" | brain2$case_control=="MSN_GBM",]
  #brain4<-brain2[brain2$case_control=="MSN_GBM",]
  brain4<-brain4[brain4$Period!="All",]
  p4<-
    ggplot(brain4,aes(x=Period,y=R,fill=case_control))+ 
    stat_boxplot(geom = "errorbar")+
    geom_boxplot(color="azure4",outlier.colour="red",
                 outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
    )+  
    
    facet_wrap(~case_control)+
    #stat_compare_means(comparisons = my_comparisons, method="t.test")+
    #stat_compare_means() +
    #scale_fill_material_d()+
    theme_bw()+  
    ylim(c(0,1.2))+
    scale_x_discrete(limits=c("Young","Middle","Late"))+
    scale_fill_manual(values = c("#bd7124",
                                 "#a64423"))+ 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Group", y="Similarity value")+ 
    theme(axis.text.x = element_text(angle = 45))+
    theme(axis.title.x =element_text(size=12,face = "bold"), axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size = 12,face = "bold"))
  ggplot_build(p4)$data[[1]] %>% select(x,ymin,ymax) -> errorbar.df
  ggsave(p4, file=paste("D:\\aging\\RESULT\\",list0[m],"_GBM_F1.pdf",sep=""), width=12, height=10)
  
  
  
  
  ggsave(p4, file=paste("D:\\aging\\RESULT\\",list0[m],"_GBM_F1.pdf",sep=""), width=12, height=10)
  
  
  
  
}