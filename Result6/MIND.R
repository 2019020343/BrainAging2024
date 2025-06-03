############### DLBS age group MIND #############
rm(list=ls())
library(reshape2)
library(tidyr)
case_control<-c("GBM","DLBS")
GBM_infor<-read.table("D:\\aging\\data\\GBM-NIHMS746836-supplement-78.txt", sep = "\t", header = T,stringsAsFactors = F,check.names=F)
DLBS_infor<-read.table("D:\\aging\\data\\DLBS\\subjects_information.txt", sep = "\t", header = T,stringsAsFactors = F,check.names=F)

for (n in 1:2) {
  if(n==1){
    infor<-GBM_infor
  }else{
    infor<-DLBS_infor
  }
  list0<-c("DKTatlas","HOAatlas","BAatlas")#Destrieux--
  for (m in 1:2) {
    Period_13<-infor[infor$Age<40,1]
    Period_14<-infor[infor$Age>=40 & infor$Age<60,1]
    Period_15<-infor[infor$Age>=60,1]
    filelist<-list.files(paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],sep="") )
    sample1<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],"\\",filelist[1],sep=""), sep = "\t", header = F,stringsAsFactors = F,check.names=F)
    
    Gi_r_13<-matrix(0,nrow = ncol(sample1),ncol = ncol(sample1))
    for (i in 1:length(Period_13)) {
      filename<-filelist[grep(Period_13[i],filelist)]  
      if(length(filename)!=0){
        si_r_13<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],"\\",filename,sep=""), sep = "\t", header = F,stringsAsFactors = F,check.names=F)
        Gi_r_13<-Gi_r_13+as.matrix(si_r_13) 
      }
      
    }
    G_r_13<-Gi_r_13/length(Period_13)
    write.table(G_r_13, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],"\\",case_control[n],"_",list0[m],"_Period_13_cor.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    
    threshold<-round(nrow(G_r_13)*nrow(G_r_13)*0.1)
    G_r_131<-cbind(colnames(G_r_13),G_r_13)
    G_r_132<-melt(as.data.frame(G_r_131), id.vars = c("V1") )
    value<-G_r_132[order(abs(as.numeric(G_r_132[,3])) ,decreasing = T)[threshold],3]
    G_r_131[G_r_131<value]<-0
    G_r_131<-G_r_131[,-1]
    write.table(G_r_131, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],"\\",case_control[n],"_",list0[m],"_Period_13_10.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    
    
    
    Gi_r_14<-matrix(0,nrow = ncol(sample1),ncol = ncol(sample1))
    for (i in 1:length(Period_14)) {
      filename<-filelist[grep(Period_14[i],filelist)]  
      if(length(filename)!=0){
        si_r_14<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],"\\",filename,sep=""), sep = "\t", header = F,stringsAsFactors = F,check.names=F)
        Gi_r_14<-Gi_r_14+as.matrix(si_r_14) 
      }}
    G_r_14<-Gi_r_14/length(Period_14)
    write.table(G_r_14, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],"\\",case_control[n],"_",list0[m],"_Period_14_cor.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    
    threshold<-round(nrow(G_r_14)*nrow(G_r_14)*0.1)
    G_r_141<-cbind(colnames(G_r_14),G_r_14)
    G_r_142<-melt(as.data.frame(G_r_141), id.vars = c("V1") )
    value<-G_r_142[order(abs(as.numeric(G_r_142[,3])) ,decreasing = T)[threshold],3]
    G_r_141[G_r_141<value]<-0
    G_r_141<-G_r_141[,-1]
    write.table(G_r_141, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],"\\",case_control[n],"_",list0[m],"_Period_14_10.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    
    
    
    Gi_r_15<-matrix(0,nrow = ncol(sample1),ncol = ncol(sample1))
    for (i in 1:length(Period_15)) {
      filename<-filelist[grep(Period_15[i],filelist)]  
      if(length(filename)!=0){
        si_r_15<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],"\\",filename,sep=""), sep = "\t", header = F,stringsAsFactors = F,check.names=F)
        Gi_r_15<-Gi_r_15+as.matrix(si_r_15) 
      }}
    G_r_15<-Gi_r_15/length(Period_15)
    write.table(G_r_15, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],"\\",case_control[n],"_",list0[m],"_Period_15_cor.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    
    threshold<-round(nrow(G_r_15)*nrow(G_r_15)*0.1)
    G_r_151<-cbind(colnames(G_r_15),G_r_15)
    G_r_152<-melt(as.data.frame(G_r_151), id.vars = c("V1") )
    value<-G_r_152[order(abs(as.numeric(G_r_152[,3])) ,decreasing = T)[threshold],3]
    G_r_151[G_r_151<value]<-0
    G_r_151<-G_r_151[,-1]
    write.table(G_r_151, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],"\\",case_control[n],"_",list0[m],"_Period_15_10.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    print(paste(case_control[n],list0[m],sep = "_"))
    
    
    
    Group_cor<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],"\\Group_cor.txt",sep=""), sep = "\t", header = F,stringsAsFactors = F,check.names=F)
    threshold<-round(nrow(Group_cor)*nrow(Group_cor)*0.1)
    Group_cor1<-cbind(colnames(Group_cor),Group_cor)
    Group_cor2<-melt(as.data.frame(Group_cor1), id.vars = c("colnames(Group_cor)") )
    value<-Group_cor2[order(abs(as.numeric(Group_cor2[,3])) ,decreasing = T)[threshold],3]
    Group_cor1[Group_cor1<value]<-0
    Group_cor1<-Group_cor1[,-1]
    write.table(Group_cor1, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MIND\\",list0[m],"\\",case_control[n],"_",list0[m],"_Group_10.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    
    
  }
  
}
