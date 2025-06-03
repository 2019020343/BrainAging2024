rm(list = ls())
library(irr)
################## os ###############
library("survival")
library("survminer")
osdata<-read.table("D:\\aging\\data\\GBM-NIHMS746836-supplement-78.txt" , sep = "\t", header = T,stringsAsFactors = F)
AverageThickness_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\lh_AverageThickness_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)
group_feature<-merge(osdata,AverageThickness_BAatlas,by.x = "Case",by.y = "lh.aparc.BAatlas.thickness")
osdata$group="Young"
osdata$group[osdata$Age>=40]="Middle"
osdata$group[osdata$Age>=60]="Late"

osdata_out<-as.data.frame(cbind(group_feature$Case,group_feature$Age,group_feature$group,group_feature$Survival..months.,group_feature$Vital.status..1.dead.)) 
osdata_out$days<-as.numeric(osdata_out$V4) *30
colnames(osdata_out)<-c("Case","Age","group","months","statu","days")

sfit <- survfit(Surv(as.numeric(Survival..months.) , as.numeric(Vital.status..1.dead.))~group, data=osdata) 
ggsurvplot(sfit, conf.int=F, pval=TRUE)
ggsurvplot(sfit,data = osdata,
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           # submain = coxdata$cerna[1], 
           legend = "top", 
           legend.title = "group",
           palette = c("#1f2a6b", "#446799","#a3ced6"))

write.table(osdata,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF8\\A.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


#################### CAS #######################

GBM_gene<-read.table("D:\\aging\\GBM\\HT_HG-U133A.txt",header=T,sep = "\t", quote = "",stringsAsFactors = F,check.names=FALSE)
case<-GBM_gene[,-grep("\\-11$",colnames(GBM_gene))]

osdata<-read.table("D:\\aging\\data\\GBM-NIHMS746836-supplement-78.txt" , sep = "\t", header = T,stringsAsFactors = F)
osdata_age<-osdata[,c(1,15,16)]
CAS_genes<-read.table("D:\\aging\\RESULT\\2-3-CAS\\B_247genes.txt",header=T,sep = "\t", quote = "")
CAS_genes_EXP<-merge(case,CAS_genes,by.x="sample",by.y="Var1")[,1:ncol(case)]


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
options(mc.cores = 10)
vis <- analyze(vis)

tsne <- getProjections(vis)[["tSNE10"]]
sigScores <- getSignatureScores(vis)[, "Common Aging Score"]
sigScores1<-cbind(names(sigScores),sigScores)
colnames(sigScores1)<-c("Sample_geo_accession1","sigScores")
sigScores1<-as.data.frame(sigScores1)
sigScores1$Case<-sapply(sigScores1$Sample_geo_accession1, function(x){substr(x,1,12) } )
sigScores_group<-merge(sigScores1,osdata_age,by= "Case" )
write.table(sigScores_group,"D:\\aging\\RESULT\\6radiomic\\sigScores_group.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

sigScores_group$group1=70
sigScores_group$group1[sigScores_group$Age < 70]="60"
sigScores_group$group1[sigScores_group$Age < 60]="50"
sigScores_group$group1[sigScores_group$Age < 50]="40"
sigScores_group$group1[sigScores_group$Age < 40]="30"
sigScores_group$group1[sigScores_group$Age < 30]="20"
write.table(sigScores_group,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF8\\C.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


P3<-ggplot(sigScores_group,aes(x=as.character(group1) ,y=as.numeric(sigScores) ,fill=as.character(group1)))+ 
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(color="azure4",outlier.colour="red",
               outlier.fill="red",outlier.size=1,outlier.alpha=0)+  
  stat_compare_means(method = "anova" )+
  theme_classic()+  
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

P3<-ggplot(sigScores_group,aes(x=group1,y=as.numeric(sigScores) ,fill=Gender))+ 
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(color="azure4",outlier.colour="red",
               outlier.fill="red",outlier.size=1,outlier.alpha=0#notch=TRUE,notchwidth = 0.8
  )+  
  theme_classic()+  
  scale_fill_manual(values = c("#D8996D",
                               "#77977A"))+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Group", y="CAS")+ 
  theme(axis.title.x =element_text(size=12,face = "bold"), axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size = 12,face = "bold"))+
  #stat_summary(fun.y = "mean", geom = "point", size = 0.5) +
  stat_summary(fun.y = "mean", geom = "line", 
               aes(group = Gender), 
               size = 0.5)

########### ICC##############
sampledata<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\T1_sample_feature\\TCGA-02-0085_result.txt" , sep = "\t", header = T,stringsAsFactors = F)
output<-c()
list<-c("mean","median","max","min","std")
for (j in 1:5) {
  data1=read.table(paste("E:\\Mask_RCNN_master\\pyradiomics\\T1_test\\F-S_matrix_",list[j],".txt",sep = "")  , sep = "\t", header = T,stringsAsFactors = F)
  for (k in 1:5) {
    data2=read.table(paste("E:\\Mask_RCNN_master\\pyradiomics\\T1_test\\F-S_matrix_",list[k],".txt",sep = "")  , sep = "\t", header = T,stringsAsFactors = F)
    for (i in 1:nrow(sampledata)) {
      featuredata<-t(rbind(data1[i,],data2[i,]))
      icc_value<-icc(featuredata, 
                     model = "twoway", 
                     type = "consistency", 
                     unit = "single")$value
      out1<-cbind(list[j],list[k], sampledata[i,1],icc_value)
      output<-rbind(output,out1)
      
    }
  }

}
output<-as.data.frame(output)
output1<-output[-which(output$V1==output$V2),]

library(ggplot2)
output1$group<-"neg"
output1$group[as.numeric(output1$icc_value)>=0]<-"0"
output1$group[as.numeric(output1$icc_value)>=0.75]<-"75"
output1$group[as.numeric(output1$icc_value)>=0.9]<-"9"

write.table(output1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF8\\H.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

ggplot(output1, aes(V1,as.numeric(icc_value))) +
  geom_jitter(aes(color = group), size = 0.001)+
  ggpubr::color_palette(palette  =c( "#E5C9CB","#C48388", "#921D22","#C4C5C6"))+
  ggpubr::theme_pubclean()

d9<-as.data.frame(unique(output1[output1$group=="9" & output1$V1=="mean",3])) 
colnames(d9)<-"features"
write.table(d9, "D:\\aging\\RESULT\\6radiomic\\icc9.txt", sep = '\t', col.names = T, quote = FALSE,row.names = F)

#################### feature selection ##############
rm(list = ls())
library(pROC)
library(randomForest)
library(rfPermute)
library(ggplot2)

d9<-read.table("D:\\aging\\RESULT\\6radiomic\\icc9.txt" , sep = "\t", header = T,stringsAsFactors = F)
mean_data=read.table("E:\\Mask_RCNN_master\\pyradiomics\\T1_test\\F-S_matrix_mean.txt", sep = "\t", header = T,stringsAsFactors = F)
sampledata<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\T1_sample_feature\\TCGA-02-0085_result.txt" , sep = "\t", header = T,stringsAsFactors = F)
mean_data1<-cbind(sampledata[,1],mean_data)
mean_icc<-merge(d9,mean_data1,by.x = "features",by.y = "sampledata[, 1]")
mean_icc1<-as.data.frame(t(mean_icc))
colnames(mean_icc1)<-mean_icc1[1,]
mean_icc1<-mean_icc1[-1,]

mean_icc_data<- matrix(as.numeric(as.matrix(mean_icc1) ),ncol = ncol(mean_icc1))

colnames(mean_icc_data)<-colnames(mean_icc1)
rownames(mean_icc_data)<-rownames(mean_icc1)

info=read.table("D:\\aging\\RESULT\\6radiomic\\GBM-info.txt", sep = "\t", header = T,stringsAsFactors = F)
sample_group<-as.data.frame(cbind(info$Case,info$Age,info$Gender))  
sample_group$group="0"
sample_group$group[ sample_group$V2>=40 & sample_group$V2<60]="1"
sample_group$group[sample_group$V2>=60]="2"

set.seed(1234)
feature_success_all<-c()
pout<-c()
aucout<-c()
seeds<-sample(1:1000000,1000,replace = F)
for (k in 1:1000) {
  set.seed(seeds[k])
  
  index<-sample(1:nrow(mean_icc_data),nrow(mean_icc_data)*0.5,replace = F)
  train_sample<-mean_icc_data[index,]
  samplename1<-sapply(rownames(train_sample),function(x){gsub('[.]','-',x)})
  train_info<-sample_group[sapply(samplename1,function(x){grep(x,sample_group$V1)}),]
  train_data<-cbind(train_sample,as.numeric(train_info$V2) )
  train_data_group<-cbind(train_sample,train_info$group )
  
  
  test_sample<-mean_icc_data[-index,]
  samplename2<-sapply(rownames(test_sample),function(x){gsub('[.]','-',x)})
  test_info<-sample_group[sapply(samplename2,function(x){grep(x,sample_group$V1)}),]
  #test_data<-cbind(test_sample,test_info[,-1])
  test_data<-cbind(test_sample,test_info$group )


  set.seed(123)
  otu_rfP <- rfPermute(train_data[,-785], train_data[,785], importance = TRUE, ntree = 500, nrep = 1000, num.cores = 1)

  importance_otu.scale <- data.frame(importance(otu_rfP, scale = TRUE), check.names = FALSE)
  importance_otu.scale$OTU_name <- rownames(importance_otu.scale)
  importance_otu.scale$OTU_name <- factor(importance_otu.scale$OTU_name, levels = importance_otu.scale$OTU_name)
  importance_otu.scale_sig<-importance_otu.scale[importance_otu.scale$`%IncMSE.pval`<=0.05,]
  if(nrow(importance_otu.scale_sig)>0){
    train_data1<-as.data.frame(train_data_group[,c(sapply(importance_otu.scale_sig$OTU_name,function(x){which(colnames(train_data_group)==x)}),785)]) 
    test_data1<-as.data.frame(test_data[,c(sapply(importance_otu.scale_sig$OTU_name,function(x){which(colnames(test_data)==x)}),785)]) 
    train_data2<- cbind(matrix(as.numeric(as.matrix(train_data1) ),ncol = ncol(train_data1)) ) 
    test_data2<- cbind(matrix(as.numeric(as.matrix(test_data1) ),ncol = ncol(test_data1)) ) 
    
    set.seed(123)
    rf_train<-randomForest(train_data2[,-ncol(train_data2)], train_data2[,ncol(train_data2)],
                           mtry=11,ntree=500,importance = T,proximity=TRUE)
    
    testpred <- predict(rf_train,newdata = test_data2[,-ncol(test_data2)])
    auc<-roc(test_data2[,ncol(test_data2)],as.numeric(as.matrix(testpred)))$auc
    if(auc>0.5){
      feature_success<-as.character(importance_otu.scale_sig$OTU_name) 
      feature_success_all<-c(feature_success_all,feature_success)
    }
  }
  #pout<-rbind(pout,importance_otu.scale[,c(5,2)])
  aucout1<-cbind(k,as.character(importance_otu.scale_sig[,5] ) ,importance_otu.scale_sig[,1],importance_otu.scale_sig[,2], auc)
  colnames(aucout1)<-c("perturbation",colnames(importance_otu.scale_sig)[c(5,1,2)],"auc")
  aucout<-rbind(aucout,aucout1)
  print(k)
}
feature_fre<-as.data.frame(table(feature_success_all))
write.table(feature_fre,"D:\\aging\\RESULT\\6radiomic\\feature_fre1_100.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
write.table(pout,"D:\\aging\\RESULT\\6radiomic\\pout1_100.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
write.table(aucout,"D:\\aging\\RESULT\\6radiomic\\aucout1_100.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



aucout<-read.table("D:\\aging\\RESULT\\6radiomic\\aucout1_100.txt" , sep = "\t", header = T,stringsAsFactors = F)
data <- data.frame(var1 = aucout$auc)
ggplot(data, aes(x=x) ) +
  # Top
  geom_density(aes(x = var1, y = after_stat(density)), fill="#921D22" )+
  theme_classic()+
  labs(x="AUC",y="Density")

write.table(data,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF8\\H_RIGHT.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

##############  RF_8 云雨图 ###################
feature_fre_auc8<-read.table("D:\\aging\\RESULT\\6radiomic\\feature_fre_auc8.txt" , sep = "\t", header = T,stringsAsFactors = F)
feature_fre_auc8_table<-as.data.frame(table(feature_fre_auc8$OTU_name))
feature_8<-as.data.frame(feature_fre_auc8_table[order(feature_fre_auc8_table$Freq,decreasing = T)[1:8], ]) 
mean_data=read.table("E:\\Mask_RCNN_master\\pyradiomics\\T1_test\\F-S_matrix_mean.txt", sep = "\t", header = T,stringsAsFactors = F,check.names=FALSE)
sampledata<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\T1_sample_feature\\TCGA-02-0085_result.txt" , sep = "\t", header = T,stringsAsFactors = F,check.names=FALSE)
mean_data1<-cbind(sampledata[,1],mean_data)
mean_feature<-merge(feature_8,mean_data1,by.x = "Var1",by.y = "sampledata[, 1]")
mean_feature1<-as.data.frame(t(mean_feature))
colnames(mean_feature1)<-mean_feature1[1,]
mean_feature1<-mean_feature1[-c(1:2),]
mean_feature1<-cbind(rownames(mean_feature1),mean_feature1)
osdata<-read.table("D:\\aging\\data\\GBM-NIHMS746836-supplement-78.txt" , sep = "\t", header = T,stringsAsFactors = F)
agedata<-osdata[,c(1,15)]
agedata<-merge(agedata,mean_feature1,by.x = "Case",by.y = "rownames(mean_feature1)")
agedata$group="Young"
agedata$group[agedata$Age>=40]="Middle"
agedata$group[agedata$Age>=60]="Late"
write.table(agedata,"D:\\aging\\RESULT\\6radiomic\\agedata.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
library(ggplot2)
library(gplots)
library(ggpubr)
library(ggsignif)
library(ggdist)
feature_list<-colnames(agedata)[3:10]
df<-c()
index<-c()
p_lst <- list()
df_ALL<-c()
agedatai_all<-c()
for (j in 3:10) {
  agedata_Y<-agedata[agedata$group=="Young",j]
  agedata_M<-agedata[agedata$group=="Middle",j]
  agedata_L<-agedata[agedata$group=="Late",j]
  df<-rbind( cbind("Young",mean(as.numeric(agedata_Y) ),sd(as.numeric(agedata_Y)) ) ,
         cbind("Middle",mean(as.numeric(agedata_M) ),sd(as.numeric(agedata_M)) ) , 
         cbind("Late",mean(as.numeric(agedata_L) ),sd(as.numeric(agedata_L)) ) )
  df<-as.data.frame(df)
  colnames(df)<-c("group","mean","se")
  df$mean<-as.numeric(df$mean)
  df$se<-as.numeric(df$se)
  agedatai<-agedata[,c(j,11)]
 p1<- ggplot()+
    geom_errorbar(data=df,aes(x = group, ymin = mean - se, ymax = mean + se),width = 0,position = position_nudge(x=-0.2),color="gray2")+
    geom_point(data=df,mapping=aes(x=df$group,y=df$mean),size=3,pch=21,position = position_nudge(x=-0.2),color="gray2",fill="white")+
    geom_jitter(data=agedatai,mapping = aes(x=group,y=as.numeric(agedatai[,1]) ,color=group),width = 0.05,alpha=0.3,size=1)+ 
    stat_halfeye(data=agedatai,mapping = aes(x=group,y=as.numeric(agedatai[,1]),fill=group),width = 0.25, justification = -0.65,.width = 0,point_colour = NA)+ 
    stat_compare_means(method = "anova",label.y.npc = "top")+ 
    scale_color_manual(values=c("#1f2a6b","#446799","#a3ced6"))+
    scale_fill_manual(values=c("#1f2a6b","#446799","#a3ced6"))+
    #geom_hline(aes(yintercept = 2.1),linetype="dashed",color="gray")+
    labs(x="", y=colnames(agedata)[j] )+
    #scale_y_continuous(expand = c(0, 0), limit = c(1.3, 2.8))+
    #annotate("rect", xmin = 2.5, xmax =3.5,  ymin = 1.3, ymax = 2.8, alpha = 0.2,fill="gray58")+
    scale_x_discrete(limits=c("Young","Middle","Late"))+
    theme(axis.text=element_text(colour='black',size=9))+theme_bw()+
    theme(panel.grid = element_blank())+
   theme(legend.position="none")

 #p_lst[[(j-2)]] <- p1
 ggsave(paste("D:\\aging\\RESULT\\6radiomic\\",colnames(agedata)[j] ,".pdf",sep=""), width = 4, height = 4 )
 df1<-cbind(df,colnames(agedata)[j])
 colnames(df1)<-c("group","mean","se","RFs" )
 df_ALL<-rbind(df_ALL,df1)
 agedatai1<- cbind(agedatai,colnames(agedata)[j])
 colnames(agedatai1)<-c("values","group","RFs" )
 agedatai_all<-rbind(agedatai_all,agedatai1)
}

write.table(df_ALL,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF10\\A_1.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
write.table(agedatai_all,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF10\\A_2.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


############## RF_CF COR ###################
rm(list = ls())
library(reshape2)
agedata<-read.table("D:\\aging\\RESULT\\6radiomic\\agedata.txt" , sep = "\t", header = T,stringsAsFactors = F)

#1
lh_AverageThickness_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\lh_AverageThickness_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35)]
rh_AverageThickness_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\rh_AverageThickness_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35)]
AverageThickness_BAatlas<-cbind(lh_AverageThickness_BAatlas,rh_AverageThickness_BAatlas[,-1])
#2
lh_curvind_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\lh_curvind_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35)]
rh_curvind_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\rh_curvind_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35)]
curvind_BAatlas<-cbind(lh_curvind_BAatlas,rh_curvind_BAatlas[,-1])
#3
lh_foldind_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\lh_foldind_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35)]
rh_foldind_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\rh_foldind_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35)]
foldind_BAatlas<-cbind(lh_foldind_BAatlas,rh_foldind_BAatlas[,-1])
#4
lh_gauscurv_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\lh_gauscurv_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35)]
rh_gauscurv_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\rh_gauscurv_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35)]
gauscurv_BAatlas<-cbind(lh_gauscurv_BAatlas,rh_gauscurv_BAatlas[,-1])
#5
lh_GrayMattervolume_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\lh_GrayMattervolume_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35)]
rh_GrayMattervolume_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\rh_GrayMattervolume_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35)]
GrayMattervolume_BAatlas<-cbind(lh_GrayMattervolume_BAatlas,rh_GrayMattervolume_BAatlas[,-1])
#6
lh_meancurv_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\lh_meancurv_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35)]
rh_meancurv_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\rh_meancurv_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35)]
meancurv_BAatlas<-cbind(lh_meancurv_BAatlas,rh_meancurv_BAatlas[,-1])
#7
lh_surfarea_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\lh_surfarea_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35,36)]
rh_surfarea_BAatlas<-read.table("D:\\aging\\RESULT\\6radiomic\\feature7\\rh_surfarea_BAatlas.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35,36)]
surfarea_BAatlas<-cbind(lh_surfarea_BAatlas,rh_surfarea_BAatlas[,-1])


group_feature<-merge(agedata,AverageThickness_BAatlas,by.x = "Case",by.y = "lh.aparc.BAatlas.thickness")
group_list<-c("Young","Middle","Late")
ri_r_mat4<-c()
for (j in 1:3) {
  ri_r_mat2<-c()
  for (i in 2:ncol(surfarea_BAatlas)) {
  index<-sapply(group_feature$Case[group_feature$group==group_list[j]] ,function(x){which(surfarea_BAatlas[,1]==x)})
  data1<-AverageThickness_BAatlas[index,]
  data2<-curvind_BAatlas[index,]
  data3<-foldind_BAatlas[index,]
  data4<-gauscurv_BAatlas[index,]
  data5<-GrayMattervolume_BAatlas[index,]
  data6<-meancurv_BAatlas[index,]
  data7<-surfarea_BAatlas[index,]
  RFdata<-  group_feature[group_feature$group==group_list[j],]

   ri_r_mat<-matrix(0, nrow = 8, ncol = 7)
    for (k in 3:10) {
      RF_value<-as.numeric(RFdata[,k])
      ri_r_mat[k-2,1]<-cor.test(RF_value,data1[,i])[[4]]
      ri_r_mat[k-2,2]<-cor.test(RF_value,data2[,i])[[4]]
      ri_r_mat[k-2,3]<-cor.test(RF_value,data3[,i])[[4]]
      ri_r_mat[k-2,4]<-cor.test(RF_value,data4[,i])[[4]]
      ri_r_mat[k-2,5]<-cor.test(RF_value,data5[,i])[[4]]
      ri_r_mat[k-2,6]<-cor.test(RF_value,data6[,i])[[4]]
      ri_r_mat[k-2,7]<-cor.test(RF_value,data7[,i])[[4]]
      
    }
   ri_r_mat<-as.data.frame(ri_r_mat)
   ri_r_mat<-cbind(colnames(RFdata)[3:10],ri_r_mat)
   colnames(ri_r_mat)<-c("Features","AverageThickness","curvind","foldind","gauscurv","GrayMattervolume","meancurv","surfarea")
   ri_r_mat1 <- melt(ri_r_mat, id.vars = "Features")
   ri_r_mat1<-cbind(colnames(surfarea_BAatlas)[i],ri_r_mat1)
   ri_r_mat2<-rbind(ri_r_mat2,ri_r_mat1)
}

colnames(ri_r_mat2)<-c("region","RFeatures","CFeatures","value")
ri_r_mat3<-cbind(group_list[j],ri_r_mat2)
ri_r_mat4<-rbind(ri_r_mat4,ri_r_mat3)
}
colnames(ri_r_mat4)<-c("group","region","RFeatures","CFeatures","value")
write.table(ri_r_mat4,"D:\\aging\\RESULT\\6radiomic\\RF_CF_cor.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
# 使用aggregate函数和mean函数计算每个分组数据的均值
Youngresult <- cbind( "Young",aggregate(value ~ region, ri_r_mat4[ri_r_mat4$group=="Young",], mean))
Middleresult <- cbind( "Middle",aggregate(value ~ region, ri_r_mat4[ri_r_mat4$group=="Middle",], mean))
Lateresult <- cbind( "Late",aggregate(value ~ region, ri_r_mat4[ri_r_mat4$group=="Late",], mean))
write.table(Youngresult,"D:\\aging\\RESULT\\6radiomic\\Youngresult_cor.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
write.table(Middleresult,"D:\\aging\\RESULT\\6radiomic\\Middleresult_cor.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
write.table(Lateresult,"D:\\aging\\RESULT\\6radiomic\\Lateresult_cor.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


################# ggsegBrodmann ###################

rm(list = ls())
library(ggseg)
#> Warning: package 'ggseg' was built under R version 4.1.1
#> Loading required package: ggplot2
library(ggseg3d)
library(ggsegBrodmann)
ri_r_mat4<-read.table("D:\\aging\\RESULT\\6radiomic\\RF_CF_cor.txt" , sep = "\t", header = T,stringsAsFactors = F)[,-c(34,35,36)]
Youngresult1<-ri_r_mat4[ri_r_mat4$group=="Young",]
Middleresult1<-ri_r_mat4[ri_r_mat4$group=="Middle",]
Lateresult1<-ri_r_mat4[ri_r_mat4$group=="Late",]
rf_list<-unique(ri_r_mat4$RFeatures)
RF_DATA_all<-c()
for (i in 1:8) {
Young_RF <- cbind( rf_list[i],"Young",aggregate(value ~ region, Youngresult1[Youngresult1$RFeatures==rf_list[i],], mean))
Middle_RF <- cbind( rf_list[i],"Middle",aggregate(value ~ region, Middleresult1[Middleresult1$RFeatures==rf_list[i],], mean))
Late_RF <- cbind( rf_list[i],"Late",aggregate(value ~ region, Lateresult1[Lateresult1$RFeatures==rf_list[i],], mean))
colnames(Young_RF)<-c("RF","group","region","value")
colnames(Middle_RF)<-c("RF","group","region","value")
colnames(Late_RF)<-c("RF","group","region","value")
brodmann1<-as.data.frame(brodmann)
label1<-sapply(brodmann1$label, function(x){gsub("BA","",x) } )
brodmann1$label2<-sapply(label1, function(x){paste(x,"_area",sep = "") } )

lh_add<-cbind(rf_list[i],"Young",
   c("lh_1_3_area","lh_26_29_30_area","lh_41_42_52_area"),   
  rbind(mean(Young_RF$value[Young_RF$region=="lh_1_area" | Young_RF$region=="lh_3_area"]),
            mean(Young_RF$value[Young_RF$region=="lh_26_area" | Young_RF$region=="lh_29_area" |Young_RF$region=="lh_30_area"]),
            mean(Young_RF$value[Young_RF$region=="lh_41_area" | Young_RF$region=="lh_42_area"]))) 
 
rh_add<-cbind(rf_list[i],"Young",
              c("rh_1_3_area","rh_26_29_30_area","rh_41_42_52_area"),   
              rbind(mean(Young_RF$value[Young_RF$region=="rh_1_area" | Young_RF$region=="rh_3_area"]),
                    mean(Young_RF$value[Young_RF$region=="rh_26_area" | Young_RF$region=="rh_29_area" |Young_RF$region=="rh_30_area"]),
                    mean(Young_RF$value[Young_RF$region=="rh_41_area" | Young_RF$region=="rh_42_area"])))
colnames(lh_add)<-colnames(Young_RF)
colnames(rh_add)<-colnames(Young_RF)
Young_RF_BA <- rbind(Young_RF,lh_add,rh_add)
lh_add<-cbind(rf_list[i],"Middle",
              c("lh_1_3_area","lh_26_29_30_area","lh_41_42_52_area"),   
              rbind(mean(Middle_RF$value[Middle_RF$region=="lh_1_area" | Middle_RF$region=="lh_3_area"]),
                    mean(Middle_RF$value[Middle_RF$region=="lh_26_area" | Middle_RF$region=="lh_29_area" |Middle_RF$region=="lh_30_area"]),
                    mean(Middle_RF$value[Middle_RF$region=="lh_41_area" | Middle_RF$region=="lh_42_area"]))) 
rh_add<-cbind(rf_list[i],"Middle",
              c("rh_1_3_area","rh_26_29_30_area","rh_41_42_52_area"),   
              rbind(mean(Middle_RF$value[Middle_RF$region=="rh_1_area" | Middle_RF$region=="rh_3_area"]),
                    mean(Middle_RF$value[Middle_RF$region=="rh_26_area" | Middle_RF$region=="rh_29_area" |Middle_RF$region=="rh_30_area"]),
                    mean(Middle_RF$value[Middle_RF$region=="rh_41_area" | Middle_RF$region=="rh_42_area"])))
colnames(lh_add)<-colnames(Middle_RF)
colnames(rh_add)<-colnames(Middle_RF)
Middle_RF_BA <- rbind(Middle_RF,lh_add,rh_add)
lh_add<-cbind(rf_list[i],"Late",
              c("lh_1_3_area","lh_26_29_30_area","lh_41_42_52_area"),   
              rbind(mean(Late_RF$value[Late_RF$region=="lh_1_area" | Late_RF$region=="lh_3_area"]),
                    mean(Late_RF$value[Late_RF$region=="lh_26_area" | Late_RF$region=="lh_29_area" |Late_RF$region=="lh_30_area"]),
                    mean(Late_RF$value[Late_RF$region=="lh_41_area" | Late_RF$region=="lh_42_area"]))) 
rh_add<-cbind(rf_list[i],"Late",
              c("rh_1_3_area","rh_26_29_30_area","rh_41_42_52_area"),   
              rbind(mean(Late_RF$value[Late_RF$region=="rh_1_area" | Late_RF$region=="rh_3_area"]),
                    mean(Late_RF$value[Late_RF$region=="rh_26_area" | Late_RF$region=="rh_29_area" |Late_RF$region=="rh_30_area"]),
                    mean(Late_RF$value[Late_RF$region=="rh_41_area" | Late_RF$region=="rh_42_area"])))
colnames(lh_add)<-colnames(Late_RF)
colnames(rh_add)<-colnames(Late_RF)
Late_RF_BA <- rbind(Late_RF,lh_add,rh_add)

Youngresult_BA1<-merge(Young_RF_BA,brodmann1,by.x = "region",by.y = "label2")
Middleresult_BA1<-merge(Middle_RF_BA,brodmann1,by.x = "region",by.y = "label2")
Lateresult_BA1<-merge(Late_RF_BA,brodmann1,by.x = "region",by.y = "label2")
Youngresult_BA2<-Youngresult_BA1[,c(1:5,7,8,10,11)]
Middleresult_BA2<-Middleresult_BA1[,c(1:5,7,8,10,11)]
Lateresult_BA2<-Lateresult_BA1[,c(1:5,7,8,10,11)]
RF_BA2<-rbind(Youngresult_BA2,Middleresult_BA2,Lateresult_BA2)
#write.table(RF_BA2,paste("D:\\aging\\RESULT\\9online tool\\3\\3-2\\",rf_list[i],"_BA.txt",sep="" ) ,col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
cols<-colorRampPalette(c("#3A264D", "#537F7C", "#FFFF06"))(200)
p1<-Youngresult_BA1 %>%
  ggplot() +
  geom_brain(atlas = brodmann, 
             position = position_brain(hemi ~ side),
             aes(fill = as.numeric(value) )) +
  theme_classic()+
  scale_fill_gradientn(limits=c(-1, 1), colours = cols)

p2<-Middleresult_BA1 %>%
  ggplot() +
  geom_brain(atlas = brodmann, 
             position = position_brain(hemi ~ side),
             aes(fill = as.numeric(value) )) +
  theme_classic()+
  scale_fill_gradientn(limits=c(-1, 1), colours = cols)

p3<-Lateresult_BA1 %>%
  ggplot() +
  geom_brain(atlas = brodmann, 
             position = position_brain(hemi ~ side),
             aes(fill = as.numeric(value) )) +
  theme_classic()+
  scale_fill_gradientn(limits=c(-1, 1), colours = cols)

p1+p2+p3

#ggsave(paste("D:\\aging\\RESULT\\6radiomic\\",rf_list[i] ,".pdf",sep=""), width = 15, height = 15 )
RF_DATA<-rbind(Youngresult_BA1,Middleresult_BA1,Lateresult_BA1)
RF_DATA_all<-rbind(RF_DATA_all,RF_DATA )
}
RF_DATA_all<-RF_DATA_all[,-12]
write.table(RF_DATA_all,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\7\\B_Bottom.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


########### GSVA ###############
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)  #install.packages("msigdbr")
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)
GBM_gene<-read.table("D:\\aging\\GBM\\HT_HG-U133A.txt",header=T,sep = "\t", quote = "",stringsAsFactors = F,check.names=FALSE)
case<-GBM_gene[,-grep("\\-11$",colnames(GBM_gene))]
agedata<-read.table("D:\\aging\\RESULT\\6radiomic\\agedata.txt" , sep = "\t", header = T,stringsAsFactors = F)
case_RF<-as.data.frame(case[,unlist(sapply(agedata$Case, function(x){grep(unlist(strsplit(x,"-"))[3]  ,colnames(case))  }))])

KEGG_df_all <-  msigdbr(species = "Homo sapiens", # Homo sapiens or Mus musculus
                        category = "C2",
                        subcategory = "CP:KEGG") 
KEGG_df <- dplyr::select(KEGG_df_all,gs_name,gs_exact_source,gene_symbol)
kegg_list <- split(KEGG_df$gene_symbol, KEGG_df$gs_name) ##按照gs_name给gene_symbol分组

##GO
GO_df_all <- msigdbr(species = "Homo sapiens",
                     category = "C5")  
GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
GO_df <- GO_df[GO_df$gs_subcat!="HPO",]
go_list <- split(GO_df$gene_symbol, GO_df$gs_name) ##按照gs_name给gene_symbol分组
dat<-as.matrix(case_RF)
rownames(dat)<-case[,1]
gsva_mat_go <- gsva(expr=dat, 
                 gset.idx.list=go_list, 
                 kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T, 
                 parallel.sz = parallel::detectCores())#调用所有核
gsva_mat_kegg <- gsva(expr=dat, 
                 gset.idx.list=kegg_list, 
                 kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T, 
                 parallel.sz = parallel::detectCores())#调用所有核

for (k in 1:2) {
  if(k==1){
    gsva_mat<-gsva_mat_go
  }else{
    gsva_mat<-gsva_mat_kegg
  }

  Young<-agedata$Case[which(agedata$group=="Young")]
  Middle<-agedata$Case[which(agedata$group=="Middle")]
  Late<-agedata$Case[which(agedata$group=="Late")]
  
  Young_index<-as.data.frame(gsva_mat[,unlist(sapply(Young, function(x){grep(unlist(strsplit(x,"-"))[3]  ,colnames(gsva_mat))  }))])
  Middle_index<-as.data.frame(gsva_mat[,unlist(sapply(Middle, function(x){grep(unlist(strsplit(x,"-"))[3]  ,colnames(gsva_mat))  }))])
  Late_index<-as.data.frame(gsva_mat[,unlist(sapply(Late, function(x){grep(unlist(strsplit(x,"-"))[3]  ,colnames(gsva_mat))  }))])
  
  re_out<-as.data.frame(rownames(gsva_mat))
  colnames(re_out)<-"rownames(nrDEG)"
  for (m in 1:3) {
    if(m==1){
      exp_cy <- cbind(Young_index, Middle_index)
      group <- c(rep('Young', ncol(Young_index)), rep('Middle', ncol(Middle_index)))
      group <- factor(group, levels = c("Young", "Middle"))
      group <- as.matrix(group)
      design <- model.matrix(~0+group)
      colnames(design)=levels(factor(group))
      rownames(design)=rownames(group)
      contrast.matrix<-makeContrasts(paste(c("Young","Middle"),collapse = "-"),levels = design)
    }else if(m==2){
      exp_cy <- cbind(Young_index, Late_index)
      group <- c(rep('Young', ncol(Young_index)), rep('Late', ncol(Late_index)))
      group <- factor(group, levels = c("Young", "Late"))
      group <- as.matrix(group)
      design <- model.matrix(~0+group)
      colnames(design)=levels(factor(group))
      rownames(design)=rownames(group)
      contrast.matrix<-makeContrasts(paste(c("Young","Late"),collapse = "-"),levels = design)
    }else if(m==3){
      exp_cy <- cbind(Middle_index, Late_index)
      group <- c(rep('Middle', ncol(Middle_index)), rep('Late', ncol(Late_index)))
      group <- factor(group, levels = c("Middle", "Late"))
      group <- as.matrix(group)
      design <- model.matrix(~0+group)
      colnames(design)=levels(factor(group))
      rownames(design)=rownames(group)
      contrast.matrix<-makeContrasts(paste(c("Middle","Late"),collapse = "-"),levels = design)
    }
    
    mydata<-matrix(as.numeric(unlist(exp_cy) ),ncol = ncol(exp_cy)) 
    rownames(mydata)<-rownames(gsva_mat) 
    ##step1
    fit <- lmFit(mydata,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput) 
    result1<-cbind(rownames(nrDEG),nrDEG[,c(1,4,5)])
    re_out<-merge(result1,re_out,by = "rownames(nrDEG)")
  }
  colnames(re_out)<-c("gene","YM_logFC","YM_P.Value","YM_adj.P.Val","YL_logFC","YL_P.Value","YL_adj.P.Val","ML_logFC","ML_P.Value","ML_adj.P.Val")
  adjpdata<-re_out[,c(3,6,9)]
  result <-re_out[which( apply(adjpdata,1,function(x){sum(x <=0.01)>=1 })),]
  if(k==1){
    result_GO<-result
    result_GO_show<-result[order(apply(result[,c(2,5,8)],1,sum),decreasing = T)[1:40] ,]
    result_GO_show_R9<-result[order(apply(result[,c(2,5,8)],1,sum),decreasing = T),]
    gsva_mat_go_orderage<-cbind(Young_index,Middle_index,Late_index)
  }else{
    result_KEGG<-result
    result_KEGG_show<-result[order(apply(result[,c(2,5,8)],1,sum),decreasing = T)[1:10] ,]
    gsva_mat_kegg_orderage<-cbind(Young_index,Middle_index,Late_index)
  }
}

result_show<-rbind(result_GO_show,result_KEGG_show)
gsva_mat_show<-rbind(gsva_mat_go_orderage,gsva_mat_kegg_orderage)
gsva_mat_show<-cbind(rownames(gsva_mat_show),gsva_mat_show)
colnames(gsva_mat_show)<-c("gene",colnames(gsva_mat_show)[-1])
gsva_mat_show_out<-merge(gsva_mat_show,result_show,by="gene") [,1:42]
rownames(gsva_mat_show_out)=gsva_mat_show_out[,1]
gsva_mat_show_out<-gsva_mat_show_out[,-1]
gsva_mat_show_out1<-matrix(as.numeric(as.matrix(gsva_mat_show_out) ),ncol = ncol(gsva_mat_show_out) )
colnames(gsva_mat_show_out1)<-colnames(gsva_mat_show_out)
rownames(gsva_mat_show_out1)<-tolower(rownames(gsva_mat_show_out))
library(pheatmap)
# 根据指定的样本顺序对数据进行重新排列
paletteLength <- 2000
myColor <- colorRampPalette(c("#154182", "#FFFFFF", "#99131D"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
annotation_col<-data.frame( c(rep("Young",length(Young)),
                              rep("Middle",length(Middle)),
                              rep("Late",length(Late)))) 
colnames(annotation_col)<-"group"
rownames(annotation_col)<-colnames(gsva_mat_show_out1)

annotation_row<-data.frame( c(rep("GO",40),
                              rep("KEGG",10))) 
colnames(annotation_row)<-"label"
rownames(annotation_row)<-rownames(gsva_mat_show_out1)
ann_colors=list(group=c(Young="#a3ced6",Middle="#446799",Late="#1f2a6b"),
                label=c(GO="#6F5A90",KEGG="#EBBF2A"))

pheatmap(gsva_mat_show_out1, 
         scale = "row",
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         cluster_rows = T,
         cluster_cols = F,
         #cellwidth = 20, cellheight = 15,
         color = myColor,
         annotation_colors = ann_colors,
         #breaks=myBreaks, 
         show_rownames = F, show_colnames = F,
         annotation_legend = T, 
         gaps_col = c(7,23))
result_all<-rbind(result_KEGG,result_GO)
gsva_mat_all<-merge(gsva_mat_show,result_all,by="gene") [,1:42]
write.table(gsva_mat_show_out1,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF8\\D.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


corout2<-c()
for (i in 3:10) {
  corout1<-c()
  for (j in 1:nrow(gsva_mat_all)) {
    r<-cor.test(as.numeric(agedata[,i]),as.numeric(gsva_mat_all[j,-1]) )[[4]]
    p<-cor.test(as.numeric(agedata[,i]),as.numeric(gsva_mat_all[j,-1]) )[[3]]
    corout<-cbind(as.character(colnames(agedata)[i]) ,as.character(gsva_mat_all[j,1] ),r,p )
    corout1<-rbind(corout1,corout)
    
    
  }
  corout2<-rbind(corout2,corout1)
}
corout2<-as.data.frame(corout2)
corout2_sig<-corout2[corout2$p<0.01,]
corout2_sig_bp<-corout2_sig[grep("GOBP",corout2_sig$V2),  ]
corout2_sig_kegg<-corout2_sig[grep("KEGG",corout2_sig$V2),  ]
corout2_sig<-rbind(corout2_sig_bp,corout2_sig_kegg)
corout2_sig$V2<-tolower(corout2_sig$V2)
node<-rbind(cbind(unique(corout2_sig$V1),"RF"),cbind(unique(corout2_sig$V2),"GOBP"))
write.table(corout2_sig,"D:\\aging\\RESULT\\6radiomic\\corout2_sig.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
write.table(node,"D:\\aging\\RESULT\\6radiomic\\node.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
GOBP<-read.table("D:\\aging\\RESULT\\6radiomic\\GOBP.txt",header=T,sep = "\t", quote = "",stringsAsFactors = F,check.names=FALSE)
GOBP$V1<-toupper(GOBP$V2)
result_show_r9<-merge(GOBP,result_GO_show_R9,by.x="V1",by.y="gene")
write.table(result_show_r9,"D:\\aging\\RESULT\\9online tool\\3\\3-2\\Table1_go_limma.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
GOBP_corout2_sig<-merge(GOBP,corout2_sig,by="V2")
write.table(GOBP_corout2_sig,"D:\\aging\\RESULT\\6radiomic\\corout2_sig.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



rm(list = ls())
library(reshape2)
library(ggplot2)
corout2_sig<-read.table("D:\\aging\\RESULT\\6radiomic\\corout2_sig.txt",header=T,sep = "\t", quote = "",stringsAsFactors = F,check.names=FALSE)
cols<-colorRampPalette(c("#45436E","#D8C7D0", "#9F4E4C"))(200)
corout2_F2<-corout2_sig[c(corout2_sig$V1=="wavelet.LH_firstorder_Variance" |corout2_sig$V1=="logarithm_firstorder_Mean"), ]

df3<-data.frame(
  x=seq(1.5,7.5,1),
  xend=seq(1.5,7.5,1),
  y=-Inf,
  yend=Inf)

df4<-data.frame(x=-Inf,
                xend=Inf,
                y=seq(1.5,45.5,1),
                yend=seq(1.5,45.5,1)
                )




ggplot(data = corout2_sig,aes(x=V1,y=V2_2))+
  geom_point(aes(color=as.numeric(r),size=0.5),
                 shape=15)+
  scale_color_gradientn(limits=c(-1, 1), colours = cols)+
  theme_bw()+
  scale_x_discrete(limits=c("logarithm_firstorder_Mean","wavelet.LH_firstorder_Variance","exponential_firstorder_10Percentile","wavelet.LL_firstorder_10Percentile","logarithm_firstorder_10Percentile",
                            "original_firstorder_10Percentile","square_firstorder_10Percentile","squareroot_firstorder_10Percentile"))+
  scale_y_discrete(limits=c("GOBP9","GOBP6","GOBP31","GOBP28","GOBP21","GOBP24","GOBP25","GOBP26",
                            "GOBP5","GOBP34","GOBP33","GOBP32","GOBP30","GOBP29","GOBP27","GOBP23",
                            "GOBP13","GOBP20","GOBP8","GOBP7","GOBP45","GOBP44","GOBP43","GOBP42",
                            "GOBP41","GOBP40","GOBP39","GOBP38","GOBP37","GOBP22","GOBP2","GOBP19",
                            "GOBP14","GOBP10","GOBP11","GOBP16","GOBP3","GOBP4","GOBP12","GOBP15",
                            "GOBP17","GOBP18","GOBP1","GOBP35","GOBP36","GOBP46"))+
  
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey"),
        axis.ticks = element_blank())+
  geom_segment(data = df3,aes(x=x,xend=xend,y=y,yend=yend),color="grey")+
  geom_segment(data = df4,aes(x=x,xend=xend,y=y,yend=yend),color="grey")+ 
  #theme(axis.text.x = element_text(angle = 45))+ 
  coord_fixed(ratio=0.8)+
 xlab("")+
  ylab("")

write.table(corout2_sig,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\7\\C.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

GOBP<-read.table("D:\\aging\\RESULT\\6radiomic\\GOBP.txt",header=T,sep = "\t", quote = "",stringsAsFactors = F,check.names=FALSE)
GOBP$V2<-toupper(GOBP$V2) 
GO_df<-read.table("D:\\aging\\RESULT\\6radiomic\\GO_df.txt",header=T,sep = "\t", quote = "",stringsAsFactors = F,check.names=FALSE)
GO_df1<-unique(GO_df[,c(1,3)]) 
GOBP1<-merge(GOBP,GO_df1,by.x="V2",by.y="gs_name")
write.table(GOBP1,"D:\\aging\\RESULT\\6radiomic\\GO_lable_term_ID.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

########### xcell #############

library(xCell)
GBM_gene<-read.table("D:\\aging\\GBM\\HT_HG-U133A.txt",header=T,sep = "\t", quote = "",stringsAsFactors = F,check.names=FALSE)
case<-GBM_gene[,-grep("\\-11$",colnames(GBM_gene))]
expmean_gene1<-case[,-1]
rownames(expmean_gene1)<- case[,1]
sigScores_group<-read.table("D:\\aging\\RESULT\\6radiomic\\sigScores_group.txt",header=T,sep = "\t", quote = "")

scores <-  xCellAnalysis(expmean_gene1,rnaseq = T)
xCell_scores<-as.data.frame(t(scores))
xCell_scores1<-cbind(substr(rownames(xCell_scores), 1,12) ,xCell_scores)
xCell_group<-merge(xCell_scores1,sigScores_group,by.x = "substr(rownames(xCell_scores), 1, 12)",by.y ="Case" )
library(tidyverse)       
library(ggsci)           
library(ggExtra)         
library(ggpmisc)        
library(palmerpenguins) 
library(ggpubr)
for (i in 2:68) {
  ggplot(xCell_group, aes(sigScores, xCell_group[,i])) +
    geom_point() +  # 添加散点图层，点的大小表示体重
    labs(x = "CAS", y = colnames(xCell_group)[i]) +  # 设置坐标轴标签
    theme_classic()  +
    stat_cor(geom = "text",label.x = 0)+
    theme(axis.title = element_text(family = "sans", 
                                    face = "bold"))
  
  
  ggsave(paste("D:\\aging\\RESULT\\6radiomic\\xcell\\",colnames(xCell_group)[i],".pdf",sep="") )
  
}

xCell_scores2 <- melt(xCell_scores1,id.vars = c("substr(rownames(xCell_scores), 1, 12)"))
xCell_group<-merge(xCell_scores2,sigScores_group,by.x = "substr(rownames(xCell_scores), 1, 12)",by.y ="Case" )
xCell_group$group[xCell_group$Age<12]="Childhood"
xCell_group$group[xCell_group$Age>=12 & xCell_group$Age<20]="Adolescence"
xCell_group$group[xCell_group$Age>=20 & xCell_group$Age<40]="Young"
xCell_group$group[xCell_group$Age>=40 & xCell_group$Age<60]="Middle"
xCell_group$group[xCell_group$Age>=60]="Late"
write.table(xCell_group,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF9\\A-C.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

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
  ggsave(paste("D:\\aging\\RESULT\\6radiomic\\xcell\\age_p\\",celltype[i],".pdf",sep=""), width = 20, height = 20, units = "cm" )
}

rm(list = ls())
library(irr)
################## os ###############
library("survival")
library("survminer")
osdata<-read.table("D:\\aging\\data\\GBM-NIHMS746836-supplement-78.txt" , sep = "\t", header = T,stringsAsFactors = F)
coordinate<-read.table("D:\\aging\\投稿\\1-NC\\修稿1\\lobe\\TCGA_120_en.txt" , sep = "\t", header = T,stringsAsFactors = F)
group_feature<-merge(osdata,coordinate,by.x = "Case",by.y = "sample")


extract_ggplotdata<-as.data.frame(table(group_feature$lobe_4C) ) 
ggplot(extract_ggplotdata, aes(x =reorder(Var1,-Freq) , y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack")+
  # geom_bar(position=position_dodge(), stat="identity",width = 0.8, color = NA, size = 0.5)+
  scale_fill_manual(values = c("#c5b5d1","#12803b","#5A7EB3","#cf1223"))+ 
  #scale_y_continuous(expand = c(0,0)) +
  #ylim(-6,15)+
  theme_classic()+
  labs(x="",y="# of samples")+　
  theme(axis.text.x = element_text(angle = 45))



group_feature<-merge(osdata,coordinate,by.x = "Case",by.y = "sample")
group_feature<-group_feature[c(group_feature$lobe_4C=="frontal" | group_feature$lobe_4C=="occipital"), ]
osdata_out<-as.data.frame(cbind(group_feature$Case,
                                group_feature$lobe_3C,
                                group_feature$Survival..months.,
                                group_feature$Vital.status..1.dead.)) 
colnames(osdata_out)<-c("Case","lobe","days","status")

sfit <- survfit(Surv(as.numeric(days) , as.numeric(status))~lobe, data=osdata_out) 
#ggsurvplot(sfit, conf.int=F, pval=TRUE)

ggsurvplot(sfit,data = osdata_out,
           pval = TRUE, #conf.int = TRUE,
           #risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           # submain = coxdata$cerna[1], 
           legend = "top", 
           legend.title = "lobe",
           palette = c("#c5b5d1", "#12803b","#cacbd0"))
write.table(osdata_out,"D:\\aging\\投稿\\1-NC\\修稿1\\fig_data\\SF8\\B.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
