rm(list=ls())
case_control<-c("GBM","DLBS")
for (n in 1:2) {
  list0<-c("a2009s","DKTatlas","HOAatlas","BAatlas")#Destrieux--
  for (m in 1:4) {
    lh_surfarea<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\lh_surfarea_",list0[m],".txt",sep = "") , sep = "\t", header = T,stringsAsFactors = F)
    lh_curvind<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\lh_curvind_",list0[m],".txt",sep = ""), sep = "\t", header = T,stringsAsFactors = F)
    lh_foldind<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\lh_foldind_",list0[m],".txt",sep = ""), sep = "\t", header = T,stringsAsFactors = F)
    lh_gauscurv<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\lh_gauscurv_",list0[m],".txt",sep = ""),sep = "\t", header = T,stringsAsFactors = F)
    lh_meancurv<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\lh_meancurv_",list0[m],".txt",sep = ""), sep = "\t", header = T,stringsAsFactors = F)
    lh_AverageThickness<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\lh_AverageThickness_",list0[m],".txt",sep = ""),sep = "\t", header = T,stringsAsFactors = F)
    lh_GrayMattervolume<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\lh_GrayMattervolume_",list0[m],".txt",sep = ""), sep = "\t", header = T,stringsAsFactors = F)
    
    
    rh_surfarea<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\rh_surfarea_",list0[m],".txt",sep = ""),sep = "\t", header = T,stringsAsFactors = F)
    rh_curvind<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\rh_curvind_",list0[m],".txt",sep = ""),  sep = "\t", header = T,stringsAsFactors = F)
    rh_foldind<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\rh_foldind_",list0[m],".txt",sep = ""), sep = "\t", header = T,stringsAsFactors = F)
    rh_gauscurv<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\rh_gauscurv_",list0[m],".txt",sep = ""),sep = "\t", header = T,stringsAsFactors = F)
    rh_meancurv<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\rh_meancurv_",list0[m],".txt",sep = ""),sep = "\t", header = T,stringsAsFactors = F)
    rh_AverageThickness<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\rh_AverageThickness_",list0[m],".txt",sep = ""), sep = "\t", header = T,stringsAsFactors = F)
    rh_GrayMattervolume<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\rh_GrayMattervolume_",list0[m],".txt",sep = ""), sep = "\t", header = T,stringsAsFactors = F)
    
    Gi_r<-matrix(0,nrow =((ncol(lh_surfarea)-4)+(ncol(rh_surfarea)-4)),ncol =((ncol(lh_surfarea)-4)+(ncol(rh_surfarea)-4)))
    
    for (i in 1:nrow(lh_surfarea)) {
      lh_sa_i<-lh_surfarea[i,2:(ncol(lh_surfarea)-3)]
      lh_curvind_i<-lh_curvind[i,2:(ncol(lh_surfarea)-3)]
      lh_foldind_i<-lh_foldind[i,2:(ncol(lh_surfarea)-3)]
      lh_gauscurv_i<-lh_gauscurv[i,2:(ncol(lh_surfarea)-3)]
      lh_meancurv_i<-lh_meancurv[i,2:(ncol(lh_surfarea)-3)]  
      lh_AverageThickness_i<-lh_AverageThickness[i,2:(ncol(lh_surfarea)-3)]
      lh_GrayMattervolume_i<-lh_GrayMattervolume[i,2:(ncol(lh_surfarea)-3)]
      
      rh_sa_i<-rh_surfarea[i,2:(ncol(rh_surfarea)-3)]
      rh_curvind_i<-rh_curvind[i,2:(ncol(rh_surfarea)-3)]
      rh_foldind_i<-rh_foldind[i,2:(ncol(rh_surfarea)-3)]
      rh_gauscurv_i<-rh_gauscurv[i,2:(ncol(rh_surfarea)-3)]
      rh_meancurv_i<-rh_meancurv[i,2:(ncol(rh_surfarea)-3)]
      rh_AverageThickness_i<-rh_AverageThickness[i,2:(ncol(rh_surfarea)-3)]
      rh_GrayMattervolume_i<-rh_GrayMattervolume[i,2:(ncol(rh_surfarea)-3)]
      if(m==1){
        lh_colname<-c("lh_G_and_S_frontomargin",	"lh_G_and_S_occipital_inf",	"lh_G_and_S_paracentral",	"lh_G_and_S_subcentral",	"lh_G_and_S_transv_frontopol",
                      "lh_G_and_S_cingul.Ant",	"lh_G_and_S_cingul.Mid.Ant",	"lh_G_and_S_cingul.Mid.Post",	"lh_G_cingul.Post.dorsal","lh_G_cingul_Post_ventral","lh_G_cuneus",
                      "lh_G_front_inf.Opercular",	"lh_G_front_inf.Orbital",	"lh_G_front_inf.Triangul",	"lh_G_front_middle",	"lh_G_front_sup",
                      "lh_G_Ins_lg_and_S_cent_ins",	"lh_G_insular_short",	"lh_G_occipital_middle",	"lh_G_occipital_sup",	"lh_G_oc.temp_lat.fusifor",
                      "lh_G_oc.temp_med.Lingual",	"lh_G_oc.temp_med.Parahip",	"lh_G_orbital",	"lh_G_pariet_inf.Angular",	"lh_G_pariet_inf.Supramar",
                      "lh_G_parietal_sup",	"lh_G_postcentral",	"lh_G_precentral",	"lh_G_precuneus",	"lh_G_rectus",	"lh_G_subcallosal",	"lh_G_temp_sup.G_T_transv",
                      "lh_G_temp_sup.Lateral",	"lh_G_temp_sup.Plan_polar",	"lh_G_temp_sup.Plan_tempo",	"lh_G_temporal_inf",	"lh_G_temporal_middle",
                      "lh_Lat_Fis.ant.Horizont",	"lh_Lat_Fis.ant.Vertical",	"lh_Lat_Fis.post",	"lh_Pole_occipital",	"lh_Pole_temporal",	"lh_S_calcarine",
                      "lh_S_central",	"lh_S_cingul.Marginalis",	"lh_S_circular_insula_ant",	"lh_S_circular_insula_inf",	"lh_S_circular_insula_sup",	"lh_S_collat_transv_ant",
                      "lh_S_collat_transv_post",	"lh_S_front_inf",	"lh_S_front_middle",	"lh_S_front_sup",	"lh_S_interm_prim.Jensen",	"lh_S_intrapariet_and_P_trans",
                      "lh_S_oc_middle_and_Lunatus",	"lh_S_oc_sup_and_transversal",	"lh_S_occipital_ant",	"lh_S_oc.temp_lat",	"lh_S_oc.temp_med_and_Lingual",	"lh_S_orbital_lateral",
                      "lh_S_orbital_med.olfact",	"lh_S_orbital.H_Shaped",	"lh_S_parieto_occipital","lh_S_pericallosal","lh_S_postcentral",	"lh_S_precentral.inf.part",	"lh_S_precentral.sup.part",
                      "lh_S_suborbital",	"lh_S_subparietal",	"lh_S_temporal_inf",	"lh_S_temporal_sup",	"lh_S_temporal_transverse")
        
        
        rh_colname<-c("rh_G_and_S_frontomargin",	"rh_G_and_S_occipital_inf",	"rh_G_and_S_paracentral",	"rh_G_and_S_subcentral",	"rh_G_and_S_transv_frontopol",
                      "rh_G_and_S_cingul.Ant",	"rh_G_and_S_cingul.Mid.Ant",	"rh_G_and_S_cingul.Mid.Post",	"rh_G_cingul.Post.dorsal","rh_G_cingul_Post_ventral",	"rh_G_cuneus",
                      "rh_G_front_inf.Opercular",	"rh_G_front_inf.Orbital",	"rh_G_front_inf.Triangul",	"rh_G_front_middle",	"rh_G_front_sup",
                      "rh_G_Ins_lg_and_S_cent_ins",	"rh_G_insular_short",	"rh_G_occipital_middle",	"rh_G_occipital_sup",	"rh_G_oc.temp_lat.fusifor",
                      "rh_G_oc.temp_med.Lingual",	"rh_G_oc.temp_med.Parahip",	"rh_G_orbital",	"rh_G_pariet_inf.Angular",	"rh_G_pariet_inf.Supramar",
                      "rh_G_parietal_sup",	"rh_G_postcentral",	"rh_G_precentral",	"rh_G_precuneus",	"rh_G_rectus",	"rh_G_subcallosal",	"rh_G_temp_sup.G_T_transv",
                      "rh_G_temp_sup.Lateral",	"rh_G_temp_sup.Plan_polar",	"rh_G_temp_sup.Plan_tempo",	"rh_G_temporal_inf",	"rh_G_temporal_middle",
                      "rh_Lat_Fis.ant.Horizont",	"rh_Lat_Fis.ant.Vertical",	"rh_Lat_Fis.post",	"rh_Pole_occipital",	"rh_Pole_temporal",	"rh_S_calcarine",
                      "rh_S_central",	"rh_S_cingul.Marginalis",	"rh_S_circular_insula_ant",	"rh_S_circular_insula_inf",	"rh_S_circular_insula_sup",	"rh_S_collat_transv_ant",
                      "rh_S_collat_transv_post",	"rh_S_front_inf",	"rh_S_front_middle",	"rh_S_front_sup",	"rh_S_interm_prim.Jensen",	"rh_S_intrapariet_and_P_trans",
                      "rh_S_oc_middle_and_Lunatus",	"rh_S_oc_sup_and_transversal",	"rh_S_occipital_ant",	"rh_S_oc.temp_lat",	"rh_S_oc.temp_med_and_Lingual",	"rh_S_orbital_lateral",
                      "rh_S_orbital_med.olfact",	"rh_S_orbital.H_Shaped",	"rh_S_parieto_occipital",	"rh_S_pericallosal","rh_S_postcentral",	"rh_S_precentral.inf.part",	"rh_S_precentral.sup.part",
                      "rh_S_suborbital","rh_S_subparietal",	"rh_S_temporal_inf",	"rh_S_temporal_sup",	"rh_S_temporal_transverse")
        
        
      }else if(m==2){
        
        lh_colname<-c("lh_caudalanteriorcingulate",	"lh_caudalmiddlefrontal",	
                      "lh_cuneus",	"lh_entorhinal",	"lh_fusiform",
                      "lh_inferiorparietal",	"lh_inferiortemporal",	
                      "lh_isthmuscingulate",	"lh_lateraloccipital",	"lh_lateralorbitofrontal",
                      "lh_lingual",	"lh_medialorbitofrontal",	"lh_middletemporal",	
                      "lh_parahippocampal","lh_paracentral",
                      "lh_parsopercularis",	"lh_parsorbitalis",	
                      "lh_parstriangularis",	"lh_pericalcarine",	"lh_postcentral",
                      "lh_posteriorcingulate",	"lh_precentral",	"lh_precuneus",	
                      "lh_rostralanteriorcingulate",	"lh_rostralmiddlefrontal",
                      "lh_superiorfrontal",	"lh_superiorparietal",	"lh_superiortemporal",	
                      "lh_supramarginal",	"lh_transversetemporal",	"lh_insula")
        
        
        rh_colname<-c("rh_caudalanteriorcingulate",	"rh_caudalmiddlefrontal",	
                      "rh_cuneus",	"rh_entorhinal",	"rh_fusiform",
                      "rh_inferiorparietal",	"rh_inferiortemporal",	"rh_isthmuscingulate",	
                      "rh_lateraloccipital",	"rh_lateralorbitofrontal",
                      "rh_lingual",	"rh_medialorbitofrontal",	"rh_middletemporal",	
                      "rh_parahippocampal",	"rh_paracentral",
                      "rh_parsopercularis",	"rh_parsorbitalis",	"rh_parstriangularis",	"rh_pericalcarine",
                      "rh_postcentral",
                      "rh_posteriorcingulate",	"rh_precentral",	"rh_precuneus",	
                      "rh_rostralanteriorcingulate",	"rh_rostralmiddlefrontal",
                      "rh_superiorfrontal",	"rh_superiorparietal",	"rh_superiortemporal",	
                      "rh_supramarginal",	"rh_transversetemporal",	"rh_insula")
      }else if(m==3){
        
        
        lh_colname<-c("FP.L",	"INS.L","F1.L",	"F2.L",	"F3t.L", "F3o.L","PRG.L","TP.L",	
                      "T1a.L","T1p.L","T2a.L","T2p.L","TO2.L","T3a.L", "T3p.L","TO3.L",
                      "POG.L","SPL.L", "SGa.L","SGp.L", "AG.L","OLs.L","OLi.L","CALC.L",	
                      "FMC.L","SMC.L","CGa.L","CN.L","TFp.L","TOF.L","OF.L","FO.L",	
                      "CO.L","PO.L","PP.L","H.L","Put.L","Pall.L","Hip.L")
        
        
        rh_colname<-c("FP.R","INS.R",	"F1.R","F2.R","F3t.R","F3o.R","PRG.R","T1a.R",	
                      "T1p.R","T2a.R", "T2p.R","TO2.R","T3a.R","T3p.R","TO3.R","POG.R",
                      "FMC.R","SMC.R", "SC.R","PAC.R", "CGa.R","CGp.R","PCN.R","CN.R",	
                      "FOC.R","TFa.R","TFp.R","TOF.R","OF.R","FO.R","CO.R","PO.R",	
                      "PP.R","H.R")
      }else{
        
        lh_colname<-c("1L",	"2L",	"3L",	"4L",	"5L",	
                      "6L",	"7L",	"8L",	"9L",	"10L",	
                      "11L",	"17L",	"18L",	"19L",	"20L",	
                      "21L",	"22L",	"23L",	"24L",	"25L",	
                      "26L",	"29L",	"30L",	"32L",	"38L",	
                      "39L",	"40L",	"41L",	"42L",	"43L",	"44L",	"45L")
        
        
        rh_colname<-c("1R",	"2R",	"3R",	"4R",	"5R",	
                      "6R",	"7R",	"8R",	"9R",	"10R",	
                      "11R",	"17R",	"18R",	"19R",	"20R",	
                      "21R",	"22R",	"23R",	"24R",	"25R",	
                      "26R",	"29R",	"30R", 	"32R",	"38R",	
                      "39R",	"40R",	"41R",	"42R",	"43R",	"44R",	"45R")
      }
      
      colnames(lh_sa_i)<-lh_colname
      colnames(lh_curvind_i)<-lh_colname
      colnames(lh_foldind_i)<-lh_colname
      colnames(lh_gauscurv_i)<-lh_colname
      colnames(lh_meancurv_i)<-lh_colname
      colnames(lh_AverageThickness_i)<-lh_colname
      colnames(lh_GrayMattervolume_i)<-lh_colname
      
      colnames(rh_sa_i)<-rh_colname
      colnames(rh_curvind_i)<-rh_colname
      colnames(rh_foldind_i)<-rh_colname
      colnames(rh_gauscurv_i)<-rh_colname
      colnames(rh_meancurv_i)<-rh_colname
      colnames(rh_AverageThickness_i)<-rh_colname
      colnames(rh_GrayMattervolume_i)<-rh_colname
      
      lh_si_m<-rbind(lh_sa_i, lh_curvind_i, lh_foldind_i, lh_gauscurv_i, lh_meancurv_i, lh_AverageThickness_i,lh_GrayMattervolume_i)
      
      rh_si_m<-rbind(rh_sa_i, rh_curvind_i, rh_foldind_i, rh_gauscurv_i, rh_meancurv_i, rh_AverageThickness_i,rh_GrayMattervolume_i)
      
      si_m1<-cbind(lh_si_m,rh_si_m)
      si_m<-t(apply(si_m1,1,function(x){x/max(x)}))
      si_r<-matrix(0,nrow = ncol(si_m),ncol = ncol(si_m))
      for (j in 1:ncol(si_m)) {
        for (k in j:ncol(si_m)) {
          re<-cor(si_m[,j],si_m[,k])
          if(is.na(re)!=1){
            si_r[j,k]<-re
            si_r[k,j]<-re
          }
          
        }
      }
      colnames(si_r)<-colnames(si_m1)
      write.table(si_r, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\",lh_surfarea[i,1],".txt",sep = "") , sep = '\t', col.names = T, quote = FALSE,row.names = F)
      Gi_r<-Gi_r+si_r
      print(paste(case_control[n],list0[m],lh_surfarea[i,1],i,nrow(lh_surfarea),sep="_"))
    }
    G_r<-Gi_r/nrow(lh_surfarea)
    write.table(G_r, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\Group_cor.txt",sep = ""), sep = '\t', col.names = T, quote = FALSE,row.names = F)
  }
  
  
}



############### DLBS age group MSN #############
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
  list0<-c("a2009s","DKTatlas","HOAatlas","BAatlas")#Destrieux--
  for (m in 1:4) {
    Period_13<-infor[infor$Age<40,1]
    Period_14<-infor[infor$Age>=40 & infor$Age<60,1]
    Period_15<-infor[infor$Age>=60,1]
    filelist<-list.files(paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],sep="") )
    sample1<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\",filelist[1],sep=""), sep = "\t", header = T,stringsAsFactors = F,check.names=F)
    
    Gi_r_13<-matrix(0,nrow = ncol(sample1),ncol = ncol(sample1))
    for (i in 1:length(Period_13)) {
      filename<-filelist[grep(Period_13[i],filelist)]  
      if(length(filename)!=0){
        si_r_13<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\",filename,sep=""), sep = "\t", header = T,stringsAsFactors = F,check.names=F)
        Gi_r_13<-Gi_r_13+as.matrix(si_r_13) 
      }
      
    }
    G_r_13<-Gi_r_13/length(Period_13)
    write.table(G_r_13, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\",case_control[n],"_",list0[m],"_Period_13_cor.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    
    threshold<-round(nrow(G_r_13)*nrow(G_r_13)*0.1)
    G_r_131<-cbind(colnames(G_r_13),G_r_13)
    G_r_132<-melt(as.data.frame(G_r_131), id.vars = c("V1") )
    value<-G_r_132[order(abs(as.numeric(G_r_132[,3])) ,decreasing = T)[threshold],3]
    G_r_131[G_r_131<value]<-0
    G_r_131<-G_r_131[,-1]
    write.table(G_r_131, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\",case_control[n],"_",list0[m],"_Period_13_10.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    
    
    
    Gi_r_14<-matrix(0,nrow = ncol(sample1),ncol = ncol(sample1))
    for (i in 1:length(Period_14)) {
      filename<-filelist[grep(Period_14[i],filelist)]  
      if(length(filename)!=0){
        si_r_14<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\",filename,sep=""), sep = "\t", header = T,stringsAsFactors = F,check.names=F)
        Gi_r_14<-Gi_r_14+as.matrix(si_r_14) 
      }}
    G_r_14<-Gi_r_14/length(Period_14)
    write.table(G_r_14, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\",case_control[n],"_",list0[m],"_Period_14_cor.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    
    threshold<-round(nrow(G_r_14)*nrow(G_r_14)*0.1)
    G_r_141<-cbind(colnames(G_r_14),G_r_14)
    G_r_142<-melt(as.data.frame(G_r_141), id.vars = c("V1") )
    value<-G_r_142[order(abs(as.numeric(G_r_142[,3])) ,decreasing = T)[threshold],3]
    G_r_141[G_r_141<value]<-0
    G_r_141<-G_r_141[,-1]
    write.table(G_r_141, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\",case_control[n],"_",list0[m],"_Period_14_10.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    
    
    
    Gi_r_15<-matrix(0,nrow = ncol(sample1),ncol = ncol(sample1))
    for (i in 1:length(Period_15)) {
      filename<-filelist[grep(Period_15[i],filelist)]  
      if(length(filename)!=0){
        si_r_15<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\",filename,sep=""), sep = "\t", header = T,stringsAsFactors = F,check.names=F)
        Gi_r_15<-Gi_r_15+as.matrix(si_r_15) 
      }}
    G_r_15<-Gi_r_15/length(Period_15)
    write.table(G_r_15, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\",case_control[n],"_",list0[m],"_Period_15_cor.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    
    threshold<-round(nrow(G_r_15)*nrow(G_r_15)*0.1)
    G_r_151<-cbind(colnames(G_r_15),G_r_15)
    G_r_152<-melt(as.data.frame(G_r_151), id.vars = c("V1") )
    value<-G_r_152[order(abs(as.numeric(G_r_152[,3])) ,decreasing = T)[threshold],3]
    G_r_151[G_r_151<value]<-0
    G_r_151<-G_r_151[,-1]
    write.table(G_r_151, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\",case_control[n],"_",list0[m],"_Period_15_10.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)
    print(paste(case_control[n],list0[m],sep = "_"))
    
    
    
    Group_cor<-read.table(paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\Group_cor.txt",sep=""), sep = "\t", header = T,stringsAsFactors = F,check.names=F)
    threshold<-round(nrow(Group_cor)*nrow(Group_cor)*0.1)
    Group_cor1<-cbind(colnames(Group_cor),Group_cor)
    Group_cor2<-melt(as.data.frame(Group_cor1), id.vars = c("colnames(Group_cor)") )
    value<-Group_cor2[order(abs(as.numeric(Group_cor2[,3])) ,decreasing = T)[threshold],3]
    Group_cor1[Group_cor1<value]<-0
    Group_cor1<-Group_cor1[,-1]
    write.table(Group_cor1, paste("D:\\aging\\sample_stats\\",case_control[n],"\\MSN\\",list0[m],"\\",case_control[n],"_",list0[m],"_Group_10.txt",sep=""), sep = '\t', col.names = F, quote = FALSE,row.names = F)

    
    }
 
}

