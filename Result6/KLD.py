import numpy as np
from scipy import stats
import pandas as pd
import os

def is_number(s):
    try:  # 如果能运行float(s)语句，返回True（字符串s是浮点数）
        float(s)
        return True
    except ValueError:  # ValueError为Python的一种标准异常，表示"传入无效的参数"
        pass  # 如果引发了ValueError这种异常，不做任何事情（pass：不做任何事情，一般用做占位语句）
    try:
        import unicodedata  # 处理ASCii码的包
        unicodedata.numeric(s)  # 把一个表示数字的字符串转换为浮点数返回的函数
        return True
    except (TypeError, ValueError):
        pass
    return False




case_control=["GBM","DLBS"]
for n in case_control:
    list0=["DKTatlas","HOAatlas","BAatlas"]#Destrieux--
    for m in list0:
        lh_surfarea = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\lh_surfarea_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        lh_curvind = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\lh_curvind_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        lh_foldind = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\lh_foldind_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        lh_gauscurv = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\lh_gauscurv_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        lh_meancurv = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\lh_meancurv_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        lh_AverageThickness = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\lh_AverageThickness_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        lh_GrayMattervolume = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\lh_GrayMattervolume_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        rh_surfarea = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\rh_surfarea_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        rh_curvind = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\rh_curvind_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        rh_foldind = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\rh_foldind_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        rh_gauscurv = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\rh_gauscurv_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        rh_meancurv = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\rh_meancurv_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        rh_AverageThickness = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\rh_AverageThickness_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        rh_GrayMattervolume = pd.read_table("D:\\aging\\sample_stats\\"+ n+"\\rh_GrayMattervolume_"+m+".txt", header=0, index_col=0, sep="\t",low_memory=False)
        n1=(lh_surfarea.shape[1]-3)+(rh_surfarea.shape[1]-3)

        Gi_r=np.zeros((n1,n1))
        n2=lh_surfarea.shape[0] # sample numbers 313
        for i in range(0,n2):
            lh_sa_i = lh_surfarea.iloc[i, 0:(lh_surfarea.shape[1]-3)]
            lh_curvind_i = lh_curvind.iloc[i, 0:(lh_surfarea.shape[1]-3)]
            lh_foldind_i = lh_foldind.iloc[i, 0:(lh_surfarea.shape[1]-3)]
            lh_gauscurv_i = lh_gauscurv.iloc[i, 0:(lh_surfarea.shape[1]-3)]
            lh_meancurv_i = lh_meancurv.iloc[i, 0:(lh_surfarea.shape[1]-3)]
            lh_AverageThickness_i = lh_AverageThickness.iloc[i, 0:(lh_surfarea.shape[1]-3)]
            lh_GrayMattervolume_i = lh_GrayMattervolume.iloc[i, 0:(lh_surfarea.shape[1]-3)]

            rh_sa_i = rh_surfarea.iloc[i, 0:(rh_surfarea.shape[1]-3)]
            rh_curvind_i = rh_curvind.iloc[i, 0:(rh_surfarea.shape[1]-3)]
            rh_foldind_i = rh_foldind.iloc[i, 0:(rh_surfarea.shape[1]-3)]
            rh_gauscurv_i = rh_gauscurv.iloc[i, 0:(rh_surfarea.shape[1]-3)]
            rh_meancurv_i = rh_meancurv.iloc[i, 0:(rh_surfarea.shape[1]-3)]
            rh_AverageThickness_i = rh_AverageThickness.iloc[i, 0:(rh_surfarea.shape[1]-3)]
            rh_GrayMattervolume_i = rh_GrayMattervolume.iloc[i, 0:(rh_surfarea.shape[1]-3)]

            lh_si_m = pd.DataFrame(list(zip(lh_sa_i, lh_curvind_i,lh_foldind_i,lh_gauscurv_i, lh_meancurv_i, lh_AverageThickness_i, lh_GrayMattervolume_i)))
            rh_si_m = pd.DataFrame(list(zip(rh_sa_i, rh_curvind_i,rh_foldind_i,rh_gauscurv_i, rh_meancurv_i, rh_AverageThickness_i, rh_GrayMattervolume_i)))
            si_m1 = lh_si_m.append(rh_si_m)
            si_m = si_m1 / si_m1.max(axis=0)
            si_r = np.zeros((si_m.shape[0], si_m.shape[0]))
            for j in range(0,si_m.shape[0]):
                for k in range(0,si_m.shape[0]):
                    px = si_m.iloc[j,:]
                    py = si_m.iloc[k, :]
                    KL = stats.entropy(px, py)
                    KL1 = stats.entropy(py, px)
                    D=1/(1+KL+KL1)
                    si_r[j,k]=D
            #si_r1 = pd.DataFrame(list(zip( colname, si_r)))
            si_r = np.nan_to_num(si_r, nan=0)
            np.savetxt("D:\\aging\\sample_stats\\"+n+"\\MIND\\"+m+"\\"+lh_surfarea.index[i]+".txt", si_r, delimiter='\t', fmt = '%s')
            Gi_r = Gi_r + si_r
        G_r = Gi_r / n2
        np.savetxt("D:\\aging\\sample_stats\\" + n + "\\MIND\\" + m + "\\Group_cor.txt", G_r, delimiter='\t', fmt = '%s')

