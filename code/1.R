### DN-heatmap

library(pheatmap)
aa<-read.table("DN-data.txt",sep = "\t",header = T,row.names = 1,check.names = F)
info<-read.table("DN-info.txt",header = T,sep = "\t",row.names = 1,check.names = F)
colnames(info)[1]<-"cancertype"
anno_colors =list(cancertype=c(ccRCC="#984EA3", MCC="#50B6BB",Melanoma="#A65628",
                               MM="#E7298A",TNBC="#66A61E",
                               SCC="#56B4E9",BCC="#0072B2",bladder.cancer="#D55E00",
                               B.ALL="#000000"),
                  'treat_type'=c("ICB"="#4A8952","CAR-T"="#C7E9B4"),
                  "PD-1"=c("PD-1"="#3871C1","NAA"="grey"),
                  "PD-L1"=c("PD-L1"="#CE0D8F","NAA"='grey'),
                  "CTLA4"=c("CTLA4"="#318698","NAA"="grey"))
bb<-aa
bb[is.na(bb)] <- 100
bb<-round(bb,2)
data_mark=bb
for(i in 1:24){
  for(j in 1:14){
    data_mark[i,j]=ifelse(bb[i,j] == 100, " ",bb[i,j] )
  }}

pheatmap(aa, display_numbers =data_mark,
         cluster_row = FALSE,
         cluster_col = FALSE,
         color = colorRampPalette(colors = c("#3980b8","white","#ef3b3c"))(100),
         cellwidth = 30,
         cellheight = 12,
         border_color = "white",
         breaks=c(seq(0.000001,0.99999,length.out=50),1,
                  seq(1.0000001,7,length.out=50)),
         annotation_col =info,
         annotation_colors =anno_colors
)       

min(bb)
max(bb)


###RN-heatmap

rm(list=ls())
library(pheatmap)
aa<-read.table("RN-data.txt",sep = "\t",header = T,row.names = 1,check.names = F)
info<-read.table("RN-info.txt",header = T,sep = "\t",row.names = 1,check.names = F)
colnames(info)[1]<-"cancertype"
anno_colors =list(cancertype=c(ccRCC="#984EA3", 
                              MCC="#50B6BB",
                              Melanoma="#A65628",
                              mUC="#7570B3",
                              MM="#E7298A",
                              TNBC="#66A61E"),
                    ICB=c(ICB="#1C3F21",NoICB="#C7E9B4"),
                  'treat_type'=c("ICB"="#4A8952","CAR-T"="#C7E9B4"),
                    "PD-1"=c("PD-1"="#3871C1","NAA"="grey"),
                    "PD-L1"=c("PD-L1"="#CE0D8F","NAA"='grey'),
                     "CTLA4"=c("CTLA4"="#318698","NAA"="grey"))
  bb<-aa
  bb[is.na(bb)] <- 0
  bb<-round(bb,2)
  data_mark=bb
for(i in 1:21){
  for(j in 1:8){
    data_mark[i,j]=ifelse(bb[i,j] == 0, " ",bb[i,j] )
  }}

pheatmap(aa, display_numbers =data_mark,
         cluster_row = FALSE,
         cluster_col = FALSE,
         color = colorRampPalette(colors = c("#3980b8","white","#ef3b3c"))(100),
         cellwidth = 30,
         cellheight = 12,
         border_color = "white",
         breaks=c(seq(0.000001,0.99999,length.out=50),1,
                  seq(1.0000001,7,length.out=50)),
         annotation_col =info,
         annotation_colors =anno_colors
         )       



###RN and DN Proportion

library(Seurat)
library(stringr)
library(do)
library(tidyverse)
pbmc<-readRDS("GSE145281.rds")
data<-data.frame(barCode = colnames(pbmc),cellType=pbmc@active.ident)
label<-read.table("GSE145281barLabel.txt",sep = "\t",
               header = T)
label<-label[,-5]
result<-merge(label,data,by="barCode")
unique(result$cellType)
write.table(result,"GSE145281barLabel_RN.txt",sep="\t",
            quote = F,row.names = F)

label<-read.table("GSE117988barLabel_Tumor_DN_RN.txt",sep = "\t",
                  header = T)
#celltype_2<-data.frame(celltype_2=Replace(label_1$cellType_1,
#                          from = "CD8[+] T cells",to="T cells"))
#celltype_2<-data.frame(celltype_2=Replace(celltype_2$celltype_2,
#                                          from = "CD8[+] T cells",to="T cells"))
#label<-cbind(label_1,celltype_2)
value<-unique(label$cellType_1)
#value

m<-c()
for (h in 1:length(value)) {
  cluster_1<-label[which(label$cellType_1==value[h]),]
  res<-cluster_1[which(cluster_1$DN=="Drug"),]
  p<-c()
  value_1<-unique(res$Sample)
  for (i in 1:length(value_1)) {
    result_2<-res[which(res$Sample==value_1[i]),]
    pi<-dim(result_2)[1]/dim(cluster_1)[1]
    p<-c(p,pi)
  }
  nores<-cluster_1[which(cluster_1$DN=="noDrug"),]
  p_1<-c()
  value_2<-unique(nores$Sample)
  for (j in 1:length(value_2)) {
    result_2<-nores[which(nores$Sample==value_2[j]),]
    pj<-dim(result_2)[1]/dim(cluster_1)[1]
    p_1<-c(p_1,pj)
  }
  m1<-mean(p)/mean(p_1)
  m1
  m<-c(m,m1)
}
value<-data.frame(value)
aaa<-cbind(value,m)
write.table(aaa,"GSE143317_DN_count.txt",row.names = F,sep = "\t",quote = F)

#FeaturePlot(pbmc, features = c("CHGA","CD68","MS4A1","CD3D","CD34"))
#VlnPlot(pbmc, features = c("CHGA","CD68","MS4A1","CD3D","CD34"))

#value[h]
#value_1
#p
#value_2
#p_1





###boxplot

rm(list=ls())
data<- read.table("RN.txt",header = T,sep="\t")
unique(data$cluster)
unique(data$GSE) 
#data_1<-data[which(data$GSE=="GSE123813_bcc"|data$GSE=="GSE153697"|
                     #data$GSE=="GSE164551"|data$GSE=="GSE143317"),]
data_1<-data[which(data$GSE=="ccRCC"),]
data_1<-data_1[which(data_1$cluster=="CD8+ T cells"),]

length(which(data_1$group=="R"))
length(which(data_1$group=="NR"))


#data_1<-data[which(data$cluster=="NK cells"),]
mycolor<-c("#0F6DB4","#E5B919")
library(ggplot2)
library(tidyverse)
library(ggpubr)
my_comparisons <- list(c("D","ND"))
ggboxplot(data_1, x="group", y="proportion", fill = "group",palette = "jco")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test",label = T)+
  scale_fill_manual(values = mycolor)# Add global p-value
label = "p.signif"

###survival

rm(list=ls())
library(ggpubr)
library(survminer)
library(survival)
a<-read.table("TCGA-READ-exp.txt",header = T,sep="\t",check.names = F,row.names = 1)
#exp_2<-read.table("TCGA-KIRP-exp.txt",header = T,sep="\t",check.names = F)
#exp_3<-read.table("TCGA-KICH-exp.txt",header = T,sep="\t",check.names = F)
#exp_4<-merge(exp_1,exp_3,by="Group.1")
#exp_2$Group.1<-rownames(exp_2)
#exp<-merge(exp_2,exp_4,by="Group.1")
#rownames(exp)<-exp$Group.1
#exp<-exp[,-1]
a<-as.matrix(a)
library(GSVA)
Mul_gene<-read.table("marker_fibroblast cell.txt",header = F,sep = "\t")
value_1<-list(Mul_gene$V1)
gsva.res<-gsva(a,value_1,method="ssgsea")
gsva.res<-t(gsva.res)
gsva.res<-cbind(rownames(gsva.res),gsva.res)
colnames(gsva.res)<-c("sample","score")
#su1<-read.table("TCGA-KIRC.survival.tsv",header=T,sep="\t")
#su2<-read.table("TCGA-KIRP.survival.tsv",header=T,sep="\t")
#su3<-read.table("TCGA-KICH.survival.tsv",header=T,sep="\t")
#su4<-rbind(su1,su2,su3)
#b<-su4
b<-read.table("TCGA-READ.survival.tsv",header=T,sep="\t")
b<-b[,-3]
colnames(b)<-c("sample","status","time")
b$time<-b$time/30
result<-merge(gsva.res,b,by="sample")
#write.table(result,"gsva.txt",row.names = F,sep = "\t",quote = F)
med<-as.numeric(result[,2])
med<-median(med)
result_med<-matrix(nrow=nrow(result),ncol=5)
result_med[,1]<-result[,1]
result_med[,2]<-result[,2]
result_med[,3]<-result[,3]
result_med[,4]<-result[,4]
result[,2]<-as.numeric(result[,2])
result_med[which(result[,2] >= med),5]<-"high"
result_med[which(result[,2]< med),5]<-"low"
colnames(result_med)<-c("Samples","Score","status","time","class")
library(ggpubr)
library(survminer)
library(survival)
result_med<-data.frame(result_med)
result_med$status<-as.numeric(result_med$status)
class(result_med$status)
result_med$time<-as.numeric(result_med$time)
class(result_med$time)
fit<-survfit(Surv(time,status)~class,data=result_med)
ggsurvplot(fit, data =result_med, risk.table = TRUE,palette = c("#f96d15","#50b6bb"),xlab="Time(Months)",
           pval =T)
library(maxstat)
library(survival)
result_med$Score<-as.numeric(result_med$Score)
#survival <- read.table("LIHC_surv.txt",sep = "\t",header = T)
mtHL <- maxstat.test(Surv(time,status) ~ Score, 
                     data=result_med,smethod="LogRank",pmethod="HL")
#plot(mtHL)
mtHL

        
result_med$group <- ifelse(result_med$Score>=1.318963,"high","low")
library(ggpubr)
library(survminer)
library(survival)
fit<-survfit(Surv(time,status)~group,data=result_med)###OS:    ?  ;status:    ??  class        ?
ggsurvplot(fit, data =result_med, risk.table = TRUE,palette = c("#f96d15","#50b6bb"),xlab="Time(Days)",
           pval =T)
result<-result[order(result$score),]
O<-dim(result)[1]*0.25
med_1<-data.frame(result[1:O,])
med_1$class<-"low"
result<-result[order(result$score,decreasing = T),]
med_2<-data.frame(result[1:O,])
med_2$class<-"high"
result_med<-rbind(med_1,med_2)
colnames(result_med)<-c("Samples","Score","status","time","class")
result_med<-data.frame(result_med)
result_med$status<-as.numeric(result_med$status)
class(result_med$status)
result_med$time<-as.numeric(result_med$time)
class(result_med$time)
fit<-survfit(Surv(time,status)~class,data=result_med)
ggsurvplot(fit, data =result_med, risk.table = TRUE,palette = c("#f96d15","#50b6bb"),xlab="Time(Months)",
           pval =T)

