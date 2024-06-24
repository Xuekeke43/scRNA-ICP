###monocle
rm(list = ls())
library(monocle)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggsci)
#data
pbmc<-readRDS("GSE115978_tumor_seurat.rds")
expr_matrix <- as.matrix(pbmc@assays$RNA@counts)
cell_meta <- data.frame(cluster=pbmc@meta.data[["seurat_clusters"]],
                        DN=pbmc@meta.data[["DN"]],
                        row.names =colnames(pbmc),
                        patient=pbmc@meta.data[["patient"]])
cell_meta<-cbind(cell_meta,D=cell_meta$cluster,ND=cell_meta$cluster)
cell_meta[which(cell_meta$DN=="noDrug"),4]<-NA
cell_meta[which(cell_meta$DN=="Drug"),5]<-NA
gene_meta <- data.frame(gene_short_name = row.names(pbmc@assays[["RNA"]]),
                        row.names = row.names(pbmc@assays[["RNA"]]))
dim(gene_meta)

pd <- new('AnnotatedDataFrame', data = cell_meta) 
fd <- new('AnnotatedDataFrame', data = gene_meta)

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData =fd,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1) 
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),
              num_cells_expressed >= dim(gene_meta)*0.05)) 

#expressed_genes <- VariableFeatures(pbmc) 
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~cluster") 
deg <- subset(diff, qval < 0.01)
deg <- deg[order(deg$qval,decreasing=F),]

write.table(deg,"DEG_GSE115978.txt",quote=F,sep="\t")

#ordergene <- row.names(deg)
ordergene <- row.names(deg)[order(deg$qval)][1:2000]
cds <- setOrderingFilter(cds, ordergene)  
plot_ordering_genes(cds)

# DDRTree
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "DN",
                     cell_name_size=cell)  
saveRDS(cds,"GSE115978_tumor_monocle2.rds")

###gsva and heatmap
rm(list=ls())
library(tidyverse)
library(Seurat)
library(dplyr)
library(ggplot2)
library(GSVA)
library(reshape)

rm(list = ls())
scRNA<-readRDS("GSE153697_tumor_seurat.rds")
DefaultAssay(scRNA)<-"RNA"
scRNA<-NormalizeData(scRNA)

#read txt
reviewpath <- "E:/6.gsva/sig"
completepath <- list.files(reviewpath, pattern = "*.txt$", full.names = TRUE)

read.txt <- function(x) {
  des <- readLines(x)  
    return(paste(des, sep = "\t"))     }
review <- lapply(completepath, read.txt)

docname <- list.files(reviewpath, pattern = "*.txt$")
docname<-data.frame(docname)
docname<-separate(docname,docname,c("name"),sep = "[.]")
names(review)<-as.vector(docname$name)

library(tidyverse)

Idents(scRNA)<-"seurat_clusters"
expr<-AverageExpression(scRNA,assays = "RNA",slot = "data")[[1]]
dim(expr)
expr<-expr[rowSums(expr)>0,] 
dim(expr)
expr<-as.matrix(expr)
gsva.res<-gsva(expr,review,method="ssgsea")
colnames(gsva.res)<-paste("c",colnames(gsva.res),sep = "")
gsva.res<-data.frame(gsva.res)
gsva.res<-data.frame(cbind(sig=rownames(gsva.res),gsva.res))
data_rowname <- c("Angiogenesis","DNA damage","EMT","Hypoxia","Invasion",
                  "Metastasis","Proliferation","Quiescence","Stemness",
                  "Apoptosis","Differentiation","DNA repair", "Inflammation",
                  "Cell Cycle")
gsva.res$sig<-factor(gsva.res$sig, levels=data_rowname, ordered=T)
gsva.res<-gsva.res[order(gsva.res$sig),]
gsva.res<-gsva.res[,-1]


##heatmap
pheatmap::pheatmap(gsva.res,show_colnames = T,
                   cluster_rows=F,cluster_cols = F,scale = "row",
                   main="GSE153697",angle_col = 0)



??pheatmap

##Violin Plot
vars<-colnames(gsva.res)
data<-melt(gsva.res,
           measure.vars =vars,
           variable.name = "cell",value.name = "x")
colnames(data)<-c("signature","cell","value")

cluster<-data.frame(cell=colnames(scRNA),cluster=scRNA@meta.data[["seurat_clusters"]],
                    RN=scRNA@meta.data[["RN"]])
data_1<-merge(data,cluster,by="cell")

data_2<-data_1[which(data_1$signature=="Proliferation"),]

ggplot(data_2,aes(x=cluster,y=value))+
  geom_violin(aes(fill = cluster))

mydata=FetchData(scRNA,vars = c("tSNE_1","tSNE_2","seurat_clusters"))
mydata<-cbind(cell=rownames(mydata),mydata)
gsva_data_1<-t(gsva.res)
gsva_data<-data.frame(apply(gsva_data_1,2,as.numeric))
rownames(gsva_data)<-rownames(gsva_data_1)
gsva_data<-cbind(cell=rownames(gsva_data),gsva_data)
mydata_1<-merge(mydata,gsva_data,by="cell")

class(mydata_1$Proliferation)
#TSNE
ggplot(data = mydata_1,aes(x=tSNE_1,y=tSNE_2,
                         colour=Proliferation))+geom_point()+
scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
  
###ccRCC
rm(list=ls())
library(Seurat)
library(dplyr)
library(magrittr)
pbmc.data <- Read10X("E:/ccRCC")
label<-read.table("ccRCCbarLabel.txt",header = T,sep = "\t")
#value<-label$barCode[which(label$cellType_1=="Tumor cells")]
#exp<-pbmc.data[,value]
pbmc<-CreateSeuratObject(counts = pbmc.data,min.cells =dim(pbmc.data)[1]*0.01,
                         min.features = dim(pbmc.data)[1]*0.01)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dim(pbmc)

library(Seurat)
library(tidyverse)
exp<-data.frame(pbmc@assays[["RNA"]]@counts)
head(exp)[1:5,1:5]
dim(exp)
label<-read.table("ccRCCbarLabel.txt",sep = "\t",header = T)
value<-unique(label$Sample)
test<-list()
for (i in 1:length(value)) {
  test[[i]]<-exp[,which(colnames(exp)%in%label[which(label$Sample==value[i]),1])]
}
test.seu<-list()
for (i in 1:length(value)) {
  test1.seu=CreateSeuratObject(counts = test[[i]])
  test1.seu <- NormalizeData(test1.seu, normalization.method = "LogNormalize", scale.factor = 10000)
  test.seu[[i]] <- FindVariableFeatures(test1.seu, selection.method = "vst", nfeatures = 2000)
}

### Integration ----
testAB.anchors <- FindIntegrationAnchors(
  object.list =test.seu , dims = 1:20)
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:20)
testAB.integrated
dim(testAB.integrated[["RNA"]]@counts)
dim(testAB.integrated[["RNA"]]@data)
dim(testAB.integrated[["integrated"]]@counts) 
dim(testAB.integrated[["integrated"]]@data)

DefaultAssay(testAB.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
png("Elbow.png")
ElbowPlot(testAB.integrated)
dev.off()
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:15)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.5)
#testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:15)
testAB.integrated <- RunTSNE(testAB.integrated, dims = 1:15)

label_1<-data.frame(row.names = label$barCode,patient=label$Sample,
                    cluster=label$cellType_1
)
testAB.integrated<-AddMetaData(testAB.integrated,metadata = label_1)
DimPlot(testAB.integrated, reduction = "tsne",label = F,group.by ="patient")
DimPlot(testAB.integrated, reduction = "tsne",label = T,group.by ="cluster")
DimPlot(testAB.integrated, reduction = "tsne",label = T)
saveRDS(testAB.integrated,"ccRCC-all.cell.rds")

library(Seurat)
testAB.integrated<-readRDS("ccRCC-all.cell.rds")
label<-read.table("ccRCCbarLabel.txt",header = T,sep = "\t")
value<-label$barCode[which(label$cellType_1=="Tumor cells")]
pbmc<-testAB.integrated[,value]
dim(pbmc)
exp<-as.matrix(pbmc@assays[["integrated"]]@data)


pbmc<-CreateSeuratObject(counts = exp,min.cells =dim(exp)[1]*0.01,
                         min.features = dim(exp)[1]*0.01)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dim(pbmc)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
png("Elbow.png")
ElbowPlot(pbmc)
dev.off()
pbmc <- FindNeighbors(pbmc, dims = 1:7)
pbmc <- FindClusters(pbmc, resolution = 0.2)
pbmc <- RunUMAP(pbmc, dims = 1:7)
pbmc <- RunTSNE(pbmc, dims = 1:7,check_duplicates = FALSE)
png("tsne.png")
DimPlot(pbmc, reduction = "tsne",label = TRUE)
dev.off()

label_1<-data.frame(row.names = label$barCode,
                    patient=label$Sample,
                    cluster=label$cellType_1,
                    DN=label$DN,RN=label$RN)
pbmc<-AddMetaData(pbmc,metadata = label_1)
DimPlot(pbmc, reduction = "tsne",label = F,group.by ="patient")
DimPlot(pbmc, reduction = "tsne",label = F,group.by ="DN")
DimPlot(pbmc, reduction = "tsne",label = F,group.by ="RN")
DimPlot(pbmc, reduction = "tsne",label = F,group.by ="cluster")

saveRDS(pbmc,"ccRCC-all.cell.rds")

pbmc<-readRDS("ccRCC-all.cell.rds")
library(Seurat)
label<-read.table("ccRCCbarLabel.txt",sep="\t",row.names = 1)


###bubble chart

library(ggplot2)
ggplot(data=Drug, mapping=aes(x=indexa,y=indexb,color=count))+
  geom_point(stat= "identity",aes(size=number),show.legend = F)+
  scale_color_gradient(low = "white", high = "#E71B58")+
  scale_size_continuous(range = c(6,9))+
  theme_grey()+
  theme(panel.grid.major = element_line(colour = "grey"),
        plot.background = element_rect(colour = "grey50"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


Drug<-read.table("Drug.txt",header = T,sep = "\t",row.names = 1)
noDrug<-read.table("noDrug.txt",header = T,sep = "\t",row.names = 1)

#unique(rownames(Drug))
#unique(rownames(noDrug))


noDrug[which(noDrug$number==0),2]<-0.1
noDrug[which(noDrug$count==0),1]<-0.01
#D_ND<-data.frame(count=Drug$count/noDrug$count,number=Drug$number/noDrug$number,
                #indexa=Drug$indexa,indexb=noDrug$indexb)

value<-union(rownames(Drug),rownames(noDrug))
res_1<-c()
for (i in 1:length(value)) {
  if (value[i]%in%rownames(noDrug)){
    index<-which(rownames(noDrug)==value[i])
    nd<-noDrug[index,]
    index_1<-which(rownames(Drug)==value[i])
    d<-Drug[index_1,]
    res<-data.frame(count=d$count/nd$count,number=d$number/nd$number,
                    indexa=d$indexa,indexb=d$indexb)
    rownames(res)<-value[i]
  }else{
    res<-Drug[which(rownames(Drug)==value[i]),]
    #res[,1]<-0
    #res[,2]<-0
  }  
  res_1<-rbind(res_1,res)
  }

D_ND<-res_1

#D_ND[which(rownames(D_ND)=="tumor_c3|tumor_c3"),1]<-2
D_ND[D_ND$count > 2,1]<-2
D_ND[D_ND$number > 3,2]<-3

ggplot(data=D_ND, mapping=aes(x=indexa,y=indexb,color=count))+
  geom_point(stat= "identity",aes(size=number),alpha=1.5,show.legend = T)+
  scale_colour_gradient2(low="blue",mid="white",high="red",midpoint=1,
                         limits = c(0,2), 
                         breaks = c(0,1,2))+
  scale_size_continuous(range = c(5,12))+
  theme_grey()+
  theme(panel.grid.major = element_line(colour = "grey"),
        plot.background = element_rect(colour = "grey50"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
theme(legend.position = "top")
 


