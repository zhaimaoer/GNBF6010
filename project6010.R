setwd('D:/project')

#check data
GSEs=list.files(path = "./GSE116256_raw/",recursive = T,full.names = T)
head(GSEs)
library(data.table)
fread(GSEs[1])[1:4,1:4]
fread(GSEs[2])[1,]

#subset healthy BM group
library(stringr)
BM=list.files("./GSE116256_raw",pattern="*_BM*",full.names = T)
BM_dem=BM[str_detect(BM,"dem")]
BM_anno=BM[str_detect(BM,"anno")]
samples=str_split(str_split(BM_dem,"_",simplify = T)[,3],"[.]",simplify = T)[,1]
res=do.call(rbind,lapply(1:4,function(i){
  dem=as.data.frame(fread(BM_dem[i]))
  rownames(dem)=dem[,1]
  dem=dem[,-1]
  anno=as.data.frame(fread(BM_anno[i]))
  return(c(dim(dem),dim(anno)))
}))
colnames(res)=c("genes","cells","cells","annotations")
rownames(res)=samples
res

#Intergrate the expression matrix
all_BM=lapply(1:length(BM_dem),function(i){
  BM <- as.data.frame(fread(BM_dem[i]))
  rownames(BM)=BM[,1]
  BM=BM[,-1]
})
all_BM.df=do.call(cbind,all_BM)
table(do.call(rbind,strsplit(colnames(all_BM.df),"_"))[,1])

#Intergrate the cell type annotation(paper)
all_BM_anno=lapply(1:length(BM_anno),function(i){
  anno <- as.data.frame(fread(BM_anno[i]))
})
all_BM_anno.df=do.call(rbind,all_BM_anno)
BM_cell_ref=all_BM_anno.df[,c("Cell","CellType")]
BM_cell_ref$sub_group=do.call(rbind,strsplit(BM_cell_ref$Cell,"_"))[,1]
BM_cell_ref=BM_cell_ref[,c(1,3,2)]
BM_cell_ref=subset(BM_cell_ref,CellType!="")
BM_cell_ref$CT4=ifelse(BM_cell_ref[,3] %in% c("HSC","Prog"),"Undifferentiated",
                       ifelse(BM_cell_ref[,3] %in% c("GMP","ProMono","Mono","cDC","pDC"), "Myeloid",
                              ifelse(BM_cell_ref[,3] %in% c("earlyEry","lateEry"), "Erythroid","Lymphoid")))
head(BM_cell_ref);dim(BM_cell_ref)
table(BM_cell_ref[,3],BM_cell_ref[,4])
save(all_BM.df,BM_cell_ref,file="./BM_input.rda")
rm(list=ls())

#Subset the aml group
library(stringr)
AML=list.files("./GSE116256_raw",pattern="*-D0",full.names = T)
AML_dem=AML[str_detect(AML,"dem")]
AML_anno=AML[str_detect(AML,"anno")]
samples=str_split(str_split(AML_dem,"_",simplify = T)[,3],"[.]",simplify = T)[,1]
res=do.call(rbind,lapply(1:16,function(i){
  dem=as.data.frame(fread(AML_dem[i]))
  rownames(dem)=dem[,1]
  dem=dem[,-1]
  anno=as.data.frame(fread(AML_anno[i]))
  return(c(dim(dem),dim(anno)))
}))
colnames(res)=c("genes","cells","cells","annotations")
rownames(res)=samples
res

#Intergrate the expression matrix
all_AML=lapply(1:length(AML_dem),function(i){
  test <- as.data.frame(fread(AML_dem[i]))
  rownames(test)=test[,1]
  test=test[,-1]
})
all_AML.df=do.call(cbind,all_AML)
table(do.call(rbind,strsplit(colnames(all_AML.df),"_"))[,1])

#Intergrate the cell type annotation(paper)
all_AML_anno=lapply(1:16,function(i){
  anno=as.data.frame(fread(AML_anno[i]))
  anno=anno[,1:28]
})
all_AML_anno.df=do.call(rbind,all_AML_anno)
AML_cell_ref=all_AML_anno.df[,c("Cell","CellType")]
AML_cell_ref$sub_group=do.call(rbind,strsplit(AML_cell_ref$Cell,"_"))[,1]
AML_cell_ref=AML_cell_ref[,c(1,3,2)]
AML_cell_ref=subset(AML_cell_ref,CellType!="")
AML_cell_ref$CT4=ifelse(AML_cell_ref[,3] %in% c("HSC","HSC-like","Prog","Prog-like"),"Undifferentiated",
                        ifelse(AML_cell_ref[,3] %in% c("GMP","GMP-like","ProMono","ProMono-like","Mono","Mono-like","cDC","cDC-like","pDC"), "Myeloid",
                               ifelse(AML_cell_ref[,3] %in% c("earlyEry","lateEry"), "Erythroid","Lymphoid")))
table(AML_cell_ref[,3],AML_cell_ref[,4])

save(all_AML.df,AML_cell_ref,file="./AML_input.rda")
rm(list=ls())

#create Seurat object and check the quality
library(Seurat)
load("./BM_input.rda")
scRNA_BM = CreateSeuratObject(counts=all_BM.df)
#filter cell cycle related genes
cc.gene=c('ASPM', 'CENPE', 'CENPF', 'DLGAP5', 'MKI67',
          'NUSAP1', 'PCLAF', 'STMN1', 'TOP2A', 'TUBB')
scRNA_BM=scRNA_BM[!(rownames(scRNA_BM) %in% cc.gene), ]
scRNA_BM
library(ggplot2)
ggplot(scRNA_BM@meta.data,aes(x=orig.ident)) +
  geom_bar()
VlnPlot(scRNA_BM,
        features = c("nFeature_RNA","nCount_RNA"),
        group.by = "orig.ident") 
FeatureScatter(scRNA_BM,
               feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
               group.by = "orig.ident",pt.size = 1.3)
#Reduce dims and cluster
scRNA_BM <- NormalizeData(scRNA_BM, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA_BM <- ScaleData(scRNA_BM, features = (rownames(scRNA_BM)))
scRNA_BM <- FindVariableFeatures(scRNA_BM, selection.method = "vst", nfeatures = 1500) 
scRNA_BM <- RunPCA(scRNA_BM, features = VariableFeatures(scRNA_BM)) 
DimPlot(scRNA_BM, reduction = "pca", group.by="orig.ident")
DimPlot(scRNA_BM, reduction = "pca", split.by="orig.ident")
ElbowPlot(scRNA_BM, ndims=50, reduction="pca")
pc.num=1:25
scRNA_BM = RunTSNE(scRNA_BM, dims = pc.num)
DimPlot(scRNA_BM, reduction = "tsne", group.by="orig.ident")
DimPlot(scRNA_BM, reduction = "tsne", split.by="orig.ident")
scRNA_BM <- FindNeighbors(scRNA_BM, dims = pc.num) 
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)) {
  scRNA_BM=FindClusters(scRNA_BM, resolution = res, algorithm = 1)
}
str(scRNA_BM@meta.data)


library(clustree)
clustree(scRNA_BM@meta.data, prefix = "RNA_snn_res.")
lapply(scRNA_BM@meta.data[,grep("snn",colnames(scRNA_BM@meta.data))],function(i){
  table(i,as.character(scRNA_BM@meta.data[,"orig.ident"]))
})
sel.clust = "RNA_snn_res.0.9"
table(scRNA_BM@meta.data[,sel.clust])
DimPlot(scRNA_BM, reduction = "tsne", group.by=sel.clust)
tb_BM_0.9=as.data.frame(table(scRNA_BM@meta.data[,sel.clust],
                              scRNA_BM@meta.data[,"orig.ident"]))
colnames(tb_BM_0.9)=c("Cluster","BM","Cells")
ggplot(tb_BM_0.9, aes(x=BM, y=Cells, fill=Cluster)) +
  geom_bar(stat = "identity",position="fill")

#find marker genes for each cluster
library(dplyr)
BM.markers <- FindAllMarkers(object = scRNA_BM, only.pos = TRUE, min.pct = 0.25)
top10 <- BM.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
DoHeatmap(scRNA_BM,features= top10$gene)+NoLegend()

#SingleR celltype annotation
library(SingleR)
#The Novershtern reference
#ref <- NovershternHematopoieticData()
ref=get(load("NovershternHematopoieticData.rda"))
clusters=scRNA_BM@meta.data$"RNA_snn_res.0.9"
testdata <- GetAssayData(scRNA_BM, slot="data")
res <- SingleR(test = testdata, ref = ref,
               labels = ref$label.main,method = "cluster",
               clusters = clusters, 
               assay.type.test = "logcounts", assay.type.ref = "logcounts")
pred.BM=data.frame(cluster=0:(nrow(res)-1),SingleR=res[,4])
scRNA_BM$SingleR=pred.BM$SingleR[match(scRNA_BM$RNA_snn_res.0.9,pred.BM$cluster)]
scRNA_BM$SingleR_CT4=ifelse(scRNA_BM$SingleR %in% c("Monocytes","Dendritic cells"),"Myeloid",
                            ifelse(scRNA_BM$SingleR %in% c("HSCs"),"Undifferentiated",
                                   ifelse(scRNA_BM$SingleR %in% c("Erythroid cells","MEPs"),
                                          "Erythroid","Lymphoid")))
table(scRNA_BM$SingleR,scRNA_BM$SingleR_CT4)

#celltype visualization(VS Paper result)
tmp=as.data.frame(table(scRNA_BM@meta.data[,"SingleR"],
                        scRNA_BM@meta.data[,"orig.ident"]))
colnames(tmp)=c("SingleR","BM","Cells")
ggplot(tmp, aes(x=BM, y=Cells, fill=SingleR)) +
  geom_bar(stat = "identity",position="fill") +
  scale_fill_brewer(palette="Set1")

tb1=as.data.frame(table(BM_cell_ref$CT4))
scRNA_BM_anno=scRNA_BM[,colnames(scRNA_BM) %in% BM_cell_ref$Cell]
tb2=as.data.frame(table(scRNA_BM_anno$SingleR_CT4))
tb1_2=rbind(tb1,tb2)
colnames(tb1_2)=c("CT4","Number")
tb1_2$group=rep(c("Ref","SingleR"),each=4)
ggplot(tb1_2, aes(x=group, y=Number, fill=CT4)) +
  geom_bar(stat = "identity",position="fill")


library(patchwork)
scRNA_BM_anno$Ref_CT4=BM_cell_ref$CT4[match(colnames(scRNA_BM_anno),BM_cell_ref$Cell)]
scRNA_BM_anno$Ref=BM_cell_ref$CellType[match(colnames(scRNA_BM_anno),BM_cell_ref$Cell)]
str(scRNA_BM_anno@meta.data)
DimPlot(scRNA_BM_anno, reduction = "tsne", group.by="SingleR_CT4") + 
  DimPlot(scRNA_BM_anno, reduction = "tsne", group.by="Ref_CT4") +
  plot_layout(guides = "collect")

DimPlot(scRNA_BM_anno, reduction = "tsne", group.by="SingleR") + 
  DimPlot(scRNA_BM_anno, reduction = "tsne", group.by="Ref")

tb=table(scRNA_BM_anno$SingleR,scRNA_BM_anno$Ref)
colnames(tb)=paste0("Ref_",colnames(tb))
rownames(tb)=paste0("SingleR_",rownames(tb))
tmp=data.frame(scRNA_BM_anno$Ref_CT4,scRNA_BM_anno$Ref)
tmp=tmp[!(duplicated(tmp[,2])),]
annotation_col=data.frame(Ref_CT4=factor(tmp[,1]))
rownames(annotation_col)=paste0("Ref_",tmp[,2])

tmp=data.frame(scRNA_BM_anno$SingleR_CT4,scRNA_BM_anno$SingleR)
tmp=tmp[!(duplicated(tmp[,2])),]
annotation_row=data.frame(SingleR_CT4=factor(tmp[,1]))
rownames(annotation_row)=paste0("SingleR_",tmp[,2])
ann_colors = list(
  Ref_CT4 = c(Erythroid="red",Lymphoid="blue",Myeloid="green",Undifferentiated="yellow"),
  SingleR_CT4 = c(Erythroid="red",Lymphoid="blue",Myeloid="green",Undifferentiated="yellow")
)
library('pheatmap')
pheatmap::pheatmap(tb,cluster_rows = F,cluster_cols = F,legend = F,
                   display_numbers = TRUE,number_format = "%.0f",
                   angle_col = "90",annotation_colors = ann_colors,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row)
write.csv(scRNA_BM@meta.data[,c(1,12,13,14)], file="./health_BM_anno.csv")
if(!(file.exists("./scRNA_BM.rda"))) {
  save(scRNA_BM,file="./scRNA_BM.rda")
}

#AML group
library(Seurat)
load("./AML_input.rda")
scRNA_AML = CreateSeuratObject(counts=all_AML.df)
#cell cycle realted genes
cc.gene=c('ASPM', 'CENPE', 'CENPF', 'DLGAP5', 'MKI67',
          'NUSAP1', 'PCLAF', 'STMN1', 'TOP2A', 'TUBB')
scRNA_AML=scRNA_AML[!(rownames(scRNA_AML) %in% cc.gene), ]
scRNA_AML
library(ggplot2)
ggplot(scRNA_AML@meta.data,aes(x=orig.ident)) +
  geom_bar()
VlnPlot(scRNA_AML,pt.size = 0,
        features = c("nFeature_RNA","nCount_RNA"),
        group.by = "orig.ident")
FeatureScatter(scRNA_AML,
               feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
               group.by = "orig.ident",pt.size = 1.3)

#Reduce dims and cluster
scRNA_AML <- NormalizeData(scRNA_AML, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA_AML <- ScaleData(scRNA_AML, features = (rownames(scRNA_AML)))
scRNA_AML <- FindVariableFeatures(scRNA_AML, selection.method = "vst", nfeatures = 1500) 
scRNA_AML <- RunPCA(scRNA_AML, features = VariableFeatures(scRNA_AML)) 
DimPlot(scRNA_AML, reduction = "pca", group.by="orig.ident")
DimPlot(scRNA_AML, reduction = "pca", split.by="orig.ident",ncol=4)
ElbowPlot(scRNA_AML, ndims=50, reduction="pca") 
pc.num=1:30
scRNA_AML = RunTSNE(scRNA_AML, dims = pc.num)
DimPlot(scRNA_AML, reduction = "tsne", group.by="orig.ident")
DimPlot(scRNA_AML, reduction = "tsne", split.by="orig.ident",ncol=4)
scRNA_AML <- FindNeighbors(scRNA_AML, dims = pc.num) 
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)) {
  scRNA_AML=FindClusters(scRNA_AML, resolution = res, algorithm = 1)
}
str(scRNA_AML@meta.data)

apply(scRNA_AML@meta.data[,grep("RNA_snn_res",colnames(scRNA_AML@meta.data))],2,table)
clustree(scRNA_AML@meta.data, prefix = "RNA_snn_res.")

sel.clust = "RNA_snn_res.0.3"
table(scRNA_AML@meta.data[,sel.clust])

DimPlot(scRNA_AML, reduction = "tsne", group.by=sel.clust)
tb_AML_0.3=as.data.frame(table(scRNA_AML@meta.data[,sel.clust],
                               scRNA_AML@meta.data[,"orig.ident"]))
colnames(tb_AML_0.3)=c("Cluster","AML","Cells")
ggplot(tb_AML_0.3, aes(x=AML, y=Cells, fill=Cluster)) +
  geom_bar(stat = "identity",position="fill")

#SingleR celltype annotation
library(SingleR)
#library(celldex)
ref <- get(load("NovershternHematopoieticData.rda"))
clusters=scRNA_AML@meta.data$"RNA_snn_res.0.3"
testdata <- GetAssayData(scRNA_AML, slot="data")
res <- SingleR(test = testdata, ref = ref,
               labels = ref$label.main,method = "cluster",
               clusters = clusters, 
               assay.type.test = "logcounts", assay.type.ref = "logcounts")
pred.AML=data.frame(cluster=0:(nrow(res)-1),SingleR=res[,4])
scRNA_AML$SingleR=pred.AML$SingleR[match(scRNA_AML$RNA_snn_res.0.3,pred.AML$cluster)]
scRNA_AML$SingleR_CT4=ifelse(scRNA_AML$SingleR %in% c("Monocytes","Granulocytes",'Dendritic cells'),"Myeloid",
                             ifelse(scRNA_AML$SingleR %in% c("HSCs"),"Undifferentiated",
                                    ifelse(scRNA_AML$SingleR %in% c("Erythroid cells","MEPs"),
                                           "Erythroid","Lymphoid")))
table(scRNA_AML$SingleR,scRNA_AML$SingleR_CT4)

#celltype visualization(VS Paper result)
tmp=as.data.frame(table(scRNA_AML@meta.data[,"SingleR"],
                        scRNA_AML@meta.data[,"orig.ident"]))
colnames(tmp)=c("SingleR","BM","Cells")
ggplot(tmp, aes(x=BM, y=Cells, fill=SingleR)) +
  geom_bar(stat = "identity",position="fill") +
  scale_fill_brewer(palette="Set1")

tb1=as.data.frame(table(AML_cell_ref$CT4))
scRNA_AML=scRNA_AML[,colnames(scRNA_AML) %in% AML_cell_ref$Cell]
tb2=as.data.frame(table(scRNA_AML$SingleR_CT4))
tb1_2=rbind(tb1,tb2)
colnames(tb1_2)=c("CT4","Number")
tb1_2$group=rep(c("Ref","SingleR"),each=4)
ggplot(tb1_2, aes(x=group, y=Number, fill=CT4)) +
  geom_bar(stat = "identity",position="fill")

library(patchwork)
scRNA_AML$Ref_CT4=AML_cell_ref$CT4[match(colnames(scRNA_AML),AML_cell_ref$Cell)]
scRNA_AML$Ref=AML_cell_ref$CellType[match(colnames(scRNA_AML),AML_cell_ref$Cell)]
str(scRNA_AML@meta.data)

DimPlot(scRNA_AML, reduction = "tsne", group.by="SingleR_CT4") + 
  DimPlot(scRNA_AML, reduction = "tsne", group.by="Ref_CT4") +
  plot_layout(guides = "collect")
DimPlot(scRNA_AML, reduction = "tsne", group.by="SingleR") + 
  DimPlot(scRNA_AML, reduction = "tsne", group.by="Ref")

tb=table(scRNA_AML$SingleR,scRNA_AML$Ref)
colnames(tb)=paste0("Ref_",colnames(tb))
rownames(tb)=paste0("SingleR_",rownames(tb))
tmp=data.frame(scRNA_AML$Ref_CT4,scRNA_AML$Ref)
tmp=tmp[!(duplicated(tmp[,2])),]
annotation_col=data.frame(Ref_CT4=factor(tmp[,1]))
rownames(annotation_col)=paste0("Ref_",tmp[,2])

tmp=data.frame(scRNA_AML$SingleR_CT4,scRNA_AML$SingleR)
tmp=tmp[!(duplicated(tmp[,2])),]
annotation_row=data.frame(SingleR_CT4=factor(tmp[,1]))
rownames(annotation_row)=paste0("SingleR_",tmp[,2])
ann_colors = list(
  Ref_CT4 = c(Erythroid="red",Lymphoid="blue",Myeloid="green",Undifferentiated="yellow"),
  SingleR_CT4 = c(Erythroid="red",Lymphoid="blue",Myeloid="green",Undifferentiated="yellow")
)
pheatmap::pheatmap(tb,cluster_rows = F,cluster_cols = F,legend = F,
                   display_numbers = TRUE,number_format = "%.0f",
                   angle_col = "90",annotation_colors = ann_colors,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row)
write.csv(scRNA_AML@meta.data[,c(1,9,13,14)], file="./AML_BM_anno.csv")
if(!(file.exists("scRNA_AML.rda"))) {
  save(scRNA_AML,file="scRNA_AML.rda")
}
rm(list=ls())

#comparison of SingleR results vs paper results by bar plot
BM=list.files("./GSE116256_raw",pattern="*_BM*",full.names = T)
BM_dem=BM[str_detect(BM,"dem")]
BM_anno=BM[str_detect(BM,"anno")]

all_BM_anno=lapply(1:length(BM_anno),function(i){
  anno <- as.data.frame(fread(BM_anno[i]))
})
all_BM_anno.df=do.call(rbind,all_BM_anno)
BM_cell_ref=all_BM_anno.df[,c("Cell","CellType")]
BM_cell_ref$sub_group=do.call(rbind,strsplit(BM_cell_ref$Cell,"_"))[,1]
BM_cell_ref=BM_cell_ref[,c(1,3,2)]
BM_cell_ref=subset(BM_cell_ref,CellType!="")
BM_cell_ref$CT4=ifelse(BM_cell_ref[,3] %in% c("HSC","Prog"),"Undifferentiated",
                       ifelse(BM_cell_ref[,3] %in% c("GMP","ProMono","Mono","cDC","pDC"), "Myeloid",
                              ifelse(BM_cell_ref[,3] %in% c("earlyEry","lateEry"), "Erythroid","Lymphoid")))
head(BM_cell_ref)

SingleR_BM=read.csv("./health_BM_anno.csv")
head(SingleR_BM)
SingleR_BM=SingleR_BM[SingleR_BM$X %in% BM_cell_ref$Cell,]

AML=list.files("./GSE116256_raw",pattern="*-D0",full.names = T)
AML_dem=AML[str_detect(AML,"dem")]
AML_anno=AML[str_detect(AML,"anno")]

library(ggplot2)
library(patchwork)
p1=ggplot(SingleR_BM,aes(x=SingleR_CT4,fill=SingleR)) + 
  geom_bar() + ylim(c(0,3000)) + xlab("4 main cell types") +
  scale_fill_discrete(name="Cell type") +
  ggtitle("4677 heathy BM cell-type annotated by SingleR") +
  theme(plot.title = element_text(hjust = 0.5))

p2=ggplot(BM_cell_ref,aes(x=CT4,fill=CellType)) + 
  geom_bar() + ylim(c(0,3000)) + xlab("4 main cell types") +
  scale_fill_discrete(name="Cell type") +
  ggtitle("4677 heathy BM cell-type annotated by Paper") +
  theme(plot.title = element_text(hjust = 0.5))
#differences of healthy group
p_1=p1 + p2


all_AML_anno=lapply(1:16,function(i){
  anno=as.data.frame(fread(AML_anno[i]))
  anno=anno[,1:28]
})
all_AML_anno.df=do.call(rbind,all_AML_anno)
AML_cell_ref=all_AML_anno.df[,c("Cell","CellType")]
AML_cell_ref$sub_group=do.call(rbind,strsplit(AML_cell_ref$Cell,"_"))[,1]
AML_cell_ref=AML_cell_ref[,c(1,3,2)]
AML_cell_ref=subset(AML_cell_ref,CellType!="")
AML_cell_ref$CT4=ifelse(AML_cell_ref[,3] %in% c("HSC","HSC-like","Prog","Prog-like"),"Undifferentiated",
                        ifelse(AML_cell_ref[,3] %in% c("GMP","GMP-like","ProMono","ProMono-like","Mono","Mono-like","cDC","cDC-like","pDC"), "Myeloid",
                               ifelse(AML_cell_ref[,3] %in% c("earlyEry","lateEry"), "Erythroid","Lymphoid")))
head(AML_cell_ref)

SingleR_AML=read.csv("./AML_BM_anno.csv")
head(SingleR_AML)
SingleR_AML=SingleR_AML[SingleR_AML$X %in% AML_cell_ref$Cell,]


p3=ggplot(SingleR_AML,aes(x=SingleR_CT4,fill=SingleR)) + 
  geom_bar() + ylim(c(0,10000)) + xlab("4 main cell types") +
  scale_fill_discrete(name="Cell type") +
  ggtitle("15685 AML cell-type annotated by SingleR") +
  theme(plot.title = element_text(hjust = 0.5))

p4=ggplot(AML_cell_ref,aes(x=CT4,fill=CellType)) + 
  geom_bar() + ylim(c(0,10000)) + xlab("4 main cell types") +
  scale_fill_discrete(name="Cell type") +
  ggtitle("15685 AML cell-type annotated by Paper") +
  theme(plot.title = element_text(hjust = 0.5))
#differences of AML group
p_2=p3 + p4

p_1
p_2

#differences between healthy and AML using SingleR
t_SingleR_BM <- as.data.frame(table(SingleR_BM[,5]))
t_SingleR_AML <- as.data.frame(table(SingleR_AML[,5]))
SinglR_4 <- rbind(t_SingleR_BM,t_SingleR_AML)
colnames(SinglR_4)<-c("CT4","Number")
SinglR_4$group=rep(c("Healthy","AML"),each=4)
ggplot(SinglR_4, aes(x=group, y=Number, fill=CT4)) +
  geom_bar(stat = "identity",position="fill")

t_SingleR_BM_all <- as.data.frame(table(SingleR_BM[,4]))
t_SingleR_AML_all <- as.data.frame(table(SingleR_AML[,4]))
SinglR_all <- rbind(t_SingleR_BM_all,t_SingleR_AML_all)
colnames(SinglR_all)<-c("cells","Number")
SinglR_all$group=rep(c("Healthy","AML"),each=8)
ggplot(SinglR_all, aes(x=group, y=Number, fill=cells)) +
  geom_bar(stat = "identity",position="fill")
