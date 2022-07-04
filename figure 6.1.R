source("http://bioconductor.org/biocLite.R")
biocLite()
library(monocle)
install.packages("devtools")
devtools::install_github("cole-trapnell-lab/monocle-release@develop")


#Reading the table, read.table default delimiter is a space
data<-read.table(file = "GSE124904_aggregate_gene_cell_matrix.txt")


library(dplyr)
library(Seurat)
library(Matrix)

### Convert my data to sparse matrix

datasparse=as(as.matrix(data),"sparseMatrix")

### Write sparse matrix
writeMM(datasparse,"matrix.mtx")
### Create barcodes and genes files trhat are required for Seurat
write.table(colnames(datasparse), 'barcodes.tsv', quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
write.table(cbind(rownames(datasparse),rownames(datasparse)), 'genes.tsv', quote=FALSE, sep='\t',row.names=FALSE, col.names=FALSE)

Whole.data <- Read10X(data.dir = "~/Desktop/Law et al SCS/8 datasets/Seurat/")
Whole <- CreateSeuratObject(counts = Whole.data, project = "WholeGB", min.cells = 3, min.features = 200)
Whole[["percent.mt"]] <- PercentageFeatureSet(Whole, pattern = "^mt-")
VlnPlot(Whole, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(Whole, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Whole, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

E16.1 <- subset(E16, idents=c("E16.1"))
VlnPlot(E16.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) ## looks like the aggregated matrix has already been QCed
E16.1 <- subset(E16.1, subset = nFeature_RNA > 1300 )
E16.2 <- subset(E16, idents=c("E16.2"))
VlnPlot(E16.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E16.2 <- subset(E16.2, subset = nFeature_RNA > 1900 & nCount_RNA < 80000 )
E16.3 <- subset(E16, idents=c("E16.3"))
VlnPlot(E16.3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Whole <- NormalizeData(Whole, normalization.method = "LogNormalize", scale.factor = 10000)

Whole <- FindVariableFeatures(Whole, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Whole), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Whole)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(Whole)
Whole <- ScaleData(Whole, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt","orig.ident"))

Whole <- RunPCA(Whole, features = VariableFeatures(object = Whole))
VizDimLoadings(Whole, dims = 1:2, reduction = "pca")
DimPlot(Whole, reduction = "pca")
DimHeatmap(Whole, dims = 1, cells = 500, balanced = TRUE)

Whole <- JackStraw(Whole, num.replicate = 100)
Whole <- ScoreJackStraw(Whole, dims = 1:40)
JackStrawPlot(Whole, dims = 1:30)
ElbowPlot(Whole)
DimHeatmap(Whole, dims = 1:32, cells = 500, balanced = TRUE)


Whole <- FindNeighbors(Whole, dims = 1:31)
Whole <- FindClusters(Whole, resolution = 0.6)
head(Idents(Whole), 5)

Whole <- RunUMAP(Whole, dims = 1:31)
DimPlot(Whole, reduction = "umap")
saveRDS(Whole, file = "~/Desktop/Law et al SCS/8 datasets/Whole.rds")

setwd("~/Desktop/Law et al SCS/")
Whole<-readRDS("8 datasets/Whole.rds")

DimPlot(Whole, reduction = "umap", label=TRUE)
FeaturePlot(Whole, features = c("Tdgf1"), pt.size=2, order=TRUE) 
FeaturePlot(Whole, features = c("Nanog"), pt.size=2, order=TRUE)
FeaturePlot(Whole, features = c("Pou5f1"), pt.size=2, order=TRUE)
FeaturePlot(Whole, features = c("Sox2"), pt.size=2, order=TRUE)
FeaturePlot(Whole, features = c("Lhx1","Id4","Gfra1","Etv5","Lhx1"), pt.size=1.5, order=TRUE)
VlnPlot(Whole, features = c("Tdgf1"))
VlnPlot(Whole, features = c("Tdgf1"))
DotPlot(Whole, features = c("Etv5","Gfra1","Id4","Lhx1","Sox2","Pou5f1","Nanog","Tdgf1"), dot.scale=25) + scale_color_viridis() +   coord_flip() 
DotPlot(Whole, features = c("Tdgf1"), dot.scale=30) + scale_color_viridis() +   coord_flip()

a<-DotPlot(object=Whole, features=c("Tdgf1"), group.by="orig.ident")
a$data
b<-subset(x = Whole, subset = Tdgf1 > 0)
FeaturePlot(b, features = c("Tdgf1"), pt.size=2, order=TRUE) 
c<-DotPlot(object=b, features=c("Tdgf1"), group.by="orig.ident")
c$data

cluster.averages <- AverageExpression(Whole)
head(cluster.averages[["RNA"]][, 1:5])

orig.levels <- levels(Whole)
Idents(Whole) <- gsub(pattern = " ", replacement = "_", x = Idents(Whole))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(Whole) <- orig.levels
cluster.averages <- AverageExpression(Whole, return.seurat = TRUE)
cluster.averages

my_levels <- c(11,4,2,6,1,5,8,7,12,3,0,9,10)
levels(Whole) <- my_levels

DoHeatmap(cluster.averagesE16, features = c("Nanog","Sox2","Tdgf1","Etv5","Id4","Lhx1","Ret","Gfra1","Nanos2","Dnmt3a","Sohlh1"), size = 3, 
          draw.lines = FALSE)

DoHeatmap(cluster.averages, features = c("Nanog","Sox2","Tdgf1","Nodal","Lefty1","Lefty2","Etv5","Id4","Lhx1","Gfra1","Ret","Dnmt3a","Sohlh1","Kit","Stra8"),  draw.lines = FALSE) + scale_fill_gradientn(colors = c("blue", "white", "red"))
DoHeatmap(cluster.averages, feature = c("Zbtb16"))
E16 <- subset(Whole, idents = c("2", "4", "6","11"))
DimPlot(E16, reduction = "umap")
DotPlot(E16, features = c("Nanog","Sox2","Tdgf1","Etv5","Id4","Lhx1","Ret","Dnmt3a","Sohlh1")) + RotatedAxis()
P0 <- subset(Whole, idents = c("5","7","8"))
cluster.averagesP0 <- AverageExpression(P0, return.seurat = TRUE)
DoHeatmap(cluster.averagesP0, features = c("Nanog","Sox2","Tdgf1","Etv5","Id4","Lhx1","Ret","Gfra1","Nanos2","Dnmt3a","Sohlh1"), size = 3, 
          draw.lines = FALSE)

cluster1.markers <- FindMarkers(Whole, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 20)
cluster0.markers <- FindMarkers(Whole, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 20)
cluster2.markers <- FindMarkers(Whole, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 20)
cluster3.markers <- FindMarkers(Whole, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 20)
cluster4.markers <- FindMarkers(Whole, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 20)
cluster5.markers <- FindMarkers(Whole, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 50)
