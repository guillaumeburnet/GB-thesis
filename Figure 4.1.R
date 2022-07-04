if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase') ##make sure the right gfortran version is installed
devtools::install_github('cole-trapnell-lab/monocle3')

install.packages("pheatmap")
library(pheatmap)
library(monocle3)
library(dplyr)
library(Seurat)
library(Matrix)
library(cellwrangler)
library(SeuratWrappers)
install.packages("systemfonts")
library(systemfonts)
ylibrary(heatmap3)
remotes::install_github("AllanCameron/geomtextpath")
library(geomtextpath)
library (viridis)

cds <- load_mm_data(mat_path = "~/Desktop/GB/Nef dataset/Seurat/matrix.mtx", 
                    feature_anno_path = "~/Desktop/GB/Nef dataset/Seurat/genes.tsv", 
                    cell_anno_path = "~/Desktop/GB/Nef dataset/Seurat/barcodes.tsv",
                    feature_metadata_column_names=c("gene_short_name")) 

cds <- preprocess_cds(cds, num_dim = 30)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds) ##Default if PCA, maybe better with ICA?
plot_cells(cds, color_cells_by="SexAge")
plot_cells(cds, color_cells_by="Batch") ## Cells cluster according to Bacth, need to correct for this..
cds = align_cds(cds, num_dim = 30, alignment_group = "Batch")## This corrects for Batch effect
cds = align_cds(cds, num_dim = 30, alignment_group = "SexAge")
cds = reduce_dimension(cds)
plot_cells(cds, color_cells_by="Batch")
plot_cells(cds, color_cells_by="SexAge")##Makes more sense in terms of clustering
cds_subset <- choose_cells(cds)  ## remove outliers
plot_cells(cds_subset, color_cells_by="SexAge", label_cell_groups = FALSE, show_trajectory_graph = FALSE, cell_size=0.25 )
plot_cells(cds, color_cells_by="pseudotime", show_trajectory_graph = FALSE) + scale_color_viridis()

male <- c("XY E10","XY E11","XY E12","XY E13","XY E16")
maleE10 <- c("XY E10")
female <- c("XX E10","XX E11","XX E12","XX E13","XX E16")
femaleE10<- c("XX E10","XX E11")
cds_male <- cds[,colData(cds)$SexAge %in% male]
cds_maleE10 <- cds_male[,colData(cds_male)$SexAge%in%maleE10]
cellsE10 <-colnames(cds_maleE10)
cds_male <- choose_cells(cds_male)
cds_male <- cluster_cells(cds_male)
cds_male <- learn_graph(cds_male, use_partition=FALSE)
cds_male <- order_cells(cds_male)
cds_maleinitial<-cds_male

cds_female <- cds[,colData(cds)$SexAge %in% female]
cds_femaleE10 <- cds_female[,colData(cds_female)$SexAge%in%femaleE10]
cellsE10 <-colnames(cds_femaleE10)
cds_female <- choose_cells(cds_female)
cds_female <- cluster_cells(cds_female, use_partition=FALSE)
cds_female <- learn_graph(cds_female, use_partition=FALSE)
cds_female <- order_cells(cds_female)
cds_female <- order_cells(cds_female, root_cells = cellsE10)
saveRDS(cds_female, file = "~/Desktop/GB/Nef dataset/cds_female.rds")
saveRDS(cds_male, file = "~/Desktop/GB/Nef dataset/cds_male.rds")
cds_male<-readRDS("~/Desktop/GB/Nef dataset/cds_male.rds")
cds_female<-readRDS("~/Desktop/GB/Nef dataset/cds_female.rds")

###Plot gene expression across pseudotime

Nanog.expression.male <- exprs(cds_male)[match(c("ENSMUSG00000012396.12"),rownames(rowData(cds_male))),order(pseudotime(cds_male))]
Nanos2.expression.male <- exprs(cds_male)[match(c("ENSMUSG00000051965.8"),rownames(rowData(cds_male))),order(pseudotime(cds_male))]
Dnmt3l.expression.male <- exprs(cds_male)[match(c("ENSMUSG00000000730.13"),rownames(rowData(cds_male))),order(pseudotime(cds_male))]
Cripto.expression.male <- exprs(cds_male)[match(c("ENSMUSG00000032494.12"),rownames(rowData(cds_male))),order(pseudotime(cds_male))]
Cfc1.expression.male <- exprs(cds_male)[match(c("ENSMUSG00000026124.4"),rownames(rowData(cds_male))),order(pseudotime(cds_male))]
Nodal.expression.male <- exprs(cds_male)[match(c("ENSMUSG00000037171.6"),rownames(rowData(cds_male))),order(pseudotime(cds_male))]
Lefty1.expression.male <- exprs(cds_male)[match(c("ENSMUSG00000038793.8"),rownames(rowData(cds_male))),order(pseudotime(cds_male))]
Lefty2.expression.male <- exprs(cds_male)[match(c("ENSMUSG00000066652.5"),rownames(rowData(cds_male))),order(pseudotime(cds_male))]


PT <- seq(0,100, length.out=4592)

x=log(Nanog.expression.male+1)
y<-rep(c("Nanog"),times=4592)
Nanog.plotframe <- data.frame (PT,x,y)

x=log(Nanos2.expression.male+1)
y<-rep(c("Nanos2"),times=4592)
Nanos2.plotframe <- data.frame (PT,x,y)

x=log(Dnmt3l.expression.male+1)
x=x/2.5
y<-rep(c("Dnmt3l"),times=4592)
Dnmt3l.plotframe <- data.frame (PT,x,y)

x=log(Cripto.expression.male+1)
y<-rep(c("Cripto"),times=4592)
Cripto.plotframe <- data.frame (PT,x,y)

x=log(Cfc1.expression.male+1)
y<-rep(c("Cryptic"),times=4592)
Cryptic.plotframe <- data.frame (PT,x,y)

x=log(Nodal.expression.male+1)
y<-rep(c("Nodal"),times=4592)
Nodal.plotframe <- data.frame (PT,x,y)

x=log(Lefty1.expression.male+1)
x=x/4
y<-rep(c("Lefty1"),times=4592)
Lefty1.plotframe <- data.frame (PT,x,y)

x=log(Lefty2.expression.male+1)
x=x/4
y<-rep(c("Lefty2"),times=4592)
Lefty2.plotframe <- data.frame (PT,x,y)

plotframe<-rbind(Cripto.plotframe,Nodal.plotframe, Lefty1.plotframe,Lefty2.plotframe)
plotframe<-rbind(Cripto.plotframe,Cryptic.plotframe)
plotframe2<-rbind(Nanog.plotframe,Nanos2.plotframe,Dnmt3l.plotframe)
colnames(plotframe)<-c("Pseudotime","Log_Expression","Gene")
colnames(plotframe2)<-c("Pseudotime","Log_Expression","Gene")
ggplot(plotframe, aes(x=Pseudotime, y=Log_Expression , color=Gene, label=Gene)) + 
  geom_textsmooth(method="gam", straight=TRUE) + scale_y_continuous(name = "", sec.axis = sec_axis( trans=~.*4, name=""))

ggplot(plotframe2, aes(x=Pseudotime, y=Log_Expression , color=Gene, label=Gene)) + 
  geom_textsmooth(method="gam", straight=TRUE) + scale_y_continuous(name = "", sec.axis = sec_axis( trans=~.*2.5, name=""))

