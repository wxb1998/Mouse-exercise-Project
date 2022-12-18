library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(SeuratData)
library(patchwork)
library(Rserve)
set.seed(1234)
pc.num=c(1:24)

###Take an organ as an example###
samples <- c('OC-Liver', 'OE-Liver', 'YC-Liver', 'YE-Liver')
for (i in c(1:4)) {
 tmp <- readRDS(file = paste0("/05_Result/05_single_tissue/Liver/01_qc/",samples[i],"/",samples[i],"_final.rds"))
 assign(samples[i],tmp)
}
int.list <- list(get(samples[1]),get(samples[2]),get(samples[3]),get(samples[4]))
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 20000 * 1024^2)
int.features <- SelectIntegrationFeatures(object.list = int.list, nfeatures = 3000)
int.list <- PrepSCTIntegration(object.list = int.list, anchor.features = int.features, verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT",
                                      anchor.features = int.features, verbose = FALSE)
saveRDS(int.anchors, '/05_Result/05_single_tissue/Liver/02_int/int_anchors.rds')
plan("multiprocess", workers = 1)
options(future.globals.maxSize = 20000 * 1024^2)
exercise <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)
DefaultAssay(exercise) <- "integrated"
saveRDS(exercise, '/05_Result/05_single_tissue/Liver/02_int/int_seu.rds')
# clustering ####
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 20000 * 1024^2)
exercise <- RunPCA(exercise, verbose = FALSE)
pdf('/05_Result/05_single_tissue/Liver/02_int/pca_heatmap.pdf', height=20, width=10)
DimHeatmap(exercise, dims = 1:50, cells = 500, balanced = TRUE)
dev.off()
pdf('/05_Result/05_single_tissue/Liver/02_int/pca_elbow.pdf', height=5, width=8)
p=ElbowPlot(exercise, ndims=50)
print(p)
dev.off()
saveRDS(exercise, file = '/05_Result/05_single_tissue/Liver/02_int/beforePC.rds')
exercise <- RunUMAP(exercise, reduction = "pca", dims = 1:50)
pdf('/05_Result/05_single_tissue/Liver/02_int/umap_sample.pdf',width=8,height=6)
DimPlot(exercise, reduction = "umap", label = TRUE)
DimPlot(exercise, reduction = "umap", pt.size = 0.4, group.by = "sample")
dev.off()

###Determine the pc###
exercise <- readRDS(file = '/05_Result/05_single_tissue/Liver/02_int/beforePC.rds')

exercise <- RunUMAP(exercise, dims = pc.num,spread=3)# min.dist=0.2,
exercise <- FindNeighbors(exercise, reduction = "pca", dims = pc.num)
exercise <- FindClusters(exercise, resolution = 2.2)#2.0

saveRDS(exercise, file = '/05_Result/05_single_tissue/Liver/02_int/umap.rds')

plot <- DimPlot(exercise, reduction = "umap", label = TRUE,pt.size = 0.1,raster = FALSE)
ggsave("/05_Result/05_single_tissue/Liver/02_int/umap.pdf", plot = plot, width = 8, height = 6)

pdf('/05_Result/05_single_tissue/Liver/02_int/umap_sample_single.pdf',width=28,height=12)
DimPlot(exercise, reduction = "umap", split.by = "sample",pt.size = 0.1,label.size = 6,ncol=4,raster = FALSE)
dev.off()

Tissue <- exercise
saveRDS(Tissue, file = '/05_Result/05_single_tissue/Liver/02_int/umap.rds')

pdf('/05_Result/05_single_tissue/Liver/02_int/umap_group.pdf',width=8,height=6)
DimPlot(Tissue, reduction = "umap", pt.size = 0.1, group.by = "sample",raster=FALSE, cols= c('#E88B88','#D8281D','#45549A','#0B3582'))
DimPlot(Tissue, reduction = "umap", pt.size = 0.1, group.by = "sample",raster=FALSE, cols= c('#E88B88',"NA","NA","NA"))
DimPlot(Tissue, reduction = "umap", pt.size = 0.1, group.by = "sample", raster=FALSE,cols= c("NA",'#D8281D',"NA","NA"))
DimPlot(Tissue, reduction = "umap", pt.size = 0.1, group.by = "sample", raster=FALSE,cols= c("NA",'NA',"#45549A","NA"))
DimPlot(Tissue, reduction = "umap", pt.size = 0.1, group.by = "sample", raster=FALSE,cols= c("NA",'NA',"NA","#0B3582"))
dev.off()


