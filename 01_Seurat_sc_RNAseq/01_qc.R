library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(SeuratData)
library(patchwork)
set.seed(2)

exercise_dir <- '/04_Raw_data/05_exercise_mice/mapping_jiace/'
read_in_sample_50 <- function(i){
  tmp <- Read10X(paste0(exercise_dir, i, '/outs/filtered_feature_bc_matrix/'))
  tmp <- CreateSeuratObject(counts = tmp, min.cells = 3, min.features = 200, project = i)
  tmp$sample <- i
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^mt-")
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_vln.pdf"), plot = p, width = 10, height = 6)
  plot1 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- plot1 + plot2
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_Scatter.pdf"),plot = plot3, width = 8, height = 4)
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_before_qc.rds"))
  
  tmp <- subset(tmp, subset = percent.mt < 50 & nFeature_RNA > 500)
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_qc_vln.pdf"), plot = p, width = 10, height = 6)
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_after_qc.rds"))
  tmp <- RenameCells(tmp, add.cell.id = i)
}

read_in_sample_20 <- function(i){
  tmp <- Read10X(paste0(exercise_dir, i, '/outs/filtered_feature_bc_matrix/'))
  tmp <- CreateSeuratObject(counts = tmp, min.cells = 3, min.features = 200, project = i)
  tmp$sample <- i
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^mt-")
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_vln.pdf"), plot = p, width = 10, height = 6)
  plot1 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- plot1 + plot2
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_Scatter.pdf"),plot = plot3, width = 8, height = 4)
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_before_qc.rds"))
  
  tmp <- subset(tmp, subset = percent.mt < 20 & nFeature_RNA > 500)
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_qc_vln.pdf"), plot = p, width = 10, height = 6)
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_after_qc.rds"))
  tmp <- RenameCells(tmp, add.cell.id = i)
}

read_in_sample_10 <- function(i){
  tmp <- Read10X(paste0(exercise_dir, i, '/outs/filtered_feature_bc_matrix/'))
  tmp <- CreateSeuratObject(counts = tmp, min.cells = 3, min.features = 200, project = i)
  tmp$sample <- i
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^mt-")
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_vln.pdf"), plot = p, width = 10, height = 6)
  plot1 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- plot1 + plot2
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_Scatter.pdf"),plot = plot3, width = 8, height = 4)
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_before_qc.rds"))
  
  tmp <- subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 500)
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_qc_vln.pdf"), plot = p, width = 10, height = 6)
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_after_qc.rds"))
  tmp <- RenameCells(tmp, add.cell.id = i)
}

read_in_sample_2_5 <- function(i){
  tmp <- Read10X_h5(paste0("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/02_cellbender/",j,i,"_remove_background_raw_feature_bc_matrix_filtered.h5"), use.names = TRUE, unique.features = TRUE)
  tmp <- CreateSeuratObject(counts = tmp, min.cells = 3, min.features = 200, project = i)
  tmp$sample <- i
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^mt-")
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_vln.pdf"), plot = p, width = 10, height = 6)
  plot1 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- plot1 + plot2
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_Scatter.pdf"),plot = plot3, width = 8, height = 4)
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_before_qc.rds"))
  
  tmp <- subset(tmp, subset = percent.mt < 2.5 & nFeature_RNA > 200)
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_qc_vln.pdf"), plot = p, width = 10, height = 6)
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_after_qc.rds"))
  tmp <- RenameCells(tmp, add.cell.id = i)
}

sample_normalize <- function(tmp, i){
  print(date())
  print(paste0(i, ': SCTransform started'))
  tmp <- SCTransform(tmp, verbose = FALSE)
  print(date())
  print(paste0(i, ': SCTransform finished'))

  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_SCT.rds"))

  tmp <- RunPCA(tmp, verbose=F)

  pdf(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_pca_heatmap.pdf"), width=10,height=20)
  DimHeatmap(tmp, dims=1:30, cells=500, balanced=T)
  dev.off()

  p<- ElbowPlot(tmp, ndims = 30)
  pdf(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_ElbowPlot.pdf"), height = 6, width = 7)
  print(p)
  dev.off()
  
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_bfPCR.rds"))

  return(tmp)
  print(paste0(i ,' completed'))
}

for (sample in c('OC-Kidney', 'OE-Kidney', 'YC-Kidney', 'YE-Kidney')){
  tmp <- read_in_sample_50(sample)
  assign(sample, tmp)
}

for (sample in c('OC-Artery', 'OE-Artery', 'YC-Artery', 'YE-Artery', 'OC-Liver', 'OE-Liver', 'YC-Liver', 'YE-Liver', 'OC-Testis', 'OE-Testis', 'YC-Testis', 'YE-Testis')){
  tmp <- read_in_sample_20(sample)
  assign(sample, tmp)
}

for (sample in c('OC-BM', 'OC-LPS-BM', 'OE-BM', 'OE-LPS-BM', 'YC-BM', 'YC-LPS-BM', 'YE-BM', 'YE-LPS-BM', 'OC-Intestine', 'OE-Intestine', 'YC-Intestine', 'YE-Intestine', 'OC-LPS-Lung', 'OC-Lung', 'OE-LPS-Lung', 'OE-Lung', 'YC-LPS-Lung', 'YC-Lung', 'YE-LPS-Lung', 'YE-Lung', 'OC-Blood', 'OC-LPS-Blood', 'OE-Blood', 'OE-LPS-Blood', 'YC-Blood', 'YC-LPS-Blood', 'YE-Blood', 'YE-LPS-Blood', 'OC-Spleen', 'OE-Spleen', 'YC-Spleen', 'YE-Spleen')){
  tmp <- read_in_sample_10(sample)
  assign(sample, tmp)
}

for (sample in c('OC-SC', 'OC-Brain', 'OC-Cerebellum', 'OC-Heart' ,'OC-Muscle', 'YC-SC' ,'YC-Brain', 'YC-Cerebellum', 'YC-Heart', 'YC-Muscle', 'YE-SC', 'YE-Brain', 'YE-Cerebellum', 'YE-Heart', 'YE-Muscle', 'OE-SC', 'OE-Brain', 'OE-Cerebellum', 'OE-Heart' ,'OE-Muscle','YC-Aorta-sn','YC-LPS-Aorta-sn','YC-Liver-sn','YC-LPS-Liver-sn','YE-Aorta-sn','YE-LPS-Aorta-sn','YE-Liver-sn','YE-LPS-Liver-sn','OC-Aorta-sn','OC-LPS-Aorta-sn','OC-Liver-sn','OC-LPS-Liver-sn','OE-Aorta-sn','OE-LPS-Aorta-sn','OE-Liver-sn','OE-LPS-Liver-sn')){
  tmp <- read_in_sample_2_5(sample)
  assign(sample, tmp)
}

samples <- ('OC-Kidney', 'OE-Kidney', 'YC-Kidney', 'YE-Kidney',
       'OC-Artery', 'OE-Artery', 'YC-Artery', 'YE-Artery', 'OC-Liver', 'OE-Liver', 'YC-Liver', 'YE-Liver', 'OC-Testis', 'OE-Testis', 'YC-Testis', 'YE-Testis',
        'OC-BM', 'OC-LPS-BM', 'OE-BM', 'OE-LPS-BM', 'YC-BM', 'YC-LPS-BM', 'YE-BM', 'YE-LPS-BM', 'OC-Intestine', 'OE-Intestine', 'YC-Intestine', 'YE-Intestine', 'OC-LPS-Lung', 'OC-Lung', 'OE-LPS-Lung', 'OE-Lung', 'YC-LPS-Lung', 'YC-Lung', 'YE-LPS-Lung', 'YE-Lung', 'OC-Blood', 'OC-LPS-Blood', 'OE-Blood', 'OE-LPS-Blood', 'YC-Blood', 'YC-LPS-Blood', 'YE-Blood', 'YE-LPS-Blood', 'OC-Spleen', 'OE-Spleen', 'YC-Spleen', 'YE-Spleen',
        'OC-SC', 'OC-Brain', 'OC-Cerebellum', 'OC-Heart' ,'OC-Muscle', 'YC-SC' ,'YC-Brain', 'YC-Cerebellum', 'YC-Heart', 'YC-Muscle', 'YE-SC', 'YE-Brain', 'YE-Cerebellum', 'YE-Heart', 'YE-Muscle', 'OE-SC', 'OE-Brain', 'OE-Cerebellum', 'OE-Heart' ,'OE-Muscle','YC-Aorta-sn','YC-LPS-Aorta-sn','YC-Liver-sn','YC-LPS-Liver-sn','YE-Aorta-sn','YE-LPS-Aorta-sn','YE-Liver-sn','YE-LPS-Liver-sn','OC-Aorta-sn','OC-LPS-Aorta-sn','OC-Liver-sn','OC-LPS-Liver-sn','OE-Aorta-sn','OE-LPS-Aorta-sn','OE-Liver-sn','OE-LPS-Liver-sn'
        )

for (sample in samples){
  tmp <- sample_normalize(get(sample), sample)
  assign(sample, tmp)
}

for (i in seq(1:84)){
  tmp <- readRDS(file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",samples[i],"/",samples[i],"_bfPCR.rds"))
  tmp <- RunUMAP(tmp, dims = dims[[i]], verbose=F)
  tmp <- FindNeighbors(tmp, reduction = "pca", dims = dims[[i]])
  tmp <- FindClusters(tmp, res=0.8)
  tmp[["cluster"]] <- Idents(tmp)
  UMAP <- DimPlot(object = tmp, reduction = "umap", label = TRUE)
  ggsave(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",samples[i],"/",samples[i],"_umap.pdf"), plot = UMAP, width = 8, height = 6)
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",samples[i],"/",samples[i],"_PCR.rds"))
  assign(samples[i], tmp)
  print(paste0(samples[i], ' completed'))
}

# Identify doublets ####
pK.df <- data.frame(matrix(nrow=0, ncol=3))
colnames(pK.df) <- c("Sample", "Optimal_pK","After_qc")

for (i in c(1:84)){
  sweep.res.list <- paramSweep_v3(get(samples[i]), PCs = dims[[i]], sct = T, num.cores=5)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- arrange(bcmvn, desc(BCmetric))$pK[1]
  tmp <- data.frame(Sample=samples[i], Optimal_pK=pK,After_qc=length(get(samples[i])@meta.data$orig.ident))
  pK.df <- rbind(pK.df, tmp)
  print(bcmvn)
  print(paste0("--------------", samples[i], " completed (", i, "/48)--------------"))
}
write.table(pK.df,file ="/05_Result/05_exercise_mice/seurat/01_qc/Optimal_pK.txt", sep = "\t" )

doublet.prop <- data.frame(matrix(nrow=0, ncol=3))
colnames(doublet.prop) <- c("Sample", "Number","Doublet_prop")

ratio  <- read.table("/05_Result/05_exercise_mice/seurat/01_qc/sample_doublet_ratio.txt")

for (i in seq(1:84)) {
  tmp <- readRDS(file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",samples[i],"/",samples[i],"_PCR.rds"))
  pK.use <- as.numeric(as.character(pK.df$Optimal_pK[i]))
  homotypic.prop <- modelHomotypic(tmp@meta.data$cluster)
  nExp_poi <- round(ratio[i]*length(tmp@meta.data$orig.ident))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = dims[[i]], pN = 0.25, pK = pK.use, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
  tmp[["doublet"]] <- tmp[[paste("DF.classifications_0.25", pK.use, nExp_poi.adj, sep="_")]]
  prop <- nExp_poi.adj/length(tmp@meta.data$cluster)
  prop.tmp <- data.frame(Sample=samples[i], Number=nExp_poi.adj, Doublet_prop=prop)
  doublet.prop <- rbind(doublet.prop, prop.tmp)
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",samples[i],"/",samples[i],"_doublets.rds"))
  assign(samples[i], tmp)
  print(paste0("--------------", samples[i], " completed (", i, "/32)--------------"))
}

write.table(doublet.prop,file ="/05_Result/05_exercise_mice/seurat/01_qc/doublet_prop_2.txt", sep = "\t" )

for (i in samples){
  doub <- DimPlot(get(i), group.by='doublet', cols=c('firebrick', 'grey90'))
  pdf(paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_UMAP_doublet.pdf"), height = 6, width = 8)
  print(doub)
  dev.off()
}




for (i in ('OC-Testis', 'OE-Testis', 'YC-Testis', 'YE-Testis','OC-Blood', 'OC-LPS-Blood', 'OE-Blood', 'OE-LPS-Blood', 'YC-Blood', 'YC-LPS-Blood', 'YE-Blood', 'YE-LPS-Blood', 'OC-Spleen', 'OE-Spleen', 'YC-Spleen', 'YE-Spleen')){
  tmp <- get(i)
  tmp <- subset(tmp, doublet=='Singlet')
  tmp <- subset(tmp, subset = nFeature_RNA < 6000 )
  tmp <- SCTransform(tmp, verbose = FALSE)
  assign(i, tmp)
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_final.rds"))
   print(paste0(i, ' completed'))
}

for (i in ('OC-Kidney', 'OE-Kidney', 'YC-Kidney', 'YE-Kidney',
       'OC-Artery', 'OE-Artery', 'YC-Artery', 'YE-Artery', 'OC-Liver', 'OE-Liver', 'YC-Liver', 'YE-Liver',
        'OC-BM', 'OC-LPS-BM', 'OE-BM', 'OE-LPS-BM', 'YC-BM', 'YC-LPS-BM', 'YE-BM', 'YE-LPS-BM', 'OC-Intestine', 'OE-Intestine', 'YC-Intestine', 'YE-Intestine', 'OC-LPS-Lung', 'OC-Lung', 'OE-LPS-Lung', 'OE-Lung', 'YC-LPS-Lung', 'YC-Lung', 'YE-LPS-Lung', 'YE-Lung',
        'OC-SC', 'OC-Brain', 'OC-Cerebellum', 'OC-Heart' ,'OC-Muscle', 'YC-SC' ,'YC-Brain', 'YC-Cerebellum', 'YC-Heart', 'YC-Muscle', 'YE-SC', 'YE-Brain', 'YE-Cerebellum', 'YE-Heart', 'YE-Muscle', 'OE-SC', 'OE-Brain', 'OE-Cerebellum', 'OE-Heart' ,'OE-Muscle','YC-Aorta-sn','YC-LPS-Aorta-sn','YC-Liver-sn','YC-LPS-Liver-sn','YE-Aorta-sn','YE-LPS-Aorta-sn','YE-Liver-sn','YE-LPS-Liver-sn','OC-Aorta-sn','OC-LPS-Aorta-sn','OC-Liver-sn','OC-LPS-Liver-sn','OE-Aorta-sn','OE-LPS-Aorta-sn','OE-Liver-sn','OE-LPS-Liver-sn'
        )){
  tmp <- get(i)
  tmp <- subset(tmp, doublet=='Singlet')
  tmp <- subset(tmp, subset = nFeature_RNA < 6000 )
  tmp <- SCTransform(tmp, verbose = FALSE)
  assign(i, tmp)
  saveRDS(tmp, file = paste0("/05_Result/05_exercise_mice/seurat/01_qc/",i,"/",i,"_final.rds"))
   print(paste0(i, ' completed'))
}
