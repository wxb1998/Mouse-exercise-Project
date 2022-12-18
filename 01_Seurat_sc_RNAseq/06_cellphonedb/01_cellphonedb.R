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

for (i in c("Blood", "BM", "Intestine", "Liver","Lung", "Spleen", "Testis", "Kidney","SC","Heart","Brain","Muscle","Cerebellum")){ 
tissue <- readRDS(file = paste0("/05_Result/05_single_tissue/",i,"/04_celltype/Combination_celltype.rds"))
###Take an group as an sample###
seu.sub <- subset(tissue, sample %in% c(paste0("YC-",i)))
aa <- as.matrix(seu.sub@assays$RNA@data)
rownames(aa) <- toupper(rownames(aa))
write.table(aa, paste0('/02_Result/01_mouse_ex/05_single_tissue/02_cell_co/',i,'/01_YC/cellphonedb_count.txt'), sep='\t', quote=F)
meta_data <- cbind(rownames(seu.sub@meta.data), seu.sub@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
#meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, paste0('/02_Result/01_mouse_ex/05_single_tissue/02_cell_co/',i,'/01_YC/cellphonedb_meta.txt'), sep='\t', quote=F, row.names=F)
}
