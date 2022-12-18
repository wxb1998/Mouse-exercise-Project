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

###Take an organ as an example###
Liver <- readRDS(file = "/05_Result/05_single_tissue/Liver/04_celltype/Combination_celltype.rds")
celltypes <- c("Hep","Hep-SC","Cho","EC","Kup","Mac","Neu","Bas","mDC",'pDC','proTC','CD4_Mem',"CD8_Naive",'CD8_Mem',"CD8_CTL",'CD4-CD8-TC',"NK","BC","Pla")
Liver@meta.data$celltype<-Liver@active.ident

Aging_DEGs <- data.frame()
Y_EX_DEGs <- data.frame()
O_EX_DEGs <- data.frame()

Idents(Liver) <- paste(Liver$celltype, Liver$sample, sep='_')
for (cell in celltypes){
  tmp <- FindMarkers(Liver, ident.1=paste0(cell, '_OC-Liver'), ident.2=paste0(cell, '_YC-Liver'))
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  tmp$compare <- 'Aging_DEGs'
  Aging_DEGs <- rbind(Aging_DEGs, tmp)
    print(paste0(cell, ' is finished'))
}
write.csv(Aging_DEGs,file=paste("/05_Result/05_single_tissue/Liver/05_DEGs/","Aging.csv",sep=""))
Aging_DEGs.deg <- subset(Aging_DEGs, p_val_adj<0.05 & abs(avg_log2FC)>0.25)
write.csv(Aging_DEGs.deg,file=paste("/05_Result/05_single_tissue/Liver/05_DEGs/","Aging_DEGs.csv",sep=""))

for (cell in celltypes){
  tmp <- FindMarkers(Liver, ident.1=paste0(cell, '_YE-Liver'), ident.2=paste0(cell, '_YC-Liver'))
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  tmp$compare <- 'Y_EX_DEGs'
  Y_EX_DEGs <- rbind(Y_EX_DEGs, tmp)
    print(paste0(cell, ' is finished'))
}
write.csv(Y_EX_DEGs,file=paste("/05_Result/05_single_tissue/Liver/05_DEGs/","Y_EX.csv",sep=""))
Y_EX_DEGs.deg <- subset(Y_EX_DEGs, p_val_adj<0.05 & abs(avg_log2FC)>0.25)
write.csv(Y_EX_DEGs.deg,file=paste("/05_Result/05_single_tissue/Liver/05_DEGs/","Y_EX_DEGs.csv",sep=""))

for (cell in celltypes){
  tmp <- FindMarkers(Liver, ident.1=paste0(cell, '_OE-Liver'), ident.2=paste0(cell, '_OC-Liver'))
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  tmp$compare <- 'O_EX_DEGs'
  O_EX_DEGs <- rbind(O_EX_DEGs, tmp)
    print(paste0(cell, ' is finished'))
}
write.csv(O_EX_DEGs,file=paste("/05_Result/05_single_tissue/Liver/05_DEGs/","O_EX.csv",sep=""))
O_EX_DEGs.deg <- subset(O_EX_DEGs, p_val_adj<0.05 & abs(avg_log2FC)>0.25)
write.csv(O_EX_DEGs.deg,file=paste("/05_Result/05_single_tissue/Liver/05_DEGs/","O_EX_DEGs.csv",sep=""))
