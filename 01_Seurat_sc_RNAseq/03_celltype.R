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
set.seed(2)


###Take an organ as an example###
Tissue <- readRDS(file = '/05_Result/05_single_tissue/Liver/02_int/umap.rds')
DefaultAssay(Tissue) <- "RNA"
Tissue <- NormalizeData(object = Tissue, normalization.method = "LogNormalize")

celltypes <- c("Hep","Hep-SC","Cho","EC","Kup","Mac","Neu","Bas","mDC",'pDC','proTC','CD4_Mem',"CD8_Naive",'CD8_Mem',"CD8_CTL",'CD4-CD8-TC',"NK","BC","Pla")
C1 <- c(29)
C2 <- c(48)
C3 <- c(26,37)
C4 <- c(3,5,6,7,8,11,13,15,16,19,20,22,23,24,27,39,45,46,52)
C5 <- c(9)
C6 <- c(14,18,21,28)

C7 <- c(4,30,47,53)
C8 <- c(49)
C9 <- c(25,50)
C10 <- c(35)
C11 <- c(41)

C12 <- c(38)
C13 <- c(2,40,42,44)
C14 <- c(12)
C15 <- c(17,31,51)
C16 <- c(32)

C17 <- c(34,36)
C18 <- c(0,1,10,43,54)
C19 <- c(33)

for (i in 1:19) {
  m <- paste0('C',i)
  tmp <- get(m)
  assign(celltypes[i],tmp)
}

new.cluster.ids <- c(rep(NA,length(unique(Idents(Tissue))  )))
cluster <- function(celltypes,new.cluster.ids){
  for (cell in celltypes) {
    tmp <- get(cell)
    for (i in tmp) {
      new.cluster.ids[i+1] <- cell
    }
  }
  return(new.cluster.ids)

}
new.cluster.ids <- cluster(celltypes,new.cluster.ids)

Tissue <- RenameIdents(Tissue,c("0"=new.cluster.ids[1], "1"=new.cluster.ids[2], "2"=new.cluster.ids[3], "3"=new.cluster.ids[4],"4"= new.cluster.ids[5], "5"=new.cluster.ids[6],"6"=new.cluster.ids[7], "7"=new.cluster.ids[8],"8"= new.cluster.ids[9],"9"= new.cluster.ids[10],"10"= new.cluster.ids[11],"11"=new.cluster.ids[12], "12"=new.cluster.ids[13], "13"=new.cluster.ids[14],"14"= new.cluster.ids[15],"15"=new.cluster.ids[16],"16"=new.cluster.ids[17], "17"=new.cluster.ids[18],"18"= new.cluster.ids[19],"19"= new.cluster.ids[20],"20"= new.cluster.ids[21],"21"=new.cluster.ids[22], "22"=new.cluster.ids[23], "23"=new.cluster.ids[24],"24"= new.cluster.ids[25],"25"=new.cluster.ids[26],"26"=new.cluster.ids[27], "27"=new.cluster.ids[28],"28"= new.cluster.ids[29],"29"= new.cluster.ids[30],"30"= new.cluster.ids[31],"31"=new.cluster.ids[32], "32"=new.cluster.ids[33], "33"=new.cluster.ids[34],"34"= new.cluster.ids[35],"35"=new.cluster.ids[36],"36"=new.cluster.ids[37], "37"=new.cluster.ids[38],"38"= new.cluster.ids[39],"39"= new.cluster.ids[40],"40"= new.cluster.ids[41],"41"=new.cluster.ids[42], "42"=new.cluster.ids[43], "43"=new.cluster.ids[44],"44"= new.cluster.ids[45], "45"=new.cluster.ids[46],"46"=new.cluster.ids[47], "47"=new.cluster.ids[48],"48"= new.cluster.ids[49],"49"= new.cluster.ids[50],"50"= new.cluster.ids[51],"51"=new.cluster.ids[52],"52"=new.cluster.ids[53],"53"=new.cluster.ids[54],"54"= new.cluster.ids[55]))

Idents(Tissue) <- factor(Idents(Tissue), levels=celltypes)
Tissue$celltype <- Idents(Tissue)
CT.col <- c('#d5695d','#b1283a',"#be9c2e",'#d3ba68',"#fb832d",'#adadad','#59a77f',"#098154",'#088158',"#5d8ca8","#65a479",'#00887d',"#098154", "#65a479","#00A087B2","#76c0c1","#5d8ca8",'#6794a7','#016392') 

UMAP <- DimPlot(Tissue, reduction = "umap", cols=CT.col, label = TRUE, pt.size = 0.1,label.size = 3.0,raster = FALSE) 
ggsave("/05_Result/05_single_tissue/Liver/04_celltype/celltype_UMAP.pdf", plot = UMAP, width = 8, height = 6)
UMAP <- DimPlot(Tissue, reduction = "umap", cols=CT.col, label = FALSE, pt.size = 0.1,raster = FALSE) 
ggsave("/05_Result/05_single_tissue/Liver/04_celltype/celltype_UMAP2.pdf", plot = UMAP, width = 8, height = 6)

markers <- c("Apoc3", "Alb","Dcn","Col1a1","Sox9","Spp1","Krt18","Stab2","Kdr","Pecam1","Adgre1","Clec4f","Cd163","Cd68","Csf1r",'Itgam',"S100a4","S100a9","S100a8","Cpa3","Ms4a2",'Clec9a','Cd74','Pld4',"Irf8","Cd3e",'Mki67','Cd4','Ccr7','Lef1','Prdm1','Pdcd1','Cd8b1','Cd8a','Gzmk',"Gzmb",'Nkg7','Prf1','Cd79a', 'Cd19','Ms4a1','Mzb1','Jchain')
Idents(Tissue) <- factor(Tissue$celltype, levels=rev(celltypes))
DotPlot<-DotPlot(Tissue,features= rev(markers), cols=c('grey90', '#C63C3C'))+RotatedAxis()
ggsave("/05_Result/05_single_tissue/Liver/04_celltype/DotPlot.pdf", plot = DotPlot, width = 13, height = 6)
saveRDS(Tissue, file = "/05_Result/05_single_tissue/Liver/04_celltype/Combination_celltype.rds")

##Analysis############################################################################
Tissue@meta.data$celltype<-Tissue@active.ident

write.csv(Tissue@meta.data,"/05_Result/05_single_tissue/Liver/04_celltype/Combination_al.csv", row.names =FALSE)

# cell proportion ####
samples <- c('OC-Liver', 'OE-Liver', 'YC-Liver', 'YE-Liver')
Tissue$sample <- factor(Tissue$sample, levels=samples)
cell_number=table(Idents(Tissue),Tissue$sample)
write.csv(cell_number,"/05_Result/05_single_tissue/Liver/04_celltype/cell_number_all_al.csv")
cell.prop <- data.frame(table(Tissue$celltype,Tissue$sample))
colnames(cell.prop) <- c('celltype', 'sample', 'number')
cell.prop <- cell.prop %>% group_by(sample) %>% mutate(total=sum(number))
cell.prop$prop <- cell.prop$number/cell.prop$total

cell.prop$sample <- factor(cell.prop$sample,
                           levels=samples, ordered=TRUE)
write.csv(cell.prop,"/05_Result/05_single_tissue/Liver/04_celltype/cell_prop_all_al.csv")
p=ggplot(cell.prop, aes(celltype, prop, fill=cell.prop$sample)) +
  geom_bar(stat='identity', position = 'dodge', width=0.9) +
  #scale_y_continuous(expand = c(0,0.005)) +
  #scale_fill_manual(values = c( 'dodgerblue1', 'orange2')) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'),
        axis.text = element_text(color='black'),axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave("/05_Result/05_single_tissue/Liver/04_celltype/cell_prop_al.pdf", plot = p, width = 24, height = 6)
