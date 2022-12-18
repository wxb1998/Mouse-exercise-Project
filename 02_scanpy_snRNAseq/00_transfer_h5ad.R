library(scater)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(patchwork)
library(hdf5r)
library(data.table)
dir.create('/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/02_cellbender/scanpy')
setwd('/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/02_cellbender/scanpy')

group <-c('OC-','OE-','YC-','YE-')
sample<- c('cerebellum','SC','heart','muscle','brain')
for (i in 1:length(sample)){
  i <-sample[i]
  for (j in 1:length(group)) {
    j <- group[j]
dir.create(paste0('/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/',j,i))
setwd(paste0('/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/',j,i))
file <- Read10X_h5(paste0("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/02_cellbender/",j,i,"_remove_background_raw_feature_bc_matrix_filtered.h5"), use.names = TRUE, unique.features = TRUE)
tmp <- CreateSeuratObject(counts = file, project =paste0(j,i))
ct=GetAssayData(object = tmp, assay = "RNA", slot = "counts") 
ct=as.data.frame(ct)
fivenum(apply(ct,2,function(x) sum(x>0)))
table(ct>0)
write.table(data.frame(rownames(ct),rownames(ct)),file = paste0(j,i,'_genes.tsv'),
            quote = F,sep = '\t',
            col.names = F,row.names = F)
write.table(colnames(ct),file = paste0(j,i,'_barcodes.tsv'),quote = F,
            col.names = F,row.names = F)
##首先写一个头信息
file=paste0(j,i,'_matrix.mtx')
sink(file)
cat("%%MatrixMarket matrix coordinate integer general\n")
cat("%\n")
cat(paste(nrow(ct),ncol(ct),sum(ct>0),"\n")) 
sink()
###再写入表达量信息
tmp1=ct[1:5,1:4]
tmp1
tmp1=do.call(rbind,lapply(1:ncol(ct),function(m){
  return(data.frame(row=1:nrow(ct),
                    col=m,
                    exp=ct[,m]))
}) )
tmp1=tmp1[tmp1$exp>0,]
head(tmp1)
write.table(tmp1,file = paste0(j,i,'_matrix.mtx'),quote = F,
            col.names = F,row.names = F,append = T )
SaveH5Seurat(tmp, filename = paste0('tmp.h5Seurat'))
Convert(paste0('tmp.h5Seurat'), dest = paste0(j,i,'.h5ad'))
unlink('tmp.h5Seurat', recursive = FALSE, force = FALSE)
  }
  }
