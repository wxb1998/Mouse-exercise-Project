library(SCENIC)
library(AUCell)
library(RcisTarget)
library(Seurat)
library(scater)
library(patchwork)
set.seed(2)

###Take an organ as an example###

#Aging#
for (i in c("Liver")){
seu.sub <- readRDS(file = paste0("/05_Result/05_single_tissue/",i,"/04_celltype/Combination_celltype.rds"))
seu.sub <- subset(seu.sub, sample %in% c(paste0("OC-",i),paste0("YC-",i)))
deg = read.csv(paste0("/05_Result/05_single_tissue/",i,"/05_DEGs/Aging_DEGs.csv"))
#tmp <- subset(seu.sub, celltype==i)
exp.mat <- GetAssayData(seu.sub, slot='data')
#deg_tmp=subset(deg, celltype==i)
genes <- unique(as.character(deg$gene))
exp.mat <- as.matrix(exp.mat)
exp.mat <- as.matrix(exp.mat[genes,])
setwd(paste0("/dellstorage02/quj_lab/pingjiale/02_Result/01_mouse_ex/05_single_tissue/all_fu_fu/01_Aging/",i,"/"))
org="mgi"
dbDir="/03_Database/02_SCENIC/mm9"
myDatasetTitle="SCENIC analysis of Aging"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)  ##nCores线程
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
# Gene filter
genesKept <- geneFiltering(exp.mat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exp.mat),
                           minSamples=ncol(exp.mat)*.01)
exprMat_filtered <- exp.mat[genesKept, ]

# Run Genie3
runCorrelation(exprMat_filtered, scenicOptions)
#judgement-------------------------------------------------------
      judgement = tryCatch({
        #correct pipeline
        runGenie3(exprMat_filtered, scenicOptions, nParts = 10)
        print(paste0('----------------------',i,'_Genie3 finished','----------------------'))
      },  error=function(e) e
      )
      if(inherits(judgement, 'error')) {
        print(paste0(i,'_error'))
        next}
      #---------------------------------------------------------------

# Run the remaining
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 15
#judgement-------------------------------------------------------
      judgement = tryCatch({
        #correct pipeline
        runSCENIC_1_coexNetwork2modules(scenicOptions)
        runSCENIC_2_createRegulons(scenicOptions)
        runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
        print(paste0('----------------------',i,'_data completed','----------------------'))
      },  error=function(e) e
      )
      if(inherits(judgement, 'error')) {
        print(paste0(i,'_error'))
        next}
      #---------------------------------------------------------------
}

#O_EX#
for (i in c("Aorta", "Blood", "BM", "Intestine", "Kidney","Liver", "Lung","Spleen", "Testis")){
seu.sub <- readRDS(file = paste0("/05_Result/05_single_tissue/",i,"/04_celltype/Combination_celltype.rds"))
seu.sub <- subset(seu.sub, sample %in% c(paste0("OE-",i),paste0("OC-",i)))
deg = read.csv(paste0("/05_Result/05_single_tissue/",i,"/05_DEGs/O_EX_DEGs.csv"))

#tmp <- subset(seu.sub, celltype==i)
exp.mat <- GetAssayData(seu.sub, slot='data')
#deg_tmp=subset(deg, celltype==i)
genes <- unique(as.character(deg$gene))
exp.mat <- as.matrix(exp.mat)
exp.mat <- as.matrix(exp.mat[genes,])

setwd(paste0("/dellstorage02/quj_lab/pingjiale/02_Result/01_mouse_ex/05_single_tissue/all_fu_fu/03_O_EX/",i,"/"))

org="mgi"
dbDir="/03_Database/02_SCENIC/mm9"
myDatasetTitle="SCENIC analysis of O_EX"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)  ##nCores线程

saveRDS(scenicOptions, file="int/scenicOptions.Rds")
# Gene filter
genesKept <- geneFiltering(exp.mat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exp.mat),
                           minSamples=ncol(exp.mat)*.01)
exprMat_filtered <- exp.mat[genesKept, ]

# Run Genie3
runCorrelation(exprMat_filtered, scenicOptions)
#judgement-------------------------------------------------------
      judgement = tryCatch({
        #correct pipeline
        runGenie3(exprMat_filtered, scenicOptions, nParts = 10)
        print(paste0('----------------------',i,'_Genie3 finished','----------------------'))
      },  error=function(e) e
      )
      if(inherits(judgement, 'error')) {
        print(paste0(i,'_error'))
        next}
      #---------------------------------------------------------------

# Run the remaining
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 15
#judgement-------------------------------------------------------
      judgement = tryCatch({
        #correct pipeline
        runSCENIC_1_coexNetwork2modules(scenicOptions)
        runSCENIC_2_createRegulons(scenicOptions)
        runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
        print(paste0('----------------------',i,'_data completed','----------------------'))
      },  error=function(e) e
      )
      if(inherits(judgement, 'error')) {
        print(paste0(i,'_error'))
        next}
      #---------------------------------------------------------------
}

#LPS#
for (i in c("Blood", "BM","Lung","Liver_lps","Aorta_lps")){
seu.sub <- readRDS(file = paste0("/05_Result/05_single_tissue/",i,"/04_celltype/Combination_celltype.rds"))
seu.sub <- subset(seu.sub, Group %in% c("YC","YC_LPS"))
deg = read.csv(paste0("/05_Result/05_single_tissue/",i,"/05_DEGs/Y_LPS_DEGs.csv"))

exp.mat <- GetAssayData(seu.sub, slot='data')
genes <- unique(as.character(deg$gene))
exp.mat <- as.matrix(exp.mat)
exp.mat <- as.matrix(exp.mat[genes,])

setwd(paste0("/dellstorage02/quj_lab/pingjiale/02_Result/01_mouse_ex/05_single_tissue/01_TFs/lps/04_Y_LPS/",i,"/"))

org="mgi"
dbDir="/03_Database/02_SCENIC/mm9"
myDatasetTitle="SCENIC analysis of Y_LPS"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)  ##nCores线程

saveRDS(scenicOptions, file="int/scenicOptions.Rds")
# Gene filter
genesKept <- geneFiltering(exp.mat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exp.mat),
                           minSamples=ncol(exp.mat)*.01)
exprMat_filtered <- exp.mat[genesKept, ]

# Run Genie3
runCorrelation(exprMat_filtered, scenicOptions)
#judgement-------------------------------------------------------
      judgement = tryCatch({
        #correct pipeline
        runGenie3(exprMat_filtered, scenicOptions, nParts = 10)
        print(paste0('----------------------',i,'_Genie3 finished','----------------------'))
      },  error=function(e) e
      )
      if(inherits(judgement, 'error')) {
        print(paste0(i,'_error'))
        next}
      #---------------------------------------------------------------

# Run the remaining
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 15
#judgement-------------------------------------------------------
      judgement = tryCatch({
        #correct pipeline
        runSCENIC_1_coexNetwork2modules(scenicOptions)
        runSCENIC_2_createRegulons(scenicOptions)
        runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
        print(paste0('----------------------',i,'_data completed','----------------------'))
      },  error=function(e) e
      )
      if(inherits(judgement, 'error')) {
        print(paste0(i,'_error'))
        next}
      #---------------------------------------------------------------
}


for (i in c("Blood", "BM","Lung","Liver_lps","Aorta_lps")){
seu.sub <- readRDS(file = paste0("/05_Result/05_single_tissue/",i,"/04_celltype/Combination_celltype.rds"))
seu.sub <- subset(seu.sub, Group %in% c("YE_LPS","YC_LPS"))
deg = read.csv(paste0("/05_Result/05_single_tissue/",i,"/05_DEGs/Y_EX_LPS_DEGs.csv"))

exp.mat <- GetAssayData(seu.sub, slot='data')

genes <- unique(as.character(deg$gene))
exp.mat <- as.matrix(exp.mat)
exp.mat <- as.matrix(exp.mat[genes,])

setwd(paste0("/dellstorage02/quj_lab/pingjiale/02_Result/01_mouse_ex/05_single_tissue/01_TFs/lps/05_Y_EX_LPS/",i,"/"))

org="mgi"
dbDir="/03_Database/02_SCENIC/mm9"
myDatasetTitle="SCENIC analysis of Y_EX_LPS"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)  ##nCores线程

saveRDS(scenicOptions, file="int/scenicOptions.Rds")
# Gene filter
genesKept <- geneFiltering(exp.mat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exp.mat),
                           minSamples=ncol(exp.mat)*.01)
exprMat_filtered <- exp.mat[genesKept, ]

# Run Genie3
runCorrelation(exprMat_filtered, scenicOptions)
#judgement-------------------------------------------------------
      judgement = tryCatch({
        #correct pipeline
        runGenie3(exprMat_filtered, scenicOptions, nParts = 10)
        print(paste0('----------------------',i,'_Genie3 finished','----------------------'))
      },  error=function(e) e
      )
      if(inherits(judgement, 'error')) {
        print(paste0(i,'_error'))
        next}
      #---------------------------------------------------------------

# Run the remaining
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 15
#judgement-------------------------------------------------------
      judgement = tryCatch({
        #correct pipeline
        runSCENIC_1_coexNetwork2modules(scenicOptions)
        runSCENIC_2_createRegulons(scenicOptions)
        runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
        print(paste0('----------------------',i,'_data completed','----------------------'))
      },  error=function(e) e
      )
      if(inherits(judgement, 'error')) {
        print(paste0(i,'_error'))
        next}
      #---------------------------------------------------------------
}
