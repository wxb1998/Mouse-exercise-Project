
library(DESeq2)
tissue <- c("Aorta","Blood","BM","Brain","Heart","Intestine","Kidney","Liver","Lung","Muscle",'Spinalcord','Spleen','Testis')
list <- c(68,68,65,68,68,72,69,68,67,69,73,69,69)

for (i in c(1:13)){
all_symbol=read.csv("E:/01_progrom/02_mouse_exercise/03_rnaseq_mouse/Mus_musculus_chr.GRCm38.91_GeneID_Symbol.txt",sep='\t')
data.path=paste0("E:/01_progrom/02_mouse_exercise/03_rnaseq_mouse/",tissue[i],"/01_pca/count/")
sample.path=paste0("E:/01_progrom/02_mouse_exercise/03_rnaseq_mouse/",tissue[i],"/06_OE/")
data.list=list.files(data.path,pattern = '*txt',full.names = T)
data.list

HTseq.handling=function(x,n){
  df=read.table(x,header = F)
  name=substr(x,n,nchar(x)-8)
  colnames(df)=c("ensembl_gene_id",name)
  df=df[-((nrow(df)-4):nrow(df)),]
  return(df)
}
##data.list[1:10] OEvsOC; data.list[11:20] YEvsYC; data.list[c(1:5,11:15)] OCvsYC
tt=lapply(data.list[c(1:6)], HTseq.handling,list[i])
head(tt[[1]])
tail(tt[[1]])
head(tt[[6]])
data.name=unique(substr(data.list,list[i],nchar(data.list[1])-10))
allmerge=function(x){
  merge.data=merge(x[[1]],x[[2]],by='ensembl_gene_id')
  for (i in 3:length(x)) {
    merge.data=merge(merge.data,x[[i]],by='ensembl_gene_id')
  }
  return(merge.data)
}
data.all=allmerge(tt)
#data.all=data.all[,c(1,5:7,2:4)]
write.csv(data.all,paste(sample.path,'merge_sample.csv',sep = ''),row.names =F)

data.combn=data.frame(c(data.name[1],data.name[2]))
all.name=as.character(data.combn[,1])
all.sample=c(paste(all.name[1],c(1,2,3),sep='-'),paste(all.name[2],c(1,2,3),sep='-'))
samplecondition=factor(c(rep(all.name[1],3),rep(all.name[2],3)),
                       levels = c(all.name[1],all.name[2]))
samplecondition 
mycol.data=data.frame(row.names = all.sample,condition=samplecondition)
mycol.data 
mycount.data=data.all
rownames(mycount.data)=mycount.data[,1]
mycount.data=mycount.data[,-1] 
class(mycount.data)

##dds
mydds=DESeqDataSetFromMatrix(countData = mycount.data,colData = mycol.data,design = ~ condition)
head(mycount.data)
head(mycol.data)
samplecondition

##??׼??
mydds=DESeq(mydds)
myres=results(mydds)
########rlog
rld=rlog(mydds)
write.csv(assay(rld),paste(sample.path,"sample_rlog.csv",sep=''),row.names = T)
myres=myres[order(myres$padj),]


all_gene.data=merge(as.data.frame(myres),as.data.frame(counts(mydds,normalize=TRUE)),by="row.names",sort=FALSE)
resdata=all_gene.data
colnames(all_symbol)[1]="gene_id"
colnames(resdata)[1]="gene_id"
resdata_anno=merge(all_symbol,resdata,by="gene_id")
write.csv(resdata_anno,paste(sample.path,"all_gene.csv",sep=''),row.names = F)
######diff_gene
diff_gene=subset(resdata_anno,padj<0.05&(log2FoldChange>=1|log2FoldChange<=-1))
write.csv(diff_gene,paste(sample.path,"diff_gene.csv",sep=''),row.names = F)
#####up_gene
up_gene=subset(diff_gene,log2FoldChange>=1)
write.csv(up_gene,paste(sample.path,"up_gene.csv",sep=''),row.names = F)
#####down_gene
down_gene=subset(diff_gene,log2FoldChange<=-1)
write.csv(down_gene,paste(sample.path,'down_gene.csv',sep=''),row.names = F)
}



library(DESeq2)
tissue <- c("Aorta","Liver","BM","Lung")
list <- c(67,67,64,66)

for (i in c(1:4)){
all_symbol=read.csv("E:/01_progrom/02_mouse_exercise/03_rnaseq_mouse/Mus_musculus_chr.GRCm38.91_GeneID_Symbol.txt",sep='\t')
data.path=paste0("E:/01_progrom/02_mouse_exercise/03_rnaseq_mouse/",tissue[i],"/02_lps/count/")
sample.path=paste0("E:/01_progrom/02_mouse_exercise/03_rnaseq_mouse/",tissue[i],"/09_O_LPS/")
data.list=list.files(data.path,pattern = '*txt',full.names = T)
data.list

HTseq.handling=function(x,n){
  df=read.table(x,header = F)
  name=substr(x,n,nchar(x)-8)
  colnames(df)=c("ensembl_gene_id",name)
  df=df[-((nrow(df)-4):nrow(df)),]
  return(df)
}
##YC:data.list[10:15] YE: data.list[13:18] OC:data.list[1:6] OE: data.list[4:9]
tt=lapply(data.list[1:6], HTseq.handling,list[i])
head(tt[[1]])
tail(tt[[1]])
head(tt[[6]])
data.name=unique(substr(data.list,list[i],nchar(data.list)-10))
allmerge=function(x){
  merge.data=merge(x[[1]],x[[2]],by='ensembl_gene_id')
  for (i in 3:length(x)) {
    merge.data=merge(merge.data,x[[i]],by='ensembl_gene_id')
  }
  return(merge.data)
}
data.all=allmerge(tt)
#data.all=data.all[,c(1,7:10,2:6)]##???? ??????ǰ ʵ???ں?
write.csv(data.all,paste(sample.path,'merge_sample.csv',sep = ''),row.names =F)

data.combn=data.frame(c(data.name[1],data.name[2]))
all.name=as.character(data.combn[,1])
all.sample=c(paste(all.name[1],c(1,2,3),sep='-'),paste(all.name[2],c(1,2,3),sep='-'))
samplecondition=factor(c(rep(all.name[1],3),rep(all.name[2],3)),
                       levels = c(all.name[1],all.name[2]))
samplecondition 
mycol.data=data.frame(row.names = all.sample,condition=samplecondition)
mycol.data 
mycount.data=data.all
rownames(mycount.data)=mycount.data[,1]
mycount.data=mycount.data[,-1] ##????????????????Ʒ??
class(mycount.data)

##????dds????
mydds=DESeqDataSetFromMatrix(countData = mycount.data,colData = mycol.data,design = ~ condition)
head(mycount.data)
head(mycol.data)
samplecondition

##??׼??
mydds=DESeq(mydds)
myres=results(mydds)
########rlog
rld=rlog(mydds)
write.csv(assay(rld),paste(sample.path,"sample_rlog.csv",sep=''),row.names = T)
myres=myres[order(myres$padj),]

##??ȡ????????????
all_gene.data=merge(as.data.frame(myres),as.data.frame(counts(mydds,normalize=TRUE)),by="row.names",sort=FALSE)
resdata=all_gene.data
colnames(all_symbol)[1]="gene_id"
colnames(resdata)[1]="gene_id"
resdata_anno=merge(all_symbol,resdata,by="gene_id")
write.csv(resdata_anno,paste(sample.path,"all_gene.csv",sep=''),row.names = F)
######diff_gene
diff_gene=subset(resdata_anno,padj<0.05&(log2FoldChange>=1|log2FoldChange<=-1))
write.csv(diff_gene,paste(sample.path,"diff_gene.csv",sep=''),row.names = F)
#####up_gene
up_gene=subset(diff_gene,log2FoldChange>=1)
write.csv(up_gene,paste(sample.path,"up_gene.csv",sep=''),row.names = F)
#####down_gene
down_gene=subset(diff_gene,log2FoldChange<=-1)
write.csv(down_gene,paste(sample.path,'down_gene.csv',sep=''),row.names = F)
}