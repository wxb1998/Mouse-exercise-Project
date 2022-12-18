library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(SeuratData)
library(patchwork)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(ggpubr)
set.seed(2)

all <- read.csv("/08_DEG/jiale/Exercise_mouse_DEGs.csv")
all$DE <- as.factor(ifelse(all$avg_logFC >=0 ,'up','down'))

Tissues <- unique(all$Tissue)

for (m in Tissues) {
all_DEGs <- subset(all, Tissue == paste0(m))
all_DEGs$ident1 <- paste0(all_DEGs$Tissue,"=",all_DEGs$Celltype,"=",all_DEGs$Gene)

for (j in c("up","down")) {

Aging<- subset(all_DEGs, Type == "Aging_DEGs" & DE==paste0(j))

OE_DEGs <- subset(all_DEGs, Type == "O_EX_DEGs")
OE_DEGs <- subset(OE_DEGs, ident1 %in% Aging$ident1)
DEGs <- rbind(Aging,OE_DEGs )
DEGs <- dcast(DEGs , ident1~Type, value.var = 'avg_logFC')
DEGs[is.na(DEGs)] <- 0
head(DEGs)
density <- melt(DEGs,id = c("ident1"))

colnames(density) <- c("ident1","group","avg_log2FC")

#####绘制aging基因和年老exercise基因的密度图####
nc=subset(density,density$group %in% "Aging_DEGs")
tb1=subset(density,density$group %in% "O_EX_DEGs")

y_peak1 <- which.max(density(nc$avg_log2FC)$y)
ncmean <- density(nc$avg_log2FC)$x[y_peak1]
y_peak2 <- which.max(density(tb1$avg_log2FC)$y)
tb1mean=density(tb1$avg_log2FC)$x[y_peak2]
wilcox.test(nc$avg_log2FC,tb1$avg_log2FC)
print(wilcox.test(nc$avg_log2FC,tb1$avg_log2FC))

p<-ggplot(density, aes(x = avg_log2FC))+ 
  geom_density(aes(fill = group),color = NA, alpha=0.7)+ 
  scale_fill_manual(values=c( "grey80","#99CD86","#9c9ece")) + #"#9c9ece","#e61c1d" "#e61c1d","#9c9ece"
  geom_vline(xintercept = 0.25,linetype="dashed",col="#e61c1d")+
  geom_vline(xintercept = -0.25,linetype="dashed",col="#9c9ece")+
  labs(title=paste0("P < 2.2e-16")) +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
ggsave(paste0("D:/7-Zip/0002_mouse_sport/08_DEG/midutu/",m,"_密度图merge-aging_",j,".pdf"), plot = p, width = 4, height = 2.5)

}
}


###计算rescue比例
OE_rescue <- read.csv("D:/7-Zip/0002_mouse_sport/08_DEG/jiale/03_OE_rescue/CT_num.csv")
OE_rescue$ident1 <- paste0(OE_rescue$tissue,'=',OE_rescue$Type,'=',OE_rescue$DE)
table(OE_rescue$ident1)
OE_rescue = OE_rescue %>% group_by(ident1) %>% mutate(sample.sum = sum(value))##统计每个总和样本的细胞数
OE_rescue <- unique(OE_rescue[,c("tissue","Type","DE","sample.sum")])
tmp_ex <- dcast(OE_rescue, tissue+DE~Type, value.var = 'sample.sum')
tmp_ex$Aging1 <- tmp_ex$Aging + tmp_ex$Rev
tmp_ex$rescue_prop <- tmp_ex$Rev / tmp_ex$Aging1 * 100
write.csv(tmp_ex,"D:/7-Zip/0002_mouse_sport/08_DEG/jiale/03_OE_rescue/OE_rescue_prop.csv")
