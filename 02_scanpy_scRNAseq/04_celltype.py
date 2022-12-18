import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as pt
# Update on Mar.11, 2021
import sys
import bbknn
from matplotlib import rcParams

results_file = '/05_Result/05_exercise_mice/02_int/Atlas_cluster.h5ad'
Atlas = sc.read(results_file)

Atlas.obs['group']=Atlas.obs['batch']

old_to_new={'OC_Artery':'OC', 'OE_Artery':'OE', 'YC_Artery':'YC', 'YE_Artery':'YE', 'OC_Blood':'OC', 'OE_Blood':'OE', 'YC_Blood':'YC', 'YE_Blood':'YE', 'OC_BM':'OC', 'OE_BM':'OE', 'YC_BM':'YC', 'YE_BM':'YE', 'OC_Intestine':'OC', 'OE_Intestine':'OE', 'YC_Intestine':'YC', 'YE_Intestine':'YE', 'OC_Kidney':'OC', 'OE_Kidney':'OE', 'YC_Kidney':'YC', 'YE_Kidney':'YE', 'OC_Liver':'OC', 'OE_Liver':'OE', 'YC_Liver':'YC', 'YE_Liver':'YE', 'OC_Lung':'OC', 'OE_Lung':'OE', 'YC_Lung':'YC', 'YE_Lung':'YE', 'OC_Spleen':'OC', 'OE_Spleen':'OE', 'YC_Spleen':'YC', 'YE_Spleen':'YE', 'OC_Testis':'OC', 'OE_Testis':'OE', 'YC_Testis':'YC', 'YE_Testis':'YE'}

Atlas.obs['group'] = (
        Atlas.obs['batch']
        .map(old_to_new)
        .astype('category')
)

Atlas.obs['tissue']=Atlas.obs['batch']

old_to_new={'OC_Artery':'Artery', 'OE_Artery':'Artery', 'YC_Artery':'Artery', 'YE_Artery':'Artery', 'OC_Blood':'Blood', 'OE_Blood':'Blood', 'YC_Blood':'Blood', 'YE_Blood':'Blood', 'OC_BM':'BM', 'OE_BM':'BM', 'YC_BM':'BM', 'YE_BM':'BM', 'OC_Intestine':'Intestine', 'OE_Intestine':'Intestine', 'YC_Intestine':'Intestine', 'YE_Intestine':'Intestine', 'OC_Kidney':'Kidney', 'OE_Kidney':'Kidney', 'YC_Kidney':'Kidney', 'YE_Kidney':'Kidney', 'OC_Liver':'Liver', 'OE_Liver':'Liver', 'YC_Liver':'Liver', 'YE_Liver':'Liver', 'OC_Lung':'Lung','OE_Lung':'Lung', 'YC_Lung':'Lung','YE_Lung':'Lung', 'OC_Spleen':'Spleen', 'OE_Spleen':'Spleen', 'YC_Spleen':'Spleen', 'YE_Spleen':'Spleen', 'OC_Testis':'Testis', 'OE_Testis':'Testis', 'YC_Testis':'Testis', 'YE_Testis':'Testis'}

Atlas.obs['tissue'] = (
        Atlas.obs['batch']
        .map(old_to_new)
        .astype('category')
)

# merge clusters, annotation
Atlas.obs['leiden_anno']=Atlas.obs['leiden']

old_to_new = {"0":"62_BC","1":"53_CD8_Naive","2":"62_BC","3":"62_BC","4":"36_SMC","5":"62_BC","6":"53_CD8_Naive","7":"62_BC","8":"20_SPC","9":"62_BC",
"10":"38_Mono","11":"15_PT-S1","12":"12_LOH","13":"62_BC","14":"54_CD8_Mem","15":"62_BC","16":"27_EC_Liver","17":"33_Fib_Aorta","18":"20_SPC","19":"55_CD8_CTL",
"20":"43_Mac1","21":"22_RS","22":"62_BC","23":"15_PT-S1","24":"27_EC_Liver","25":"33_Fib_Aorta","26":"62_BC","27":"43_Mac1","28":"57_CD4-CD8-TC","29":"43_Mac1",
"30":"27_EC_Liver","31":"62_BC","32":"59_NKT","33":"29_EC_Aorta","34":"62_BC","35":"62_BC","36":"27_EC_Liver","37":"62_BC","38":"22_RS","39":"46_Neu",
"40":"52_CD4_Mem","41":"62_BC","42":"59_NKT","43":"51_CD4_Naive","44":"20_SPC","45":"15_PT-S1","46":"46_Neu","47":"29_EC_Aorta","48":"53_CD8_Naive","49":"52_CD4_Mem",
"50":"21_Early_RS","51":"52_CD4_Mem","52":"17_PT-S3","53":"38_Mono","54":"30_EC_Kidney","55":"53_CD8_Naive","56":"62_BC","57":"46_Neu","58":"53_CD8_Naive","59":"43_Mac1",
"60":"46_Neu","61":"46_Neu","62":"15_PT-S1","63":"63_Pla","64":"31_EC","65":"62_BC","66":"56_CD4+CD8+TC","67":"48_mDC","68":"22_RS","69":"46_Neu",
"70":"43_Mac1","71":"44_Mac2","72":"22_RS","73":"62_BC","74":"17_PT-S3","75":"62_BC","76":"53_CD8_Naive","77":"62_BC","78":"43_Mac1","79":"62_BC",
"80":"46_Neu","81":"46_Neu","82":"46_Neu","83":"15_PT-S1","84":"09_CD-PC","85":"59_NKT","86":"53_CD8_Naive","87":"20_SPC","88":"58_NK","89":"62_BC",
"90":"37_Mast","91":"33_Fib_Aorta","92":"28_EC_Lung","93":"43_Mac1","94":"27_EC_Liver","95":"20_SPC","96":"59_NKT","97":"15_PT-S1","98":"43_Mac1","99":"51_CD4_Naive",
"100":"23_ES","101":"63_Pla","102":"58_NK","103":"36_SMC","104":"33_Fib_Aorta","105":"16_PT-S2","106":"21_Early_RS","107":"23_ES","108":"33_Fib_Aorta","109":"36_SMC",
"110":"36_SMC","111":"23_ES","112":"36_SMC","113":"43_Mac1","114":"27_EC_Liver","115":"49_pDC","116":"61_LProBC","117":"44_Mac2","118":"43_Mac1","119":"63_Pla",
"120":"17_PT-S3","121":"12_LOH","122":"36_SMC","123":"14_DLOH","124":"46_Neu","125":"29_EC_Aorta","126":"06_Pericyte","127":"23_ES","128":"15_PT-S1","129":"34_Fib_Intestine",
"130":"37_Mast","131":"08_CD-IC","132":"41_Kup","133":"31_EC","134":"15_PT-S1","135":"42_AMac","136":"58_NK","137":"23_ES","138":"12_LOH","139":"06_Pericyte",
"140":"12_LOH","141":"43_Mac1","142":"47_Bas","143":"26_Epi_Intestine","144":"43_Mac1","145":"21_Early_RS","146":"29_EC_Aorta","147":"50_ProTC","148":"24_Sertoil","149":"33_Fib_Aorta",
"150":"59_NKT","151":"45_ProNeu","152":"41_Kup","153":"02_Cho","154":"24_Sertoil","155":"29_EC_Aorta","156":"62_BC","157":"39_Mono_BM","158":"12_LOH","159":"46_Neu",
"160":"43_Mac1","161":"43_Mac1","162":"60_ProBC","163":"23_ES","164":"38_Mono","165":"43_Mac1","166":"59_NKT","167":"01_Hep","168":"62_BC","169":"18_SPG",
"170":"07_Progenitor","171":"29_EC_Aorta","172":"05_Aorl_Neuron","173":"62_BC","174":"10_CD-Trans","175":"43_Mac1","176":"39_Mono_BM","177":"11_DCT","178":"32_Fib_Lung","179":"31_EC",
"180":"35_Fib_Testis","181":"33_Fib_Aorta","182":"19_Early_SPC","183":"04_AT2","184":"03_AT1","185":"40_Mega","186":"53_CD8_Naive","187":"32_Fib_Lung","188":"43_Mac1","189":"22_RS",
"190":"29_EC_Aorta","191":"59_NKT","192":"62_BC","193":"20_SPC","194":"34_Fib_Intestine","195":"25_Epi_Lung","196":"59_NKT","197":"13_Podocytes"}
 
Atlas.obs['leiden_anno'] = (
    Atlas.obs['leiden']
    .map(old_to_new)
    .astype('category')
)

sc.pl.tsne(Atlas, color='leiden_anno', legend_loc='on data',legend_fontsize=6, title='')
pt.savefig('/05_Result/05_exercise_mice/04_celltype/cell_type_umap.pdf')
sc.pl.tsne(Atlas, color='leiden_anno', title='')
pt.savefig('/05_Result/05_exercise_mice/04_celltype/cell_type_umap_fu.pdf', bbox_inches='tight')

marker_genes={'01_Hep':['Alb','Apoc3'],'02_Cho':['Krt18','Spp1','Sox9'],'03_AT1':['Hopx','Clic5','Cav1'],'04_AT2':['Sftpc', 'Sftpb', 'Sftpd'],'05_Aorl_Neuron':['Plp1','Kcna1'],'06_Pericyte':['Higd1b'],'07_Progenitor':['Sox9','Muc4'],
'08_CD-IC':['Atp6v0d2','Atp6v1g3','Slc26a4'],'09_CD-PC':['Aqp2','Hsd11b2'],'10_CD-Trans':['Insrr','Rhbg'],'11_DCT':['Pvalb','Slc12a3','Slc8a1','Calb1'],'12_LOH':['Umod','Slc12a1'],'13_Podocytes':['Nphs1','Nphs2'],'14_DLOH':['Aqp1'],'15_PT-S1':['Cubn','Lrp2','Slc5a2','Slc5a12'],'16_PT-S2_PT-S3':['Atp11a','Slc13a3'],
'18_SPG':['Uchl1','Dmrt1','Prdm9','Stra8'],'19_Early_SPC':['Sycp1', 'Sycp2', 'Sycp3'],'20_SPC':['Tbpl1','Mns1','Lyar'],'21_Early_RS':['Chd5','Calr3'],'22_RS':['Acrv1','Spag6','Tekt1','Txndc8'],'23_ES':['Prm2','Tnp2','Prm1','Tnp1'],'24_Sertoil':['Ctsl','Clu'],
'24_Epi':['Cdh1','Epcam'],'26_EC':['Pecam1','Esam','Ctla2a'],'31_Fib':['Lum','Dcn'],'36_SMC':['Tagln','Postn'],'IC':['Ptprc','Spi1'],
'37_Mast':['Kit'],'38_Mono':['Cd14','Fcgr3','Vcan'],'39_Mono_BM':['Ctsg','Prtn3'],'40_Mega':['Pf4','Ppbp'],'41_Kup':['Cd163','Clec4f','Adgre1'],'42_AMac':['Ear2'],'43_Mac1':['Cd68', 'Csf1r'],'44_Mac2':['Cd163'],'45_ProNeu':['Mki67'],'46_Neu':['S100a9','S100a8','Ngp','Cxcr2'],'47_Bas':['Ms4a2', 'Cpa3'],'48_mDC':['Irf8','Pld4','Clec9a','Fcer1a'],'49_pDC':['Siglech'],
'49_TC':['Cd3d', 'Cd3g','Cd4','Ccr7','Lef1','Prdm1','Pdcd1','Cd8a','Gzmk',"Gzmb"],'58_NK':['Nkg7','Prf1','Klrd1'],'60_ProBC':['Ezh2'],'61_LProBC':['Rag1'],'62_BC':['Cd79a', 'Cd19','Ms4a1'],'63_Pla':['Mzb1', 'Jchain']}

sc.pl.dotplot(Atlas, marker_genes, groupby='leiden_anno')
pt.savefig('/05_Result/05_exercise_mice/04_celltype/dotplot.pdf', bbox_inches='tight')
sc.pl.stacked_violin(Atlas, marker_genes, groupby='leiden_anno', rotation=90)
pt.savefig('/05_Result/05_exercise_mice/04_celltype/violin.pdf', bbox_inches='tight')

results_file = '/05_Result/05_exercise_mice/04_celltype/cell_type.h5ad'
Atlas.write(results_file)

Atlas.obs[['batch', 'leiden_anno','tissue','group']].to_csv('/05_Result/05_exercise_mice/04_celltype/celltype_prop.csv')
