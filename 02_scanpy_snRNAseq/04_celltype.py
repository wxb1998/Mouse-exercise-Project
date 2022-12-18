import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as pt
# Update on Mar.11, 2021
import sys
import bbknn
from matplotlib import rcParams

results_file = '/data/home/quj_lab/wangxuebao/00_tmp/sn_pc50_Atlas_cluster.h5ad'
Atlas = sc.read(results_file)

Atlas.obs['group']=Atlas.obs['batch']

old_to_new={'OC_spinal_cord':'OC','OE_spinal_cord':'OE', 'YC_spinal_cord':'YC', 'YE_spinal_cord':'YE', 'OC_cerebellum':'OC', 'OE_cerebellum':'OE', 'YC_cerebellum':'YC', 'YE_cerebellum':'YE','OC_brain':'OC','OE_brain':'OE', 'YC_brain':'YC', 'YE_brain':'YE', 'OC_heart':'OC', 'OE_heart':'OE', 'YC_heart':'YC', 'YE_heart':'YE', 'OC_muscle':'OC', 'OE_muscle':'OE', 'YC_muscle':'YC', 'YE_muscle':'YE'}

Atlas.obs['group'] = (
        Atlas.obs['batch']
        .map(old_to_new)
        .astype('category'))

Atlas.obs['tissue']=Atlas.obs['batch']

old_to_new={'OC_spinal_cord':'spinal_cord','OE_spinal_cord':'spinal_cord', 'YC_spinal_cord':'spinal_cord', 'YE_spinal_cord':'spinal_cord', 'OC_cerebellum':'cerebellum', 'OE_cerebellum':'cerebellum', 'YC_cerebellum':'cerebellum', 'YE_cerebellum':'cerebellum','OC_brain':'brain','OE_brain':'brain', 'YC_brain':'brain', 'YE_brain':'brain', 'OC_heart':'heart', 'OE_heart':'heart', 'YC_heart':'heart', 'YE_heart':'heart', 'OC_muscle':'muscle', 'OE_muscle':'muscle', 'YC_muscle':'muscle', 'YE_muscle':'muscle'}

Atlas.obs['tissue'] = (
        Atlas.obs['batch']
        .map(old_to_new)
        .astype('category')
)


# merge clusters, annotation
Atlas.obs['leiden_anno']=Atlas.obs['leiden']

#01_Ast  02_Ast_Brain
#03_OPC 04_OL
#05_OL_Brain 06_EC
#07_EC_CNS 
#08_Fib_Muscle 09_Fib_Heart 10_Mic 11_Schwann 12_Men 13_Epe 14_Neuron_SC 15_Neuron_CNS 16_Neuron_Brain 17_Per_Muscle
#18_Per_Heart 19_Per_CNS 20_SMC 21_Granule 22_MLI 23_ExN_Brain 24_InN_Brain
#25_Car 26_Fast_IIB
#27_Fast_IIX 28_Fast_IIA 29_Adi 30_NMJ_post 31_NMJ_pre
#32_Tendon 33_MJ 34_Satellite_cell 35_Mac1
#36_Mac2
#37_TC 38_Purkinje

old_to_new = {"0":"05_OL_Brain","1":"08_Fib_Muscle","2":"09_Fib_Heart","3":"26_Fast_IIB","4":"06_EC","5":"04_OL","6":"26_Fast_IIB","7":"16_Neuron_Brain","8":"04_OL","9":"15_Neuron_CNS",
"10":"28_Fast_IIA","11":"21_Granule","12":"25_Car","13":"05_OL_Brain","14":"27_Fast_IIX","15":"26_Fast_IIB","16":"25_Car","17":"21_Granule","18":"36_Mac2","19":"21_Granule",
"20":"21_Granule","21":"21_Granule","22":"21_Granule","23":"04_OL","24":"26_Fast_IIB","25":"26_Fast_IIB","26":"21_Granule","27":"04_OL","28":"25_Car","29":"38_Purkinje",
"30":"17_Per_Muscle","31":"27_Fast_IIX","32":"04_OL","33":"01_Ast","34":"25_Car","35":"05_OL_Brain","36":"03_OPC","37":"21_Granule","38":"18_Per_Heart","39":"23_ExN_Brain",
"40":"33_MJ","41":"36_Mac2","42":"26_Fast_IIB","43":"25_Car","44":"04_OL","45":"06_EC","46":"21_Granule","47":"10_Mic","48":"01_Ast","49":"04_OL",
"50":"21_Granule","51":"10_Mic","52":"04_OL","53":"05_OL_Brain","54":"26_Fast_IIB","55":"26_Fast_IIB","56":"21_Granule","57":"02_Ast_Brain","58":"24_InN_Brain","59":"16_Neuron_Brain",
"60":"01_Ast","61":"15_Neuron_CNS","62":"25_Car","63":"25_Car","64":"35_Mac1","65":"12_Men","66":"14_Neuron_SC","67":"22_MLI","68":"10_Mic","69":"02_Ast_Brain",
"70":"32_Tendon","71":"08_Fib_Muscle","72":"29_Adi","73":"06_EC","74":"07_EC_CNS ","75":"21_Granule","76":"11_Schwann","77":"10_Mic","78":"14_Neuron_SC","79":"21_Granule",
"80":"23_ExN_Brain","81":"23_ExN_Brain","82":"25_Car","83":"20_SMC","84":"14_Neuron_SC","85":"12_Men","86":"06_EC","87":"31_NMJ_pre","88":"19_Per_CNS","89":"23_ExN_Brain",
"90":"04_OL","91":"23_ExN_Brain","92":"37_TC","93":"13_Epe","94":"34_Satellite_cell","95":"23_ExN_Brain","96":"14_Neuron_SC","97":"30_NMJ_post","98":"10_Mic","99":"27_Fast_IIX",
"100":"25_Car","101":"36_Mac2","102":"09_Fib_Heart","103":"34_Satellite_cell","104":"18_Per_Heart","105":"25_Car","106":"17_Per_Muscle","107":"01_Ast"}
 
Atlas.obs['leiden_anno'] = (
    Atlas.obs['leiden']
    .map(old_to_new)
    .astype('category')
)

sc.pl.tsne(Atlas, color='leiden_anno', legend_loc='on data',legend_fontsize=6, title='')
pt.savefig('/data/home/quj_lab/wangxuebao/00_tmp/0004_celltype/sn_pc50_cell_type_umap.pdf')
sc.pl.tsne(Atlas, color='leiden_anno', title='')
pt.savefig('/data/home/quj_lab/wangxuebao/00_tmp/0004_celltype/sn_pc50_cell_type_umap_fu.pdf', bbox_inches='tight')

marker_genes={'01_Ast':["Aqp4","Gja1","Slc1a2"],'03_OPC':['Ppfibp1','Pdgfra'],'04_OL':["Mobp",'Mog','Plp1'],'06_EC':['Flt1','Pecam1',"Tek"],
'08_Fib_Muscle':['Pdgfra',"Dcn"],'10_Mic':["Ctss","Apbb1ip","Ptprc"],'11_Schwann':["Pmp22","Prx"],'12_Men':['Slc6a13',"Col3a1","Dcn"],'13_Epe':['Cfap54','Dnah6'],
'14_Neuron_SC':["Syt1",'Syp', "Snap25"],'17_Per_Muscle':['Mcam','Rgs5',"Pdgfrb"],
'20_SMC':["Tagln","Myh11"],'21_Granule':['Gabra6','Camk4','Etv1'],'22_MLI':["Lypd6"],'23_ExN_Brain':["Slc17a7","Satb2"],'24_InN_Brain':["Gad1",'Gad2'],
'25_Car':['Tnnt2','Myh7'],'26_Fast_IIB':['Myh4'],'27_Fast_IIX':["Myh1"],'28_Fast_IIA':['Myh2'],'29_Adi':['Adipoq'],
'30_NMJ_post':["Chrne"],'31_NMJ_pre':["Grik2", "Cdh19"],'32_Tendon':["Thbs4","Mkx"],'33_MJ':["Col22a1"],'34_Satellite_cell':["Pax7","Calcr"],'35_Mac1':["Mrc1","Cd86"],'36_Mac2':["Cd163"],'37_TC':["Cd247"],'38_Purkinje':['Pcp4', "Klhl1"]}

sc.pl.dotplot(Atlas, marker_genes, groupby='leiden_anno')
pt.savefig('/data/home/quj_lab/wangxuebao/00_tmp/0004_celltype/sn_pc50_dotplot.pdf', bbox_inches='tight')
sc.pl.stacked_violin(Atlas, marker_genes, groupby='leiden_anno', rotation=90)
pt.savefig('/data/home/quj_lab/wangxuebao/00_tmp/0004_celltype/sn_pc50_violin.pdf', bbox_inches='tight')

results_file = '/data/home/quj_lab/wangxuebao/00_tmp/0004_celltype/sn_pc50_cell_type.h5ad'
Atlas.write(results_file)

Atlas.obs[['batch', 'leiden_anno','tissue','group']].to_csv('/data/home/quj_lab/wangxuebao/00_tmp/0004_celltype/sn_pc50_celltype_prop.csv')
