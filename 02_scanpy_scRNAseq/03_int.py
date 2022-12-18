import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as pt
# Update on Mar.6, 2021
import sys
import bbknn
from matplotlib import rcParams

OC_Artery = sc.read("/05_Result/05_exercise_mice/01_qc/OC-Artery/qc_after.h5ad")
OE_Artery = sc.read("/05_Result/05_exercise_mice/01_qc/OE-Artery/qc_after.h5ad")
YC_Artery = sc.read("/05_Result/05_exercise_mice/01_qc/YC-Artery/qc_after.h5ad")
YE_Artery = sc.read("/05_Result/05_exercise_mice/01_qc/YE-Artery/qc_after.h5ad")

OC_Blood = sc.read("/05_Result/05_exercise_mice/01_qc/OC-Blood/qc_after.h5ad")
OE_Blood = sc.read("/05_Result/05_exercise_mice/01_qc/OE-Blood/qc_after.h5ad")
YC_Blood = sc.read("/05_Result/05_exercise_mice/01_qc/YC-Blood/qc_after.h5ad")
YE_Blood = sc.read("/05_Result/05_exercise_mice/01_qc/YE-Blood/qc_after.h5ad")

OC_BM = sc.read("/05_Result/05_exercise_mice/01_qc/OC-BM/qc_after.h5ad")
OE_BM = sc.read("/05_Result/05_exercise_mice/01_qc/OE-BM/qc_after.h5ad")
YC_BM = sc.read("/05_Result/05_exercise_mice/01_qc/YC-BM/qc_after.h5ad")
YE_BM = sc.read("/05_Result/05_exercise_mice/01_qc/YE-BM/qc_after.h5ad")

OC_Intestine = sc.read("/05_Result/05_exercise_mice/01_qc/OC-Intestine/qc_after.h5ad")
OE_Intestine = sc.read("/05_Result/05_exercise_mice/01_qc/OE-Intestine/qc_after.h5ad")
YC_Intestine = sc.read("/05_Result/05_exercise_mice/01_qc/YC-Intestine/qc_after.h5ad")
YE_Intestine = sc.read("/05_Result/05_exercise_mice/01_qc/YE-Intestine/qc_after.h5ad")

OC_Kidney = sc.read("/05_Result/05_exercise_mice/01_qc/OC-Kidney/qc_after.h5ad")
OE_Kidney = sc.read("/05_Result/05_exercise_mice/01_qc/OE-Kidney/qc_after.h5ad")
YC_Kidney = sc.read("/05_Result/05_exercise_mice/01_qc/YC-Kidney/qc_after.h5ad")
YE_Kidney = sc.read("/05_Result/05_exercise_mice/01_qc/YE-Kidney/qc_after.h5ad")

OC_Liver = sc.read("/05_Result/05_exercise_mice/01_qc/OC-Liver/qc_after.h5ad")
OE_Liver = sc.read("/05_Result/05_exercise_mice/01_qc/OE-Liver/qc_after.h5ad")
YC_Liver = sc.read("/05_Result/05_exercise_mice/01_qc/YC-Liver/qc_after.h5ad")
YE_Liver = sc.read("/05_Result/05_exercise_mice/01_qc/YE-Liver/qc_after.h5ad")

OC_Lung = sc.read("/05_Result/05_exercise_mice/01_qc/OC-Lung/qc_after.h5ad")
OE_Lung = sc.read("/05_Result/05_exercise_mice/01_qc/OE-Lung/qc_after.h5ad")
YC_Lung = sc.read("/05_Result/05_exercise_mice/01_qc/YC-Lung/qc_after.h5ad")
YE_Lung = sc.read("/05_Result/05_exercise_mice/01_qc/YE-Lung/qc_after.h5ad")

OC_Spleen = sc.read("/05_Result/05_exercise_mice/01_qc/OC-Spleen/qc_after.h5ad")
OE_Spleen = sc.read("/05_Result/05_exercise_mice/01_qc/OE-Spleen/qc_after.h5ad")
YC_Spleen = sc.read("/05_Result/05_exercise_mice/01_qc/YC-Spleen/qc_after.h5ad")
YE_Spleen = sc.read("/05_Result/05_exercise_mice/01_qc/YE-Spleen/qc_after.h5ad")

OC_Testis = sc.read("/05_Result/05_exercise_mice/01_qc/OC-Testis/qc_after.h5ad")
OE_Testis = sc.read("/05_Result/05_exercise_mice/01_qc/OE-Testis/qc_after.h5ad")
YC_Testis = sc.read("/05_Result/05_exercise_mice/01_qc/YC-Testis/qc_after.h5ad")
YE_Testis = sc.read("/05_Result/05_exercise_mice/01_qc/YE-Testis/qc_after.h5ad")

Atlas = OC_Artery.concatenate(OE_Artery, YC_Artery, YE_Artery, OC_Blood, OE_Blood, YC_Blood, YE_Blood, OC_BM, OE_BM, YC_BM, YE_BM, OC_Intestine, OE_Intestine, YC_Intestine, YE_Intestine, OC_Kidney, OE_Kidney, YC_Kidney, YE_Kidney, OC_Liver, OE_Liver, YC_Liver, YE_Liver, OC_Lung, OE_Lung, YC_Lung, YE_Lung, OC_Spleen, OE_Spleen, YC_Spleen, YE_Spleen, OC_Testis, OE_Testis, YC_Testis, YE_Testis, join = 'outer', batch_categories=['OC_Artery', 'OE_Artery', 'YC_Artery', 'YE_Artery', 'OC_Blood', 'OE_Blood', 'YC_Blood', 'YE_Blood', 'OC_BM', 'OE_BM', 'YC_BM', 'YE_BM', 'OC_Intestine', 'OE_Intestine', 'YC_Intestine', 'YE_Intestine', 'OC_Kidney', 'OE_Kidney', 'YC_Kidney', 'YE_Kidney', 'OC_Liver', 'OE_Liver', 'YC_Liver', 'YE_Liver', 'OC_Lung', 'OE_Lung', 'YC_Lung', 'YE_Lung', 'OC_Spleen', 'OE_Spleen', 'YC_Spleen', 'YE_Spleen', 'OC_Testis', 'OE_Testis', 'YC_Testis', 'YE_Testis'])

results_file = '/05_Result/05_exercise_mice/02_int/int_only.h5ad'
Atlas.write(results_file)

sc.pp.normalize_per_cell(Atlas, counts_per_cell_after=1e4)

Atlas = sc.pp.filter_genes_dispersion(Atlas, subset = False, min_disp=.5, max_disp=None, 
                              min_mean=.0125, max_mean=10, n_bins=20, n_top_genes=None, 
                              log=True, copy=True)
                              
sc.pp.log1p(Atlas)
results_file2 = '/05_Result/05_exercise_mice/02_int/int.h5ad'
Atlas.write(results_file2)

Atlas.raw = Atlas
sc.pp.highly_variable_genes(Atlas, min_mean=0.0125, max_mean=3, min_disp=0.5,inplace=False)
sc.pp.regress_out(Atlas, ['n_counts', 'pct_counts_MT','n_genes_by_counts'])
sc.pp.scale(Atlas, max_value=10)	

results_file3 = '/05_Result/05_exercise_mice/02_int/scale.h5ad'
Atlas.write(results_file3)

sc.tl.pca(Atlas,use_highly_variable=True,n_comps=100)

sc.pl.pca_variance_ratio(Atlas, log=True,n_pcs=100)
pt.savefig('/05_Result/05_exercise_mice/02_int/pca.pdf')

sc.external.pp.bbknn(Atlas, batch_key='batch',n_pcs=100)

sc.tl.tsne(Atlas)
sc.tl.diffmap(Atlas,n_comps=100)
sc.pp.neighbors(Atlas, n_neighbors=10, use_rep='X_diffmap',n_pcs=100)

sc.tl.leiden(Atlas,resolution = 2.2)
sc.pl.tsne(Atlas, color=['leiden'],size=0.5, legend_loc='on data', legend_fontsize=5])
pt.savefig('/05_Result/05_exercise_mice/02_int/leiden_umap.pdf', bbox_inches='tight')

sc.pl.tsne(Atlas, color=['leiden'],size=0.5, legend_fontsize=5)
pt.savefig('/05_Result/05_exercise_mice/02_int/leiden_umap2.pdf', bbox_inches='tight')
sc.tl.rank_genes_groups(Atlas, 'leiden', method='wilcoxon')

results_file5 = '/05_Result/05_exercise_mice/02_int/Atlas_cluster.h5ad'
Atlas.write(results_file5)
