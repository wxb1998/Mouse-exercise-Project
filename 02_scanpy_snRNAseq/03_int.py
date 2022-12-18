import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as pt
# Update on Mar.6, 2021
import sys
import bbknn
from matplotlib import rcParams



####读取数据
OC_spinal_cord = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/OC-SC/final/qc_after.h5ad")
OE_spinal_cord = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/OE-SC/final/qc_after.h5ad")
YC_spinal_cord = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/YC-SC/final/qc_after.h5ad")
YE_spinal_cord = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/YE-SC/final/qc_after.h5ad")

OC_cerebellum = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/OC-cerebellum/final/qc_after.h5ad")
OE_cerebellum = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/OE-cerebellum/final/qc_after.h5ad")
YC_cerebellum = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/YC-cerebellum/final/qc_after.h5ad")
YE_cerebellum = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/YE-cerebellum/final/qc_after.h5ad")

OC_brain = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/OC-brain/final/qc_after.h5ad")
OE_brain = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/OE-brain/final/qc_after.h5ad")
YC_brain = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/YC-brain/final/qc_after.h5ad")
YE_brain = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/YE-brain/final/qc_after.h5ad")

OC_heart = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/OC-heart/final/qc_after.h5ad")
OE_heart = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/OE-heart/final/qc_after.h5ad")
YC_heart = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/YC-heart/final/qc_after.h5ad")
YE_heart = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/YE-heart/final/qc_after.h5ad")

OC_muscle = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/OC-muscle/final/qc_after.h5ad")
OE_muscle = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/OE-muscle/final/qc_after.h5ad")
YC_muscle = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/YC-muscle/final/qc_after.h5ad")
YE_muscle = sc.read("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/YE-muscle/final/qc_after.h5ad")


Atlas =OC_spinal_cord.concatenate(OE_spinal_cord, YC_spinal_cord, YE_spinal_cord, OC_cerebellum, OE_cerebellum, YC_cerebellum, YE_cerebellum,OC_brain, OE_brain, YC_brain, YE_brain, OC_heart, OE_heart, YC_heart, YE_heart, OC_muscle, OE_muscle, YC_muscle, YE_muscle, join = 'outer', batch_categories=['OC_spinal_cord','OE_spinal_cord', 'YC_spinal_cord', 'YE_spinal_cord', 'OC_cerebellum', 'OE_cerebellum', 'YC_cerebellum', 'YE_cerebellum','OC_brain','OE_brain', 'YC_brain', 'YE_brain', 'OC_heart', 'OE_heart', 'YC_heart', 'YE_heart', 'OC_muscle', 'OE_muscle', 'YC_muscle', 'YE_muscle'])


sc.pp.normalize_per_cell(Atlas, counts_per_cell_after=1e4)

Atlas = sc.pp.filter_genes_dispersion(Atlas, subset = False, min_disp=.5, max_disp=None, 
                              min_mean=.0125, max_mean=10, n_bins=20, n_top_genes=None, 
                              log=True, copy=True)

sc.pp.log1p(Atlas)
results_file2 = '/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/04_int/000_scanpy/sn_int2_only.h5ad'
Atlas.write(results_file2)

Atlas.raw = Atlas
sc.pp.highly_variable_genes(Atlas, min_mean=0.0125, max_mean=3, min_disp=0.5,inplace=False)
sc.pp.regress_out(Atlas, ['n_counts', 'pct_counts_MT','n_genes_by_counts'])
sc.pp.scale(Atlas, max_value=10)	


results_file3 = '/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/04_int/000_scanpy/sn_scale.h5ad'
Atlas.write(results_file3)

sc.tl.pca(Atlas,use_highly_variable=True,n_comps=50)

sc.pl.pca_variance_ratio(Atlas, log=True,n_pcs=50)
pt.savefig('/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/04_int/000_scanpy/sn_pca.pdf')

results_file4 = '/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/04_int/000_scanpy/sn_Atlas_pca.h5ad'
Atlas.write(results_file4)
