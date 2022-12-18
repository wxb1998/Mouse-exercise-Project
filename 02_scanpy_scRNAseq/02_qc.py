import numpy as np
import pandas as pd
import scanpy as sc
# Update on Apr.23, 2021
import sys
import matplotlib.pyplot as pt


samples = open("/01_Script/05_exercise_mice/02_yuchuli/sample10.txt","r")
lines = samples.readlines()
samples.close()

for i in range(len(lines)):
	sc.settings.verbosity = 3
	sc.logging.print_header()
	sample = lines[i].strip("\n")
	in_dir = "/04_Raw_data/05_exercise_mice/mapping_jiace/"+sample+"/outs/filtered_feature_bc_matrix/final/"
	adata = sc.read_10x_mtx(
		in_dir,
		var_names='gene_symbols',
		cache=True)
	adata.var_names_make_unique()

	sc.pp.filter_cells(adata, min_genes=200) 
	sc.pp.filter_genes(adata, min_cells=3)
	outdir = "/05_Result/05_exercise_mice/01_qc/"+sample+"/qc_before.h5ad"
	adata.write(outdir)
	
	adata.var['MT'] = adata.var_names.str.startswith('mt-')
	sc.pp.calculate_qc_metrics(adata, qc_vars=['MT'], percent_top=None, log1p=False, inplace=True)

	sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4, multi_panel=True)
	pt.savefig("/05_Result/05_exercise_mice/01_qc/"+sample+"/n_genes_by_counts.pdf")
	sc.pl.violin(adata, ['total_counts'], jitter=0.4, multi_panel=True)
	pt.savefig("/05_Result/05_exercise_mice/01_qc/"+sample+"/total_counts.pdf")
	sc.pl.violin(adata, ['pct_counts_MT'],jitter=0.4, multi_panel=True)
	pt.savefig("/05_Result/05_exercise_mice/01_qc/"+sample+"/pct_counts_MT.pdf")

	sc.pl.scatter(adata, x='total_counts', y='pct_counts_MT')
	pt.savefig("/05_Result/05_exercise_mice/01_qc/"+sample+"/scatter_pct_counts_MT.pdf")
	sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
	pt.savefig("/05_Result/05_exercise_mice/01_qc/"+sample+"/scatter_n_genes_by_counts.pdf")
	adata = adata[adata.obs.n_genes_by_counts > 500, :]
#The cells with more than 4,000 genes detected were excluded for peripheral blood, spleen and testis were excluded; The cells with more than 6,000 genes detected were excluded for bone marrow, small intestine, lung, aorta, liver, and kidney were excluded. 
	adata = adata[adata.obs.n_genes_by_counts < 6000, :]
#The cells with more than 10% mitochondrial gene ratio for were excluded; The cells with more than 20% mitochondrial gene ratio for aorta, liver and testis were excluded; The cells with more than 50% mitochondrial gene ratio for kidney were excluded. 
	adata = adata[adata.obs.pct_counts_MT < 10, :]

	sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4, multi_panel=True)
	pt.savefig("/05_Result/05_exercise_mice/01_qc/"+sample+"/qc_n_genes_by_counts.pdf")
	sc.pl.violin(adata, ['total_counts'], jitter=0.4, multi_panel=True)
	pt.savefig("/05_Result/05_exercise_mice/01_qc/"+sample+"/qc_total_counts.pdf")
	sc.pl.violin(adata, ['pct_counts_MT'],jitter=0.4, multi_panel=True)
	pt.savefig("/05_Result/05_exercise_mice/01_qc/"+sample+"/qc_pct_counts_MT.pdf")

	adata.write('/05_Result/05_exercise_mice/01_qc/'+sample+'/qc_after.h5ad')
