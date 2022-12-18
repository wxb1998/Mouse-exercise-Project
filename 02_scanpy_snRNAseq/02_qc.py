import numpy as np
import pandas as pd
import scanpy as sc
# Update on Apr.23, 2021
import sys
import matplotlib.pyplot as pt


samples = open("/data/home/quj_lab/wangxuebao/03_script/03_mouse_exerise/11_scanpy/02_yuchuli/sample2.5.txt","r")
lines = samples.readlines()
samples.close()

for i in range(len(lines)):
	sc.settings.verbosity = 3
	sc.logging.print_header()
	sample = lines[i].strip("\n")
	in_dir = "/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/"+ sample+'/final/'
	adata = sc.read_10x_mtx(
		in_dir,
		var_names='gene_symbols',
		cache=True)
	adata.var_names_make_unique()

	sc.pp.filter_cells(adata, min_genes=200) 
	sc.pp.filter_genes(adata, min_cells=3)
	outdir = "/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/"+ sample+'/final/'+"qc_before.h5ad"
	adata.write(outdir)
	
	adata.var['MT'] = adata.var_names.str.startswith('mt-')
	sc.pp.calculate_qc_metrics(adata, qc_vars=['MT'], percent_top=None, log1p=False, inplace=True)

	sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4, multi_panel=True)
	pt.savefig("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/"+ sample+'/final'+"/n_genes_by_counts.pdf")
	sc.pl.violin(adata, ['total_counts'], jitter=0.4, multi_panel=True)
	pt.savefig("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/"+ sample+'/final'+"/total_counts.pdf")
	sc.pl.violin(adata, ['pct_counts_MT'],jitter=0.4, multi_panel=True)
	pt.savefig("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/"+ sample+'/final'+"/pct_counts_MT.pdf")

	sc.pl.scatter(adata, x='total_counts', y='pct_counts_MT')
	pt.savefig("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/"+ sample+'/final'+"/scatter_pct_counts_MT.pdf")
	sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
	pt.savefig("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/"+ sample+'/final'+"/scatter_n_genes_by_counts.pdf")

	adata = adata[adata.obs.n_genes_by_counts > 200, :]
	adata = adata[adata.obs.n_genes_by_counts < 6000, :]
	adata = adata[adata.obs.pct_counts_MT < 2.5, :]

	sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4, multi_panel=True)
	pt.savefig("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/"+ sample+'/final'+"/qc_n_genes_by_counts.pdf")
	sc.pl.violin(adata, ['total_counts'], jitter=0.4, multi_panel=True)
	pt.savefig("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/"+ sample+'/final'+"/qc_total_counts.pdf")
	sc.pl.violin(adata, ['pct_counts_MT'],jitter=0.4, multi_panel=True)
	pt.savefig("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/"+ sample+'/final'+"/qc_pct_counts_MT.pdf")

	adata.write("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/"+ sample+'/final'+'/qc_after.h5ad')
