import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

samples = open("/data/home/quj_lab/wangxuebao/03_script/03_mouse_exerise/11_scanpy/samples.txt","r")
lines = samples.readlines()
samples.close()

for i in range(len(lines)):
        sample = lines[i].strip("\n")
        input_dir = '/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/03_qc/000_scanpy/001_scrublet/' + sample + '/'
        counts_matrix = scipy.io.mmread(input_dir + sample +'_matrix.mtx').T.tocsc()
        genes = np.array(scr.load_genes(input_dir + sample + '_genes.tsv', delimiter='\t', column=1))
        out_df = pd.read_csv(input_dir + sample + '_barcodes.tsv', header = None, index_col=None, names=['barcode'])

        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.069)

        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=85,
                                                          n_prin_comps=30)
        scrub.call_doublets(threshold=0.28)
        scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

        print (scrub.detected_doublet_rate_)
        out_df['doublet_rate'] = scrub.detected_doublet_rate_
        out_df['doublet_scores'] = doublet_scores
        out_df['predicted_doublets'] = predicted_doublets
        out_df.to_csv(input_dir + sample + '_doublet.txt', index=False,header=True,sep = "\t")
