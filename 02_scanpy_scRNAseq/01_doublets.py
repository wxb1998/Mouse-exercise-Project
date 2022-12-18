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

samples = open("/01_Script/05_exercise_mice/01_doublets/sample.txt","r")
lines = samples.readlines()
samples.close()

for i in range(len(lines)):
        sample = lines[i].strip("\n")
        input_dir = '/04_Raw_data/05_exercise_mice/mapping_jiace/'+sample+'/outs/filtered_feature_bc_matrix'

        counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
        genes = np.array(scr.load_genes(input_dir + '/features2.tsv', delimiter='\t', column=1))
        out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])

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
        out_df.to_csv(input_dir + '/doublet.txt', index=False,header=True,sep = "\t")
