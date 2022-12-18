#!/bin/bash

for tissue in Blood BM Intestine Kidney Liver Lung Spleen Testis;

do
cd /02_Result/01_mouse_ex/05_single_tissue/02_cell_co/$tissue/01_YC/
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --counts-data=gene_name

cd /02_Result/01_mouse_ex/05_single_tissue/02_cell_co/$tissue/01_YC/
cellphonedb plot dot_plot 
cellphonedb plot heatmap_plot cellphonedb_meta.txt

done
