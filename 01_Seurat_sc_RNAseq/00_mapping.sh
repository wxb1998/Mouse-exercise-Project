##PBS -q long
#PBS -l nodes=node16:ppn=25
#count=/04_Raw_data/05_exercise_mice/buce/
count=/01_rawdata/01_mouse_sn
cd $count
 
for sample in YC-Testis;
do
cellranger count \
--id=$sample \
--transcriptome=/03_Database/01_cellranger_4/refdata-gex-mm10-2020-A \
--fastqs=/01_rawdata/04_Exercise_mice/02_jiace/$sample \
--sample=${sample}-1,${sample}-2,${sample}-3 \
--localcores=25 \
--localmem=64 \
--mempercore=150 \
--chemistry=SC3Pv3
rm $count/$sample/outs/*.bam
done


###just for sn
#PBS -q long
#PBS -l nodes=1:ppn=28
#PBS -N mice-brain_cb
/data/home/quj_lab/wangxuebao/miniconda3/envs/cellbender/bin/cellbender remove-background \
--input /dellstorage02/quj_lab/pingjiale/01_rawdata/02_mouse_sn/OC-brain/outs/raw_feature_bc_matrix.h5 \
--output /data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/02_cellbender/OC_brain_remove_background_raw_feature_bc_matrix.h5 \
--expected-cells 8209 \
--total-droplets-included 10453 \
--fpr 0.01 \
--epochs 100
