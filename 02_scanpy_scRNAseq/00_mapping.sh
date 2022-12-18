##PBS -q long
#PBS -l nodes=node16:ppn=25
#count=/04_Raw_data/05_exercise_mice/buce/

###Take an organ as an example###
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
