#!/bin/bash

tar=/05_Result/05_rnaseq_mouse
raw=/01_rawdata/04_Exercise_mice/04_bulkRNA/00_new/Rawdata
map=$tar/03_map

for sample in `ls $raw`
do

echo -e '#!/bin/bash'>$map/scripts/${sample}_map.sh
echo "#PBS -N map_$sample">>$map/scripts/${sample}_map.sh
echo "#PBS -o $sample.out">>$map/scripts/${sample}_map.sh
echo "#PBS -e $sample.err">>$map/scripts/${sample}_map.sh
echo "#PBS -l nodes=1:ppn=12">>$map/scripts/${sample}_map.sh
echo "#PBS -q batch">>$map/scripts/${sample}_map.sh
echo '#PBS -p 1023'>>$map/scripts/${sample}_map.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'
echo Mapping:$sample is Starting
trim=$tar/02_trim
map=$tar/03_map

clean1=$trim/$sample/*_R1_val_1.fq.gz
clean2=$trim/$sample/*_R2_val_2.fq.gz

index=/03_Database/05_rna_seq/01_hisat2_index/mm10

result=$map/$sample
log=$map/logs
mkdir $result
hisat2=/anaconda3/bin/hisat2
$hisat2 -p 24 -x $index -1 $clean1 -2 $clean2 --dta -S $result/${sample}.sam 2>$log/${sample}_map.log

echo Mapping:$sample is Done

bam=$result/${sample}.bam
sort=$result/${sample}.sort.bam
sort1=$result/${sample}.sort.name.bam

samtools=/anaconda3/bin/samtools
$samtools view -bS -@ 24 -q 10 $result/${sample}.sam>$bam &&
$samtools sort -@ 24 $bam -o $sort &&
$samtools sort -@ 24 -n $bam -o $sort1 &&
$samtools index -@ 24 $sort


'>>$map/scripts/${sample}_map.sh

done


echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`
do
qsub $i &
done'>$map/run_map.sh

