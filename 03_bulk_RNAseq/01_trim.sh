#!/bin/bash

tar=/05_Result/05_rnaseq_mouse
raw=/01_rawdata/04_Exercise_mice/04_bulkRNA/00_new/Rawdata
trim=$tar/02_trim

for sample in `ls $raw`
do

echo -e '#!/bin/bash'>$trim/scripts/${sample}_trim.sh
echo "#PBS -N Trim_$sample">>$trim/scripts/${sample}_trim.sh
echo "#PBS -o $sample.out">>$trim/scripts/${sample}_trim.sh
echo "#PBS -e $sample.err">>$trim/scripts/${sample}_trim.sh
echo "#PBS -l nodes=1:ppn=28">>$trim/scripts/${sample}_trim.sh
echo "#PBS -q batch">>$trim/scripts/${sample}_trim.sh
echo '#PBS -p 1023'>>$trim/scripts/${sample}_trim.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Trimming:$sample is Starting

trim=$tar/02_trim
log=$trim/logs

fq1=$raw/$sample/*_R1.fq.gz
fq2=$raw/$sample/*_R2.fq.gz

result=$trim/$sample
mkdir $result

trim_galore=/anaconda3/envs/py2/bin/trim_galore
$trim_galore --fastqc --path_to_cutadapt /anaconda3/envs/py2/bin/cutadapt --stringency 3 --paired --output_dir $result $fq1 $fq2 2>>$log/${sample}.log

echo Trimming has been Done'>>$trim/scripts/${sample}_trim.sh

done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`
do
qsub $i &
done'>$trim/run_trim.sh
