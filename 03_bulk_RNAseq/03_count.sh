#!/bin/bash
tar=/05_Result/05_rnaseq_mouse
raw=/01_rawdata/04_Exercise_mice/04_bulkRNA/00_new/Rawdata
count=$tar/04_count

for sample in `ls $raw`
do

echo -e '#!/bin/bash'>$count/scripts/${sample}_count.sh
echo "#PBS -N count_$sample">>$count/scripts/${sample}_count.sh
echo "#PBS -o $sample.out">>$count/scripts/${sample}_count.sh
echo "#PBS -e $sample.err">>$count/scripts/${sample}_count.sh
echo "#PBS -l nodes=1:ppn=12">>$count/scripts/${sample}_count.sh
echo "#PBS -q batch">>$count/scripts/${sample}_count.sh
echo '#PBS -p 1023'>>$count/scripts/${sample}_count.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

gtf=/03_Database/05_rna_seq/Mus_musculus.GRCm38.91.chr_final.gtf

uni=$tar/03_map
srt_bam=$uni/$sample/${sample}.sort.bam
count=$tar/04_count
result=$count/$sample
log=$count/logs
mkdir $result

echo HTseq:$sample is Starting
htseq=/anaconda3/bin/htseq-count
$htseq -f bam -r name -s no -a 10 $srt_bam $gtf > $result/${sample}_all.txt 2>$log/${sample}_count.log

echo HTseq-count has been Done'>>$count/scripts/${sample}_count.sh

done


echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`
do
qsub $i &
done'>$count/run_count.sh

