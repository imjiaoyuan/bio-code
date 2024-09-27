#!/bin/bash

set -e
set -o pipefail

############# environmental preparation #############

conda activate ATAC-seq # 使用conda环境

mkdir project
wkd='./project/'
projectName='project'
cd $wkd
mkdir $projectName && cd $projectName

# SRR_Acc_List.txt文件放在project根目录便于调用
mkdir sra fastq_gz fastqc_reports clean sorted genome_index compared rdup sam_index
 
############# data preparation #############

nohup prefetch -O . $(<SRR_Acc_List.txt) &
wait # 等待prefetch完成

cat ./SRR_Acc_List.txt | while read id; do
    mv -f "${id}/${id}.sra" ./sra/
    rm -rf "${id}"
done

############# sra to fastq #############

nohup fastq-dump --gzip --split-3 ./sra/*.sra --outdir ./fastqgz &

############# FastQc quality control #############

ls ./fastq_gz/*.fastq.gz | xargs fastqc -t 2 -O ./fastqc_reports
multiqc ./fastqc_reports  # 合并质量检测报告

############# Fastp filtering sequence #############

cat ./SRR_Acc_List.txt | while read id; do 
    fastp \
    -i ./fastq_gz/${id}_1.fastq.gz \
    -o ./clean/${id}_1.fastq.gz \
    -I ./fastq_gz/${id}_2.fastq.gz \
    -O ./clean/${id}_2.fastq.gz \
    -f 16 \
    -t 2 \
    -L
done

############# Reply to reference genome #############

genome='../oryza_sativa.fa'   # 基因组文件路径（需提前下载解压）
species='oryza_sativa'  # 物种名称（与基因组文件名一致）
cd ./genome_index
cp ${genome} ./

bwa index -a bwtsw ${genome}
cd ..

cat ./SRR_Acc_List.txt | while read id;
do
    bwa mem -v 3 -t 4 ./genome_index/${species}.fa ./fastp/${id}_1.fastq.gz ./fastp/${id}_2.fastq.gz -o ./compared/${id}.sam
done

############# sam to bam #############

cat ./SRR_Acc_List.txt | while read id ;do
    samtools \
    sort \
    -n -@ 5 \
    ./compared/${id}.sam \
    -o ./sorted/${id}.bam
done

############# Remove PCR duplicates #############

cat ./SRR_Acc_List.txt | while read id;
do
    sambamba markdup -r \
    -t 4 \
    ./sorted/${id}.bam \
    ./rdup/${id}_rdup.bam
done

############# Create sam index #############

cat ./SRR_Acc_List.txt | while read id;
do
    samtools index \
    -@ 4 \
    ./rdup/${id}_rdup.bam
done

mv ./rdup/*.bai ./sam_index/

############# Create bw files #####

cat ./SRR_Acc_List.txt | while read id;
do
    BAMscale scale \
    -t 4 \
    --bam ./rdup/${id}_rdup.bam
done