#!/bin/bash

## http://www.bio-info-trainee.com/2218.html
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50177
conda activate rna3 # 使用python3 环境
# new stepup project: GSE50177
wkd='/home/fleame/project/'
projectName='GSE50177'
cd $wkd
mkdir $projectName && cd $projectName
 
mkdir sra fastq reports clean aligned 
 
############# sra download #####
cd sra
# download the SRR_Acc_List.txt ## 
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA217298&o=acc_s%3Aa
# download the data 
cat SRR_Acc_List.txt | while read id; do (prefetch ${id});done
ps -ef | grep prefetch | awk '{print $2}' | while read id; do kill ${id}; done ## kill download processes
 
# 下载好之后的每个sra放在单独一个文件夹，可以自行手动处理，一个个copy到sra。这里尝试用shell脚本批处理。过程：从SRR_Acc_List.txt读取，移动到sra
cat SRR_Acc_List.txt | while read id; do mv -f $id/$id.sra  ./; done 
cat SRR_Acc_List.txt | while read id; do rm -rf $id; done #分两步走，检出全部移动完毕再删除文件夹
 
############# sra to fastq #####
for i in *sra
do
	fasterq-dump  --split-3 $i -f ../fastq 
	#fastq-dump --split-3 --skip-technical --clip --gzip $i ../fastq # python2环境老是报错，切换到python3的环境OK！
done
 
############# fastqc multiqc #####
cd ../fastq
ls *fastq | xargs fastqc -t 2 -O ../reports
cd ../reports && multiqc ./
    
############# data trim and clean #####
# paired
cd ../clean
#ls ../fastq/*_1.fastq.gz >1
#ls ../fastq/*_2.fastq.gz >2
#paste 1 2 > config
 
bin_trim_galore=trim_galore
#cat config |while read id
#do
  #arr=(${id})
  #fq1=${arr[0]}
  #fq2=${arr[1]} 
  #$bin_trim_galore -q 25 --phred33 --length 36 --stringency 3 --paired -o ./  $fq1 $fq2 
#done 
 
# single sequence
ls ../fastq/*.fastq|while read id;do ($bin_trim_galore -q 25 --phred33 --length 36 --stringency 3  -o ./  $id);done
 
############# align #####
# index has been build
# hisat2-build -p 4 '/home/fleame/public/references/genome/GRCh38.p13.genome.fa' genome
#ls *gz|cut -d"_" -f 1 | sort -u| while read id;
#do
  #hisat2 -p 10 -x /home/fleame/public/references/index/hisat/grch38/genome -1 ${id}_1_val_1.fq.gz   -2 ${id}_2_val_2.fq.gz  -S ../aligned/${id}.hisat.sam
#done
 
#unpaired
ls *fq|cut -d"." -f 1 | sort -u| while read id; 
do 
hisat2 -p 5 -x /home/fleame/public/references/index/hisat/grch38/genome -U ${id}.sra_trimmed.fq -S ../aligned/${id}.hisat.sam
done
 
############# sam to bam #####
cd ../aligned
ls *.sam| while read id ;do (samtools sort -O bam -@ 5 -o $(basename ${id} ".sam").bam ${id});done
rm *.sam
# 为bam文件建立索引
ls *.bam |xargs -i samtools index {}
# 比对结果统计
ls *.bam |while read id ;do ( samtools flagstat -@ 1 $id >  $(basename ${id} ".bam").flagstat  );done
 
############# featureCounts #####
gtf='/home/fleame/public/references/gtf/gencode.v32.annotation.gtf.gz'   
featureCounts -T 5 -p -t exon -g gene_id  -a $gtf -o  all.id.txt  *.bam  1>counts.id.log 2>&1 &
# 这样得到的  all.id.txt  文件就是表达矩阵，这个 featureCounts有非常多的参数可以调整。
 
#ls *.bam |while read id;do (nohup samtools view $id | htseq-count -f sam -s no -i gene_name - $gtf 1>${id%%.*}.geneCounts 2>${id%%.*}.HTseq.log&);done