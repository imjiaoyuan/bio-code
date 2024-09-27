wkd='project'
SRRfileshead='SRR109915'

conda env create -f env.yml
conda activate genomics

cd ..
mkdir $wkd
cd $wkd
mkdir SRR fastqgz fastqc_reports clean sorted
mv SRR_Acc_List.txt SRR/

cd SRR
nohup prefetch -O . $(<SRR_Acc_List.txt) &
cat SRR_Acc_List.txt | while read id; do mv -f $id/$id.sra  ./; done

echo "原始数据下载完成" >> ../report.txt

for i in *sra
do
	fastq-dump --gzip --split-3 $i -f ../fastqgz
done

echo "sra数据解压完成" >> ../report.txt

cd ../fastq_gz
ls *.fastq.gz | xargs fastqc -t 2 -O ../fastqc_reports
cd ../fastqc_reports && multiqc ./

echo "质量检测完成" >> ../report.txt

#!/bin/bash
for i in {78..82}
do
    hisat2 -x Paeoniaostii/Paeoniaostii -p 5 -1 SRR179997${i}_1.fastq.gz -2 SRR179997${i}_2.fastq.gz -S Paeoniaostii_${i}.sam
done

echo "将reads比对到参考基因组完成" >> ../report.txt

for i in {78..82}
do
    sambamba sort -n -@ 5 Paeoniaostii_${i}.sam -o Paeoniaostii_${i}.bam
done

echo "sam文件排序转bam文件完成" >> ../report.txt

for i in {73..84}
do
    featureCounts -T 5 -t exon -g gene_id -a IRGSP-1.0_representative_transcript_exon_2024-01-11.gtf -o  gene.counts -p IRGSP_genome_57.bam IRGSP_genome_58.bam IRGSP_genome_59.bam IRGSP_genome_60.bam IRGSP_genome_61.bam IRGSP_genome_62.bam IRGSP_genome_63.bam IRGSP_genome_64.bam
done

featureCounts -T 5 -t exon -g gene_id -a IRGSP-1.0_representative_transcript_exon_2024-01-11.gtf -o  gene.counts -p IRGSP_genome_73 IRGSP_genome_74 IRGSP_genome_75 IRGSP_genome_76 IRGSP_genome_77 IRGSP_genome_78 IRGSP_genome_79 IRGSP_genome_80 IRGSP_genome_81 IRGSP_genome_82 IRGSP_genome_83 IRGSP_genome_84 
