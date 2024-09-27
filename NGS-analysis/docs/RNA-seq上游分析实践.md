# RNA Seq 上游分析实践

之前那一篇文章主要讲的是一些知识与工具的用法，这次用六组数据进行分析，得到基因表达矩阵。

<!--more-->

还是先从 NCBI 搜索数据，查看数据是双端测序的数据之后，多选条目，直接 Accession list

![](https://images.yuanj.top/2023090820339.png)

我使用的是下列数据

```txt
SRR25907783
SRR25907784
SRR25907785
SRR25907786
SRR25907787
SRR25907788
```

使用 sra-tools 批量下载数据，并且使用 nohup 把程序挂在后台下载

```bash
nohup prefetch --option-file acc_list.txt 
```

再使用脚本进行拆分

```bash
#!/bin/bash
mkdir SRR
mv ./SRR*/*.sra ./SRR
cd SRR
nohup fastq-dump --gzip --split-3 SRR*.sra 
```

这里千万不要在程序没有运行完成的时候就进行下一步操作，使用 top 可以查看后台进程

![](https://images.yuanj.top/20230908204934.png)

当程序运行完成后便可看到后台任务中已经没有前面运行的程序了

![](https://images.yuanj.top/20230908205137.png)

同一目录下的 nohup.out 文件中是后台进程的运行记录

拆分完成后就需要进行质量检测，使用通配符批量检测，并且将检测报告放在单独一个文件夹以便后面进行压缩

```bash
mkdir fastqc_report
fastqc SRR2590778*.fastq.gz -o ./fastqc_report
```

完成之后，将检测报告进行压缩，以便下载查看

```bash
zip -r fastqc_report.zip fastqc_report/*.html
```

依据检测报告对序列进行过滤，参数之前已经讲过，这里数据比较多，写一个 bash 脚本的 for 循环

```bash
#!/bin/bash
for i in {3..8}
do
    trimmomatic PE -threads 1 -phred33 SRR2590778${i}_1.fastq.gz SRR2590778${i}_2.fastq.gz -summary oryza_sativa_${i}.summary -baseout SRR2590778${i}.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 HEADCROP:15 MINLEN:36
done
```

下载水稻的参考基因组和注释文件进行 hisat2 比对

```bash
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/gff3/oryza_sativa/Oryza_sativa.IRGSP-1.0.57.gff3.gz
```

把下载的文件解压后重命名，方便使用

```bash
gzip -d Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz
mv Oryza_sativa.IRGSP-1.0.dna.toplevel.fa oryza_sativa.fa

gzip -d Oryza_sativa.IRGSP-1.0.57.gff3.gz
mv Oryza_sativa.IRGSP-1.0.57.gff3 oryza_sativa.gff3
```

先构建 hisat2 索引

```bash
hisat2-build oryza_sativa.fa oryza_sativa
```

构建完成之后开始将 reads 比对到参考基因组，还是写成脚本进行

```bash
#!/bin/bash
for i in {3..8}
do
    hisat2 -x hisat2_index/oryza_sativa -p 5 -1 SRR2590778${i}_1P.fastq.gz -2 SRR2590778${i}_2P.fastq.gz -S oryza_sativa_${i}.sam
done
```

还是写成 bash 脚本一键进行，这里我将每次排序的结构按照 1-6 分别命名

```bash
#!/bin/bash
for i in {3..8}
do
    samtools sort -n -@ 5 oryza_sativa_${i}.sam -o oryza_sativa_${i}
done
```

最后一步 featureCounts 生成基因计数表

```bash
#!/bin/bash
bam_files=(*.bam)

if [ ${#bam_files[@]} -gt 0 ]; then
    featureCounts -T 5 -t exon -g Name -a Lolium_perenne.gff3 -o  gene.counts -p "${bam_files[@]}"
else
    echo "No BAM files found in the current directory."
fi
```

至此完成，拿到矩阵

个人习惯完成分析后把文件进行一下归类

![](https://images.yuanj.top/20230908205327.png)

最后将数据下载到本地进行下游分析

![](https://images.yuanj.top/2023090820577.png)

全流程分析代码可以在我的 [GitHub](https://github.com/imjiaoyuan/Omics) 获取。