for i in {73..84}
do
    trim_galore -q 25 --phred33 --stringency 3 --length 36  --paired SRR109915${i}_1.fastq.gz SRR109915${i}_2.fastq.gz --gzip -o ./clean/
done
