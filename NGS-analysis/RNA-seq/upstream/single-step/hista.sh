#!/bin/bash
for i in {73..84}
do
    hisat2 -x IRGSP_genome/IRGSP_genome -p 5 -1 SRR109915${i}_1_val_1.fq.gz -2 SRR109915${i}_2_val_2.fq.gz -S  SRR109915${i}.sam
done
