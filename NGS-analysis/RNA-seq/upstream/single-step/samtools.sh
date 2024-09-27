#!/bin/bash
for i in {73..84}
do
    samtools sort -n -@ 5  SRR109915${i}.sam -o SRR109915${i}.bam
done
