#!/bin/bash

sample_info=sample_summary.txt
sample=($(cut -f 1 "$sample_info" | tail -n +2))
for fastq in ${sample[@]}
do
echo -e "gene.id\t$fastq" | cat - ${fastq}_counts.txt > /tmp/out && mv /tmp/out ${fastq}_counts.txt
join -1 1 -2 1 SRR1238715_counts.txt ${fast1}_counts.txt > /tmp/out && mv /tmp/out SRR1238715_counts.txt
done