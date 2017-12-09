#!/bin/bash
sample_info=sample_summary.txt
sample=($(cut -f 1 "$sample_info"))

for fastq in ${sample[@]}
do
	echo -e "gene.id\t$fastq" | cat - ${fastq}_counts.txt > /tmp/out && mv /tmp/out ${fastq}_counts.txt
	join -1 1 -2 1 SRR1238715_counts.txt ${fastq}_counts.txt > /tmp/out && mv /tmp/out SRR1238715_counts.txt
done
cp SRR1238715_counts.txt counts.txt
cp SRR1238715_counts.bak SRR1238715_counts.txt
cut -f2 -d' ' --complement counts.txt > /tmp/out && mv /tmp/out counts.txt
