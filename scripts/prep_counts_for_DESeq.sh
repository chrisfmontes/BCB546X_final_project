#!/bin/bash

sample_info=sample_summary.txt
# Specify the input sample file, where the SRR accession number is the first column
sample=($(cut -f 1 "$sample_info"))
# Create a bash array from the first column of $sample_info

for fastq in ${sample[@]}
# Go through each element of the $sample bash array
do
	echo -e "gene.id\t$fastq" | cat - ${fastq}_counts.txt > /tmp/out && mv /tmp/out ${fastq}_counts.txt
	# prints the 'gene.id' string on top of the first column and the corresponding SRR accession number
	# on top of the second column.

	join -1 1 -2 1 SRR1238715_counts.txt ${fastq}_counts.txt > /tmp/out && mv /tmp/out SRR1238715_counts.txt
	# Concatenates recursively the file SRR1238715_counts.txt with the each element of the array by column 1 as the common
	# data (gene.id). 

done

cp SRR1238715_counts.txt counts.txt
# Changes the name for SRR1238715_counts.txt to counts.txt

cp SRR1238715_counts.bak SRR1238715_counts.txt
# Restores SRR1238715_counts.txt as it was on the beginning of the script

cut -f2 -d' ' --complement counts.txt > /tmp/out && mv /tmp/out counts.txt
# Removes second column, which has duplicated data (because of how the script works)