#!/bin/bash

for file in chr*.fasta # For every file that start with chr and finishes with .fasta do as follows:
do
  if [ "$file" == "chrMt.fasta" ] 
  then
    sed 's/.*mito.*/>Mt/g' $file > $(basename $file .fasta).temp 
    # If filename is chrMt.fasta replace the fasta header with ">Mt"
    # fasta header need to have the srting `mito` on it.
  elif [ "$file" == "chrPt.fasta" ] 
  then
    sed 's/.*chloro.*/>Pt/g' $file > $(basename $file .fasta).temp
    # If filename is chrPt.fasta replace the fasta header with ">Pt"
    # fasta header need to have the srting `chloro` on it.
  else sed 's/chromosome:AGPv2://g' $file | sed 's/:1:[1-9][0-9].*$//g' > $(basename $file .fasta).temp
  	# else get rid of whatever is before and after the chromosome number on the fasta header.
  fi
mv *.temp $file
# replace the original files with the formated ones.
done
