#!/bin/bash
for file in chr*.fasta
do
  if [ "$file" == "chrMt.fasta" ]
  then
    sed 's/.*mito.*/>Mt/g' $file > $(basename $file .fasta).temp
  elif [ "$file" == "chrPt.fasta" ]
  then
    sed 's/.*chloro.*/>Pt/g' $file > $(basename $file .fasta).temp
  else sed 's/chromosome:AGPv2://g' $file | sed 's/:1:[1-9][0-9].*$//g' > $(basename $file .fasta).temp
  fi
mv *.temp $file
done
