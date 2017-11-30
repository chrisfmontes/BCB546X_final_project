#!/bin/sh
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4:00:00
#SBATCH --mail-user=chrisfm@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --error=SRR.%J.err
#SBATCH --output=SRR.%J.out

module load sratoolkit
#module load parallel

# parallel --jobs 3 "fastq-dump --split-files --origfmt --gzip {}" ::: SRR.numbers
cat SRR.numbers | xargs -P 3 fastq-dump --split-files --gzip