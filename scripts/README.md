# Script list and description:  
  
## About file extensions  
1. \*.sh files: This is an unix bash scripts.  
2. \*.sub files: This is a [Slurm Workload Manager](https://slurm.schedmd.com) script. Used to submit job on hpc-class.  
3. \*.R files: This is a R script.  
4. \*.pl files: This is a Perl script.  

## Script description

`fetch_fastq_names.sh`: fetch all the SRR accesions associated with a NCBI Bioproject from the SRA. You need to provide the Bioproject accession number and the output file name. Borrowed from [ISU GIF site](https://gif.biotech.iastate.edu/downloading-all-sra-files-related-bioprojectstudy). You need to install the [`edirect` software from NCBI](https://www.ncbi.nlm.nih.gov/books/NBK179288/) in order for this script to work.

`getSRR.sub`: 

build_bowtie2DB.sub
DESeq2.R

run_fastqc.sub
run_tophat.sub
transcript_index.sub
deinterleave_fastq.sh
download_SRR.sh
rename_chr.sh
run_tophat_rev.sub
splitsrrfiles.pl