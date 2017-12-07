# Script list and description:  
  
## About file extensions  
1. \*.sh files: This is an unix bash scripts.  
2. \*.sub files: This is a [Slurm Workload Manager](https://slurm.schedmd.com) script. Used to submit job on hpc-class.  
3. \*.R files: This is a R script.  
4. \*.pl files: This is a Perl script.  

## Script description

`fetch_fastq_names.sh`: fetch all the SRR accesions associated with a NCBI Bioproject from the SRA. You need to provide the Bioproject accession number and the output file name. Code borrowed from [ISU GIF site](https://gif.biotech.iastate.edu/downloading-all-sra-files-related-bioprojectstudy). You need to install the [`edirect` software from NCBI](https://www.ncbi.nlm.nih.gov/books/NBK179288/) in order for this script to work.

`getSRR.sub`: Download the requested SRR accesions in compressed fastq format (fastq.gz) and split them into forward and reverse reads (provided that they are paired-end reads). The list of accessions to be downloaded is provided in the `SRR.numbers` file. The job is going to ask for 1 node and 16 threads and is going to run for 20 hours or until complete, whatever happens first.  
  
`deinterleave_fastq.sh`: This script split forward and reverse reads in an interleaved paired-end fastq file. Borrowed from [nathanhaigh](https://gist.github.com/nathanhaigh/3521724).  
  
`splitsrrfiles.pl`: This script split forward and reverse reads in a concatenated paired-end fastq file (such as those in the NCBI SRA). Borrowed from [David Langenberger](https://www.biostars.org/p/19446/).  
  
`run_fastqc.sub`: Perform quality check on all the fastq.gz format file in the `input_files` folder, reading 16 files at the same time. The job is going to ask for 1 node and 16 threads and is going to run for 4 hours or until complete, whatever happens first.  
  
`build_bowtie2DB.sub`: Build bowtie2 index for the provded reference genome fasta file(s). The job is going to ask for 1 node and 16 threads and is going to run for 2 hours or until complete, whatever happens first.  
  
`rename_chr.sh`: This script changes the fasta header for each fasta file and replace it with the chromosome number or either "Mt", "Pt" or "UNKNOWN" for those chromosomes. This information is taken from the original fasta header. It is specific for maizeGDB reference V2 genomic fasta files.  
  
`transcript_index.sub`: Build tophat index for the provided reference transcriptome gff file. The job is going to ask for 1 node and 16 threads and is going to run for 4 hours or until complete, whatever happens first.  
  
`run_tophat.sub`: Runs tophat aligner individually for each SRR accession. The job is going to ask for 5 nodes and 16 threads per node and is going to run for 20 hours or until complete, whatever happens first.  
  
`run_tophat_rev.sub`: same as `run_tophat.sub` but starting from the last SRR accession, in order to run both script in parallel. The job is going to ask for 5 nodes and 16 threads per node and is going to run for 20 hours or until complete, whatever happens first. 

`DESeq2.R`: Script to format the obtained raw counts table and get the differentially expressed genes using DESeq2.