**A Reproduction Of: Transposable Elements Contribute to Activation of Maize Genes in Response to Abiotic Stress**

**By: Maxwell McReynolds, Christian Montes-Serey, and Brian Zebosi**

**Original Paper Authors: Irina Makarevitch, Amanda J. Waters, Patrick T. West, Michelle Stitzer, Candice N. Hirsch, Jeffrey Ross-Ibarra, Nathan M. Springer**

**Paper Published: January 8, 2015, PLOS Genetics**

**Paper Introduction**

Transposable elements (TEs) have been shown to make up a large portion of many eukaryotic genomes, especially that of *Zea Mays*.  Commonly known for their deleterious effects on host DNA, it has been shown that TEs also can contribute to regulatory variation in gene expression via complex interactions with host DNA. Previous work in *Nicotiana tabacum, Arabidopsis thaliana, *and *Oryzae sativa* has shown that abiotic stresses activate the TEs and genes they insert near.  The rationale of their research was to fill the gap in knowledge as to what maize genes and nearby TEs are being activated in response to abiotic stress and what the proposed mechanism(s) of this process might be. 

Profiling of the gene and TE transcript levels was performed to analyze the effect of TEs on gene expression response to environmental stress using RNA-sequencing methods. The maize genotypes B74, Oh43, and Mo17 were the three maize inbred lines that were used for transciptomic profiling. There was a mix of genes up- and down-regulated in the three genotypes in response to the individually-applied abiotic stresses: cold (5°C for 16 hours), heat (50°C for 4 hours), high salt (300mM NaCl 20 hours prior to collection) or UV stress (2 hours). There was similar transcriptional response to each individual stress between the genotypes along with evidence of genotype-specific responses.  Analysis of TE families located within 1KB upstream of the transcription start site of the up-regulated abiotic stress responsive genes (identified by their RNA-seq experiment) revealed four to nine different TE families associated with induced expression in each stress condition. TEs were suggested to provide local enhancer activities that stimulate expression as revealed by analysis of stress-induced transcripts and TE proximity to response genes. This form of "TE expression enhancement" may drive expression of cryptic promoters in non-coding regions of the genome thus producing stress-responsive transcripts that may not produce functional proteins. They used whole-genome shotgun re-sequencing data from Mo17 and Oh43 to identify novel insertions of elements from the TE families they identified earlier.  For each of the genes they identified they found that the alleles that lack TE insertions did not exhibit stress responsive expression.  This suggests that insertion polymorphisms for the TE families they identified can generate novel expression responses for genes nearby. 

The research presented in this paper suggests that TEs are an important source of allelic regulatory variation in individual gene response to abiotic stresses in maize.  The exact mechanism to which TEs provide novel insertions and whether those inserts function as enhancer elements, novel cis-regulatory binding sites, or even sites of post-transcriptional modification still remains unclear.

**The A-MAIZE-ING Trio’s Objective**

For our reproduction we decided to focus primarily on the initial RNA-sequencing analysis outlined in the paper.  This analysis looks at whole transcriptome response to individual stresses of cold, heat, salt, and UV stress.  The authors used three agronomically important maize genotypes (Oh43, B73, and Mo17) each exposed to the stress conditions.  The RNA collected was from three biological replicates of each genotype exposed to heat and cold, and one biological replicate for the UV and salt stresses.  We aim to replicate this major portion of the paper and generate differential expression data for each genotype exposed to each individual stress.  We used two datasets for this project, one provided by the authors from the Expression Atlas from EMBL-EBI and one from the SRA database.  

**Technical Details **

1. Transcriptomic Data Set Preparation (not performed by the Trio)

RNA samples prepared into TruSeq libraries and sequenced on an Illumina HiSeq2000.  Reads were deposited on the NCBI Short Read Archive as SRA format files

2. Data Inspection

Following the "Data Availability" table we accessed the short read archive under study identity PRJNA244661 on the NCBI BioProject data collection.  Inspection of the page showed there was 85 gb of data with 33 paired end read files equaling 66 files in total. 

3. Obtaining Read Names

To obtain read names from the NCBI Bioproject accession we used fetch_fastq_names.sh unix bash script to fetch all of the SRR accessions from the NCBI BioProject.  We provided the accession number and the output file in our script while running the edirect software that’s available from NCBI. 

4. Downloading the Reads

Due to the large size of the files we developed a script to download SRR accessions in a compressed fastq format and split the paired end reads using get SRR.sub on a high performance computing (HPC-class) cluster.  We used 1 node and 16 threads for the download (appx. 20 hours in length).

5. Quality Checking the Reads 

Quality checks of the fastq.gz format files was performed using the run_fasttqc.subscript on HPC-class in the unix environment.  This job analyzed 16 files at a time and was performed using 1 node and 16 threads and took 4 hours to complete. 

6. Read Alignment

The splice junction mapper TopHat was used was used on the downloaded files and aligned the RNA-Seq reads to the B73 RefGen_V4 using the high-throughput short read aligner Bowtie. build_bowtie2DB.sub built the bowtie2 index on HPC-class, rename_chr.sh was a bash unix script that renamed each fasta header with the chromosome number, transcript_index.sub build a tophat index for the provided reference transcriptome gff file.  Running the tophat aligner was used in both forward and reverse using the run_tophat.sub (1 node, 16 threads, 4 hours) and run_tophat_rev.sub (5 nodes and 16 threads per node, 20 hours) scripts in HPC-class in parallel. We had to infer the insert size, standard deviations, and the arguments.  Each aligned SRR accession was collected as an accepted_hits.bam file.  

7. Count File Generation 

Accepted bam files were run using run_htseq.sub on HPC-class to generate HTseq-counts for each of the aligned SRR accessions (1 node, 16 threads, 20 hours or to completion)

8.  Differential Expression Analysis 

The differential expression analysis was performed locally on R version 3.4.0 under the DESeq2 package.  The raw counts table from the HTSeq analysis were formatted via the DESeq2.R script to obtain the differentially expressed genes

9.  Gene Differential Expression Results Figure Generation

The figures outlining the whole-transcriptome response of the maize genotypes were generated through the plotting function in R using gplots.

**Summary**

Our group reproduced a large-scale transcriptomic analysis of three maize genotypes exposed to individual stresses.  Utilizing a prototypical RNA-sequencing analysis workflow with modification set to work with a very large data set we replicated to the best of our abilities the data presented in Makarevitch et al.  Our final representation displayed as a heat map had general trends that matched the map as displayed by the original authors. 

1. Stress had substantial influence on the maize transcriptome in both maps

2. The three inbred lines in both maps showed similar transcriptional responses to the stress conditions 

3. Our map had slightly different clustering of the genotypes (ie. B73 clustered more closely with Mo17 during heat stress compared to Oh43)

You can find the figures in the figure directory. 

**Conclusions**

**Things We Learned**

1. How to modify and adapt scripts to permit the manipulation of large sets of data in order to obtain RNA-seq read information.

2. What conditions are necessary and what must be changed when using the splice junction mapper TopHat and it’s alignment partner Bowtie.

3. The ability to use use HTSeq-count to return raw count values from bam files generated via TopHat.

4. Figured out how to appropriately use DESeq software to generate the differential expression values.

**Obstacles Encountered**

1. There were key gaps in the analysis outlined in the paper versus the actual analysis. 

2. Little to no documentation was available resulting in lowered reproducibility by outside groups.

3. There were errors in uploading/downloading the data in which we encountered errors and corrupted files when analyzing data.

4. Our counting values were normalized while there’s was based on RPKM values which is something DESeq recommends not doing

5. Versioning of the packages (DESeq1 vs DESeq2) and the reference genome version seems to have effects on the genes annotated.

