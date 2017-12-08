**Materials and Methods**

**Plant growth and stress conditions**

Maize seeding of three genotypes B73, Mo17, and Oh43 were grown under natural light conditions and used for the study. All maize genotypes were subjected to four abiotic stress conditions namely heat stress, cold stress, salt stress and Ultraviolet  (UV) stress. 

For cold stress, seedlings were incubated at 5°C for 16 hours. For heat stress, seedlings were incubated at 50°C for 4 hours. For salt stress, plants were watered with 300 mM NaCl 20 hours prior to tissue collection. For UV stress, plants were subjected to UV  in the growth chamber conditions using UV-B lamps for 2 hours prior to tissue collection. Light conditions were the same for all stress and control conditions. 

14 days after germination, six seedlings were pooled together for each sample and then tissue collected for RNA extraction. Three biological replicates were used for heat, cold treated genotypes and controls for the genotypes.  one biological replicate was used for UV , high salt stressed and c. control genotypes

**RNA isolation and RNAseq analysis**

RNA was isolated using Trizol and purified with LiCl. All RNA samples were prepared in accordance with the TruSeq library creation protocol (Illumina, San Diego, CA). Samples were sequenced on the HiSeq 2000 to develop 10–20 million reads per sample. 

Transcript abundance was calculated by mapping reads to the combined transcript models of the maize reference genome (AGP v2) using TopHat software. Reads were filtered to allow for only uniquely mapped reads. RPKM values were developed using ‘BAM to Counts' across the exon space of the maize genome reference working gene set (ZmB73_5a) within the iPlant Discovery Environment.

Genes were only considered expressed if RPKM>1 and differentially expressed if log2(stress/control)> 1 or log2(stress/control) <-1. Statistical significance of expression differences was determined using DeSeq package for all fully replicated samples.  

**Group Documentation**

Every step of our data analysis was documented and files deposited on our group github repository [https://github.com/dormilon/BCB546X_final_project](https://github.com/dormilon/BCB546X_final_project) .  The repository was subdivided into directories or folder namely; - figures, paper,  fastqc_out, original_data, raw_tables, scripts and general readme. 

Directory figures contained the all visual summaries made using R software, paper directory contained our group documentation (md), fastqc_out contained quality outputs, original_data contained the research article being studied and raw data used (*however, dataset too big to be uploaded in the github repository (~86GB), so we have our raw data stored on an external hard drive*), raw_tables contained processed data summaries, scripts contained scripts ( with comments)  used for file inspection, processing and running the analysis. Above, we put a readme in every directory describing what the directory is composed of. 



