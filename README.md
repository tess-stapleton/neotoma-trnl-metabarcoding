# This repository contains code associated with "Plant Metabarcoding as a Tool for Dietary Composition Analysis: Successes and Limitations"


File contents are as follows: 

QIIME_trnL_processing.txt: Initial processing of plant trnL g,h  amplicon sequences, initial filtering, quality control, clustering into OTUs, and taxonomic assignment. Amplicon sequences are available from the NCBI sequence archive under BioProject PRJNA766427. Inputs are fastq.gz files.

glmms_script.R: R script for all statistical analyses conducted for this manuscript. 

Constructing_trnL_database_and_taxonomy_assignment_with_BLAST+: all shell commands, Python scripts, and resulting DADA2 FASTA for constructing trnL reference database used in manuscript. Also includes shell commands and script used to assign sequences (assign_taxonomy.py) using BLAST+. 

R_trnL_processing.R: Initial processing of plant trnL g,h  amplicon sequences, initial filtering, quality control, clustering into OTUs in R . Amplicon sequences are available from the NCBI sequence archive under BioProject PRJNA766427. Inputs are fastq.gz files.
