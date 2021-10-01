# TrnL Metabarcoding Sequence Analysis: 
  # Stapleton, Tess: Last updated 10/22/20 
  # this script is adapted from the following tutorials: 
  # https://benjjneb.github.io/dada2/tutorial.html
  # https://joey711.github.io/phyloseq/

#packages
library(Rcpp)
library(ggplot2)
library(ape)
library(gridExtra)
library(dada2)
library(S4Vectors)
library(stats4)
library(IRanges)
library(XVector)
library(RSQLite)
library(ShortRead)
library(Biostrings)
library(insect)
library(readxl)

# load fastq files into R
sequencepath <- "location/of/fastq/files"

#Sort to ensure fwd and reverse reads in same order
fnFs <- sort(list.files(sequencepath, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(sequencepath, pattern="_R2_001.fastq"))

# specify full path to the fnFs and fnRs
fnFs <- file.path(sequencepath, fnFs)
fnRs <- file.path(sequencepath, fnRs)
fnFs[1:3]
fnRs[1:3]

#################### removing primers before further processing #################
# CutAdapt: removing primers before further processing in dada2
# tell R where to find cutadapt
# https://cutadapt.readthedocs.io/en/stable/installation.html'
# code from: https://benjjneb.github.io/dada2/ITS_workflow.html
cutadapt <- "cutadapt/location" # version 2.4 
system2(cutadapt, args = "--version") # Run shell commands from R

FWD <- "GGGCAATCCTGAGCCAA" 
REV <- "CCATTGAGTCTCTGCACCTATC" #
# verify primers 
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# remove sequences with ambigous nucleotides: 
fnFs.filtN <- file.path(sequencepath, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(sequencepath, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# Counting the # of primers that appear in the fwd and rev reads
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# cutadapt removal of primers 
path.cut <- file.path(sequencepath, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "--discard-untrimmed", "--minimum-length", 8, 
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# check that the primers were removed: 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

# check quality of reads after primer removal: 
plotQualityProfile(fnFs[1:4]) # before primer removal
plotQualityProfile(cutFs[1:4]) # fwd after primer removal 
plotQualityProfile(cutRs[1:4]) # reverse after primer removal 

########################### Read Filtering and Trimming ##########################
# Filtering and Trimming of Reads for Quality Control 
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

f_path <- file.path(sequencepath, "filtered") 
if(!file_test("-d", f_path)) dir.create(f_path)
filtFs <- file.path(f_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(f_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Defining filtering parameters:  
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN=0, maxEE=c(2,2), rm.phix=TRUE, minQ = 20, truncQ = 0, minLen = 20, 
                     compress=TRUE, verbose = TRUE, multithread=TRUE) 

head(out)

# plot quality after filtering 
# make sure quality drop improves  
plotQualityProfile(filtFs[1:4])
plotQualityProfile(filtRs[1:4])

#Dereplication: 
exists <- file.exists(filtRs)
derepFs <- derepFastq(filtFs[exists], verbose=TRUE)
derepRs <- derepFastq(filtRs[exists], verbose=TRUE)

# call the derep-class objects by their sample names: 
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]

# estimate error rates 
errF <- learnErrors(filtFs[exists], multithread=TRUE, verbose = TRUE)
errR <- learnErrors(filtRs[exists], multithread=TRUE, verbose = TRUE)

# plot errors: 
plotErrors(errF)
plotErrors(errR)

#infer sequence variants at 100% identity
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, minOverlap = 10) 

# make a sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# remove artificial chimeras: 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#save sequence table 
saveRDS(seqtab.nochim, "<file_name>.rds")