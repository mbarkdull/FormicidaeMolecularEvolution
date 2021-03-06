#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

# Make a directory for outputs:
dir.create("./6_2_FilteredCDS")

# Load packages:
library(phylotools)
library(plyr)
library(tidyverse)

# Read in the input data:
speciesInfo <- read.table(file = args[1], sep = ",")
#speciesInfo <- read.table(file = "./scripts/inputurls_full.txt", sep = ",")

# Get a vector from the abbreviations column:
abbreviations <- speciesInfo$V4

# I've written this function that will produce a subsetted nucleotide sequence file:
# I should also add error messages into this function. 
# The cdsFile argument is the path to the coding sequences output by Transdecoder, something like: ./4_2_TransdecoderCodingSequences/cds_acol_filteredTranscripts.fasta
# The mdsFile argument is the path to the recombined MSA file output by DataMSA.R, something like: ./6_1_SpeciesMSA/acol_msa.fasta
# The msaPrefix argument is 
cdsFiltering <- function(cdsFile, msaFile, msaPrefix, filteredOutput, msaOutput) {
  # Load packages:
  library(phylotools)
  library(plyr)
  library(tidyverse)
  # Read in coding sequences and split the sequence name column by the space:
  CodingSequences <- read.fasta(file = cdsFile)
  CodingSequences <- separate(data = CodingSequences, col = seq.name, into = c("seq.name", "extraInfo"), sep = " ")
  CodingSequences <- select(CodingSequences, -c("extraInfo"))
  # Read in protein sequences and remove the string "mpha_transcripts_" from the sequence names:
  MSAProteinSequences <- read.fasta(file = msaFile)
  MSAProteinSequences$seq.name <- gsub(msaPrefix, '', MSAProteinSequences$seq.name)
  # Create a column that checks if the coding sequence found in the protein sequence file?
  CodingSequences$inMSA <- CodingSequences$seq.name %in% MSAProteinSequences$seq.name
  # Subset the coding sequences based on the value in that column:
  FilteredCodingSequences <- subset(CodingSequences, inMSA == "TRUE", select = c("seq.name", "seq.text"))
  # Save the output
  dat2fasta(FilteredCodingSequences, outfile = filteredOutput)
  dat2fasta(MSAProteinSequences, outfile = msaOutput)
}

# This for loop iterates the subsetting function on all species:
for (i in abbreviations)
{
  print(i)
  cdsFile <- (paste("./4_2_TransdecoderCodingSequences/cds_", i, "_filteredTranscripts.fasta", sep = ""))
  msaFile <- (paste("./6_1_SpeciesMSA/", i, "_msa.fasta", sep = ""))
  msaPrefix <- (paste(i, "_filteredTranscripts_", sep = ""))
  filteredOutput <- (paste("./6_2_FilteredCDS/filtered_", i, "_cds.fasta", sep = ""))
  print(filteredOutput)
  msaOutput <- (paste("./6_1_SpeciesMSA/proteins_", i, ".fasta", sep = ""))
  print(msaOutput)
  print("__________________________________________________")
  
  cdsFiltering(cdsFile = cdsFile, msaFile = msaFile, msaPrefix = msaPrefix, filteredOutput = filteredOutput, msaOutput = msaOutput)
  }

