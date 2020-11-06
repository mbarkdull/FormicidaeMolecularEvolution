#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

# Make a directory for outputs:
dir.create("./FilteredCDS")

# Load packages:
library(phylotools)
library(plyr)
library(tidyverse)

# Read in the input data:
# The only thing left to do is figure out how to assign this at the command line:
speciesInfo <- read.table(file = args[1], sep = ",")
# Split the second column to get a column with only abbreviations:
speciesInfo <- separate(data = speciesInfo, col = V2, into = c("abbrev", "transcript"), sep = "_")
# Get a vector from that column:
abbreviations <- speciesInfo$abbrev

# I've written this function that will produce a subsetted nucleotide sequence file:
# I should also add error messages into this function. 
cdsFiltering <- function(cdsFile, msaFile, msaPrefix, output) {
  # Load packages:
  library(phylotools)
  library(plyr)
  library(tidyverse)
  # Read in coding sequences and split the sequence name column by the space:
  CodingSequences <- read.fasta(file = "./CodingSequences/cds_mpha_transcripts.fasta")
  CodingSequences <- separate(data = CodingSequences, col = seq.name, into = c("seq.name", "extraInfo"), sep = " ")
  CodingSequences <- select(CodingSequences, -c("extraInfo"))
  # Read in protein sequences and remove the string "mpha_transcripts_" from the sequence names:
  MSAProteinSequences <- read.fasta(file = "./Proteins/proteins_mpha.fasta")
  MSAProteinSequences$seq.name <- gsub("mpha_transcripts_", '', MSAProteinSequences$seq.name)
  # Create a column that checks if the coding sequence found in the protein sequence file?
  CodingSequences$inMSA <- CodingSequences$seq.name %in% MSAProteinSequences$seq.name
  # Subset the coding sequences based on the value in that column:
  FilteredCodingSequences <- subset(CodingSequences, inMSA == "TRUE", select = c("seq.name", "seq.text"))
  # Save the output
  dat2fasta(FilteredCodingSequences, outfile = "test.fasta")
}

# This for loop iterates the subsetting function on all species:
for (i in abbreviations)
{
  print(i)
  cdsFile <- (paste("./CodingSequences/cds_", i, "_transcripts.fasta", sep = ""))
  msaFile <- (paste("./Proteins/proteins_", i, ".fasta", sep = ""))
  msaPrefix <- (paste(i, "_transcripts_", sep = ""))
  output <- (paste("./FilteredCDS/filtered_", i, "_cds.fasta", sep = ""))
  print(output)
  print("__________________________________________________")
  
  cdsFiltering(cdsFile = cdsFile, msaFile = msaFile, msaPrefix = msaPrefix, output = output)
  }

