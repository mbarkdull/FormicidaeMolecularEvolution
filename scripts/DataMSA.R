#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "/workdir/mb2337/FormicidaeMolecularEvolution/OrthoFinder/fasta/OrthoFinder/Results_*/MultipleSequenceAlignments"
}

# The command to run this script is `Rscript ./scripts/DataMSA.R [inputurls file] [path to Orthofinder MSA files]`, for example: `./scripts/DataMSA.R ./scripts/inputurls_partial /workdir/mb2337/FormicidaeMolecularEvolution/OrthoFinder/fasta/OrthoFinder/Results_Oct26/MultipleSequenceAlignments`

# This script will reshuffle the Orthofinder MSA files so that instead of one file per orthogroup with all species present, we have one file per specise with all orthogroups present. 

library(phylotools)
library(plyr)
library(tidyverse)

# First make a working directory and copy the folder with the MSA files there (THIS WILL NEED TO BE A VALUE SET ON THE COMMAND LINE). 
dir.create("./SpeciesMSA")
file.copy(args[2], "./SpeciesMSA", recursive = TRUE)
# Concatenate all of the MSA files into a single file:
setwd("./SpeciesMSA/MultipleSequenceAlignments")
msaFiles <- list.files(pattern = "*.fa")
allMSAFiles <- bind_rows(lapply(msaFiles, read.fasta))

# Create a function to create a subsetted fasta file for each species:
speciesChecking <- function(abbreviation, outFile){
    # Create a column in the concatenated file that tells us if it matches the species abbreviation:
    allMSAFiles$SpeciesMatch <- ifelse(grepl(abbreviation, allMSAFiles$seq.name), "True", "False")
    # Make a dataframe that contains only the rows where $SpeciesMatch == True. 
    SpeciesMSA <- subset(allMSAFiles, SpeciesMatch == "True", select = c("seq.name", "seq.text"))
    # Save that to a .fasta file: ../test.fasta
    dat2fasta(SpeciesMSA, outfile = outFile)
}

# Read in the inputurl file:
speciesInfo <- read.table(file = args[1], sep = ",")
# Split the second column to get a column with only abbreviations:
speciesInfo <- separate(data = speciesInfo, col = V2, into = c("abbrev", "transcript"), sep = "_")
# Get a vector from that column:
abbreviations <- speciesInfo$abbrev

# Write a for loop to run the function on all species:
for (i in abbreviations)
{
  print(i)
  abbreviation <- i
  print(abbreviation)
  outFile <- (paste("../", i, "_msa.fasta", sep = ""))
  print(outFile)
  print("__________________________________________________")
  
  speciesChecking(abbreviation = abbreviation, outFile = outFile)
}

setwd("../")
setwd("../")

