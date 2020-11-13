#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[1] = "/workdir/mb2337/FormicidaeMolecularEvolution/OrthoFinder/fasta/OrthoFinder/Results_*/MultipleSequenceAlignments"
}

# The command to run this script is `Rscript ./scripts/DataMSA.R [path to Orthofinder MSA files]`, for example: `./scripts/DataMSA.R /workdir/mb2337/FormicidaeMolecularEvolution/OrthoFinder/fasta/OrthoFinder/Results_Oct26/MultipleSequenceAlignments`

# I need to break up nucleotide sequence files into many small files, like the orthogroup MSA files. 

library(phylotools)
library(plyr)
library(tidyverse)

# First make a working directory and copy the folder with the cds files there. 
dir.create("./CDSOrthogroups")
file.copy("./TransdecoderCodingSequences", "./CDSOrthogroups", recursive = TRUE)
file.copy(args[1], "./CDSOrthogroups", recursive = TRUE)
setwd("./CDSOrthogroups/TransdecoderCodingSequences")

# Concatenate all of the cds files into a single file:
cdsFiles <- list.files(pattern = "*.fasta")
allCDSFiles <- bind_rows(lapply(cdsFiles, read.fasta))
allCDSFiles <- separate(data = allCDSFiles, col = seq.name, into = c("seq.name", "extraInfo"), sep = " ")
allCDSFiles <- select(allCDSFiles, -c("extraInfo"))

setwd("../")

cdsSubsetting <- function(orthogroup, outfile){
  # Then check if the gene names in an MSA file are in the CDS file. 
  OrthogroupSequences <- read.fasta(orthogroup)
  allCDSFiles$OrthogroupMatch <- allCDSFiles$seq.name %in% OrthogroupSequences$seq.name
  
  # If they are, subset the cds sequences to a new file. 
  CDSOrthogroup <- subset(allCDSFiles, OrthogroupMatch == "TRUE", select = c("seq.name", "seq.text"))
  dat2fasta(CDSOrthogroup, outfile = outfile)
  # Iterate that over all of the MSA files. 
}

# Now run this function over all orthogroups with a for loop. 
# Construct a list of all orthogroup files:
OrthogroupList <- list.files(path = "./OrthogroupSequences", pattern = "*.fa(?!\S)", full.names = TRUE)
for (i in OrthogroupList)
{
  print(i)
  orthogroup <- i
  # Split up the value of i so it's just the orthogroup name, not with .fa 
  orthogoupName <- sapply(strsplit(i, "\\."), `[`, 2)
  print(orthogoupName)
  outfile <- (paste(".", orthogoupName, "_cds.fasta", sep = ""))
  print(outfile)
  cdsSubsetting(orthogroup = orthogroup, outfile = outfile)
}

#############################
# This script will reshuffle the Orthofinder MSA files so that instead of one file per orthogroup with all species present, we have one file per specise with all orthogroups present. 

#library(phylotools)
#library(plyr)
#library(tidyverse)

# First make a working directory and copy the folder with the MSA files there. 
#dir.create("./SpeciesMSA")
#file.copy(args[2], "./SpeciesMSA", recursive = TRUE)
# Concatenate all of the MSA files into a single file:
#setwd("./SpeciesMSA/MultipleSequenceAlignments")
#msaFiles <- list.files(pattern = "*.fa")
#allMSAFiles <- bind_rows(lapply(msaFiles, read.fasta))

# Create a function to create a subsetted fasta file for each species:
#speciesChecking <- function(abbreviation, outFile){
  # Create a column in the concatenated file that tells us if it matches the species abbreviation:
  #allMSAFiles$SpeciesMatch <- ifelse(grepl(abbreviation, allMSAFiles$seq.name), "True", "False")
  # Make a dataframe that contains only the rows where $SpeciesMatch == True. 
  #SpeciesMSA <- subset(allMSAFiles, SpeciesMatch == "True", select = c("seq.name", "seq.text"))
  # Save that to a .fasta file: ../test.fasta
  #dat2fasta(SpeciesMSA, outfile = outFile)
#}

# Read in the inputurl file:
#speciesInfo <- read.table(file = args[1], sep = ",")
# Split the second column to get a column with only abbreviations:
#speciesInfo <- separate(data = speciesInfo, col = V2, into = c("abbrev", "transcript"), sep = "_")
# Get a vector from that column:
#abbreviations <- speciesInfo$abbrev

# Write a for loop to run the function on all species:
#for (i in abbreviations)
#{
  #print(i)
  #abbreviation <- i
  #print(abbreviation)
  #outFile <- (paste("../", i, "_msa.fasta", sep = ""))
  #print(outFile)
  #print("__________________________________________________")
  
  #speciesChecking(abbreviation = abbreviation, outFile = outFile)
#}

#setwd("../")
#setwd("../")


