#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[1] = "/workdir/mb2337/FormicidaeMolecularEvolution/4_OrthoFinder/fasta/OrthoFinder/Results_Nov12/MultipleSequenceAlignments/"
}

# The command to run this script is `Rscript ./scripts/DataSubsetCDS.R ./scripts/inputurls_partial [path to Orthofinder MSA files]`, for example: `Rscript ./scripts/DataSubsetCDS.R ./scripts/inputurls_partial ./4_OrthoFinder/fasta/OrthoFinder/Results_Nov12/MultipleSequenceAlignments/`

# I need to break up nucleotide sequence files into many small files, like the orthogroup MSA files. 

library(phylotools)
library(plyr)
library(tidyverse)

# Read in the abbreviation data:
speciesInfo <- read.table(file = args[1], sep = ",")
#speciesInfo <- read.table(file = "./scripts/inputurls_partial", sep = ",")
# Split the second column to get a column with only abbreviations:
speciesInfo <- separate(data = speciesInfo, col = V2, into = c("abbrev", "transcript"), sep = "_")
# Get a vector from that column:
abbreviations <- speciesInfo$abbrev

# Make a working directory:
dir.create("./7_1_CDSOrthogroups")
# Copy the nucleotide alignments to the working directory:
file.copy("./6_PAL2NALOutput", "./7_1_CDSOrthogroups", recursive = TRUE)
# Copy the multiple sequence alignments to the directory:
file.copy(args[2], "./7_1_CDSOrthogroups", recursive = TRUE)
#file.copy("./MultipleSequenceAlignments", "./7_1_CDSOrthogroups", recursive = TRUE)
setwd("./7_1_CDSOrthogroups/6_PAL2NALOutput")

# Concatenate all of the cds files into a single file:
cdsFiles <- list.files(pattern = "*.fasta")
allCDSFiles <- bind_rows(lapply(cdsFiles, read.fasta))

# Fix the prefixes of the genes names in the coding sequences:
for (i in abbreviations)
{
  print(i)
  currentPrefix <- (paste(i, "_", sep = ""))
  print(currentPrefix)
  newPrefix <- (paste(i, "_transcripts_", i, "_", sep = ""))
  print (newPrefix)
  allCDSFiles$seq.name <- gsub(currentPrefix, newPrefix, allCDSFiles$seq.name)
}

setwd("../")

cdsSubsetting <- function(orthogroup, outfile){
  # Then check if the gene names in an MSA file are in the CDS file. 
  OrthogroupSequences <- read.fasta(orthogroup)
  allCDSFiles$OrthogroupMatch <- allCDSFiles$seq.name %in% OrthogroupSequences$seq.name
  # If they are, subset the cds sequences to a new file. 
  CDSOrthogroup <- subset(allCDSFiles, OrthogroupMatch == "TRUE", select = c("seq.name", "seq.text"))
  dat2fasta(CDSOrthogroup, outfile = outfile)
}

# Make a directory for output:
dir.create("./MultipleSequenceAlignments/Output")

# Now run this function over all orthogroups with a for loop. 
# Construct a list of all orthogroup files:
OrthogroupList <- list.files(path = "./MultipleSequenceAlignments", pattern = "*.fa", full.names = TRUE)
for (i in OrthogroupList)
{
  print(i)
  orthogroup <- i
  # Split up the value of i so it's just the orthogroup name, not with .fa 
  orthogoupName <- sapply(strsplit(i, "\\."), `[`, 2)
  print(orthogoupName)
  orthogoupName <- sapply(strsplit(orthogoupName, "\\/"), `[`, 3)
  print(orthogoupName)
  outfile <- (paste("./MultipleSequenceAlignments/Output/", orthogoupName, "_cds.fasta", sep = ""))
  print(outfile)
  cdsSubsetting(orthogroup = orthogroup, outfile = outfile)
}
