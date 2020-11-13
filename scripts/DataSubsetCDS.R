#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[1] = "/workdir/mb2337/FormicidaeMolecularEvolution/OrthoFinder/fasta/OrthoFinder/Results_*/OrthogroupSequences"
}

# The command to run this script is `Rscript ./scripts/DataMSA.R [path to Orthofinder MSA files]`, for example: `./scripts/DataMSA.R /workdir/mb2337/FormicidaeMolecularEvolution/OrthoFinder/fasta/OrthoFinder/Results_Oct26/MultipleSequenceAlignments`

# I need to break up nucleotide sequence files into many small files, like the orthogroup MSA files. 

library(phylotools)
library(plyr)
library(tidyverse)

# First make a working directory and copy the folder with the cds files there. 
dir.create("./CDSOrthogroups")
file.copy("./PAL2NALOutput", "./CDSOrthogroups", recursive = TRUE)
file.copy(args[1], "./CDSOrthogroups", recursive = TRUE)
setwd("./CDSOrthogroups/PAL2NALOutput")

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