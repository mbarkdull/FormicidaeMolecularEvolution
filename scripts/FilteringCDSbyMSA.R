library(phylotools)
library(tidyverse)
library(plyr)

# I've written a function to do what I want:
cdsFiltering <- function(cdsFile, msaFile, msaPrefix, output) {
  # Load packages:
  library(phylotools)
  library(plyr)
  library(tidyverse)
  # Read in coding sequences and split the sequence name column by the space:
  CodingSequences <- read.fasta(file = cdsFile)
  CodingSequences <- separate(data = CodingSequences, col = seq.name, into = c("name", "extraInfo"), sep = " ")
  # Read in protein sequences and remove the string "mpha_transcripts_" from the sequence names:
  MSAProteinSequences <- read.fasta(file = msaFile)
  MSAProteinSequences$seq.name <- gsub(msaPrefix, '', MSAProteinSequences$seq.name)
  # Create a column that checks if the coding sequence found in the protein sequence file?
  CodingSequences$inMSA <- CodingSequences$name %in% MSAProteinSequences$seq.name
  # Subset the coding sequences based on the value in that column:
  FilteredCodingSequences <- subset(CodingSequences, inMSA == "TRUE", select = c("name", "seq.text"))
  # Save the output
  dat2fasta(FilteredCodingSequences, outfile = output)
}

# Now I need to run this iteratively based on command line input:
