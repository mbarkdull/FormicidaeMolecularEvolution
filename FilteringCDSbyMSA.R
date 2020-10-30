library(phylotools)
library(tidyverse)
library(plyr)

# Read in coding sequences:
CodingSequences <- read.fasta(file = "./cds_mpha_transcripts.fasta")
# split the sequence name column by the space:
CodingSequences <- separate(data = CodingSequences, col = seq.name, into = c("name", "extraInfo"), sep = " ")

# Read in protein sequences:
MSAProteinSequences <- read.fasta(file = "./proteins_mpha.fasta")
# remove the string "mpha_transcripts_" from the sequence names:
MSAProteinSequences$seq.name <- gsub('mpha_transcripts_', '', MSAProteinSequences$seq.name)

# Is the coding sequence found in the protein sequence file?
CodingSequences$inMSA <- CodingSequences$name %in% MSAProteinSequences$seq.name

# Subset the coding sequences:
FilteredCodingSequences <- subset(CodingSequences, inMSA == "TRUE", select = c("name", "seq.text"))
