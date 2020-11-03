# Read in the input data:
# The only thing left to do is figure out how to assign this at the command line:
speciesInfo <- read.table(file = "./scripts/inputurls_partial", sep = ",")
# Split the second column to get a column with only abbreviations:
speciesInfo <- separate(data = speciesInfo, col = V2, into = c("abbrev", "transcript"), sep = "_")
# Get a vector from that column:
abbreviations <- speciesInfo$abbrev

# I've written this function that will produce a subsetted nucleotide sequence file:
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

# This while loop iterates the subsetting function on all species:
for (i in abbreviations)
{
  print(i)
  # The last thing I'll need to change is to make sure the paths to these things are correct in the context of the whole directory.
  cdsFile <- (paste("cds_", i, "_transcripts.fasta", sep = ""))
  print(cdsFile)
  msaFile <- (paste("proteins_", i, ".fasta", sep = ""))
  print(msaFile)
  msaPrefix <- (paste(i, "_transcripts_", sep = ""))
  print(msaPrefix)
  output <- (paste("filtered_", i, "_cds.fasta", sep = ""))
  print(output)
  
  cdsFiltering(cdsFile = cdsFile, msaFile = msaFile, msaPrefix = msaPrefix, output = output)
  }

