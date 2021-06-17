#!/usr/bin/env Rscript

# This script will identify the longest isoforms based on proteome and annotation (GFF) files, and will then filtered transcript files to include only genes that match the longest isoforms. 

# Read in the command line arguments; here, only one argument, the path to the input url/data file, needs to be provided.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
install.packages("orthologr", "biomartr")
library(orthologr)
library(tidyverse)
library(biomartr)
library(phylotools)
library(data.table)

# Read in the input urls file based on the command line argument:
speciesInfo <- read.table(file = args[1], sep = ",")
# If testing locally, can use this line:
# speciesInfo <- read.table(file = "./inputurls_full.txt", sep = ",")
# Filter to only include species where the "FilterIsoforms" attribute is yes. 
speciesInfo <- filter(speciesInfo, speciesInfo$V5 == "yes")
# Get a list of species abbreviations:
species <- speciesInfo$V4

# For each species i in the list of species, species:
for (i in species) {
  print(i)
  # Get the name of the proteome file that has been downloaded:
  proteomeFile <- paste("./1_RawData/",i,"_proteins.faa", sep = "")
  print(proteomeFile)
  # Get the name of the GFF annotation file that has been downloaded:
  annotationFile <- paste("./1_RawData/",i,"_GFF.gff", sep = "")
  print(annotationFile)
  # Create the name of the file to which the longest isoforms file should be saved:
  newFile <- paste("./1_RawData/",i,"_longestIsoforms.fasta", sep = "")
  print(newFile)
  # Retrieve the longest isoforms for each species:
  retrieve_longest_isoforms(proteome_file = proteomeFile,
                            annotation_file = annotationFile,
                            new_file = newFile,
                            annotation_format = "gff")
  
  # Get the name of the file to which the longest isoforms file should be saved:
  longestIsoformsFile <- paste("./1_RawData/",i,"_longestIsoforms.fasta", sep = "")
  print(longestIsoformsFile)
  # Get the name of the transcripts file that has been downloaded:
  transcriptsFile <- paste("./1_RawData/",i,"_transcripts.fasta", sep = "")
  print(transcriptsFile)
  # Read in the isoform and transcript files:
  isoforms <- phylotools::read.fasta(longestIsoformsFile)
  transcripts <- phylotools::read.fasta(transcriptsFile)
  
  # Do some processing to turn the transcript file gene names into something that will match the isoform gene names. 
    # This expects transcript gene names to look like this:
      # >lcl|NW_012130065.1_cds_XP_012054525.1_1 [gene=LOC105617575] [db_xref=GeneID:105617575] [protein=proton-coupled amino acid transporter 4-like] [frame=2] [partial=5'] [protein_id=XP_012054525.1] [location=join(<125..532,740..1100,1848..2438)] [gbkey=CDS]
    # And processes them to look like:
      # XP_012054525.1
  transcripts <- separate(transcripts, 
                          col = seq.name,
                          into = c("seq.name", "extra"),
                          sep = " ")
  transcripts <- separate(transcripts, 
                          col = seq.name,
                          into = c("prefix", "seq.name"),
                          sep = "cds_")
  transcripts <- separate(transcripts, 
                          col = seq.name,
                          into = c("seq.name", "extra"),
                          sep = "_")
  transcripts$seq.name <- paste(transcripts$seq.name, transcripts$extra, sep = "_")
  
  # Check if the transcript gene name is in the longest isoforms file:
  transcripts$check <- transcripts$seq.name %in% isoforms$seq.name
  # Filter to include only matches:
  longestTranscripts <- filter(transcripts, check == TRUE)
  # Select only the gene name and gene sequence columns:
  longestTranscripts <- select(longestTranscripts, seq.name, seq.text)
  # Create the file name to which the filtered transcripts should be saved:
  filteredTranscriptsOutput <- paste("./1_RawData/",i,"_filteredTranscripts.fasta", sep = "")
  # Output the filtered transcripts in fasta format. 
  dat2fasta(longestTranscripts, outfile = filteredTranscriptsOutput)
}
