#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

library(orthologr)
library(tidyverse)
library(biomartr)
library(phylotools)
library(data.table)

# Read in the input urls file:
speciesInfo <- read.table(file = args[1], sep = ",")
# Get a list of species abbreviations:
species <- speciesInfo$V4

for (i in species) {
  print(i)
  proteomeFile <- paste("./1_RawData/",i,"_proteins.faa", sep = "")
  print(proteomeFile)
  annotationFile <- paste("./1_RawData/",i,"_GFF.gff", sep = "")
  print(annotationFile)
  newFile <- paste("./1_RawData/",i,"_longestIsoforms.fasta", sep = "")
  print(newFile)
  # Retrieve the longest isoforms for each species:
  retrieve_longest_isoforms(proteome_file = proteomeFile,
                            annotation_file = annotationFile,
                            new_file = newFile,
                            annotation_format = "gff")
}


