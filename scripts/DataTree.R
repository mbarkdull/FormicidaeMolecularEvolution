#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

# This script will generate a tree for Orthofinder and BUSTED based on an input tree and a set of tips to retain. 
# The command to run this is `Rscript ./scripts/DataTree.R [path to full tree] [path to inputurls file]`

library(ape)
library(phylotools)

# Read in the large tree file:
FullTree <- read.tree(file = args[1])

# Construct a vector of species to retain:
speciesInfo <- read.table(file = args[2], sep = ",")
StudySpecies <- speciesInfo$V3

#Construct a vector of file names corresponding to the species (for use with Orthofinder, and potentially with BUSTED)
# Split the second column by period (.) to get a column with just abbrev_transcripts:
speciesInfo <- separate(data = speciesInfo, col = V2, into = c("transcripts", ".fasta"), sep = "\\.")
filenames <- speciesInfo$transcripts

# Construct the species pruned tree:
StudyTree <- keep.tip(FullTree, StudySpecies)
# Export the species pruned tree:
write.tree(StudyTree, file = "./Tree/StudyTree.tre")

# Construct the filename pruned tree:
rename.tips <- function(phy, old_names, new_names) {
  mpos <- match(old_names,phy$tip.label)
  phy$tip.label[mpos] <- new_names
  return(phy)
}

FilenameTree <-  StudyTree
FilenameTree <-  rename.tips(phy = FilenameTree, old_names = StudySpecies, new_names = filenames)
# Export the filename pruned tree:
write.tree(FilenameTree, file = "./Tree/FilenameTree.tree")

