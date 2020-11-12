#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

# This script will generate a tree for BUSTED based on an input tree and a set of tips to retain. 
# The command to run this is `Rscript ./scripts/DataTree.R [path to full tree] [path to inputurls file]`

library(ape)
library(phylotools)

# Read in the large tree file:
FullTree <- read.tree(file = args[1])

# Construct a vector of species to retain:
speciesInfo <- read.table(file = args[2], sep = ",")
StudySpecies <- speciesInfo$V3

# Construct the pruned tree:
StudyTree <- keep.tip(FullTree, StudySpecies)

# Export the pruned tree:
write.tree(StudyTree, file = "./Tree/StudyTree.tre")



