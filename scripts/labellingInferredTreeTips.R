# Load packages
library(tidyverse)
library(ape)
library(purrr)

# List the inferred trees:
# treeFiles <- list.files(path = "./inferredTrees", full.names = TRUE)
treeFiles <- list.files(path = args[1], full.names = TRUE)

multiTreeLabelling <- function(i) {
  tree <- ape::read.tree(i)
  
  # Write a function that runs over the tip labels, splits them to extract the species abbreviation prefix, and adds [abbreviation]_filteredTranscripts_. 
  labellingFunction <- function(i) {
    # i looks like: cvar_filteredTranscripts_cvar_CVAR_01478RA_p1
    # we are extracting the species abbreviation, for example here cvar
    species <- sapply(strsplit(i, "\\_"), `[`, 1)
    print(i)
    prefix <- paste(species, "_filteredTranscripts_", i, sep = "")
  }
  
  # Apply that to all of the tip labels with purrr:map.
  tree[["tip.label"]] <- purrr::map(tree[["tip.label"]], labellingFunction)
  tree[["tip.label"]]
  
  ape::write.tree(tree, file = i)
}

map(treeFiles, multiTreeLabelling)
