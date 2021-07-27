#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least three arguments must be supplied (the full path to Orthofinder's gene trees, a list of species abbreviations to label as foreground, and a prefix for the output files.)", call.=FALSE)
} else if (length(args)==1) {
  # default label:
  args[3] = "labelled"
}

# I need to be able to label all of my phylogenies in order to run aBSREL and RELAX. 

library(ape)
library(tidyverse)

# List all of the unlabelled tree files:
treeFiles <- list.files(path = args[1], full.names = TRUE)
#tree <- ape::read.tree("/Users/meganbarkdull/mb2337/FormicidaeMolecularEvolution/5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Resolved_Gene_Trees/OG0001224_tree.txt")
#interest <- read_csv(file = "WorkerPolymorphismLabelling.txt", col_names = FALSE)

interest <- read_csv(file = args[2], col_names = FALSE)

# Write a function that can relabel any tree with any vector of species of interest:
multiTreeLabelling <- function(tree, speciesOfInterest, exportFile) {
  # Read in a phylogeny:
  tree <- ape::read.tree(tree)
  plot(tree)
  
  # Assign the list of species of interest:
  interest <- speciesOfInterest
  
  # Write a function that runs over the tip labels, splits them to extract the species abbreviation prefix, checks if that prefix is in the species of interest vector, and if it is, appends {Foreground} to the end of the tip label. 
  labellingFunction <- function(i) {
    # i looks like: cvar_filteredTranscripts_cvar_CVAR_01478RA_p1
    # we are extracting the species abbreviation, for example here cvar
    species <- sapply(strsplit(i, "\\_"), `[`, 1)
    if (species %in% interest) {
      new <- paste(i, "{Foreground}", sep = "")
      print(new)
    } else {
      print(i)
    }
  }
  
  # Apply that to all of the tip labels with purrr:map.
  tree[["tip.label"]] <- map(tree[["tip.label"]], labellingFunction)
  plot(tree)
  
  
  selected_tips <- grep("\\{Foreground\\}", tree$tip.label)
  
  ## Finding the direct ancestor for each of these tips
  descendants <- tree$edge[tree$edge[, 2] %in% selected_tips, 1]
  
  ## Adding the term "{Foreground}" to the selected descendants (the -Ntip(tree) part is because the node counting in the $edge table starts at the value Ntip(tree)).
  tree$node.label[descendants-Ntip(tree)] <- paste(tree$node.label[descendants-Ntip(tree)], "{Foreground}")
  ## Replacing all the non selected node labels by nothing ("")
  #tree$node.label[-c(descendants-Ntip(tree))] <- ""
  
  ## Plotting the results
  plot(tree, show.node.label = TRUE)
  ape::write.tree(tree, file = exportFile)
  return(tree)
}

test <- multiTreeLabelling(tree = "/Users/meganbarkdull/mb2337/FormicidaeMolecularEvolution/5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Resolved_Gene_Trees/OG0001224_tree.txt", speciesOfInterest = interest$X1, exportFile = "test.txt")
plot(test, show.node.label = TRUE)

for (i in treeFiles) {
  print(i)
  orthogoupName <- sapply(strsplit(i, "\\/"), `[`, 11)
  orthogoupName <- paste(args[3], "Labelled_", orthogoupName, sep = "")
  print(orthogoupName)
  multiTreeLabelling(tree = i, speciesOfInterest = interest, exportFile = orthogoupName)
}
