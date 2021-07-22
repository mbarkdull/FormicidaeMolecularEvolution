# I need to be able to label all of my phylogenies in order to run aBSREL and RELAX. 

library(ape)
library(tidyverse)

# List all of the unlabelled tree files:
treeFiles <- list.files(path = "/Users/meganbarkdull/mb2337/FormicidaeMolecularEvolution/5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Resolved_Gene_Trees", full.names = TRUE)

# Write a function that can relabel any tree with any vector of species of interest:
multiTreeLabelling <- function(tree, speciesOfInterest, exportFile) {
  # Read in a phylogeny:
  tree <- ape::read.tree(tree)
  plot(tree)
  
  # Make your list of species of interest:
  interest <- speciesOfInterest
  
  # Write a function that runs over the tip labels, splits them to extract the species abbreviation prefix, checks if that prefix is in the species of interest vector, and if it is, appends {Foreground} to the end of the tip label. 
  labellingFunction <- function(i) {
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
  ape::write.tree(tree, file = exportFile)
  return(tree)
}

test <- multiTreeLabelling(tree = "/Users/meganbarkdull/mb2337/FormicidaeMolecularEvolution/5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Resolved_Gene_Trees/OG0012798_tree.txt", speciesOfInterest = c("pbar", "cvar"), exportFile = "test.txt")
plot(test)







