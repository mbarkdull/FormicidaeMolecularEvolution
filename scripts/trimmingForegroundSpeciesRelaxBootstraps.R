#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

# I need to trim out the previously foreground species for each bootstrapping analysis:
library(ape)
library(tidyverse)
library(phylotools)

# List the randomly selected tree files from boostrappingRelax.sh
treesToTrim <- list.files("./bootstrappingTreeFiles",
                          full.names = TRUE)

# Read in the set of foreground species that should now be removed:
#interest <- read_csv(file = "./ForegroundGroupings/workerReproductionQueens", col_names = FALSE)
interest <- read_csv(file = args[1], col_names = FALSE)

# Trim them out:
dir.create("./bootstrappingTreeFiles/trimmed/")
dir.create("./8_2_RemovedStops/sequences/")
trimmingForegroundSpecies <- function(treeFile) {
  orthogroup <- str_split_i(string = treeFile,
                            pattern = "/",
                            i = 3) %>%
    str_split_i(pattern = "_",
                i = 1)
  tree <- ape::read.tree(file = treeFile)
  tree[["tip.label"]]
  
  findingTipsToDrop <- function(species) {
    selectedTips <- grep(species, tree[["tip.label"]])
  }
  
  tipsToDrop <- purrr::map(interest$X1, findingTipsToDrop)
  tipsToDrop <- (do.call(c, tipsToDrop))
  tipsToDrop
  trimmedTree <- drop.tip(tree, tip = tipsToDrop)
  
  if (is.null(trimmedTree)) {
    print("This tree had only foreground tips")
  } else {
    print("Trimming tree.")
    ape::write.tree(trimmedTree, 
                    file = paste("./bootstrappingTreeFiles/trimmed/",
                                 orthogroup,
                                 "_tree.txt",
                                 sep = "")) 
    
    sequences <- read.fasta(file = paste("/workdir/mb2337/FormicidaeMolecularEvolution/8_2_RemovedStops/cleaned_",
                                         orthogroup,
                                         "_cds.fasta",
                                         sep = ""))
    sequences$species <- str_split_i(sequences$seq.name,
                                     pattern = "_",
                                     i = 1)
    sequences <- sequences[ ! sequences$species %in% interest$X1, ]
    dat2fasta(sequences, outfile = paste("./8_2_RemovedStops/sequences/cleaned_",
                                         orthogroup,
                                         "_cds.fasta",
                                         sep = ""))
  }
}

possiblyTrimmingForegroundSpecies <- possibly(trimmingForegroundSpecies,
                                              otherwise = "Error")
purrr::map(treesToTrim, 
           possiblyTrimmingForegroundSpecies)
