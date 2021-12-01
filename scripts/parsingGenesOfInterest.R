# Load necessary packages (I'm using RJSONIO here because rjson threw errors for one of my results files):
library(RJSONIO)
library(tidyverse)
library(janitor)


# Write a function to read in the aBSREL JSON result file and transform it into a dataframe:
absrelJSONProcessing <- function(file) {
  # Read in the JSON file:
  aBSRELJSON <- RJSONIO::fromJSON(content = file)
  
  # Create a data frame with the information for each branch tested:
  aBSRELDataframe <- data.table::rbindlist(aBSRELJSON[["branch attributes"]][["0"]], fill = TRUE, idcol = TRUE)
  aBSRELDataframe <- subset(aBSRELDataframe, select = -c(`Rate Distributions`))
  aBSRELDataframe <- distinct (aBSRELDataframe)
  
  #Add a column to hold the orthogroup name:
  orthogroup <- sapply(strsplit(aBSRELJSON[["input"]][["file name"]],"/"), tail, 1)
  orthogroup <- sapply(strsplit(orthogroup, "_"), `[`, 2)
  aBSRELDataframe$Orthogroup <- orthogroup

  return(aBSRELDataframe)
}
aBSRELJSON <- RJSONIO::fromJSON(content = "./12_genesOfInterest/absrel/OG0000796_absrel.json")

# List all of the results files:
aBSRELResults <- list.files(path = "./12_genesOfInterest/absrel/", 
                            pattern = "*.json", 
                            full.names = TRUE)
# Drop any files with file size of zero:
aBSRELResults <- aBSRELResults[sapply(aBSRELResults, file.size) > 0]

# Make a "Safe" version of the function that will return an error if there is a problem with a particular file:
possiblyabsrelJSONProcessing <- possibly(absrelJSONProcessing, 
                                         otherwise = "File empty.")
# Map the function over all foreground files:
aBSRELResults <- map_dfr(aBSRELResults, 
                         possiblyabsrelJSONProcessing)
significantGenes <- filter(aBSRELResults, 
                           `Corrected P-value` < 0.05)

# Read in the gene info:
info <- read_delim(file = "orthogroupsOfInterestFullInfo.txt", 
                   delim = "\t",
                   col_names = FALSE)
info$X2 <- gsub("^.*/", "", info$X2)

# Join it to the aBSREL results:
aBSRELResults <- inner_join(aBSRELResults, 
                            info, 
                            by = c("Orthogroup" = "X1"))

ggplot(data = filter(aBSRELResults, !is.na(`original name`))) +
  geom_point(mapping = aes(x = .id,
                           y = `Corrected P-value`)) + 
  geom_hline(yintercept = 0.05) +
  facet_wrap("X2")


# Plot the trees:
library(ggtree)
library(ape)
# I have to edit these functions from ggtree, because there's a bug that's making it not work with my data. 
  # The issue is in `dd <- dplyr::left_join(d1, d2, by="node")`, which I have changed to `dd <- dplyr::left_join(d1, d2, by=c("label" = "node"))`
`%<+%` <- function(pg, data) {
  if (! is.data.frame(data)) {
    stop("input should be a data.frame...")
  }
  pg %add% data
}
`%add%` <- function(p, data) {
  p$data <- p$data %add2% data
  return(p)
}
`%add2%` <- function(d1, d2) {
  if ("node" %in% colnames(d2)) {
    cn <- colnames(d2)
    ii <- which(cn %in% c("node", cn[!cn %in% colnames(d1)]))
    d2 <- d2[, ..ii]
    dd <- dplyr::left_join(d1, d2, by=c("label" = "node"))
  } else {
    d2[,1] <- as.character(unlist(d2[,1])) ## `unlist` to work with tbl_df
    d2 <- dplyr::rename(d2, label = 1) ## rename first column name to 'label'
    dd <- dplyr::left_join(d1, d2, by="label")
  }
  dd <- dd[match(d1$node, dd$node),]
  return(dd)
}

# Write a function to plot phylogeny, colored by p-value from aBSREL:
treeWithpValues <- function(orthogroup, aBSRELResults) {
  # Read in the tree file:
  tree <- read.tree(file = paste("5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Resolved_Gene_Trees/", 
                                 orthogroup, 
                                 "_tree.txt", 
                                 sep = ""))
  # Create a column in aBSRELResults called node, so that these data can be joined to the tree:
  aBSRELResults$node <- aBSRELResults$.id
  OG0004124 <- filter(aBSRELResults, 
                      Orthogroup == orthogroup) 
  OG0004124 <- mutate(OG0004124, 
                      node = as.character(node))
  OG0004124 <- dplyr::select(OG0004124, 
                             node, 
                             `Corrected P-value`, 
                             X2)
  OG0004124 <- mutate(OG0004124, 
                      Selected = if_else(`Corrected P-value` <= 0.05, 
                                         "Yes", 
                                         "No"))
  gene <- base::unique(OG0004124$X2)
  test <- ggtree(tree)
  tree <- test %<+% OG0004124
  tree[["data"]] <- full_join(tree[["data"]], OG0004124, by = c("label" = "node"))
  tree[["data"]][["label"]] <- sapply(strsplit(tree[["data"]][["label"]], "_"), `[`, 1)
  tree <- tree + 
    geom_tree(aes(color = `Corrected P-value.x`)) +
    geom_tiplab() +  
    scale_color_gradientn(colours=c("red", 
                                    'orange', 
                                    'green', 
                                    'cyan', 
                                    'blue')) +
    labs(title = paste("Pattern of positive selection in the gene family", 
                       gene, 
                       sep = " "))
  plot(tree)
}

orthogroups <- base::unique(aBSRELResults$Orthogroup)

plots <- map(orthogroups, ~treeWithpValues(.x, aBSRELResults = aBSRELResults))
plots[[1]]
cowplot::plot_grid(plotlist = plots)
ggplot2::ggsave(filename = "./genesOfInterestaBSRELTrees.pdf",
                width = 20,
                height = 20,
                units = "in")

treeWithBinarySelection <- function(orthogroup) {
  tree <- read.tree(file = paste("5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Resolved_Gene_Trees/", 
                                 orthogroup, 
                                 "_tree.txt", 
                                 sep = ""))
  aBSRELResults$node <- aBSRELResults$.id
  OG0004124 <- filter(aBSRELResults, 
                      Orthogroup == orthogroup & !is.na(`original name`)) 
  OG0004124 <- mutate(OG0004124, 
                      node = as.character(node))
  OG0004124 <- dplyr::select(OG0004124, 
                             node, 
                             `Corrected P-value`, 
                             X2)
  OG0004124 <- mutate(OG0004124, 
                      Selected = ifelse(OG0004124$`Corrected P-value` <= 0.05, 
                                         "Yes", 
                                         "No"))
  gene <- base::unique(OG0004124$X2)
  test <- ggtree(tree)
  tree <- test %<+% OG0004124
  tree[["data"]][["label"]] <- sapply(strsplit(tree[["data"]][["label"]], "_"), `[`, 1)
  tree <- tree + 
    geom_tree(aes(color = Selected)) +
    geom_tiplab() +
    labs(title = paste("Pattern of positive selection in the gene family", 
                       gene, 
                       sep = " "))
  plot(tree)
}

binaryPlots <- map(orthogroups, ~treeWithBinarySelection(.x))
binaryPlots[[1]]
cowplot::plot_grid(plotlist = binaryPlots)
ggplot2::ggsave(filename = "./genesOfInterestBinaryaBSRELTrees.pdf",
                width = 20,
                height = 20,
                units = "in")




tree <- read.tree(file = "5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Resolved_Gene_Trees/OG0000565_tree.txt")
aBSRELResults$node <- aBSRELResults$.id
OG0004124 <- filter(aBSRELResults, 
                    Orthogroup == "OG0000565") 
OG0004124 <- mutate(OG0004124, 
                    node = as.character(node))
OG0004124 <- dplyr::select(OG0004124, 
                           node, 
                           `Corrected P-value`, 
                           X2)
OG0004124 <- mutate(OG0004124, 
                    Selected = if_else(X2 <= 0.05, 
                                       "Yes", 
                                       "No"))
gene <- base::unique(OG0004124$X2)
test <- ggtree(tree)
tree <- test %<+% OG0004124

a <- tree[["data"]]

a <- inner_join(tree[["data"]], OG0004124, by = c("label" = "node"))
tree[["data"]] <- a  

tree[["data"]] <- mutate(tree[["data"]], 
                    Selected = ifelse(tree[["data"]]$`Corrected P-value.y` <= 0.05, 
                                      "Yes", 
                                      "No"))


x <- tree + 
  geom_tree(aes(color = Selected)) +
  geom_tiplab()
plot(x) 




