#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Must provide the path to BUSTED results as a command line argument.)", call.=FALSE)
}

###############################################
####### Convert JSON files to spreadsheet #####
###############################################

# This script reads in the results JSON files produced by BUSTED, and creates a list of orthogroups with evidence for positive selection somewhere in the tree. 
library(plyr)
library(tidyverse)
library(rjson)

# Construct a list of all of the json files:
jsonFiles <- list.files(path = args[1], pattern = "*.json", full.names = TRUE)
#jsonFiles <- list.files(path = "/Users/meganbarkdull/mb2337/FormicidaeMolecularEvolution/8_3_BustedResults", pattern = "*.json", full.names = TRUE)
jsonFiles <- sort(jsonFiles, decreasing = TRUE)
# bustedResult <- fromJSON(file = "/Users/meganbarkdull/mb2337/FormicidaeMolecularEvolution/8_3_BustedResults/OG0013390_busted.json")

# Write a function that will process each individual json file and extract the file name, orthogroup number, p-value, and return a text description of the p-value:
bustedJSONProcessing <- function(i) {
  bustedResults <- fromJSON(file = i)
  # Now run my if else statement:
  # If the p value is less than 0.05, that means there is positive selection. 
  if (bustedResults[["test results"]][["p-value"]] < 0.05) {
    print("There is evidence for positive selection.")
    orthogoupName <- sapply(strsplit(i, "\\/"), `[`, 7)
    orthogoupName <- sapply(strsplit(orthogoupName, "\\_"), `[`, 1)
    
    # Construct a vector of data containing the file name, the orthogroup number, the p-value, and the text "yes, evidence for positive selection":
    data <- c(bustedResults[["input"]][["file name"]], orthogoupName, bustedResults[["test results"]][["p-value"]], "yes, BUSTED found evidence for positive selection")
    return(data)
    
  } else {
    print("No positive selection.")
    orthogoupName <- sapply(strsplit(i, "\\/"), `[`, 7)
    orthogoupName <- sapply(strsplit(orthogoupName, "\\_"), `[`, 1)
    
    # Construct a vector of data containing the file name, the orthogroup number, the p-value, and the text "no evidence for positive selection":
    data <- c(bustedResults[["input"]][["file name"]], orthogoupName, bustedResults[["test results"]][["p-value"]], "no evidence for positive selection from BUSTED")
    return(data)
  }
}
# Create a version of the function that returns an error if there's an empty file (from https://www.r-bloggers.com/2017/12/skip-errors-in-r-loops-by-not-writing-loops/):
possiblyBustedJSONProcessing <- possibly(bustedJSONProcessing, otherwise = "File empty.")
# Run this function with purrr:map so as to avoid for loops (from https://www.r-bloggers.com/2017/12/skip-errors-in-r-loops-by-not-writing-loops/ and https://jennybc.github.io/purrr-tutorial/ls01_map-name-position-shortcuts.html):
bustedResults <- map(jsonFiles, possiblyBustedJSONProcessing)

# Convert the results to a dataframe:
bustedResults <- as.data.frame(do.call(rbind, bustedResults))   
bustedResults$V3 <- as.numeric(as.character(bustedResults$V3), scientific = FALSE)
colnames(bustedResults) <- c("file", "orthogroup", "test results p-value", "selectionOn")

# Get the number of genes with and without positive selection:
numberPositive <- sum(bustedResults$V4 == "yes, BUSTED found evidence for positive selection")
numberNone <- sum(bustedResults$V4 == "no evidence for positive selection from BUSTED")
percentPositive <- (numberPositive / (numberPositive + numberNone))*100

# Create an output directory:
dir.create("./Results")
# Export the results:
write_csv(bustedResults, "./Results/bustedResults.csv")

###############################################
####### Check for GO term enrichment ##########
###############################################

# Load packages:
library(topGO)
library(tidyverse)
library(RJSONIO)

# Read in the GO term annotations for each orthogroup, which we obtained from KinFin:
GOannotations <- read_delim("cluster_domain_annotation.GO.txt", delim = "\t")
# Subset so we just have orthogroup name and GO domain IDs:
longAnnotations <- dplyr::select(GOannotations, `#cluster_id`, domain_id)
# Take out genes without GO terms
longAnnotations <- longAnnotations[which(longAnnotations$domain_id != ""),] 
# Rename the column #cluster_id to orthogroup
longAnnotations <- longAnnotations %>%
  dplyr::rename(orthogroup = `#cluster_id`)

# Create list with element for each gene, containing vectors with all terms for each gene
wideListAnnotations <- tapply(longAnnotations$domain_id, longAnnotations$orthogroup, function(x)x)

# Define vector that is 1 if gene is significantly DE (`test results p-value` < 0.05) and 0 otherwise:
significanceInfo <- dplyr::select(bustedResults, orthogroup, `test results p-value`)
# Set each gene to 1 if adjP < cutoff, 0, otherwise
pcutoff <- 0.05 
tmp <- ifelse(significanceInfo$`test results p-value` < pcutoff , 1, 0)
geneList <- tmp

# Give geneList names:
names(geneList) <- significanceInfo$orthogroup

# Create the GOdata object:
GOdataBP <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
# Run Fisher's exact test to check for enrichment:
resultFisherBP <- runTest(GOdataBP, algorithm = "elim", statistic = "fisher")
resultFisherBP
resultsFisherBPTable <- GenTable(GOdataBP, raw.p.value = resultFisherBP, topNodes = length(resultFisherBP@score),
                                 numChar = 120)
GOdataMF <- new("topGOdata",
                ontology = "MF",
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
# Run Fisher's exact test to check for enrichment:
resultFisherMF <- runTest(GOdataMF, algorithm = "elim", statistic = "fisher")
resultFisherMF
resultsFisherMFTable <- GenTable(GOdataMF, raw.p.value = resultFisherMF, topNodes = length(resultFisherMF@score),
                                 numChar = 120)

GOdataCC <- new("topGOdata",
                ontology = "CC",
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
# Run Fisher's exact test to check for enrichment:
resultFisherCC <- runTest(GOdataCC, algorithm = "elim", statistic = "fisher")
resultFisherCC
resultsFisherCCTable <- GenTable(GOdataCC, raw.p.value = resultFisherCC, topNodes = length(resultFisherCC@score),
                                 numChar = 120)
# Use the Kolomogorov-Smirnov test:
geneList <- as.numeric(as.character(bustedResults$`test results p-value`))
names(geneList) <- bustedResults$orthogroup
# Create topGOData object
GOdataBP <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSelectionFun = function(x)x,
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
resultKSBP <- runTest(GOdataBP, algorithm = "weight01", statistic = "ks")
resultKSBP
resultKSBPTable <- GenTable(GOdataBP, raw.p.value = resultKSBP, topNodes = length(resultKSBP@score), numChar = 120)

# Create topGOData object
GOdataMF <- new("topGOdata",
                ontology = "MF",
                allGenes = geneList,
                geneSelectionFun = function(x)x,
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
resultKSMF <- runTest(GOdataMF, algorithm = "weight01", statistic = "ks")
resultKSMF
resultKSMFTable <- GenTable(GOdataMF, raw.p.value = resultKSMF, topNodes = length(resultKSMF@score), numChar = 120)


# Create topGOData object
GOdataCC <- new("topGOdata",
                ontology = "CC",
                allGenes = geneList,
                geneSelectionFun = function(x)x,
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
resultKSCC <- runTest(GOdataCC, algorithm = "weight01", statistic = "ks")
resultKSCC
resultKSCCTable <- GenTable(GOdataCC, raw.p.value = resultKSCC, topNodes = length(resultKSCC@score), numChar = 120)

goEnrichmentSummaries <- capture.output(print(resultFisherBP), 
                                        print(resultFisherMF), 
                                        print(resultFisherCC), 
                                        print(resultKSBP), 
                                        print(resultKSMF), 
                                        print(resultKSCC))

writeLines(goEnrichmentSummaries, con = file("./Results/bustedGOSummaries.csv"))




