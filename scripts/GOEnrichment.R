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

# Read in BUSTED-PH results:
# Get a list of results files, sort by decreasing size, and remove any that are empty:
files <- list.files("./9_3_BustedPHResults/polygyny", pattern = "*.json", full.names = TRUE)
files <- sort(files, decreasing = TRUE)
files <- files[sapply(files, file.size) > 0]

#Write a function to read a single BUSTEDPH result file and extract relevant information from the JSON:
parsingBustedPH <- function(content) {
  result <- RJSONIO::fromJSON(content = content)
  if (result[["test results"]][["p-value"]] <= 0.05) {
    if (result[["test results background"]][["p-value"]] > 0.05) {
      if (result[["test results shared distributions"]][["p-value"]] <= 0.05) {
        print("Selection is associated with the phenotype / trait")
        textResult <- "Selection is associated with the phenotype / trait"
        data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "Yes")
        return(data)
      } else {
        print("Selection is associated with the phenotype / trait, but there is no significant difference between test and background branches in terms of selective pressure")
        textResult <- "Selection is associated with the phenotype / trait, but there is no significant difference between test and background branches in terms of selective pressure"
        data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "No")
      }
    } else {
      #print("Selection is acting on the branches with the phenotype / trait, but is **also** acting on background branches.")
      #textResult <- "Selection is acting on the branches with the phenotype / trait, but is **also** acting on background branches."
      #data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]])
      if (result[["test results shared distributions"]][["p-value"]] <= 0.05) {
        print("There is a significant difference between test and background branches in terms of selective pressure")
        textResult <- "There is a significant difference between test and background branches in terms of selective pressure"
        data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "Yes")
      } else {
        print("There is no significant difference between test and background branches in terms of selective pressure")
        textResult <- "There is no significant difference between test and background branches in terms of selective pressure"
        data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "No")
      }
    }
  } else {
    print("There is **no evidence** of episodic diversifying selection on test branches; selection is not associated with phenotype/trait")
    textResult <- "There is **no evidence** of episodic diversifying selection on test branches; selection is not associated with phenotype/trait"
    data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "No")
  }
  
}
# Make a safer version of that function with `possibly`:
possiblyparsingBustedPH <- possibly(parsingBustedPH, otherwise = "File empty.")
# Map the function over all results files to construct a master dataframe:
workerPolymorphismBustedPHResults <- map(files, possiblyparsingBustedPH)

# Write a function to process and fix the column names of the results dataframe:
bustedPHDataframeProcessing <- function(resultsDataframe) {
  BustedPHResults <- as.data.frame(do.call(rbind, resultsDataframe))   
  colnames(BustedPHResults) <- c("file", "textResult", "test results p-value", "test results background p-value", "test results shared distributions p-value", "differenceInSelection")
  orthogroup <- sapply(strsplit(as.character(BustedPHResults$file),"/"), tail, 1)
  orthogroup <- sapply(strsplit(orthogroup, "_"), `[`, 2)
  BustedPHResults$orthogroup <- orthogroup
  trait <- sapply(strsplit(as.character(BustedPHResults$file),"/"), tail, 1)
  trait <- sapply(strsplit(trait, "_"), `[`, 1)
  BustedPHResults$trait <- trait
  return(BustedPHResults)
}
workerPolymorphismBustedPHResults <- bustedPHDataframeProcessing(workerPolymorphismBustedPHResults)

workerPolymorphismBustedPHResults <- workerPolymorphismBustedPHResults %>% mutate(selectionOn =
                                                                                    case_when(as.numeric(as.character(`test results p-value`)) <= 0.05 & 
                                                                                                as.numeric(as.character(`test results background p-value`)) > 0.05 &
                                                                                                as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 ~ "ForegroundOnly", 
                                                                                              as.numeric(as.character(`test results p-value`)) > 0.05 &
                                                                                                as.numeric(as.character(`test results background p-value`)) <= 0.05 &
                                                                                                as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 ~ "BackgroundOnly", 
                                                                                              as.numeric(as.character(`test results p-value`)) > 0.05 &
                                                                                                as.numeric(as.character(`test results background p-value`)) > 0.05  ~ "NoSelection",
                                                                                              as.numeric(as.character(`test results shared distributions p-value`)) > 0.05 ~ "NoSignificantDifferenceBetweenForegroundAndBackground"))

### FIX THIS: we don't just care if the test p value is signficant. 
# Define vector that is 1 if gene is significantly DE (`test results p-value` < 0.05) and 0 otherwise:
significanceInfo <- dplyr::select(workerPolymorphismBustedPHResults, orthogroup, `test results p-value`, `test results background p-value`, `test results shared distributions p-value`)
# Set each gene to 1 if adjP < cutoff, 0, otherwise
pcutoff <- 0.05 
tmp <- ifelse(significanceInfo$`test results p-value` < pcutoff & 
                significanceInfo$`test results background p-value` > pcutoff & 
                significanceInfo$`test results shared distributions p-value` < pcutoff, 1, 0)
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
geneList <- as.numeric(as.character(workerPolymorphismBustedPHResults$`test results p-value`))
names(geneList) <- workerPolymorphismBustedPHResults$orthogroup
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

########################################################
# Check for GO term enrichment of plain BUSTED results:
########################################################

jsonFiles <- list.files(path = "./8_3_BustedResults", pattern = "*.json", full.names = TRUE)
jsonFiles <- sort(jsonFiles, decreasing = TRUE)

i <- "./8_3_BustedResults/OG0014518_busted.json"

# Write a function that will process each individual json file and extract the file name, orthogroup number, p-value, and return a text description of the p-value:
bustedJSONProcessing <- function(i) {
  bustedResults <- RJSONIO::fromJSON(content = i)
  # Now run my if else statement:
  # If the p value is less than 0.05, that means there is positive selection. 
  if (bustedResults[["test results"]][["p-value"]] < 0.05) {
    print("There is evidence for positive selection.")
    orthogoupName <- sapply(strsplit(i, "\\/"), tail, 1)
    orthogoupName <- sapply(strsplit(orthogoupName, "\\_"), `[`, 1)
    
    # Construct a vector of data containing the file name, the orthogroup number, the p-value, and the text "yes, evidence for positive selection":
    data <- c(bustedResults[["input"]][["file name"]], orthogoupName, bustedResults[["test results"]][["p-value"]], "yes, BUSTED found evidence for positive selection")
    return(data)
    
  } else {
    print("No positive selection.")
    orthogoupName <- sapply(strsplit(i, "\\/"), tail, 1)
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

# Define vector that is 1 if gene is significantly DE (`test results p-value` < 0.05) and 0 otherwise:
significanceInfo <- dplyr::select(bustedResults, V2, V3)
# Set each gene to 1 if adjP < cutoff, 0, otherwise
pcutoff <- 0.05 
tmp <- ifelse(significanceInfo$V3 < pcutoff , 1, 0)
geneList <- tmp

# Give geneList names:
names(geneList) <- significanceInfo$V2

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
geneList <- as.numeric(as.character(bustedResults$V3))
names(geneList) <- bustedResults$V2
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



