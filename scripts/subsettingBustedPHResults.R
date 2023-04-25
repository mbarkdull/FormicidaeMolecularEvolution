# Load required packages:
library(RJSONIO)
library(tidyverse)
library(stringr)

# Get a list of results files, sort by decreasing size, and remove any that are empty:
files <- list.files(path = "./9_3_BustedPHResults/workerReproductionQueens", pattern = "*.json", full.names = TRUE)
files <- sort(files, decreasing = TRUE)
files <- files[sapply(files, file.size) > 0]

#### Check if an orthogroup contain obir genes: ####
subsettingBySpeciesPresence <- function(content) {
  # Read in the results files:
  result <- RJSONIO::fromJSON(content = content)
  
  # Check if the result contains sequences from O. biroi:
  anyObirMatches <- str_subset(names(result[["tested"]][["0"]]), "obir")
  if(length(anyObirMatches) < 1){
    return(result)
  }
}

subsetFiles <- purrr::map(files, subsettingBySpeciesPresence)

subsetFiles <- subsetFiles[-which(sapply(subsetFiles, is.null))]

extractResults <- function(input) {
  results <- c(`test results`= input$`test results`[2],
               `test results background`= input$`test results background`[2],
               `test results shared distributions`= input$`test results shared distributions`[2],
               tree = input$input$trees,
               file = input$input$`file name`)
}

subsetFiles2 <- purrr::map(subsetFiles, extractResults)

# Turn the results into a dataframe:
subsetFiles2 <- as.data.frame(do.call(rbind, subsetFiles2))   
subsetFiles <- subsetFiles2

# Get a column for the orthogroup:
subsetFiles <- separate(subsetFiles, col = file, into = c("junk1", "junk2", "junk3", "orthogroup"), sep = "_")

# Read in the hyphy results with FDR correction applied:
allHyphy <- read_csv(file = "./relaxAndBustedPH.csv")
workerReproductionBustedPH <- filter(allHyphy, 
                                     trait == "workerReproductionQueens", 
                                     orthogroup %in% subsetFiles$orthogroup)

# Create a dataframe with genes only under selection in the background:
backgroundSelected <- workerReproductionBustedPH %>% 
  filter(as.numeric(as.character(testResultspValueFDR)) > 0.05 &
           as.numeric(as.character(testResultsBackgroundpValueFDR)) <= 0.05 & 
           as.numeric(as.character(testResultsSharedDistributionspValueFDR)) <= 0.05)
nBackgroundSelected <- nrow(backgroundSelected)

# Create a dataframe with genes only under selection in the foreground:
foregroundSelected <- workerReproductionBustedPH %>% 
  filter(as.numeric(as.character(testResultspValueFDR)) <= 0.05 & 
           as.numeric(as.character(testResultsBackgroundpValueFDR)) > 0.05 & 
           as.numeric(as.character(testResultsSharedDistributionspValueFDR)) <= 0.05)
nForegroundSelected <- nrow(foregroundSelected)

nTotal <- workerReproductionBustedPH %>% 
  nrow()

# Construct a frequency table with that information:
selected <- c(nForegroundSelected, nBackgroundSelected)

# Run a test of equal proportions:
proportionTest <- chisq.test(selected)

# Print the resulting p-value:
proportionTest[["p.value"]]

#### Check the subsetted results for GO term enrichment: ####
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

goEnrichmentBUSTEDPH <- function(specificTrait) {
  significanceInfo <- dplyr::select(subsetFiles, 
                                    orthogroup, 
                                    testResultspValueFDR, 
                                    testResultsBackgroundpValueFDR, 
                                    testResultsSharedDistributionspValueFDR) 
  pcutoff <- 0.05 
  tmp <- ifelse(significanceInfo$testResultspValueFDR < pcutoff & 
                  significanceInfo$testResultsBackgroundpValueFDR > pcutoff & 
                  significanceInfo$testResultsSharedDistributionspValueFDR < pcutoff, 1, 0)
  geneList <- tmp
  names(geneList) <- significanceInfo$orthogroup
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
  
  goTermsForeground <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
    dplyr::filter(raw.p.value <= 0.01)
  goTermsForeground$selectionCategory <- "Positive selection in species with the trait"
  
  head(goTermsForeground)
  
  # Check GO term enrichment of genes under positive selection in background species:
  significanceInfo <- dplyr::select(bustedPHResults, 
                                    orthogroup, 
                                    testResultspValueFDR, 
                                    testResultsBackgroundpValueFDR, 
                                    testResultsSharedDistributionspValueFDR)
  significanceInfo <- dplyr::filter(significanceInfo)
  
  # Set each gene to 1 if adjP < cutoff, 0, otherwise
  pcutoff <- 0.05 
  tmp <- ifelse(significanceInfo$testResultspValueFDR > pcutoff & 
                  significanceInfo$testResultsBackgroundpValueFDR < pcutoff & 
                  significanceInfo$testResultsSharedDistributionspValueFDR < pcutoff, 1, 0)
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
  
  # Combine results of all the tests:
  goTermsBackground <- rbind(resultsFisherBPTable, 
                             resultsFisherMFTable, 
                             resultsFisherCCTable) %>%
    dplyr::filter(raw.p.value <= 0.01)
  goTermsBackground$selectionCategory <- "Positive selection in species lacking the trait"
  allGOTerms <- rbind(goTermsForeground,
                      goTermsBackground)
  allGOTerms$trait <- specificTrait
  allGOTerms
}

# Get a vector of traits:
traits <- unique(bustedPHResults$trait)

# Map the GO enrichment function over the traits:
goTermsBustedPH <- purrr::map(traits, goEnrichmentBUSTEDPH)








##### Create a function that will subset only orthogroups with equal representation from fore and background: ####
library(ape)
files <- list.files(path = "./9_1_LabelledPhylogenies/workerPolymorphism", pattern = "*.txt", full.names = TRUE)
files <- sort(files, decreasing = TRUE)
files <- files[sapply(files, file.size) > 0]

checkingBalancedRepresentation <- function(treeFile) {
  # Read in the tree:
  tree <- read.tree(file = treeFile)
  # Look for "{Foreground}" in the tip labels:
  foregroundPresent <- grepl("Foreground", tree[["tip.label"]])
  # Count the number of tips with and without Foreground:
  numberForeground <- length(which(foregroundPresent == TRUE))
  numberBackground <- length(which(foregroundPresent == FALSE))
  # Check if those two numbers are equal:
  checkingEquality <- numberForeground == numberBackground
  # If they are equal, get the file name:
  if(checkingEquality == TRUE){
    treeFileWithoutPath <- sapply(strsplit(as.character(treeFile),"/"), tail, 1)
    orthogroupTested <- sapply(strsplit(as.character(treeFileWithoutPath), "_tree"), head, 1)
    orthogroupTested <- str_replace(orthogroupTested, "Labelled", "")
    trait <- sapply(strsplit(as.character(orthogroupTested), "_"), head, 1)
    resultFile <- paste("./9_3_BustedPHResults/", trait, "/", orthogroupTested, "_BUSTEDPH.fas.BUSTED.json", sep = "")
    return(resultFile)
  }
}

balancedOrthogroupResults <- purrr::map(files, checkingBalancedRepresentation)
balancedOrthogroupResults <- balancedOrthogroupResults[-which(sapply(balancedOrthogroupResults, is.null))]
balancedOrthogroupResults <- as.data.frame(do.call(rbind, balancedOrthogroupResults))   
balancedOrthogroupResults <- as.character(balancedOrthogroupResults$V1)

#Write a function to read a single BUSTEDPH result file and extract relevant information from the JSON:
parsingBustedPH <- function(content) {
  result <- RJSONIO::fromJSON(content = content)
  if (result[["test results"]][["p-value"]] <= 0.05) {
    if (result[["test results background"]][["p-value"]] > 0.05) {
      if (result[["test results shared distributions"]][["p-value"]] <= 0.05) {
        textResult <- "Selection is associated with the phenotype / trait"
        data <- c(result[["input"]][["file name"]], 
                  textResult, 
                  result[["test results"]][["p-value"]], 
                  result[["test results background"]][["p-value"]], 
                  result[["test results shared distributions"]][["p-value"]], 
                  "Yes",
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["proportion"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["proportion"]])
        return(data)
      } else {
        textResult <- "Selection is associated with the phenotype / trait, but there is no significant difference between test and background branches in terms of selective pressure"
        data <- c(result[["input"]][["file name"]], 
                  textResult, 
                  result[["test results"]][["p-value"]], 
                  result[["test results background"]][["p-value"]], 
                  result[["test results shared distributions"]][["p-value"]], 
                  "No",
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["proportion"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["proportion"]])
      }
    } else {
      #print("Selection is acting on the branches with the phenotype / trait, but is **also** acting on background branches.")
      #textResult <- "Selection is acting on the branches with the phenotype / trait, but is **also** acting on background branches."
      #data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]])
      if (result[["test results shared distributions"]][["p-value"]] <= 0.05) {
        textResult <- "There is a significant difference between test and background branches in terms of selective pressure"
        data <- c(result[["input"]][["file name"]], 
                  textResult, 
                  result[["test results"]][["p-value"]], 
                  result[["test results background"]][["p-value"]], 
                  result[["test results shared distributions"]][["p-value"]], 
                  "Yes",
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["proportion"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["proportion"]])
      } else {
        textResult <- "There is no significant difference between test and background branches in terms of selective pressure"
        data <- c(result[["input"]][["file name"]], 
                  textResult, 
                  result[["test results"]][["p-value"]], 
                  result[["test results background"]][["p-value"]], 
                  result[["test results shared distributions"]][["p-value"]], 
                  "No",
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["proportion"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["proportion"]])
      }
    }
  } else {
    textResult <- "There is **no evidence** of episodic diversifying selection on test branches; selection is not associated with phenotype/trait"
    data <- c(result[["input"]][["file name"]], 
              textResult, 
              result[["test results"]][["p-value"]], 
              result[["test results background"]][["p-value"]], 
              result[["test results shared distributions"]][["p-value"]], 
              "No",
              result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["omega"]],
              result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["proportion"]],
              result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["omega"]],
              result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["proportion"]])
  }
  
}
# Make a safer version of that function with `possibly`:
possiblyparsingBustedPH <- possibly(parsingBustedPH, otherwise = content)
# Map the function over all results files to construct a master dataframe:
balancedBustedPHResults <- purrr::map(balancedOrthogroupResults, possiblyparsingBustedPH)

# Write a function to process and fix the column names of the results dataframe:
bustedPHDataframeProcessing <- function(resultsDataframe) {
  BustedPHResults <- as.data.frame(do.call(rbind, resultsDataframe))   
  colnames(BustedPHResults) <- c("file", 
                                 "textResult", 
                                 "test results p-value", 
                                 "test results background p-value", 
                                 "test results shared distributions p-value", 
                                 "differenceInSelection",
                                 "unconstrainedTestOmega",
                                 "unconstrainedTestProportion",
                                 "unconstrainedBackgroundOmega",
                                 "unconstrainedBackgroundProportion")
  orthogroup <- sapply(strsplit(as.character(BustedPHResults$file),"/"), tail, 1)
  orthogroup <- sapply(strsplit(orthogroup, "_"), `[`, 2)
  BustedPHResults$orthogroup <- orthogroup
  trait <- sapply(strsplit(as.character(BustedPHResults$file),"/"), tail, 1)
  trait <- sapply(strsplit(trait, "_"), `[`, 1)
  BustedPHResults$trait <- trait
  return(BustedPHResults)
}
balancedBustedPHResults <- bustedPHDataframeProcessing(balancedBustedPHResults)

# Create a dataframe with genes only under selection in the background:
backgroundSelected <- balancedBustedPHResults %>% 
  filter(as.numeric(as.character(`test results p-value`)) > 0.05 &
           as.numeric(as.character(`test results background p-value`)) <= 0.05 & 
           as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05)
nBackgroundSelected <- nrow(backgroundSelected)

# Create a dataframe with genes only under selection in the foreground:
foregroundSelected <- balancedBustedPHResults %>% 
  filter(as.numeric(as.character(`test results p-value`)) <= 0.05 & 
           as.numeric(as.character(`test results background p-value`)) > 0.05 & 
           as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05)
nForegroundSelected <- nrow(foregroundSelected)

# Construct a frequency table with that information:
frequencyTableWithNoCorrection <- matrix(c(nBackgroundSelected, nForegroundSelected, 
                                           (nrow(balancedBustedPHResults) - nBackgroundSelected),
                                           (nrow(balancedBustedPHResults) - nForegroundSelected)),
                                         nrow = 2,
                                         byrow = TRUE,
                                         dimnames = list("category" = c("selected", "noSelection") ,
                                                         "selection" = c("backgroundSelected", "foregroundSelected")))
# Run a fisher's exact test to see if there is a difference in the proportion of genes under selection in the fore- vs. background:
fishersExactTest <- fisher.test(frequencyTableWithNoCorrection)
# Print the resulting p-value:
fishersExactTest[["p.value"]]
frequencyTableWithNoCorrection






