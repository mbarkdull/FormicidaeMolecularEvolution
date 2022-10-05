# Load in packages:
library(tidyverse)
library(ape)
library(stringr)
library(purrr)

# Write a function that checks a tree to see if it has equal number of fore- and background tips, and returns the tree file only if it does.
subsettingBalancedTrees <- function(inputFile, equality) {
  # Read in a labelled tree:
  tree <- ape::read.tree(inputFile)
  # Count the number of foreground tips and background tips:
  # Create a function that checks each tip to see if it is labelled as foreground:
  countingTips <- function(tip) {
    foreground <- "oreground"
    grepl(foreground, tip)
  }

  # Make a safe version of this function with possibly:
  possiblyCountingTips <- possibly(countingTips, otherwise = "Error")
  
  # Map this function over all the tip labels, checking each one for the foreground label:
  tips <- purrr::map(tree[["tip.label"]], possiblyCountingTips)
  tips <- as.data.frame(do.call(rbind, tips))   
  # Count the number of tips that DO and DO NOT contain the foreground label:
  numberForeground <- length(which(tips$V1 == "TRUE"))
  numberBackground <- length(which(tips$V1 == "FALSE"))
  # See if those quantities are equal to one another:
  areEqual <- setequal(numberForeground, numberBackground) 
  # Depending on the input, check if the number of foreground tips is equal to, or greater than/equal to, the number of background tips:
  if (equality == "equal") {
    print("Checking for equal representation.")
    # If they are equal, return the input file name:
    if(numberForeground == numberBackground){
      print("Equal fore- and background")
      return(inputFile)
    } else {
      print("Unbalanced.")
    }
  } else if (equality == "greater") {
    print("Checking for equal or greater representation")
    # If they are equal, return the input file name:
    if(numberForeground >= numberBackground){
      print("Equal fore- and background")
      return(inputFile)
    } else {
      print("Unbalanced.")
    }
  } else {
    print("Error.")
  }
}

# Make a safe version of the tree selection function:
possiblySubsettingBalancedTrees <- possibly(subsettingBalancedTrees, otherwise = "Error")

# List all tree files:
workerPolymorphismTreeFiles <- list.files(path = "9_1_LabelledPhylogenies/workerPolymorphism", pattern = "*_tree.txt", full.names = TRUE)
workerReproductionTreeFiles <- list.files(path = "9_1_LabelledPhylogenies/workerReproductionQueens", pattern = "*_tree.txt", full.names = TRUE)

# Map it over all tree files:
workerPolymorphismBalancedTrees <- workerPolymorphismTreeFiles %>% 
  purrr::map(~ possiblySubsettingBalancedTrees(.x, "greater"))
workerReproductionBalancedTrees <- workerReproductionTreeFiles %>% 
  purrr::map(~ possiblySubsettingBalancedTrees(.x, "greater"))

workerPolymorphismBalancedTrees <- as.data.frame(do.call(rbind, workerPolymorphismBalancedTrees)) %>%
  dplyr::filter(V1 != "Unbalanced.")
workerPolymorphismBalancedTrees$trait <- "workerPolymorphism"
workerReproductionBalancedTrees <- as.data.frame(do.call(rbind, workerReproductionBalancedTrees)) %>%
  dplyr::filter(V1 != "Unbalanced.")
workerReproductionBalancedTrees$trait <- "workerReproductionQueens"

# Combine the lists of balanced trees for each trait:
allBalancedTrees <- rbind(workerPolymorphismBalancedTrees, workerReproductionBalancedTrees)

# Extract the orthogroup number from the tree file name:
getOrthogroup <- function(file) {
  filePieces <- strsplit(file, split = "_")
  orthogroup <- tail(filePieces[[1]], 2)[1]
  return(orthogroup)
}  
allBalancedTrees$orthogroups <- purrr::map(allBalancedTrees$V1, getOrthogroup)
allBalancedTrees$orthogroups <- as.character(allBalancedTrees$orthogroups)

# Read in the results:
selectionResults <- read_csv("./Results/relaxAndBustedPH.csv")
traits <- unique(selectionResults$trait)

# Subset to only balanced orthogroups:
allBalancedResults <- left_join(allBalancedTrees,
                                selectionResults,
                                by = c("trait" = "trait",
                                       "orthogroups" = "orthogroup"))

# Calculate p-values for each trait, on the FDR adjusted p-values:
# Get the RELAX p-values:
pValueByTraitRelax <- function(specificTrait, inputData) {
  relaxResultsTrait <- filter(inputData, 
                              trait == specificTrait & !is.na(kValue))
  nSelectionIntensified <- relaxResultsTrait %>% 
    filter(pValueFDR <= 0.05,
           kValue > 1) %>%
    nrow()
  
  nSelectionRelaxed <- relaxResultsTrait %>% 
    filter(pValueFDR <= 0.05,
           kValue < 1) %>%
    nrow()
  nTotal <- relaxResultsTrait %>% 
    nrow()
  # Construct a frequency table with that information:
  shiftInIntensity <- c(nSelectionIntensified, 
                        nSelectionRelaxed)
  population <- c((nrow(relaxResultsTrait) - nSelectionIntensified),
                  (nrow(relaxResultsTrait) - nSelectionRelaxed))
  proportionTest <- prop.test(shiftInIntensity,
                              population)
  # Print the resulting p-value:
  proportionTest[["p.value"]]
  pValue <- if (proportionTest[["p.value"]] < 0.000001) {
    formatC(proportionTest[["p.value"]], format = "e", digits = 2)
  } else {
    round(proportionTest[["p.value"]], digits = 4)
  }
  textHeight <- as.numeric(max(nSelectionIntensified, nSelectionRelaxed))
  return(c(specificTrait, pValue, textHeight, nSelectionIntensified, nSelectionRelaxed, nTotal))
}
possiblypValueByTraitRelax <- possibly(pValueByTraitRelax, otherwise = "Error")
pValuesRelax <- traits %>% 
  purrr::map(~ possiblypValueByTraitRelax(.x, allBalancedResults))
pValuesRelax <- as.data.frame(do.call(rbind, pValuesRelax))   
colnames(pValuesRelax) <- c("trait", "pValue", "maxHeight", "nSelectionIntensified", "nSelectionRelaxed", "nTotal")
pValuesRelax$test <- "relax"

# Get the BUSTED-PH p-values:
pValueByTraitBUSTEDPH <- function(specificTrait, inputData) {
  bustedPHResultsTrait <- filter(inputData, 
                                 trait == specificTrait & !is.na(`test results p-value`))
  nSelectionForeground <- bustedPHResultsTrait %>% 
    filter(testResultspValueFDR <= 0.05,
           testResultsBackgroundpValueFDR > 0.05,
           testResultsSharedDistributionspValueFDR <= 0.05) %>%
    nrow()
  
  nSelectionBackground <- bustedPHResultsTrait %>% 
    filter(testResultspValueFDR > 0.05,
           testResultsBackgroundpValueFDR <= 0.05,
           testResultsSharedDistributionspValueFDR <= 0.05) %>%
    nrow()
  nTotal <- bustedPHResultsTrait %>% 
    nrow()
  
  # Construct a frequency table with that information:
  selected <- c(nSelectionForeground, nSelectionBackground)
  population <- c((nrow(bustedPHResultsTrait) - nSelectionForeground),
                  (nrow(bustedPHResultsTrait) - nSelectionBackground))
  # Run a test of equal proportions:
  proportionTest <- prop.test(selected,
                              population)
  # Print the resulting p-value:
  proportionTest[["p.value"]]
  pValue <- if (proportionTest[["p.value"]] < 0.000001) {
    formatC(proportionTest[["p.value"]], format = "e", digits = 2)
  } else {
    round(proportionTest[["p.value"]], digits = 4)
  }
  textHeight <- as.numeric(max(nSelectionForeground, nSelectionBackground))
  return(c(specificTrait, pValue, textHeight, nSelectionForeground, nSelectionBackground, nTotal))
}
possiblypValueByTraitBUSTEDPH <- possibly(pValueByTraitBUSTEDPH, otherwise = "Error")
pValuesBUSTEDPH <- traits %>% 
  map(~ possiblypValueByTraitBUSTEDPH(.x, allBalancedResults))
pValuesBUSTEDPH <- as.data.frame(do.call(rbind, pValuesBUSTEDPH))   
colnames(pValuesBUSTEDPH) <- c("trait", "pValue", "maxHeight", "nSelectionForeground", "nSelectionBackground", "nTotal")
pValuesBUSTEDPH$test <- "bustedPH"

# Get all the p-values:
pValuesAll <- bind_rows(pValuesRelax, pValuesBUSTEDPH)
