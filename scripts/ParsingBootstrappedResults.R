# Load packages:
library(tidyverse)

# Function to read in all results from a single run, do FDR corrections, and get the number and proportion of genes under positive selection in the fore- and backgrounds. 
parsingBootstrapping <- function(i) {
  # Get a list of results files, sort by decreasing size, and remove any that are empty:
  files <- list.files(path = i, pattern = "*.json", full.names = TRUE)
  files <- sort(files, decreasing = TRUE)
  files <- files[sapply(files, file.size) > 0]
  
  #Write a function to read a single BUSTEDPH result file and extract relevant information from the JSON:
  parsingBustedPH <- function(content) {
    result <- RJSONIO::fromJSON(content = content)
    if (result[["test results"]][["p-value"]] <= 0.05) {
      if (result[["test results background"]][["p-value"]] > 0.05) {
        if (result[["test results shared distributions"]][["p-value"]] <= 0.05) {
          textResult <- "Selection is associated with the phenotype / trait"
          data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "Yes")
          return(data)
        } else {
          textResult <- "Selection is associated with the phenotype / trait, but there is no significant difference between test and background branches in terms of selective pressure"
          data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "No")
        }
      } else {
        if (result[["test results shared distributions"]][["p-value"]] <= 0.05) {
          textResult <- "There is a significant difference between test and background branches in terms of selective pressure"
          data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "Yes")
        } else {
          textResult <- "There is no significant difference between test and background branches in terms of selective pressure"
          data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "No")
        }
      }
    } else {
      textResult <- "There is **no evidence** of episodic diversifying selection on test branches; selection is not associated with phenotype/trait"
      data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "No")
    }
    
  }
  # Make a safer version of that function with `possibly`:
  possiblyparsingBustedPH <- purrr::possibly(parsingBustedPH, otherwise = "File empty.")
  # Map the function over all results files to construct a master dataframe:
  bootstrappingResults <- purrr::map(files, possiblyparsingBustedPH)
  
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
  bootstrappingResults <- bustedPHDataframeProcessing(bootstrappingResults)
  
  # Do FDR correction:
  bootstrappingResultsAdjusted <- bootstrappingResults %>%
    mutate(testResultspValueFDR = p.adjust(`test results p-value`, method='BH')) %>% 
    mutate(testResultsBackgroundpValueFDR = p.adjust(`test results background p-value`, method='BH')) %>% 
    mutate(testResultsSharedDistributionspValueFDR = p.adjust(`test results shared distributions p-value`, method='BH'))
  
  # Calculate the number of orthogroups under positive selection only in the foreground and only in the background:
  foreground <- length(which(bootstrappingResults$`test results p-value` <= 0.05 &
                               bootstrappingResults$`test results background p-value` > 0.05 &
                               bootstrappingResults$`test results shared distributions p-value` <= 0.05))
  background <- length(which(bootstrappingResults$`test results p-value` > 0.05 &
                               bootstrappingResults$`test results background p-value` <= 0.05 &
                               bootstrappingResults$`test results shared distributions p-value` <= 0.05))
  
  foregroundAdjusted <- length(which(bootstrappingResultsAdjusted$testResultspValueFDR <= 0.05 &
                                       bootstrappingResultsAdjusted$testResultsBackgroundpValueFDR > 0.05 &
                                       bootstrappingResultsAdjusted$testResultsSharedDistributionspValueFDR <= 0.05))
  backgroundAdjusted <- length(which(bootstrappingResultsAdjusted$testResultspValueFDR > 0.05 &
                                       bootstrappingResultsAdjusted$testResultsBackgroundpValueFDR <= 0.05 &
                                       bootstrappingResultsAdjusted$testResultsSharedDistributionspValueFDR <= 0.05))
  percentForeground <- foregroundAdjusted/length(bootstrappingResultsAdjusted$file)
  percentBackground <- backgroundAdjusted/length(bootstrappingResultsAdjusted$file)
  proportionForeToBack <- foregroundAdjusted/backgroundAdjusted
  results <- c(foregroundAdjusted, 
               backgroundAdjusted, 
               percentForeground,
               percentBackground,
               proportionForeToBack)
  names(results) <- c("foregroundAdjusted", 
                      "backgroundAdjusted", 
                      "percentForeground",
                      "percentBackground",
                      "proportionForeToBack")
  return(results)
}

# Make a safe version with possibly:
possiblyparsingBootstrapping <- possibly(parsingBootstrapping, otherwise = "Error")

# List all the bootstrapping runs:
runs <- list.files(path = "/workdir/mb2337/FormicidaeMolecularEvolution/9_3_BustedPHResults", 
                   pattern = "*bootstrapping",
                   full.names = TRUE)

# Iterate over all your runs:
resultsBootstrapping <- purrr::map(runs, possiblyparsingBootstrapping)
resultsBootstrapping <- as.data.frame(do.call(rbind, resultsBootstrapping))   
write.csv(resultsBootstrapping, 
          "./Results/threeSpeciesBootstrappingResults.csv")

# Plot the results as a barplot, marking the value for the real result:
ggplot(resultsBootstrapping) +
  geom_histogram(mapping = aes(x = proportionForeToBack),
                 binwidth = 0.02) +
  geom_vline(xintercept = 0.1288, 
             color = "red") +
  ggtitle("Summary of boostrap iterations\nwith three random foreground species") +
  xlab("Ratio of orthogroups positively selected\nin the foreground vs. background") +
  ylab("Count of bootstrap iterations") + 
  annotate("text", 
           x = 0.14, 
           y = 6.75, 
           label = "Ratio in actual results",
           color = "red",
           hjust = 0)

# Plot as a scatterplot:
ggplot(resultsBootstrapping) +
  geom_point(mapping = aes(x = percentForeground,
                           y = percentBackground))
