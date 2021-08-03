#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Must provide the full path to Relax results, and the prefix used when generation labelled phylogenies, as command line arguments.)", call.=FALSE)
}

library(plyr)
library(tidyverse)
library(rjson)


# A significant result of k>1 indicates that selection strength has been intensified along the test branches, and a significant result of k<1 indicates that selection strength has been relaxed along the test branches. (https://stevenweaver.github.io/hyphy-site/methods/selection-methods/)

jsonFiles <- list.files(path = args[1], pattern = "*.json", full.names = TRUE)
#jsonFiles <- list.files(path = "./9_3_RelaxResults/workerPolymorphism", pattern = "*.json", full.names = TRUE)

jsonFiles <- sort(jsonFiles, decreasing = TRUE)

relaxJSONProcessing <- function(i) {
  relaxResult <- fromJSON(file = i)
  # Now run my if else statement:
  # If the p value is less than 0.05, that means there is some kind of difference between the foreground and background. 
  if (relaxResult[["test results"]][["p-value"]] < 0.05) {
    print(i)
    print("Evidence for a difference in selective regime between foreground and background branches.")
    orthogoupName <- sapply(strsplit(i, "\\/"), `[`, 7)
    orthogoupName <- sapply(strsplit(orthogoupName, "\\_"), `[`, 1)
    
    # Get the K value, which tells us if selection is intensified or relaxed along the foreground:
    kValue <- relaxResult[["test results"]][["relaxation or intensification parameter"]]
    
    if (kValue > 1) {
      kValue <- relaxResult[["test results"]][["relaxation or intensification parameter"]]
      kValueDescriptive <- "Selection strength has been intensified along the test branches."
      # Construct a vector of data containing the file name, the orthogroup number, the p-value, and the text "yes, evidence for positive selection":
      data <- c(relaxResult[["input"]][["file name"]], orthogoupName, relaxResult[["test results"]][["p-value"]], "Evidence for a difference in selective regime between foreground and background branches.", kValue, kValueDescriptive)
      return(data)
    } else {
      kValue <- relaxResult[["test results"]][["relaxation or intensification parameter"]]
      kValueDescriptive <- "Selection strength has been relaxed along the test branches"
      # Construct a vector of data containing the file name, the orthogroup number, the p-value, and the text "yes, evidence for positive selection":
      data <- c(relaxResult[["input"]][["file name"]], orthogoupName, relaxResult[["test results"]][["p-value"]], "Evidence for a difference in selective regime between foreground and background branches.", kValue)
      return(data)
    }
  } else {
    print(i)
    print("No difference between foreground and background")
    orthogoupName <- sapply(strsplit(i, "\\/"), `[`, 7)
    orthogoupName <- sapply(strsplit(orthogoupName, "\\_"), `[`, 1)
    
    if (kValue > 1) {
      kValue <- relaxResult[["test results"]][["relaxation or intensification parameter"]]
      kValueDescriptive <- "Nonsignificant increase in selection intensity on the test branches."
      data <- c(relaxResult[["input"]][["file name"]], orthogoupName, relaxResult[["test results"]][["p-value"]], "No difference between foreground and background", kValue, kValueDescriptive)
      return(data)
    } else {
      kValue <- relaxResult[["test results"]][["relaxation or intensification parameter"]]
      kValueDescriptive <- "Nonsignificant relaxation of selection on the test branches."
      data <- c(relaxResult[["input"]][["file name"]], orthogoupName, relaxResult[["test results"]][["p-value"]], "No difference between foreground and background", kValue, kValueDescriptive)
      return(data)
    }
  }
}

possiblyRelaxJSONProcessing <- possibly(relaxJSONProcessing, otherwise = "File empty.")
# Run this function with purrr:map so as to avoid for loops (from https://www.r-bloggers.com/2017/12/skip-errors-in-r-loops-by-not-writing-loops/ and https://jennybc.github.io/purrr-tutorial/ls01_map-name-position-shortcuts.html):
relaxResults <- map(jsonFiles, possiblyRelaxJSONProcessing)

# Convert the results to a dataframe:
relaxResults <- as.data.frame(do.call(rbind, relaxResults))   
relaxResults$V3 <- as.numeric(as.character(relaxResults$V3), scientific = FALSE)

numberRelax <- sum(relaxResults$V6 == "Selection strength has been relaxed along the test branches")
numberNonsignficant <- sum(relaxResults$V4 == "No difference between foreground and background")
numberIntensified <- sum(relaxResults$V6 == "Selection strength has been intensified along the test branches.")
percentRelaxed <- (numberRelax / (numberRelax + numberNonsignficant + numberIntensified))*100
percentIntensified <- (numberIntensified / (numberRelax + numberNonsignficant + numberIntensified))*100

# Create an output directory:
dir.create("./Results")
outputDirectory <- paste("./Results/", args[2], sep = "")
dir.create(outputDirectory)
outputFile <- paste("./Results/", args[2], "/relaxResults.csv", sep = "")
print(outputFile)
# Export the results:
write_csv(relaxResults, outputFile)
