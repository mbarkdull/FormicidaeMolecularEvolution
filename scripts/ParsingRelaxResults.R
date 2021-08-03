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
  if (relaxResult[["test results"]][["p-value"]] < 0.05) {
    if (relaxResult[["test results"]][["relaxation or intensification parameter"]] > 1) {
      data <- c(relaxResult[["input"]][["file name"]], relaxResult[["test results"]][["p-value"]], relaxResult[["test results"]][["relaxation or intensification parameter"]], "Significant difference in selective regime between foreground and background branches", "Intensification of selection along foreground branches")
      return(data)
    } else {
      data <- c(relaxResult[["input"]][["file name"]], relaxResult[["test results"]][["p-value"]], relaxResult[["test results"]][["relaxation or intensification parameter"]], "Significant difference in selective regime between foreground and background branches", "Relaxation of selection along foreground branches")
      return(data)
    }
  } else {
    if (relaxResult[["test results"]][["relaxation or intensification parameter"]] > 1) {
      data <- c(relaxResult[["input"]][["file name"]], relaxResult[["test results"]][["p-value"]], relaxResult[["test results"]][["relaxation or intensification parameter"]], "No significant difference in selective regime between foreground and background branches", "Nonsignificant intensification")
      return(data)
    } else {
      data <- c(relaxResult[["input"]][["file name"]], relaxResult[["test results"]][["p-value"]], relaxResult[["test results"]][["relaxation or intensification parameter"]], "No significant difference in selective regime between foreground and background branches", "Nonsignificant relaxation")
      return(data)
    }
  }
}


possiblyRelaxJSONProcessing <- possibly(relaxJSONProcessing, otherwise = "File empty.")
# Run this function with purrr:map so as to avoid for loops (from https://www.r-bloggers.com/2017/12/skip-errors-in-r-loops-by-not-writing-loops/ and https://jennybc.github.io/purrr-tutorial/ls01_map-name-position-shortcuts.html):
relaxResults <- map(jsonFiles, possiblyRelaxJSONProcessing)

# Convert the results to a dataframe:
relaxResults <- as.data.frame(do.call(rbind, relaxResults))   
relaxResults$V2 <- as.numeric(as.character(relaxResults$V2), scientific = FALSE)
relaxResults$V3 <- as.numeric(as.character(relaxResults$V3), scientific = FALSE)


numberRelax <- sum(relaxResults$V5 == "Relaxation of selection along foreground branches")
numberNonsignficant <- sum(relaxResults$V4 == "No significant difference in selective regime between foreground and background branches")
numberIntensified <- sum(relaxResults$V5 == "Intensification of selection along foreground branches")
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






