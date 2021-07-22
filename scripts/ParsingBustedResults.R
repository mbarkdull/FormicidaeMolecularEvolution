# This script reads in the results JSON files produced by BUSTED, and creates a list of orthogroups with evidence for positive selection somewhere in the tree. 
library(plyr)
library(tidyverse)
library(rjson)

# Construct a list of all of the json files:
jsonFiles <- list.files(path = "/Users/meganbarkdull/mb2337/FormicidaeMolecularEvolution/8_3_BustedResults", pattern = "*.json", full.names = TRUE)
jsonFiles <- sort(jsonFiles, decreasing = TRUE)
# bustedResults <- fromJSON(file = "/Users/meganbarkdull/mb2337/FormicidaeMolecularEvolution/8_3_BustedResults/OG0011321_busted.json")

# If the value of bustedResults[["test results"]][["p-value"]] is less than 0.05, append the value of bustedResults[["input"]][["file name"]] to a list
# Write a function that will process the json file:
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
# Run this function with map so as to avoid for loops (from https://www.r-bloggers.com/2017/12/skip-errors-in-r-loops-by-not-writing-loops/ and https://jennybc.github.io/purrr-tutorial/ls01_map-name-position-shortcuts.html):
result <- map(jsonFiles, possiblyBustedJSONProcessing)
