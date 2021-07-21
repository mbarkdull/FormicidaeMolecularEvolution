# This script reads in the results JSON files produced by BUSTED, and creates a list of orthogroups with evidence for positive selection somewhere in the tree. 
library(plyr)
library(tidyverse)
library(rjson)

# Construct a list of all of the json files:
jsonFiles <- list.files(path = "/Users/meganbarkdull/mb2337/FormicidaeMolecularEvolution/8_3_BustedResults", pattern = "*.json", full.names = TRUE)
# bustedResults <- fromJSON(file = "/Users/meganbarkdull/mb2337/FormicidaeMolecularEvolution/8_3_BustedResults/OG0013829_busted.json")

# If the value of bustedResults[["test results"]][["p-value"]] is less than 0.05, append the value of bustedResults[["input"]][["file name"]] to a list
# Create an empty dataframe:
listOfResults <- data.frame(file = character(), orthogroup = character(), pValue = integer(), result = character(), stringsAsFactors = FALSE)

# Write a for loop to iterate over the json files:
for (i in jsonFiles) {
  print(i)
  # Read in the json file:
  bustedResults <- fromJSON(file = i)

  # Now run my if else statement:
  # If the p value is less than 0.05, that means there is positive selection. 
  if (bustedResults[["test results"]][["p-value"]] < 0.05) {
    print("There is evidence for positive selection")
    orthogoupName <- sapply(strsplit(i, "\\/"), `[`, 7)
    orthogoupName <- sapply(strsplit(orthogoupName, "\\_"), `[`, 1)
    
    # Construct a vector of data containing the file name, the orthogroup number, the p-value, and the text "yes, evidence for positive selection":
    data <- c(bustedResults[["input"]][["file name"]], orthogoupName, bustedResults[["test results"]][["p-value"]], "yes, BUSTED found evidence for positive selection")
    
    # Add that vector to the end of the results dataframe constructed earlier. 
    listOfResults[nrow(listOfResults) + 1, ] <- data  
  } else {
    print("no positive selection.")
    orthogoupName <- sapply(strsplit(i, "\\/"), `[`, 7)
    orthogoupName <- sapply(strsplit(orthogoupName, "\\_"), `[`, 1)
   
    # Construct a vector of data containing the file name, the orthogroup number, the p-value, and the text "no evidence for positive selection":
    data <- c(bustedResults[["input"]][["file name"]], orthogoupName, bustedResults[["test results"]][["p-value"]], "no evidence for positive selection from BUSTED")
    
    # Add that vector to the end of the results dataframe constructed earlier. 
    listOfResults[nrow(listOfResults) + 1, ] <- data   }
}

listOfResults$pValue <- as.numeric(listOfResults$pValue)
