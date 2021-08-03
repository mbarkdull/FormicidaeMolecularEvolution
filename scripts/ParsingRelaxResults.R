library(plyr)
library(tidyverse)
library(rjson)


# A significant result of k>1 indicates that selection strength has been intensified along the test branches, and a significant result of k<1 indicates that selection strength has been relaxed along the test branches. (https://stevenweaver.github.io/hyphy-site/methods/selection-methods/)


relaxResult <- fromJSON(file = "./9_3_RelaxResults/workerPolymorphism/OG0013353_relax.json")

relaxJSONProcessing <- function(i) {
  relaxResult <- fromJSON(file = i)
  # Now run my if else statement:
  # If the p value is less than 0.05, that means there is some kind of difference between the foreground and background. 
  if (relaxResult[["test results"]][["p-value"]] < 0.05) {
    print("Evidence for a difference in selective regime between foreground and background branches.")
    orthogoupName <- sapply(strsplit(i, "\\/"), `[`, 7)
    orthogoupName <- sapply(strsplit(orthogoupName, "\\_"), `[`, 1)
    
    # Get the K value, which tells us if selection is intensified or relaxed along the foreground:
    kValue <- relaxResult[["test results"]][["relaxation or intensification parameter"]]
    
    if (kValue > 1) {
      kValueDescriptive <- "Selection strength has been intensified along the test branches."
      # Construct a vector of data containing the file name, the orthogroup number, the p-value, and the text "yes, evidence for positive selection":
      data <- c(relaxResult[["input"]][["file name"]], orthogoupName, relaxResult[["test results"]][["p-value"]], "Evidence for a difference in selective regime between foreground and background branches.", kValue, kValueDescriptive)
      return(data)
    } else {
      kValueDescriptive <- "Selection strength has been relaxed along the test branches"
      # Construct a vector of data containing the file name, the orthogroup number, the p-value, and the text "yes, evidence for positive selection":
      data <- c(relaxResult[["input"]][["file name"]], orthogoupName, relaxResult[["test results"]][["p-value"]], "Evidence for a difference in selective regime between foreground and background branches.", kValue)
      return(data)
    }
  } else {
    print("No difference between foreground and background.")
    orthogoupName <- sapply(strsplit(i, "\\/"), `[`, 7)
    orthogoupName <- sapply(strsplit(orthogoupName, "\\_"), `[`, 1)
    
    if (kValue > 1) {
      kValueDescriptive <- "Nonsignificant increase in selection intensity on the test branches."
      data <- c(relaxResult[["input"]][["file name"]], orthogoupName, relaxResult[["test results"]][["p-value"]], "No difference between foreground and background", kValue, kValueDescriptive)
      return(data)
    } else {
      kValueDescriptive <- "Nonsignificant relaxation of selection on the test branches."
      data <- c(relaxResult[["input"]][["file name"]], orthogoupName, relaxResult[["test results"]][["p-value"]], "No difference between foreground and background", kValue, kValueDescriptive)
      return(data)
    }
  }
}
test <- relaxJSONProcessing("./9_3_RelaxResults/workerPolymorphism/OG0013353_relax.json")
print(test)
