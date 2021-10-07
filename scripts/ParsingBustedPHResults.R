library(RJSONIO)
library(tidyverse)
library(ggthemes)

###############################################
####### Convert JSON files to spreadsheet #####
###############################################

# Get a list of results files, sort by decreasing size, and remove any that are empty:
files <- list.files(path = "./BUSTEDPH", pattern = "*.json", full.names = TRUE)
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
# Create an output directory:
dir.create("./Results")
# Export the results:
write_csv(workerPolymorphismBustedPHResults, "./Results/bustedPHResults.csv")

###############################################
####### Fisher's exact test ###################
###############################################

backgroundSelected <- workerPolymorphismBustedPHResults %>% 
  filter(as.numeric(as.character(`test results p-value`)) > 0.05 & 
           as.numeric(as.character(`test results background p-value`)) <= 0.05 & 
           as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05)
nBackgroundSelected <- nrow(backgroundSelected)

foregroundSelected <- workerPolymorphismBustedPHResults %>% 
  filter(as.numeric(as.character(`test results p-value`)) <= 0.05 & 
           as.numeric(as.character(`test results background p-value`)) > 0.05 & 
           as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05)
nForegroundSelected <- nrow(foregroundSelected)

frequencyTableWithNoCorrection <- matrix(c(nBackgroundSelected, nForegroundSelected, 
                                           (nrow(workerPolymorphismBustedPHResults) - nBackgroundSelected),
                                           (nrow(workerPolymorphismBustedPHResults) - nForegroundSelected)),
                                         nrow = 2,
                                         byrow = TRUE,
                                         dimnames = list("category" = c("selected", "noSelection") ,
                                                         "selection" = c("backgroundSelected", "foregroundSelected")))
fishersExactTest <- fisher.test(frequencyTableWithNoCorrection)

fishersExactTest[["p.value"]]
if (fishersExactTest[["p.value"]] <= 0.05) {
  print("Signficant difference between foreground and background") 
} else {
  print("No signficant difference.")
}

fishersExactTest <- capture.output(print(fishersExactTest))

writeLines(fishersExactTest, con = file("./Results/bustedPHFisher.csv"))














