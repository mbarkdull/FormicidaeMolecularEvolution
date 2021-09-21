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
#jsonFiles <- list.files(path = "./10_1_RelaxResults/workerPolymorphism", pattern = "*.json", full.names = TRUE)

jsonFiles <- sort(jsonFiles, decreasing = TRUE)

relaxJSONProcessing <- function(i) {
  relaxResult <- fromJSON(file = i)
  if (relaxResult[["test results"]][["p-value"]] < 0.05) {
    if (relaxResult[["test results"]][["relaxation or intensification parameter"]] > 1) {
      data <- c(relaxResult[["input"]][["file name"]], relaxResult[["test results"]][["p-value"]], relaxResult[["test results"]][["relaxation or intensification parameter"]], "Significant difference in selective regime between foreground and background branches", "Intensification of selection along foreground branches", args[2])
      return(data)
    } else {
      data <- c(relaxResult[["input"]][["file name"]], relaxResult[["test results"]][["p-value"]], relaxResult[["test results"]][["relaxation or intensification parameter"]], "Significant difference in selective regime between foreground and background branches", "Relaxation of selection along foreground branches", args[2])
      return(data)
    }
  } else {
    if (relaxResult[["test results"]][["relaxation or intensification parameter"]] > 1) {
      data <- c(relaxResult[["input"]][["file name"]], relaxResult[["test results"]][["p-value"]], relaxResult[["test results"]][["relaxation or intensification parameter"]], "No significant difference in selective regime between foreground and background branches", "Nonsignificant intensification", args[2])
      return(data)
    } else {
      data <- c(relaxResult[["input"]][["file name"]], relaxResult[["test results"]][["p-value"]], relaxResult[["test results"]][["relaxation or intensification parameter"]], "No significant difference in selective regime between foreground and background branches", "Nonsignificant relaxation", args[2])
      return(data)
    }
  }
}


possiblyRelaxJSONProcessing <- possibly(relaxJSONProcessing, otherwise = "File empty.")
# Run this function with purrr:map so as to avoid for loops (from https://www.r-bloggers.com/2017/12/skip-errors-in-r-loops-by-not-writing-loops/ and https://jennybc.github.io/purrr-tutorial/ls01_map-name-position-shortcuts.html):
relaxResults <- map(jsonFiles, possiblyRelaxJSONProcessing)

# Convert the results to a dataframe:
relaxResults <- as.data.frame(do.call(rbind, relaxResults))   
colnames(relaxResults) <- c("fileName", "pValue", "kValue", "description", "shortDescription", "traitTested")
#relaxResults$V2 <- as.numeric(as.character(relaxResults$V2), scientific = FALSE)
#relaxResults$V3 <- as.numeric(as.character(relaxResults$V3), scientific = FALSE)


numberRelax <- sum(relaxResults$V5 == "Relaxation of selection along foreground branches")
numberNonsignficant <- sum(relaxResults$V4 == "No significant difference in selective regime between foreground and background branches")
numberIntensified <- sum(relaxResults$V5 == "Intensification of selection along foreground branches")
percentRelaxed <- (numberRelax / (numberRelax + numberNonsignficant + numberIntensified))*100
percentIntensified <- (numberIntensified / (numberRelax + numberNonsignficant + numberIntensified))*100

# Create an output directory:
dir.create("./Results")
outputDirectory <- paste("./Results/", args[2], sep = "")
dir.create(outputDirectory)
outputFile <- paste("./Results/", args[2], "/", args[2], "RelaxResults.csv", sep = "")
print(outputFile)
# Export the results:
write_csv(relaxResults, outputFile)

outputFile <- read_csv("./workerPolymorphismRelaxResults.csv")

# Assessing significance:
  # "The number of relaxed orthogroups in Coregonus and Salvelinus was significantly higher compared to the number of intensified orthogroups (one-sided Fisher’s Exact Test, p =0.035)." From Schneider, Kevin, Colin E. Adams, and Kathryn R. Elmer. "Parallel selection on ecologically relevant gene functions in the transcriptomes of highly diversifying salmonids." BMC genomics 20.1 (2019): 1-23. 	

### Kevin Schneider; last modified: 07.10.2019 ###
### KevinSchneider@gmx.at ###

library(ggplot2)

relaxframe <- outputFile
relaxframe <- subset(relaxframe, fileName != "File empty.")
relaxframe$qval <- NA
relaxframe$sig <- NA

relaxframe$qval <- p.adjust(relaxframe$pValue, method = "fdr", n = nrow(relaxframe))
Bonfcutoff <- (0.05 / nrow(relaxframe))

# We want to count the number of genes that are under relaxed and intensified selection. This person chose to do this with a for loop, so first create values set at zero. Each time the for loop encounters a row that would fit into one of these categories, it will increase the counter by 1. 
relaxedWithNoCorrection <- 0
relaxedWithFDR <- 0
relaxedWithBonferroni <- 0
intensifiedWithNoCorrection <- 0
intensifiedWithFDR <- 0
intensifiedWithBonferroni <- 0
for (i in 1:nrow(relaxframe)) {
  # If the p-value is significant and the k-value is less than 1, meaning selection has been relaxed, increase the relaxedWithNoCorrection counter by 1. 
  if ((relaxframe$pValue[i] < 0.05) && (relaxframe$kValue[i] < 1)) {
    relaxedWithNoCorrection <- relaxedWithNoCorrection + 1
    relaxframe$sig[i] <- "relaxedWithNoCorrection"
  } 
  if ((relaxframe$qval[i] < 0.10) && (relaxframe$kValue[i] < 1)) {
    relaxedWithFDR <- relaxedWithFDR + 1
    relaxframe$sig[i] <- "relaxedWithFDR"
  } 
  if ((relaxframe$pValue[i] < Bonfcutoff) && (relaxframe$kValue[i] < 1)) {
    relaxedWithBonferroni <- relaxedWithBonferroni + 1
    relaxframe$sig[i] <- "relaxedWithBonferroni"
  } 
  if ((relaxframe$pValue[i] < 0.05) && (relaxframe$kValue[i] > 1)) {
    intensifiedWithNoCorrection <- intensifiedWithNoCorrection + 1
    relaxframe$sig[i] <- "intensifiedWithNoCorrection"
  } 
  if ((relaxframe$qval[i] < 0.10) && (relaxframe$kValue[i] > 1)) {
    intensifiedWithFDR <- intensifiedWithFDR + 1
    relaxframe$sig[i] <- "intensifiedWithFDR"
  } 
  if ((relaxframe$pValue[i] < Bonfcutoff) && (relaxframe$kValue[i] > 1)) {
    intensifiedWithBonferroni <- intensifiedWithBonferroni + 1
    relaxframe$sig[i] <- "intensifiedWithBonferroni"
  }
}

frequencyTableWithNoCorrection <- matrix(c(relaxedWithNoCorrection, intensifiedWithNoCorrection, 
                                           (nrow(relaxframe) - relaxedWithNoCorrection),
                                           (nrow(relaxframe) - intensifiedWithNoCorrection)),
                                         nrow = 2,
                                         byrow = TRUE,
                                         dimnames = list("category" = c("outliers", "non-outliers") ,
                                                         "selection" = c("relaxed", "intensified")))
frequencyTableWithFDR <- matrix(c(relaxedWithFDR, intensifiedWithFDR, 
                                  (nrow(relaxframe) - relaxedWithFDR),
                                  (nrow(relaxframe) - intensifiedWithFDR)),
                                nrow = 2,
                                byrow = TRUE,
                                dimnames = list("category" = c("outliers", "non-outliers") ,
                                                "selection" = c("relaxed", "intensified")))

frequencyTableWithBonferroni <- matrix(c(relaxedWithBonferroni, intensifiedWithBonferroni, 
                                         (nrow(relaxframe) - relaxedWithBonferroni),
                                         (nrow(relaxframe) - intensifiedWithBonferroni)),
                                       nrow = 2,
                                       byrow = TRUE,
                                       dimnames = list("category" = c("outliers", "non-outliers") ,
                                                       "selection" = c("relaxed", "intensified")))

below1 <- relaxframe$kValue[which(relaxframe$kValue < 1)]
above1 <- relaxframe$kValue[which(relaxframe$kValue > 1)]

freq_table_below_or_above_1 <- matrix(c(length(below1), length(above1), 
                                        1347.5, 1347.5),
                                      nrow = 2,
                                      byrow = TRUE,
                                      dimnames = list("category" = c("outliers", "non-outliers") ,
                                                      "selection" = c("relaxed", "intensified")))


fisher.test(frequencyTableWithNoCorrection)
fisher.test(frequencyTableWithFDR)
fisher.test(frequencyTableWithBonferroni)
fisher.test(freq_table_below_or_above_1)

relaxedWithNoCorrection
relaxedWithFDR
relaxedWithBonferroni
intensifiedWithNoCorrection
intensifiedWithFDR
intensifiedWithBonferroni

write.csv(relaxframe, "RELAX_R_dataframe.csv")

median(relaxframe$kValue)
mean(relaxframe$kValue)

freq_table_selparam <- matrix(c((length(relaxframe$kValue[relaxframe$kValue != 1]) / 2), 
                                (length(relaxframe$kValue[relaxframe$kValue != 1]) / 2), 
                                length(relaxframe$kValue[relaxframe$kValue < 1]),
                                length(relaxframe$kValue[relaxframe$kValue > 1])),
                              nrow = 2,
                              byrow = TRUE,
                              dimnames = list("category" = c("expected", "observed") ,
                                              "selection" = c("relaxed", "intensified")))

fisher.test(freq_table_selparam)
