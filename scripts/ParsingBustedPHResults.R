#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("At least two arguments required: the path to BUSTED-PH results, and the corresponding trait prefix.", call.=FALSE)
} else if (length(args)==1) {
  # default label:
  args[2] = "labelled"
}


library(RJSONIO)
library(tidyverse)
library(ggthemes)

###############################################
####### Convert JSON files to spreadsheet #####
###############################################

# Get a list of results files, sort by decreasing size, and remove any that are empty:
files <- list.files(path = args[1], pattern = "*.json", full.names = TRUE)
#files <- list.files(path = "./9_3_BustedPHResults/polygyny", pattern = "*.json", full.names = TRUE)
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
      #print("Selection is acting on the branches with the phenotype / trait, but is **also** acting on background branches.")
      #textResult <- "Selection is acting on the branches with the phenotype / trait, but is **also** acting on background branches."
      #data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]])
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
dir.create("./Results/")
dir.create(base::paste("./Results/", args[2], sep = ""))
# Export the results:
write_csv(workerPolymorphismBustedPHResults, base::paste("./Results/", args[2], "/bustedPHResults.csv", sep = ""))

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

writeLines(fishersExactTest, con = file(base::paste("./Results/", args[2], "/bustedPHFisher.csv", sep = "")))

###############################################
####### Make some nice plots ##################
###############################################

workerPolymorphismBustedPHResults$selectionOn <-
  factor(workerPolymorphismBustedPHResults$selectionOn,
         levels = c("ForegroundOnly",
                    "BackgroundOnly",
                    "NoSignificantDifferenceBetweenForegroundAndBackground",
                    "NoSelection",
                    NA))

plot <- ggplot(data = workerPolymorphismBustedPHResults)

plot + 
  geom_bar(mapping = aes(x = selectionOn)) + 
  theme_bw() +
  labs(x = "Selective regime", 
       y = "Count of orthogroups", 
       title = "Distribution of selective regimes") +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) + 
  scale_x_discrete(labels=c("ForegroundOnly" = "Selection on\nforeground only",
                            "BackgroundOnly" = "Selection on\nbackground only",
                            "NoSignificantDifferenceBetweenForegroundAndBackground" = "No significant\ndifference\nbetween fore- \nand background",
                            "NoSelection" = "No evidence\nfor selection\non foreground"))

###############################################
####### Check for GO term enrichment ##########
###############################################

# Load packages:
library(topGO)
library(tidyverse)
library(RJSONIO)

print("####################################################################")
print("#### GO term enrichment on genes under  selection in foreground ####")
print("####################################################################")

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

# Define vector that is 1 if gene is significantly DE (`test results p-value` < 0.05) and 0 otherwise:
significanceInfo <- dplyr::select(workerPolymorphismBustedPHResults, orthogroup, `test results p-value`, `test results background p-value`, `test results shared distributions p-value`)
# Set each gene to 1 if adjP < cutoff, 0, otherwise
pcutoff <- 0.05 
tmp <- ifelse(significanceInfo$`test results p-value` < pcutoff & 
                significanceInfo$`test results background p-value` > pcutoff & 
                significanceInfo$`test results shared distributions p-value` < pcutoff, 1, 0)
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
head(resultsFisherBPTable)
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
head(resultsFisherMFTable)

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
head(resultsFisherCCTable)

# Use the Kolomogorov-Smirnov test:

# Define vector that is 1 if gene is significantly DE (`test results p-value` < 0.05) and 0 otherwise:
significanceInfo <- dplyr::select(workerPolymorphismBustedPHResults, orthogroup, `test results p-value`, `test results background p-value`, `test results shared distributions p-value`)
# Set each gene to 1 if adjP < cutoff, 0, otherwise
significanceInfo <- filter(significanceInfo, significanceInfo$`test results p-value` < 0.05 & 
                             significanceInfo$`test results background p-value` > 0.05 & 
                             significanceInfo$`test results shared distributions p-value` < 0.05)

geneList <- as.numeric(as.character(significanceInfo$`test results p-value`))

# Give geneList names:
names(geneList) <- significanceInfo$orthogroup

# Create topGOData object
GOdataBP <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSelectionFun = function(x)x,
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
resultKSBP <- runTest(GOdataBP, algorithm = "weight01", statistic = "ks")
resultKSBP
resultKSBPTable <- GenTable(GOdataBP, raw.p.value = resultKSBP, topNodes = length(resultKSBP@score), numChar = 120)
head(resultKSBPTable)

# Create topGOData object
GOdataMF <- new("topGOdata",
                ontology = "MF",
                allGenes = geneList,
                geneSelectionFun = function(x)x,
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
resultKSMF <- runTest(GOdataMF, algorithm = "weight01", statistic = "ks")
resultKSMF
resultKSMFTable <- GenTable(GOdataMF, raw.p.value = resultKSMF, topNodes = length(resultKSMF@score), numChar = 120)
head(resultKSMFTable)

# Create topGOData object
GOdataCC <- new("topGOdata",
                ontology = "CC",
                allGenes = geneList,
                geneSelectionFun = function(x)x,
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
resultKSCC <- runTest(GOdataCC, algorithm = "weight01", statistic = "ks")
resultKSCC
resultKSCCTable <- GenTable(GOdataCC, raw.p.value = resultKSCC, topNodes = length(resultKSCC@score), numChar = 120)
head(resultKSCCTable)

goEnrichmentSummaries <- capture.output(print(resultFisherBP), 
                                        print(resultFisherMF), 
                                        print(resultFisherCC), 
                                        print(resultKSBP), 
                                        print(resultKSMF), 
                                        print(resultKSCC))

writeLines(goEnrichmentSummaries, con = file(base::paste("./Results/", args[2], "/bustedPHGOSummariesForeground.csv", sep = "")))


print("####################################################################")
print("#### GO term enrichment on genes under  selection in background ####")
print("####################################################################")

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

# Define vector that is 1 if gene is significantly DE (`test results p-value` < 0.05) and 0 otherwise:
significanceInfo <- dplyr::select(workerPolymorphismBustedPHResults, orthogroup, `test results p-value`, `test results background p-value`, `test results shared distributions p-value`)
# Set each gene to 1 if adjP < cutoff, 0, otherwise
pcutoff <- 0.05 
tmp <- ifelse(significanceInfo$`test results p-value` > pcutoff & 
                significanceInfo$`test results background p-value` < pcutoff & 
                significanceInfo$`test results shared distributions p-value` < pcutoff, 1, 0)
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
head(resultsFisherBPTable)
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
head(resultsFisherMFTable)

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
head(resultsFisherCCTable)

# Use the Kolomogorov-Smirnov test:

# Define vector that is 1 if gene is significantly DE (`test results p-value` < 0.05) and 0 otherwise:
significanceInfo <- dplyr::select(workerPolymorphismBustedPHResults, orthogroup, `test results p-value`, `test results background p-value`, `test results shared distributions p-value`)
# Set each gene to 1 if adjP < cutoff, 0, otherwise
significanceInfo <- filter(significanceInfo, significanceInfo$`test results p-value` > 0.05 & 
                             significanceInfo$`test results background p-value` < 0.05 & 
                             significanceInfo$`test results shared distributions p-value` < 0.05)

geneList <- as.numeric(as.character(significanceInfo$`test results background p-value`))

# Give geneList names:
names(geneList) <- significanceInfo$orthogroup

# Create topGOData object
GOdataBP <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSelectionFun = function(x)x,
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
resultKSBP <- runTest(GOdataBP, algorithm = "weight01", statistic = "ks")
resultKSBP
resultKSBPTable <- GenTable(GOdataBP, raw.p.value = resultKSBP, topNodes = length(resultKSBP@score), numChar = 120)
head(resultKSBPTable)

# Create topGOData object
GOdataMF <- new("topGOdata",
                ontology = "MF",
                allGenes = geneList,
                geneSelectionFun = function(x)x,
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
resultKSMF <- runTest(GOdataMF, algorithm = "weight01", statistic = "ks")
resultKSMF
resultKSMFTable <- GenTable(GOdataMF, raw.p.value = resultKSMF, topNodes = length(resultKSMF@score), numChar = 120)
head(resultKSMFTable)

# Create topGOData object
GOdataCC <- new("topGOdata",
                ontology = "CC",
                allGenes = geneList,
                geneSelectionFun = function(x)x,
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
resultKSCC <- runTest(GOdataCC, algorithm = "weight01", statistic = "ks")
resultKSCC
resultKSCCTable <- GenTable(GOdataCC, raw.p.value = resultKSCC, topNodes = length(resultKSCC@score), numChar = 120)
head(resultKSCCTable)

goEnrichmentSummaries <- capture.output(print(resultFisherBP), 
                                        print(resultFisherMF), 
                                        print(resultFisherCC), 
                                        print(resultKSBP), 
                                        print(resultKSMF), 
                                        print(resultKSCC))

writeLines(goEnrichmentSummaries, con = file(base::paste("./Results/", args[2], "/bustedPHGOSummariesBackground.csv", sep = "")))



