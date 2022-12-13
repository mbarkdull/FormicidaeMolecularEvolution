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
#jsonFiles <- list.files(path = "./10_1_RelaxResults/polygyny", pattern = "*.json", full.names = TRUE)
#jsonFiles <- list.files(path = "./10_1_RelaxResults/workerPolymorphism", pattern = "*.json", full.names = TRUE)
#jsonFiles <- list.files(path = "./10_1_RelaxResults/workerReproductionQueens", pattern = "*.json", full.names = TRUE)
#jsonFiles <- list.files(path = "./10_1_RelaxResults/multilineage", pattern = "*.json", full.names = TRUE)
#jsonFiles <- list.files(path = "./10_1_RelaxResults/polyandry", pattern = "*.json", full.names = TRUE)

jsonFiles <- sort(jsonFiles, decreasing = TRUE)
jsonFiles <- jsonFiles[sapply(jsonFiles, file.size) > 0]

relaxJSONProcessing <- function(i) {
  relaxResult <- rjson::fromJSON(file = i)
  if (relaxResult[["test results"]][["p-value"]] <= 0.05) {
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
colnames(relaxResults) <- c("fileName", "pValue", "kValue", "description", "shortDescription")
orthogroup <- sapply(strsplit(as.character(relaxResults$fileName),"/"), tail, 1)
orthogroup <- sapply(strsplit(orthogroup, "_"), `[`, 2)
relaxResults$orthogroup <- orthogroup
relaxResults$pValue <- as.numeric(as.character(relaxResults$pValue))
relaxResults$kValue <- as.numeric(as.character(relaxResults$kValue))

numberRelax <- sum(relaxResults$shortDescription == "Relaxation of selection along foreground branches")
numberNonsignficant <- sum(relaxResults$description == "No significant difference in selective regime between foreground and background branches")
numberIntensified <- sum(relaxResults$shortDescription == "Intensification of selection along foreground branches")
percentRelaxed <- (numberRelax / (numberRelax + numberNonsignficant + numberIntensified))*100
percentIntensified <- (numberIntensified / (numberRelax + numberNonsignficant + numberIntensified))*100

relaxResults %>%
  dplyr::select(orthogroup, kValue, pValue) %>%
  dplyr::filter(pValue < 0.055) %>%
  dplyr::filter(pValue > 0.045) %>%
  head() %>%
  gt::gt()

# Create an output directory:
dir.create("./Results")
outputDirectory <- paste("./Results/", args[2], sep = "")
dir.create(outputDirectory)
outputFile <- paste("./Results/", args[2], "/RelaxResults.csv", sep = "")
print(outputFile)
# Export the results:
write_csv(relaxResults, outputFile)

###############################################
####### Fisher's exact test ###################
###############################################
  # "The number of relaxed orthogroups in Coregonus and Salvelinus was significantly higher compared to the number of intensified orthogroups (one-sided Fisher’s Exact Test, p =0.035)." From Schneider, Kevin, Colin E. Adams, and Kathryn R. Elmer. "Parallel selection on ecologically relevant gene functions in the transcriptomes of highly diversifying salmonids." BMC genomics 20.1 (2019): 1-23. 	

### Kevin Schneider; last modified: 07.10.2019 ###
### KevinSchneider@gmx.at ###

library(ggplot2)

relaxframe <- relaxResults
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

fishersExactTest <- fisher.test(frequencyTableWithNoCorrection)
fisher.test(frequencyTableWithNoCorrection)
fisher.test(frequencyTableWithFDR)
fisher.test(frequencyTableWithBonferroni)
fisher.test(freq_table_below_or_above_1)

freq_table_selparam <- matrix(c((length(relaxframe$kValue[relaxframe$kValue != 1]) / 2), 
                                (length(relaxframe$kValue[relaxframe$kValue != 1]) / 2), 
                                length(relaxframe$kValue[relaxframe$kValue < 1]),
                                length(relaxframe$kValue[relaxframe$kValue > 1])),
                              nrow = 2,
                              byrow = TRUE,
                              dimnames = list("category" = c("expected", "observed") ,
                                              "selection" = c("relaxed", "intensified")))

fisher.test(freq_table_selparam)
fishersExactTestPrint <- capture.output(print(fishersExactTest))

writeLines(fishersExactTestPrint, con = file(base::paste("./Results/", args[2], "/RelaxFisher.csv", sep = "")))

###############################################
####### Make some nice plots ##################
###############################################
print(fishersExactTest[["p.value"]])

pValue <- if (fishersExactTest[["p.value"]] < 0.000001) {
  formatC(fishersExactTest[["p.value"]], format = "e", digits = 2)
} else {
  round(fishersExactTest[["p.value"]], digits = 4)
}
print(pValue)
pValueLabel <- paste("p-value = ", pValue, sep = "")
trait <- args[2]

relaxResults$shortDescription <-
  factor(relaxResults$shortDescription,
         levels = c("Intensification of selection along foreground branches", 
                    "Relaxation of selection along foreground branches",
                    "Nonsignificant intensification",
                    "Nonsignificant relaxation",
                    "File empty."))
relaxResults$kValue <- as.numeric(as.character(relaxResults$kValue))

# Create a base plot, filtering out relax files that are empty:
plot <- ggplot(data = filter(relaxResults, 
                             shortDescription != "File empty."))

# Combine that base plot with a bar plot to show the number of genes intensified, relaxed, and NS:
plot + 
  geom_bar(mapping = aes(x = shortDescription)) +
  labs(x = "Selective regime", 
       y = "Count of orthogroups", 
       title = "Distribution of selective regimes") +
  theme(axis.text.x = element_text(angle = 0,
                                 hjust = 0.5,
                                 vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) + 
  scale_x_discrete(labels=c("Nonsignificant relaxation"  = "Nonsignificant\nrelaxation" ,
                            "Intensification of selection along foreground branches" = "Intensification of\nselection along\nforeground branches",
                            "Nonsignificant intensification" = "Nonsignificant\nintensification",
                            "Relaxation of selection along foreground branches" = "Relaxation of\nselection along\nforeground branches",
                            "File empty." = "File empty.")) + 
  geom_linerange(aes(xmin = 1, 
                     xmax = 2, 
                     y = (1.75 * length(which(shortDescription == "Intensification of selection along foreground branches")))), 
                 color = "grey26", 
                 size = 0.3) +
  geom_text(aes(x = 1.5, 
                y = (1.9 * length(which(shortDescription == "Intensification of selection along foreground branches"))),
                label = pValueLabel),
            stat = "unique") + 
  theme_bw() 

# Create a plot that shows the distribution of k-values:
ggplot(data = filter(relaxResults, 
                     shortDescription != "File empty.", 
                     kValue < 2)) + 
  geom_histogram(mapping = aes(x = kValue),
                 binwidth = 0.05,
                 center = 0,
                 color = "grey26",
                 alpha = 0) + 
  geom_vline(xintercept = 1,
             color = "darkred",
             size = 0.75) +
  geom_text(aes(x = 1.05, y = 175,
                label = "Neutral expectation, k = 1"),
            stat = "unique",
            hjust = 0, 
            color = "darkred") +
  theme_bw() +
  labs(x = "Relaxation parameter, k", 
       y = "Count of orthogroups", 
       title = paste("Distribution of selective regimes associated with", trait, sep = " ")) +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

relaxResults$pValue <- as.numeric(relaxResults$pValue)
relaxResults$kValue <- as.numeric(relaxResults$kValue)

ggplot(data = relaxResults) +
  geom_hex(mapping = aes(x = pValue,
                         y = kValue),
           bins = 50) +
  geom_vline(xintercept = 0.05)

###############################################
####### Check for GO term enrichment ##########
###############################################

# Load packages:
library(topGO)
library(tidyverse)
library(RJSONIO)
library(gt)

print("####################################################################")
print("###### GO term enrichment on genes under relaxed selection #####")
print("####################################################################")

# Define a function to generate violin plots:
violinPlotting <- function(fisherResultsTable, inferredGeneAnnotations, plotTitle) {
  violinData <- full_join(fisherResultsTable, inferredGeneAnnotations, by = c("GO.ID" = "goTerm"))
  violinData <- full_join(violinData, relaxResults, by = c("orthogroup" = "orthogroup"))
  
  ggplot(data = filter(violinData, 
                       raw.p.value < 0.01),
         mapping = aes(x = Term,
                       y = as.numeric(as.character(kValue)))) +
    geom_violin() +
    geom_jitter(height = 0, 
                width = 0.1) +
    stat_summary(fun = median, 
                 geom = "crossbar", 
                 size = 0.2, 
                 color = "darkred") +
    geom_hline(yintercept = 1) +
    theme_bw() +
    labs(title = plotTitle,
         x = "Significantly enriched GO terms",
         y = "Distribution of relaxtion parameters (k values)") +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     vjust = 1),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())
  
}

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
significanceInfo <- dplyr::select(relaxResults, orthogroup, pValue, kValue) 
# Set each gene to 1 if adjP < cutoff, 0, otherwise
pcutoff <- 0.05 
tmp <- ifelse(significanceInfo$pValue < pcutoff & significanceInfo$kValue < 1, 1, 0)
geneList <- tmp
rm(tmp)

# Give geneList names:
names(geneList) <- significanceInfo$orthogroup

# Create the GOdata object:
GOdataBP <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
# Run Fisher's exact test to check for enrichment:
resultsFisherBP <- runTest(GOdataBP, algorithm = "elim", statistic = "fisher")
resultsFisherBP
resultsFisherBPTable <- GenTable(GOdataBP, raw.p.value = resultsFisherBP, topNodes = length(resultsFisherBP@score),
                                 numChar = 120)
head(resultsFisherBPTable)

# topGO infers additional annotations not present in your input. We want to map these back to their orthogroups. 
inferredGeneAnnotationsBP <- genesInTerm(GOdataBP)
inferredGeneAnnotationsBP <- map_dfr(.x = inferredGeneAnnotationsBP,
        .f = ~ enframe(x = .x,
                       name = NULL,
                       value = "orthogroup"),
        .id = "goTerm")


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
inferredGeneAnnotationsMF <- genesInTerm(GOdataMF)
inferredGeneAnnotationsMF <- map_dfr(.x = inferredGeneAnnotationsMF,
                                     .f = ~ enframe(x = .x,
                                                    name = NULL,
                                                    value = "orthogroup"),
                                     .id = "goTerm")


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
inferredGeneAnnotationsCC <- genesInTerm(GOdataCC)
inferredGeneAnnotationsCC <- map_dfr(.x = inferredGeneAnnotationsCC,
                                     .f = ~ enframe(x = .x,
                                                    name = NULL,
                                                    value = "orthogroup"),
                                     .id = "goTerm")


allResultsFisherRelaxed <- bind_rows(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable)
# Combine all of the results for GO terms enriched in the foreground:
goTermsRelaxed <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
  dplyr::filter(raw.p.value <= 0.01)
# Generate a nice table:
goTermsRelaxedTable <- goTermsRelaxed %>% 
  gt() %>%
  tab_header(title = "GO terms enriched for relaxed selection") %>%
  cols_label(GO.ID = "Go term ID",
             Term = "GO Term",
             Annotated = "Orthogroups annotated",
             Significant = "Orthogroups under relaxed selection",
             Expected = "Expected number of orthogroups under relaxed selection",
             raw.p.value = "p-value") %>%
  cols_move_to_start(columns = c(Term)) %>%
  cols_align(align = c("left"),
             columns = everything())
gtsave(goTermsRelaxedTable, filename = "goTermsRelaxed.png", path = base::paste("./Results/", args[2], "/", sep = ""))

allInferredAnnotationsRelaxed <- bind_rows(inferredGeneAnnotationsBP, inferredGeneAnnotationsMF, inferredGeneAnnotationsCC)

if (length(which(allResultsFisherRelaxed$raw.p.value <= 0.01)) > 0) {
  violinPlotting(allResultsFisherRelaxed, allInferredAnnotationsRelaxed, plotTitle = "Relaxation parameter distributions across significantly enriched GO terms")
  ggsave(filename = paste(args[2], "violinRelaxedSelection.png", sep = ""), path = base::paste("./Results/", args[2], "/", sep = ""))
} else {
  print("no significant GO terms")
}

goEnrichmentSummaries <- capture.output(print(resultsFisherBP), 
                                        print(resultFisherMF), 
                                        print(resultFisherCC))

writeLines(goEnrichmentSummaries, con = file(base::paste("./Results/", args[2], "/", args[2], "RelaxGOSummariesRelaxed.csv", sep = "")))

############# GO Term enrichment for genes under intensified selection #################

print("####################################################################")
print("###### GO term enrichment on genes under intensified selection #####")
print("####################################################################")

# Define vector that is 1 if gene is significantly DE (`test results p-value` < 0.05) and 0 otherwise:
significanceInfo <- dplyr::select(relaxResults, orthogroup, pValue, kValue) 
# Set each gene to 1 if adjP < cutoff, 0, otherwise
pcutoff <- 0.05 
tmp <- ifelse(significanceInfo$pValue < pcutoff & significanceInfo$kValue > 1, 1, 0)
geneList <- tmp
rm(tmp)

# Give geneList names:
names(geneList) <- significanceInfo$orthogroup

# Create the GOdata object:
GOdataBP <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
# Run Fisher's exact test to check for enrichment:
resultsFisherBP <- runTest(GOdataBP, algorithm = "elim", statistic = "fisher")
resultsFisherBP
resultsFisherBPTable <- GenTable(GOdataBP, raw.p.value = resultsFisherBP, topNodes = length(resultsFisherBP@score),
                                 numChar = 120)
head(resultsFisherBPTable)
inferredGeneAnnotationsBP <- genesInTerm(GOdataBP)
inferredGeneAnnotationsBP <- map_dfr(.x = inferredGeneAnnotationsBP,
                                     .f = ~ enframe(x = .x,
                                                    name = NULL,
                                                    value = "orthogroup"),
                                     .id = "goTerm")

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
inferredGeneAnnotationsMF <- genesInTerm(GOdataMF)
inferredGeneAnnotationsMF <- map_dfr(.x = inferredGeneAnnotationsMF,
                                     .f = ~ enframe(x = .x,
                                                    name = NULL,
                                                    value = "orthogroup"),
                                     .id = "goTerm")

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
inferredGeneAnnotationsCC <- genesInTerm(GOdataCC)
inferredGeneAnnotationsCC <- map_dfr(.x = inferredGeneAnnotationsCC,
                                     .f = ~ enframe(x = .x,
                                                    name = NULL,
                                                    value = "orthogroup"),
                                     .id = "goTerm")

allResultsFisherIntensified <- bind_rows(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable)
goTermsBackground <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
  dplyr::filter(raw.p.value <= 0.01)
# Generate a nice table:
goTermsBackgroundTable <- goTermsBackground %>% 
  gt() %>%
  tab_header(title = "GO terms enriched for intensified selection") %>%
  cols_label(GO.ID = "Go term ID",
             Term = "GO Term",
             Annotated = "Orthogroups annotated",
             Significant = "Orthogroups under intensified selection",
             Expected = "Expected number of orthogroups under intensified selection",
             raw.p.value = "p-value") %>%
  cols_move_to_start(columns = c(Term)) %>%
  cols_align(align = c("left"),
             columns = everything())
gtsave(goTermsBackgroundTable, filename = "goTermsIntensified.png", path = base::paste("./Results/", args[2], "/", sep = ""))

allInferredAnnotationsIntensified <- bind_rows(inferredGeneAnnotationsBP, inferredGeneAnnotationsMF, inferredGeneAnnotationsCC)

if (length(which(allResultsFisherIntensified$raw.p.value <= 0.01)) > 0) {
  violinPlotting(allResultsFisherIntensified, allInferredAnnotationsIntensified, plotTitle = "Relaxation parameter distributions across significantly enriched GO terms")
  ggsave(filename = paste(args[2], "violinIntensifiedSelection.png", sep = ""), path = base::paste("./Results/", args[2], "/", sep = ""))
} else {
  print("no significant GO terms")
}



##### FIX!!!!!!! #####
# Use the Kolomogorov-Smirnov test:
#significanceInfo <- dplyr::select(relaxResults, orthogroup, pValue, kValue) %>%
  #dplyr::filter(kValue > 1)
#geneList <- as.numeric(as.character(significanceInfo$pValue))
#names(geneList) <- significanceInfo$orthogroup
# Create topGOData object
#GOdataBP <- new("topGOdata",
                #ontology = "BP",
                #allGenes = geneList,
                #geneSelectionFun = function(x)x,
                #annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
#resultKSBP <- runTest(GOdataBP, algorithm = "weight01", statistic = "ks")
#resultKSBP
#resultKSBPTable <- GenTable(GOdataBP, raw.p.value = resultKSBP, topNodes = length(resultKSBP@score), numChar = 120)
#head(resultKSBPTable)

# Create topGOData object
#GOdataMF <- new("topGOdata",
                #ontology = "MF",
                #allGenes = geneList,
                #geneSelectionFun = function(x)x,
                #annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
#resultKSMF <- runTest(GOdataMF, algorithm = "weight01", statistic = "ks")
#resultKSMF
#resultKSMFTable <- GenTable(GOdataMF, raw.p.value = resultKSMF, topNodes = length(resultKSMF@score), numChar = 120)
#head(resultKSMFTable)

# Create topGOData object
#GOdataCC <- new("topGOdata",
                #ontology = "CC",
                #allGenes = geneList,
                #geneSelectionFun = function(x)x,
                #annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
#resultKSCC <- runTest(GOdataCC, algorithm = "weight01", statistic = "ks")
#resultKSCC
#resultKSCCTable <- GenTable(GOdataCC, raw.p.value = resultKSCC, topNodes = length(resultKSCC@score), numChar = 120)
#head(resultKSCCTable)

goEnrichmentSummaries <- capture.output(print(resultsFisherBP), 
                                        print(resultFisherMF), 
                                        print(resultFisherCC))

writeLines(goEnrichmentSummaries, con = file(base::paste("./Results/", args[2], "/RelaxGOSummariesIntensified.csv", sep = "")))

# Look at particularly interesting genes:
krf1 <- filter(relaxResults, orthogroup == "OG0004124")
