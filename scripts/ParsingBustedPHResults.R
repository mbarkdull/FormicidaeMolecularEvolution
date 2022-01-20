#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("At least two arguments required: the path to BUSTED-PH results, and the corresponding trait prefix.", call.=FALSE)
} else if (length(args)==1) {
  # default label:
  args[2] = "labelled"
}

# Load required packages:
library(RJSONIO)
library(tidyverse)
library(ggthemes)
library(gt)
#May need to install webshot to save gt outputs:
#webshot::install_phantomjs()

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
        data <- c(result[["input"]][["file name"]], 
                  textResult, 
                  result[["test results"]][["p-value"]], 
                  result[["test results background"]][["p-value"]], 
                  result[["test results shared distributions"]][["p-value"]], 
                  "Yes",
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["proportion"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["proportion"]])
        return(data)
      } else {
        textResult <- "Selection is associated with the phenotype / trait, but there is no significant difference between test and background branches in terms of selective pressure"
        data <- c(result[["input"]][["file name"]], 
                  textResult, 
                  result[["test results"]][["p-value"]], 
                  result[["test results background"]][["p-value"]], 
                  result[["test results shared distributions"]][["p-value"]], 
                  "No",
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["proportion"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["proportion"]])
      }
    } else {
      #print("Selection is acting on the branches with the phenotype / trait, but is **also** acting on background branches.")
      #textResult <- "Selection is acting on the branches with the phenotype / trait, but is **also** acting on background branches."
      #data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]])
      if (result[["test results shared distributions"]][["p-value"]] <= 0.05) {
        textResult <- "There is a significant difference between test and background branches in terms of selective pressure"
        data <- c(result[["input"]][["file name"]], 
                  textResult, 
                  result[["test results"]][["p-value"]], 
                  result[["test results background"]][["p-value"]], 
                  result[["test results shared distributions"]][["p-value"]], 
                  "Yes",
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["proportion"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["proportion"]])
      } else {
        textResult <- "There is no significant difference between test and background branches in terms of selective pressure"
        data <- c(result[["input"]][["file name"]], 
                  textResult, 
                  result[["test results"]][["p-value"]], 
                  result[["test results background"]][["p-value"]], 
                  result[["test results shared distributions"]][["p-value"]], 
                  "No",
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["proportion"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["omega"]],
                  result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["proportion"]])
      }
    }
  } else {
    textResult <- "There is **no evidence** of episodic diversifying selection on test branches; selection is not associated with phenotype/trait"
    data <- c(result[["input"]][["file name"]], 
              textResult, 
              result[["test results"]][["p-value"]], 
              result[["test results background"]][["p-value"]], 
              result[["test results shared distributions"]][["p-value"]], 
              "No",
              result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["omega"]],
              result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["proportion"]],
              result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["omega"]],
              result[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["proportion"]])
  }
  
}
# Make a safer version of that function with `possibly`:
possiblyparsingBustedPH <- possibly(parsingBustedPH, otherwise = "File empty.")
# Map the function over all results files to construct a master dataframe:
bustedPHResults <- purrr::map(files, possiblyparsingBustedPH)

# Write a function to process and fix the column names of the results dataframe:
bustedPHDataframeProcessing <- function(resultsDataframe) {
  BustedPHResults <- as.data.frame(do.call(rbind, resultsDataframe))   
  colnames(BustedPHResults) <- c("file", 
                                 "textResult", 
                                 "test results p-value", 
                                 "test results background p-value", 
                                 "test results shared distributions p-value", 
                                 "differenceInSelection",
                                 "unconstrainedTestOmega",
                                 "unconstrainedTestProportion",
                                 "unconstrainedBackgroundOmega",
                                 "unconstrainedBackgroundProportion")
  orthogroup <- sapply(strsplit(as.character(BustedPHResults$file),"/"), tail, 1)
  orthogroup <- sapply(strsplit(orthogroup, "_"), `[`, 2)
  BustedPHResults$orthogroup <- orthogroup
  trait <- sapply(strsplit(as.character(BustedPHResults$file),"/"), tail, 1)
  trait <- sapply(strsplit(trait, "_"), `[`, 1)
  BustedPHResults$trait <- trait
  return(BustedPHResults)
}
bustedPHResults <- bustedPHDataframeProcessing(bustedPHResults)

# Create a column that classifies orthogroups based on selective regime:
bustedPHResults <- bustedPHResults %>% mutate(selectionOn =
                     case_when(as.numeric(as.character(`test results p-value`)) <= 0.05 & 
                                 as.numeric(as.character(`test results background p-value`)) > 0.05 &
                                 as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 ~ "ForegroundOnly",
                               
                               as.numeric(as.character(`test results p-value`)) <= 0.05 & 
                                 as.numeric(as.character(`test results background p-value`)) <= 0.05 &
                                 as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 ~ "SelectionOnBothButDifferent",
                               
                               as.numeric(as.character(`test results p-value`)) <= 0.05 & 
                                 as.numeric(as.character(`test results background p-value`)) <= 0.05 &
                                 as.numeric(as.character(`test results shared distributions p-value`)) > 0.05 ~ "SelectionOnBothButNoSignificantDifference",
                               
                               as.numeric(as.character(`test results p-value`)) <= 0.05 & 
                                 as.numeric(as.character(`test results background p-value`)) > 0.05 &
                                 as.numeric(as.character(`test results shared distributions p-value`)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithTraitButNS",
                               
                               as.numeric(as.character(`test results p-value`)) > 0.05 & 
                                 as.numeric(as.character(`test results background p-value`)) <= 0.05 &
                                 as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 ~ "BackgroundOnly",
                               
                               as.numeric(as.character(`test results p-value`)) > 0.05 & 
                                 as.numeric(as.character(`test results background p-value`)) <= 0.05 &
                                 as.numeric(as.character(`test results shared distributions p-value`)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithLackOfTraitButNS",
                               
                               as.numeric(as.character(`test results p-value`)) > 0.05 & 
                                 as.numeric(as.character(`test results background p-value`)) > 0.05 &
                                 as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 ~ "NoEvidenceOfSelection",
                               
                               as.numeric(as.character(`test results p-value`)) > 0.05 & 
                                 as.numeric(as.character(`test results background p-value`)) > 0.05 &
                                 as.numeric(as.character(`test results shared distributions p-value`)) > 0.05 ~ "NoEvidenceOfSelection"))

# Convert the p-value and omega columns to numeric, not character:

bustedPHResults$`test results p-value` <- as.numeric(as.character(bustedPHResults$`test results p-value`))
bustedPHResults$`test results background p-value` <- as.numeric(as.character(bustedPHResults$`test results background p-value`))
bustedPHResults$`test results shared distributions p-value` <- as.numeric(as.character(bustedPHResults$`test results shared distributions p-value`))
bustedPHResults$unconstrainedTestOmega <- as.numeric(as.character(bustedPHResults$unconstrainedTestOmega))
bustedPHResults$unconstrainedTestProportion <- as.numeric(as.character(bustedPHResults$unconstrainedTestProportion))
bustedPHResults$unconstrainedBackgroundOmega <- as.numeric(as.character(bustedPHResults$unconstrainedBackgroundOmega))
bustedPHResults$unconstrainedBackgroundProportion <- as.numeric(as.character(bustedPHResults$unconstrainedBackgroundProportion))

# Create a nice table of some of the results, using gt:
bustedPHResults %>%
  dplyr::select(orthogroup, `test results p-value`, `test results background p-value`, `test results shared distributions p-value`) %>%
  dplyr::arrange(`test results p-value`) %>%
  dplyr::filter(`test results p-value` > 0.001) %>%
  dplyr::filter(`test results background p-value` < 0.06) %>%
  dplyr::filter(`test results background p-value` > 0.04) %>%
  dplyr::filter(`test results shared distributions p-value` < 0.06) %>%
  head() %>%
  gt::gt()

# Create an output directory:
dir.create("./Results/")
dir.create(base::paste("./Results/", args[2], sep = ""))
# Export the results as a csv:
write_csv(bustedPHResults, base::paste("./Results/", args[2], "/bustedPHResults.csv", sep = ""))

###############################################
####### Fisher's exact test ###################
###############################################

# Create a dataframe with genes only under selection in the background:
backgroundSelected <- bustedPHResults %>% 
  filter(as.numeric(as.character(`test results p-value`)) > 0.05 &
           as.numeric(as.character(`test results background p-value`)) <= 0.05 & 
           as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05)
nBackgroundSelected <- nrow(backgroundSelected)

# Create a dataframe with genes only under selection in the foreground:
foregroundSelected <- bustedPHResults %>% 
  filter(as.numeric(as.character(`test results p-value`)) <= 0.05 & 
           as.numeric(as.character(`test results background p-value`)) > 0.05 & 
           as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05)
nForegroundSelected <- nrow(foregroundSelected)

# Construct a frequency table with that information:
frequencyTableWithNoCorrection <- matrix(c(nBackgroundSelected, nForegroundSelected, 
                                           (nrow(bustedPHResults) - nBackgroundSelected),
                                           (nrow(bustedPHResults) - nForegroundSelected)),
                                         nrow = 2,
                                         byrow = TRUE,
                                         dimnames = list("category" = c("selected", "noSelection") ,
                                                         "selection" = c("backgroundSelected", "foregroundSelected")))
# Run a fisher's exact test to see if there is a difference in the proportion of genes under selection in the fore- vs. background:
fishersExactTest <- fisher.test(frequencyTableWithNoCorrection)
# Print the resulting p-value:
fishersExactTest[["p.value"]]

# Print text explaining the result:
if (fishersExactTest[["p.value"]] <= 0.05) {
  print("Signficant difference between foreground and background") 
} else {
  print("No signficant difference.")
}

# Export the result of the fisher's exact test:
fishersExactTestPrint <- capture.output(print(fishersExactTest))
writeLines(fishersExactTestPrint, con = file(base::paste("./Results/", args[2], "/bustedPHFisher.csv", sep = "")))

###############################################
####### Make some nice plots ##################
###############################################

print(fishersExactTest[["p.value"]])

pValue <- if (fishersExactTest[["p.value"]] < 0.000001) {
  formatC(fishersExactTest[["p.value"]], format = "e", digits = 2)
} else {
  round(fishersExactTest[["p.value"]], digits = 3)
}
print(pValue)
trait <- args[2]

bustedPHResults$selectionOn <-
  factor(bustedPHResults$selectionOn,
         levels = c("ForegroundOnly",
                    "BackgroundOnly",
                    "EvidenceOfSelectionAssociatedWithTraitButNS",
                    "EvidenceOfSelectionAssociatedWithLackOfTraitButNS",
                    "SelectionOnBothButDifferent",
                    "SelectionOnBothButNoSignificantDifference",
                    "NoEvidenceOfSelection"))

plot <- ggplot(data = dplyr::filter(bustedPHResults, !is.na(bustedPHResults$selectionOn)))

plot + 
  geom_bar(mapping = aes(x = selectionOn)) + 
  geom_linerange(aes(xmin = 1, 
                     xmax = 2, 
                     y = (1.75 * length(which(selectionOn == "ForegroundOnly")))), 
                 color = "grey26", 
                 size = 0.3) +
  geom_text(aes(x = 1.5, 
                y = (2.25 * length(which(selectionOn == "ForegroundOnly"))),
                label = paste("p-value = ", pValue, sep = "")),
            stat = "unique") +
  theme_bw() +
  labs(x = "Selective regime", 
       y = "Count of orthogroups", 
       title = paste("Distribution of positive selection regimes \nassociated with", trait, sep = " ")) +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) + 
  scale_x_discrete(labels=c("ForegroundOnly" = "Selection on\nforeground only",
                            "BackgroundOnly" = "Selection on\nbackground only",
                            "EvidenceOfSelectionAssociatedWithTraitButNS" = "Nonsignificant\nevidence of\nselection associated\nwith trait",
                            "EvidenceOfSelectionAssociatedWithLackOfTraitButNS" = "Nonsignificant\nevidence of\nselection associated\nwith lacking trait",
                            "SelectionOnBothButNoSignificantDifference" = "Selection on both fore-\nand background,\nno signficant difference",
                            "SelectionOnBothButDifferent" = "Selection on\nfore- and background\nwith significant differences",
                            "NoEvidenceOfSelection" = "No evidence\nfor selection"))

ggsave(filename = "barPlotSelectionBustedPH.png", 
       path = base::paste("./Results/", args[2], "/", sep = ""), 
       width = 7,
       height = 4,
       units = "in",
       dpi = 600)

plot + 
  geom_point(mapping = aes(x = `test results p-value`,
                           y = `test results background p-value`,
                           color = `test results shared distributions p-value`)) +
  geom_vline(xintercept = 0.05) +
  geom_hline(yintercept = 0.05)

ggplot(data = filter(bustedPHResults, 
                     `test results shared distributions p-value` <= 0.05), 
       mapping = aes(x = `test results p-value`,
                     y = `test results background p-value`)) +
  geom_hex(bins = 50) +
  scale_fill_continuous(trans = "log10") +
  geom_vline(xintercept = 0.05) +
  geom_hline(yintercept = 0.05) + 
  labs(x = "Significance (p-value) for evidence\nof selection on species with the trait",
       y = "Significance (p-value) for evidence\nof selection on species without the trait",
       title = "Distribution of selective regimes on\ngenes with a difference between species with and\nwithout the trait")

maxOmega <- max(c(bustedPHResults$unconstrainedBackgroundOmega, bustedPHResults$unconstrainedTestOmega))
maxProportion <- max(c(bustedPHResults$unconstrainedBackgroundProportion, bustedPHResults$unconstrainedTestProportion))

p1 <- ggplot(filter(bustedPHResults,
              `test results p-value` <= 0.05 &
                `test results shared distributions p-value` <= 0.05)) +
  geom_hex(mapping = aes(x = unconstrainedTestOmega,
                         y = unconstrainedTestProportion),
           bins = 50) +
  labs(x = "Selection strength (omega values)",
       y = "Proportion of sites in gene under selection",
       title = "Selective regimes on genes \nunder selection in species \nwith the trait")

p2 <- ggplot(filter(bustedPHResults,
              `test results background p-value` <= 0.05 &
                `test results shared distributions p-value` <= 0.05)) +
  geom_hex(mapping = aes(x = unconstrainedBackgroundOmega,
                         y = unconstrainedBackgroundProportion),
           bins = 50) +
  labs(x = "Selective strength (omega values)",
       y = "Proportion of sites in gene under selection",
       title = "Selection regimes on genes \nunder selection in species \nwithout the trait")


ggpubr::ggarrange(p1, 
                  p2,
                  ncol = 2, 
                  nrow = 1) 
rm(p1, p2)
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
significanceInfo <- dplyr::select(bustedPHResults, orthogroup, `test results p-value`, `test results background p-value`, `test results shared distributions p-value`)
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


# Combine all of the results for GO terms enriched in the foreground:
goTermsForeground <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
  dplyr::filter(raw.p.value <= 0.01)
# Generate a nice table:
goTermsForegroundTable <- goTermsForeground %>% 
  gt() %>%
  tab_header(title = "GO terms enriched for positive selection in foreground species") %>%
  cols_label(GO.ID = "Go term ID",
             Term = "GO Term",
             Annotated = "Orthogroups annotated",
             Significant = "Orthogroups under positive selection",
             Expected = "Expected number of orthogroups",
             raw.p.value = "p-value") %>%
  cols_move_to_start(columns = c(Term)) %>%
  cols_align(align = c("left"),
    columns = everything())
gtsave(goTermsForegroundTable, filename = "goTermsForeground.png", path = base::paste("./Results/", args[2], "/", sep = ""))

goEnrichmentSummaries <- capture.output(print(resultFisherBP), 
                                        print(resultFisherMF), 
                                        print(resultFisherCC))

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
significanceInfo <- dplyr::select(bustedPHResults, orthogroup, `test results p-value`, `test results background p-value`, `test results shared distributions p-value`)
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

# Combine results of all the tests:
goTermsBackground <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
  dplyr::filter(raw.p.value <= 0.01)
# Generate a nice table:
goTermsBackgroundTable <- goTermsBackground %>% 
  gt() %>%
  tab_header(title = "GO terms enriched for positive selection in background species") %>%
  cols_label(GO.ID = "Go term ID",
             Term = "GO Term",
             Annotated = "Orthogroups annotated",
             Significant = "Orthogroups under positive selection",
             Expected = "Expected number of orthogroups",
             raw.p.value = "p-value") %>%
  cols_move_to_start(columns = c(Term)) %>%
  cols_align(align = c("left"),
             columns = everything())
gtsave(goTermsBackgroundTable, filename = "goTermsBackground.png", path = base::paste("./Results/", args[2], "/", sep = ""))

goEnrichmentSummaries <- capture.output(print(resultFisherBP), 
                                        print(resultFisherMF), 
                                        print(resultFisherCC))

writeLines(goEnrichmentSummaries, con = file(base::paste("./Results/", args[2], "/bustedPHGOSummariesBackground.csv", sep = "")))

# Look at particular, interesting genes:
# krf1 <- filter(bustedPHResults, orthogroup == "OG0004124")

