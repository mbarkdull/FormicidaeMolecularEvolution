#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Must provide the path to BUSTED results as a command line argument.)", call.=FALSE)
}

###############################################
####### Convert JSON files to spreadsheet #####
###############################################

# This script reads in the results JSON files produced by BUSTED, and creates a list of orthogroups with evidence for positive selection somewhere in the tree. 
library(plyr)
library(tidyverse)
library(rjson)

# Check how many orthogroups should have been tested. Some orthogroups contain only a single gene, and so weren't tested by BUSTED[S]. 
orthogroupMembers <- readr::read_table2("./5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Orthogroups/Orthogroups.txt", col_names = FALSE) 
singletons <- subset(orthogroupMembers, is.na(orthogroupMembers$X3)) 
nonCobsSingletons <- singletons %>% 
  filter(!str_detect(X2, 'cobs')) %>%
  distinct()
multiGeneGroups <- orthogroupMembers %>%
  drop_na(X3)
multiGeneGroups$X1 <- gsub(':', '', multiGeneGroups$X1) 
multiGeneGroups <- multiGeneGroups%>%
  dplyr::select(c(X1))
  
# Construct a list of all of the json files:
jsonFiles <- list.files(path = args[1], pattern = "*.json", full.names = TRUE)
jsonFiles <- sort(jsonFiles, decreasing = TRUE)

#jsonFiles <- list.files(path = "./8_3_BustedResults/", pattern = "*.json", full.names = TRUE)
# Get a list of all orthogroups tested:
getOrthogroup <- function(i) {
  orthogoupName <- sapply(strsplit(as.character(i),"/"), tail, 1)
  orthogoupName <- sapply(strsplit(orthogoupName, "\\_"), `[`, 1)
}
possiblygetOrthogroup <- possibly(getOrthogroup, otherwise = "Error.")
orthogroupsFiles <- list.files(path = "5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Orthogroup_Sequences/", full.names = TRUE)
allOrthogroups <- map(jsonFiles, possiblygetOrthogroup)
testedOrthogroups <- map(jsonFiles, possiblygetOrthogroup)
testedOrthogroups <- as.data.frame(do.call(rbind, testedOrthogroups))   
# See which orthogroups weren't tested:
multiGeneGroups$tested <- multiGeneGroups$X1 %in% testedOrthogroups$V1  

# bustedResult <- rjson::fromJSON(file = "./8_3_BustedResults/OG0000067_busted.json")
# Write a function that will process each individual json file and extract the file name, orthogroup number, p-value, and return a text description of the p-value:
bustedJSONProcessing <- function(i) {
  bustedResults <- rjson::fromJSON(file = i)
  # Now run my if else statement:
  # If the p value is less than 0.05, that means there is positive selection. 

    orthogoupName <- sapply(strsplit(as.character(i),"/"), tail, 1)
    orthogoupName <- sapply(strsplit(orthogoupName, "\\_"), `[`, 1)
    
    # Construct a vector of data containing the file name, the orthogroup number, the p-value, and the text "yes, evidence for positive selection":
    data <- c(bustedResults[["input"]][["file name"]], 
              orthogoupName, 
              bustedResults[["test results"]][["p-value"]], 
              bustedResults[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["omega"]],
              bustedResults[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["proportion"]])
    return(data)
}

# Create a version of the function that returns an error if there's an empty file (from https://www.r-bloggers.com/2017/12/skip-errors-in-r-loops-by-not-writing-loops/):
possiblyBustedJSONProcessing <- possibly(bustedJSONProcessing, otherwise = "File empty.")
# Run this function with purrr:map so as to avoid for loops (from https://www.r-bloggers.com/2017/12/skip-errors-in-r-loops-by-not-writing-loops/ and https://jennybc.github.io/purrr-tutorial/ls01_map-name-position-shortcuts.html):
bustedResults <- map(jsonFiles, possiblyBustedJSONProcessing)

# Convert the results to a dataframe:
bustedResults <- as.data.frame(do.call(rbind, bustedResults))   
colnames(bustedResults) <- c("file", "orthogroup", "test results p-value", "omega3", "proportionOfSitesUnderSelection")
bustedResults$`test results p-value` <- as.numeric(as.character(bustedResults$`test results p-value`))
bustedResults$omega3 <- as.numeric(as.character(bustedResults$omega3))
bustedResults$proportionOfSitesUnderSelection <- as.numeric(as.character(bustedResults$proportionOfSitesUnderSelection))

# Adjust the raw p-values with Benjamini-Hochberg FDR:
bustedResults <- bustedResults %>%
  mutate(pValueFDR = p.adjust(`test results p-value`, method='BH')) 

# Get the number of genes with and without positive selection:
numberPositive <- length(which(bustedResults$pValueFDR <= 0.05))
numberNone <- length(which(bustedResults$pValueFDR > 0.05))
percentPositive <- (numberPositive / (numberPositive + numberNone))*100

# Create an output directory:
dir.create("./Results")
# Export the results:
write_csv(bustedResults, "./Results/bustedResults.csv")

#### Check result for just single-copy orthogroups: ####
singleCopyOrthogroups <- read_csv("./Results/singleCopyOrthogroups.csv")
singleCopyOrthogroups <- singleCopyOrthogroups$x

singleCopyResults <- dplyr::filter(bustedResults, 
                                   orthogroup %in% singleCopyOrthogroups)
# Get the number of genes with and without positive selection:
numberPositive <- length(which(singleCopyResults$pValueFDR <= 0.05))
numberNone <- length(which(singleCopyResults$pValueFDR > 0.05))
percentPositive <- (numberPositive / (numberPositive + numberNone))*100


#### Make a heat map of p-values vs. omega3 values, a la https://github.com/veg/hyphy/issues/737 ###

d <- ggplot(filter(bustedResults, `test results p-value` <= 0.05), 
            mapping = aes(x = proportionOfSitesUnderSelection, 
                          y = omega3))

d + 
  geom_hex() + 
  scale_fill_gradient(trans = "log10")

d + 
  geom_hex(bins = 40) + 
  scale_fill_gradient(trans = "log10",
                      low = "#ccf0ed",
                      high = "#014a44") +
  scale_x_log10() + 
  scale_y_log10() +
  labs(x = "Proportion of sites in the gene that are under selection",
       y = "Strength of selection (omega)",
       fill = "Count of genes") +
  theme_bw()

###############################################
####### Check for GO term enrichment ##########
###############################################

# Load packages:
library(topGO)
library(tidyverse)
library(RJSONIO)

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
significanceInfo <- dplyr::select(bustedResults, orthogroup, pValueFDR)
# Set each gene to 1 if adjP < cutoff, 0, otherwise
pcutoff <- 0.05 
tmp <- ifelse(significanceInfo$pValueFDR < pcutoff , 1, 0)
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
significantresultsFisherBPTable <- dplyr::filter(resultsFisherBPTable, 
                                          as.numeric(raw.p.value) <= 0.01)

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
geneList <- as.numeric(as.character(bustedResults$pValueFDR))
names(geneList) <- bustedResults$orthogroup
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

resultKSBPTable %>%
  head(n = 7) %>%
  dplyr::select(GO.ID, Term, raw.p.value) %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Biological Process GO terms enriched across ants")

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

writeLines(goEnrichmentSummaries, con = file("./Results/bustedGOSummaries.csv"))
allResults <- rbind(resultsFisherBPTable,
                    resultsFisherMFTable,
                    resultsFisherCCTable,
                    resultKSBPTable,
                    resultKSMFTable,
                    resultKSCCTable)
significantGOResults <- dplyr::filter(allResults,
                                      as.numeric(raw.p.value) <= 0.01) %>%
  distinct(Term, .keep_all = TRUE)
significantGOResults$raw.p.value <- as.numeric(significantGOResults$raw.p.value)
significantGOResults <- arrange(significantGOResults, 
                                dplyr::desc(raw.p.value))
goTermsTable <- significantGOResults %>% 
  dplyr::select(-c(Expected)) %>%
  dplyr::select(c(`p-value` = raw.p.value,
                  `GO ID` = GO.ID,
                  `Total orthogroups annotated` = Annotated,
                  `Significant orthogroups` = Significant,
                  Term = Term)) %>%
  gt(rowname_col = "Term") %>%
  tab_header(title = paste("Gene ontology terms under\npositive selection across the ant phylogeny")) %>% 
  tab_options(table.width = pct(100),
              container.width = pct(100),
              container.height = pct(100),
              row_group.background.color = "#f2f2fa",
              heading.title.font.size = 30,
              heading.title.font.weight = "bolder") %>%
  cols_width(everything() ~ pct(25))

# Get a dataframe wtih just significant Go terms and the columns you want:
printResults <- significantGOResults %>% 
  dplyr::select(-c(Expected)) %>%
  dplyr::select(c(Term = Term,
                  `GO ID` = GO.ID,
                  `p-value` = raw.p.value,
                  `Total orthogroups annotated` = Annotated,
                  `Significant orthogroups` = Significant)) %>%
  arrange(`p-value`)

library(flextable)

ft <- flextable(printResults)

ft <- theme_vanilla(ft)

ft <- set_caption(ft, caption = "Gene ontology terms under\npositive selection across the ant phylogeny")

save_as_docx(ft,
             path = "./Plots/significantGOTermsBustedS.docx")

#%>%
 # gtsave(paste("GOTerms_BUSTED.png", 
  #             sep = ""), 
   #      path = "./Plots/",
    #     vwidth = 1000)








################## Generate some violin plots ###############################


#fullData <- inner_join(resultKSBPTable, GOannotations, by = c("GO.ID" = "domain_id"))
#fullData <- inner_join(fullData, bustedResults, by = c("#cluster_id" = "orthogroup"))

#ggplot(data = filter(fullData, raw.p.value < 0.01)) +
  #geom_violin(mapping = aes(x = Term,
                            #y = `test results p-value`))


