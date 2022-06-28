# Load packages:
library(topGO)
library(tidyverse)
library(RJSONIO)
library(gt)
library(cowplot)

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

# Read in BUSTED-PH results:
bustedPHResults <- read_csv("./Results/allBustedPHResults.csv", col_names = TRUE)

# For each trait, check for GO term enrichment:
goEnrichmentBUSTEDPH <- function(specificTrait) {
  significanceInfo <- dplyr::select(bustedPHResults, 
                                    orthogroup, 
                                    `test results p-value`, 
                                    `test results background p-value`, 
                                    `test results shared distributions p-value`, 
                                    trait) 
  significanceInfo <- dplyr::filter(significanceInfo, trait == specificTrait)
  pcutoff <- 0.05 
  tmp <- ifelse(significanceInfo$`test results p-value` < pcutoff & 
                  significanceInfo$`test results background p-value` > pcutoff & 
                  significanceInfo$`test results shared distributions p-value` < pcutoff, 1, 0)
  geneList <- tmp
  names(geneList) <- significanceInfo$orthogroup
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

  goTermsForeground <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
    dplyr::filter(raw.p.value <= 0.01)
  goTermsForeground$selectionCategory <- "Positive selection in species with the trait"
  
  head(goTermsForeground)
  
  # Check GO term enrichment of genes under positive selection in background species:
  significanceInfo <- dplyr::select(bustedPHResults, 
                                    orthogroup, 
                                    `test results p-value`, 
                                    `test results background p-value`, 
                                    `test results shared distributions p-value`,
                                    trait)
  significanceInfo <- dplyr::filter(significanceInfo, trait == specificTrait)
  
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

  # Combine results of all the tests:
  goTermsBackground <- rbind(resultsFisherBPTable, 
                             resultsFisherMFTable, 
                             resultsFisherCCTable) %>%
    dplyr::filter(raw.p.value <= 0.01)
  goTermsBackground$selectionCategory <- "Positive selection in species lacking the trait"
  allGOTerms <- rbind(goTermsForeground,
                      goTermsBackground)
  allGOTerms$trait <- specificTrait
  allGOTerms
}

# Get a vector of traits:
traits <- unique(bustedPHResults$trait)

# Map the GO enrichment function over the traits:
goTermsBustedPH <- purrr::map(traits, goEnrichmentBUSTEDPH)

# Read in the RELAX results:
relaxResults <- read_csv("./Results/allRelaxResults.csv", col_names = TRUE)

# For each trait, check for GO term enrichment:
goEnrichmentRELAX <- function(specificTrait) {
  # Define vector that is 1 if gene is under relaxed selection and 0 otherwise:
  significanceInfo <- dplyr::select(relaxResults, 
                                    orthogroup, 
                                    pValue, 
                                    kValue,
                                    trait) 
  significanceInfo <- dplyr::filter(significanceInfo, 
                                    trait == specificTrait)
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

  # Combine all of the results for GO terms enriched in the foreground:
  goTermsRelaxed <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
    dplyr::filter(raw.p.value <= 0.01)
  goTermsRelaxed$selectionCategory <- "Relaxed selection associated with the trait"
  
  goTermsRelaxed
  
  significanceInfo <- dplyr::select(relaxResults, 
                                    orthogroup, 
                                    pValue, 
                                    kValue,
                                    trait) 
  significanceInfo <- dplyr::filter(significanceInfo, 
                                    trait == specificTrait)
  # Set each gene to 1 if k >1, 0, otherwise
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
  
  # Combine all of the results for GO terms enriched in the foreground:
  goTermsIntensified <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
    dplyr::filter(raw.p.value <= 0.01)
  goTermsIntensified$selectionCategory <- "Intensified selection associated with the trait"
  goTermsIntensified
  goTermsRELAX <- rbind(goTermsIntensified, 
                        goTermsRelaxed)
  goTermsRELAX$trait <- specificTrait
  goTermsRELAX
}

goTermsRELAX <- purrr::map(traits, goEnrichmentRELAX)

# Convert the outputs to dataframes:
goTermsBustedPH <- as.data.frame(do.call(rbind, goTermsBustedPH))   
goTermsRELAX <- as.data.frame(do.call(rbind, goTermsRELAX))   
goTermsAll <- rbind(goTermsBustedPH,
                    goTermsRELAX)

# Add names to the vector of traits, so that they can be nicely formatted for the tables:
names(traits) <- c("multilineage colonies", 
                   "polyandry",
                   "polygyny",
                   "worker reproduction",
                   "worker polymorphism")

# Make a nice table:
subTables <- function(specificTrait, formattedTrait) {
  goTermsTable <- goTermsAll %>% 
    filter(trait == specificTrait) %>%
    dplyr::select(-c(trait, 
                     Expected)) %>%
    dplyr::select(c(`p-value` = raw.p.value,
                    `GO ID` = GO.ID,
                    `Total orthogroups annotated` = Annotated,
                    `Significant orthogroups` = Significant,
                    `Type of selection` = selectionCategory,
                    Term = Term)) %>%
    gt(rowname_col = "Term",
       groupname_col = "Type of selection") %>%
    tab_header(title = paste("GO terms associated with", formattedTrait, sep = " ")) %>% 
    tab_options(table.width = pct(100),
                container.width = pct(100),
                container.height = pct(100),
                row_group.background.color = "#f2f2fa",
                heading.title.font.size = 30,
                heading.title.font.weight = "bolder") %>%
    cols_width(everything() ~ pct(25)) %>%
    gtsave(paste("GOTerms_", 
                 specificTrait, 
                 ".png", 
                 sep = ""), 
           path = "./Plots/",
           vwidth = 1000)
  
}

goTermsAllTable <- purrr::map2(traits, names(traits), subTables)

listOfGOTables <- list.files(path = "./Plots", 
                             pattern = "GOTerms_*", 
                             full.names = TRUE)

plotGTTable <- function(file) {
  plot <- ggdraw() + 
    draw_image(file, 
               scale = 0.8) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  return(plot)
}

testPlotList <- purrr::map(listOfGOTables, 
                   plotGTTable)

test <- plot_grid(plotlist = testPlotList, 
                  ncol = 1)
test
ggsave(filename = "goTermsAll.pdf", path = "./Plots/")



