# Load required packages:
library(tidyverse)
library(ggthemes)
library(gt)
library(plyr)
library(rjson)

traits <- c("polyandry", "polygyny", "workerReproductionQueens", "workerPolymorphism", "multilineage")

processingRELAXFiles <- function(trait) {
  jsonFiles <- list.files(path = paste("./10_1_RelaxResults/", trait, sep = ""), pattern = "*.json", full.names = TRUE)
  #jsonFiles <- list.files(path = "./10_1_RelaxResults/workerPolymorphism", pattern = "*.json", full.names = TRUE)
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
  relaxResults$trait <- trait
  return(relaxResults)
}

# Write a safe version with possibly:
possiblyprocessingRELAXFiles <- possibly(processingRELAXFiles, otherwise = "File empty.")
# Run this function with purrr:map so as to avoid for loops (from https://www.r-bloggers.com/2017/12/skip-errors-in-r-loops-by-not-writing-loops/ and https://jennybc.github.io/purrr-tutorial/ls01_map-name-position-shortcuts.html):
relaxResults <- purrr::map(traits, possiblyprocessingRELAXFiles)
relaxResults <- as.data.frame(do.call(rbind, relaxResults))  

# Create a column that accurate categories genes by selective regime:
relaxResults <- relaxResults %>% mutate(selectionCategory =
                                          case_when(pValue <= 0.05 & kValue > 1 ~ "signficantIntensification",
                                                    pValue <= 0.05 & kValue < 1 ~ "signficantRelaxation",
                                                    pValue > 0.05 & kValue > 1 ~ "nonsignficantIntensification",
                                                    pValue > 0.05 & kValue < 1 ~ "nonsignficantRelaxation")
                                        )

#Calculate p-values for each trait:
pValueByTrait <- function(specificTrait) {
  relaxResultsTrait <- filter(relaxResults, 
                                 trait == specificTrait)
  nSelectionIntensified <- relaxResultsTrait %>% 
    filter(pValue <= 0.05,
           kValue > 1) %>%
    nrow()
  
  nSelectionRelaxed <- relaxResultsTrait %>% 
    filter(pValue <= 0.05,
           kValue < 1) %>%
    nrow()
  # Construct a frequency table with that information:
  frequencyTableWithNoCorrection <- matrix(c(nSelectionIntensified, nSelectionRelaxed, 
                                             (nrow(relaxResultsTrait) - nSelectionIntensified),
                                             (nrow(relaxResultsTrait) - nSelectionRelaxed)),
                                           nrow = 2,
                                           byrow = TRUE,
                                           dimnames = list("category" = c("evidenceForSelection", "nonsignificantResult") ,
                                                           "selection" = c("intensified", "relaxed")))
  # Run a fisher's exact test to see if there is a difference in the proportion of genes under selection in the fore- vs. background:
  fishersExactTest <- fisher.test(frequencyTableWithNoCorrection)
  # Print the resulting p-value:
  fishersExactTest[["p.value"]]
  pValue <- if (fishersExactTest[["p.value"]] < 0.000001) {
    formatC(fishersExactTest[["p.value"]], format = "e", digits = 2)
  } else {
    round(fishersExactTest[["p.value"]], digits = 4)
  }
  textHeight <- as.numeric(max(nSelectionIntensified, nSelectionRelaxed))
  return(c(specificTrait, pValue, textHeight))
}
possiblypValueByTrait <- possibly(pValueByTrait, otherwise = "Error")
pValues <- purrr::map(traits, possiblypValueByTrait)
pValues <- as.data.frame(do.call(rbind, pValues))   
colnames(pValues) <- c("trait", "pValue", "maxHeight")

# Plot:
relaxResults$selectionCategory <- factor(relaxResults$selectionCategory,
                                      levels = c("signficantIntensification",
                                                 "signficantRelaxation",
                                                 "nonsignficantIntensification",
                                                 "nonsignficantRelaxation"))

plot <- ggplot(data = filter(relaxResults,
                             !is.na(selectionCategory)))
plot + 
  geom_bar(mapping = aes(x = selectionCategory),
           position = "dodge") +
  labs(x = "Selective regime", 
       y = "Count of orthogroups", 
       title = "Relationship between selective intensity and focal traits") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  geom_text(data = pValues, 
            aes(x = 1.5,
                y = as.numeric(maxHeight) + 300,
                label = pValue)) + 
  geom_linerange(data = pValues,
                 aes(xmin = 1, 
                     xmax = 2, 
                     y = (as.numeric(maxHeight) + 200)), 
                 color = "grey26", 
                 size = 0.3) + 
  scale_x_discrete(labels=c("nonsignficantRelaxation" = "Nonsignificant\nrelaxation",
                            "signficantIntensification" = "Intensification of\nselection in species\nwith focal trait",
                            "nonsignficantIntensification" = "Nonsignificant\nintensification",
                            "signficantRelaxation" = "Relaxation of\nselection in species\nwith focal trait")) + 
  facet_wrap(~trait,
             nrow = 3) 

# Plot distributions of k-values:
ggplot(data = filter(relaxResults,
                     fileName != "File empty." & pValue <= 0.05)) +
  geom_histogram(mapping = aes(x = log10(as.numeric(as.character(kValue)))),
           bins = 150) + 
  facet_wrap(~trait) 


ggplot(filter(relaxResults, 
              pValue <= 0.05), 
       aes(x = trait, 
           y = log(as.numeric(as.character(kValue))))) + 
  geom_violin()

# Read in the orthogroup membership data:
orthogroupMembership <- read_delim("5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Orthogroups/Orthogroups.tsv",
                                   delim = "\t")
# Pivot longer:
orthogroupMembershipLong <- pivot_longer(orthogroupMembership,
                                         cols = acol_filteredTranscripts:waur_filteredTranscripts,
                                         names_to = "Species",
                                         values_to = "Gene")
orthogroupMembershipLong <- orthogroupMembershipLong %>% 
  mutate(new = Gene) %>%
  group_by(Orthogroup) %>% 
  nest() %>% 
  mutate(
    temp_col = map(
      data, 
      ~ stringr::str_split(.x$new, pattern = ",") %>% 
        flatten() %>% 
        map_chr(~return(.x)) %>% 
        as_tibble()
    )
  ) %>% 
  unnest(temp_col) %>% 
  separate(value, into = c("Gene", "count"), sep = ", ") %>% 
  dplyr::select(Gene, count, Orthogroup) 

orthogroupMembershipLong$Gene <- trimws(orthogroupMembershipLong$Gene)
orthogroupMembershipLong$Gene <-gsub("\"",
                                     "",
                                     as.character(orthogroupMembershipLong$Gene))
orthogroupMembershipLong <- orthogroupMembershipLong %>% 
  drop_na(Gene) %>% 
  dplyr::select(Gene, Orthogroup) 

selectionOnGenes <- right_join(relaxResults,
                               orthogroupMembershipLong,
                               by = c("orthogroup" = "Orthogroup")) %>% 
  drop_na(fileName)
selectionOnGenes$Species <- sapply(strsplit(as.character(selectionOnGenes$Gene),
                                            "_"),
                                   head,
                                   1)

selectionOnGenes$shortDescription <- factor(selectionOnGenes$shortDescription,
                                        levels = c("Intensification of selection along foreground branches",
                                                   "Relaxation of selection along foreground branches",
                                                   "Nonsignificant intensification",
                                                   "Nonsignificant relaxation"))
plot2 <- ggplot(data = filter(selectionOnGenes,
                             fileName != "File empty."))
plot2 + 
  geom_bar(mapping = aes(x = shortDescription),
           position = "dodge") +
  labs(x = "Selective regime", 
       y = "Count of orthogroups", 
       title = "Distribution of selective regimes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) + 
  facet_wrap(~trait) 

# Looking at only single copy orthogroups:
singleCopyOrthogroups <- read_delim("5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Orthogroups/Orthogroups_SingleCopyOrthologues.txt",
                                    delim = "\t",
                                    col_names = FALSE)
                                   
singleCopyOrthogroupResults <- dplyr::filter(relaxResults, 
                                             relaxResults$orthogroup %in% singleCopyOrthogroups$X1)

pValueByTraitSingleCopy <- function(specificTrait) {
  relaxResultsTrait <- filter(singleCopyOrthogroupResults, 
                              trait == specificTrait)
  nSelectionIntensified <- relaxResultsTrait %>% 
    filter(pValue <= 0.05,
           kValue > 1) %>%
    nrow()
  
  nSelectionRelaxed <- relaxResultsTrait %>% 
    filter(pValue <= 0.05,
           kValue < 1) %>%
    nrow()
  # Construct a frequency table with that information:
  frequencyTableWithNoCorrection <- matrix(c(nSelectionIntensified, nSelectionRelaxed, 
                                             (nrow(relaxResultsTrait) - nSelectionIntensified),
                                             (nrow(relaxResultsTrait) - nSelectionRelaxed)),
                                           nrow = 2,
                                           byrow = TRUE,
                                           dimnames = list("category" = c("evidenceForSelection", "nonsignificantResult") ,
                                                           "selection" = c("intensified", "relaxed")))
  # Run a fisher's exact test to see if there is a difference in the proportion of genes under selection in the fore- vs. background:
  fishersExactTest <- fisher.test(frequencyTableWithNoCorrection)
  # Print the resulting p-value:
  fishersExactTest[["p.value"]]
  pValue <- if (fishersExactTest[["p.value"]] < 0.000001) {
    formatC(fishersExactTest[["p.value"]], format = "e", digits = 2)
  } else {
    round(fishersExactTest[["p.value"]], digits = 4)
  }
  textHeight <- as.numeric(max(nSelectionIntensified, nSelectionRelaxed))
  return(c(specificTrait, pValue, textHeight))
}
possiblypValueByTraitSingleCopy <- possibly(pValueByTraitSingleCopy, otherwise = "Error")
pValues <- purrr::map(traits, possiblypValueByTraitSingleCopy)
pValues <- as.data.frame(do.call(rbind, pValues))   
colnames(pValues) <- c("trait", "pValue", "maxHeight")


# Plot:
singleCopyOrthogroupResults$selectionCategory <- factor(singleCopyOrthogroupResults$selectionCategory,
                                                        levels = c("signficantIntensification",
                                                                   "signficantRelaxation",
                                                                   "nonsignficantIntensification",
                                                                   "nonsignficantRelaxation"))

plot <- ggplot(data = filter(singleCopyOrthogroupResults,
                             !is.na(selectionCategory)))
plot + 
  geom_bar(mapping = aes(x = selectionCategory),
           position = "dodge") +
  labs(x = "Selective regime", 
       y = "Count of orthogroups", 
       title = "Relationship between selective intensity and focal traits\nfor single-copy orthogroups only") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  geom_text(data = pValues, 
            aes(x = 1.5,
                y = as.numeric(maxHeight) + 300,
                label = pValue)) + 
  geom_linerange(data = pValues,
                 aes(xmin = 1, 
                     xmax = 2, 
                     y = (as.numeric(maxHeight) + 200)), 
                 color = "grey26", 
                 size = 0.3) + 
  scale_x_discrete(labels=c("nonsignficantRelaxation" = "Nonsignificant\nrelaxation",
                            "signficantIntensification" = "Intensification of\nselection in species\nwith focal trait",
                            "nonsignficantIntensification" = "Nonsignificant\nintensification",
                            "signficantRelaxation" = "Relaxation of\nselection in species\nwith focal trait")) + 
  facet_wrap(~trait,
             nrow = 3) 
