# Load required packages:
library(RJSONIO)
library(tidyverse)
library(ggthemes)
library(gt)

# Write a function to read in RELAX results for a single trait:
processingBUSTEDPHFiles <- function(trait) {
  print("Listing files")
  jsonFiles <- list.files(path = paste("./9_3_BustedPHResults/", trait, sep = ""), pattern = "*.json", full.names = TRUE)
  jsonFiles <- sort(jsonFiles, decreasing = TRUE)
  jsonFiles <- jsonFiles[sapply(jsonFiles, file.size) > 0]
  
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
  bustedPHResults <- purrr::map(jsonFiles, possiblyparsingBustedPH)
  
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
  
  return(bustedPHResults)
}

# Get a list of traits:
traits <- c("multilineage", "polyandry", "polygyny", "workerReproductionQueens", "workerPolymorphism")
# Write a safe version with possibly:
possiblyprocessingBUSTEDPHFiles <- possibly(processingBUSTEDPHFiles, otherwise = "File empty.")
# Run this function with purrr:map so as to avoid for loops (from https://www.r-bloggers.com/2017/12/skip-errors-in-r-loops-by-not-writing-loops/ and https://jennybc.github.io/purrr-tutorial/ls01_map-name-position-shortcuts.html):
bustedPHResults <- purrr::map(traits, possiblyprocessingBUSTEDPHFiles)
bustedPHResults <- as.data.frame(do.call(rbind, bustedPHResults))   
#colnames(bustedPHResults) <- c("fileName", "pValue", "kValue", "description", "shortDescription", "orthogroup", "trait")

# Calculate p-values for each trait:
pValueByTrait <- function(specificTrait) {
  bustedPHResultsTrait <- filter(bustedPHResults, 
                                 trait == specificTrait)
  nSelectionForeground <- bustedPHResultsTrait %>% 
    filter(`test results p-value` <= 0.05,
           `test results background p-value` > 0.05,
           `test results shared distributions p-value` <= 0.05) %>%
    nrow()
  
  nSelectionBackground <- bustedPHResultsTrait %>% 
    filter(`test results p-value` > 0.05,
           `test results background p-value` <= 0.05,
           `test results shared distributions p-value` <= 0.05) %>%
    nrow()
  # Construct a frequency table with that information:
  frequencyTableWithNoCorrection <- matrix(c(nSelectionForeground, nSelectionBackground, 
                                             (nrow(bustedPHResultsTrait) - nSelectionForeground),
                                             (nrow(bustedPHResultsTrait) - nSelectionBackground)),
                                           nrow = 2,
                                           byrow = TRUE,
                                           dimnames = list("category" = c("selected", "noSelection") ,
                                                           "selection" = c("foreground", "background")))
  # Run a fisher's exact test to see if there is a difference in the proportion of genes under selection in the fore- vs. background:
  fishersExactTest <- fisher.test(frequencyTableWithNoCorrection)
  # Print the resulting p-value:
  fishersExactTest[["p.value"]]
  pValue <- if (fishersExactTest[["p.value"]] < 0.000001) {
    formatC(fishersExactTest[["p.value"]], format = "e", digits = 2)
  } else {
    round(fishersExactTest[["p.value"]], digits = 4)
  }
  textHeight <- as.numeric(max(nSelectionForeground, nSelectionBackground))
  return(c(specificTrait, pValue, textHeight, nSelectionForeground, nSelectionBackground))
}
possiblypValueByTrait <- possibly(pValueByTrait, otherwise = "Error")
pValues <- purrr::map(traits, possiblypValueByTrait)
pValues <- as.data.frame(do.call(rbind, pValues))   
colnames(pValues) <- c("trait", "pValue", "maxHeight", "nSelectionForeground", "nSelectionBackground")

#### Plotting
bustedPHResults$selectionOn <- factor(bustedPHResults$selectionOn,
                                      levels = c("ForegroundOnly",
                                                 "BackgroundOnly",
                                                 "EvidenceOfSelectionAssociatedWithTraitButNS",
                                                 "EvidenceOfSelectionAssociatedWithLackOfTraitButNS",
                                                 "SelectionOnBothButDifferent",
                                                 "SelectionOnBothButNoSignificantDifference",
                                                 "NoEvidenceOfSelection"))

plot <- ggplot(data = filter(bustedPHResults,
                             file != "File empty.")) + 
  geom_bar(mapping = aes(x = selectionOn),
           position = "dodge") +
  labs(x = "Selective regime", 
       y = "Count of orthogroups", 
       title = "Association of positive selection and focal traits") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) + 
  scale_x_discrete(labels=c("ForegroundOnly" = "Selection on species\nwith the trait only",
                            "BackgroundOnly" = "Selection on species\nwithout the trait only",
                            "EvidenceOfSelectionAssociatedWithTraitButNS" = "Selection on species\nwith the trait,\nbut difference with\nspecies lacking the trait\nis non-signficant",
                            "EvidenceOfSelectionAssociatedWithLackOfTraitButNS"= "Selection on species\nwithout the trait,\nbut difference with species\nlacking the trait\nis non-signficant",
                            "SelectionOnBothButDifferent" = "Signficant evidence for\nselection in species with\nand without the trait\nwith signficant differences\nin selective regime",
                            "SelectionOnBothButNoSignificantDifference" = "Signficant evidence for\nselection in species with\nand without the trait\nbut with no\nsignficant differences\nin selective regime",
                            "NoEvidenceOfSelection" = "No evidence\nfor any selection")) +
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
  facet_wrap(~trait,
             nrow = 3) 

plot

# Expand out to the individual gene level:
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

selectionOnGenes <- right_join(bustedPHResults,
                               orthogroupMembershipLong,
                               by = c("orthogroup" = "Orthogroup")) %>% 
  drop_na(file)
selectionOnGenes$Species <- sapply(strsplit(as.character(selectionOnGenes$Gene),
                                            "_"),
                                   head,
                                   1)


plot <- ggplot(data = filter(selectionOnGenes,
                             file != "File empty."))
plot + 
  geom_bar(mapping = aes(x = selectionOn),
           position = "dodge") +
  labs(x = "Selective regime", 
       y = "Count of individual genes", 
       title = "Distribution of selective regimes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) + 
  scale_x_discrete(labels=c("ForegroundOnly" = "Selection on species\nwith the trait only",
                            "BackgroundOnly" = "Selection on species\nwithout the trait only",
                            "EvidenceOfSelectionAssociatedWithTraitButNS" = "Selection on species\nwith the trait,\nbut difference with\nspecies lacking the trait\nis non-signficant",
                            "EvidenceOfSelectionAssociatedWithLackOfTraitButNS"= "Selection on species\nwithout the trait,\nbut difference with species\nlacking the trait\nis non-signficant",
                            "SelectionOnBothButDifferent" = "Signficant evidence for\nselection in species with\nand without the trait\nwith signficant differences\nin selective regime",
                            "SelectionOnBothButNoSignificantDifference" = "Signficant evidence for\nselection in species with\nand without the trait\nbut with no\nsignficant differences\nin selective regime",
                            "NoEvidenceOfSelection" = "No evidence\nfor any selection")) + 
  facet_wrap(~trait) 
