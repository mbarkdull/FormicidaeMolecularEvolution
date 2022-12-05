# Script to create supplemental tables for the manuscript. 
library(flextable)
library(googlesheets4)
library(tidyverse)

# Read in data for the phenotype sources:
gs4_auth()
phenotypeData <- read_sheet("https://docs.google.com/spreadsheets/d/1O0aVh6LkTpS4MvJjDegTTy1IrIueJrY65s9xwb8DML8/edit?usp=sharing") %>%
  filter(Using == "yes")

# Select the columns you want to use:
phenotypeDataColumns <- dplyr::select(phenotypeData,
                                  Species,
                                  Subfamily,
                                  WorkerReproductionOffspringType,
                                  WorkerReproductionOffspringTypeCitation,
                                  DegreeOfWorkerPolymorphism,
                                  DegreeOfWorkerPolymorphismCitation,
                                  BioProject,
                                  GenomeCitation) %>%
  arrange(desc(Subfamily)) 

# Rename the columns:
phenotypeDataColumns <- phenotypeDataColumns %>%
  dplyr::rename("Species" = "Species",
         "Subfamily" = "Subfamily",
         "Worker reproductive capacity" = "WorkerReproductionOffspringType",
         "Citation for worker reproductive capacity" = "WorkerReproductionOffspringTypeCitation",
         "Worker subcaste type" = "DegreeOfWorkerPolymorphism",
         "Citation for worker subcaste" = "DegreeOfWorkerPolymorphismCitation",
         "BioProject accession number" = "BioProject",
         "Genome citation" = "GenomeCitation")

# Group the rows by subfamily:
phenotypeDataColumns <- as_grouped_data(x = phenotypeDataColumns,
                                        groups = c("Subfamily"))

# Make and save a .docx flextable:
phenotypeDataColumns <- flextable(phenotypeDataColumns)
phenotypeDataColumns <- theme_vanilla(phenotypeDataColumns)
phenotypeDataColumns
save_as_docx(phenotypeDataColumns,
             path = "./Plots/supplementalGenomeDataTable.docx")

# Do the same for the candidate genes:
candidateGenesSupplement <- read_sheet("https://docs.google.com/spreadsheets/d/1nOIuzA7f6yJbdXLH1A4LgXRo-HSKM8xRsXDCMr8PYZQ/edit?usp=sharing")

# Make a nice table for the publication:
publicationTable <- candidateGenesSupplement %>%
  dplyr::select(...1, 
         `candidate genes`,
         `gene symbol`,
         selectionOn,
         SelectionIntensity,
         pValue,
         pValueFDR,
         `test results p-value`,
         `test results background p-value`,
         `test results shared distributions p-value`,
         testResultspValueFDR,
         testResultsBackgroundpValueFDR,
         testResultsSharedDistributionspValueFDR,
         `paper citation`) %>%
  filter(...1 != "multilineage") %>%
  arrange(desc(...1))

publicationTable <- publicationTable %>% 
  mutate(trait = case_when(...1 == "workerReproductionQueens" ~ "Worker reproduction",
                           ...1 == "workerPolymorphism" ~ "Worker polymorphism"))

publicationTable <- publicationTable %>% dplyr::rename("Focal trait" = "trait",
                                                "Candidate gene" = "candidate genes",
                                                "NCBI gene symbol" = "gene symbol",
                                                "Positive selection" = "selectionOn",
                                                "Selection intensity" = "SelectionIntensity",
                                                "unadjusted RELAX p-value" = "pValue",
                                                "FDR-adjusted RELAX p-value" = "pValueFDR",
                                                "unadjusted BUSTED-PH p-value, foreground selection" = "test results p-value",
                                                "unadjusted BUSTED-PH p-value, background selection" = "test results background p-value",
                                                "unadjusted BUSTED-PH p-value, difference in selective regimes" = "test results shared distributions p-value",
                                                "FDR-adjusted BUSTED-PH p-value, foreground selection" = "testResultspValueFDR",
                                                "FDR-adjusted BUSTED-PH p-value, background selection" = "testResultsBackgroundpValueFDR",
                                                "FDR-adjusted BUSTED-PH p-value, difference in selective regimes" = "testResultsSharedDistributionspValueFDR",
                                                "Originally discussed in" = "paper citation") %>%
  dplyr::select(-c("...1"))

library(flextable)

table <- as_grouped_data(x = publicationTable,
                                    groups = c("Focal trait"))

table <- flextable(table)
table <- theme_vanilla(table)
table <- set_caption(table,
                     caption = "Selective regimes for candidate genes")
table
save_as_docx(table,
             path = "./Plots/suppplementalCandidateGenesPublicationTable.docx")


