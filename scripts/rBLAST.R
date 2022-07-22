#devtools::install_github("mhahsler/rBLAST")
library(rBLAST)
library(tidyverse)
library(googlesheets4)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/usr/local/ncbi/blast/bin/", sep= .Platform$path.sep))

# Read in concatenated genomes:
database <- Biostrings::readDNAStringSet('./allGenomes.fasta', format = 'fasta')

# Read in the orthogroups file:
orthogroupMembers <- readr::read_table2("./5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Orthogroups/Orthogroups.txt", col_names = FALSE) %>%
  tidyr::pivot_longer(cols = -c(X1),
                      names_to = "junk") %>%
  dplyr::select(-junk) %>%
  dplyr::distinct()
orthogroupMembers$X1 <- gsub(':', '', orthogroupMembers$X1)

# Construct the database:
makeblastdb(file = './allGenomes.fasta', dbtype = "nucl")

# List files:
files <- list.files(path = "./genesOfInterest", pattern = "*.fasta", full.names = TRUE)

# Write a function that will blast the candidate gene sequences against my ant orthogroup sequences:
fetchingOrthogroups <- function(i) {
  # Read in my query sequence:
  querySequence <- readDNAStringSet(i, format = 'fasta')
  
  # Prep the search:
  blastSearch <- blast(db = './allGenomes.fasta', type = 'blastn')
  # Run the search:
  searchResults <- predict(blastSearch, querySequence, BLAST_args = "-max_target_seqs 1")
  # Filter the search to include only matches with 100% identity:
  searchResults <- dplyr::filter(searchResults, Perc.Ident == max(Perc.Ident))
  
  # Extract the percent identity of the match:
  percentID <- max(searchResults$Perc.Ident)
  
  # Get out the matching sequence name:
  matchingSequence <- searchResults$SubjectID
  
  # Get the orthogroup that contains the sequence from our BLAST:
  orthogroup <- filter(orthogroupMembers, matchingSequence == value)
  orthogroup <- orthogroup$X1
  
  #Construct an object that contains that orthogroup, and the file it matches to, so you know which gene of interest goes with each orthogroup:
  results <- c(orthogroup, i, percentID)
  return(results)
}

# Make a safer version of that function with `possibly`:
possiblyfetchingOrthogroups <- purrr::possibly(fetchingOrthogroups, otherwise = "error")
# Map the function over all results files to construct a master dataframe:
orthogroupsToTest <- purrr::map(files, possiblyfetchingOrthogroups)
orthogroupsToTest <- as.data.frame(do.call(rbind, orthogroupsToTest))   

# Make a file to export to feed to a BASH script for running HYPHY tests:
write.table(orthogroupsToTest$V1, 
            file = "orthogroupsOfInterest.txt", 
            sep = "\t",
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)
write.table(orthogroupsToTest, 
            file = "orthogroupsOfInterestFullInfo.txt", 
            sep = "\t",
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

# Check the genes of interest for their selective regimes:
allHyphyResults <- read_csv("relaxAndBustedPH.csv", col_names = TRUE)

selectionOnGenesOfInterest <- filter(allHyphyResults,
                                     orthogroup %in% orthogroupsToTest$V1)
selectionOnGenesOfInterest <- full_join(orthogroupsToTest,
                                        selectionOnGenesOfInterest,
                                        by = c("V1" = "orthogroup"))

test <- dplyr::select(selectionOnGenesOfInterest,
                      c("orthogroup" = "V1",
                        "V2" = "V2",
                        "positive selection" = "selectionOn",
                        "selection intensity" = "shortDescription",
                        "trait" = "trait")) %>%
  #pivot_wider(names_from = c(trait),                     #No longer pivoting wider because I want just results for the relevant trait
              #values_from = c(`positive selection`, `selection intensity`)) %>%  
  dplyr::select(-ends_with(c("polyandry",
                             "polygyny",
                             "NA"))) %>%
  mutate(across(everything(), as.character))

test <- test %>%
  mutate(`positive selection` = case_when(`positive selection` == "EvidenceOfSelectionAssociatedWithTraitButNS" ~ "Nonsignificant evidence of positive selection associated with trait",
                                          `positive selection` == "SelectionOnBothButDifferent" ~ "Positive selection associated with trait presence and absence",
                                          `positive selection` == "EvidenceOfSelectionAssociatedWithLackOfTraitButNS" ~ "Nonsignificant evidence of positive selection associated with trait absence",
                                          `positive selection` == "BackgroundOnly" ~ "Positive selection associated with trait absence",
                                          `positive selection` == "SelectionOnBothButNoSignificantDifference" ~ "Positive selection associated with trait presence and absence",
                                          `positive selection` == "ForegroundOnly" ~ "Positive selection associated with trait presence",
                                          `positive selection` == "NoEvidenceOfSelection" ~ "No evidence of selection on orthogroup",
                                          `positive selection` == NA ~ "Orthogroup not tested for selection associated with this trait"))

test <- test %>%
  mutate(`selection intensity` = case_when(`selection intensity` == "Nonsignificant relaxation" ~ "Nonsignificant relaxation associated with trait",
                                          `selection intensity` == "Intensification of selection along foreground branches" ~ "Intensified selection significantly associated with trait",
                                          `selection intensity` == "Nonsignificant intensification" ~ "Nonsignificant intensification associated with trait",
                                          `selection intensity` == "Relaxation of selection along foreground branches" ~ "Relaxed selection significantly associated with trait"))

literatureReview <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1lVe5BsNmvGiF4jfC6ypwHSMDuEziXUQhsPqk7Xvf458/edit?usp=sharing")

# Join by trait and file name from literature review, and by trait and file name in the selection results. 
allResults <- full_join(test,
                        literatureReview,
                        by = c("V2" = "blast gene file name",
                               "trait" = "trait")) %>%
  filter(!is.na(`database query`),
         `candidate genes` != "-") %>%
  arrange(trait)
googlesheets4::write_sheet(allResults,
                           "https://docs.google.com/spreadsheets/d/1nOIuzA7f6yJbdXLH1A4LgXRo-HSKM8xRsXDCMr8PYZQ/edit?usp=sharing", 
                           sheet = "allResults")

# Make a nice table for the publication:
publicationTable <- allResults %>%
  select(trait, 
         `candidate genes`,
         `gene symbol`,
         orthogroup,
         `positive selection`,
         `selection intensity`) %>%
  filter(!is.na(orthogroup)) %>%
  arrange(desc(trait))

publicationTable <- publicationTable %>% 
  mutate(trait = case_when(trait == "workerReproductionQueens" ~ "Worker reproduction",
                           trait == "workerPolymorphism" ~ "Worker polymorphism"))

publicationTable <- publicationTable %>% 
  mutate(`Differential evolution` = case_when(`positive selection` == "Positive selection associated with trait absence"  ~ "Differentially evolving",
                                              `positive selection` == "Positive selection associated with trait presence" ~ "Differentially evolving",
                                              `selection intensity` == "Intensified selection significantly associated with trait" ~ "Differentially evolving",
                                              `selection intensity` == "Relaxed selection significantly associated with trait" ~ "Differentially evolving",
                                              TRUE ~ "None"))
publicationTable <- publicationTable %>% rename("Focal trait" = "trait",
                                                "Candidate gene" = "candidate genes",
                                                "NCBI gene symbol" = "gene symbol",
                                                "Ant orthogroup" = "orthogroup",
                                                "Positive selection" = "positive selection",
                                                "Selection intensity" = "selection intensity",
                                                "Differential evolution" = "Differential evolution")

library(flextable)

publicationTable <- as_grouped_data(x = publicationTable,
                                    groups = c("Focal trait"))

publicationTable <- flextable(publicationTable)
publicationTable <- theme_vanilla(publicationTable)
publicationTable <- set_caption(publicationTable,
                                caption = "Selective regimes for candidate genes")
publicationTable
save_as_docx(publicationTable,
             path = "./Plots/candidateGenesPublicationTable.docx")

