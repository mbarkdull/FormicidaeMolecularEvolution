#devtools::install_github("mhahsler/rBLAST")
library(rBLAST)
library(tidyverse)
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
literatureReview <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1lVe5BsNmvGiF4jfC6ypwHSMDuEziXUQhsPqk7Xvf458/edit?usp=sharing",
                                              sheet = "workerPolymorphism",
                                              col_names = TRUE) %>%
  filter(`including?` == "yes")
files <- literatureReview$`blast gene file name`

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

# Add back in filenames for the genes with no BLAST match:
  # Get all the filenames we tested as a dataframe:
  filenames <- as.data.frame(files)
  allOrthogroupsToTest <- full_join(orthogroupsToTest, filenames, by = c("V2" = "files"))
  # Drop any row where V2 isn't a valid filename (these are the rows corresponding to BLAST errors):
  allOrthogroupsToTest <- allOrthogroupsToTest %>%
    filter(str_detect(V2, "./"))
  
# Read in the literature search data from Google sheets:
literatureReview <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1lVe5BsNmvGiF4jfC6ypwHSMDuEziXUQhsPqk7Xvf458/edit?usp=sharing")

# Combine the BLASt results with the literature review info so I know which under which trait I'm testing each gene:
allOrthogroupsToTestPlusTraits <- right_join(literatureReview,
                                            allOrthogroupsToTest,
                                            by = c("blast gene file name" = "V2"))

# Check the genes of interest for their selective regimes:
allHyphyResults <- read_csv("relaxAndBustedPH.csv", col_names = TRUE) %>%
  dplyr::filter(trait != "polyandry",
                trait != "polygyny",
                trait != "multilineage")

genesOfInterest <- filter(allHyphyResults,
                          orthogroup %in% allOrthogroupsToTestPlusTraits$V1)
selectionOnGenesOfInterest <- full_join(allOrthogroupsToTestPlusTraits,
                                        genesOfInterest,
                                        by = c("V1" = "orthogroup", 
                                               "...1" = "trait")) %>%
  dplyr::select(-c(X1)) %>%
  drop_na(`paper number`) %>%
  drop_na(V1) %>%
  #drop_na(pValue) %>%
  #drop_na(`test results p-value`) %>%
  distinct() 

# Correct the p-values by the number of genes being tested for each trait:
correctingPValues <- selectionOnGenesOfInterest %>% group_by(`...1`) %>% 
  mutate(pValueFDR = p.adjust(pValue, method='BH')) 

correctingPValues <- correctingPValues %>%
  group_by(`...1`) %>% 
  mutate(testResultspValueFDR = p.adjust(`test results p-value`, method='BH')) %>% 
  mutate(testResultsBackgroundpValueFDR = p.adjust(`test results background p-value`, method='BH')) %>% 
  mutate(testResultsSharedDistributionspValueFDR = p.adjust(`test results shared distributions p-value`, method='BH'))

allResults <- correctingPValues %>%
  dplyr::select(c(`...1`,
           `candidate genes`,
           `gene symbol`,
           `gene function in original study`,
           pValue,
           pValueFDR,
           kValue,
           `test results p-value`,
           `test results background p-value`,
           `test results shared distributions p-value`,
           testResultspValueFDR,
           testResultsBackgroundpValueFDR,
           testResultsSharedDistributionspValueFDR,
           `paper citation`))

test <- allResults %>% mutate(selectionOn =
                                case_when(as.numeric(as.character(`test results p-value`)) <= 0.05 & 
                                            as.numeric(as.character(`test results background p-value`)) > 0.05 &
                                            as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 &
                                            as.numeric(as.character(pValueFDR)) >= 0.05 ~ "Positive selection in the foreground only, prior to FDR correction",
                                          as.numeric(as.character(`test results p-value`)) > 0.05 & 
                                            as.numeric(as.character(`test results background p-value`)) <= 0.05 &
                                            as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 &
                                            as.numeric(as.character(testResultsBackgroundpValueFDR)) >= 0.05 ~ "Positive selection in the background only prior to FDR correction",
                                          as.numeric(as.character(`testResultspValueFDR`)) <= 0.05 & 
                                            as.numeric(as.character(`testResultsBackgroundpValueFDR`)) > 0.05 &
                                            as.numeric(as.character(`testResultsSharedDistributionspValueFDR`)) <= 0.05 ~ "Positive selection in the foreground only, after FDR correction",
                                          as.numeric(as.character(`testResultspValueFDR`)) > 0.05 & 
                                            as.numeric(as.character(`testResultsBackgroundpValueFDR`)) <= 0.05 &
                                            as.numeric(as.character(`testResultsSharedDistributionspValueFDR`)) <= 0.05 ~ "Positive selection in the background only, after FDR correction",
                                           TRUE ~ "Trait presence or absence is not associated with a shift in selective regime."))

allResults <- test %>% mutate(SelectionIntensity =
                                case_when(as.numeric(as.character(`pValue`)) <= 0.05 &
                                            as.numeric(as.character(kValue)) <= 1 &
                                            as.numeric(as.character(`pValueFDR`)) >= 0.05 ~ "Relaxed selection is associated with the trait, prior to FDR correction",
                                          as.numeric(as.character(`pValue`)) <= 0.05 &
                                            as.numeric(as.character(kValue)) >= 1 &
                                            as.numeric(as.character(`pValueFDR`)) >= 0.05 ~ "Intensified selection is associated with the trait, prior to FDR correction",
                                          as.numeric(as.character(`pValueFDR`)) <= 0.05 &
                                            as.numeric(as.character(kValue)) <= 1 ~ "Relaxed selection is associated with the trait after FDR correction",
                                          as.numeric(as.character(`pValueFDR`)) <= 0.05 &
                                            as.numeric(as.character(kValue)) >= 1 ~ "Intensified selection is associated with the trait after FDR correction",
                                          TRUE ~ "Trait presence or absence is not associated with a shift in selective regime."))

allResults <- allResults %>% mutate(`Differential evolution` =
                                case_when(selectionOn == "Trait presence or absence is not associated with a shift in selective regime." &
                                          SelectionIntensity == "Trait presence or absence is not associated with a shift in selective regime." ~ "None",
                                          TRUE ~ "Some differential evolution."))

differentiallyEvolvingGenes <- allResults %>%
  filter(`Differential evolution` != "None") %>%
  distinct(...1, `gene symbol`, .keep_all = TRUE)

googlesheets4::write_sheet(allResults,
                           "https://docs.google.com/spreadsheets/d/1nOIuzA7f6yJbdXLH1A4LgXRo-HSKM8xRsXDCMr8PYZQ/edit?usp=sharing", 
                           sheet = "allResults")

# Make publication quality tables for these results:
library(flextable)
significantCandidateGenes <- allResults %>%
  drop_na(pValue) %>%
  drop_na(`test results p-value`) %>%
  dplyr::select(c(`...1`,
                  `candidate genes`,
                  `gene symbol`,
                  `gene function in original study`,
                  selectionOn,
                  SelectionIntensity,
                  `Differential evolution`,
                  pValue,
                  pValueFDR,
                  `test results p-value`,
                  `test results background p-value`,
                  `test results shared distributions p-value`,
                  testResultspValueFDR,
                  testResultsBackgroundpValueFDR,
                  testResultsSharedDistributionspValueFDR)) %>%
  dplyr::select(c(Trait = `...1`,
                  `Candidate gene` = `candidate genes`,
                  `NCBI gene symbol` = `gene symbol`,
                  `Gene function in original study` = `gene function in original study`,
                  `Positive selection` = selectionOn,
                  `Selection intensity` = SelectionIntensity,
                  `Differential evolution` = `Differential evolution`,
                  `unadjusted RELAX p-value` = pValue,
                  `FDR-adjusted RELAX p-value` = pValueFDR,
                  `unadjusted BUSTED-PH p-value, foreground selection` = `test results p-value`,
                  `unadjusted BUSTED-PH p-value, background selection` = `test results background p-value`,
                  `unadjusted BUSTED-PH p-value, difference in selective regimes` = `test results shared distributions p-value`,
                  `FDR-adjusted BUSTED-PH p-value, foreground selection` = testResultspValueFDR,
                  `FDR-adjusted BUSTED-PH p-value, background selection` = testResultsBackgroundpValueFDR,
                  `FDR-adjusted BUSTED-PH p-value, difference in selective regimes` = testResultsSharedDistributionspValueFDR)) %>%
  arrange(desc(Trait)) %>% 
  mutate(Trait = case_when(Trait == "workerReproductionQueens" ~ "Worker reproduction",
                           Trait == "workerPolymorphism" ~ "Worker polymorphism"))
significantCandidateGenes <- as_grouped_data(x = significantCandidateGenes, 
                                               groups = c("Trait")) 


significantCandidateGenes <- flextable(significantCandidateGenes)
significantCandidateGenes <- theme_vanilla(significantCandidateGenes)

save_as_docx(significantCandidateGenes,
             path = "./Plots/fullResultsSignificantCandidateGenes.docx")

significantCandidateGenes <- allResults %>%
  dplyr::select(c(`...1`,
                  `candidate genes`,
                  `gene symbol`,
                  `gene function in original study`,
                  selectionOn,
                  SelectionIntensity,
                  `Differential evolution`)) %>%
  dplyr::select(c(Trait = `...1`,
                  `Candidate gene` = `candidate genes`,
                  `NCBI gene symbol` = `gene symbol`,
                  `Gene function in original study` = `gene function in original study`,
                  `Positive selection` = selectionOn,
                  `Selection intensity` = SelectionIntensity,
                  `Differential evolution` = `Differential evolution`)) %>%
  arrange(desc(Trait)) %>% 
  mutate(Trait = case_when(Trait == "workerReproductionQueens" ~ "Worker reproduction",
                           Trait == "workerPolymorphism" ~ "Worker polymorphism"))
significantCandidateGenes <- as_grouped_data(x = significantCandidateGenes, 
                                             groups = c("Trait")) 


significantCandidateGenes <- flextable(significantCandidateGenes)
significantCandidateGenes <- theme_vanilla(significantCandidateGenes)

save_as_docx(significantCandidateGenes,
             path = "./Plots/significantCandidateGenes.docx")


