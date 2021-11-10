#devtools::install_github("mhahsler/rBLAST")
library(rBLAST)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/programs/ncbi-blast-2.9.0+/bin/", sep= .Platform$path.sep))

# Read in concatenated genomes:
database <- Biostrings::readDNAStringSet('./blast/allGenomes.fasta', format = 'fasta')

# Read in the orthogroups file:
orthogroupMembers <- readr::read_table2("./5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Orthogroups/Orthogroups.txt", col_names = FALSE) %>%
  tidyr::pivot_longer(cols = -c(X1),
                      names_to = "junk") %>%
  dplyr::select(-junk) %>%
  dplyr::distinct()
orthogroupMembers$X1 <- gsub(':', '', orthogroupMembers$X1)

# Construct the database:
makeblastdb(file = './blast/allGenomes.fasta', dbtype = "nucl")

# List files:
files <- list.files(path = "./blast/genesOfInterest", pattern = "*.fasta", full.names = TRUE)

fetchingOrthogroups <- function(i) {
  # Read in my query sequence:
  querySequence <- readDNAStringSet(i, format = 'fasta')
  
  # Prep the search:
  blastSearch <- blast(db = './blast/allGenomes.fasta', type = 'blastn')
  # Run the search:
  searchResults <- predict(blastSearch, querySequence, BLAST_args = "-max_target_seqs 1")
  # Filter the search to include only matches with 100% identity:
  searchResults <- dplyr::filter(searchResults, Perc.Ident == 100)
  
  # Get out the matching sequence name:
  matchingSequence <- searchResults$SubjectID
  
  # Get the orthogroup that contains the sequence from our BLAST:
  orthogroup <- filter(orthogroupMembers, matchingSequence == value)
  orthogroup <- orthogroup$X1
  
  #Construct an object that contains that orthogroup, and the file it matches to, so you know which gene of interest goes with each orthogroup:
  results <- c(orthogroup, i)
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
