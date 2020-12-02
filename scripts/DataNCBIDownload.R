library(devtools)
install_github("ropensci/rentrez")
library(rentrez)

# Create a directory in which to store data retrieved from NCBI:
dir.create("./1_RawData")

# Set the API key so that you can get requests slightly faster:
set_entrez_key(args[1])

# Use rentrez to figure out the complete list of Formicidae species that have fairly complete gene sets available on NCBI?


# Write a function to pull all gene sequences for a species:
# Species looks like: "Monomorium pharaonis"
getNCBIsequences <- function(species, output) {
  # Get all gene IDs for the species:
  Sys.sleep(10)
  query <- paste(species, "[organism] AND (alive[prop])", sep = "")
  print("Query is:")
  print(query)
  GeneIDs <- entrez_search(db = "gene", term = query, use_history = TRUE, retmax = 30000)
  print("Gene IDs retrieved")
  # Now link those gene IDs to the nucleotide sequence:
  Sys.sleep(10)
  IDLinks <- entrez_link(dbfrom = "gene", web_history = GeneIDs$web_history, db = "nuccore", cmd = "neighbor_history")
  print("Gene IDs linked to sequence IDs")
  # Now pull down those nucleotide sequences:
  Sys.sleep(10)
  CDS <- entrez_fetch(db = "nuccore", web_history = IDLinks$web_histories$gene_nuccore, rettype = "fasta")
  print("Sequences fetched")
  # Now export the sequences:
  print("Exporting sequences")
  write(CDS, file = output)
}

# Function on one species:
getNCBIsequences(species = "Monomorium pharaonis", output = "mpha_transcripts.fasta")

# Read in the species data:
speciesInfo <- read.table(file = args[1], sep = ",")
speciesInfo <- read.table(file = "./scripts/inputurls_partial", sep = ",")
# Split the second column to get a column with only abbreviations:
#speciesInfo$V3 <- gsub("_", " ", speciesInfo$V3)
# Get a vector from that column:
species <- speciesInfo$V3

# For loop to iterate NCBI download over all species:
for (i in species)
{
  print(i)
  speciesName <- gsub("_", " ", i)
  print(speciesName)
  output <- paste(i, "_cds.fasta", sep = "")
  print(output)
  Sys.sleep(5)
  getNCBIsequences(species = speciesName, output = output)
  Sys.sleep(5)
  
}


# Code for a single species:
# Get all gene IDs for the species:
Sys.sleep(5)
mphaGeneIDs <- entrez_search(db = "gene", term = "Monomorium pharaonis[organism] AND (alive[prop])", use_history = TRUE, retmax = 30000)
Sys.sleep(5)
#mphaGeneIDs_ids <- entrez_search(db = "gene", term = "Monomorium pharaonis[organism] AND (alive[prop])", retmax = 30000)

# Now link those gene IDs to the nucleotide sequence:
mphaIDLinks <- entrez_link(dbfrom = "gene", web_history = mphaGeneIDs$web_history, db = "nuccore", cmd = "neighbor_history")
Sys.sleep(5)
# Now pull down those nucleotide sequences:
mphaCDS <- entrez_fetch(db = "nuccore", web_history = mphaIDLinks$web_histories$gene_nuccore, rettype = "fasta")
Sys.sleep(5)
# Now export the sequences:
write(mphaCDS, file = "./1_RawData/mpha_cds.fasta")
