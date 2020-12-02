library(devtools)
install_github("ropensci/rentrez")
library(rentrez)

dir.create("./1_RawData")


# Write a function to pull gene sequences for a species:
# Query looks like: "Monomorium pharaonis[organism] AND (alive[prop])"
getNCBIsequences <- function(species, output) {
  # Get all gene IDs for the species:
  Sys.sleep(0.5)
  query <- paste(species, "[organism] AND (alive[prop])", sep = "")
  print("Query is:")
  print(query)
  print("Getting gene IDs")
  GeneIDs <- entrez_search(db = "gene", term = query, use_history = TRUE, retmax = 30000)
  # Now link those gene IDs to the nucleotide sequence:
  Sys.sleep(0.5)
  print("Linking gene IDs to sequence IDs")
  IDLinks <- entrez_link(dbfrom = "gene", web_history = GeneIDs$web_history, db = "nuccore", cmd = "neighbor_history")
  # Now pull down those nucleotide sequences:
  Sys.sleep(0.5)
  print("Fetching sequences")
  CDS <- entrez_fetch(db = "nuccore", web_history = IDLinks$web_histories$gene_nuccore, rettype = "fasta")
  # Now export the sequences:
  print("Exporting sequences")
  write(CDS, file = output)
}

getNCBIsequences(species = "Monomorium pharaonis", output = "mpha_transcripts.fasta")



# Read in the input data:
speciesInfo <- read.table(file = args[1], sep = ",")
#speciesInfo <- read.table(file = "./scripts/inputurls_partial", sep = ",")
# Split the second column to get a column with only abbreviations:
speciesInfo$V3 <- gsub("_", " ", speciesInfo$V3)
# Get a vector from that column:
species <- speciesInfo$V3


# Code for a single species:
# Get all gene IDs for the species:
mphaGeneIDs <- entrez_search(db = "gene", term = "Monomorium pharaonis[organism] AND (alive[prop])", use_history = TRUE, retmax = 30000)
mphaGeneIDs_ids <- entrez_search(db = "gene", term = "Monomorium pharaonis[organism] AND (alive[prop])", retmax = 30000)

# Now link those gene IDs to the nucleotide sequence:
mphaIDLinks <- entrez_link(dbfrom = "gene", web_history = mphaGeneIDs$web_history, db = "nuccore", cmd = "neighbor_history")

# Now pull down those nucleotide sequences:
mphaCDS <- entrez_fetch(db = "nuccore", web_history = mphaIDLinks$web_histories$gene_nuccore, rettype = "fasta")

# Now export the sequences:
write(mphaCDS, file = "./1_RawData/mpha_cds.fasta")
