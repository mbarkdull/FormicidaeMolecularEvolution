library(ape)
library(phylotools)
library(ggtree)
library(tidyverse)

# Read in the large tree file:
#FullTree <- read.tree(file = args[1])
fullTree <- read.nexus(file = "./Tree/BlanchardMoreau_ConcatenatedMatrix_MCCTree_Combined")

# Get list of tips to keep:
tipsToKeep <- c("Lasius_californicus",
                "Monomorium_pharaonis",
                "Acromyrmex_versicolor",
                "Atta_texana",
                "Camponotus_maritimus",
                "Cardiocondyla_mauritanica",
                "Eciton_vagans",
                "Harpegnathos_saltator",
                "Linepithema_humile",
                "Pogonomyrmex_angustus",
                "Solenopsis_molesta",
                "Cerapachys_dohertyi_Dohertyi_grp",
                "Formica_moki",
                "Vollenhovia_emeryi",
                "Wasmannia_auropunctata",
                "Cephalotes_texanus",
                "Pseudomyrmex_gracilis",
                "Nylanderia_hystrix",
                "Cyphomyrmex_cornutus",
                "Temnothorax_poeyi",
                "Dinoponera_australis",
                "Trachymyrmex_arizonensis",
                "Odontomachus_coquereli"
)
# Trim the tree to retain only those tips:
trimmedTree <- keep.tip(fullTree, tip = tipsToKeep)
# Plot the trimmed tree:
ggtree(trimmedTree) + 
  geom_tiplab(size = 5) +
  xlim(-100, 1000)

# Switch the tip labels to reflect the actual species in my study:
# Note it is CRUCIAL that they are in the same order as in the data structure for the first tree! In the first tree, they will be alphabetical, so if any genus names are changing, you'll need to take that into account (that is why Oocera here is out of place, because in the original tree it is Cerapachys).
newTipLabels <- c("Acromyrmex_echinatior",
                  "Atta_colombica",
                  "Camponotus_floridanus",
                  "Cardiocondyla_obscurior",
                  "Cephalotes_varians",
                  "Ooceraea_biroi",
                  "Cyphomyrmex_costatus",
                  "Dinoponera_quadriceps",
                  "Eciton_burchellii",
                  "Formica_exsecta",
                  "Harpegnathos_saltator",
                  "Lasius_niger",
                  "Linepithema_humile",
                  "Monomorium_pharaonis",
                  "Nylanderia_fulva",
                  "Odontomachus_brunneus",
                  "Pogonomyrmex_barbatus",
                  "Pseudomyrmex_gracilis",
                  "Solenopsis_invicta",
                  "Temnothorax_curvispinosus",
                  "Trachymyrmex_septentrionalis",
                  "Vollenhovia_emeryi",
                  "Wasmannia_auropunctata"
)
# Create a new tree object so you don't overwrite the original one:
trimmedTreeMySpecies <- trimmedTree
# Replace the original tip labels with my new vector of tip labels:
trimmedTreeMySpecies$tip.label <- newTipLabels
# Plot the new tree:
ggtree(trimmedTreeMySpecies) + 
  geom_tiplab(size = 5) +
  xlim(-100, 500)

# Export the species pruned tree:
write.tree(trimmedTreeMySpecies, file = "./Tree/StudyTree.tre")

#Construct a vector of file names corresponding to the species (for use with Orthofinder, and potentially with BUSTED)
#speciesInfo <- read.table(file = args[2], sep = ",")
speciesInfo <- read.table(file = "./scripts/inputurls_full.txt", sep = ",")

# Construct the filenames:
speciesInfo$transcripts <- paste(speciesInfo$V4, "_filteredTranscripts.fasta", sep = "")
# Sort the dataframe by the species so that we can use the filenames as new tip labels:
speciesInfo <- plyr::arrange(speciesInfo, V6)
filenames <- speciesInfo$transcripts
filenameTree <- trimmedTreeMySpecies

filenameTree$tip.label <- tolower(filenameTree$tip.label)
filenameTree$tip.label <- paste(str_sub(filenameTree$tip.label, 1, 1), substring(filenameTree$tip.label, regexpr("_", filenameTree$tip.label) + 1, regexpr("_", filenameTree$tip.label) + 3), "_filteredTranscripts.fasta", sep = "")

ggtree(filenameTree) + 
  geom_tiplab(size = 5) +
  xlim(-100, 1000)
write.tree(FilenameTree, file = "./Tree/FilenameTree.tree")
