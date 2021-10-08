library(ape)
library(geiger)
library(nlme)
library(phytools)
library(googlesheets4)
library(motmot)

absrelJSONProcessing <- function(file, analysisType) {
  # Read in the JSON file:
  aBSRELJSON <- RJSONIO::fromJSON(content = file)
  
  # Create a data frame with the information for each branch tested:
  aBSRELDataframe <- data.table::rbindlist(aBSRELJSON[["branch attributes"]][["0"]], fill = TRUE, idcol = TRUE)
  aBSRELDataframe <- subset(aBSRELDataframe, select = -c(`Rate Distributions`))
  aBSRELDataframe <- distinct (aBSRELDataframe)
  
  #Add a column to hold the orthogroup name:
  orthogroup <- sapply(strsplit(aBSRELJSON[["input"]][["file name"]],"/"), tail, 1)
  orthogroup <- sapply(strsplit(orthogroup, "_"), `[`, 2)
  aBSRELDataframe$Orthogroup <- orthogroup
  
  # Add a column to label genes as having been part of the foreground or background
  aBSRELDataframe$ForegroundBackground <- analysisType
  
  return(aBSRELDataframe)
}

# List all of the foreground analysis results files:
foregroundFiles <- list.files(path = "/Users/meganbarkdull/Projects/FormicidaeMolecularEvolution/9_2_ABSRELResults/multilineage/foreground", pattern = "*.json", full.names = TRUE)
# Drop any files with file size of zero:
foregroundFiles <- foregroundFiles[sapply(foregroundFiles, file.size) > 0]

# Make a "Safe" version of the function that will return an error if there is a problem with a particular file:
possiblyabsrelJSONProcessing <- possibly(absrelJSONProcessing, otherwise = "File empty.")
# Map the function over all foreground files:
foregroundaBSRELResults <- map_dfr(foregroundFiles, analysisType = "foreground", possiblyabsrelJSONProcessing)

# Do the same for background files:
backgroundFiles <- list.files(path = "/Users/meganbarkdull/Projects/FormicidaeMolecularEvolution/9_2_ABSRELResults/multilineage/background", pattern = "*.json", full.names = TRUE)
# Drop any files with file size of zero:
backgroundFiles <- backgroundFiles[sapply(backgroundFiles, file.size) > 0]

# Map the function over all background results:
backgroundaBSRELResults <- map_dfr(backgroundFiles, analysisType = "background", possiblyabsrelJSONProcessing)

# Combine the foreground and background tables:
allResults <- rbind(backgroundaBSRELResults, foregroundaBSRELResults)

# Remove rows with NAs in the corrected p-value columns
allResults <- drop_na(allResults, `Corrected P-value`)

# Remove results with NAs in the original name column; these are internal branches. 
allResults <- drop_na(allResults, `original name`)

# Create a column for species:
allResults$Species <- sapply(strsplit(allResults$.id, "_"), head, 1)

# Group by species, get the count of genes for each species, the number of genes with evidence for positive selection, and turn into a dataframe with one row per species:
groupedResults <- allResults %>%
  group_by(Species) %>%
  mutate(allGenes = n(),
         significantGenes = sum(`Corrected P-value` <= 0.05)) %>% 
  ungroup() %>%
  select(ForegroundBackground, Species, allGenes, significantGenes) %>% 
  distinct()

# What proportion of genes are under selection in each species?
groupedResults$realProportion <- groupedResults$significantGenes / groupedResults$allGenes


# Read in tree:
tree <- read.tree("./Tree/StudyTree.tre")
plot(tree)

# Read in phenotype data:
phenotypesRaw <- read_sheet("https://docs.google.com/spreadsheets/d/1O0aVh6LkTpS4MvJjDegTTy1IrIueJrY65s9xwb8DML8/edit?usp=sharing") 

# Select only columns of interest:
phenotypes <- phenotypesRaw %>%
  filter(Using == "yes") %>%
  select(SpeciesUnderscore, SpeciesAbbreviation, Polygyny) %>%
  as.data.frame()

# Combine the phenotype data and the positive selection data:
completeData <- merge(phenotypes, groupedResults, by.x = "SpeciesAbbreviation", by.y = "Species")
completeData <- completeData %>%
  mutate(polygynyQuantitative = case_when(
    Polygyny == "Yes" ~ 1,
    Polygyny == "Sometimes" ~ 1,
    Polygyny == "Two different morphs" ~ 1,
    Polygyny == "No" ~ 0,))

# Assign row names:
row.names(completeData) <- completeData$SpeciesUnderscore

# Make sure the species names match between the tree and the data:
test <- name.check(tree, completeData)

# Get vectors for our two traits 
polygynyQuantitative <- completeData$polygynyQuantitative
positiveSelection <- completeData$realProportion

# Give them names
names(polygynyQuantitative) <- names(positiveSelection) <- rownames(completeData)

# Create a vector of species names:
species <- rownames(completeData)

# Calculate PICs
polygynyPic <- pic(polygynyQuantitative, tree)
positiveSelectionPic <- pic(positiveSelection, tree)

# Make a model
picModel <- lm(polygynyPic ~ positiveSelectionPic - 1)

# No, not significant
summary(picModel)

# plot results
plot(polygynyPic ~ positiveSelectionPic)
abline(a = 0, b = coef(picModel))

# With pgls
pglsModel <- gls(positiveSelection ~ polygynyQuantitative, 
                 correlation = corBrownian(phy = tree, 
                                           form = ~species),
                 data = completeData, 
                 method = "ML")
summary(pglsModel)
plot(positiveSelection ~ polygynyQuantitative)
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])


# Compute Pagel's lambda to see if there's even phylogenetic signal:
lambda <- phylosig(tree, polygynyQuantitative, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
         control=list())

test <- as.matrix(polygynyQuantitative)

lambda.ml <- transformPhylo.ML(phy = tree, y = test, 
                               model = "lambda")
lambda.ml
