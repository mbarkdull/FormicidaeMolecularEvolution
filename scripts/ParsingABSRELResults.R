# Parsing aBSREL results

# I am testing for positive selection on foreground and background branches separately, because this improves power. I want to know if there are more genes under positive selection in the foreground than you would expect by chance.

# I could use:
  # Fisher's exact test
  # for each species, calculate the proportion of genes under positive selection, then see if this proportion is higher for foreground species (https://academic.oup.com/gbe/article/13/7/evab103/6275684?login=true#267264559). 
    # "We coded species according to their diet (carnivore, omnivore, or herbivore), their microhabitat (terrestrial, arboreal, or semi-aquatic), and their reproductive output, (based on the number of mammae for each species, supplementary table S3, Supplementary Material online). Using the time-calibrated species tree inferred in MCMCtree, we performed phylogenetic generalized least squares (PGLS) regression and phylogenetic ANOVA with a Bonferroni correction in phytools (Revell 2012), to test the effects of diet, microhabitat, reproductive output, and body size on genome-wide positive selection."


library(RJSONIO)
library(boot)

# Write a function to read in the aBSREL JSON result file and transform it into a dataframe:
absrelJSONProcessing <- function(file, analysisType) {
  OG0008911_absrel <- RJSONIO::fromJSON(content = file)
  
  # Create a data frame with the information for each branch tested:
  OG0008911_absrelResult <- data.table::rbindlist(OG0008911_absrel[["branch attributes"]][["0"]], fill = TRUE, idcol = TRUE)
  OG0008911_absrelResult <- subset(OG0008911_absrelResult, select = -c(`Rate Distributions`))
  OG0008911_absrelResult <- distinct (OG0008911_absrelResult)
  
  #Add a column to hold the orthogroup name:
  orthogroup <- sapply(strsplit(OG0008911_absrel[["input"]][["file name"]],"/"), tail, 1)
  orthogroup <- sapply(strsplit(orthogroup, "_"), `[`, 2)
  OG0008911_absrelResult$Orthogroup <- orthogroup
  
  # Add a column to label genes as having been part of the foreground or background
  OG0008911_absrelResult$ForegroundBackground <- analysisType

  return(OG0008911_absrelResult)
}

# List all of the foreground analysis results files:
foregroundFiles <- list.files(path = "/Users/meganbarkdull/Projects/FormicidaeMolecularEvolution/9_2_ABSRELResults/polygyny/foreground", pattern = "*.json", full.names = TRUE)
# Drop any files with file size of zero:
foregroundFiles <- foregroundFiles[sapply(foregroundFiles, file.size) > 0]

# Make a "Safe" version of the function that will return an error if there is a problem with a particular file:
possiblyabsrelJSONProcessing <- possibly(absrelJSONProcessing, otherwise = "File empty.")
# Map the function over all foreground files:
foregroundaBSRELResults <- map_dfr(foregroundFiles, analysisType = "foreground", possiblyabsrelJSONProcessing)

# Do the same for background files:
backgroundFiles <- list.files(path = "/Users/meganbarkdull/Projects/FormicidaeMolecularEvolution/9_2_ABSRELResults/polygyny/background", pattern = "*.json", full.names = TRUE)
# Drop any files with file size of zero:
backgroundFiles <- backgroundFiles[sapply(backgroundFiles, file.size) > 0]
backgroundaBSRELResults <- map_dfr(backgroundFiles, analysisType = "background", possiblyabsrelJSONProcessing)

# Combine the foreground and background tables:
allResults <- rbind(backgroundaBSRELResults, foregroundaBSRELResults)

# Remove rows with NAs in the corrected p-value columns
allResults <- drop_na(allResults, `Corrected P-value`)

hist(allResults$`Corrected P-value`)

# How many significant p-values are on foreground branches?
ForegroundSignificant <- nrow(allResults[ForegroundBackground == "foreground" & `Corrected P-value` < 0.05, ])

# If I reshuffle the labels foreground/background, how many?
countForeground <- function(allResults) {
  shuffled <- transform(allResults, ForegroundBackground = sample(ForegroundBackground) )
  shuffledForegroundSignificant <- nrow(shuffled[ForegroundBackground == "foreground" & `Corrected P-value` < 0.05, ])
}

# If I shuffle 100 times, what is my distribution of counts?
replicates <- replicate(500, {
  countForeground(allResults = allResults)
})
hist(replicates)
mean <- mean(replicates)
sd(replicates)

