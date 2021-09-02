# Parsing aBSREL results

# I am testing for positive selection on foreground and background branches separately, because this improves power. I want to know if there are more genes under positive selection in the foreground than you would expect by chance.

# I could use:
  # Fisher's exact test
  # for each species, calculate the proportion of genes under positive selection, then see if this proportion is higher for foreground species (https://academic.oup.com/gbe/article/13/7/evab103/6275684?login=true#267264559). 
    # "We coded species according to their diet (carnivore, omnivore, or herbivore), their microhabitat (terrestrial, arboreal, or semi-aquatic), and their reproductive output, (based on the number of mammae for each species, supplementary table S3, Supplementary Material online). Using the time-calibrated species tree inferred in MCMCtree, we performed phylogenetic generalized least squares (PGLS) regression and phylogenetic ANOVA with a Bonferroni correction in phytools (Revell 2012), to test the effects of diet, microhabitat, reproductive output, and body size on genome-wide positive selection."

# Load necessary packages (I'm using RJSONIO here because rjson threw errors for one of my results files):
library(RJSONIO)
library(tidyverse)


# Write a function to read in the aBSREL JSON result file and transform it into a dataframe:
absrelJSONProcessing <- function(file, analysisType) {
  # Read in the JSON file:
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

# Make a quick histogram of the p-values, mostly for curiosity's sake:
hist(allResults$`Corrected P-value`)

# How many significant p-values are on foreground branches?
ForegroundSignificant <- nrow(filter(allResults, ForegroundBackground == "foreground" & `Corrected P-value` <= 0.05))

# Now I'll use a permutation test to see how extreme my number of significant foreground genes is:
# If I reshuffle the labels foreground/background, how many significant foreground genes?shuffleForeground <- function(allResults) {
  shuffled <- transform(allResults, ForegroundBackground = sample(ForegroundBackground, replace = FALSE))
  shuffledForegroundSignificant <- nrow(filter(shuffled, ForegroundBackground == "foreground" & `Corrected P-value` <= 0.05))
}

# If I shuffle 10000 times, what is my distribution of counts?
replicates <- replicate(1000, {
  shuffleForeground(allResults = allResults)
})
replicates <- sort(replicates)

# Plot the distribution of replicates:
hist(replicates)

# Get the mean and standard deviation of the distribution of replicates:
mean <- mean(replicates)
standardDeviation <- sd(replicates)
oneStandardDeviationHigher <- mean + standardDeviation
oneStandardDeviationLower <- mean - standardDeviation
twoStandardDeviationsHigher <- oneStandardDeviationHigher - standardDeviation
twoStandardDeviationsLower <- oneStandardDeviationLower - standardDeviation

# How many times did the replication process generate values more extreme than my real data:
ForegroundSignificant <- as.integer(ForegroundSignificant)
replicatesGreater <- sum(replicates > ForegroundSignificant)
replicatesLess <- sum(replicates < ForegroundSignificant)

# Group the data by foreground/background:
groupedResults <- group_by(allResults, ForegroundBackground)
groupedStats <- summarize(groupedResults,
          meanpValue = mean(`Corrected P-value`),
          standardDeviationpValue = sd(`Corrected P-value`))
