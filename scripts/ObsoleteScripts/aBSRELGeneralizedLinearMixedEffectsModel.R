# Load necessary packages (I'm using RJSONIO here because rjson threw errors for one of my results files):
library(RJSONIO)
library(tidyverse)
library(janitor)
library(tidytext)
library(lme4)
library(emmeans)

# Write a function to read in the aBSREL JSON result file and transform it into a dataframe:
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

# Get a binary column that classifies each gene as selected or not. 
allResults <- allResults %>%
  mutate(binarySelection = case_when(
    `Corrected P-value` >= 0.05 ~ "NotSelected",
    `Corrected P-value` <= 0.05 ~ "Selected"
  ))

# Convert that column to a factor, so that glmer will be happy with it:
allResults$binarySelection <- as.factor(allResults$binarySelection)
levels(allResults$binarySelection)

allResultsExport <- select(allResults, c(`original name`, `Corrected P-value`, binarySelection, Orthogroup, Species, ForegroundBackground))
write.csv(allResultsExport, "./allResults.csv", row.names = FALSE)

# Run a generalized linear mixed effects model:
  # binarySelection is my response
  # the predictors are on the right of the tilde; the terms with a vertical bar are random-effects terms, because we need to account for the effect of species and orthogroup (in other words, I don't truly have nrows(allResults) truly independent observations)
  # (1|group)	is the term for the random group intercept
Model <- glmer(binarySelection ~ ForegroundBackground + (1|Species) + (1|Orthogroup), 
               data = allResults, 
               family =  "binomial")

# Look under fixed effects. "Estimate" for ForegroundBackgroundforeground is -0.17206, which tells me foreground species are less likely to be selected, because of the negative value. The Pr(>|z|)  tells me my result is significant. https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-glm-and-lme4/
summary(Model)


contingency <- table(allResults$ForegroundBackground, allResults$binarySelection)
# Gives the row proportions of the contingency table. 
prop.table(contingency, margin = 1)
# Take the model, figure out predicted values for groups, and back transform to the prob
emmeans(Model, ~ForegroundBackground, type = "response")

plot(emmeans(Model, ~ForegroundBackground, type = "response"))





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

a <- by(groupedResults, groupedResults$ForegroundBackground, summary)
a <-prop.table(table(groupedResults$ForegroundBackground, groupedResults$realProportion))







