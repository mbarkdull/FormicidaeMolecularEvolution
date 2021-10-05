library(RJSONIO)
library(tidyverse)
library(ggthemes)

# Get a list of results files, sort by decreasing size, and remove any that are empty:
files <- list.files(path = "./BUSTEDPH", pattern = "*.json", full.names = TRUE)
files <- sort(files, decreasing = TRUE)
files <- files[sapply(files, file.size) > 0]

#Write a function to read a single BUSTEDPH result file and extract relevant information from the JSON:
parsingBustedPH <- function(content) {
  result <- RJSONIO::fromJSON(content = content)
  if (result[["test results"]][["p-value"]] <= 0.05) {
    if (result[["test results background"]][["p-value"]] > 0.05) {
      if (result[["test results shared distributions"]][["p-value"]] <= 0.05) {
        print("Selection is associated with the phenotype / trait")
        textResult <- "Selection is associated with the phenotype / trait"
        data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "Yes")
        return(data)
      } else {
        print("Selection is associated with the phenotype / trait, but there is no significant difference between test and background branches in terms of selective pressure")
        textResult <- "Selection is associated with the phenotype / trait, but there is no significant difference between test and background branches in terms of selective pressure"
        data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "No")
      }
    } else {
      #print("Selection is acting on the branches with the phenotype / trait, but is **also** acting on background branches.")
      #textResult <- "Selection is acting on the branches with the phenotype / trait, but is **also** acting on background branches."
      #data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]])
      if (result[["test results shared distributions"]][["p-value"]] <= 0.05) {
        print("There is a significant difference between test and background branches in terms of selective pressure")
        textResult <- "There is a significant difference between test and background branches in terms of selective pressure"
        data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "Yes")
      } else {
        print("There is no significant difference between test and background branches in terms of selective pressure")
        textResult <- "There is no significant difference between test and background branches in terms of selective pressure"
        data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "No")
      }
    }
  } else {
    print("There is **no evidence** of episodic diversifying selection on test branches; selection is not associated with phenotype/trait")
    textResult <- "There is **no evidence** of episodic diversifying selection on test branches; selection is not associated with phenotype/trait"
    data <- c(result[["input"]][["file name"]], textResult, result[["test results"]][["p-value"]], result[["test results background"]][["p-value"]], result[["test results shared distributions"]][["p-value"]], "No")
  }
  
}
# Make a safer version of that function with `possibly`:
possiblyparsingBustedPH <- possibly(parsingBustedPH, otherwise = "File empty.")
# Map the function over all results files to construct a master dataframe:
workerPolymorphismBustedPHResults <- map(files, possiblyparsingBustedPH)

# Write a function to process and fix the column names of the results dataframe:
bustedPHDataframeProcessing <- function(resultsDataframe) {
  BustedPHResults <- as.data.frame(do.call(rbind, resultsDataframe))   
  colnames(BustedPHResults) <- c("file", "textResult", "test results p-value", "test results background p-value", "test results shared distributions p-value", "differenceInSelection")
  orthogroup <- sapply(strsplit(as.character(BustedPHResults$file),"/"), tail, 1)
  orthogroup <- sapply(strsplit(orthogroup, "_"), `[`, 2)
  BustedPHResults$orthogroup <- orthogroup
  trait <- sapply(strsplit(as.character(BustedPHResults$file),"/"), tail, 1)
  trait <- sapply(strsplit(trait, "_"), `[`, 1)
  BustedPHResults$trait <- trait
  return(BustedPHResults)
}
workerPolymorphismBustedPHResults <- bustedPHDataframeProcessing(workerPolymorphismBustedPHResults)

workerPolymorphismBustedPHResults <- workerPolymorphismBustedPHResults %>% mutate(selectionOn =
                     case_when(as.numeric(as.character(`test results p-value`)) <= 0.05 & 
                                 as.numeric(as.character(`test results background p-value`)) > 0.05 &
                                 as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 ~ "ForegroundOnly", 
                                as.numeric(as.character(`test results p-value`)) > 0.05 &
                                   as.numeric(as.character(`test results background p-value`)) <= 0.05 &
                                   as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 ~ "BackgroundOnly", 
                               as.numeric(as.character(`test results p-value`)) > 0.05 &
                                 as.numeric(as.character(`test results background p-value`)) > 0.05  ~ "NoSelection",
                               as.numeric(as.character(`test results shared distributions p-value`)) > 0.05 ~ "NoSignificantDifferenceBetweenForegroundAndBackground"))

# Write a function to merge the BUSTEDPH results and the annotations from InterProScan; add all the annotations:
mergeAnnotations <- function(annotationFile, data) {
  annotations <- read_delim(annotationFile, delim = "\t")
  # Convert them to wide format:
  annotations <- select(annotations, "#cluster_id", "domain_id", "domain_description") %>%
    group_by(`#cluster_id`) %>%
    mutate(uniqueID = row_number()) %>%
    tidyr::pivot_wider(names_from = uniqueID, values_from = c(domain_description, domain_id)) 
  data <- left_join(data, annotations, by = c("orthogroup" = "#cluster_id"))
}
# Read in annotations:
data <- mergeAnnotations(annotationFile = "cluster_domain_annotation.Pfam.txt", data = workerPolymorphismBustedPHResults)
data <- mergeAnnotations(annotationFile = "cluster_domain_annotation.IPR.txt", data = data)
data <- mergeAnnotations(annotationFile = "cluster_domain_annotation.GO.txt", data = data)

# Filter to only positively selected genes and then pivot longer for readability:
positivelySelectedGenes <- filter(data, differenceInSelection == "Yes") %>%
  pivot_longer(cols = starts_with('domain'), 
                     names_to = c('domain', 'annotation'), 
                     names_sep = "_") %>%
  filter(!is.na(value))


ggplot(data, 
       aes(x = selectionOn)) +
  geom_bar() + 
  theme(axis.text.x = element_text(angle = 0)) + 
  scale_x_discrete(labels =c("NoSelection" = "No selection \non foreground \nor background",
                             "NoSignificantDifferenceBetweenForegroundAndBackground" = "Foreground and \nbackground regimes \nare not \nsignificantly \ndifferent",
                             "BackgroundOnly" = "Selection \non background \nonly",
                             "ForegroundOnly" = "Selection \non foreground \nonly")) +
  labs(x = "Selective regime",
       y = "Count of orthogroups",
       title = "Distribution of selective regimes associated with trait")  +
  theme_hc()


backgroundSelected <- data %>% 
  filter(as.numeric(as.character(`test results p-value`)) > 0.05 & 
           as.numeric(as.character(`test results background p-value`)) <= 0.05 & 
           as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05)
nBackgroundSelected <- nrow(backgroundSelected)

foregroundSelected <- data %>% 
  filter(as.numeric(as.character(`test results p-value`)) <= 0.05 & 
           as.numeric(as.character(`test results background p-value`)) > 0.05 & 
           as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05)
nForegroundSelected <- nrow(foregroundSelected)

frequencyTableWithNoCorrection <- matrix(c(nBackgroundSelected, nForegroundSelected, 
                                           (nrow(data) - nBackgroundSelected),
                                           (nrow(data) - nForegroundSelected)),
                                         nrow = 2,
                                         byrow = TRUE,
                                         dimnames = list("category" = c("selected", "noSelection") ,
                                                         "selection" = c("backgroundSelected", "foregroundSelected")))
fishersExactTest <- fisher.test(frequencyTableWithNoCorrection)

fishersExactTest[["p.value"]]
if (fishersExactTest[["p.value"]] <= 0.05) {
  print("Signficant difference between foreground and background") 
} else {
  print("No signficant difference.")
}




