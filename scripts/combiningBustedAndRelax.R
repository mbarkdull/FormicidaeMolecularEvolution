# Load required packages:
library(tidyverse)
library(ggthemes)
library(gt)
library(plyr)
library(rjson)
library(MetBrewer)
library(SuperExactTest)
library(topGO)

##### Analyze results for all orthogroups: ##### 
# Read in results:
relaxResults <- read_csv("./Results/allRelaxResults.csv", col_names = TRUE)
bustedPHResults <- read_csv("./Results/allBustedPHResults.csv", col_names = TRUE)

# Combine relax and busted results:
relaxAndBustedPH <- full_join(relaxResults, 
                  bustedPHResults, 
                  by = c("orthogroup" = "orthogroup", 
                         "trait" = "trait"))

relaxAndBustedPH$combo <- paste(relaxAndBustedPH$selectionOn, relaxAndBustedPH$selectionCategory)
write.csv(relaxAndBustedPH, file = "relaxAndBustedPH.csv")

relaxAndBustedPHPlot <- ggplot(data = filter(relaxAndBustedPH,
                                             file != "File empty." & !is.na(trait) & !is.na(file))) + 
  geom_bar(mapping = aes(x = combo),
           position = "dodge") +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  coord_flip() + 
  scale_x_discrete(labels=c("EvidenceOfSelectionAssociatedWithLackOfTraitButNS nonsignficantRelaxation" = "nonsignficant background positive, nonsignificant foreground relaxation",      
                            "NoEvidenceOfSelection nonsignficantRelaxation" = "no positive selection, nonsignficant foreground relaxation",                                 
                            "EvidenceOfSelectionAssociatedWithLackOfTraitButNS signficantIntensification" = "nonsignificant background positive, significant foreground intensification",
                            "BackgroundOnly signficantIntensification" = "significant background positive selection, significant foreground intensification",
                            "NoEvidenceOfSelection signficantIntensification" = "no positive selection, significant foreground intensification",
                            "SelectionOnBothButNoSignificantDifference nonsignficantIntensification" = "similar positive selection on both, nonsignificant foreground intensification",
                            "EvidenceOfSelectionAssociatedWithTraitButNS signficantIntensification" = "nonsignificant foreground positive selection, significant foreground intensification",
                            "EvidenceOfSelectionAssociatedWithTraitButNS nonsignficantRelaxation" = "nonsignificant foreground positive selection, nonsignificant foreground relaxation",
                            "NoEvidenceOfSelection nonsignficantIntensification" = "no positive selection, nonsignificant foreground intensification",
                            "EvidenceOfSelectionAssociatedWithLackOfTraitButNS nonsignficantIntensification" = "nonsignificant foreground positive selection, nonsignificant foreground intensification",
                            "ForegroundOnly nonsignficantRelaxation" = "significant foreground positive selection, nonsignificant foreground relaxation",
                            "NoEvidenceOfSelection signficantRelaxation" = "no positive selection, significant foreground relaxation",
                            "EvidenceOfSelectionAssociatedWithTraitButNS nonsignficantIntensification" = "nonsignificant foreground positive selection, nonsignificant foreground relaxation",
                            "BackgroundOnly nonsignficantIntensification" = "significant background positive selection, nonsignificant foreground intensification",
                            "BackgroundOnly nonsignficantRelaxation" = "significant background positive selection, nonsignificant foreground relaxation",                                         
                            "BackgroundOnly signficantRelaxation" = "significant background positive selection, significant foreground relaxation",                                           
                            "ForegroundOnly signficantRelaxation" = "significant foreground positive selection, significant foreground relaxation",                                            
                            "EvidenceOfSelectionAssociatedWithLackOfTraitButNS signficantRelaxation" = "nonsignificant background positive selection, significant foreground relaxation",        
                            "EvidenceOfSelectionAssociatedWithTraitButNS signficantRelaxation"  = "nonsignificant foreground positive selection, significant foreground relaxation",                
                            "ForegroundOnly signficantIntensification" = "significant foreground positive selection, significant foreground intensification",                                      
                            "NoEvidenceOfSelection NA" = "no positive selection, not tested with RELAX",                                                       
                            "SelectionOnBothButDifferent nonsignficantIntensification" = "positive and different selection on fore- and background, nonsignificant foreground intensification",                      
                            "ForegroundOnly nonsignficantIntensification"  = "significant foreground positive selection, nonsignficant foreground intensification",                                    
                            "SelectionOnBothButNoSignificantDifference signficantIntensification" = "positive and not different selection on fore- and background, significant foreground intensification",           
                            "SelectionOnBothButDifferent signficantIntensification" = "positive and different selection on fore- and background, significant foreground intensification",                          
                            "NA nonsignficantIntensification" = "not tested with BUSTED-PH, nonsignificant foreground intensification",                                               
                            "SelectionOnBothButNoSignificantDifference nonsignficantRelaxation" = "positive and not different selection on fore- and background, nonsignificant foreground relaxation",              
                            "SelectionOnBothButDifferent signficantRelaxation" = "positive and different selection on fore- and background, significant foreground relaxation",                              
                            "SelectionOnBothButNoSignificantDifference signficantRelaxation"  = "positive and not different selection on fore- and background, significant foreground relaxation",                  
                            "EvidenceOfSelectionAssociatedWithTraitButNS NA" = "nonsignificant foreground positive selection, not tested with RELAX",                                 
                            "SelectionOnBothButDifferent nonsignficantRelaxation"  = "positive and different selection on fore- and background, nonsignificant foreground relaxation",                             
                            "EvidenceOfSelectionAssociatedWithLackOfTraitButNS NA" = "nonsignificant background positive selection, not tested with RELAX",                           
                            "BackgroundOnly NA" = "signficant background positive selection, not tested with RELAX",                                                              
                            "SelectionOnBothButDifferent NA" = "positive and different selection on fore- and background, not tested with RELAX",                                                
                            "SelectionOnBothButNoSignificantDifference NA"  = "positive and not different selection on fore- and background, not tested with RELAX",                                    
                            "ForegroundOnly NA"  = "significant foreground positive selection, not tested with RELAX"  )) + 
  facet_wrap(~trait.x,
             nrow = 3) 
  
relaxAndBustedPHPlot

traits <- unique(relaxAndBustedPH$trait)

# Calculate p-values for each trait:
pValueByTraitRelax <- function(specificTrait, inputData) {
  relaxResultsTrait <- filter(inputData, 
                              trait == specificTrait & !is.na(kValue))
  nSelectionIntensified <- relaxResultsTrait %>% 
    filter(pValue <= 0.05,
           kValue > 1) %>%
    nrow()
  
  nSelectionRelaxed <- relaxResultsTrait %>% 
    filter(pValue <= 0.05,
           kValue < 1) %>%
    nrow()
  # Construct a frequency table with that information:
  frequencyTableWithNoCorrection <- matrix(c(nSelectionIntensified, nSelectionRelaxed, 
                                             (nrow(relaxResultsTrait) - nSelectionIntensified),
                                             (nrow(relaxResultsTrait) - nSelectionRelaxed)),
                                           nrow = 2,
                                           byrow = TRUE,
                                           dimnames = list("category" = c("evidenceForSelection", "nonsignificantResult") ,
                                                           "selection" = c("intensified", "relaxed")))
  # Run a fisher's exact test to see if there is a difference in the proportion of genes under selection in the fore- vs. background:
  fishersExactTest <- fisher.test(frequencyTableWithNoCorrection)
  # Print the resulting p-value:
  fishersExactTest[["p.value"]]
  pValue <- if (fishersExactTest[["p.value"]] < 0.000001) {
    formatC(fishersExactTest[["p.value"]], format = "e", digits = 2)
  } else {
    round(fishersExactTest[["p.value"]], digits = 4)
  }
  textHeight <- as.numeric(max(nSelectionIntensified, nSelectionRelaxed))
  return(c(specificTrait, pValue, textHeight))
}
possiblypValueByTraitRelax <- possibly(pValueByTraitRelax, otherwise = "Error")
pValuesRelax <- traits %>% 
  purrr::map(~ possiblypValueByTraitRelax(.x, relaxAndBustedPH))
pValuesRelax <- as.data.frame(do.call(rbind, pValuesRelax))   
colnames(pValuesRelax) <- c("trait", "pValue", "maxHeight")
pValuesRelax$test <- "relax"

pValueByTraitBUSTEDPH <- function(specificTrait, inputData) {
  bustedPHResultsTrait <- filter(inputData, 
                                 trait == specificTrait & !is.na(`test results p-value`))
  nSelectionForeground <- bustedPHResultsTrait %>% 
    filter(`test results p-value` <= 0.05,
           `test results background p-value` > 0.05,
           `test results shared distributions p-value` <= 0.05) %>%
    nrow()
  
  nSelectionBackground <- bustedPHResultsTrait %>% 
    filter(`test results p-value` > 0.05,
           `test results background p-value` <= 0.05,
           `test results shared distributions p-value` <= 0.05) %>%
    nrow()
  # Construct a frequency table with that information:
  frequencyTableWithNoCorrection <- matrix(c(nSelectionForeground, nSelectionBackground, 
                                             (nrow(bustedPHResultsTrait) - nSelectionForeground),
                                             (nrow(bustedPHResultsTrait) - nSelectionBackground)),
                                           nrow = 2,
                                           byrow = TRUE,
                                           dimnames = list("category" = c("selected", "noSelection") ,
                                                           "selection" = c("foreground", "background")))
  # Run a fisher's exact test to see if there is a difference in the proportion of genes under selection in the fore- vs. background:
  fishersExactTest <- fisher.test(frequencyTableWithNoCorrection)
  # Print the resulting p-value:
  fishersExactTest[["p.value"]]
  pValue <- if (fishersExactTest[["p.value"]] < 0.000001) {
    formatC(fishersExactTest[["p.value"]], format = "e", digits = 2)
  } else {
    round(fishersExactTest[["p.value"]], digits = 4)
  }
  textHeight <- as.numeric(max(nSelectionForeground, nSelectionBackground))
  return(c(specificTrait, pValue, textHeight))
}
possiblypValueByTraitBUSTEDPH <- possibly(pValueByTraitBUSTEDPH, otherwise = "Error")
pValuesBUSTEDPH <- traits %>% 
  map(~ possiblypValueByTraitBUSTEDPH(.x, relaxAndBustedPH))
pValuesBUSTEDPH <- as.data.frame(do.call(rbind, pValuesBUSTEDPH))   
colnames(pValuesBUSTEDPH) <- c("trait", "pValue", "maxHeight")
pValuesBUSTEDPH$test <- "bustedPH"

pValuesAll <- rbind(pValuesRelax, pValuesBUSTEDPH)

# Make a lollipop plot:
# Create data
summaryData <- function(specificTrait) {
  relaxResultsTrait <- filter(relaxResults, 
                              trait == specificTrait)
  nSelectionIntensified <- relaxResultsTrait %>% 
    filter(pValue <= 0.05,
           kValue > 1) %>%
    nrow()
  
  nSelectionRelaxed <- relaxResultsTrait %>% 
    filter(pValue <= 0.05,
           kValue < 1) %>%
    nrow()
  bustedPHResultsTrait <- filter(bustedPHResults, 
                                 trait == specificTrait)
  nSelectionForeground <- bustedPHResultsTrait %>% 
    filter(`test results p-value` <= 0.05,
           `test results background p-value` > 0.05,
           `test results shared distributions p-value` <= 0.05) %>%
    nrow()
  
  nSelectionBackground <- bustedPHResultsTrait %>% 
    filter(`test results p-value` > 0.05,
           `test results background p-value` <= 0.05,
           `test results shared distributions p-value` <= 0.05) %>%
    nrow()
  
  nSelectionBoth <- bustedPHResultsTrait %>% 
    filter(`test results p-value` <= 0.05,
           `test results background p-value` <= 0.05) %>%
    nrow()
  # Create a dataframe that lists these data:
  test <- c("relax", "relax", "bustedPH", "bustedPH", "bustedPH")
  selection <- c(paste(specificTrait, "nSelectionIntensified"), 
                 paste(specificTrait, "nSelectionRelaxed"), 
                 paste(specificTrait, "nSelectionForeground"), 
                 paste(specificTrait, "nSelectionBackground"),
                 paste(specificTrait, "nSelectionBoth"))
  number <- c(nSelectionIntensified, nSelectionRelaxed, nSelectionForeground, nSelectionBackground, nSelectionBoth)
  countsOfSelectionRegimes <- data.frame(test, selection, number)
  
  return(c(countsOfSelectionRegimes))
}
possiblysummaryData <- possibly(summaryData, otherwise = "Error")
summaryDataframe <- purrr::map_dfr(traits, possiblysummaryData)
summaryDataframe <- separate(summaryDataframe,
                             selection, 
                             into = c("trait", "selection"),
                             sep = " ")


summaryDataframe <- dplyr::full_join(summaryDataframe, pValuesAll, by = c("trait", "test")) 

# Create a color scheme and nicely formatted trait category labels:
colorScheme <- c("nSelectionIntensified" = "#819fe3",
                 "nSelectionRelaxed" = "#391463",
                 "nSelectionForeground" = "#fc9086",
                 "nSelectionBackground" = "#ba1f19",
                 "nSelectionBoth" = "#e0d8d7")
traitLabels <- c("polyandry" = "polyandry",
                 "polygyny" = "polygyny",
                 "workerReproductionQueens" = "reproductive workers",
                 "workerPolymorphism" = "polymorphic workers",
                 "multilineage" = "multilineage colonies")
traitLabels <- data.frame(traitLabels, trait = names(traitLabels))
summaryDataframe <- left_join(summaryDataframe, traitLabels)
rm(traitLabels)

traits <- unique(summaryDataframe$trait)
summaryDataframe$selection <- factor(summaryDataframe$selection,
                                     levels = c("nSelectionForeground",
                                                "nSelectionBackground",
                                                "nSelectionBoth",
                                                "nSelectionRelaxed",
                                                "nSelectionIntensified"))
#traits <- c("workerReproductionQueens", "workerPolymorphism")

#### Plot all results with lollipop plots: ####
lollipopPlots <- function(specificTrait, inputData) {
  testLabels <- c("bustedPH" = "Distribution of positive selection\n(BUSTED-PH)",
                  "relax" = "Distribution of selection intensity\n(RELAX)")
  plotData <- dplyr::filter(inputData, 
                            trait == specificTrait,
                            selection != "nSelectionBoth")
  plot <- ggplot(plotData, 
                 aes(x = selection,
                     y = number,
                     fill = selection)) +
    geom_point(mapping = aes(colour = selection),
               position = position_dodge2(width = 0.75),
               size = 4.5) + 
    geom_linerange(aes(x = selection,
                       ymin = 0,
                       ymax = number,
                       colour = selection),
                   position = position_dodge2(width = 0.75),
                   size = 2) +
    geom_text(data = dplyr::filter(plotData, 
                                   trait == specificTrait), 
              aes(x = 1.5,
                  y = as.numeric(as.character(maxHeight)) + 300,
                  label = pValue),
              colour = "grey50",
              size = 3) + 
    geom_linerange(data = dplyr::filter(plotData, 
                                        trait == specificTrait),
                   aes(xmin = 1, 
                       xmax = 2, 
                       y = (as.numeric(as.character(maxHeight)) + 200)), 
                   color = "grey50", 
                   size = 0.3) +
    theme_bw() +
    scale_color_manual(values = colorScheme,
                       labels = c("nSelectionIntensified" = "Selection intensified\non foreground",
                                  "nSelectionRelaxed" = "Selection relaxed\non foreground",
                                  "nSelectionForeground" = "Positive selection\non foreground",
                                  "nSelectionBackground" = "Positive selection\non background")) +
    labs(y = "Count of orthogroups", 
         title = paste(str_to_sentence(plotData$traitLabels), 
                       sep = " ")) +
    scale_x_discrete(labels = c("bustedPH" = "Positive selection\non fore- or background",
                                "relax" = "Relaxed or\nintensified selection")) + 
    facet_wrap(~test,
               nrow = 2,
               labeller = labeller(test = testLabels),
               scales = 'free_x', 
               strip.position = "top") + 
    theme(axis.text.x = element_text(angle = 0,
                                     hjust = 0.5,
                                     vjust = 1,
                                     size = 8),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.spacing = unit(0.5, "cm"),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "gray98"),
          legend.position = "none",
          panel.spacing.x = unit(0, "pt"), 
          panel.spacing.y = unit(30, "pt"),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(), 
          rect = element_rect(fill = "transparent"),
          plot.title = element_text(size = 15)) +
    scale_x_discrete(labels = c("nSelectionBackground" = "Positive selection\non species without trait",
                                "nSelectionForeground" = "Positive selection\non species with trait",
                                "nSelectionIntensified" = "Intensified selection\nassociated with trait",
                                "nSelectionRelaxed"  = "Relaxed selection\nassociated with trait"))
}
subsetTraits <- c("workerReproductionQueens", "workerPolymorphism", "multilineage")
test <- purrr::map(subsetTraits, ~lollipopPlots(.x, inputData = summaryDataframe))
cowplot::plot_grid(plotlist = test,
                   ncol = 3)
patchworkTest <- patchwork::wrap_plots(test, 
                                       ncol = 3) + 
  patchwork::plot_annotation(title = 'Impact of sociality-associated traits on molecular evolution',
                             theme = theme(plot.title = element_text(size = 19)))
patchworkTest[[1]] <- patchworkTest[[1]] +
  ylim(0, 2700)
patchworkTest[[2]] <- patchworkTest[[2]] + 
  ylim(0, 2700)
patchworkTest[[3]] <- patchworkTest[[3]] +
  ylim(0, 2700)
#patchworkTest[[4]] <- patchworkTest[[4]] + 
  #ylim(0, 2700)
#patchworkTest[[5]] <- patchworkTest[[5]] + 
  #ylim(0, 2700)
patchworkTest[[1]] <- patchworkTest[[1]] +
  theme(strip.text.x = element_text(color = "white")) 

patchworkTest[[2]] <- patchworkTest[[2]] + 
  theme(axis.text.y = element_blank(),           
        axis.ticks.y = element_blank(),        
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 11)) 

patchworkTest[[3]] <- patchworkTest[[3]] + 
  theme(axis.text.y = element_blank(),           
        axis.ticks.y = element_blank(),        
        axis.title.y = element_blank(),
        strip.text.x = element_text(color = "white")) 

patchworkTest

ggsave(filename = "comparativeGenomicsResultsLollipop.png", 
       device = "png",
       path = "./Plots/", 
       bg = "transparent", 
       width = 9, 
       height = 8)

for (i in test) {
  trait <- unique(i[["data"]][["trait"]])
  filename <- paste(trait, "ComparativeGenomicsResultsLollipop.png", sep = "")
  plot(i)
  print(filename)
  ggsave(filename = filename, path = "./Plots/", bg = "transparent", width = 7, height = 5)
}

#### Construct euler diagrams (like Venn diagrams but doesn't have to show all possible intersections) ####
library(eulerr)

eulerPlot <- function(specificTrait) {
  intensifiedOrthogroups <- relaxAndBustedPH %>%
    filter(trait  == specificTrait) %>%
    filter(pValue <= 0.05 & kValue > 1) %>%
    pull(orthogroup)
  
  relaxedOrthogroups <- relaxAndBustedPH %>% 
    filter(trait  == specificTrait) %>%
    filter(pValue <= 0.05 & kValue < 1) %>%
    pull(orthogroup)
  
  foregroundPositiveOrthogroups <- relaxAndBustedPH %>% 
    filter(trait  == specificTrait) %>%
    filter(selectionOn == "ForegroundOnly") %>%
    pull(orthogroup)
  
  backgroundPositiveOrthogroups <- relaxAndBustedPH %>% 
    filter(trait  == specificTrait) %>%
    filter(selectionOn == "BackgroundOnly") %>%
    pull(orthogroup)
  
  allOrthogroups <- relaxAndBustedPH %>% 
    filter(trait  == specificTrait) %>%
    pull(orthogroup) %>%
    unique()
  
  selectionSets <- list(`Intensified selection` = intensifiedOrthogroups,
                        `Relaxed selection` = relaxedOrthogroups,
                        `Positive selection\nassociated\nwith trait` = foregroundPositiveOrthogroups,
                        `Positive selection\nassociated\nwith lack of trait` = backgroundPositiveOrthogroups)
  
  intersectionsOfSelection <- euler(selectionSets,
                                    shape = "ellipse")
  #library(Hmisc)
  # If you want to capitalize the first letter of something, use Hmisc::capitalize()
  #eulerr_options(padding = grid::unit(0, "lines"))
  test <- plot(intersectionsOfSelection,
          quantities = list(fontsize = 5,
                            labels = paste(intersectionsOfSelection$original, "\northogroups", sep = "")),
          labels = list(fontsize = 8),
          adjust_labels = TRUE,
          main = list(label = paste("Combinations of selective regimes\nrelated to", tolower(gsub("([A-Z])", " \\1", specificTrait)), sep = " "),
                      fontsize = 7),
          edges = list(col = "grey50"),
          fills = c("#dce7f5", 
                    "#b6a5c9", 
                    "#f5ccc9", 
                    "#d1726f")) 
  test
}

eulerPlotsList <- purrr::map(traits, ~eulerPlot(.x))
eulerPatchwork <- patchwork::wrap_plots(eulerPlotsList,
                                        nrow = 2)
eulerPatchwork
ggsave(filename = "eulerPatchwork.svg", 
       device = "svg",
       path = "./Plots/", 
       bg = "transparent", 
       width = 12, 
       height = 18)

#### Test for differences in the number of orthogroups under BOTH intensified and positive selection in the foreground vs the background: ####
differencesInPurifyingSelection <- function(specificTrait, inputData) {
  resultsTrait <- filter(inputData,
                        trait == specificTrait & !is.na(`test results p-value`) & !is.na(pValue))
  nPurifyingForeground <- resultsTrait %>% 
    filter(`test results p-value` <= 0.05,
           `test results background p-value` > 0.05,
           `test results shared distributions p-value` <= 0.05,
           kValue > 1, 
           pValue <= 0.05) %>%
    nrow()
  
  nPurifyingBackground <- resultsTrait %>% 
    filter(`test results p-value` > 0.05,
           `test results background p-value` <= 0.05,
           `test results shared distributions p-value` <= 0.05,
           kValue < 1,
           pValue <= 0.05) %>%
    nrow()
  # Construct a frequency table with that information:
  frequencyTableWithNoCorrection <- matrix(c(nPurifyingForeground, nPurifyingBackground, 
                                             (nrow(resultsTrait) - nPurifyingForeground),
                                             (nrow(resultsTrait) - nPurifyingBackground)),
                                           nrow = 2,
                                           byrow = TRUE,
                                           dimnames = list("category" = c("purifying", "noPurifying") ,
                                                           "selection" = c("foreground", "background")))
  # Run a fisher's exact test to see if there is a difference in the proportion of genes under selection in the fore- vs. background:
  fishersExactTest <- fisher.test(frequencyTableWithNoCorrection)
  # Print the resulting p-value:
  fishersExactTest[["p.value"]]
  pValue <- if (fishersExactTest[["p.value"]] < 0.000001) {
    formatC(fishersExactTest[["p.value"]], format = "e", digits = 2)
  } else {
    round(fishersExactTest[["p.value"]], digits = 4)
  }
  textHeight <- as.numeric(max(nPurifyingForeground, nPurifyingBackground))
  return(c(specificTrait, pValue, textHeight, nPurifyingForeground, nPurifyingBackground))
}

possiblyDifferencesInPurifyingSelection <- possibly(differencesInPurifyingSelection, otherwise = "Error")

pValuesPurifying <- traits %>% 
  map(~ possiblyDifferencesInPurifyingSelection(.x, relaxAndBustedPH))

pValuesPurifying <- as.data.frame(do.call(rbind, pValuesPurifying))   

colnames(pValuesPurifying) <- c("trait", "pValue", "maxHeight", "nPurifyingForeground", "nPurifyingBackground")

# Construct a list of orthogroups under positive selection in the foreground:
foregroundPositiveOrthogroups <- relaxAndBustedPH %>% 
  filter(trait  == specificTrait) %>%
  filter(selectionOn == "ForegroundOnly") %>%
  pull(orthogroup)

# Construct a list of orthogroups under positive selection in the background:
backgroundPositiveOrthogroups <- relaxAndBustedPH %>% 
  filter(trait  == specificTrait) %>%
  filter(selectionOn == "BackgroundOnly") %>%
  pull(orthogroup)

# Construct a list of orthogroups under intensified selection in the foreground:
foregroundIntensifiedOrthogroups <- relaxAndBustedPH %>% 
  filter(trait  == specificTrait) %>%
  filter(shortDescription == "Intensification of selection along foreground branches") %>%
  pull(orthogroup)

# Construct a list of orthogroups under intensified selection in the background:
backgroundIntensifiedOrthogroups <- relaxAndBustedPH %>% 
  filter(trait  == specificTrait) %>%
  filter(shortDescription == "Relaxation of selection along foreground branches") %>%
  pull(orthogroup)

selectionSets <- list(`Foreground positive selection` = foregroundPositiveOrthogroups,
                      `Background positive selection` = backgroundPositiveOrthogroups,
                      `Foreground intensified selection` = foregroundIntensifiedOrthogroups,
                      `Background intensified selection` = backgroundIntensifiedOrthogroups)
totalNumberOfGenes <- length(which(relaxAndBustedPH$trait == specificTrait))

test <- supertest(selectionSets, 
          n = totalNumberOfGenes,
          degree = c(2))
summary(test)

##### Do these same analyses but only for the single-copy orthogroups: #####
# Get a list of the single copy orthogroups:
singleCopyOrthogroups <- list.files(path = paste("./5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Single_Copy_Orthologue_Sequences/", sep = ""), full.names = FALSE)
singleCopyOrthogroups <- sapply(strsplit(as.character(singleCopyOrthogroups), "[.]"), head, 1)
singleCopyOrthogroups <- read_csv("./Results/singleCopyOrthogroups.csv")
singleCopyOrthogroups <- singleCopyOrthogroups$x

singleCopyResults <- dplyr::filter(relaxAndBustedPH, 
                                   orthogroup %in% singleCopyOrthogroups)


# Get p-values for the relax results for each trait:
pValuesRelaxSingle <- traits %>% 
  purrr::map(~ possiblypValueByTraitRelax(.x, singleCopyResults))
pValuesRelaxSingle <- as.data.frame(do.call(rbind, pValuesRelaxSingle))   
colnames(pValuesRelaxSingle) <- c("trait", "pValue", "maxHeight")
pValuesRelaxSingle$test <- "relax"

# Get p-values for the BUSTED-PH results for each trait:
pValuesBUSTEDPHSingle <- traits %>% 
  map(~ possiblypValueByTraitBUSTEDPH(.x, singleCopyResults))
pValuesBUSTEDPHSingle <- as.data.frame(do.call(rbind, pValuesBUSTEDPHSingle))   
colnames(pValuesBUSTEDPHSingle) <- c("trait", "pValue", "maxHeight")
pValuesBUSTEDPHSingle$test <- "bustedPH"

# Combine the p-values:
pValuesAllSingle <- rbind(pValuesRelaxSingle, pValuesBUSTEDPHSingle)

# Get summary data for the single copy orthogroups:
summaryData <- function(specificTrait, inputData) {
  relaxResultsTrait <- filter(inputData, 
                              trait == specificTrait)
  nSelectionIntensified <- relaxResultsTrait %>% 
    filter(pValue <= 0.05,
           kValue > 1) %>%
    nrow()
  
  nSelectionRelaxed <- relaxResultsTrait %>% 
    filter(pValue <= 0.05,
           kValue < 1) %>%
    nrow()
  bustedPHResultsTrait <- filter(inputData, 
                                 trait == specificTrait)
  nSelectionForeground <- bustedPHResultsTrait %>% 
    filter(`test results p-value` <= 0.05,
           `test results background p-value` > 0.05,
           `test results shared distributions p-value` <= 0.05) %>%
    nrow()
  
  nSelectionBackground <- bustedPHResultsTrait %>% 
    filter(`test results p-value` > 0.05,
           `test results background p-value` <= 0.05,
           `test results shared distributions p-value` <= 0.05) %>%
    nrow()
  # Create a dataframe that lists these data:
  test <- c("relax", "relax", "bustedPH", "bustedPH")
  selection <- c(paste(specificTrait, "nSelectionIntensified"), 
                 paste(specificTrait, "nSelectionRelaxed"), 
                 paste(specificTrait, "nSelectionForeground"), 
                 paste(specificTrait, "nSelectionBackground"))
  number <- c(nSelectionIntensified, nSelectionRelaxed, nSelectionForeground, nSelectionBackground)
  countsOfSelectionRegimes <- data.frame(test, selection, number)
  
  return(c(countsOfSelectionRegimes))
}
possiblysummaryData <- possibly(summaryData, otherwise = "Error")

summaryDataframeSingle <- traits %>% 
  map_dfr(~ possiblysummaryData(.x, singleCopyResults))

summaryDataframeSingle <- separate(summaryDataframeSingle,
                             selection, 
                             into = c("trait", "selection"),
                             sep = " ")

summaryDataframeSingle <- dplyr::full_join(summaryDataframeSingle, pValuesAllSingle, by = c("trait", "test")) 
traitLabels <- c("polyandry" = "polyandry",
                 "polygyny" = "polygyny",
                 "workerReproductionQueens" = "reproductive workers",
                 "workerPolymorphism" = "polymorphic workers",
                 "multilineage" = "multilineage colonies")
traitLabels <- data.frame(traitLabels, trait = names(traitLabels))
summaryDataframeSingle <- left_join(summaryDataframeSingle, traitLabels)
rm(traitLabels)

# Make lollipop plots:
test <- purrr::map(traits, ~lollipopPlots(.x, inputData = summaryDataframeSingle))
cowplot::plot_grid(plotlist = test,
                   nrow = 2)

#### Plot the distributons of k-values: ####
# Define a new axis transformation that does log10(x + 1) to the k values, so that values of zero are not dropped.
log10plusone_trans <- function() scales::trans_new("one_over", function(x) log10(x + 1), function(x) (10^x) -1)

# Add a column with human-readable trait labels to the results data:
traitLabels <- factor(c("polyandry" = "Polyandry",
                        "polygyny" = "Polygyny",
                        "workerReproductionQueens" = "Reproductive workers",
                        "workerPolymorphism" = "Polymorphic workers",
                        "multilineage" = "Multilineage colonies"))
traitLabels <- data.frame(traitLabels, trait = names(traitLabels))
relaxResults <- left_join(relaxResults, traitLabels)
rm(traitLabels)

kValuePlot <- ggplot(data = filter(relaxResults,
                                   !is.na(kValue),
                                   trait %in% subsetTraits,
                                   pValue <= 0.05)) +
  geom_histogram(mapping = aes(x = as.numeric(kValue)),
                 bins = 200,
                 fill = "#655569") +
  theme_bw() +
  geom_vline(xintercept = 1) + 
  labs(title = paste("Distribution of selection intensity parameters associated with each focal trait"),
       x = "k-value",
       y = "Count of orthogroups") +
  scale_x_continuous(trans = "log10",
                     breaks = scales::log_breaks(n = 10),
                     labels = function(x) format(x, scientific = FALSE)) +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1,
                                   vjust = 1,
                                   size = 6),
        panel.spacing = unit(0.5, "cm"),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "gray98"),
        legend.position = "none",
        panel.spacing.x = unit(0, "pt"), 
        panel.spacing.y = unit(0, "pt"),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        rect = element_rect(fill = "transparent"),
        plot.title = element_text(size = 10),
        strip.text = element_text(hjust = 0)) +
  facet_wrap(~factor(traitLabels, 
                     levels = c('Reproductive workers',
                                'Polymorphic workers',
                                'Multilineage colonies')),
             nrow = 3)

kValuePlot

ggsave(filename = "kValueHistogramsPVALUE.png", 
       device = "png",
       path = "./Plots/", 
       bg = "transparent", 
       width = 8, 
       height = 8)

kValueViolinPlot <- ggplot(data = filter(relaxResults,
                                   !is.na(kValue),
                                   trait %in% subsetTraits,
                                   pValue <= 0.05)) +
  geom_violin(mapping = aes(x = traitLabels, 
                            y = as.numeric(kValue)),
              fill = "#efebf0") +
  geom_hline(yintercept = 1) + 
  labs(title = paste("Distribution of selection intensity parameters associated with each focal trait"),
       x = "Focal sociality-associated trait",
       y = "K-values") +
  scale_y_continuous(trans = "log10",
                     breaks = scales::log_breaks(n = 10),
                     labels = function(x) format(x, scientific = FALSE)) + 
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   vjust = 1,
                                   size = 6),
        panel.spacing = unit(0.5, "cm"),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "gray98"),
        legend.position = "none",
        panel.spacing.x = unit(0, "pt"), 
        panel.spacing.y = unit(0, "pt"),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        rect = element_rect(fill = "transparent"),
        plot.title = element_text(size = 13),
        strip.text = element_text(hjust = 0))

kValueViolinPlot

ggsave(filename = "kValueViolins.png", 
       device = "png",
       path = "./Plots/", 
       bg = "transparent", 
       width = 8, 
       height = 8)

# Make a flipped violin plot:
kValueViolinPlotFlipped <- ggplot(data = filter(relaxResults,
                                         !is.na(kValue),
                                         trait %in% subsetTraits,
                                         pValue <= 0.05)) +
  geom_violin(mapping = aes(x = traitLabels, 
                            y = as.numeric(kValue)),
              fill = "#efebf0") +
  geom_hline(yintercept = 1) + 
  labs(title = paste("Distribution of selection intensity parameters associated with each focal trait"),
       x = "Focal sociality-associated trait",
       y = "K-values") +
  scale_y_continuous(trans = "log10",
                     breaks = scales::log_breaks(n = 10),
                     labels = function(x) format(x, scientific = FALSE)) + 
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0.5, "cm"),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "gray98"),
        legend.position = "none",
        panel.spacing.x = unit(0, "pt"), 
        panel.spacing.y = unit(0, "pt"),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_blank(), 
        rect = element_rect(fill = "transparent"),
        plot.title = element_text(size = 13),
        strip.text = element_text(hjust = 0)) +
  facet_wrap(~factor(traitLabels, 
                     levels = c('Reproductive workers',
                                'Polymorphic workers',
                                'Multilineage colonies')),
             ncol = 3,
             scales = 'free_y', 
             strip.position = "top")

kValueViolinPlotFlipped

#### Plot lollipop and violin plots together: ####
lollipopViolin <- function(specificTrait) {
  plotData <- dplyr::filter(summaryDataframe, 
                            trait == specificTrait,
                            selection != "nSelectionBoth",
                            test == "relax")
  lollipop <- ggplot(data = dplyr::filter(summaryDataframe, 
                                          selection != "nSelectionBoth",
                                          test == "relax",
                                          trait == specificTrait),
                     mapping = aes(x = trait,
                                   y = number,
                                   fill = selection)) +
    geom_point(mapping = aes(colour = selection),
               position = position_dodge2(width = 0.75),
               size = 4.5) + 
    geom_linerange(aes(x = trait,
                       ymin = 0,
                       ymax = number,
                       colour = selection),
                   position = position_dodge2(width = 0.75),
                   size = 2) +
    geom_text(data = dplyr::filter(plotData,
                                   trait == specificTrait), 
              aes(x = 1,
                  y = as.numeric(as.character(maxHeight)) + 300,
                  label = pValue),
              colour = "grey50",
              size = 3) + 
    geom_linerange(data = dplyr::filter(plotData,
                                        trait == specificTrait),
                   aes(xmin = 0.75, 
                       xmax = 1.25, 
                       y = (as.numeric(as.character(maxHeight)) + 200)), 
                   color = "grey50", 
                   size = 0.3) +
    scale_color_manual(values = colorScheme,
                       labels = c("nSelectionIntensified" = "Selection intensified\non foreground",
                                  "nSelectionRelaxed" = "Selection relaxed\non foreground",
                                  "nSelectionForeground" = "Positive selection\non foreground",
                                  "nSelectionBackground" = "Positive selection\non background")) +
    labs(y = "Count of orthogroups", 
         title = "Ratios of relaxed to intensified selection") + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.spacing = unit(0.5, "cm"),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "gray98"),
          legend.position = "none",
          panel.spacing.x = unit(0, "pt"), 
          panel.spacing.y = unit(30, "pt"),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(), 
          rect = element_rect(fill = "transparent"),
          plot.title = element_text(size = 15)) +
    facet_wrap(~factor(traitLabels, 
                       levels = c('reproductive workers',
                                  'polymorphic workers',
                                  'multilineage colonies')),
               ncol = 3,
               scales = 'free_y', 
               strip.position = "top")
  
  violin <- ggplot(data = filter(relaxResults,
                                 !is.na(kValue),
                                 trait == specificTrait,
                                 pValue <= 0.05)) +
    geom_violin(mapping = aes(x = traitLabels, 
                              y = as.numeric(kValue)),
                fill = "#efebf0") +
    geom_hline(yintercept = 1) + 
    labs(title = paste("Distribution of selection intensity parameters\nassociated with each focal trait"),
         x = "Focal sociality-associated trait",
         y = "K-values") +
    scale_y_continuous(trans = "log10",
                       breaks = scales::log_breaks(n = 10),
                       labels = function(x) format(x, scientific = FALSE),
                       limits = c(0.001, 400)) + 
    coord_flip() +
    theme(axis.text.y = element_text(color = "white",
                                     angle = 90),
          axis.ticks.y = element_line(color = "white"),
          panel.spacing = unit(0.5, "cm"),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "gray98"),
          legend.position = "none",
          panel.spacing.x = unit(0, "pt"), 
          panel.spacing.y = unit(0, "pt"),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.y = element_text(color = "white"), 
          rect = element_rect(fill = "transparent"),
          plot.title = element_text(size = 13),
          strip.text = element_text(hjust = 0)) +
    facet_wrap(~factor(traitLabels, 
                       levels = c('Reproductive workers',
                                  'Polymorphic workers',
                                  'Multilineage colonies')),
               ncol = 3,
               scales = 'free_y', 
               strip.position = "top") 
  
  ggarrange(lollipop,
            violin,
            ncol = 1)
}

lollipopViolin(specificTrait = "multilineage")
test <- purrr::map(subsetTraits, lollipopViolin)
patchworkTest <- cowplot::plot_grid(plotlist = test,
                   ncol = 3)

patchworkTest
# Plot the distribution of k-values and color them by whether they are under positive selection in foreground or background species:
ggplot(data = filter(relaxAndBustedPH,
                     trait == "polyandry")) +
  geom_histogram(mapping = aes(x = as.numeric(kValue),
                               fill = selectionOn),
                 bins = 60) +
  theme_bw() +
  labs(title = "Distribution of k-values associated with",
       x = "k-value",
       y = "Count of orthogroups") +
  scale_y_log10() +
  scale_x_log10()

# Get the mean k-value for each trait:
relaxAndBustedPH$kValue <- as.numeric(relaxAndBustedPH$kValue)

traitMeanKValues <- relaxAndBustedPH %>%
  filter(trait == "polyandry") %>%
  group_by(selectionOn) %>%
  summarise_at(vars(kValue), funs(mean(., na.rm=TRUE)))


ggplot(data = filter(relaxAndBustedPH,
                     trait == "polyandry")) +
  geom_histogram(mapping = aes(x = (as.numeric(kValue))),
                 bins = 60) +
  theme_bw() +
  labs(title = "Distribution of k-values associated with",
       x = "k-value",
       y = "Count of orthogroups") +
  facet_wrap(vars(selectionOn),
             ncol = 1) +
  geom_vline(data = traitMeanKValues,
               mapping = aes(xintercept = kValue),
             col = "red") 

ggplot(data = filter(relaxAndBustedPH,
                     trait == "polyandry")) +
  geom_boxplot(mapping = aes(x = selectionOn, 
                            y = kValue)) +
  scale_y_log10()

##### Look at whether common genes are under positive selection in multiple traits (specifically polyandry, polygyny, and multilineage): #####
polyandryOrthogroups <- bustedPHResults %>%
  filter(trait  == "polyandry") %>%
  filter(selectionOn == "ForegroundOnly") %>%
  pull(orthogroup)

polygynyOrthogroups <- bustedPHResults %>% 
  filter(trait  == "polygyny") %>%
  filter(selectionOn == "ForegroundOnly") %>%
  pull(orthogroup)

multilineageOrthogroups <- bustedPHResults %>% 
  filter(trait  == "multilineage") %>%
  filter(selectionOn == "ForegroundOnly") %>%
  pull(orthogroup)

allOrthogroups <- bustedPHResults %>%
  pull(orthogroup) %>%
  unique()

selectionSets <- list(`Positive selection associated with\npolyandry` = polyandryOrthogroups,
                      `Positive selection associated with\npolygyny` = polygynyOrthogroups,
                      `Positive selection associated with\nmultilineage` = multilineageOrthogroups)

intersectionsOfSelection <- euler(selectionSets,
                                  shape = "ellipse")
#library(Hmisc)
# If you want to capitalize the first letter of something, use Hmisc::capitalize()
eulerr_options(padding = grid::unit(0.1, "lines"))
test <- plot(intersectionsOfSelection,
             quantities = list(fontsize = 15,
                               labels = paste(intersectionsOfSelection$original, "\northogroups", sep = "")),
             labels = list(fontsize = 10),
             adjust_labels = TRUE,
             #main = paste("Combinations of selective regimes related to", tolower(gsub("([A-Z])", " \\1", specificTrait)), sep = " "),
             edges = list(col = "grey50"),
             fills = c("#dce7f5", 
                       "#b6a5c9", 
                       "#f5ccc9")) 
test

#### GO Term enrichment for the intersection of selection categories: ####
# The four selection categories I'm interested in are purifying selection in the foreground, purifying selection in the background, relaxed-then-positive selection in the foreground, and relaxed-then-positive selection in the background. 

# Read in the GO term annotations for each orthogroup, which we obtained from KinFin:
GOannotations <- read_delim("cluster_domain_annotation.GO.txt", delim = "\t")
# Subset so we just have orthogroup name and GO domain IDs:
longAnnotations <- dplyr::select(GOannotations, `#cluster_id`, domain_id)
# Take out genes without GO terms
longAnnotations <- longAnnotations[which(longAnnotations$domain_id != ""),] 
# Rename the column #cluster_id to orthogroup
longAnnotations <- longAnnotations %>%
  dplyr::rename(orthogroup = `#cluster_id`)

# Create list with element for each gene, containing vectors with all terms for each gene
wideListAnnotations <- tapply(longAnnotations$domain_id, longAnnotations$orthogroup, function(x)x)

goTermsPurifyingForeground <- function(specificTrait) {
   # Define vector that is 1 if gene is under purifying selection in the foreground and 0 otherwise:
  significanceInfo <- dplyr::select(relaxAndBustedPH, orthogroup, pValue, kValue, `test results p-value`, `test results background p-value`, `test results shared distributions p-value`, trait) %>%
    filter(trait == specificTrait)
  # Set each gene to 1 if under purifying selection in the foreground and 0 otherwise:
  tmp <- ifelse(significanceInfo$pValue <= 0.05 &
                  significanceInfo$kValue > 1 &
                  significanceInfo$`test results p-value` <= 0.05 &
                  significanceInfo$`test results background p-value` > 0.05 &
                  significanceInfo$`test results shared distributions p-value` <= 0.05,
                1, 
                0)
  geneList <- tmp
  rm(tmp)
  
  # Give geneList names:
  names(geneList) <- significanceInfo$orthogroup
  
  # Create the GOdata object:
  GOdataBP <- new("topGOdata",
                  ontology = "BP",
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
  # Run Fisher's exact test to check for enrichment:
  resultsFisherBP <- runTest(GOdataBP, algorithm = "elim", statistic = "fisher")
  resultsFisherBP
  resultsFisherBPTable <- GenTable(GOdataBP, raw.p.value = resultsFisherBP, topNodes = length(resultsFisherBP@score),
                                   numChar = 120)
  head(resultsFisherBPTable)
  
  GOdataMF <- new("topGOdata",
                  ontology = "MF",
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
  # Run Fisher's exact test to check for enrichment:
  resultFisherMF <- runTest(GOdataMF, algorithm = "elim", statistic = "fisher")
  resultFisherMF
  resultsFisherMFTable <- GenTable(GOdataMF, raw.p.value = resultFisherMF, topNodes = length(resultFisherMF@score),
                                   numChar = 120)
  head(resultsFisherMFTable)
  
  GOdataCC <- new("topGOdata",
                  ontology = "CC",
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
  # Run Fisher's exact test to check for enrichment:
  resultFisherCC <- runTest(GOdataCC, algorithm = "elim", statistic = "fisher")
  resultFisherCC
  resultsFisherCCTable <- GenTable(GOdataCC, raw.p.value = resultFisherCC, topNodes = length(resultFisherCC@score),
                                   numChar = 120)
  head(resultsFisherCCTable)
  
  # Combine all of the results for GO terms enriched in the foreground:
  goTermsPurifyingFore <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
    dplyr::filter(raw.p.value <= 0.01)
  
  goTermsPurifyingFore$trait <- specificTrait
  goTermsPurifyingFore
}

test <- purrr::map(traits, goTermsPurifyingForeground)
