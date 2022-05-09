# Load required packages:
library(tidyverse)
library(ggthemes)
library(gt)
library(plyr)
library(rjson)
library(MetBrewer)

relaxResults <- read_csv("./Results/allRelaxResults.csv", col_names = TRUE)
bustedPHResults <- read_csv("./Results/allBustedPHResults.csv", col_names = TRUE)


relaxAndBustedPH <- full_join(relaxResults, 
                  bustedPHResults, 
                  by = c("orthogroup" = "orthogroup", 
                         "trait" = "trait"))

relaxAndBustedPH$combo <- paste(relaxAndBustedPH$selectionOn, relaxAndBustedPH$selectionCategory)

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

#Calculate p-values for each trait:
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
                 "nSelectionBackground" = "#ba1f19")
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
                                                "nSelectionRelaxed",
                                                "nSelectionIntensified"))
#traits <- c("workerReproductionQueens", "workerPolymorphism")
lollipopPlots <- function(specificTrait, inputData) {
  testLabels <- c("bustedPH" = "Distribution of\npositive selection\n(BUSTED-PH)",
                  "relax" = "Distribution of\nselection intensity\n(RELAX)")
  plotData <- dplyr::filter(inputData, 
                            trait == specificTrait)
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
    geom_text(data = dplyr::filter(inputData, 
                                   trait == specificTrait), 
              aes(x = 1.5,
                  y = as.numeric(maxHeight) + 400,
                  label = pValue),
              colour = "grey50") + 
    geom_linerange(data = dplyr::filter(inputData, 
                                        trait == specificTrait),
                   aes(xmin = 1, 
                       xmax = 2, 
                       y = (as.numeric(maxHeight) + 200)), 
                   color = "grey50", 
                   size = 0.3) +
    theme_bw() +
    scale_color_manual(values = colorScheme,
                       labels = c("nSelectionIntensified" = "Selection intensified\non foreground",
                                  "nSelectionRelaxed" = "Selection relaxed\non foreground",
                                  "nSelectionForeground" = "Positive selection\non foreground",
                                  "nSelectionBackground" = "Positive selection\non background")) +
    labs(y = "Count of orthogroups", 
         title = paste("Patterns of molecular evolution associated with", 
                       plotData$traitLabels, 
                       sep = " ")) +
    scale_x_discrete(labels = c("bustedPH" = "Positive selection\non fore- or background",
                                "relax" = "Relaxed or\nintensified selection")) + 
    facet_wrap(~test,
               labeller = labeller(test = testLabels),
               scales = 'free_x', 
               strip.position = "bottom") + 
    theme(axis.text.x = element_text(angle = 0,
                                     hjust = 0.5,
                                     vjust = 1),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.spacing = unit(0.5, "cm"),
          panel.border = element_rect(fill = "NA", 
                                      color = "gray90", 
                                      size = 0.5, 
                                      linetype = "solid"),
          panel.background = element_rect(fill = "gray95"),
          legend.position = "none",
          panel.spacing.x = unit(0, "pt"), 
          panel.spacing.y = unit(0, "pt"),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(), 
          rect = element_rect(fill = "transparent")) +
    scale_x_discrete(labels = c("nSelectionBackground" = "Positive selection\non species\nwithout trait",
                                "nSelectionForeground" = "Positive selection\non species\nwith trait",
                                "nSelectionIntensified" = "Intensified selection\nassociated\nwith trait",
                                "nSelectionRelaxed"  = "Relaxed selection\nassociated\nwith trait"))
}

test <- purrr::map(traits, ~lollipopPlots(.x, inputData = summaryDataframe))
cowplot::plot_grid(plotlist = test,
                   nrow = 1)
patchworkTest <- patchwork::wrap_plots(test, 
                                       nrow = 1)
patchworkTest[[1]] <- patchworkTest[[1]] +
  ylim(0, 2700)
patchworkTest[[2]] <- patchworkTest[[2]] + 
  ylim(0, 2700)

patchworkTest[[2]] <- patchworkTest[[2]] + 
  theme(axis.text.y = element_blank(),           
        axis.ticks.y = element_blank(),        
        axis.title.y = element_blank()) 
  
patchworkTest[[4]] = patchworkTest[[4]] + theme(axis.text.y = element_blank(),
                                                axis.ticks.y = element_blank(),
                                                axis.title.y = element_blank())


patchworkTest
ggsave(filename = "posterComparativeGenomicsResultsLollipop1.svg", 
       device = "svg",
       path = "./Plots/", 
       bg = "transparent", 
       width = 12, 
       height = 4.5)

for (i in test) {
  trait <- unique(i[["data"]][["trait"]])
  filename <- paste(trait, "ComparativeGenomicsResultsLollipop.png", sep = "")
  plot(i)
  print(filename)
  ggsave(filename = filename, path = "./Plots/", bg = "transparent", width = 7, height = 5)
}

# Construct euler diagrams (like Venn diagrams but doesn't have to show all possible intersections)
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
  eulerr_options(padding = grid::unit(0.1, "lines"))
  test <- plot(intersectionsOfSelection,
          quantities = list(fontsize = 15,
                            labels = paste(intersectionsOfSelection$original, "\northogroups", sep = "")),
          labels = list(fontsize = 0),
          adjust_labels = TRUE,
          #main = paste("Combinations of selective regimes related to", tolower(gsub("([A-Z])", " \\1", specificTrait)), sep = " "),
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

# Do these same analyses but only for the single-copy orthogroups:
# Get a list of the single copy orthogroups:
singleCopyOrthogroups <- list.files(path = paste("./5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Single_Copy_Orthologue_Sequences/", sep = ""), full.names = FALSE)
singleCopyOrthogroups <- sapply(strsplit(as.character(singleCopyOrthogroups), "[.]"), head, 1)

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

# Make lollipop plots:
test <- purrr::map(traits, ~lollipopPlots(.x, inputData = summaryDataframeSingle))
cowplot::plot_grid(plotlist = test,
                   nrow = 1)








