# Load required packages:
library(tidyverse)
library(ggthemes)
library(gt)
library(plyr)
library(rjson)
library(MetBrewer)
library(SuperExactTest)
library(topGO)
library(ggpubr)

##### Analyze results for all orthogroups: ##### 
# Read in results:
relaxResults <- read_csv("./Results/allRelaxResults.csv", col_names = TRUE)
relaxResults$pValue <- as.numeric(as.character(relaxResults$pValue))
relaxResults$kValue <- as.numeric(as.character(relaxResults$kValue))

bustedPHResults <- read_csv("./Results/allBustedPHResults.csv", col_names = TRUE)
bustedPHResults$`test results background p-value` <- as.numeric(as.character(bustedPHResults$`test results background p-value`))
bustedPHResults$`test results p-value` <- as.numeric(as.character(bustedPHResults$`test results p-value`))
bustedPHResults$`test results shared distributions p-value` <- as.numeric(as.character(bustedPHResults$`test results shared distributions p-value`))

# Adjust the raw p-values with Benjamini-Hochberg FDR:
relaxResultsAdjusted <- relaxResults %>%
  group_by(trait) %>% 
  mutate(pValueFDR = p.adjust(pValue, method='BH')) 
bustedPHResultsAdjusted <- bustedPHResults %>%
  group_by(trait) %>% 
  mutate(testResultspValueFDR = p.adjust(`test results p-value`, method='BH')) %>% 
  mutate(testResultsBackgroundpValueFDR = p.adjust(`test results background p-value`, method='BH')) %>% 
  mutate(testResultsSharedDistributionspValueFDR = p.adjust(`test results shared distributions p-value`, method='BH'))

# Combine relax and busted results:
relaxAndBustedPH <- full_join(relaxResultsAdjusted,
                              bustedPHResultsAdjusted,
                              by = c("orthogroup" = "orthogroup",
                                     "trait" = "trait"))

relaxAndBustedPH$combo <- paste(relaxAndBustedPH$selectionOn, relaxAndBustedPH$selectionCategory)
write.csv(relaxAndBustedPH, file = "relaxAndBustedPH.csv")

traits <- unique(relaxAndBustedPH$trait)

# Calculate p-values for each trait, on the FDR adjusted p-values:
pValueByTraitRelax <- function(specificTrait, inputData) {
  relaxResultsTrait <- filter(inputData, 
                              trait == specificTrait & !is.na(kValue))
  nSelectionIntensified <- relaxResultsTrait %>% 
    filter(pValueFDR <= 0.05,
           kValue > 1) %>%
    nrow()
  
  nSelectionRelaxed <- relaxResultsTrait %>% 
    filter(pValueFDR <= 0.05,
           kValue < 1) %>%
    nrow()
  nTotal <- relaxResultsTrait %>% 
    nrow()
  # Construct a frequency table with that information:
  shiftInIntensity <- c(nSelectionIntensified, 
                        nSelectionRelaxed)
  population <- c((nrow(relaxResultsTrait) - nSelectionIntensified),
                  (nrow(relaxResultsTrait) - nSelectionRelaxed))
  proportionTest <- prop.test(shiftInIntensity,
                              population)
  # Print the resulting p-value:
  proportionTest[["p.value"]]
  pValue <- if (proportionTest[["p.value"]] < 0.000001) {
    formatC(proportionTest[["p.value"]], format = "e", digits = 2)
  } else {
    round(proportionTest[["p.value"]], digits = 4)
  }
  textHeight <- as.numeric(max(nSelectionIntensified, nSelectionRelaxed))
  return(c(specificTrait, pValue, textHeight, nSelectionIntensified, nSelectionRelaxed, nTotal))
}

possiblypValueByTraitRelax <- possibly(pValueByTraitRelax, otherwise = "Error")

pValuesRelax <- traits %>% 
  purrr::map(~ possiblypValueByTraitRelax(.x, relaxAndBustedPH))

pValuesRelax <- as.data.frame(do.call(rbind, pValuesRelax))   

colnames(pValuesRelax) <- c("trait", "pValue", "maxHeight", "nSelectionIntensified", "nSelectionRelaxed", "nTotal")

pValuesRelax$test <- "relax"

pValueByTraitBUSTEDPH <- function(specificTrait, inputData) {
  bustedPHResultsTrait <- filter(inputData, 
                                 trait == specificTrait & !is.na(`test results p-value`))
  nSelectionForeground <- bustedPHResultsTrait %>% 
    filter(testResultspValueFDR <= 0.05,
           testResultsBackgroundpValueFDR > 0.05,
           testResultsSharedDistributionspValueFDR <= 0.05) %>%
    nrow()
  
  nSelectionBackground <- bustedPHResultsTrait %>% 
    filter(testResultspValueFDR > 0.05,
           testResultsBackgroundpValueFDR <= 0.05,
           testResultsSharedDistributionspValueFDR <= 0.05) %>%
    nrow()
  nTotal <- bustedPHResultsTrait %>% 
    nrow()
  
  # Construct a frequency table with that information:
  selected <- c(nSelectionForeground, nSelectionBackground)
  population <- c((nrow(bustedPHResultsTrait) - nSelectionForeground),
                  (nrow(bustedPHResultsTrait) - nSelectionBackground))
  # Run a test of equal proportions:
  proportionTest <- prop.test(selected,
                              population)
  # Print the resulting p-value:
  proportionTest[["p.value"]]
  pValue <- if (proportionTest[["p.value"]] < 0.000001) {
    formatC(proportionTest[["p.value"]], format = "e", digits = 2)
  } else {
    round(proportionTest[["p.value"]], digits = 4)
  }
  textHeight <- as.numeric(max(nSelectionForeground, nSelectionBackground))
  return(c(specificTrait, pValue, textHeight, nSelectionForeground, nSelectionBackground, nTotal))
}
possiblypValueByTraitBUSTEDPH <- possibly(pValueByTraitBUSTEDPH, otherwise = "Error")
pValuesBUSTEDPH <- traits %>% 
  map(~ possiblypValueByTraitBUSTEDPH(.x, relaxAndBustedPH))
pValuesBUSTEDPH <- as.data.frame(do.call(rbind, pValuesBUSTEDPH))   
colnames(pValuesBUSTEDPH) <- c("trait", "pValue", "maxHeight", "nSelectionForeground", "nSelectionBackground", "nTotal")
pValuesBUSTEDPH$test <- "bustedPH"

pValuesAll <- bind_rows(pValuesRelax, pValuesBUSTEDPH)

# Make a lollipop plot:
# Create a summary of the data
summaryData <- function(specificTrait) {
  relaxResultsTrait <- filter(relaxResultsAdjusted, 
                              trait == specificTrait)
  nSelectionIntensified <- relaxResultsTrait %>% 
    filter(pValueFDR <= 0.05,
           kValue > 1) %>%
    nrow()
  
  nSelectionRelaxed <- relaxResultsTrait %>% 
    filter(pValueFDR <= 0.05,
           kValue < 1) %>%
    nrow()
  bustedPHResultsTrait <- filter(bustedPHResultsAdjusted, 
                                 trait == specificTrait)
  nSelectionForeground <- bustedPHResultsTrait %>% 
    filter(testResultspValueFDR <= 0.05,
           testResultsBackgroundpValueFDR > 0.05,
           testResultsSharedDistributionspValueFDR <= 0.05) %>%
    nrow()
  
  nSelectionBackground <- bustedPHResultsTrait %>% 
    filter(testResultspValueFDR > 0.05,
           testResultsBackgroundpValueFDR <= 0.05,
           testResultsSharedDistributionspValueFDR <= 0.05) %>%
    nrow()
  
  nSelectionBoth <- bustedPHResultsTrait %>% 
    filter(testResultspValueFDR <= 0.05,
           testResultsBackgroundpValueFDR <= 0.05) %>%
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
colorScheme <- c("nSelectionIntensified" = "#48266e",
                 "nSelectionRelaxed" = "#819fe3",
                 "nSelectionForeground" = "#bf3934",
                 "nSelectionBackground" = "#e3a59f",
                 "nSelectionBoth" = "#e0d8d7")
traitLabels <- c("polyandry" = "polyandry",
                 "polygyny" = "polygyny",
                 "workerReproductionQueens" = "reproductive workers",
                 "workerPolymorphism" = "polymorphic workers",
                 "multilineage" = "multilineage colonies")
traitLabels <- data.frame(traitLabels, trait = names(traitLabels))

# Add the trait labels to the summary data:
summaryDataframe <- left_join(summaryDataframe, traitLabels)

traits <- unique(summaryDataframe$trait)
summaryDataframe$selection <- factor(summaryDataframe$selection,
                                     levels = c("nSelectionForeground",
                                                "nSelectionBackground",
                                                "nSelectionBoth",
                                                "nSelectionRelaxed",
                                                "nSelectionIntensified"))
summaryDataframe$maxHeight <- as.numeric(as.character(summaryDataframe$maxHeight))
summaryDataframe$nSelectionIntensified <- as.numeric(as.character(summaryDataframe$nSelectionIntensified))
summaryDataframe$nSelectionRelaxed <- as.numeric(as.character(summaryDataframe$nSelectionRelaxed))
summaryDataframe$nSelectionForeground <- as.numeric(as.character(summaryDataframe$nSelectionForeground))
summaryDataframe$nSelectionBackground <- as.numeric(as.character(summaryDataframe$nSelectionBackground))

#### Plot only BUSTED-PH results with bar plots: ####
lollipopPlots <- function(specificTrait, inputData) {
  testLabels <- c("bustedPH" = "Distribution of positive selection\n(BUSTED-PH)",
                  "relax" = "Distribution of selection intensity\n(RELAX)")
  plotData <- dplyr::filter(inputData, 
                            trait == specificTrait,
                            selection != "nSelectionBoth",
                            test == "bustedPH")
  plot <- ggplot(plotData, 
                 aes(x = selection,
                     y = number,
                     fill = selection)) +
    geom_col(plotData,
             mapping = aes(fill = selection),
             position = position_dodge2(width = 0.75)) +
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
    scale_fill_manual(values = colorScheme,
                      labels = c("nSelectionIntensified" = "Selection intensified\non foreground",
                                 "nSelectionRelaxed" = "Selection relaxed\non foreground",
                                 "nSelectionForeground" = "Positive selection\non foreground",
                                 "nSelectionBackground" = "Positive selection\non background")) +
    labs(y = "Count of orthogroups", 
         title = paste(str_to_sentence(plotData$traitLabels), 
                       sep = " ")) +
    scale_x_discrete(labels = c("bustedPH" = "Positive selection\non fore- or background",
                                "relax" = "Relaxed or\nintensified selection")) + 
    theme(axis.text.x = element_text(angle = 0,
                                     hjust = 0.5,
                                     vjust = 1,
                                     size = 10),
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
  plot
}
subsetTraits <- c("workerReproductionQueens", "workerPolymorphism")
test <- purrr::map(subsetTraits, ~lollipopPlots(.x, inputData = summaryDataframe))
cowplot::plot_grid(plotlist = test,
                   ncol = 2)

# Render the plots together in a patchwork:
patchworkTest <- patchwork::wrap_plots(test, 
                                       ncol = 2) + 
  patchwork::plot_annotation(title = 'Impact of sociality-associated traits on molecular evolution',
                             theme = theme(plot.title = element_text(size = 19)))
# Set the individual plots to have the same y limits:
patchworkTest[[1]] <- patchworkTest[[1]] +
  ylim(0, 2100)
patchworkTest[[2]] <- patchworkTest[[2]] + 
  ylim(0, 2100)

patchworkTest[[1]] <- patchworkTest[[1]] +
  theme(strip.text.x = element_text(color = "white")) 

patchworkTest[[2]] <- patchworkTest[[2]] + 
  theme(axis.text.y = element_blank(),           
        axis.ticks.y = element_blank(),        
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 11)) 

patchworkTest

ggsave(filename = "bustedPHResultsBarPlot.png", 
       device = "png",
       path = "./Plots/", 
       bg = "transparent", 
       width = 8, 
       height = 4)

##### RELAX results ####
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
relaxResultsAdjusted <- left_join(relaxResultsAdjusted, traitLabels)
rm(traitLabels)

kValuePlot <- ggplot(data = filter(relaxResultsAdjusted,
                                   !is.na(kValue),
                                   trait %in% subsetTraits,
                                   pValueFDR <= 0.05)) +
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

kValueViolinPlot <- ggplot(data = filter(relaxResultsAdjusted,
                                         !is.na(kValue),
                                         trait %in% subsetTraits,
                                         pValueFDR <= 0.05)) +
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
kValueViolinPlotFlipped <- ggplot(data = filter(relaxResultsAdjusted,
                                                !is.na(kValue),
                                                trait %in% subsetTraits,
                                                pValueFDR <= 0.05)) +
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
  inputData <- summaryDataframe
  testLabels <- c("bustedPH" = "Distribution of positive selection\n(BUSTED-PH)",
                  "relax" = "Distribution of selection intensity")
  plotData <- dplyr::filter(inputData, 
                            trait == specificTrait,
                            selection != "nSelectionBoth")
  barPlot <- ggplot(filter(plotData,
                           test != "bustedPH"), 
                    aes(x = selection,
                        y = number,
                        fill = selection)) +
    geom_col(filter(plotData,
                    test != "bustedPH"),
             mapping = aes(fill = selection),
             position = position_dodge2(width = 0.75)) +
    geom_text(data = dplyr::filter(plotData, 
                                   trait == specificTrait,
                                   test != "bustedPH"), 
              aes(x = 1.5,
                  y = as.numeric(as.character(maxHeight)) + 300,
                  label = pValue),
              colour = "grey50",
              size = 3) + 
    geom_linerange(data = dplyr::filter(plotData, 
                                        trait == specificTrait,
                                        test != "bustedPH"),
                   aes(xmin = 1, 
                       xmax = 2, 
                       y = (as.numeric(as.character(maxHeight)) + 200)), 
                   color = "grey50", 
                   size = 0.3) +
    scale_fill_manual(values = colorScheme,
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
          strip.text = element_text(hjust = 0),
          axis.title.x = element_blank(), 
          rect = element_rect(fill = "transparent"),
          plot.title = element_text(size = 15)) +
    scale_x_discrete(labels = c("nSelectionBackground" = "Positive selection\non species without trait",
                                "nSelectionForeground" = "Positive selection\non species with trait",
                                "nSelectionIntensified" = "Intensified selection\nassociated with trait",
                                "nSelectionRelaxed"  = "Relaxed selection\nassociated with trait")) +
    ylim(0, 3200)
  
  violin <- ggplot(data = filter(relaxResultsAdjusted,
                                 !is.na(kValue),
                                 trait == specificTrait,
                                 pValueFDR <= 0.05)) +
    geom_violin(mapping = aes(x = traitLabels, 
                              y = as.numeric(kValue)),
                fill = "#efebf0") +
    geom_hline(yintercept = 1) + 
    labs(title = paste("Distribution of selection intensity parameters\nassociated with each focal trait"),
         y = "K-values") +
    scale_y_continuous(trans = "log10",
                       breaks = scales::log_breaks(n = 10),
                       labels = function(x) sprintf("%g", x),
                       limits = c(0.001, 400)) + 
    coord_flip() +
    theme(axis.text.y = element_text(color = "white",
                                     angle = 90),
          axis.ticks.y = element_line(color = "white"),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1),
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
          plot.title = element_blank(),
          strip.text = element_blank()) +
    facet_wrap(~factor(traitLabels, 
                       levels = c('Reproductive workers',
                                  'Polymorphic workers',
                                  'Multilineage colonies')),
               ncol = 3,
               scales = 'free_y', 
               strip.position = "top") 
  
  ggarrange(barPlot,
            violin,
            ncol = 1,
            heights = c(2, 1))
}

lollipopViolin(specificTrait = "workerPolymorphism")
test <- purrr::map(subsetTraits, lollipopViolin)
patchworkTest <- cowplot::plot_grid(plotlist = test,
                                    ncol = 2)

patchworkTest
ggsave(filename = "relaxResultsBarAndViolin.png", 
       device = "png",
       path = "./Plots/", 
       bg = "transparent", 
       width = 8, 
       height = 5)

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

#### Look at only single-copy orthogroups: ####
singleCopyOrthogroups <- read_delim("5_OrthoFinder/fasta/OrthoFinder/Results_Jul13/Orthogroups/Orthogroups_SingleCopyOrthologues.txt",
                                    delim = "\t",
                                    col_names = FALSE)

relaxResultsAdjustedUngrouped <- ungroup(relaxResultsAdjusted)

singleCopyOrthogroupResults <- dplyr::filter(relaxResultsAdjustedUngrouped, 
                                             relaxResultsAdjustedUngrouped$orthogroup %in% singleCopyOrthogroups$X1) %>%
  group_by(trait)

pValuesRelaxSingleCopy <- traits %>% 
  purrr::map(~ possiblypValueByTraitRelax(.x, singleCopyOrthogroupResults))
pValuesRelaxSingleCopy <- as.data.frame(do.call(rbind, pValuesRelaxSingleCopy))   
colnames(pValuesRelaxSingleCopy) <- c("trait", "pValue", "maxHeight", "nSelectionIntensified", "nSelectionRelaxed", "nTotal")
pValuesRelaxSingleCopy$test <- "relax"

bustedPHResultsAdjustedUngrouped <- ungroup(bustedPHResultsAdjusted)

singleCopyOrthogroupResults <- dplyr::filter(bustedPHResultsAdjustedUngrouped, 
                                             bustedPHResultsAdjustedUngrouped$orthogroup %in% singleCopyOrthogroups$X1) %>%
  group_by(trait)

pValuesBUSTEDPHSingleCopy <- traits %>% 
  map(~ possiblypValueByTraitBUSTEDPH(.x, singleCopyOrthogroupResults))
pValuesBUSTEDPHSingleCopy <- as.data.frame(do.call(rbind, pValuesBUSTEDPHSingleCopy))   
colnames(pValuesBUSTEDPHSingleCopy) <- c("trait", "pValue", "maxHeight", "nSelectionForeground", "nSelectionBackground", "nTotal")
pValuesBUSTEDPHSingleCopy$test <- "bustedPH"

pValuesAllSingleCopy <- bind_rows(pValuesRelaxSingleCopy, pValuesBUSTEDPHSingleCopy)

#### Plot the single-copy RELAX results ####


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
