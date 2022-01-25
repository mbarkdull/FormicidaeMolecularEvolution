library(gt)
library(tidyverse)
library(googlesheets4)

# Read in phenotype data:
rawFocalPhenotypes <- read_sheet("https://docs.google.com/spreadsheets/d/1O0aVh6LkTpS4MvJjDegTTy1IrIueJrY65s9xwb8DML8/edit?usp=sharing")
# Prune to columns and rows you want:
dataSourceTable <- rawFocalPhenotypes %>%
  filter(Using == "yes")%>%
  select(Species, GenomeCitation, TranscriptDownloadLink, ProteinsLink, GFFLink) %>%
  gt() %>% 
  cols_label(Species = md("Species"),
             TranscriptDownloadLink = md("Source for coding sequences"),
             ProteinsLink = md("Source for proteome sequences"),
             GFFLink = md("Source for GFF annotations"),
             GenomeCitation = md("Genome citation") 
  ) %>% 
  tab_style(style = list(cell_text(style = "italic")),
            locations = cells_body(columns = Species)) %>%
  opt_row_striping() %>% 
  tab_header(title = "Sources for genomic data used in this paper") %>%
  tab_options(container.width = c(10000),
              table.width = c(1000))

gtsave(dataSourceTable, filename = "dataSourceTable.png")

phenotypeDataTable <- rawFocalPhenotypes %>%
  filter(Using == "yes")%>%
  select(Species,
         Subfamily,
         Polyandry, 
         Polygyny, 
         WorkerReproductionOffspringType, 
         DegreeOfWorkerPolymorphism) %>%
  arrange(Subfamily) %>%
  gt() %>% 
  cols_label(Species = md("Species"),
             Subfamily = md("Subfamily"),
             Polyandry = md("Polyandry"),
             Polygyny = md("Polygyny"),
             WorkerReproductionOffspringType = md("Worker reproductive capacity"),
             DegreeOfWorkerPolymorphism = md("Degree of worker polymorphism") 
  ) %>% 
  tab_style(style = list(cell_text(style = "italic")),
            locations = cells_body(columns = Species)) %>%
  opt_row_striping() %>% 
  tab_header(title = "Trait data") %>%
  tab_options(container.width = c(10000),
              table.width = c(1000))
phenotypeDataTable
gtsave(phenotypeDataTable, filename = "phenotypeDataTable.png")



phenotypeSourcesDataTable <- rawFocalPhenotypes %>%
  filter(Using == "yes")%>%
  select(Species,
         Subfamily,
         Polyandry, 
         PolyandryCitation,
         Polygyny, 
         PolygynyCitation,
         WorkerReproductionOffspringType, 
         WorkerReproductionOffspringTypeCitation,
         DegreeOfWorkerPolymorphism,
         DegreeOfWorkerPolymorphismCitation) %>%
  arrange(Subfamily) %>%
  gt() %>% 
  cols_label(Species = md("Species"),
             Subfamily = md("Subfamily"),
             Polyandry = md("Polyandry"),
             PolyandryCitation = md("Citation for polyandry"),
             Polygyny = md("Polygyny"),
             PolygynyCitation = md("Citation for polygyny"),
             WorkerReproductionOffspringType = md("Worker reproductive capacity"),
             WorkerReproductionOffspringTypeCitation = md("Citation for worker reproduction"),
             DegreeOfWorkerPolymorphism = md("Degree of worker polymorphism"),
             DegreeOfWorkerPolymorphismCitation = md("Citation for worker polymorphism")
  ) %>% 
  tab_style(style = list(cell_text(style = "italic")),
            locations = cells_body(columns = Species)) %>%
  opt_row_striping() %>% 
  tab_header(title = "Trait data and sources") %>%
  tab_options(container.width = c(10000),
              table.width = c(1000))
phenotypeSourcesDataTable
gtsave(phenotypeSourcesDataTable, filename = "phenotypeSourcesDataTable.png")
