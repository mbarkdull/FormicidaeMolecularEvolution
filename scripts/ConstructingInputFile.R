library(googlesheets4)
inputdata <- read_sheet("https://docs.google.com/spreadsheets/d/1O0aVh6LkTpS4MvJjDegTTy1IrIueJrY65s9xwb8DML8/edit?usp=sharing")
inputdata <- filter(inputdata, !is.na(inputdata$SpeciesAbbreviation)) %>%
  select(TranscriptDownloadLink, ProteinsLink, GFFLink, SpeciesAbbreviation, FilterIsoforms)
write.table(inputdata, file = "inputurls_full.txt", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
