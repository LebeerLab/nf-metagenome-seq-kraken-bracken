library(tidyverse)
library(tidytacos)

RN <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")

M <- read.csv("count_matrix.tsv", sep="\t", row.names=1)
tt <- create_tidytacos(as.matrix(M), taxa_are_columns=F)
tt$taxa <- tt$taxa %>% separate(taxon,RN, sep="\\|") 
tt %>% 
  set_rank_names(RN) %>%
  write_tidytacos("tidytacos")
