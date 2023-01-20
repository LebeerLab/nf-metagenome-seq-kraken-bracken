#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(tidyamplicons))

ta <- read_tidyamplicons('results/tidyamplicons')
t_abundances <- ta %>% 
    filter_samples(sample == "SRR20285099") %>%
    add_rel_abundance() %>%
    abundances() %>%
    arrange(desc(rel_abundance)) %>%
    head(10)

t_abundances %>% left_join(ta$taxa, by="taxon_id") %>% select(rel_abundance, genus, species)
