#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(tidyamplicons))

args = commandArgs(trailingOnly=TRUE)

ta <- read_tidyamplicons(args[1])
t_abundances <- ta %>% 
    filter_samples(sample == "ERR5621728") %>%
    add_rel_abundance() %>%
    abundances() %>%
    arrange(desc(rel_abundance)) %>%
    head(10)

t_abundances %>% left_join(ta$taxa, by="taxon_id") %>% select(rel_abundance, genus, species)
