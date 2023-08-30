#!/usr/bin/env Rscript
# Author: Stijn Wittouck
# Last modified: 06/01/2022

library(tidyverse)
library(tidyamplicons)

args = commandArgs(trailingOnly=TRUE)
run <- Sys.Date()
pipeline = args[1]

# converts per-sample kraken2 results to a taxonomy table
kraken2taxtable <- function(din_krakensamples, fout_taxtable, file_pattern=".mpa$") {
  readcount_M <- list.files(din_krakensamples, pattern = file_pattern) %>%
    {names(.) <- str_extract(., "^[^_]+"); .} %>%
    map(~ read_tsv(
      skip=1,
      str_c(din_krakensamples, ., sep = "/"), 
      col_names = c("classification", "readcount", "taxid"),
      col_types = cols(classification = col_character(), readcount = col_double(), taxid = col_character()),
      lazy = F
    )) %>%
    # remove empty tables (they'll give trouble otherwise)
    keep(~ nrow(.) != 0) %>%
    map(calculate_unclassified) %>%
    #map2(., names(.), ~ {names(.x) <- c("classification", .y); .x}) %>%
    map2(., names(.), ~ {names(.x) <- c("classification", "taxid", "readcount"); .x}) %>%
    reduce(., full_join, by = c("classification", "taxid")) %>%
    replace(., is.na(.), 0) 
    
    readcount_M %>%
    select(-taxid) %>%
    write_tsv(fout_taxtable)

    readcount_M %>%
    select(-classification) %>%
    write_tsv(paste0("taxid_", fout_taxtable))

}

# needed by kraken2taxtable
calculate_unclassified <- function(taxa) {
  
  taxa <- 
    taxa %>%
    # remove the kingdom level 
    # (because: not present for bacteria, in krakendb)
    filter(classification != "d__Eukaryota|k__Metazoa") %>%
    mutate_at("classification", str_remove, fixed("|k__Metazoa")) %>% 
    rename(total = readcount)
  
  taxa_classified <-
    taxa %>%
    filter(str_detect(classification, "\\|")) %>%
    mutate_at("classification", str_remove, "\\|[^|]+$") %>%
    group_by(classification) %>%
    summarize(classified = sum(total))
  
  taxa %>%
    left_join(taxa_classified, by = "classification") %>%
    mutate(unclassified = total - classified) %>%
    mutate(readcount = if_else(is.na(unclassified), total, unclassified)) %>%
    select(classification, taxid, readcount)
  
}

# imports a taxonomy table as a tidyamplicons object
import_tidyamplicons <- function(fin_taxtable, domains_present = T) {
  
  if (domains_present) {
    
    taxlevels <- 
      c("domain", "phylum", "class", "order", "family", "genus", "species")
    
  } else {
    
    taxlevels <- 
      c("phylum", "class", "order", "family", "genus", "species")
    
  }
  
  fin_taxtable %>%
    read.table(
      header = T, sep = "\t", row.names = 1, stringsAsFactors = F, quote = ""
    ) %>%
    as.matrix() %>%
    create_tidyamplicons(taxa_are_columns = F) %>%
    modify_at(
      "taxa", separate, taxon, into = taxlevels, sep = "\\|", 
      fill = "right"
    ) %>%
    set_rank_names(taxlevels)
  
}


# convert kraken2 results to a nice taxonomy table
kraken2taxtable(".", "taxtable")

ta <- import_tidyamplicons("taxtable")
ta_taxid <- import_tidyamplicons("taxid_taxtable") %>% 
  taxa() %>% 
  mutate(taxid = domain) %>%
  select(taxon_id, taxid)

ta$taxa <- ta$taxa %>% left_join(ta_taxid, by=c("taxon_id"="taxon_id"))

# add the run and pipeline names to the tidyamplicons object
ta$samples$run <- run
ta$samples$pipeline <- pipeline

if ("domain" %in% colnames(ta$taxa)) {
  ta$taxa <- ta$taxa %>% 
    mutate(kingdom = domain) %>% 
    select(-domain) %>%
    select("kingdom", everything())
}

# save the tidyamplicons object as three tidy tables
ta %>% write_tidyamplicons("tidyamplicons")
