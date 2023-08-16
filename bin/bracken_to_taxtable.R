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
  list.files(din_krakensamples, pattern = file_pattern) %>%
    {names(.) <- str_extract(., "^[^_]+"); .} %>%
    map(~ read_tsv(
      skip=1,
      str_c(din_krakensamples, ., sep = "/"), 
      col_names = c("classification", "readcount"),
      col_types = cols(classification = col_character(), readcount = col_double()),
      lazy = F
    )) %>%
    # remove empty tables (they'll give trouble otherwise)
    keep(~ nrow(.) != 0) %>%
    map(calculate_unclassified) %>%
    map2(., names(.), ~ {names(.x) <- c("classification", .y); .x}) %>%
    reduce(., full_join, by = "classification") %>%
    replace(., is.na(.), 0) %>%
    write_tsv(fout_taxtable)
    
}

# util to reduce mpa output to taxon-coverage
select_coverage <- function(din){

	for (file in list.files(din, pattern=".mpa$")){
		f <- read_tsv(file, skip=1, 
			col_names = c("taxon","readcount","genomesize","coverage"),
			col_types = ("cddd")
			)
 
		f %>% select(taxon, coverage) %>% rename(readcount=coverage) %>%
			write_tsv(paste0(file ,"_cov"))
	} 
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
    select(classification, readcount)
  
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

# add the run and pipeline names to the tidyamplicons object
ta$samples$run <- run
ta$samples$pipeline <- pipeline

# add genomesizes and norm_readcount to abundances table
## build new ta object with coverage

ta %>% write_tidyamplicons("ta_interm")

#ta$abundances <- ta$abundances %>% 
#left_join(ta_cov$abundances, by=c("taxon_id", "sample_id"))

# save the tidyamplicons object as three tidy tables
ta %>% write_tidyamplicons("tidyamplicons")
