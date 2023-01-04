# Author: Stijn Wittouck
# Last modified: 02/11/2021

# runs kraken2 on a set of forward and reverse sample pairs
run_kraken <- function(din_samples, din_database, dout) {
  
  if (! dir.exists(dout)) dir.create(dout)
  
  samples <- 
    tibble(
      file_f = list.files(din_samples, pattern = "_R1_001.fastq") %>% sort(),
      file_r = list.files(din_samples, pattern = "_R2_001.fastq") %>% sort()
    ) %>%
    separate(
      file_f, c("description", "sample", "lane", "direction", "appendix"), 
      sep = "_", remove = F
    ) %>%
    select(sample, description, file_f, file_r) %>%
    filter(description != "Undetermined") %>%
    mutate(file_out = str_c(description, "_counts.tsv", sep = ""))
  
  samples %>%
    mutate(
      file_f = str_c(din_samples, file_f, sep = "/"),
      file_r = str_c(din_samples, file_r, sep = "/"),
      file_out = str_c(dout, file_out, sep = "/")
    ) %>%
    select(file_f, file_r, file_out) %>%
    pwalk(function(file_f, file_r, file_out) {
      run_kraken_onesample(
        freadsfin = file_f, rreadsfin = file_r, db = din_database,
        fout = file_out
      )
    })
  
}

# needed by run_kraken
run_kraken_onesample <- function(freadsfin, rreadsfin, db, fout) {
  
  if (file.exists(fout)) {
    message("kraken output already exists for file")
  } else {
    system(glue::glue(
      "kraken2 --db {db} --report {fout} --use-mpa-style --threads 16 --paired ",
      "{freadsfin} {rreadsfin} > /dev/null"
    ))
  }
  
}

# converts per-sample kraken2 results to a taxonomy table
kraken2taxtable <- function(din_krakensamples, fout_taxtable) {
  
  list.files(din_krakensamples, pattern = ".mpa") %>%
    {names(.) <- str_extract(., "^[^_]+"); .} %>%
    map(~ read_tsv(
      str_c(din_krakensamples, ., sep = "/"), 
      col_names = c("classification", "count"),
      col_types = cols(classification = col_character(), count = col_integer()),
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

# needed by kraken2taxtable
calculate_unclassified <- function(taxa) {
  
  taxa <- 
    taxa %>%
    # remove the kingdom level 
    # (because: not present for bacteria, in krakendb)
    filter(classification != "d__Eukaryota|k__Metazoa") %>%
    mutate_at("classification", str_remove, fixed("|k__Metazoa")) %>% 
    rename(total = count)
  
  taxa_classified <-
    taxa %>%
    filter(str_detect(classification, "\\|")) %>%
    mutate_at("classification", str_remove, "\\|[^|]+$") %>%
    group_by(classification) %>%
    summarize(classified = sum(total))
  
  taxa %>%
    left_join(taxa_classified, by = "classification") %>%
    mutate(unclassified = total - classified) %>%
    mutate(count = if_else(is.na(unclassified), total, unclassified)) %>%
    select(classification, count)
  
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
