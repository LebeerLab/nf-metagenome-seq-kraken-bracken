// PARAMETERS
params.CONTAINER = "tehoaphidian/kraken-bracken"
params.OUTPUT = "results"

process CREATE_TIDYAMPLICONS {
    publishDir "${params.outdir}",  mode:  'copy'
    container params.CONTAINER

    input:
    path("*")

    output:
    tuple path("tidyamplicons/samples.csv"), path("tidyamplicons/taxa.csv"), path("tidyamplicons/abundances.csv")


    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(tidyamplicons)

    source('$baseDir/bin/functions.R')

    run <- Sys.Date()
    pipeline <- "${launchDir.getName()}"

    # convert kraken2 results to a nice taxonomy table
    kraken2taxtable(".", "taxtable")

    ta <- import_tidyamplicons("taxtable")

    # add the run and pipeline names to the tidyamplicons object
    ta\$samples\$run <- run
    ta\$samples\$pipeline <- pipeline

    # save the tidyamplicons object as three tidy tables
    ta %>% write_tidyamplicons("tidyamplicons")

    """
}

process PRINT_TOP10 {
    container params.CONTAINER
    input:
    path(ta)
    
    output:
    stdout

    script:
    """
    top10.R ${ta}
    """
}