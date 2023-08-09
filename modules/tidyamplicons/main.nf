// PARAMETERS
params.CONTAINER = "theoaphidian/kraken-bracken"
params.OUTPUT = "results"
params.SKIP_NORM = true
params.GENOMESIZES = null
params.TEST_PIPELINE = false
params.REPORT_TYPE = "bracken"

process CONVERT_MPA {
    tag "${pair_id}"

    input:
    tuple val(pair_id), path(brck_rpt), val(readlen) 

    output:
    tuple val(pair_id), path("${pair_id}_bracken.report.mpa"), val(readlen) 

    script:
    """
    kreport2mpa.py -r ${brck_rpt} -o "${pair_id}_bracken.report.mpa"
    """
}

process NORMALIZE_READCOUNT {
    tag "${pair_id}"

    input:
    tuple val(pair_id), path(mpa_rpt), val(readlen) 
    path(genomesizes)

    output:
    path("${pair_id}_normalized_rc.mpa")
    
    script:
    """
    normalize_abundances.py ${mpa_rpt} "${genomesizes}" "${pair_id}_normalized_rc.mpa" --factor ${readlen}
    """
}

process CREATE_TIDYAMPLICONS {
    publishDir "${params.OUTPUT}",  mode:  'copy'
    container params.CONTAINER

    input:
    path("*")
    val PREFIX

    output:
    tuple path("${PREFIX}_tidyamplicons/samples.csv"), path("${PREFIX}_tidyamplicons/taxa.csv"), path("${PREFIX}_tidyamplicons/abundances.csv")

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

    # add genomesizes and norm_readcount to abundances table
    ## build new ta object with coverage
    select_coverage(".")
    kraken2taxtable(".", "taxtable_cov", file_pattern=".mpa_cov")
    ta_cov <- import_tidyamplicons("taxtable_cov")

    ## left join to original 
    ## (! this table has different taxon_ids, as it keeps track of higher orders of lower assigned reads)
    taxid_ref <- ta\$taxa %>% 
	left_join(ta_cov\$taxa, by=c("family", "genus", "species"), suffix=c(".true", ".false")) 
    taxid_ref <- taxid_ref %>%
	select(taxon_id.true, taxon_id.false)
    # Keep only taxa from original
    ta_cov\$abundances <- ta_cov\$abundances %>% 
	filter(taxon_id %in% taxid_ref\$taxon_id.false) %>% 
	left_join(taxid_ref, by=c("taxon_id" = "taxon_id.false")) %>% 
	select(-taxon_id) %>% rename(taxon_id = taxon_id.true, coverage = abundance)

    ta\$abundances <- ta\$abundances %>% left_join(ta_cov\$abundances)

    # save the tidyamplicons object as three tidy tables
    ta %>% write_tidyamplicons("${PREFIX}_tidyamplicons")

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

process EXTRACT_PROPORTIONS {
    publishDir "${params.OUTPUT}",  mode:  'copy'
    container params.CONTAINER

    input:
    path("*")

    output:
    path("bact_proportions.tsv")


    script:
    """
    extract_bacterial_proportion.py . 
    """
}

workflow CONVERT_REPORT_TO_TA {
    take:
        report_ch

    main:
        // Convert to mpa format
        CONVERT_MPA( report_ch )
            .set { mpa_reports }

        // Normalize using genome size
        if (params.GENOMESIZES) {
        ch_genomesizes = Channel.value(file ("${params.GENOMESIZES}"))    

        NORMALIZE_READCOUNT( mpa_reports, ch_genomesizes )
            .collect()
            .set{ norm_rc }

        } else { 
            mpa_reports
                .collect{it[1]}
                .set {norm_rc}
        }

	EXTRACT_PROPORTIONS(norm_rc) 
        CREATE_TIDYAMPLICONS(norm_rc, params.REPORT_TYPE)
            .map {it.first().getParent()}
            .set { ta }

        if (params.TEST_PIPELINE){
            PRINT_TOP10(ta) | view {"$it"}
        }
}
