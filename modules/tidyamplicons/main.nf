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

workflow CONVERT_REPORT_TO_TA {
    take:
        report_ch

    main:
        // Convert to mpa format
        CONVERT_MPA( report_ch )
            .set { mpa_reports }

        // Normalize using genome size
        if (!params.SKIP_NORM){
        ch_genomesizes = Channel.value(file ("${params.GENOMESIZES}"))    

        NORMALIZE_READCOUNT( mpa_reports, ch_genomesizes )
            .collect()
            .set{ norm_rc }

        } else { 
            mpa_reports
                .collect{it[1]}
                .set {norm_rc} 
        }

        CREATE_TIDYAMPLICONS(norm_rc, params.REPORT_TYPE)
            .map {it.first().getParent()}
            .set { ta }

        if (params.TEST_PIPELINE){
            PRINT_TOP10(ta) | view {"$it"}
        }
}