// PARAMETERS
params.CONTAINER = "theoaphidian/kraken-bracken"
params.OUTPUT = "results"
params.SKIP_NORM = true
params.GENOMESIZES = null
params.TEST_PIPELINE = false
params.REPORT_TYPE = "bracken"
params.readlen = 50

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

    input:
    tuple path("tidyamplicons/samples.csv"), path("tidyamplicons/taxa.csv"), path("tidyamplicons/abundances.csv") 
    path(genomesizes)

    output:
    tuple path("tidyamplicons_coverage/samples.csv"), path("tidyamplicons_coverage/taxa.csv"), path("tidyamplicons_coverage/abundances.csv") 

    
    script:
    """
    normalize_abundances.py tidyamplicons "${genomesizes}" tidyamplicons --factor ${params.readlen}
    """
}

process CREATE_TIDYAMPLICONS {
    publishDir "${params.OUTPUT}",  mode:  'copy'
    container params.CONTAINER

    input:
    path("*.mpa")

    output:
    tuple path("tidyamplicons/samples.csv"), path("tidyamplicons/taxa.csv"), path("tidyamplicons/abundances.csv")

    script:
    """
    bracken_to_taxtable.R "${launchDir.getName()}"
    """
}


workflow CONVERT_REPORT_TO_TA {
    take:
        report_ch

    main:
        // Convert to mpa format
        CONVERT_MPA( report_ch )
            .collect{ it[1] }
            .set { mpa_reports }
        
        CREATE_TIDYAMPLICONS( mpa_reports )

        // Normalize using genome size
        if (params.GENOMESIZES) {
        ch_genomesizes = Channel.value(file ("${params.GENOMESIZES}"))    

        NORMALIZE_READCOUNT( CREATE_TIDYAMPLICONS.out, ch_genomesizes )
            .collect()
            .set{ norm_rc }

        } else { 
        mpa_reports
            .collect{it[1]}
            .set {norm_rc}
        }

}
