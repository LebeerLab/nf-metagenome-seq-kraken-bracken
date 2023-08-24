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
    //publishDir "${params.OUTPUT}/mpa"

    input:
    tuple val(pair_id), path(brck_rpt)

    output:
    tuple val(pair_id), path("${pair_id}_bracken.report.mpa")

    def converter = (params.profiler == "kraken") ? "kreport2mpa.py" : "mbuli2mpa.py"

    script:
    """
    ${converter} -r ${brck_rpt} -o "${pair_id}_bracken.report.mpa"
    """
}

process NORMALIZE_READCOUNT {
    publishDir "${params.OUTPUT}",  mode:  'copy'

    input:
    tuple path("tidyamplicons/samples.csv"), path("tidyamplicons/taxa.csv"), path("tidyamplicons/abundances.csv") 
    path(genomesizes)
    val(readlen)

    output:
    tuple path("tidyamplicons_coverage/samples.csv"), path("tidyamplicons_coverage/taxa.csv"), path("tidyamplicons_coverage/abundances.csv") 

    
    script:
    """
    normalize_abundances.py tidyamplicons "${genomesizes}" tidyamplicons --factor ${readlen}
    """
}

process CREATE_TIDYAMPLICONS {
    //publishDir "${params.OUTPUT}/temp",  mode:  'copy'
    container params.CONTAINER

    input:
    path(mpas)

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
        readlen

    main:
        // Convert to mpa format
        CONVERT_MPA( report_ch )
        CONVERT_MPA.out
            .collect{ it[1] }
            .set { mpa_reports }
        
        CREATE_TIDYAMPLICONS( mpa_reports )

        // Normalize using genome size

        if (params.GENOMESIZES) {
            ch_genomesizes = Channel.value(file ("${params.GENOMESIZES}"))    

            NORMALIZE_READCOUNT( CREATE_TIDYAMPLICONS.out, ch_genomesizes, readlen )

        }
}

