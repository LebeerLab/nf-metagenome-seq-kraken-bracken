// PARAMETERS
params.CONTAINER = "theoaphidian/kraken-bracken"
params.OUTPUT = "results"
params.SKIP_NORM = true
params.GENOMESIZES = null
params.TEST_PIPELINE = false
params.readlen = 50

process CONVERT_MPA {
    tag "${pair_id}"
    //publishDir "${params.OUTPUT}/mpa"

    input:
    tuple val(pair_id), path(brck_rpt), val(readlen)

    output:
    tuple val(pair_id), path("${pair_id}_bracken.report.mpa"), val(readlen)

    def converter = (params.profiler == "kraken") ? "kreport2mpa.py" : "mbuli2mpa.py"

    script:
    """
    ${converter} -r ${brck_rpt} -o "${pair_id}_bracken.report.mpa"
    """
}

process NORMALIZE_READCOUNT {
    publishDir "${params.OUTPUT}",  mode:  'copy'

    input:
    tuple path("tidyamplicons/samples.csv"), path("tidyamplicons/taxa.csv"), path("tidyamplicons/abundances.csv"), emit: ta
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
    publishDir "${params.OUTPUT}",  mode:  'copy', pattern: "tidyamplicons/*"
    //publishDir "${params.OUTPUT}/tidyamplicons",  mode:  'copy', pattern: "taxtable"
    container params.CONTAINER
    
    input:
    path(mpas)

    output:
    tuple path("tidyamplicons/samples.csv"), path("tidyamplicons/taxa.csv"), path("tidyamplicons/abundances.csv"), emit: ta
    path("taxtable"), emit: taxtable

    script:
    """
    bracken_to_taxtable.R "${launchDir.getName()}"
    """
}

process ADD_ASVS{
    publishDir "${params.OUTPUT}/tidyamplicons",  mode:  'copy'
    
    input:
    path(taxa)
    path(sequences)

    output:
    path("taxa.csv")
    
    script:
    """
    cat *_sequences > all_sequences
    add_sequences_to_ta.py $taxa all_sequences
    """

}


workflow CONVERT_REPORT_TO_TA {
    take:
        reports
        min_len

    main:
        // Convert to mpa format

        reports.join(min_len)
            .set{reports_len}
        
        CONVERT_MPA( reports_len )

        CONVERT_MPA.out
            .collect{ it[1] }
            .set { mpa_reports }
        
        CREATE_TIDYAMPLICONS( mpa_reports )

        // Normalize using genome size

        if (!params.SKIP_NORM) {
            ch_genomesizes = Channel.value(file ("${params.GENOMESIZES}"))    

            NORMALIZE_READCOUNT( CREATE_TIDYAMPLICONS.out, ch_genomesizes, min_len )
            tidyamplicons = NORMALIZE_READCOUNT.out.ta
        } else {
            tidyamplicons = CREATE_TIDYAMPLICONS.out.ta
        }
    
    emit:
    ta = tidyamplicons

}

workflow MERGE_ASV_SEQUENCES {
    take:
    tidyamplicons
    sequences

    main:
    taxa = tidyamplicons.map{ it[1] }
    ADD_ASVS(taxa, sequences)

}
