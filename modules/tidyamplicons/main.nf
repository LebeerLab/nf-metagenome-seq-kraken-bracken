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
    tuple val(pair_id), path(brck_rpt), val(readlen) 

    output:
    tuple val(pair_id), path("${pair_id}_bracken.report.mpa"), val(readlen) 

    script:
    """
    kreport2mpa.py -r ${brck_rpt} -o "${pair_id}_bracken.report.mpa"
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
    cat * > all_sequences
    add_sequences_to_ta.py $taxa all_sequences
    """

}


workflow CONVERT_REPORT_TO_TA {
    take:
        reports
        min_len
        sequences

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

            readlen = CONVERT_MPA.out
                        .map{it[2].toInteger()}.min()

            NORMALIZE_READCOUNT( CREATE_TIDYAMPLICONS.out.ta, ch_genomesizes, readlen )

        } else {
            taxa = CREATE_TIDYAMPLICONS.out.ta.map{ it[1] }
            seqs = sequences.collect{ it[1] }
            
            ADD_ASVS(taxa, seqs)
        }
}

