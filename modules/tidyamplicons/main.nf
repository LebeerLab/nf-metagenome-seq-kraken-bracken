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
    tuple path("tidyamplicons/samples.csv"), path("tidyamplicons/taxa.csv"), path("tidyamplicons/abundances.csv")
    path(genomesizes)
    val(readlen)

    output:
    tuple path("tidyamplicons_coverage/samples.csv"), path("tidyamplicons_coverage/taxa.csv"), path("tidyamplicons_coverage/abundances.csv"), emit: ta 

    
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
    tuple path("tidyamplicons/samples.csv"), path("tidyamplicons/taxa.csv"), path("tidyamplicons/abundances.csv"), emit: ta
    path("versions.yml"), emit: versions

    script:
    """
    bracken_to_taxtable.R "${launchDir.getName()}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | grep -Po "version \\d\\.\\d\\.\\d" | sed 's/version //';)
    END_VERSIONS  
    """
}


workflow CONVERT_REPORT_TO_TA {
    take:
        reports
        min_len

    main:
        ch_versions = Channel.empty()
        // Convert to mpa format
        reports.join(min_len)
            .set{ report_ch }
        CONVERT_MPA( report_ch )

        CONVERT_MPA.out
            .collect{ it[1] }
            .set { mpa_reports }
        
        CREATE_TIDYAMPLICONS( mpa_reports )
        ch_versions = ch_versions.mix(
            CREATE_TIDYAMPLICONS.out.versions
        )

        // Normalize using genome size

        if (params.GENOMESIZES) {
            ch_genomesizes = Channel.value(file ("${params.GENOMESIZES}"))    

            readlen = CONVERT_MPA.out
                        .map{it[2].toInteger()}.min()

            NORMALIZE_READCOUNT( CREATE_TIDYAMPLICONS.out.ta, ch_genomesizes, readlen )
            tidyamplicons = NORMALIZE_READCOUNT.out.ta

        } else {
            tidyamplicons = CREATE_TIDYAMPLICONS.out.ta
        }
    
    emit:
    ta = tidyamplicons
    versions = ch_versions

}

