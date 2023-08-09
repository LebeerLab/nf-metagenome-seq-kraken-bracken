params.b_treshold = 10
params.confidence = 0
params.min_hit_groups = 2
params.bracken_treshold = 10

process KRAKEN {
    tag "${pair_id}"
    publishDir "${params.outdir}/kraken", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)
    path db

    output:
    tuple val(pair_id), path("${pair_id}.kraken2.report")
    

    script:
    def single = reads instanceof Path

    def read1 = !single ? "${reads[0]}" : "'${reads}'"
    def read2 = !single ? "${reads[1]}" : ''
    def mode = !single ? "--paired" : "" 

    def report = pair_id + ".kraken2.report"

    """
    kraken2 --db "${db}" --report "${report}" --threads ${task.cpus} \
    --minimum-base-quality ${params.b_treshold} --confidence ${params.confidence} \
    --minimum-hit-groups ${params.min_hit_groups} --memory-mapping ${mode} "${read1}" "${read2}" > /dev/null
    """

}

process BRACKEN {
    tag "${pair_id}"
    publishDir "${params.outdir}/bracken", mode: 'copy'

    input:
    tuple val(pair_id) , path(kraken_rpt), path(reads), val(min_len)
    path db

    output:
    tuple val(pair_id), path("${pair_id}_bracken.report"), val(min_len) 
    //path("${pair_id}_bracken.out")
    
    script:
    def minLen = params.test_pipeline ? 100 : min_len
    """
    bracken -d ${db} -i ${kraken_rpt}  -w "${pair_id}_bracken.report" \
    -o "${pair_id}_bracken.out" -t ${params.bracken_treshold} -r ${minLen}     
    """

}

process DETERMINE_MIN_LENGTH {
    tag "${pair_id}"

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path(reads), stdout


    script:
    def single = reads instanceof Path

    def read1 = !single ? "${reads[0]}" : "${reads}"
    def read2 = !single ? "${reads[1]}" : '' 
    """
    determine_minlen.py ${read1} ${read2}
    """
    
}

workflow KRACKEN_BRACKEN {
    take: reads

    main:
        ch_KrakenDB = Channel.value(file ("${params.krakendb}"))

        KRAKEN(reads, ch_KrakenDB)
            .branch {
                success: it[1].countLines() > 1
                failed: it[1].countLines() == 1
            }
            .set { kraken_reports }

        kraken_reports.failed.subscribe { println "Sample ${it[0]} only contained unclassified reads and is excluded for further bracken processing."}
        // Run bracken on appropriate readlength sizes   

        // Determine % reads of sizes 200, 100, 50
        DETERMINE_MIN_LENGTH(reads)
            // Filter out empty reads (readlength=0)
            .branch { 
                success : it[2] as int > 0
                failed : it[2] as int == 0
            }
            .set { max_length }
        kraken_reports.success
            .join(max_length.success)
            .set {kraken_reportsLength}

        kraken_reportsLength
	        .map{ tuple(it[0], it[1], it[3] )}
    	    .set{ kreports }

        BRACKEN(kraken_reportsLength, ch_KrakenDB)

    emit: BRACKEN.out
}
