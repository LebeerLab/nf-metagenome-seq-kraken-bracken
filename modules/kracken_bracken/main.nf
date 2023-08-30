params.b_treshold = 10
params.confidence = 0
params.min_hit_groups = 2
params.bracken_treshold = 10
params.level = "S"

process KRAKEN {
    tag "${pair_id}"
    publishDir "${params.outdir}/kraken", mode: 'copy', pattern: "*.kraken2.report"
    //publishDir "${params.outdir}/kraken/raw", mode: 'copy', pattern: "*.kraken2.out"
    //publishDir "${params.outdir}/kraken/reads", mode: 'copy', pattern: "*.fq"
    

    input:
    tuple val(pair_id), path(reads)
    path db

    output:
    tuple val(pair_id), path("${pair_id}.kraken2.report"), emit: reports
    tuple val(pair_id), path("${pair_id}_classified_*.fq"), emit: classified_reads
    tuple val(pair_id), path("${pair_id}.kraken2.out"), emit: raw_output

    script:
    def single = reads instanceof Path

    def read1 = !single ? "${reads[0]}" : "'${reads}'"
    def read2 = !single ? "${reads[1]}" : ''
    def mode = !single ? "--paired" : "" 

    def report = pair_id + ".kraken2.report"
    def classif = pair_id + "_classified_#.fq"
    def outp = pair_id + ".kraken2.out"   

    """
    kraken2 --db "${db}" --report "${report}" --threads ${task.cpus} \
    --minimum-base-quality ${params.b_treshold} --confidence ${params.confidence} \
    --classified-out "$classif" --minimum-hit-groups ${params.min_hit_groups} \
    --memory-mapping ${mode} "${read1}" "${read2}" > $outp
    """

}

process BRACKEN {
    tag "${pair_id}"
    publishDir "${params.outdir}/bracken", mode: 'copy', pattern: "*.bracken.report"
    //publishDir "${params.outdir}/bracken/raw", mode: 'copy', pattern: "*.bracken.out"
    

    input:
    tuple val(pair_id) , path(kraken_rpt), path(reads), val(min_len)
    path db

    output:
    tuple val(pair_id), path("${pair_id}.bracken.report"), emit: reports 
    tuple val(pair_id), path("${pair_id}.bracken.out"), emit: raw_output
    tuple val(pair_id), val(min_len), emit: min_len

    script:
    def minLen = params.test_pipeline ? 100 : min_len
    """
    bracken -d ${db} -i ${kraken_rpt}  -w "${pair_id}.bracken.report" \
    -o "${pair_id}.bracken.out" -l ${params.level} -t ${params.bracken_treshold} -r ${minLen}   
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


process EXTRACT_SEQUENCES {
    tag "${pair_id}"
    publishDir "${params.outdir}/bracken/sequences", mode: 'copy'

    input:
    tuple val(pair_id), path(kraken2_out), path(bracken_report), path(classified_reads)

    output:
    tuple val(pair_id), path("*_sequences")


    script:
    
    def single = classified_reads instanceof Path

    def read1 = !single ? "${classified_reads[0]}" : "${classified_reads}"
    
    """
    tail +2 $bracken_report | cut -f 5 > taxids
    kraken-extract.py --kraken $kraken2_out --taxids taxids $read1 > "${pair_id}_sequences"
    """
}



workflow KRACKEN_BRACKEN {
    take: reads

    main:
        ch_KrakenDB = Channel.value(file ("${params.krakendb}"))

        KRAKEN(reads, ch_KrakenDB)
	    KRAKEN.out.reports
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

        KRAKEN.out.raw_output.join(BRACKEN.out.reports).join(KRAKEN.out.classified_reads)
            .set{ sequences_input }

        EXTRACT_SEQUENCES(sequences_input)

    emit: 
    reports = BRACKEN.out.reports
    min_len = BRACKEN.out.min_len
    sequences = EXTRACT_SEQUENCES.out
}
