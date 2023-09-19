BASE_QUAL= 10
CONF = 0
MIN_HIT_GROUP = 2
BRACKEN_TRESH = 10
LEVEL = "S"

process KRAKEN {
    tag "${pair_id}"
    publishDir "${params.outdir}/kraken", mode: 'copy', pattern: '*[!.yaml]'

    input:
    tuple val(pair_id), path(reads)
    path db

    output:
    tuple val(pair_id), path("${pair_id}.kraken2.report"), emit: reports
    tuple val(pair_id), path("${pair_id}_classified_*.fq"), emit: classified_reads
    tuple val(pair_id), path("${pair_id}.kraken2.out"), emit: raw_output
    path("versions.yaml"), emit: versions
    
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
    --minimum-base-quality $BASE_QUAL --confidence $CONF \
    --classified-out "${classif}" --minimum-hit-groups $MIN_HIT_GROUP \
    --memory-mapping ${mode} "${read1}" "${read2}" > $outp
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version | grep version | sed 's/kraken2 version //';))
    END_VERSIONS 
    """
}

process BRACKEN {
    tag "${pair_id}"
    publishDir "${params.outdir}/bracken", mode: 'copy', pattern: '*[!.yaml]'

    input:
    tuple val(pair_id) , path(kraken_rpt), path(reads), val(min_len)
    path db

    output:
    tuple val(pair_id), path("${pair_id}.bracken.report"), emit: reports
    tuple val(pair_id), path("${pair_id}.bracken.out")
    tuple val(pair_id), val(min_len), emit: min_len 
    path("versions.yaml"), emit: versions
    //path("${pair_id}_bracken.out")
    
    script:
    def minLen = params.test_pipeline ? 100 : min_len
    """
    bracken -d ${db} -i ${kraken_rpt}  -w "${pair_id}.bracken.report" \
    -o "${pair_id}.bracken.out" -l $LEVEL -t $BRACKEN_TRESH -r ${minLen}   
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: 2.8.0
    END_VERSIONS 
    """
    // UPDATE THIS ONCE VERSION ECHO VERSION OF BRACKEN IS RELEASED...

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
        ch_versions = Channel.empty()
        KRAKEN(reads, ch_KrakenDB)
        ch_versions = ch_versions.mix(
            KRAKEN.out.versions.first()
        )
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
        ch_versions = ch_versions.mix(
            BRACKEN.out.versions.first()
        )

    emit: 
    reports = BRACKEN.out.reports
    min_len = BRACKEN.out.min_len
    versions = ch_versions
}
