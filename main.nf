// PARAMS ======================================================================
params.reads = "${projectDir}data/samples/*_R{1,2}_001.fastq.gz"
params.krakendb = "/mnt/ramdisk/krakendb"
params.outdir = "results"

params.debug = false
params.skip_fastp = false

params.debug = false
params.skip_fastp = false

params.pairedEnd = true
params.min_reads=800

params.truncLen = 0
params.trimLeft = 0
params.trimRight = 0
params.minLen = 50
params.maxN = 2

params.b_treshold = 10
params.confidence = 0
params.min_hit_groups = 2
params.bracken_treshold = 10

params.genomesizes = null

// INCLUDE MODULES ===============================================================
include { FASTP; MULTIQC } from './modules/qc' addParams(OUTPUT: "${params.outdir}")
include { CREATE_TIDYAMPLICONS; PRINT_TOP10 } from './modules/tidyamplicons' addParams(
    OUTPUT: "${params.outdir}", MINLEN: "${params.minLen}", TRIMLEFT: "${params.trimLeft}", 
    TRIMRIGHT: "${params.trimRight}", TRUNCLEN: "${params.truncLen}" , MAXN: "${params.maxN}")

//======= INFO ===================================================================
def helpMessage() {
    log.info"""
     Name: nf-kraken2-bracken
     Author: LAMB (UAntwerp)
    =========================================
    Required arguments:
      --reads                   Path to directory with input samples. If using paired reads 
                                they need to be captured using a glob expression such as the following:
                                data/samples/*_R{1,2}_001.fastq.gz

      --krakendb                Path to kraken database.
    Optional arguments:

      --help  --h               Shows this help page
      --test_pipeline           Run a test of the pipeline on SRR2085099 and print the 10 most abundant taxa at the end of the pipeline.
      --debug                   Run on a small subset of samples, for debugging purposes.
      --outdir                  The output directory where the results will be saved. Defaults to ${params.outdir}

      --pairedEnd               Specifies if reads are paired-end (true | false). Default = ${params.pairedEnd}
      --min_reads               Minimum amount of reads needed for analysis. Default = ${params.min_size}
      
      --truncLen                Truncation length used by fastp. Default = ${params.truncLen}
      --trimLeft --trimRight    Trimming on left or right side of reads by fastp. Default = ${params.trimLeft}
      --minLen                  Minimum length of reads kept by fastp. Default = ${params.minLen}
      --maxN                    Maximum amount of uncalled bases N to be kept by fastp. Default = ${params.maxN}

      --b_treshold              Minimum base quality used in classification with Kraken2. Default = ${params.b_treshold}
      --confidence              The confidence used in Kraken2 classfication. Default = ${params.confidence}
      --min_hit_groups          The minimum number of hit groups needed to make a classification call. Default = ${params.min_hit_groups}

      --bracken_treshold        The minimum number of reads required for a classification at a specified rank. Default = ${params.bracken_treshold} 
      --genomesizes             A tsv containing genomesizes per taxon. Default = ${params.genomesizes}

    Usage example:
        nextflow run main.nf --reads '/path/to/reads' \
        --krakendb '/path/to/krakendb/' 
    """.stripIndent()
}

def paramsUsed() {
    log.info"""
    N F - K R A K E N 2 - B R A C K E N
    =========================================
    reads:            ${params.reads}
    krakendb:         ${params.krakendb}
    outdir:           ${params.outdir}
    
    KRAKEN: 
    b_treshold:       ${params.b_treshold}
    confidence:       ${params.confidence}
    min_hit_groups:   ${params.min_hit_groups}
    
    BRACKEN:
    bracken_treshold: ${params.bracken_treshold}

    ABUNDANCE NORMALIZATION:
    genomesizes:      ${params.genomesizes}
    """.stripIndent()
}

if (params.help  || params.h){
    helpMessage()
    exit 0
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

process WRITE_READCOUNTS {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path readcounts

    output:
    path "readcountsFiltered.csv"

    script:
    """
    awk '(NR == 1) || (FNR > 1)' ${readcounts} > readcountsFiltered.csv
    """
}

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

process CONVERT_MPA {
    tag "${pair_id}"
    //publishDir "${params.outdir}/mpa", mode: 'copy'

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
    //publishDir "${params.outdir}/norm", mode: 'copy'

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

workflow {

    paramsUsed()

    // Collect all fastq files
    Channel
        .fromFilePairs(params.reads, size: params.pairedEnd ? 2 : 1)
        .ifEmpty { error "Could not find any reads matching the pattern ${params.reads}"}
        .take( params.debug ? 3 : -1 )
        //remove 'empty' samples
        .branch {
            success : params.pairedEnd ? it[1][1].countFastq() >= params.min_reads &&  it[1][0].countFastq() >= params.min_reads : it[1][0].countFastq() >= params.min_reads 
            failed : params.pairedEnd ? it[1][1].countFastq() < params.min_reads &&  it[1][0].countFastq() < params.min_reads : it[1][0].countFastq() < params.min_reads
        }
        .set { reads }

    reads.failed.subscribe { println "Sample ${it[0]} did not meet minimum reads requirement of ${params.min_reads} reads and is excluded."}

    if (!params.skip_fastp){

        // Filter and trim using fastp
        FASTP(reads.success)

        FASTP.out.filteredReads
            .ifEmpty { error "No reads to filter"}
            .set { filteredReads }

        FASTP.out.fastp
            .collect()
            .set{fastp}

        MULTIQC(fastp)
    } else {
        filteredReads = reads.success
    }

    // Run kraken on the samples with readlength > 0
    ch_KrakenDB = Channel.value(file ("${params.krakendb}"))    

    KRAKEN(filteredReads, ch_KrakenDB)
        .branch {
            success: it[1].countLines() > 1
            failed: it[1].countLines() == 1
        }
        .set { kraken_reports }
    
    kraken_reports.failed.subscribe { println "Sample ${it[0]} only contained unclassified reads and is excluded for further bracken processing."}
    // Run bracken on appropriate readlength sizes   

    // Determine % reads of sizes 200, 100, 50
    DETERMINE_MIN_LENGTH(reads.success)
    // Filter out empty reads (readlength=0)
        .branch { 
            success : it[2] as int > 0
            failed : it[2] as int == 0
        }
        .set { max_length }

    kraken_reports.success
        .join(max_length.success)
        .set {kraken_reportsLength}

    
    BRACKEN(kraken_reportsLength, ch_KrakenDB)
        .set{ brck_reports }

    // Convert to mpa format
    CONVERT_MPA( brck_reports )
        .set { mpa_reports }

    CONVERT_MPA( kraken_reports.success)
        .set { mpa_reports_kraken }

    // Normalize using genome size
    if (!params.skip_norm){
    ch_genomesizes = Channel.value(file ("${params.genomesizes}"))    

    NORMALIZE_READCOUNT( mpa_reports, ch_genomesizes )
        .collect()
        .set{ norm_rc }

    NORMALIZE_READCOUNT( mpa_reports_kraken, ch_genomesizes)
        .collect()
        .set{ norm_rc_kr }

    } else { 
        mpa_reports
            .collect{it[1]}
            .set {norm_rc} 
        mpa_reports_kraken
            .collect{it[1]}
            .set {norm_rc_kr}
    }

    CREATE_TIDYAMPLICONS(norm_rc, "bracken")
        .map {it.first().getParent()}
        .set { ta }

    CREATE_TIDYAMPLICONS(norm_rc_kr, "kracken")
    
    if (params.test_pipeline){
        PRINT_TOP10(ta) | view {"$it"}
    }

}
