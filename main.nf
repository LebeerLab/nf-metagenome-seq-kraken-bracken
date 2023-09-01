// PARAMS ======================================================================
params.reads = "${projectDir}data/samples/*_R{1,2}_001.fastq.gz"
params.profiler = "kraken"
params.krakendb = "/mnt/ramdisk/krakendb"
params.outdir = "results"

params.debug = false
params.skip_fastp = false

params.pairedEnd = true
params.min_reads=800

params.truncLen = 0
params.trimLeft = 0
params.trimRight = 0
params.minLen = 50
params.maxN = 2
params.windowFront = 4
params.windowTail = 5
params.genomesizes = null


// INCLUDE WORKFLOW ==============================================================
include { KRACKEN_BRACKEN } from './modules/kracken_bracken'

// INCLUDE MODULES ===============================================================
include { METABULI_CLASSIFY } from './modules/metabuli/classify' addParams(
    OUTPUT : "${params.outdir}"
)
include { DETERMINE_MIN_LENGTH as GET_MINLEN} from './modules/kracken_bracken'
include { FASTP; MULTIQC } from './modules/qc' addParams(
    OUTPUT: "${params.outdir}"    
)
include { CONVERT_REPORT_TO_TA; MERGE_ASV_SEQUENCES} from './modules/tidyamplicons' addParams(
    OUTPUT: "${params.outdir}", SKIP_NORM : "${params.skip_norm}", 
    GENOMESIZES : "${params.genomesizes}", TEST_PIPELINE : "${params.test_pipeline}"
)

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
      --cutFront --cutTail      Remove low quality bases from front or end of the reads.
      --windowFront/Tail        Size of the sliding window for cutting low quality bases. 
                                Set to 1 for trailing method. Default = ${params.windowFront} / ${params.windowTail}
      --tailQual --frontQual    Set a mean quality treshold for cutting low quality bases. Default = ${params.tailQual} / ${params.frontQual}

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


workflow PROFILING {
    take: reads

    main:

    if (params.profiler == "kraken") {
        KRACKEN_BRACKEN(reads)
        CONVERT_REPORT_TO_TA(KRACKEN_BRACKEN.out.reports, KRACKEN_BRACKEN.out.min_len)
        KRACKEN_BRACKEN.out.sequences
            .map{it -> it[1]}
            .collect()
            .set{ all_sequences }
	
        MERGE_ASV_SEQUENCES(CONVERT_REPORT_TO_TA.out.ta, all_sequences)

    } else if (params.profiler == "metabuli") {
        ch_MetabuliDB = Channel.value(file ("${params.metabulidb}"))
        METABULI_CLASSIFY(reads, ch_MetabuliDB)
        profile = METABULI_CLASSIFY.out.report
        minlen = GET_MINLEN(reads).map{it[2].toInteger()}.min()
        CONVERT_REPORT_TO_TA(profile, minlen)

    } else {
        error "Not a valid profiler: ${params.profiler}"
    }

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

    PROFILING(filteredReads)
    
}
