params.reads = "${projectDir}data/samples/*_R{1,2}_001.fastq.gz"
params.krakendb = "/mnt/ramdisk/krakendb"
params.pairedEnd = true
params.min_reads=800
params.b_treshold = 10
params.truncLen = 0
params.trimLeft = 0
params.trimRight = 0
params.minLen = 50
params.maxN = 2
params.maxEE = 2

params.debug = false

def helpMessage() {
    log.info"""
     Name: nf-kraken2-bracken
     Author: LAMB (UAntwerp)
    =========================================
    Required arguments:
      --reads                       Path to directory with input samples. If using paired reads 
                                    they need to be captured using a glob expression such as the following:
                                    data/samples/*_R{1,2}_001.fastq.gz

      --krakendb                    Path to kraken database.
    Optional arguments:

      --pairedEnd                   Specifies if reads are paired-end (true | false). Default = ${params.pairedEnd}
      --min_reads                   Minimum amount of reads needed for analysis. Default = ${params.min_size}
      --b_treshold                  
      --outdir                      The output directory where the results will be saved. Defaults to ${params.outdir}
      --help  --h                   Shows this help page

      --truncLen                    Truncation length used by dada2 FilterandTrim algorithm. Default = ${params.truncLen}
      --trimLeft --trimRight        Trimming on left or right side of reads by dada2 FilterandTrim algorithm. Default = ${params.trimLeft}
      --minLen                      Minimum length of reads kept by dada2 FilterandTrim algorithm. Default = ${params.minLen}
      --maxN                        Maximum amount of uncalled bases N to be kept by dada2 FilterandTrim algorithm. Default = ${params.maxN}
      --maxEE                       Maximum number of expected errors allowed in a read by dada2 FilterandTrim algorithm. Default = ${params.maxEE}

      --debug

    Usage example:
        nextflow run main.nf --reads '/path/to/reads' \
        --krakendb '/path/to/krakendb/' 
    """.stripIndent()
}

def paramsUsed() {
    log.info"""
    N F - K R A K E N 2 - B R A C K E N
    =========================================
    reads: ${params.reads}
    krakendb: ${params.krakendb}
    outdir: ${params.outdir}
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

process READLENGTH_DISTRIBUTION {
    tag "${pair_id}"

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("readlengths")


    script:
    def single = reads instanceof Path

    def read1 = !single ? "${reads[0]}" : "${reads}"
    def read2 = !single ? "${reads[1]}" : '' 
    """
    zcat ${read1} ${read2} | awk '{if(NR%4==2) print length(\$1)}'| sort -n | uniq -c > readlengths
    """

}

process DETERMINE_MAX_LENGTH {
    tag "${pair_id}"

    input:
    tuple val(pair_id), path(readLengths)

    output:
    tuple val(pair_id), path(readLengths), stdout


    script:
    """
    determine_kmer.py ${readLengths[0]}
    """
    
}

process FILTER_TRIM {
    tag "${pair_id}"

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("filtered_${pair_id}_*"), emit: filteredReads
    path "readCounts_${pair_id}", emit: readCounts 
    
    script:
    def single = reads instanceof Path

    def read1 = !single ? "fwd='${reads[0]}'" : "'${reads}'"
    def read2 = !single ? "rev='${reads[1]}', filt.rev='filtered_${pair_id}_rev.fastq.gz'," : ''
    def filts = !single ? "'filtered_fwd', 'filtered_rev'" : "'filtered_fwd'" 

    def truncLen = !single ? "c(${params.truncLen}, ${params.truncLen})" : ${params.truncLen}
    def trimLeft = !single ? "c(${params.trimLeft}, ${params.trimLeft})" : ${params.trimLeft}
    def trimRight = !single ? "c(${params.trimRight}, ${params.trimRight})" : ${params.trimRight}
    def minLen = !single ? "c(${params.minLen}, ${params.minLen})" : ${params.minLen}
    def maxN = !single ? "c(${params.maxN}, ${params.maxN})" : ${params.maxN}
    def maxEE = !single ? "c(${params.maxEE}, ${params.maxEE})" : ${params.maxEE}
    """
    #!/usr/bin/env Rscript

    library(dada2)
    
    # create files in case of empty reads
    file.create(${filts}, showWarnings=F)
    
    readCounts <- filterAndTrim(
        ${read1},${read2} filt='filtered_${pair_id}_fwd.fastq.gz',
        truncLen = ${truncLen},
        trimLeft = ${trimLeft},
        trimRight = ${trimRight},
        minLen = ${minLen},
        maxN = ${maxN},
        maxEE = ${maxEE}
    )

    write.csv(readCounts, "readCounts_${pair_id}")

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
    output:
    tuple val(pair_id), path("${pair_id}.kraken2.report")
    

    script:
    def single = reads instanceof Path

    def read1 = !single ? "${reads[0]}" : "'${reads}'"
    def read2 = !single ? "${reads[1]}" : ''
    def mode = !single ? "--paired" : "" 

    def report = pair_id + ".kraken2.report"
    //def out = pair_id + ".kraken.out"

    """
    kraken2 --db "${params.krakendb}" --report "${report}" --threads ${task.cpus} \
    --minimum-base-quality ${params.b_treshold} --confidence ${params.confidence} \
    --memory-mapping ${mode} "${read1}" "${read2}" > /dev/null
    """

}

process BRACKEN {
    tag "${pair_id}"
    publishDir "${params.outdir}/bracken", mode: 'copy'

    input:
    tuple val(pair_id) , path(kraken_rpt), path(readLengths), val(rlen)

    output:
    tuple val(pair_id), path("${pair_id}_bracken.report")
    //path("${pair_id}_bracken.out")
    
    script:
    """
    bracken -d ${params.krakendb} -i ${kraken_rpt} -w "${pair_id}_bracken.report" \
    -o "${pair_id}_bracken.out" -r ${rlen}    
    """

}

// process KRONA_VISUALIZATION {
//     tag "${pair_id}"
//     publishDir "${params.outdir}/krona", mode: 'copy'

//     input:
//     tuple val(pair_id), path(brk_rpt)

//     output:
//     path("${pair_id}_krona.html")

//     script:
//     """
//     ktImportTaxonomy -t 5 -m 2 -o "${pair_id}_krona.html" ${brk_rpt}
//     """
// }

process CONVERT_MPA {
    tag "${pair_id}"
    //publishDir "${params.outdir}/mpa", mode: 'link'

    input:
    tuple val(pair_id), path(brck_rpt)

    output:
    path("${pair_id}_bracken.report.mpa")

    script:
    """
    kreport2mpa.py -r ${brck_rpt} -o "${pair_id}_bracken.report.mpa"
    """
}

process CREATE_TIDYAMPLICONS {
    publishDir "${params.outdir}",  mode:  'copy'

    input:
    path("*")

    output:
    tuple path("tidyamplicons/samples.csv"), path("tidyamplicons/taxa.csv"), path("tidyamplicons/abundances.csv")


    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(tidyamplicons)

    source('$baseDir/bin/functions.R')

    run <- Sys.Date()
    pipeline <- "${launchDir.getName()}"

    # convert kraken2 results to a nice taxonomy table
    kraken2taxtable(".", "taxtable")

    ta <- import_tidyamplicons("taxtable")

    # add the run and pipeline names to the tidyamplicons object
    ta\$samples\$run <- run
    ta\$samples\$pipeline <- pipeline

    # save the tidyamplicons object as three tidy tables
    ta %>% write_tidyamplicons("tidyamplicons")

    """
}

process FASTQC {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path(reads)

    output:
    path("fastqc")

    script:
    """
    mkdir fastqc && fastqc -t ${task.cpus} -f fastq -o fastqc ${reads}
    """
}

process MULTIQC {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path('*')

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc .
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
            succes : params.pairedEnd ? it[1][1].countFastq() >= params.min_reads &&  it[1][0].countFastq() >= params.min_reads : it[1][0].countFastq() >= params.min_reads 
            failed : params.pairedEnd ? it[1][1].countFastq() < params.min_reads &&  it[1][0].countFastq() < params.min_reads : it[1][0].countFastq() < params.min_reads
        }
        .set { reads }

    reads.failed.subscribe { println "Sample ${it[0]} did not meet minimum reads requirement of ${params.min_reads} reads and is excluded."}

    // Filter and trim using dada 
    FILTER_TRIM(reads.succes)
        //.ifEmpty { error "No reads to filter" }
    
    FILTER_TRIM.out.filteredReads
        .ifEmpty { error "No reads to filter" }
        .set { filteredReads }

    WRITE_READCOUNTS(
        FILTER_TRIM.out.readCounts.collect()    
    )
        .set{ readCounts }

    // Summarize readlengths
    READLENGTH_DISTRIBUTION( filteredReads )
        .set { readLengths }
    
    // QC
    filteredReads
        .mix(reads.succes)
        .mix(reads.failed)
        .collect{ it[1] }
        .set { allReads }

    
    FASTQC(allReads)
        .set { fastqc }

    MULTIQC(fastqc)

    // Determine % reads of sizes 200, 100, 50
    DETERMINE_MAX_LENGTH(readLengths)
    // Filter out empty reads (kmer=0)
        .branch { 
            success : it[2] as int > 0
            failed : it[2] as int == 0
        }
        .set { max_length }

    // Run kraken on the samples with kmer > 0
    KRAKEN(filteredReads)
        .branch {
            success: it[1].countLines() > 1
            failed: it[1].countLines() == 1
        }
        .set { kraken_reports }
    
    kraken_reports.failed.subscribe { println "Sample ${it[0]} only contained unclassified reads and is excluded for further bracken processing."}

    // join max length info with reports
    kraken_reports.success
        .join(max_length.success)
        .set { reportsAndLengths }
    // Run bracken on appropriate kmer sizes   
    BRACKEN(reportsAndLengths)
        .set{ brck_reports }

    // Convert to mpa format
    CONVERT_MPA( brck_reports )
        .collect()
        .set { mpa_reports }

    CREATE_TIDYAMPLICONS(mpa_reports)

    // viz with krona
    //KRONA_VISUALIZATION(brck_reports)

}