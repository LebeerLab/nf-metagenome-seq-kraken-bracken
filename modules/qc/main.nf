process FASTP {
    tag "${pair_id}"    
    publishDir "${params.outdir}/fastp", mode: 'copy'
    
    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("filtered_${pair_id}*.fastq.gz"), emit: filteredReads
    path("${pair_id}_fastp.json"), emit: fastp

    script:
    def single = reads instanceof Path

    def input = !single ? "-i '${reads[0]}' -I '${reads[1]}'" : "-i '${reads}'"
    def output = !single ? "-o 'filtered_${pair_id}_fwd.fastq.gz' -O 'filtered_${pair_id}_rev.fastq.gz'" : "-o 'filtered_${pair_id}.fastq.gz'"

    """
    fastp ${input} ${output} --json ${pair_id}_fastp.json \\
     --length_required ${params.minLen} --trim_front1 ${params.trimLeft} \\
     --trim_tail1 ${params.trimRight} --max_len1 ${params.truncLen} \\
     --n_base_limit ${params.maxN} 
    """
}

process MULTIQC {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path('fastp/*')

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc -m fastp .
    """
}