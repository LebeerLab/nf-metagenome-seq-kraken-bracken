params.OUTPUT = "results"

process METABULI_CLASSIFY {
    tag "$meta"
    label 'process_medium'

    cpus 20

    publishDir "${params.OUTPUT}/metabuli/class", pattern: "*/*_classifications.tsv", mode: "copy"
    publishDir "${params.OUTPUT}/metabuli/reports", pattern: "*/*report.tsv", mode: "copy"
    publishDir "${params.OUTPUT}/metabuli/krona", pattern: "*/*_krona.html", mode: "copy"
    

    conda "bioconda::metabuli=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabuli:1.0.0--pl5321hf1761c0_0':
        'quay.io/biocontainers/metabuli:1.0.1--pl5321h6a68c12_0' }"

    input:
    tuple val(meta), path(fastas)
    path(db)

    output:
    tuple val(meta), path("*/*_classifications.tsv"), emit: classification
    tuple val(meta), path("*/*_report.tsv"), emit: report
    tuple val(meta), path("*/*_krona.html"), emit: krona
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def single_end = ! params.pairedEnd
    def is_compressed = single_end ? fastas.getName().endsWith(".gz") : fastas[0].getName().endsWith(".gz")
    def input = single_end ? "--seq-mode 1 ${fastas}" : "${fastas[0]} ${fastas[1]}"
    if (is_compressed && single_end) {
      input = "--seq-mode 1 ${fastas.baseName}"
    } else if (is_compressed) {
      input =  "${fastas[0].baseName} ${fastas[1].baseName}"
    }
    
    """
    if [ "$is_compressed" == "true" ]; then
    gzip -d *.gz
    fi

    metabuli \\
        classify \\
        $args \\
        ${input} \\
        ${db} \\
        ${prefix}_out \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabuli: \$(metabuli | grep Version | sed 's/^metabuli Version: //';))
    END_VERSIONS
    """
}
