params.OUTPUT = "results"

process SINGLEM_PROFILE {
    tag "$meta"
    label 'process_medium'

    cpus 10

    //publishDir "${params.OUTPUT}/singleM/class", pattern: "**_classifications.tsv", mode: "copy"

    //conda "bioconda::metabuli=1.0.0"
    container "wwood/singlem:1.0.0beta7"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*_profile.csv"), emit: profile
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def single_end = ! params.pairedEnd
    
    def input = single_end ? "--sequences `pwd`/${fastas}" : "-1 `pwd`/${fastas[0]} -2 ´pwd´/${fastas[1]}"
    
    """

    pipe \\
      $args \\
      ${input} \\
      -p ${prefix}_profile.csv \\
      --threads $task.cpus
      
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        singleM: 1.0.0beta7
    END_VERSIONS
    """
}
