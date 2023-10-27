process FILTER_HOST_READS {
    tag "${id}"
    container "quay.io/biocontainers/hostile:0.1.0--pyhdfd78af_0"

    cpus 6

    input:
    tuple val(id), path(fastq)

    output:
    tuple val(id), path("*.clean*.fastq.gz"), emit: host_removed
    path("versions.yml"), emit: versions

    script:
    def single = fastq instanceof Path
    def input = single ? "--fastq1 ${fastq}" : "--fastq1 ${fastq[0]} --fastq2 ${fastq[1]}"
    """
    hostile clean \\
     $input \\
     --threads $task.cpus 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version ;)
    END_VERSIONS    
    """

}