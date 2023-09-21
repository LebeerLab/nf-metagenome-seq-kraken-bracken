process FILTER_HOST_READS {
    tag "${id}"
    container null

    cpus 6

    input:
    tuple val(id), path(fastq)
    path(index)

    output:
    tuple val(id), path("${id}_host_removed*"), emit: host_removed
    path("versions.yml"), emit: versions

    script:
    def single = fastq instanceof Path
    def db = "$index/$index.baseName"
    def input = single ? "-U ${fastq}" : "-1 ${fastq[0]} -2 ${fastq[1]}"
    // gz doesn't work due to symlinks?
    def output = single ? "${id}_host_removed.fastq" : "-1 ${id}_host_removed_1.fastq -2 ${id}_host_removed_2.fastq"
    """
    bowtie2 \\
     -x $db \\
     $input \\
     --threads $task.cpus | \\
     samtools view -bS -f 12 -F 256 > in.bam
     
    samtools fastq -@ $task.cpus $output -0 /dev/null -s /dev/null -n in.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | grep -oP '\\d\\.\\d\\.\\d\\.\\d\$';)
    END_VERSIONS    
    """

}