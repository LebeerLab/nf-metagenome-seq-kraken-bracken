process SOURMASH_SKETCH {
    tag "$id"
    conda: 'bioconda:sourmash'
    container null

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("*.sig")
    //path("versions.yml")

    script:
    def input = "${reads[0]} ${reads[1]}"
    """
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund \
      --merge $id -o ${id}.sig $input
    """
}

process SOURMASH_GATHER {
    tag "$id"
    conda: 'bioconda:sourmash'
    container null

    input:
    tuple val(id), path(sign)
    path(db)

    output:
    tuple val(id), path("*.csv")
    //path("versions.yml")

    script:
    """
    sourmash gather $sign $db \
      -o ${id}_gather.csv
    """
}

process SOURMASH_TAX {
    tag "$id"
    conda: 'bioconda:sourmash'
    container null

    input:
    tuple val(id), path(gather)
    path(sql)

    output:
    tuple val(id), path("*.csv")
    //path("versions.yml")

    script:
    """
    sourmash tax annotate -g $gather -t $sql 
    """
}



