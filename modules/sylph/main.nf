process SYLPH_PROFILE {

    cpus 12
    conda 'bioconda::sylph'
    container 'staphb/sylph:latest'
    publishDir "${params.OUTDIR}/sylph", mode: 'copy'

    input:
    path reads
    path db

    output:
    path("results.tsv"), emit: profile
    //path("versions.yml"), emit: versions
    

    script:
    """
    sylph profile $db \
    -1 *_1.fastq.gz -2 *_2.fastq.gz \
    -t $task.cpus \
    -o results.tsv

    sed -s 's/filtered_//' results.tsv > results.tsv.tmp
    sed -s 's/_fwd.clean_1.fastq.gz//' results.tsv.tmp > results.tsv
    """
}

process SYLP_TO_MPA {
    publishDir "${params.OUTDIR}/sylph", mode: 'copy'
    container null

    input:
    path report

    output:
    path("mpa-results.tsv"), emit: mpa
    path("tidytacos"), emit: taco

    script:
    """
    python3 $moduleDir/bin/sylph-utils/sylph_to_taxprof.py -s $report -m $moduleDir/bin/sylph-utils/prokaryote/gtdb_r220_metadata.tsv.gz
    $moduleDir/bin/extract_strain_counts.sh .
    Rscript $moduleDir/bin/to_tidytacos.R
    python3 $moduleDir/bin/sylph-utils/merge_sylph_taxprof.py *.sylphmpa -o mpa-results.tsv --column relative_abundance
    """
}