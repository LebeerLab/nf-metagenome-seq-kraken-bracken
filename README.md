# Kraken-Bracken Pipeline

Metagenomics analysis pipeline using kraken for OTU detection and bracken for correction of detected abundance values.

## Input
 
You can check all input arguments by using the --help flag of the pipeline.
All inputs and options can be modified either from the command line or directly by changing their default value in the nextflow.config file.  
The following arguments are required:

### reads  
Location of the input fastq files. Example:  

```--reads '/data/samples/*_R{1,2}_001.fastq.gz' ```
for paired end reads
```--reads 'data/samples/*.fastq.gz'```
for single end reads

The name of the path is provided in quotes and a * glob pattern is used to find all fastq files.  
It is possible to specify whether the reads are paired or with the parameter pairedEnd. 
If they are paired, it is necessary to use **{1,2}** notation to specify read pairs.  

### krakendb  
Path to Kraken database. Before execution of the pipeline it is wise to copy the database to a ramdisk to improve read/write speed.

## Output

### ReadCountsFiltered

### Kraken

### Bracken

## Help
```
$ nextflow run main.nf --help
N E X T F L O W  ~  version 22.04.5

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
      --outdir                      The output directory where the results will be saved. Defaults to ${params.outdir}
      --help  --h                   Shows this help page

      --truncLen                    Truncation length used by dada2 FilterandTrim algorithm.
      --trimLeft --trimRight        Trimming on left or right side of reads by dada2 FilterandTrim algorithm.
      --minLen                      Minimum length of reads kept by dada2 FilterandTrim algorithm.
      --maxN                        Maximum amount of uncalled bases N to be kept by dada2 FilterandTrim algorithm.
      --maxEE                       Maximum number of expected errors allowed in a read by dada2 FilterandTrim algorithm. 

    Usage example:
        nextflow run main.nf --reads '/path/to/reads' \
        --krakendb '/path/to/krakendb/' 

```
