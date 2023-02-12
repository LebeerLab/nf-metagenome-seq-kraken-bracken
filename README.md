# Kraken-Bracken Pipeline

Metagenomics analysis pipeline using kraken for OTU detection and bracken for correction of detected abundance values.

## Input

You can check all input arguments by using the --help flag of the pipeline.
All inputs and options can be modified either from the command line or directly by changing their default value in the nextflow.config file.  
The following arguments are required:

### reads

Location of the input fastq files. Example:

`--reads '/data/samples/*_R{1,2}_001.fastq.gz' `
for paired end reads
`--reads 'data/samples/*.fastq.gz'`
for single end reads

The name of the path is provided in quotes and a \* glob pattern is used to find all fastq files.  
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
N E X T F L O W  ~  version 22.10.3
Launching `main.nf` [golden_wiles] DSL2 - revision: e005bf7aee
WARN: Access to undefined parameter `min_size` -- Initialise it to a default value eg. `params.min_size = some_value`

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
  --outdir                  The output directory where the results will be saved. Defaults to ./results

  --pairedEnd               Specifies if reads are paired-end (true | false). Default = true
  --min_reads               Minimum amount of reads needed for analysis. Default = null

  --truncLen                Truncation length used by fastp. Default = 0
  --trimLeft --trimRight    Trimming on left or right side of reads by fastp. Default = 0
  --minLen                  Minimum length of reads kept by fastp. Default = 50
  --maxN                    Maximum amount of uncalled bases N to be kept by fastp. Default = 2

  --b_treshold              Minimum base quality used in classification with Kraken2.
  --confidence              The confidence used in Kraken2 classfication. Default = 0

  --bracken_treshold        The minimum number of reads required for a classification at a specified rank.

Usage example:
    nextflow run main.nf --reads '/path/to/reads' --krakendb '/path/to/krakendb/'


```
