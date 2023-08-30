import sys
import pandas as pd
from Bio import SeqIO
import skbio.io

test_reads_r = "/mnt/c/seqdata_illumina_amplicon/runs_2023/10_20230804_bonobo-3/results-mgs-pipeline/kraken/reads/BVF-065-4_S330_classified"
test_reads_1 = f'{test_reads_r}__1.fq'
test_reads_2 = f'{test_reads_r}__2.fq'
test_report = "/mnt/c/seqdata_illumina_amplicon/runs_2023/10_20230804_bonobo-3/results-mgs-pipeline/kraken/BVF-065-4_S330.kraken2.report"

taxonomy_file = "/media/ramdisk/krakendb/taxonomy/nodes.dmp"

test_report_bracken = "/mnt/c/seqdata_illumina_amplicon/runs_2023/10_20230804_bonobo-3/results-mgs-pipeline/bracken/BVF-065-4_S330.bracken.report"
test_output_bracken = "/mnt/c/seqdata_illumina_amplicon/runs_2023/10_20230804_bonobo-3/results-mgs-pipeline/bracken/raw/BVF-065-4_S330.bracken.out"

classif = pd.read_table(test_report_bracken, names=["total_perc", "reads_tot", "reads_spec", "level", "taxID", "taxName"])
classif["sequence"] = ""

with open(taxonomy_file) as fh:
    nodes = skbio.io.read(fh, format="taxdump", into=pd.DataFrame, scheme="nodes_slim")

def extract_asvs_from_reads_file(reads_file, classif, taxid_field="taxID"):
    
    records = SeqIO.to_dict(SeqIO.parse(reads_file, "fastq"))

    for _, item in records.items():
        taxid = int(item.description.split("kraken:taxid|")[1])
        sequence = item.seq
        
        if len(sequence) > len(classif.loc[classif[taxid_field] == taxid, ["sequence"]]):
            classif.loc[classif[taxid_field] == taxid, ["sequence"]] = str(sequence)
    
    return classif

for f in [test_reads_1, test_reads_2]:
    classif = extract_asvs_from_reads_file(f, classif)
    
classif = classif.merge(nodes["parent_tax_id"], how='left', left_on='taxID', right_index=True)

classif.to_csv("temp.csv", index=False)

for f in [test_reads_1, test_reads_2]:
    classif = extract_asvs_from_reads_file(f, classif, taxid_field="parent_tax_id")
    
print(classif)