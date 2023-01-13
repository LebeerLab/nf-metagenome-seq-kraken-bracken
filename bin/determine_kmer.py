#!/usr/bin/env python3
import os
import sys
import gzip
import pandas as pd
from Bio import SeqIO

KMER_SIZES = (200, 100, 50)
readlengths = {}

def fetch_readlengths(fastq_f:str) -> pd.DataFrame:

    if not os.path.exists(fastq_f):
        raise FileNotFoundError(f"Could not locate {fastq_f}.")

    with gzip.open(sys.argv[1], "rt") as f:
        for record in SeqIO.parse(f, format="fastq"):
            r_len = len(record.seq)

            if r_len not in readlengths.keys():
                readlengths[r_len] = 1
            else: readlengths[r_len] += 1

    return pd.DataFrame(readlengths, index=[0]).transpose()


def determine_kmer():
    for kmer in KMER_SIZES:
        reads_this_size = reads[reads[0] >= kmer][0].sum()
        read_perc = reads_this_size / total_reads
        if read_perc >= 0.8:
            return kmer
    return kmer

shortest_kmer = 200
# loop over args (in case of paired reads)
for fastq_f in sys.argv[1:]:
    reads = fetch_readlengths(fastq_f)
    total_reads = reads[0].sum()
    kmer = determine_kmer()
    if kmer < shortest_kmer:
        shortest_kmer = kmer
# output to stdout
print(kmer)
