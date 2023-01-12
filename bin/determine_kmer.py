#!/usr/bin/env python3
import sys
import gzip
import pandas as pd
from Bio import SeqIO

KMER_SIZES = (200, 100, 50)
readlengths = {}

with gzip.open(sys.argv[1], "rt") as f:
    for record in SeqIO.parse(f, format="fastq"):
        r_len = len(record.seq)

        if r_len not in readlengths.keys():
            readlengths[r_len] = 1
        else: readlengths[r_len] += 1

reads = pd.DataFrame(readlengths, index=[0]).transpose()
total_reads = reads[0].sum()

def determine_kmer():
    for kmer in KMER_SIZES:
        reads_this_size = reads[reads[0] >= kmer][0].sum()
        read_perc = reads_this_size / total_reads
        if read_perc >= 0.8:
            return kmer
    return kmer


kmer = determine_kmer()
print(kmer)
