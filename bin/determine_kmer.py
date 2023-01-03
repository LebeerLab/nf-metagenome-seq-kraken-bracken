#!/usr/bin/env python3
import sys
import pandas as pd

try:
    reads = pd.read_csv(sys.argv[1], header=None, sep="\s+")
except pd.errors.EmptyDataError:
    print(0)
    exit()

total_reads = reads[0].sum()

kmer_sizes = (200, 100, 50)


def determine_kmer():
    for kmer in kmer_sizes:
        reads_this_size = reads[reads[1] >= kmer][0].sum()
        read_perc = reads_this_size / total_reads
        if read_perc >= 0.8:
            return kmer
    return kmer


kmer = determine_kmer()
print(kmer)
