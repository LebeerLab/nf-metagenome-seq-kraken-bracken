#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="aggregates taxids towards a desired rank")
parser.add_argument("sequences", metavar="SEQ", help="a tsv with taxids and read sequences")
parser.add_argument("mpa", metavar="MPA", help="the kraken report in MPA format")
parser.add_argument("--level", default="S")

args = parser.parse_args()

seqs = pd.read_table(args.sequences, names=["taxid", "sequence"])
mpa = pd.read_table(args.mpa, names = ["taxonomy","full_counts","taxid"])

RANKS = ["D", "P", "C", "O", "F", "G", "S", "S1"]
try:
    desired_rank = str(args.level)[0].upper()
    desired_i = RANKS.index(desired_rank)
except:
    raise IndexError(f"Rank name {args.level} not found in {RANKS}")

mpa[RANKS] = mpa["taxonomy"].str.split("|", expand=True)
ranks_to_keep = RANKS[:desired_i+1]
ranks_to_keep.extend(["taxonomy", "taxid"])

seqs = seqs[["taxid", "sequence"]].merge(mpa[ranks_to_keep], how='left', on="taxid").drop(columns="taxonomy")
seqs.to_csv(args.sequences, sep="\t", index=False)