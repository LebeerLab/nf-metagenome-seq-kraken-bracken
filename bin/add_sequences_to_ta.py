#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Join sequences to report by taxid')
parser.add_argument('taxa', type=str,
                    help='path to taxa.csv table')
parser.add_argument('sequences', type=str,
                    help='path to sequences tsv')

args = parser.parse_args()

RANKS = ["kingdom", "phyllum","class","order","family","genus", "species"]
snames = ["taxid","sequence"]
snames.extend(RANKS)

seqs = pd.read_table(args.sequences, names=snames, skiprows=1)
seq_per_taxid = seqs.drop(columns=RANKS).groupby("taxid").first()

seq_per_species = seqs.drop(columns=RANKS[:-1]).groupby("species").first()
print(seqs.groupby("species", dropna=False).first())
taxa = pd.read_csv(args.taxa, dtype={"taxid":"str"})

# S1 assigned reads won't find their taxid --> use mpa format, groupby species and match on species
taxa_seqs_taxid = taxa.merge(seq_per_taxid, left_on="taxid", right_on="taxid", how="inner")
taxa_seqs_species = taxa.merge(seq_per_species, left_on="species", right_on="species", how="inner")
taxa_seqs = pd.concat([taxa_seqs_species, taxa_seqs_taxid])

taxa_seqs.to_csv(args.taxa, index=False)
