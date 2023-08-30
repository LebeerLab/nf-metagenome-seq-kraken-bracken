#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Join sequences to report by taxid')
parser.add_argument('report', type=str,
                    help='path to bracken report')
parser.add_argument('sequences', type=str,
                    help='path to sequences tsv')

args = parser.parse_args()


seqs = pd.read_table(args.sequences, names=["taxid","sequence"])

seq_per_taxid = seqs.groupby("taxid").first()

taxa = pd.read_csv(args.report, dtype={"taxid":"str"})
taxa['taxid'] = taxa['taxid']

taxa_seqs = taxa.merge(seq_per_taxid, left_on="taxid", right_on="taxid")
taxa_seqs.to_csv(args.report, index=False)
