#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd

TAXA_LVLS = ["domain", "phyllum", "class", "order", "family", "genus", "species"]


def add_aggregated_genome_sizes(genomesize_df):
    def split_taxon(taxon):
        tx = taxon.split("|")
        while True:
            if len(tx) < len(TAXA_LVLS):
                tx.append("")
            else:
                return tx

    genomesize_df[TAXA_LVLS] = genomesize_df["taxon"].apply(split_taxon)


def normalize_readcount(abundances, genomesizes, factor):

    factor = int(factor)

    # Read the abundances from the MPA format table
    df_ab = pd.read_table(abundances, header=None, names=["taxon", "readcount"])
    
    df_gs = pd.read_table(
        genomesizes, header=None, skiprows=[0], names=["genome_size", "taxon"]
    )
    
    # Adapt to MPA style
    df_gs["taxon"] = df_gs["taxon"].apply(
        lambda x: x.replace(";", "|").replace("__", "_")
    )
    df_m = df_ab.merge(df_gs, how="left", on="taxon")

    df_m["norm_readcount"] = factor * df_m["readcount"] / df_m["genome_size"]
    # Give missing genomesizes average size of bact genome... Unless human, then use 3*10e9
    df_m.loc[df_m.taxon.str.endswith("Homo sapiens"), "norm_readcount"] = factor * df_m.readcount / 3*10e9
    df_m["norm_readcount"] = df_m["norm_readcount"].fillna(factor * df_m["readcount"] /38*10e5)

    return df_m


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("abundances")
    parser.add_argument("genomesizes")
    parser.add_argument("output")
    parser.add_argument("-f", "--factor", default=200)

    args = parser.parse_args()

    df = normalize_readcount(args.abundances, args.genomesizes, args.factor)
    outf = args.output
    df.to_csv(outf, sep="\t", index=False)
