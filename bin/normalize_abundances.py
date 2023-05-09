#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import numpy as np

TAXA_LVLS = ["domain", "phyllum", "class", "order", "family", "genus", "species"]


def add_aggregated_genome_sizes(genomesize_df):
    def split_taxon(taxon):
        tx = taxon.split("|")
        while len(tx) < len(TAXA_LVLS):
            tx.append("")
        return tx
    genomesize_df[TAXA_LVLS] = genomesize_df.apply(
        lambda x: split_taxon(x["taxon"]), axis=1, result_type='expand')
    return genomesize_df

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
    # Alternative
    df_ab_spl = add_aggregated_genome_sizes(df_ab)
    df_gs_spl = add_aggregated_genome_sizes(df_gs)
    TAXA_LVLS.reverse()
    indexes = pd.Index([])
    df_m = pd.DataFrame()
    for level in TAXA_LVLS:
        df_class = df_ab_spl.loc[df_ab_spl[level].str.len() != 0]
        df_class = df_class.drop(indexes, axis=0)
        corr_i = df_class.index
        df_class = df_class.merge(df_gs_spl[[level, "genome_size"]], how="left", on=level)
        df_class = df_class.groupby(level).mean()
        df_class.index = corr_i
        out_df = df_ab.join(df_class["genome_size"])[["taxon","readcount","genome_size"]]
        if (df_m.empty):
            df_m = out_df
        else:
            df_m["genome_size"] = df_m["genome_size"].fillna(out_df["genome_size"])
        indexes = indexes.union(df_class.index)




    # Give missing genomesizes average size of bact genome... Unless human, then use 3*10e9
    df_m.loc[df_m.taxon.str.endswith("Homo sapiens"), "genome_size"] = 3*10e9
    df_m["genome_size"] = df_m["genome_size"].fillna(38*10e5)
    df_m["norm_readcount"] = factor * df_m["readcount"] / df_m["genome_size"]

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
