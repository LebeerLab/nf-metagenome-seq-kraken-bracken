#!/usr/bin/env python3

import pandas as pd
import sys


TAXA_LVLS = ["domain", "phyllum", "class", "order", "family", "genus", "species"]


def aggregate_genome_sizes(genomesize):

    gs_df = pd.read_table(
        genomesize, header=None, skiprows=[0], names=["genome_size", "taxon"]
    )

    gs_out = pd.DataFrame(columns=gs_df.columns)

    gs_df[TAXA_LVLS] = gs_df["taxon"].str.split(";", expand=True)

    gs_agg = {}
    gs_agg["genome_size"] = "mean" 
    
    for level in TAXA_LVLS:
        gs_agg[level] = "first"     
        gs_lvl = gs_df.groupby(level).agg(gs_agg)
        #.mean(numeric_only=True)
        gs_lvl["taxon"] = gs_lvl.loc[:, gs_lvl.columns != "genome_size"].agg(";".join, axis=1)
        gs_select = gs_lvl[["genome_size", "taxon"]].reset_index()
        gs_select = gs_select[["genome_size", "taxon"]]
        
        gs_select["genome_size"] = gs_select["genome_size"].round()
        gs_out = pd.concat([gs_out, gs_select], axis=0)
    
    return gs_out
    

if __name__ == "__main__":
    
    infile = sys.argv[1]
    outfile = infile.split(".")[0] + "_aggregated.tsv"
    df = aggregate_genome_sizes(sys.argv[1])
    
    df.to_csv(outfile, sep="\t", index=False)