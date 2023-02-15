#!/usr/bin/env python3

import sys
import pandas as pd


TAXA_LVLS = ["domain", "phyllum", "class", "order", "family", "genus", "species"]
TAXA_LVLS_REV = TAXA_LVLS.copy()
TAXA_LVLS_REV.reverse()


def split_taxon_id(gs_df: pd.DataFrame, levels: list) -> pd.DataFrame:
    gs_out = gs_df.copy()

    gs_out[levels] = gs_out["taxon"].str.split(";", expand=True)
    # Use index to fetch full taxon name later
    gs_out = gs_out.set_index("taxon")
    return gs_out


def aggregate_genome_sizes(genomesize):

    gs_df = pd.read_table(
        genomesize, header=None, skiprows=[0], names=["genome_size", "taxon"]
    )
    gs_out = pd.DataFrame(columns=gs_df.columns)

    gs_df_working = split_taxon_id(gs_df, TAXA_LVLS)
    gs_df_working["index"] = gs_df_working.index

    agg_methods = {}
    agg_methods["genome_size"] = "mean"
    agg_methods["index"] = "first"

    for lvl in TAXA_LVLS:
        agg_methods[lvl] = "first"

    # Start from species, aggregate and then use aggregated means for the next level
    for i, level in enumerate(TAXA_LVLS_REV):
        # Aggregate on this level
        gs_lvl = gs_df_working.groupby(level).agg(agg_methods)
        gs_lvl = gs_lvl.set_index("index")

        gs_lvl["taxon"] = gs_lvl.index
        # Cleanup and add to results
        gs_select = gs_lvl[["genome_size", "taxon"]].copy()
        gs_select["genome_size"] = gs_select["genome_size"].round()
        gs_out = pd.concat([gs_out, gs_select], axis=0)

        if level != "domain":
            # Zoom out taxa level for next round
            del agg_methods[level]
            # Split off last used level
            gs_select[["taxon", "garbage"]] = gs_select["taxon"].str.rsplit(
                ";", n=1, expand=True
            )

            gs_df_working = split_taxon_id(
                gs_select[["genome_size", "taxon"]].copy(),
                TAXA_LVLS[: len(TAXA_LVLS) - i - 1],
            )
            gs_df_working["index"] = gs_df_working.index

    return gs_out


if __name__ == "__main__":

    infile = sys.argv[1]
    outfile = infile.split(".")[0] + "_aggregated.tsv"
    df = aggregate_genome_sizes(sys.argv[1])
    df.reset_index()

    df[["genome_size", "taxon"]].to_csv(outfile, sep="\t", index=False)
