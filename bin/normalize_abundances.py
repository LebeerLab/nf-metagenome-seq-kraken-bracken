#!/usr/bin/env python3
import os
import shutil
import argparse
import pandas as pd
import numpy as np
from itertools import chain

# Take ta object and add coverages

TAXA_LVLS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
output_folder = "tidyamplicons_coverage"
human_gs = 3*10e9

def add_aggregated_genome_sizes(genomesize_df):
    def split_taxon(taxon):
        tx = taxon.split("|")
        while len(tx) < len(TAXA_LVLS):
            tx.append("")
        return tx
    genomesize_df[TAXA_LVLS] = genomesize_df.apply(
        lambda x: split_taxon(x["taxon"]), axis=1, result_type='expand')
    return genomesize_df

def normalize_readcount(taxa_f, genomesizes_f, factor):

    # Read the abundances from the MPA format table
    df_tax = pd.read_csv(taxa_f)
    df_gs = pd.read_table(
        genomesizes_f, header=None, skiprows=[0], names=["genome_size", "taxon"]
    )
    
    # Adapt to MPA style
    df_gs["taxon"] = df_gs["taxon"].apply(
        lambda x: x.replace(";", "|").replace("__", "_")
    )
    df_gs[TAXA_LVLS] = df_gs["taxon"].str.split("|", expand=True)

    taxids_assigned = []
    taxes = []
    for tax in reversed(TAXA_LVLS):
        tax_class = df_tax[~df_tax[tax].isna() & ~df_tax["taxon_id"].isin(taxids_assigned)]
        df_m = pd.merge(left=tax_class, right=df_gs[["genome_size", tax]], how="left", on=tax)
        df_m = df_m.groupby(tax).agg({"genome_size": ['mean'], "taxon_id":  ['first']}).dropna()
        taxi = list(chain.from_iterable(df_m.taxon_id.values))
        taxids_assigned.extend(taxi)
        taxes.append(df_m)
      
    df_m = pd.concat(taxes).reset_index(drop=True)
    df_m.columns = df_m.columns.droplevel(1)

    df_tax_gs = pd.merge(left=df_tax, right=df_m, how="left", left_on="taxon_id", right_on="taxon_id")
    # Give missing genomesizes average size of bact genome... Unless human, then use 3*10e9
    df_tax_gs.loc[df_tax_gs["domain"] == "d_Eukarya", "genome_size"] = human_gs
    df_tax_gs["genome_size"] = df_tax_gs["genome_size"].fillna(38*10e5)
    df_tax_gs[TAXA_LVLS] = df_tax_gs[TAXA_LVLS].fillna("")
    return df_tax_gs

def normalize_tidyamplicons(ta_folder, genomesizes_folder, factor, save=True):
    
    taxa_f = os.path.join(ta_folder, "taxa.csv")
    abund_f = os.path.join(ta_folder, "abundances.csv")
    taxa_fout = os.path.join(output_folder, "taxa.csv")
    abund_fout = os.path.join(output_folder, "abundances.csv")
    
    taxa = normalize_readcount(taxa_f,genomesizes_folder,factor)

    ab = pd.read_csv(abund_f)
    ab_m = pd.merge(left=ab, right=taxa[["taxon_id", "genome_size"]], on="taxon_id")
    ab_m["coverage"] = factor * ab_m["abundance"].astype(np.float32) / ab_m["genome_size"].astype(np.float32)
    
    if save:
        taxa.to_csv(taxa_fout, index=False)
        ab_m.to_csv(abund_fout, index=False)    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("tidyamplicons")
    parser.add_argument("genomesizes")
    parser.add_argument("output")
    parser.add_argument("-f", "--factor", default=200)

    args = parser.parse_args()
    os.mkdir(output_folder)
    shutil.copy("tidyamplicons/samples.csv", os.path.join(output_folder, "samples.csv"))
    df = normalize_tidyamplicons(args.tidyamplicons, args.genomesizes, np.float32(args.factor), save=True)
