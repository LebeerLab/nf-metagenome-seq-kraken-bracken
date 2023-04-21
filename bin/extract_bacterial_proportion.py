#!/usr/bin/env python3

import sys
import os
import pandas as pd

proportions = []
unassigned = []
samples = []

def sane_perc(fl):
	try:
		return round(float(fl)/100, 3)
	except:
		return 0

working_dir = sys.argv[1]

for f in os.listdir(working_dir):
	if f.endswith('.kraken2.report'):
		filepath = os.path.join(working_dir, f)
		content = pd.read_table(filepath, names=["proportion", "reads_assigned",
					  "reads_assigned_lowest","class_level", 
					  "bs", "taxon"])
		# Damn whitespace...
		content.taxon = content.taxon.str.strip()
		proportion = content.loc[content.taxon == "Bacteria", "proportion"]
		sample = f.split(".")[0]
		unassigned_reads = content.loc[content.taxon == "unclassified", "proportion"]
		unassigned.append(sane_perc(unassigned_reads))
		proportions.append(sane_perc(proportion))
		samples.append(sample)

overview_bact_proportions = pd.DataFrame({
	"sample":samples, 
	"bacterial_proportion":proportions, 
	"unassigned_proportion":unassigned
})
overview_bact_proportions.to_csv("bact_proportions.tsv", sep="\t", index=False)
