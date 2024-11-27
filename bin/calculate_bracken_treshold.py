import sys
import pandas as pd

kraken_report = sys.argv[1]
import os
rootdir = "/mnt/b/seqdata/illumina_mgs/16_20231116_Bonobos-Isala/results-profile/human/kraken"
for f in os.listdir(rootdir):
    kraken_report = os.path.join(rootdir, f)

df = pd.read_table(kraken_report, names=["coverage", "fragments_clade", "fragments", "minimizers","distinct_minimizers","level","ncbi_tax", "taxonomy"])
df_specs = df[df["level"] == "S"]
print(df_specs.describe())
df_specs[df_specs["distinct_minimizers"] > 200].to_csv("test", sep="\t", index=False)