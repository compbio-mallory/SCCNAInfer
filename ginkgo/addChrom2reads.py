import pandas as pd
import sys

reads = sys.argv[1]
CN = sys.argv[2]
out = sys.argv[3]

df_reads = pd.read_csv(reads, sep="\t")
df_CN = pd.read_csv(CN, sep="\t")
df_chrom = df_CN[['CHROM', 'START', 'END']]
df_reads = pd.concat([df_chrom, df_reads], axis=1)
df_reads.to_csv(out, sep="\t", index=False)