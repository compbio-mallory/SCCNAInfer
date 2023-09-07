import pandas as pd
import numpy as np
import sys
import re

def split_row(row, step):
	chrom, start, end, cn = row
	rows = []
	for i in range(start, end, step):
		new_end = min(i+step, end)
		rows.append((chrom, i, new_end, cn))
	return rows


cnvFile = sys.argv[1]
f = open(cnvFile, "r")
lines = f.readlines()
f.close()
cell = re.search('for (.+?)\.bam', lines[0]).group().split()[1].split(".")[0]

curr_cn = [['CHROM', 'START', 'END', cell]]
cell2cnv = {}
for i in range(1, len(lines)):
	
	while lines[i].startswith("chr"):
		temp = lines[i].rstrip().split()
		curr_cn.append([temp[0], int(temp[1]), int(temp[2]), int(temp[3].split("-")[0])])
		i += 1
		if i > len(lines) -1:
			break	
	cell2cnv[cell] = curr_cn
	if i > len(lines) -1:
		break	
	cell = re.search('for (.+?)\.bam', lines[i]).group().split()[1].split(".")[0]
	curr_cn = [['CHROM', 'START', 'END', cell]]
		


cell2df = {}
for c in cell2cnv:
	df = pd.DataFrame(cell2cnv[c][1:], columns=cell2cnv[c][0])
	new_rows = []
	for _, row in df.iterrows():
		new_rows += split_row(row, int(sys.argv[3]))

				
# create a new dataframe with the split rows
	new_df = pd.DataFrame(new_rows, columns=['CHROM', 'START', 'END', c])
	new_df['START'] = new_df['START']  + 1
	
	cell2df[c] = new_df

merged_df = None
for key, df in cell2df.items():
	df.columns = ['CHROM', 'START', 'END', key]
	if merged_df is None:
		merged_df = df
	else:
		merged_df = pd.merge(merged_df, df, on=['CHROM', 'START', 'END'], how='outer')

# Sort the columns in the desired order
cols = ['CHROM', 'START', 'END'] + list(cell2df.keys())
merged_df = merged_df[cols]

#merged_df.dropna(inplace=True)
merged_df.to_csv(sys.argv[2], index=False, sep="\t")
