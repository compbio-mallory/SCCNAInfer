# This file calculate the gc and mappability given bins
import os
import sys
import pandas as pd
import linecache

work_path = sys.argv[1]
ref = sys.argv[2]
ref_type = sys.argv[3]
bed_file = sys.argv[4]
'''
Bed_file can be the final tab-separated copy number file
with CHROM, START, and END as first three columns. 
'''
def chrom_key(chrom):
  if chrom == 'chrX':
    return 23
  elif chrom == 'chrY':
    return 24
  else:
    return int(chrom[3:])

def extract_bins(bed_file):
  df = pd.read_csv(bed_file, sep = "\t")
  df = df.iloc[:, :3]
  df.to_csv(os.path.join(work_path, 'consecutive_bins.tsv'), sep='\t', header=False, index=False)

# Sort the DataFrame
#  df['CHROM'] = pd.Categorical(df['CHROM'], sorted(df['CHROM'].unique(), key=chrom_key), ordered=True)

# Sort the DataFrame by 'CHROM' and then by 'START'
#  df = df.sort_values(by=['CHROM', 'START'])
  df['count'] = range(1, len(df) + 1)
  df.to_csv(os.path.join(work_path, 'consecutive_bins_add.tsv'), sep='\t', header=False, index=False)

def cal_gc():
  print("Calculating GC content...")
  genome_gc = os.path.join(work_path, "genome_gc.bed")
  bed_file = os.path.join(work_path, 'consecutive_bins.tsv')
  os.system("bedtools nuc -fi %s -bed %s | cut -f 1-3,5 > %s"%(ref, bed_file, genome_gc))

def cal_map():
  print("Calculating mappability...")
  genome_map = os.path.join(work_path, "genome_mappability.tab")
  bed_file = os.path.join(work_path, 'consecutive_bins_add.tsv')
  if ref_type == "hg19":
    os.system("bigWigAverageOverBed hg19_mappability.bigWig %s %s"%(bed_file, genome_map))
  elif ref_type == "hg38":
    os.system("bigWigAverageOverBed hg38_mappability.bigWig %s %s"%(bed_file, genome_map))

def merge_gc_map():
  map_file = os.path.join(work_path, "genome_mappability.tab")
  gc_file = os.path.join(work_path, "genome_gc.bed")
  bin_num = len(open(map_file, 'r').readlines())
  f_gc_map = open(os.path.join(work_path,"gc_map.tsv"), "w")
  f_gc_map.write("CHROM\tSTART\tEND\tgc\tmap\n")
  for i in range(bin_num):
      gc_line = linecache.getline(gc_file, i+2).strip("\n").split("\t")
      mappability = float(linecache.getline(map_file, i+1).strip("\n").split("\t")[5])
      f_gc_map.write(gc_line[0] + "\t" + gc_line[1] + "\t" + gc_line[2] + "\t" +
                          gc_line[3] + "\t" + str(mappability) + "\n")
  f_gc_map.close()  

extract_bins(bed_file)
cal_gc()
cal_map()
merge_gc_map()
