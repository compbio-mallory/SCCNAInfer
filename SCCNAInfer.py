# this is a wrapper script to run any of Aneufinder, Ginkgo, SCOPE, and SeCNV 
# followed by the correction
import os
import argparse
import sys
import glob

src_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(src_path)
print(src_path)

bdir = ""
wd = ""
ref = "hg19"
fa = ""
pattern = "*dedup.bam"
binw = "500000"
out = ""
minP = "1.5"
maxP = "5"
K = "auto_set"
s = "auto_set"
norm = "None"
def run_secnv(minP, maxP, K, s, norm):
  src_dir = src_path + "/scripts"	
  # create output directory
  if not os.path.isdir(os.path.join(wd, "SeCNV")):
    os.system("mkdir " + os.path.join(wd, "SeCNV"))
  secnv_dir = os.path.join(wd, "SeCNV")
  # preprocess
  print("SeCNV: processing bam files and  getting coverage ...")
  cmd = "python " + src_dir + "/preprocess.py " + " ".join([secnv_dir,  fa, bdir, pattern, binw, ref])
  #print(cmd)
  os.system(cmd) 
  print("SeCNV: getting CNV and running correction ......")
  # run secnv and correction
  cmd = "python " + src_dir + "/SCCNAInfer_main.py " +  " -cov genome_cov.bed -path " + secnv_dir  + " -gc gc_map.tsv -out " + out + " -ref " + ref 
  cmd = cmd + " -minP " + minP + " -maxP " + maxP + " -K " + K + " -s " + s + " -norm "  + norm
  os.system(cmd)
  # pass test on crc1

def run_scope(pa):
  src_dir = src_path + "/scope"
  scope_dir = os.path.join(wd, "SCOPE")
  if not os.path.isdir(scope_dir):
    os.system("mkdir " + scope_dir)
  paired = ""
  if not pa:
    # single end
    paired = "-s" 
  # run scope
  cmd = "Rscript " + src_dir + "/run_scope.R -b " + bdir + " -r " + ref + " -o out " + paired
  os.system(cmd)  
  # run correction
  cmd = "python " + src_path + "/scripts/SCCNAInfer_main.py " +  " -cov " + out + "_SCOPE_reads.tsv -CN " + out + "_SCOPE_cnv.tsv -path " + scope_dir  + " -gc  gc_map.tsv  -out " + out + " -ref " + ref
  cmd = cmd + " -minP " + minP + " -maxP " + maxP + " -K " + K + " -s " + s + " -norm "  + norm 
  os.system(cmd)
  # pass test on crc1

def run_aneufinder(bw, binw):
  # aneufinder/run_aneufinder.R -b $bamdir -d $wd
  src_dir = src_path + "/aneufinder"
  an_dir = os.path.join(wd, "Aneufinder")
  if not os.path.isdir(an_dir):
    os.system("mkdir " + an_dir)
  cmd = "Rscript " + src_dir + "/run_aneufinder.R -b " + bdir + " -d " + an_dir + " -w " + str(bw)
  os.system(cmd)
  # extract CNV
  print("Getting Aneufinder CNV ...")
  cmd = "python " +  src_dir + "/get_aneufinderCN.py " +  an_dir + "/BROWSERFILES/method-edivisive/binsize_" + str(bw) + "_stepsize_" + str(bw) +"_CNV.bed.gz " + an_dir + "/" +  out + "_Aneufinder_cnv.tsv " + str(binw) 
  #print(cmd)
  os.system(cmd)
  # extract reads
  print("Extracting Raw Reads ...")
  cmd = "Rscript " + src_dir + "/getRawReads.R " + an_dir + "/binned " + an_dir + "/" + out + "_Aneufinder_reads.tsv"

  #print(cmd)
  os.system(cmd)
  # run correction 
  os.system("cp " + src_path + "/scripts/gc_map_hg19.tsv " + an_dir) 
  cmd = cmd = "python " + src_path + "/scripts/SCCNAInfer_main.py " +  " -cov " + out + "_Aneufinder_reads.tsv -CN " + out + "_Aneufinder_cnv.tsv -path " + an_dir  + " -gc gc_map_hg19.tsv  -out " + out + " -ref " + ref
  cmd = cmd + " -minP " + minP + " -maxP " + maxP + " -K " + K + " -s " + s + " -norm "  + norm 
  os.system(cmd) 
  # past test on crc2 data

def run_ginkgo(gp, gbin):
  # while read line; do
  # do something with the line
  # bedtools bamtobed -i $line > $line.bed
  # done < bam.list
  ginkgo_dir = os.path.join(wd, "Ginkgo/")
  if not os.path.isdir(ginkgo_dir):
    os.system("mkdir " + ginkgo_dir)
  files = glob.glob(os.path.join(bdir, pattern))
  print("Running Ginkgo ......")
  
  f = open(os.path.join(ginkgo_dir, "cell.list"), "w")
  print(files)
  
  for bam in files:
    base_name = os.path.basename(bam).split(".")[0]
    print("Getting bed file for bam " + base_name)
    os.system("bedtools bamtobed -i " + bam + " > " + ginkgo_dir + "/" + base_name + ".bed")
    f.write(base_name + ".bed\n")
  f.close() 
  cli/ginkgo.sh --input /gpfs/research/fangroup/lz20w/detect-CNA/T10/Ginkgo/  --genome hg19 --binning variable_500000_48_bwa  --cells cell.list
  cmd = "cd " + ginkgo_dir + "; bash " + os.path.join(gp, "cli/ginkgo.sh") + " --input " + ginkgo_dir + " --genome hg19 --binning " + gbin + " --cells " + os.path.join(ginkgo_dir, "cell.list")
  os.system(cmd)
  
  print("Running Correction ...")
  os.system("python " + src_path + "/ginkgo/addChrom2reads.py " + ginkgo_dir + "data " + ginkgo_dir + "SegCopy " + ginkgo_dir + out + "_Ginkgo_reads.tsv")
  os.system("cp " + src_path + "/scripts/gc_map_hg19.tsv " + ginkgo_dir)
  cmd = cmd = "python " + src_path + "/scripts/SCCNAInfer_main.py " +  " -cov " + ginkgo_dir + out + "_Ginkgo_reads.tsv " +  " -CN SegCopy -path " + ginkgo_dir  + " -gc gc_map_hg19.tsv  -out " + out + " -ref " + ref
  cmd = cmd + " -minP " + minP + " -maxP " + maxP + " -K " + K + " -s " + s + " -norm "  + norm
  os.system(cmd)
  # test pass on crc2 data

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  # all parameters
  parser.add_argument("-m", help = "Method: Aneufinder, Ginkgo, SCOPE, SeCNV")
  parser.add_argument("-bdir", help ="bam directory")
  parser.add_argument("-wd", help = "working directory")
  parser.add_argument("-ref", help = "reference type, hg19 or hg38")
  parser.add_argument("-fa", help = "absolute reference fasta file")
  parser.add_argument("-out", help="Output file name")
  parser.add_argument("-gc", help = "Path to the gc and maapability file", type=str)
  parser.add_argument("-binw", help = "Bin width, default is 500000", type = str, default = "500000")
  parser.add_argument("-pa", help = "Bam files patterns, default *.dedup.bam", default = "*dedup.bam", type = str)
  parser.add_argument("-n", help = "Normal cell file, default None", type = str, default = "None")
  # secnv parameters
  parser.add_argument("-minP", help="Miminum ploidy", type=str, default="1.5")
  parser.add_argument("-maxP", help="Maximum ploidy", type=str, default="5")
  parser.add_argument("-K", help = "SeCNV: The K largest distances used to construct adjacency matrix", type = str, default="auto_set")
  parser.add_argument("-s", help = "SeCNV: The standard deviation of the Gaussian kernel function", type = str, default="auto_set")
  
  # scope parameters
  parser.add_argument("-pair", help="paired end.",  action='store_true')

  # aneufinder parameters
  
  # ginkgo parameter
  parser.add_argument("-gp", help = "Abosolute path for ginkgo main directory")
  parser.add_argument("-gbin", help = "Ginkgo binning option, default is fixed_500000", default = "fixed_500000")
  args = parser.parse_args()
  # set global parameters
  args = parser.parse_args() 
  bdir = args.bdir
  wd = args.wd
  ref = args.ref
  fa = args.fa
  binw = args.binw
  pattern = args.pa
  out = args.out
  norm = args.n
  minP = args.minP
  maxP = args.maxP
  K = args.K
  s = args.s
  if args.m == "SeCNV":
    run_secnv(args.minP, args.maxP, args.K, args.s, args.n)
  if args.m == "SCOPE":
    run_scope(args.pair)
  if args.m == "Aneufinder":
    binformated = "{:.0e}".format(int(binw))
    run_aneufinder(binformated, binw)
  if args.m == "Ginkgo":
    run_ginkgo(args.gp, args.gbin)
