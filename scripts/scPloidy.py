
import os
import sys
import numpy as np
import pandas as pd
from secnv_segmentation import *
#import secnv_segmentation
from sklearn.neighbors import KernelDensity
from sklearn import metrics
import warnings
import copy
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import pairwise_distances
from collections import Counter
import numpy as np
import argparse
from collections import Counter

# find the common breakpoints in non-cross sample segmentaiton
# select the breakpoint as common breakpoints if
# 1. more than 3 cells have it, or
# 2. 
# avg_reads take care of non-cross sample
# get_CN use avg_reads, so non-cross sample is taken care as well
builtin = True

src_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(src_path)

def read_matrix(file_name):
  f = open(file_name)
  matrix = []
  chr_name = []
  bin_list = []
  count = 0
  for line in f:
    line = line.strip("\n").rstrip("\t").split("\t")
    if len(line) < 4:
      continue
    count += 1
    if not count == 1:
      matrix.append(line[3:])
      chr_name.append(line[0])
      #bin_list.append(line[0]+":"+line[1]+"-"+line[2])
      bin_list.append([line[0], int(float(line[1])), int(float(line[2]))])
    else:
      sample_list = line
  matrix = np.array(matrix)
  #bin_list = np.array(bin_list)
  bin_list_df = pd.DataFrame(bin_list, columns = ['CHROM', 'START', 'END'])
  sample_list = np.array(sample_list[3:])
  matrix = matrix.astype(float)
  matrix = matrix.astype(int)
  for i in range(len(sample_list)):
    sample_list[i] = str(sample_list[i]).split("/")[-1]
  chr_name = np.array(chr_name)
  f.close()
  cov_df = pd.read_csv(file_name, sep = "\t")
  print(cov_df.head())
  cov_df = cov_df.dropna(axis=1, how='all')
  cov_df.dropna(inplace=True)
  sample_list = cov_df.columns.to_list()[3:]
  return matrix, cov_df, chr_name, bin_list_df, sample_list

def save_matrix(matrix, bin_list_df, sample_list, cnv_name, meta_name):
  ploidy = np.average(matrix, axis=1)
  ploidy = np.round(ploidy, 2)
  #df = pd.DataFrame(matrix.T, columns=sample_list, index=bin_list)
  df = pd.DataFrame(matrix.T, columns=sample_list)
  #ind = df.index.to_series().str.split(r':|-', expand = True)
  #ind = ind.rename(columns={0: "CHROM", 1: "START", 2:"END"})
  #df = pd.concat([ind, df], axis=1, join='inner')
  df = pd.concat([bin_list_df, df], axis=1, join='inner')
  df.to_csv(cnv_name, index=False, sep = "\t")
  df_2 = pd.DataFrame(np.vstack((sample_list, ploidy)).T, columns=["cell_id", "c_ploidy"])
  df_2.to_csv(meta_name, index=False)
  return df


# this function average the read counts to the median for each segment
# take care of 2d list of bps
def avg_reads(Y, BPs):
  Y_cp = copy.deepcopy(Y)
  if not isinstance(BPs[0], list):
    for j in range(len(BPs)):
      if j == 0:
        Y_cp[0:BPs[j]+1, ] = np.median(Y_cp[0:BPs[j]+1, ], axis = 0)
      else:
        Y_cp[BPs[j-1]+1:BPs[j]+1, ] = np.median(Y_cp[BPs[j-1]+1:BPs[j]+1, ], axis = 0)
  else:
    for c in range(Y_cp.shape[1]):
      bps = BPs[c]
      for j in range(len(bps)):
        if j == 0:
          Y_cp[0:bps[j]+1, c] = np.median(Y_cp[0:bps[j]+1, c])
        else:
          Y_cp[bps[j-1]+1:bps[j]+1, c] = np.median(Y_cp[bps[j-1]+1:bps[j]+1, c])
  return Y_cp
  
# get non interger copy number
# @Y: read count matrix
# @ploidy: list of ploidy
# @breakpoints: list of breakpoint for finding median read counts for each seg
def get_CN(Y, ploidy, BPs):
  # Y: row: bins, col: cells
  # get median read count for each segments
  reads = avg_reads(Y, BPs)
  res = reads*ploidy/np.mean(reads, axis=0)
  res[res > 10] = 10
  return res

# get pairwise ploidy difference
def ploidy_dist(ploidy):
  n = len(ploidy)
  dist = np.empty([n,n])
  for i in range(n):
    for j in range(n):
      dist[i,j] = abs(ploidy[i] - ploidy[j])
  return dist


def get_entropy(CN_nonInt,  ploidy, cluster, BPs):
  print("get entropy")
  k = max(cluster) + 1
  CN = np.around(CN_nonInt).T
  absCNsum = np.sum(np.abs(np.around(CN_nonInt) - CN_nonInt)) 
  #print(cluster)
  intraCluster = 0
  pseudo = []
  M_P = []
  interCluster = 0

  for c in range(max(cluster) + 1):
    cells = np.where(cluster == c)
    # get median cp number
    M_CN = np.median(CN[cells], axis = 0)
    M_P.append(ploidy[cells[0][0]])
    pseudo.append(M_CN)

    for cell in cells:
      intraCluster += np.sum(np.abs(M_CN - CN[cell]))
  # get inter cluster distance
  pseudo = np.array(pseudo)
  interCluster = pairwise_distances((pseudo), metric="manhattan")
  if max(cluster) == 0:
    return np.log( absCNsum + intraCluster)
  distSum = 0
  for c in range(max(cluster) + 1):
    cells = np.where(cluster == c)
    interCluster[c][c] = np.Inf
    #print(interCluster[c])
    c1 = interCluster[c].argmin()

    distSum += np.sum(np.abs(pseudo[c1] - pseudo[c]))
  distSum = distSum
  E = absCNsum*(1+0.05*(k))  - distSum +  intraCluster*(1+0.05*k)
  if E <= 0:
    E = 10**-20
  return np.log( E)

def init_ploidy(Y, BPs, minP, maxP):
  reads = avg_reads(Y, BPs)
  Y_cp = reads.T # row = cells, col = bins
  X = np.arange(minP, maxP, 0.05)
  initP = []
  for cell in Y_cp:
    RCNP = cell  / (np.mean(cell))
    sose_list = []
    for i in X:
      scnp = RCNP * i
      sose_list.append(np.sum(np.square(np.round(scnp) - scnp)))
    ploidy = X[sose_list.index(min(sose_list))]
    #print(sample_list[j], "mean", np.mean(cell), "median", np.median(cell),"std", np.std(cell), ploidy)
    initP.append(ploidy)
  return np.array(initP)

# get breakpoints from all chromosome
# @bins_len, a list contains the number of bins in each chromosome
# @breakpoints, a list of breakpoints index for each chromosome
# return a list of breakpoint index with respect to all bins across all chromosome
def getBPs(bins_len, breakpoints):
  bps = breakpoints[0]
  # merge breakpoints into a single list and change index
  total_len = bins_len[0]
  for i in range(1, len(breakpoints)):
    bp = np.array(breakpoints[i]) + total_len
    bps = bps + bp.tolist()
    total_len += bins_len[i]
  return bps

# matrix: seg * cell
def getSim(matrix, BPs, CN, ploidy_list):
  dist, R = getDist(matrix, BPs, CN, ploidy_list)
  return 1 - dist

# get the distance between each pair of cells
# @matrix, read count matrix. bins*cells
# @breakpoints, breakpoints index. Returned from getBPs
# @sample_list: list of sample name
# @ploidy_list: list of ploidy for all cells
# define the reads ratio as abs(R_x - R_y) / mean(R)
def getDist(matrix, BPs, CN_, ploidy_list):
  initP = ploidy_list
  if CN_ is None:
    new_reads = avg_reads(matrix, BPs)
    CN = get_CN(matrix, initP, BPs)
    CN = np.around(CN)
    BPs_ = BPs
  else:
    # BPs is a 2d list of bp
    new_reads = avg_reads(matrix, BPs)
    CN = CN_
    BPs_ = [item for sublist in BPs for item in sublist]
    BPs_ = sorted(list(set(BPs_)))
   # get read count ratio
  R_list = []

  for i in range(matrix.shape[1]):
    r = []
    for j in (BPs_):
      if j + 1 > len(CN) - 1:
        break
      if CN[j-1,i] == CN[j+1,i]:
        r.append(0)
      else:
        r.append(abs(new_reads[j+1,i]-new_reads[j-1,i]) / np.mean(new_reads[:,i]))
    R_list.append(r)
  R_list = np.array(R_list)
  dist = metrics.pairwise.manhattan_distances(R_list)
  # print(R_list[7])
  # print(R_list[12])
  # for d in dist:
  #   print(d)
  dist_ = (dist - np.amin(dist))/(np.amax(dist) - np.amin(dist) + 10**(-16))
  return dist_, R_list

# this function update the reads counts by gc and map for each cluster
# @reads: read count matrix, bins*cell
# @ploidy_list: list of ploidy for each cell
# @df: dataframe of gc and mappability 3/29		
def update_reads(reads, cluster, df, ploidy, BPs):
	CN = get_CN(reads, ploidy, BPs)
	CN = np.round(CN.T) # now is cell by bins
	reads = reads.T # now is cell by bins
	used_bins = []
	print("Updating reads", reads.shape)
	M_CN = {}
	# get median copy number for each
	for c in range(max(cluster) + 1):
		cells = np.where(cluster == c)[0]
		M_CN[c] = np.round(np.median(CN[cells],axis = 0))
	for b in range(reads.shape[1]):
		# find bins with same gc and map
		if b in used_bins:
			continue
		curr_gc = df.iloc[b,3]
		curr_map = df.iloc[b,4]
		same_gc_map = df[(df['gc'] == curr_gc) & (df['map'] == curr_map)]
		ind = (same_gc_map.index.to_list())
		used_bins = used_bins+ind
		if len(ind) >1:
			for c in range(max(cluster) + 1):
			# for each cluster
			# get cell index in this cluster
				cells = np.where(cluster == c)[0]
				# get median CN for this gc map groups
				CNs = set(M_CN[c][ind])
				for cn in CNs:
					ind_ = np.where(M_CN[c] == cn)[0]
					ind1 = []
					for x in ind_:
						if x in ind:
							ind1.append(x)
					if len(ind1) > 1:
						for cell in cells:
							med_read = np.median(reads[cell][ind1], axis = None)
							reads[cell][ind1] =  med_read
	return reads.T


# search ploidy that minimize the cost function
# reads: updated reads
# BPs: breakpoints
# clusters: list of lables
# ploidy: initial ploidy
def searchP( reads, BPs, initCN, clusters, ploidy):
  print("in search p")
  print(clusters)
  resP = copy.deepcopy(ploidy)
  sim = getSim(reads, BPs, initCN,  ploidy.reshape((1,-1)))
  for c in range(max(clusters) + 1):
    cells = np.where(np.array(clusters) == c)
    #minP = np.min(resP[cells[0]])
    #maxP = np.max(resP[cells[0]])
    minP = np.percentile(resP[cells[0]], 25)
    maxP = np.percentile(resP[cells[0]], 75)
    if maxP - minP <= 0.05:
       X = [maxP]
    else:
      X = np.arange(minP, maxP, 0.05)
    #X = np.arange(1.5, 5, 0.05)
    reads_ = reads.T[cells]
    bestE = np.Infinity
    sim_ = sim[cells[0]][:,cells[0]]
    weight = np.sum(sim_, axis = 1)/len(cells[0])
    bestE = np.Infinity
    bestP = 0
    for p in X:
      CN_nonInt = get_CN(reads_.T, [p], BPs).T # cells * bins
      # work good on both  and CRC2
      E = np.sum(np.abs(np.round(CN_nonInt) - CN_nonInt)*weight.reshape((-1,1))) + np.sum(pairwise_distances(np.round(CN_nonInt))*weight.reshape((-1,1)))/2
      print("current p", p,"E", E)
      if E < bestE:
        bestE = E
        bestP = p
    print("cluster", c, "best P", bestP)
    for cell in cells:
      #resP[cell] = np.median(resP[cells[0]])
      resP[cell] = bestP
  return resP
def scale_gc_map(gc_map_df, bin_list_df):
  #df = pd.DataFrame(bin_list)
  #df = df[0].str.split(r':|-', expand = True)
  #df = df.rename(columns={0: "CHROM", 1: "START", 2:"END"})
  #df['START'] = df['START'].astype(int)
  #df['END'] = df['END'].astype(int)
  df = bin_list_df
  gc = gc_map_df
  df = df.merge(gc)
  gc_min = df['gc'].min()
  gc_max = df['gc'].max()
  map_min = df['map'].min()
  map_max = df['map'].max()
  df['gc'] = round((100/(gc_max - gc_min))*(df['gc'] - gc_max) + 100)
  df['map'] = round((100/(map_max - map_min))*(df['map'] - map_max) + 100)
  df['gc'] = df['gc'].astype(int)
  df['map'] = df['map'].astype(int)
  return df


# find the best clusters
# return the clustering with min cost
# if a cluster has < 2 cells, this clustering will not be considered
# BPs, breakpoints
# cov_matrix, coverage matrix
# sample, sample name list
# 4/15 maximize inter-clusters distance
def getClusters(BPs, initCN, cov_matrix, initP,  maxK = 8):
  # remove outliner from clustering
  outliers = []

  dist, R = getDist(cov_matrix, BPs, initCN, initP)
  q2 = np.percentile(np.array(dist), 50)
  print("q2", q2)
  for i in range(len(initP)):
    dist_ = dist[i]
    dist_ = np.delete(dist_,[i])
    if np.all(dist_ > q2):
      outliers.append(i)
  # get remaining matrix, sample_ and get Ratio matrix etc.
  matrix_ = np.delete(cov_matrix, outliers, axis = 1)
  ploidy_list = initP
  ploidy_list = np.delete(ploidy_list, outliers, axis = 0)

  R = np.delete(R, outliers, axis = 0)
  print(outliers)
  if isinstance(BPs[0], list):
    BPs_ = [sublist for i, sublist in enumerate(BPs) if i not in outliers]

  else:
    BPs_ = BPs
  cells_count = matrix_.shape[1]
  bestE = sys.maxsize 
  bestK = 1
  # when cluster size = 1
  medPloidy = np.median(ploidy_list)
  bestPloidy = [medPloidy] * cells_count
  CN_nonInt = get_CN(matrix_, bestPloidy, BPs_)
  bestE = get_entropy(CN_nonInt, bestPloidy, np.array([0]*len(bestPloidy)), BPs)		
  bestCluster = [0] * cells_count
  print("cluters size 1", bestE)
  for c in range(2, maxK + 1):
    #clustering = AgglomerativeClustering(n_clusters = c, metric="precomputed",linkage="complete" ).fit(dist)
    clustering = AgglomerativeClustering(n_clusters = c, affinity="l1",linkage="complete" ).fit(R)    
    cluster2ploidy = {}
    # get median ploidy of each cluster
    for i in range(c):
      cells = np.where(clustering.labels_ == i)
      medPloidy = np.median(ploidy_list[cells[0]])
      print(ploidy_list[cells[0]])
      print("cluster", i, "median ploidy", medPloidy)
      cluster2ploidy[i] = medPloidy
    print(clustering.n_clusters_, clustering.labels_)
    newPloidy_list = []
    for i in range(len(ploidy_list)):
      cluster = clustering.labels_[i]
      newPloidy_list.append(round(cluster2ploidy[cluster],3))
    CN_nonInt = get_CN(matrix_, newPloidy_list, BPs_)
    E = get_entropy(CN_nonInt, newPloidy_list, clustering.labels_, BPs)
    print("clusters size", c, "E", E)
    freq = Counter(clustering.labels_)
    flag = False
    for value in freq.values():
      if value < 2:
        flag = True
        break
    if E < bestE and not flag:
      bestE = E
      bestPloidy = newPloidy_list
      bestCluster = clustering.labels_
      bestK = c
  return bestK, bestPloidy, bestCluster, outliers

# use SeCNV's segmentation method
def segmentation(path, cov, gc, norm_file, ref, K, s):
  Y, cov_df, chr_name, bin_list_df, sample_list = read_matrix(os.path.join(path, cov))
  gc_map_df = pd.read_csv(os.path.join(path,gc), sep = "\t")
  bin_list_df, gc_map_df, Y, cov_df = filterBins(cov_df, gc_map_df)
  print(cov_df.head())
  chr_name = np.array(cov_df['CHROM'].values.tolist())
  print("Y shape after filtering bin", Y.shape)
  var, norm_cell_index, abnorm_cell_index = get_norm_cell(Y, sample_list, norm_file )
  cov_matrix = Bias_norm(Y, norm_cell_index, ref)
  temp = pd.DataFrame(cov_matrix, columns = sample_list)
  temp['CHROM'] = bin_list_df['CHROM']
  temp['START'] = bin_list_df['START']
  temp['END'] = bin_list_df['END']
  
  temp.to_csv("normalized_reads.tsv", sep = "\t", index = False)
  chosen_chr = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
  breakpoints = []
  bins_len = []
  for chrom in chosen_chr:
    print("process %s..."%chrom)
    indexes = np.where(chr_name == chrom)
    chrom_matrix = copy.deepcopy(cov_matrix[indexes])
    # matrix = chrom_matrix[:, abnorm_cell_index]
    matrix = chrom_matrix
    adj_matrix = make_adj_matrix(matrix, K, s)
    print(adj_matrix.shape)
    print("Segmenting %s ... "%chrom)
    dp_process = DP_process(adj_matrix, 20)
    dp_process.dp_init()
    dp_process.dp_process()
    dp_process.find_k()
    boundaries = dp_process.find_boundaries()
    breakpoints.append(boundaries)
    bins_len.append(adj_matrix.shape[0])
  BPs = getBPs(bins_len, breakpoints)
  print(len(BPs), BPs)
  cov_df = pd.DataFrame(cov_matrix, columns = cov_df.columns.tolist()[3:])
  cov_df = pd.concat([bin_list_df, cov_df], axis = 1)
  return BPs, cov_matrix, cov_df, bin_list_df, sample_list

# filter bins fo other input reads/CNs
# cov: np array, cov from read_matrix
# gc_df: pd df, from scale_gs_map
def filterBins(cov_df, gc_df):
  print("Filtering bins based on gc and mappability")
  df_filtered = gc_df[(gc_df['gc'] > 0.2) & (gc_df['gc'] < 0.8) & (gc_df['map'] > 0.9)]
  print(df_filtered.shape)
  print(cov_df.shape)
  df_filtered = df_filtered.merge(cov_df)
  bin_list_df = df_filtered[['CHROM', 'START', 'END']]
  gc = df_filtered[['CHROM', 'START', 'END', "gc", "map"]]
  cov_ = df_filtered.drop(['CHROM', 'START', 'END', "gc", "map"], axis = 1).to_numpy()
  cov_df = df_filtered.drop(["gc", "map"], axis = 1)
  return bin_list_df, gc, cov_, cov_df

# normalize read counts for other methods
# when no normal cell detected
# has to be 500k. Need to do some filtering as other methods other than secnv are likely to have different bed as the bias. 
def Bias_norm_others(cov_df, gc_map_df, ref):
  print("Normalizting reads using built-in bias")
  if ref == "hg19":
    ave_bias = pd.read_csv(os.path.join(src_path, "hg19_bias.txt"), sep = "\t")
  else:
    ave_bias = pd.read_csv(os.path.join(src_path, "hg38_bias.txt"), sep = "\t")
  # for i in range(ave_bias.shape[0]):
  #   chrom, start, end = ave_bias.loc[i, ["CHROM", "START", "END"]]
  #   if len(cov_df[(cov_df['CHROM'] == chrom) & (cov_df['START'] == start) & (cov_df['END'] == end)]) == 0:
  #     print( chrom, start, end)
  merged_df = ave_bias.merge(cov_df)
  print(merged_df.shape)
  
  merged_df = merged_df.merge(gc_map_df)
  merged_df.dropna(inplace=True)
  if merged_df.shape[0] == 0:
    print("Reads are not in 500kb. Normalization is not performed")
    Y = cov_df.drop(["CHROM", "START", "END"], axis = 1).to_numpy()
    return Y, gc_map_df, gc_map_df[["CHROM", "START", "END"]] 
  cov_df_ = merged_df.drop(["Bias","gc", "map"], axis=1)
  gc_map_df_ = merged_df[["CHROM", "START", "END", "gc", "map"]]
  ave_bias = merged_df['Bias'].to_numpy()
  Y = cov_df_.drop(["CHROM", "START", "END"], axis = 1).to_numpy()
  #Y = Y.dropna().to_numpy() # 5071* #cells
  gc_nor_Y = Y.T / ave_bias # ( #cells, 5071) / (5071,)
  return gc_nor_Y.T, gc_map_df_, gc_map_df_[["CHROM", "START", "END"]]  
  
# find breakpoints given copy number profile from other method
# return breakpoints and cov matrix
def findBPs(path, covF, CNF, gc, ref):
  df = pd.read_csv(os.path.join(path, CNF), sep = "\t")
  cov, cov_df, chr_name, bin_list_df, sample_list = read_matrix(os.path.join(path, covF))
  gc_map_df = pd.read_csv(os.path.join(path,gc), sep = "\t")
  gc_map_df = gc_map_df.merge(bin_list_df)
  print(cov_df.head())
  bin_list_df, gc_map_df, cov, cov_df = filterBins(cov_df, gc_map_df)
  print("shape after filtering", cov.shape)
  var, norm_cell_index, abnorm_cell_index = get_norm_cell(cov, sample_list, "None" )
  if len(norm_cell_index) > 0:
    cov = Bias_norm(cov, norm_cell_index, ref)
  else:
    cov, gc_map_df, bin_list_df = Bias_norm_others(cov_df, gc_map_df, ref)
  gc_map_df = scale_gc_map(gc_map_df, bin_list_df)
  print("shape after normalization", cov.shape)
  BPlist = []
  df = df.merge(gc_map_df, how="inner").drop(['gc', 'map'],axis=1)
  df = df.reset_index(drop=True)
  print(df.head())
  CN = df.iloc[:, 3:]
  # Loop through the columns and find the row indices where values change
  initP = []
  # get breakpoint for each sample
  BPlists = []
  for column in CN.columns:
    change_indices = df[column].ne(df[column].shift()).iloc[1:].index[df[column].ne(df[column].shift()).iloc[1:]].tolist()
    # if selected_columns.shape[0] - 1 not in change_indices:
    #    change_indices.append(selected_columns.shape[0] - 1)
    # if len(BPlists) > 0 and crossSample and change_indices not in BPlists:
    #   crossSample = False
    BPlist.extend(change_indices)
    BPlists.append(change_indices)
  count_dict = Counter(BPlist)
  count_dict = dict(sorted(count_dict.items()))
  BPs = []
  for key, value in count_dict.items():
    if value >= 3:
      if key - 1 in BPs:
        continue
      BPs.append(key)
  print(len(BPs),BPs)  
#   # get common breakpoints for cross sample data
#   print(crossSample)
# # Compare each row with the previous row
#   all_bps = selected_columns.ne(selected_columns.shift()).any(axis=1)
# # Get the row indices where the current row is different from the previous row
#   allBPlist = all_bps[all_bps].index.tolist()
#   allBPlist = allBPlist[1:-1]
  print(CN.shape)
  initP = CN.mean().tolist()
#   allBPlist = sorted(allBPlist)
  return BPlists, cov, CN.to_numpy(), bin_list_df, gc_map_df, sample_list, initP


def main(covfile, CNfile, path, minP, maxP, K, s, gc, outfile, norm_file, ref, secnv = True):
  BPs = None
  cov_matrix = None
  initP = None
  bin_list = None
  sample = None
  BPlists = None
  initCN = None
  gc_map_df = None
  if CNfile == "None":
    BPs, cov_matrix, cov_df, bin_list, sample = segmentation(path, covfile, gc, norm_file, ref, K, s)
    initP =init_ploidy(cov_matrix, BPs, minP, maxP)
    gc_map_df = pd.read_csv(os.path.join(path, gc), sep = "\t")
    gc_map_df = scale_gc_map(gc_map_df, bin_list)
    # return secnv's result as intermediate result
    seCNV_matrix = get_CN(cov_matrix, initP, BPs)
    df = save_matrix(np.round(seCNV_matrix.T), bin_list, sample, os.path.join(path,  "SeCNV_" + outfile + "_cnv.tsv"),
            os.path.join(path,"SeCNV_"+outfile + "_cnv_meta.tsv"))
    print(initP)
  else:
    print("Use breakpoints from other method")
    global builtin
    builtin = False
    BPs, cov_matrix, initCN, bin_list, gc_map_df, sample, initP = findBPs(path, covfile, CNfile, gc, ref)
    print(initP)
    initP = init_ploidy(cov_matrix, BPs, minP, maxP)
  similarity = getSim(cov_matrix, BPs, initCN, initP)
  matrix = cov_matrix
  bestK, bestPloidy, bestCluster, outliers = getClusters(BPs, initCN, cov_matrix, initP,  maxK = 8)
  #print("best number of clusters is ", bestK)
  #print(bestCluster)
  smallMat = np.delete(matrix, outliers, axis = 1)
  smallInitP = np.delete(initP, outliers, axis=0)
  #print(smallInitP)
  bigCluster = copy.deepcopy(bestCluster)
  if isinstance(BPs[0], list):
    BPs_ = [sublist for i, sublist in enumerate(BPs) if i not in outliers]
  else:
    BPs_ = BPs
  #print("before search P")
  smallP = searchP(smallMat, BPs_, initCN, bigCluster, np.array(smallInitP).reshape(-1,1))
  smallP = smallP.flatten()
  bigMat = copy.deepcopy(smallMat)
  bigP = copy.deepcopy(smallP)
  for c in outliers:
    bigMat = np.insert(bigMat, c, matrix[:,c], axis = 1)
    bigP= np.insert(bigP, c, -1)
    bigCluster = np.insert(bigCluster, c, -1)
  for c in outliers:
    indexed_list = list(enumerate(similarity[c]))
    sorted_list = sorted(indexed_list, key=lambda x: x[1], reverse=True)
    for d in sorted_list:
      if d[0] != c and bigCluster[d[0]] != -1:
        bigCluster[c] = bigCluster[d[0]]
        bigP[c] = bigP[d[0]]
        break

  bigMat = update_reads(bigMat, bigCluster, gc_map_df, bigP, BPs)
  seCNV_matrix = get_CN(bigMat, bigP.reshape((1,-1)), BPs)
  df = save_matrix(np.round(seCNV_matrix.T), bin_list, sample, os.path.join(path, outfile+ "_cnv.tsv"),
                  os.path.join(path, outfile + "_cnv_meta.tsv"))

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-cov", help="Coverage matrix", type=str)
  parser.add_argument("-CN", help="Copy number profile if not use built-in segmentaiton method", default="None")
  parser.add_argument("-path", help="work path", type=str)
  parser.add_argument("-minP", help="Miminum ploidy", type=float, default=1.5)
  parser.add_argument("-maxP", help="Maximum ploidy", type=float, default=5)
  parser.add_argument("-K", help = "SeCNV: The K largest distances used to construct adjacency matrix", default="auto_set")
  parser.add_argument("-s", help = "SeCNV: The standard deviation of the Gaussian kernel function", default="auto_set")
  parser.add_argument("-gc", help = "Path to the gc and maapability file", type=str)
  parser.add_argument("-out", help="Output file name")
  parser.add_argument("-norm", help="Normal cell file", default="None")
  parser.add_argument("-ref", help= "reference version hg38 or hg19", default ="hg19")
  args = parser.parse_args()
  main(args.cov, args.CN, args.path, args.minP, args.maxP, args.K, args.s, args.gc, args.out, args.norm, args.ref)
