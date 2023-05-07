
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

def read_matrix(file_name):
  f = open(file_name)
  matrix = []
  chr_name = []
  bin_list = []
  count = 0
  for line in f:
    line = line.strip("\n").rstrip("\t").split("\t")
    count += 1
    if not count == 1:
      matrix.append(line[3:])
      chr_name.append(line[0])
      bin_list.append(line[0]+":"+line[1]+"-"+line[2])
    else:
      sample_list = line
  matrix = np.array(matrix)
  bin_list = np.array(bin_list)
  sample_list = np.array(sample_list[3:])
  matrix = matrix.astype(int)
  for i in range(len(sample_list)):
    sample_list[i] = str(sample_list[i]).split("/")[-1]
  chr_name = np.array(chr_name)
  f.close()
  return matrix, chr_name, bin_list, sample_list

def save_matrix(matrix, bin_list, sample_list, cnv_name, meta_name):
  ploidy = np.average(matrix, axis=1)
  ploidy = np.round(ploidy, 2)
  df = pd.DataFrame(matrix.T, columns=sample_list, index=bin_list)
  ind = df.index.to_series().str.split(r':|-', expand = True)
  ind = ind.rename(columns={0: "CHROM", 1: "START", 2:"END"})
  df = pd.concat([ind, df], axis=1, join='inner')
  df.to_csv(cnv_name, index=False, sep = "\t")
  df_2 = pd.DataFrame(np.vstack((sample_list, ploidy)).T, columns=["cell_id", "c_ploidy"])
  df_2.to_csv(meta_name, index=False)
  return df


# this function average the read counts to the median for each segment
def avg_reads(Y, BPs):
  Y_cp = copy.deepcopy(Y)
  for j in range(len(BPs)):
    if j == 0:
      Y_cp[0:BPs[j]+1,] = np.median(Y_cp[0:BPs[j]+1,], axis = 0)
    else:
      Y_cp[BPs[j-1]+1:BPs[j]+1, ] = np.median(Y_cp[BPs[j-1]+1:BPs[j]+1, ], axis = 0)
  return Y_cp
  
# get non interger copy number
# @Y: read count matrix
# @ploidy: list of ploidy
# @breakpoints: list of breakpoint for finding median read counts for each seg
def get_CN(Y, ploidy, BPs):
  # Y: row: bins, col: cells
  # get median read count for each segments
  reads = avg_reads(Y, BPs)
  #print("get_CN avg reads", np.array_equal(reads[:,0], reads[:,91]), np.array_equal(reads[:,0], reads[:,98]))
  res = reads*ploidy/np.mean(reads, axis=0)
  res[res > 10] = 10
  #print("get_CN CN", np.array_equal(np.round(res[:,0]), np.round(res[:,91])), np.array_equal(np.round(res[:,0]), np.round(res[:,98])))
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
  sim = getSim(pseudo.T, BPs, np.array(M_P))
  distSum = 0
  for c in range(max(cluster) + 1):
    cells = np.where(cluster == c)
    interCluster[c][c] = np.Inf
    #print(interCluster[c])
    c1 = interCluster[c].argmin()

    distSum += np.sum(np.abs(pseudo[c1] - pseudo[c]))
  distSum = distSum
  return np.log(  absCNsum*(1+0.05*(k))  - distSum +  intraCluster*(1+0.05*k))

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

# get read count ratio
def getR(matrix, BPs, ploidy_list):
  initP = ploidy_list
  new_reads = avg_reads(matrix, BPs)
  Y = get_CN(matrix, initP, BPs)
  Y = np.around(Y)
  # get read count ratio
  R_list = []
  for i in range(matrix.shape[1]):
    r = []
    for j in (BPs):
      if j + 1 > len(Y) - 1:
        break
      if Y[j-1,i] == Y[j+1,i]:
        r.append(0)
      else:
        r.append(abs(new_reads[j+1,i]-new_reads[j-1,i]) / np.mean(new_reads[:,i]))
    R_list.append(r)
  R_list = np.array(R_list)
  return R_list

# matrix: seg * cell
def getSim(matrix, BPs,  ploidy_list):
  dist = getDist(matrix, BPs,   ploidy_list)
  return 1 - dist

# get the distance between each pair of cells
# @matrix, read count matrix. bins*cells
# @breakpoints, breakpoints index. Returned from getBPs
# @sample_list: list of sample name
# @ploidy_list: list of ploidy for all cells
# define the reads ratio as abs(R_x - R_y) / mean(R)
def getDist(matrix, BPs,   ploidy_list):
  initP = ploidy_list
  new_reads = avg_reads(matrix, BPs)
  Y = get_CN(matrix, initP, BPs)
  Y = np.around(Y)
   # get read count ratio
  R_list = []
  for i in range(matrix.shape[1]):
    r = []
    for j in (BPs):
      if j + 1 > len(Y) - 1:
        break
      if Y[j-1,i] == Y[j+1,i]:
        r.append(0)
      else:
        r.append(abs(new_reads[j+1,i]-new_reads[j-1,i]) / np.mean(new_reads[:,i]))
    R_list.append(r)
  R_list = np.array(R_list)
  dist = metrics.pairwise.manhattan_distances(R_list)
  dist_ = (dist - np.amin(dist))/(np.amax(dist) - np.amin(dist))
  return dist_

# this function update the reads counts by gc and map for each cluster
# @reads: read count matrix, bins*cell
# @ploidy_list: list of ploidy for each cell
# @df: dataframe of gc and mappability 3/29		
def update_reads(reads, cluster, df, ploidy, BPs):
	CN = get_CN(reads, ploidy, BPs)
	CN = np.round(CN.T) # now is cell by bins
	reads = reads.T # now is cell by bins
	used_bins = []
	print(reads.shape)
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
def searchP( reads, BPs, clusters, ploidy):
	resP = copy.deepcopy(ploidy)
	sim = getSim(reads, BPs, ploidy.reshape((1,-1)))
	for c in range(max(clusters) + 1):
		cells = np.where(clusters == c)
		minP = np.min(resP[cells[0]])
		maxP = np.max(resP[cells[0]])
		X = np.arange(minP, maxP, 0.05)
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
			#print("current p", p,"E", E)
			if E < bestE:
				bestE = E
				bestP = p
		print("cluster", c, "best P", bestP)
		for cell in cells:
			resP[cell] = bestP
	return resP
def scale_gc_map(gc_map, bin_list):
  df = pd.DataFrame(bin_list)
  df = df[0].str.split(r':|-', expand = True)
  df = df.rename(columns={0: "CHROM", 1: "START", 2:"END"})
  df['START'] = df['START'].astype(int)
  df['END'] = df['END'].astype(int)
  gc = pd.read_csv(gc_map, sep = "\t")
  df = df.merge(gc)
  gc_min = df['gc'].min()
  gc_max = df['gc'].max()
  
  map_min = df['map'].min()
  map_max = df['map'].max()
  df['gc'] = round((100/(gc_max - gc_min))*(df['gc'] - gc_max) + 100)
  df['map'] = round((100/(map_max - map_min))*(df['map'] - map_max) + 100)
  return df


# find the best clusters
# return the clustering with min cost
# if a cluster has < 2 cells, this clustering will not be considered
# BPs, breakpoints
# cov_matrix, coverage matrix
# sample, sample name list
# 4/15 maximize inter-clusters distance
def getClusters(BPs, cov_matrix, initP,  maxK = 8):
  # remove outliner from clustering
  outliers = []
  dist = getDist(cov_matrix, BPs, initP)
  for i in range(len(initP)):
    dist_ = dist[i]
    dist_ = np.delete(dist_,[i])
    if np.all(dist_ > 0.5):
      outliers.append(i)
  # get remaining matrix, sample_ and get Ratio matrix etc.
  matrix_ = np.delete(cov_matrix, outliers, axis = 1)
  ploidy_list = initP
  ploidy_list = np.delete(ploidy_list, outliers, axis = 0)
  cells_count = matrix_.shape[1]
  R = getR(matrix_, BPs,  ploidy_list)
  print(R.shape)
  bestE = sys.maxsize 
  bestK = 0
  # when cluster size = 1
  medPloidy = np.median(ploidy_list)
  bestPloidy = [medPloidy] * cells_count
  CN_nonInt = get_CN(matrix_, bestPloidy, BPs)
  bestE = get_entropy(CN_nonInt, bestPloidy, np.array([0]*len(bestPloidy)), BPs)		
  bestCluster = [0] * cells_count
  print("cluters size 1", bestE)
  for c in range(2, maxK):
    clustering = AgglomerativeClustering(n_clusters = c, metric="l1",linkage="complete" ).fit(R)
    cluster2ploidy = {}
    # get median ploidy of each cluster
    for i in range(c):
      cells = np.where(clustering.labels_ == i)
      medPloidy = np.median(ploidy_list[cells])
      #print("cluster", i, "cells", sample_[cells], "median ploidy", medPloidy)
      cluster2ploidy[i] = medPloidy
    #print(clustering.n_clusters_, clustering.labels_)
    newPloidy_list = []
    for i in range(len(ploidy_list)):
      cluster = clustering.labels_[i]
      newPloidy_list.append(round(cluster2ploidy[cluster],3))
    CN_nonInt = get_CN(matrix_, newPloidy_list, BPs)
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

def main(cov, p, minP, maxP, K, s, gc, outfile, norm_file):
  Y, chr_name, bin_list, sample_list = read_matrix(os.path.join(p, cov))
  var, norm_cell_index, abnorm_cell_index = get_norm_cell(Y, sample_list,norm_file )
  gc_map_df = scale_gc_map(gc, bin_list)
  cov_matrix = Bias_norm(Y, norm_cell_index)

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

  matrix= copy.deepcopy(cov_matrix)
  sample = copy.deepcopy(sample_list)
  BPs = getBPs(bins_len, breakpoints)

  initP =init_ploidy(matrix, BPs, minP, maxP)

  similarity = getSim(matrix, BPs, initP)

  bestK, bestPloidy, bestCluster, outliers = getClusters(BPs, matrix, initP, maxK = 15)
  smallMat = np.delete(matrix, outliers, axis = 1)
  smallInitP = np.delete(initP, outliers, axis=0)
  bigCluster = copy.deepcopy(bestCluster)
  smallP = searchP(smallMat, BPs, bigCluster, np.array(smallInitP).reshape(-1,1))
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
  df = save_matrix(np.round(seCNV_matrix.T), bin_list, sample,   outfile+ "_cnv.csv",  outfile + "_cnv_meta.csv")

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-cov", help="Coverage matrix", type=str)
  parser.add_argument("-p", help="work path", type=str)
  parser.add_argument("-minP", help="Miminum ploidy", type=float, default=1.5)
  parser.add_argument("-maxP", help="Maximum ploidy", type=float, default=5)
  parser.add_argument("-K", help = "SeCNV: The K largest distances used to construct adjacency matrix", default="auto_set")
  parser.add_argument("-s", help = "SeCNV: The standard deviation of the Gaussian kernel function", default="auto_set")
  parser.add_argument("-gc", help = "Path to the gc and maapability file", type=str)
  parser.add_argument("-out", help="Output file name")
  parser.add_argument("-norm", help="Normal cell file", default="None")
  args = parser.parse_args()
  main(args.cov, args.p, args.minP, args.maxP, args.K, args.s, args.gc, args.out, args.norm)
