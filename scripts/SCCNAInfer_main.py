
import os
import sys
import numpy as np
import pandas as pd
from secnv_segmentation import *
from extract_BP import *
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

# version2: take segmentation as input
# assume input segmentation file is chrom, start, end

src_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(src_path)

# read in read count matrix
# @file_name, read file name
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
      bin_list.append([line[0], int(float(line[1])), int(float(line[2]))])
    else:
      sample_list = line
  matrix = np.array(matrix)
  bin_list_df = pd.DataFrame(bin_list, columns = ['CHROM', 'START', 'END'])
  sample_list = np.array(sample_list[3:])
  matrix = matrix.astype(float)
  matrix = matrix.astype(int)
  for i in range(len(sample_list)):
    sample_list[i] = str(sample_list[i]).split("/")[-1]
  chr_name = np.array(chr_name)
  f.close()
  cov_df = pd.read_csv(file_name, sep = "\t")
  # drop last empty column if any
  cov_df = cov_df.dropna(axis=1, how='all')
  cov_df.dropna(inplace=True)
  sample_list = cov_df.columns.to_list()[3:]
  return matrix, cov_df, chr_name, bin_list_df, sample_list

# write out matrix
# @matrix: copy number matrix
# @bin_list_df: df with bin list, chrom, start, end
# @sample_list: list of sample name
# @cnv_name: output cnv file name
# @meta_name: output meta data file name
def save_matrix(matrix, bin_list_df, sample_list, cnv_name, meta_name):
  ploidy = np.average(matrix, axis=1)
  ploidy = np.round(ploidy, 2)
  df = pd.DataFrame(matrix.T, columns=sample_list)
  df = pd.concat([bin_list_df, df], axis=1, join='inner')
  df.to_csv(cnv_name, index=False, sep = "\t")
  df_2 = pd.DataFrame(np.vstack((sample_list, ploidy)).T, columns=["cell_id", "c_ploidy"])
  df_2.to_csv(meta_name, index=False)
  return df


# this function average the read counts to the median for each segment
# @Y: read count matrix
# @BPs: break point list
def avg_reads(Y, BPs):
  Y_cp = copy.deepcopy(Y)
  for j in range(len(BPs)):
    if j == 0:
      Y_cp[0:BPs[j]+1, ] = np.median(Y_cp[0:BPs[j]+1, ], axis = 0)
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
  res = reads*ploidy/np.mean(reads, axis=0)
  res[res > 10] = 10
  return res


# get the cost for clustering
# @CN_nonInt: non interger copy number matrix
# @cluster: cluster list
# @return: return the cost
def get_entropy(CN_nonInt,   cluster):
  k = max(cluster) + 1
  CN = np.around(CN_nonInt).T
  absCNsum = np.sum(np.abs(np.around(CN_nonInt) - CN_nonInt)) 
  #print(cluster)
  intraCluster = 0
  pseudo = []
  interCluster = 0

  for c in range(max(cluster) + 1):
    cells = np.where(cluster == c)
    # get median cp number
    M_CN = np.median(CN[cells], axis = 0)
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

# initialize the ploidy
# @Y: read counts
# @BPs: break points list
# @minP: minimum ploidy
# @maxP: maximum ploidy
# @return: a list of ploidy 
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
    initP.append(ploidy)
  return np.array(initP)

# get breakpoints from all chromosomes, merged into one list
# @bins_len, a list contains the number of bins in each chromosome
# @breakpoints, a list of breakpoints index for each chromosome
# @return a list of breakpoint index with respect to all bins across all chromosome
def mergeBPs(bins_len, breakpoints):
  bps = breakpoints[0]
  # merge breakpoints into a single list and change index
  total_len = bins_len[0]
  for i in range(1, len(breakpoints)):
    bp = np.array(breakpoints[i]) + total_len
    bps = bps + bp.tolist()
    total_len += bins_len[i]
  return bps

# get the similarity matrix between cells
# @matrix: read count matrix
# @ploidy: list ploidy list
# @return: return similarity
def getSim(matrix, BPs, ploidy_list):
  dist, R = getDist(matrix, BPs, ploidy_list)
  return 1 - dist

# get the distance between each pair of cells
# define the reads ratio as abs(R_x - R_y) / mean(R)
# @matrix, read count matrix. bins*cells
# @breakpoints, breakpoints index
# @sample_list: list of sample name
# @ploidy_list: list of ploidy for all cells
# @return: return the distant matrix, and read count ratio list 
def getDist(matrix, BPs, ploidy_list):
  initP = ploidy_list
  new_reads = avg_reads(matrix, BPs)
  CN = get_CN(matrix, initP, BPs)
  CN = np.around(CN)
  BPs_ = BPs
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
  dist_ = (dist - np.amin(dist))/(np.amax(dist) - np.amin(dist) + 10**(-16))
  return dist_, R_list

# this function update the reads counts by gc and map for each cluster
# @reads: read count matrix, bins*cell
# @clusters: cluster list
# @ploidy: list of ploidy for each cell
# @gc_map_df: dataframe of gc and mappability 
# @BPs: breakpoint list	
# @return update reads
def update_reads(reads, cluster, gc_map_df, ploidy, BPs):
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
		curr_gc = gc_map_df.iloc[b,3]
		curr_map = gc_map_df.iloc[b,4]
		same_gc_map = gc_map_df[(gc_map_df['gc'] == curr_gc) & (gc_map_df['map'] == curr_map)]
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
# @reads: updated reads
# @BPs: breakpoints
# @clusters: list of lables
# @ploidy: initial ploidy
# @return: return  a list of optimal copy number
def searchP( reads, BPs,  clusters, ploidy):
  print("in search p")
  print(clusters)
  resP = copy.deepcopy(ploidy)
  sim = getSim(reads, BPs,  ploidy.reshape((1,-1)))
  for c in range(max(clusters) + 1):
    cells = np.where(np.array(clusters) == c)
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
      E = np.sum(np.abs(np.round(CN_nonInt) - CN_nonInt)*weight.reshape((-1,1))) + np.sum(pairwise_distances(np.round(CN_nonInt))*weight.reshape((-1,1)))/2
      #print("current p", p,"E", E)
      if E < bestE:
        bestE = E
        bestP = p
    print("cluster", c, "best P", bestP)
    for cell in cells:
      #resP[cell] = np.median(resP[cells[0]])
      resP[cell] = bestP
  return resP

# scale gc and mapp to integer from 0 - 100
# @gc_map_df: gc and mapp df
# @bin_list_df: a df of bin list
def scale_gc_map(gc_map_df, bin_list_df):
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
# if a cluster has < 2 cells, this clustering will not be considered
# @BPs, breakpoints
# @cov_matrix, coverage matrix
# @sample, sample name list
# @initP: initial ploidy
# @maxK: maximum number of cluster being considered
# @return the clustering with min cost
def getClusters(BPs,  cov_matrix, initP, samples, perc = 50, maxK = 8):
  # remove outliner from clustering
  outliers = []
  dist, R = getDist(cov_matrix, BPs,  initP)
  q2 = np.percentile(np.array(dist), perc)
  print("q2", q2)
  sim = 1 - dist
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
  dist_ = pd.DataFrame(dist)
  dist_.drop(outliers, inplace = True, axis = 1)
  dist_ = dist_.values
  dist_ = np.delete(dist_, outliers, axis = 0)
  BPs_ = BPs
  cells_count = matrix_.shape[1]
  bestE = sys.maxsize 
  bestK = 1
  # when cluster size = 1
  medPloidy = np.median(ploidy_list)
  bestPloidy = [medPloidy] * cells_count
  CN_nonInt = get_CN(matrix_, bestPloidy, BPs_)
  bestE = get_entropy(CN_nonInt, np.array([0]*len(bestPloidy)))		
  bestCluster = [0] * cells_count
  print("cluters size 1", bestE)
  for c in range(2, maxK + 1):
    clustering = AgglomerativeClustering(n_clusters = c, metric="l1",linkage="complete" ).fit(R)    
    #clustering = AgglomerativeClustering(n_clusters = c, metric="precomputed",linkage="complete" ).fit(1- dist_)  
    cluster2ploidy = {}
    # get median ploidy of each cluster
    for i in range(c):
      cells = np.where(clustering.labels_ == i)
      medPloidy = np.median(ploidy_list[cells[0]])
      print(ploidy_list[cells[0]])
      #print(cells[0])
      #print(np.array(samples)[cells[0]])
      print("cluster", i, "median ploidy", medPloidy)
      cluster2ploidy[i] = medPloidy
    print(clustering.n_clusters_, clustering.labels_)
    newPloidy_list = []
    for i in range(len(ploidy_list)):
      cluster = clustering.labels_[i]
      newPloidy_list.append(round(cluster2ploidy[cluster],3))
    CN_nonInt = get_CN(matrix_, newPloidy_list, BPs_)
    E = get_entropy(CN_nonInt, clustering.labels_)
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
# @path: working path
# @cov: coverage matrix
# @gcF: gc and mapp file
# @norm_file: file with nomral cells
# @ref: reference type, hg19, hg38 
# @K: top K to consider when segmenting
# @return: BPs, cov_matrix, cov_df, bin_list_df, sample_list
def segmentation(path, cov, gcF, norm_file, ref, K, s):
  Y, cov_df, chr_name, bin_list_df, sample_list = read_matrix(os.path.join(path, cov))
  gc_map_df = pd.read_csv(os.path.join(path,gcF), sep = "\t")
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
  BPs = mergeBPs(bins_len, breakpoints)
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
  df_filtered = df_filtered.reset_index(drop=True)
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
# @df: cnv df
# @gc_map_df: df for gc and mapp 
# return breakpoints 
def findBPs(df, gc_map_df):
  df = df.merge(gc_map_df, how="inner").drop(['gc', 'map'],axis=1)
  df = df.reset_index(drop=True)
  BPlist = []
  CN = df.iloc[:, 3:]
  # Loop through the columns and find the row indices where values change
  for column in CN.columns:
    change_indices = df[column].ne(df[column].shift()).iloc[1:].index[df[column].ne(df[column].shift()).iloc[1:]].tolist()
    BPlist.extend(change_indices)
  count_dict = Counter(BPlist)
  count_dict = dict(sorted(count_dict.items()))
  BPs = []
  for key, value in count_dict.items():
    if value >= 3:
      if key - 1 in BPs and value + 3< count_dict[key - 1]:
        continue
      BPs.append(key)
  BPs = sorted(BPs)
  return BPs

# get BPs from CNV matrix or directly from BP files, filter bins, normalize reads
# @path: working path
# @cov: coverage matrix
# @gcF: gc and mapp file
# @norm_file: file with nomral cells
# @ref: reference type, hg19, hg38 
# return breakpoints and cov matrix
def getBPs(path, covF, CNF, gc, ref):
  df = pd.read_csv(os.path.join(path, CNF), sep = "\t")
  cov, cov_df, chr_name, bin_list_df, sample_list = read_matrix(os.path.join(path, covF))
  gc_map_df = pd.read_csv(os.path.join(path,gc), sep = "\t")
  gc_map_df = gc_map_df.merge(bin_list_df)
  # filter bins
  bin_list_df, gc_map_df, cov, cov_df = filterBins(cov_df, gc_map_df)
  # normalize reads
  var, norm_cell_index, abnorm_cell_index = get_norm_cell(cov, sample_list, "None" )
  if len(norm_cell_index) > 0:
    cov = Bias_norm(cov, norm_cell_index, ref)
  else:
    cov, gc_map_df, bin_list_df = Bias_norm_others(cov_df, gc_map_df, ref)
  gc_map_df = scale_gc_map(gc_map_df, bin_list_df)
  print("shape after normalization", cov.shape)
  BPlist = []
  # only contains chrom, start, end
  # input is a breakpoint matrix
  if df.shape[1] == 3: 
    # sort and deduplicates the bp df
    df = df.sort_values(by=['CHROM', 'START', 'END'])
    df.drop_duplicates(inplace=True)
    df = df.reset_index(drop=True)
    # add the end of the chromosome to the BPlist
    for i in range(1, 23):
      temp = bin_list_df[(bin_list_df['CHROM'] == 'chr' + str(i))]
      BPlist.append(temp.index.tolist()[-1])
    for idx, row in df.iterrows():
      temp = bin_list_df[(bin_list_df['CHROM'] == row["CHROM"])& (bin_list_df['END'] == row["END"])]
      if len(temp.index.tolist()) == 0:
        continue
      i = temp.index.tolist()[0]
      BPlist.append(i)
      BPlist = set(BPlist)
      BPlist = sorted(BPlist)
  else:
    BPlist = findBPs(df, gc_map_df)
  return BPlist, cov, gc_map_df, bin_list_df, sample_list

# this function find the cloest cluster
def findCloestCluster(clusters, ploidy, outliers, similarity):
  for o in outliers:
    sim = similarity[o]
    max_sim = 0
    p = 0
    g = -1
    for c in range(max(clusters) + 1):
      cells = np.where(clusters == c)[0]
      avg_sim = np.mean(sim[cells])
      print("outlier", o, "cluster", c, "sim", avg_sim, ploidy[cells[0]])
      if avg_sim > max_sim:
        max_sim = avg_sim
        p = ploidy[cells[0]]
        g = c
    ploidy[o] = p
    clusters[o] = g
  return ploidy, clusters
  
   
def main(covfile, CNfile, path, minP, maxP, K, s, gc, outfile, norm_file, ref, perc):
  BPs = None
  cov_matrix = None
  initP = None
  bin_list = None
  sample = None
  gc_map_df = None
  if CNfile == "None":
    BPs, cov_matrix, cov_df, bin_list, sample = segmentation(path, covfile, gc, norm_file, ref, K, s)
    gc_map_df = pd.read_csv(os.path.join(path, gc), sep = "\t")
    gc_map_df = scale_gc_map(gc_map_df, bin_list)
    # return secnv's result as intermediate result
    initP =init_ploidy(cov_matrix, BPs, minP, maxP) 
    seCNV_matrix = get_CN(cov_matrix, initP, BPs)
    df = save_matrix(np.round(seCNV_matrix.T), bin_list, sample, os.path.join(path,  "SeCNV_" + outfile + "_cnv.tsv"),
            os.path.join(path,"SeCNV_"+outfile + "_cnv_meta.tsv"))  
  else:
    print("Use breakpoints from other method")
    BPs, cov_matrix, gc_map_df, bin_list, sample = getBPs(path, covfile, CNfile, gc, ref)

  initP =init_ploidy(cov_matrix, BPs, minP, maxP) 
  similarity = getSim(cov_matrix, BPs, initP)
  matrix = cov_matrix
  bestK, bestPloidy, bestCluster, outliers = getClusters(BPs, cov_matrix, initP,  sample, perc, maxK = 8)
  print("best K is ", bestK)
  print(bestCluster)
  smallMat = np.delete(matrix, outliers, axis = 1)
  smallInitP = np.delete(initP, outliers, axis=0)
  print(smallInitP)
  bigCluster = copy.deepcopy(bestCluster)

  BPs_ = BPs
  print("before search P")
  smallP = searchP(smallMat, BPs_,  bigCluster, np.array(smallInitP).reshape(-1,1))
  smallP = smallP.flatten()
  bigMat = copy.deepcopy(smallMat)
  bigP = copy.deepcopy(smallP)
  for c in outliers:
    bigMat = np.insert(bigMat, c, matrix[:,c], axis = 1)
    bigP= np.insert(bigP, c, -1)
    bigCluster = np.insert(bigCluster, c, -1)
  bigP, bigCluster = findCloestCluster(bigCluster, bigP, outliers, similarity)
  for c in outliers:
    print(bigP[c], bigCluster[c])
  # for c in outliers:
  #   indexed_list = list(enumerate(similarity[c]))
  #   sorted_list = sorted(indexed_list, key=lambda x: x[1], reverse=True)
  #   for d in sorted_list:
  #     if d[0] != c and bigCluster[d[0]] != -1:
  #       bigCluster[c] = bigCluster[d[0]]
  #       bigP[c] = bigP[d[0]]
  #       break

  bigMat = update_reads(bigMat, bigCluster, gc_map_df, bigP, BPs)
  seCNV_matrix = get_CN(bigMat, bigP.reshape((1,-1)), BPs)
  df = save_matrix(np.round(seCNV_matrix.T), bin_list, sample, os.path.join(path, outfile+ "_cnv.tsv"),
                  os.path.join(path, outfile + "_cnv_meta.tsv"))

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-cov", help="Coverage matrix", type=str)
  parser.add_argument("-CN", help="Copy number profile or segmetation file if not use built-in segmentaiton method", default="None")
  parser.add_argument("-path", help="work path", type=str)
  parser.add_argument("-minP", help="Miminum ploidy", type=float, default=1.5)
  parser.add_argument("-maxP", help="Maximum ploidy", type=float, default=5)
  parser.add_argument("-K", help = "SeCNV: The K largest distances used to construct adjacency matrix", default="auto_set")
  parser.add_argument("-s", help = "SeCNV: The standard deviation of the Gaussian kernel function", default="auto_set")
  parser.add_argument("-gc", help = "Path to the gc and maapability file", type=str)
  parser.add_argument("-out", help="Output file name")
  parser.add_argument("-norm", help="Normal cell file", default="None")
  parser.add_argument("-ref", help= "reference version hg38 or hg19", default ="hg19")
  parser.add_argument("-perc", help= "percentile to use to remove outlier", type = int, default =50)
  args = parser.parse_args()
  main(args.cov, args.CN, args.path, args.minP, args.maxP, args.K, args.s, args.gc, args.out, args.norm, args.ref, args.perc)
