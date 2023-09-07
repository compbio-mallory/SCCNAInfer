#! /usr/bin/env python

# adapted from SeCNV

import numpy as np
from scipy import stats
import os 
import sys
src_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(src_path)
import KR_norm_juicer
from sklearn.neighbors import KernelDensity
from scipy import signal, stats
import numpy as np
import pandas as pd

def Gaussian_kernel(data,mu,sig):
  temp = np.array(-(data - mu) ** 2 / (2 * sig ** 2), dtype="float")
  density = np.exp(temp)
  return density

# this function get the sd of coverage and return idx for normal and abnormal cells
def get_norm_cell(cov_data, sample_list, normal_cell_file):
  cov_std = stats.variation(cov_data, axis=0)
  if normal_cell_file == "None":
    kde = KernelDensity(kernel="gaussian", bandwidth=0.01).fit(cov_std.reshape(-1, 1))
    dens = np.exp(kde.score_samples(np.arange(0, 1, 0.001).reshape(-1, 1)))
    peaks_pos = signal.argrelextrema(dens, np.greater)[0]
    first_peak = np.arange(0, 1, 0.001)[peaks_pos][0]
    print(first_peak)
    if first_peak < 0.2:
      density_list = Gaussian_kernel(cov_std, first_peak, 0.01)
      norm_index = []
      abnorm_index = []
      for i in range(len(cov_std)):
        if density_list[i] > 1e-3:
          norm_index.append(i)
        else:
          abnorm_index.append(i)
    else:
      print("Warning: No normal cell detected! Use reference to correct bias. Please make sure the bin size is 500 kb.")
      norm_index = []
      abnorm_index = []
      for i in range(len(cov_std)):
        abnorm_index.append(i)
  else:
    print("Use the cells in %s as the normal cells."%normal_cell_file)
    normal_cell_id = np.loadtxt(normal_cell_file, dtype=str)
    norm_index = []
    abnorm_index = []
    for i in range(len(sample_list)):
      if sample_list[i] in normal_cell_id:
        norm_index.append(i)
      else:
        abnorm_index.append(i)

  return cov_std, norm_index, abnorm_index

def make_adj_matrix(matrix, K, sigma):
  if K == "auto_set":
    K = 5
  else:
    K = int(K)
  if sigma == "auto_set":
    sigma = np.median(matrix)
  else:
    sigma = float(sigma)

  adj_matrix = np.zeros((matrix.shape[0], matrix.shape[0]))
  for i in range(matrix.shape[0]):
    for j in range(i+1, matrix.shape[0]):
      temp = (matrix[i]-matrix[j]) ** 2
      temp.sort()
      temp = temp[-K:]
      adj_matrix[i][j] = np.exp(-np.sum(temp)/sigma**2)
      adj_matrix[j][i] = np.exp(-np.sum(temp)/sigma**2)
  print(adj_matrix)
  print(np.sum(adj_matrix))
  adj_matrix = KR_norm_juicer.KR_norm(adj_matrix)

  return adj_matrix
class DP_process:
  def __init__(self,input_matrix,K_max):
    self.matrix = input_matrix
    self.n = len(input_matrix)
    self.max_k = K_max
    self.dp_table = np.zeros((self.n, self.max_k))
    self.index_table = np.zeros((self.n, self.max_k))
    self.pre_table = np.zeros((self.n, self.n))
    for start in range(0, self.n):
      for end in range(start, self.n):
        current_matrix = self.matrix[start:end + 1, start:end + 1]
        self.pre_table[start][end] = 2 * np.sum(np.triu(current_matrix,1))
        self.pre_table[end][start] = np.sum(self.matrix[0:start, start:end + 1]) \
                       + np.sum(self.matrix[end + 1:self.n, start:end + 1])
    self.sum = self.pre_table[0][self.n - 1]
    self.dlogd_sum = np.zeros(self.n)
    self.dlogd_sum[0] = self.pre_table[0][0] * np.log2(self.pre_table[0][0])
    for bin in range(1, self.n):
      self.dlogd_sum[bin] = self.dlogd_sum[bin - 1] + self.pre_table[bin][bin] * np.log2(self.pre_table[bin][bin]+1e-7)

  def calculate_se(self, g, V_p, V):
    if (V == 0):
      return 0
    else:
      return g / self.sum * np.log2(V_p / V)

  #record the initial state of dynamic programming (K = 0)
  def dp_init(self):
    for i in range(0, self.n):
      if i == 0:
        current_volume = self.pre_table[0][0]
      else:
        current_volume = self.pre_table[0][i] + self.pre_table[i][0]
      intra_se = (current_volume * np.log2(current_volume) - self.dlogd_sum[i])/self.sum
      inter_se = self.calculate_se(self.pre_table[i][0], self.sum, current_volume)
      self.dp_table[i][0] = inter_se + intra_se

  #dynamic process
  def dp_process(self):
    for temp_k in range(1, self.max_k):
      for temp_n in range(0, self.n):
        min_se = np.inf
        min_index = 0
        for i in range(0,temp_n):
          se_1 =  self.dp_table[i][temp_k - 1]  # T(i,k-1)
          if i + 1 == temp_n:
            current_volume = self.pre_table[temp_n][temp_n]
          else:
            current_volume = self.pre_table[i+1][temp_n] + self.pre_table[temp_n][i+1]
          intra_se = (current_volume * np.log2(current_volume) - (self.dlogd_sum[temp_n]-self.dlogd_sum[i])) / self.sum
          inter_se = self.calculate_se(self.pre_table[temp_n][i + 1], self.sum, current_volume)
          se_2 = intra_se + inter_se  # H(i+1;n)
          temp_se = se_1 + se_2
          if temp_se < min_se:
            min_se = temp_se
            min_index = i
        self.dp_table[temp_n][temp_k] = min_se
        self.index_table[temp_n][temp_k] = min_index

  #find the best K
  def find_k(self):
    k_list = self.dp_table[-1].tolist()
    self.k = k_list.index(np.min(k_list)) + 1
    self.min_entropy  = np.min(k_list)


  #find the boundaries:
  def find_boundaries(self):
    self.boundaries = []
    self.boundaries.append(int(self.n - 1))
    for i in range(1, self.k):
      self.boundaries.append(int(self.index_table[self.boundaries[i-1]][self.k-i]))
    self.boundaries.reverse()
    return self.boundaries

def Bias_norm(Y, norm_cell_index, ref):
  Y = Y.T
  # clip the extreme reads
  # Calculate the mean and standard deviation of each row
  #cell_means = np.mean(Y, axis=1)
  #cell_stdevs = np.std(Y, axis=1)
  # Calculate the threshold value for each row
  #thresholds = cell_means + 2 * cell_stdevs
  # Replace values greater than the threshold with the threshold value
  for i in range(Y.shape[0]):
    #Y[i][Y[i] > thresholds[i]] = thresholds[i]
    Y[i][Y[i] < 10] = 10
    # use the normal cells to normalize the data
  if len(norm_cell_index) > 0:
    norm_cell_Y = Y[norm_cell_index]
    bias_matrix = []
    for cell in norm_cell_Y:
      bias_list = []
      median = np.median(cell)
      for bin in cell:
        bias = bin/median
        bias_list.append(bias)
      bias_list = np.array(bias_list)
      bias_matrix.append(bias_list)
    bias_matrix = np.array(bias_matrix)
    ave_bias = bias_matrix.mean(axis=0)
    ave_bias = np.where(ave_bias==0, 1, ave_bias)
    temp = pd.DataFrame(ave_bias, columns = ["bias"])
    temp.to_csv("inferred_bias.tsv", sep = "\t", index = False, header = False)  
  else:
    if ref == "hg19":
      ave_bias = pd.read_csv(os.path.join(src_path, "hg19_bias.txt"), sep = "\t")
      ave_bias = ave_bias['Bias'].values
    else:
      ave_bias = pd.read_csv(os.path.join(src_path, "hg38_bias.txt"), sep = "\t")
      ave_bias = ave_bias['Bias'].values
  gc_nor_Y = Y / ave_bias
  return gc_nor_Y.T

def Bias_norm_secnv(Y, norm_cell_index, ref):
  Y = Y.T

    # use the normal cells to normalize the data
  if len(norm_cell_index) > 0:
    norm_cell_Y = Y[norm_cell_index]
    bias_matrix = []
    for cell in norm_cell_Y:
      bias_list = []
      median = np.median(cell)
      for bin in cell:
        bias = bin/median
        bias_list.append(bias)
      bias_list = np.array(bias_list)
      bias_matrix.append(bias_list)
    bias_matrix = np.array(bias_matrix)
    ave_bias = bias_matrix.mean(axis=0)
    ave_bias = np.where(ave_bias==0, 1, ave_bias)
  else:
    if ref == "hg19":
        ave_bias = pd.read_csv(os.path.join(src_path, "hg19_bias.txt"), sep = "\t")
        ave_bias = ave_bias['Bias'].values
    else:
        ave_bias = pd.read_csv(os.path.join(src_path, "hg38_bias.txt"), sep = "\t")
        ave_bias = ave_bias['Bias'].values
  gc_nor_Y = Y / ave_bias
  return gc_nor_Y.T
