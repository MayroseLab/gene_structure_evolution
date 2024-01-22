"""
Given a GFF3 file containing intron features,
calculate various stats regarding the introns
and print to a stats file
"""

import sys
import os
import numpy as np
import pandas as pd
import scipy


### FUNCTIONS

def attributes_to_dict(attr):
  return {a.split('=')[0]: a.split('=')[1] for a in attr.split(';')}

def attributes_to_column(row, attr):
  attr_dict = attributes_to_dict(row['attributes'])
  if attr in attr_dict:
    return attr_dict[attr]
  else:
    return None
  
def gff_to_df(gff, type=None):
  headers = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
  gff_df = pd.read_csv(gff, sep='\t', names=headers, comment='#')
  if type:
    gff_df = gff_df.query('type == @type')
  gff_df['ID'] = gff_df.apply(attributes_to_column, axis=1, args=('ID',))
  gff_df['Parent'] = gff_df.apply(attributes_to_column, axis=1, args=('Parent',))
  return gff_df
  
def kde(v, steps=1000):
  kde = scipy.stats.gaussian_kde(v)
  vmin = v.min()
  vmax = v.max()
  step = (vmax-vmin)/steps
  x = np.arange(vmin,vmax,step)
  y = kde.evaluate(x)
  return x,y

def kde_peaks(x, y, frac=0.1, min_dist=20):
  """
  frac: height fraction of peak from previous highest peak to be keps
  min_dist: minimal distance on X axis between peaks
  """
  # find all peaks
  yrange = y.max() - y.min()
  peaks = scipy.signal.find_peaks(y, distance=min_dist, prominence=yrange/20)
  peaks_x = x[peaks[0]]
  peaks_y = y[peaks[0]]
  # convert peaks to (x,y pairs)
  peaks_pairs = list(zip(peaks_x, peaks_y))
  # sort peaks by descending y
  peaks_pairs_sort = sorted(peaks_pairs, key=lambda p: p[1], reverse=True)
  # iterate on peaks until peak y is smaller than frac of the previous peak
  i = 1
  while i < len(peaks_pairs_sort):
    if peaks_pairs_sort[i][1]/peaks_pairs_sort[i-1][1] >= frac:
      i += 1
    else:
      break
  return peaks_pairs_sort[:i]

def intron_len_ratios(trans_df, column):
  """
  Expects a DF of introns from a
  single transcript (result of groupby).
  Calculates several intron length ratios
  """
  colnames = ['first_from_total_ratio', 'longest_intron', 'first_to_max_nonfirst_ratio', 'max_to_min_ratio']
  # if single intron - return NA
  if trans_df.shape[0] == 1:
    return pd.Series([np.nan]*len(colnames), index=colnames)
  first_intron_len = trans_df.query('intron_index == 1')[column].iloc[0]
  total_intron_len = trans_df[column].sum()
  max_nonfirst_intron_len = trans_df.query('intron_index > 1')[column].max()
  max_intron_length = max(first_intron_len,max_nonfirst_intron_len)
  longest_intron = trans_df.query(f'{column} == {max_intron_length}')['intron_index'].values[0]
  min_intron_length = trans_df[column].min()
  first_from_total_ratio = (first_intron_len/total_intron_len)
  first_to_max_nonfirst_ratio = (first_intron_len/max_nonfirst_intron_len)
  max_to_min_ratio = max_intron_length/min_intron_length
  return pd.Series([first_from_total_ratio, longest_intron, first_to_max_nonfirst_ratio, max_to_min_ratio], index=colnames)

def simulate(df):
  """
  Simulate a DF of introns by
  permutations of the true one,
  so the number of introns of
  each rank is kept
  """
  n_trans = len(df['Parent'].unique())
  sim_trans_names = pd.Series([f'transSim{i}' for i in range(n_trans)])
  sim_df = df.copy()
  sim_df.sort_values(by='intron_index', inplace=True)
  rank_counts = sim_df['intron_index'].value_counts()
  parents = [sim_trans_names[:c].sample(frac=1) for c in rank_counts]
  parents = pd.concat(parents)
  parents.index = sim_df.index
  sim_df['Parent'] = parents
  return sim_df

def calc_stats(gff_df, column, name, query=None):
  
  stats_list = ['Min', 'Max','Mean', 'STD', 'Q10','Q25','Q50','Q75','Q90','Modes', 'M1_x', 'M1_y', 'M2_x', 'M2_y', 'Introns_count', 'Total_intron_length', 'Mean_exon', 'Transcripts_count', 'Total_transcript_length', 'Transcripts_containing_introns', 'Mean_per_transcript', 'Mean_total_intron_length_per_transcript', 'Mean_total_exon_length_per_transcript', 'Total_intron_fraction', 'Mean_intron_fraction', 'Mean_intron_ratio', 'mean_first_from_total_intronic_ratio', 'mean_first_to_max_nonfirst_intron_ratio', 'perc_first_intron_longest', 'mean_max_to_min_intron_ratio', 'mean_max_to_min_intron_ratio_perm_p', 'Dataset']

  # per-transcript stats
  trans_df = gff_df.query('type == "mRNA"')
  df = gff_df.query('type == "intron"')
  exon_df = gff_df.query('type == "exon"')

  df.sort_values(by=['seqid','start'], inplace=True)
  df_plus = df.query('strand == "+"')
  df_plus['intron_index'] = df_plus.groupby('attributes').cumcount()
  df_minus = df.query('strand == "-"')
  df_minus['intron_index'] = df_minus.groupby('attributes').cumcount(ascending=False)
  df = pd.concat([df_plus, df_minus])
  df['intron_index'] = df['intron_index'] + 1

  if query:
    df = df.query(query)
  if df.shape[0] == 0:
    stats = pd.Series([np.nan]*len(stats_list), index=stats_list)
    stats = pd.DataFrame(stats).transpose()
    return stats

  trans_count = trans_df.shape[0]
  trans_with_introns = df['Parent'].nunique()
  introns_per_trans = df.groupby('Parent')['intron_index'].max()
  single_exon_trans = set(trans_df['ID']) - set(introns_per_trans.index)
  single_exon_trans_introns = pd.Series(0, index=single_exon_trans)
  introns_per_trans = pd.concat([introns_per_trans, single_exon_trans_introns])
  mean_introns_per_trans = introns_per_trans.mean()

  trans_len = trans_df['length']
  trans_len.index = trans_df['ID']
  trans_len.name = 'transcript_length' 
  total_trans_len = trans_len.sum()
  
  intron_len_per_trans = df.groupby('Parent')[column].sum()
  mean_intron_len_per_trans = intron_len_per_trans.mean()
  exon_len_per_trans = exon_df.groupby('Parent')[column].sum()
  mean_exon_len_per_trans = exon_len_per_trans.mean()
  single_exon_trans_introns_len = pd.Series(0, index=single_exon_trans)
  intron_len_per_trans = pd.concat([intron_len_per_trans, single_exon_trans_introns_len])
  intron_len_per_trans.name = column
  intron_len_per_trans = pd.concat([intron_len_per_trans, trans_len], axis=1)
  intron_len_per_trans['intron_fraction'] = intron_len_per_trans[column]/intron_len_per_trans['transcript_length']
  mean_intron_frac = intron_len_per_trans['intron_fraction'].mean()
  intron_len_per_trans['exon_length'] = intron_len_per_trans['transcript_length'] - intron_len_per_trans[column]
  intron_len_per_trans['intron_ratio'] = intron_len_per_trans[column]/intron_len_per_trans['exon_length']
  mean_intron_ratio = intron_len_per_trans['intron_ratio'].mean()

  # introns stats

  vec = df[column]
  q10 = vec.quantile(0.1)
  q25 = vec.quantile(0.25)
  q50 = vec.quantile(0.5)
  q75 = vec.quantile(0.75)
  q90 = vec.quantile(0.9)
  min_ = vec.min()
  max_ = vec.max()
  mean = vec.mean()
  std = vec.std()

  # exons stats
  mean_exon = exon_df[column].mean()
  
  kde_x, kde_y = kde(vec)
  peaks = kde_peaks(kde_x, kde_y)

  n_modes = len(peaks)
  if n_modes > 0:
    m1_x, m1_y = peaks[0]
  else:
    m1_x, m1_y = (np.nan, np.nan)
  if n_modes > 1:
    m2_x, m2_y = peaks[1]
  else:
    m2_x, m2_y = (np.nan, np.nan)

  count = df.shape[0]
  mean_per_transcript = df.groupby('Parent')['intron_index'].max().mean()

  total_len = vec.sum()
  total_intron_fraction = total_len/total_trans_len

  intron_len_ratios_df = df.groupby('Parent').apply(intron_len_ratios, column=column)
  mean_first_from_toal_intronic_ratio = intron_len_ratios_df['first_from_total_ratio'].mean()
  mean_first_to_max_nonfirst_intron_ratio = intron_len_ratios_df['first_to_max_nonfirst_ratio'].mean()
  mean_max_to_min_intron_ratio = intron_len_ratios_df['max_to_min_ratio'].mean()
  try:
    perc_first_intron_longest = intron_len_ratios_df.query('longest_intron == 1').shape[0] / intron_len_ratios_df.query('longest_intron >= 1').shape[0] * 100
  except ZeroDivisionError:
    perc_first_intron_longest = np.nan

  nsim = 1000
  sim_stats = []
  for s in range(nsim):
    sim_df = simulate(df)
    max_intron = sim_df.groupby('Parent')[column].max()
    min_intron = sim_df.groupby('Parent')[column].min()
    max2min_ratio = max_intron/min_intron
    sim_stat = max2min_ratio.loc[max2min_ratio != 1].mean()
    sim_stats.append(sim_stat)
  sim_stats = pd.Series(sim_stats)
  mean_max_to_min_intron_ratio_perm_p = len(sim_stats.loc[sim_stats <= mean_max_to_min_intron_ratio])/nsim


  stats = pd.Series([min_, max_, mean, std, q10, q25, q50, q75, q90, n_modes, m1_x, m1_y, m2_x, m2_y, count, total_len, mean_exon,  trans_count, total_trans_len, trans_with_introns, mean_introns_per_trans, mean_intron_len_per_trans, mean_exon_len_per_trans, total_intron_fraction, mean_intron_frac, mean_intron_ratio, mean_first_from_toal_intronic_ratio, mean_first_to_max_nonfirst_intron_ratio, perc_first_intron_longest, mean_max_to_min_intron_ratio, mean_max_to_min_intron_ratio_perm_p, name], index=stats_list)

  stats = pd.DataFrame(stats).transpose()
  return stats

### MAIN

if __name__ == "__main__":

  in_gff = sys.argv[1]
  out_stats = sys.argv[2]
  species = sys.argv[3]

  gff_df = gff_to_df(in_gff)
  gff_df['length'] = gff_df['end'] - gff_df['start'] + 1
  gff_df['log_length'] = np.log10(gff_df['length'])

  all_introns_stats = calc_stats(gff_df, 'length', 'all')
  all_introns_log_stats = calc_stats(gff_df, 'log_length', 'all_log')
  #first_introns_stats = calc_stats(gff_df, 'length', 'first', query='intron_index == 1')
  #first_introns_log_stats = calc_stats(gff_df, 'log_length', 'first_log', query='intron_index == 1')
  #nonfirst_introns_stats = calc_stats(gff_df, 'length', 'nonfirst', query='intron_index > 1')
  #nonfirst_introns_log_stats = calc_stats(gff_df, 'log_length', 'nonfirst_log', query='intron_index > 1')

  #stats_df = pd.concat([all_introns_stats,all_introns_log_stats,first_introns_stats,first_introns_log_stats,nonfirst_introns_stats,nonfirst_introns_log_stats])
  stats_df = pd.concat([all_introns_stats,all_introns_log_stats])

  stats_df.index = [species]*stats_df.shape[0]

  stats_df.to_csv(out_stats, sep='\t') 
