"""
Calculate pairwise distances
between species, based on
intron length vectors.
We use the KS statistic as
the distance metric.
"""

import sys
import os
from itertools import combinations
import pandas as pd
from scipy import stats

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
  gff_df['length'] = gff_df['end'] - gff_df['start'] + 1
  return gff_df

def get_stats(gff_df):
  introns_df = gff_df.query('type == "intron"')
  trans_df = gff_df.query('type == "mRNA"')
  exon_df = gff_df.query('type == "exon"')

  # intron lengths
  intron_lengths = introns_df['length']

  # intron counts
  trans_with_introns = introns_df['Parent'].nunique()
  introns_per_trans = introns_df.groupby('Parent')['type'].count()
  single_exon_trans = set(trans_df['ID']) - set(introns_per_trans.index)
  single_exon_trans_introns = pd.Series(0, index=single_exon_trans)
  introns_per_trans = pd.concat([introns_per_trans, single_exon_trans_introns])
    
  # intron fractions
  intron_len_per_trans = introns_df.groupby('Parent')['length'].sum()
  exon_len_per_trans = exon_df.groupby('Parent')['length'].sum()
  single_exon_trans_introns_len = pd.Series(0, index=single_exon_trans)
  intron_len_per_trans = pd.concat([intron_len_per_trans, single_exon_trans_introns_len])
  intron_len_per_trans.name = 'length'
  trans_len = trans_df['length']
  trans_len.index = trans_df['ID']
  trans_len.name = 'transcript_length'
  intron_len_per_trans = pd.concat([intron_len_per_trans, trans_len], axis=1)
  intron_len_per_trans['intron_fraction'] = intron_len_per_trans['length']/intron_len_per_trans['transcript_length']

  return {'intron_lengths': intron_lengths, 'intron_counts': introns_per_trans, 'intron_fractions': intron_len_per_trans['intron_fraction']}

if __name__ == "__main__":

  out_dir = sys.argv[1]
  gff_list = sys.argv[2:]

  print('Parsing GFF files...')
  stat_vectors = {}
  for gff in gff_list:
    sp = os.path.basename(os.path.dirname(gff))
    print(sp)
    gff_df = gff_to_df(gff)
    stat_vectors[sp] = get_stats(gff_df)

  print('Calculating pairwise distances...')
  stat_names = ['intron_lengths', 'intron_counts', 'intron_fractions']
  for st in stat_names:
    print('STAT: %s' % st)
    out_matrix = os.path.join(out_dir, st+'_KS_dist.tsv')
    ks_similarity = []
    for pair in combinations(stat_vectors.keys(),2):
      sp1, sp2 = pair
      print(sp1, sp2)
      ks = stats.kstest(stat_vectors[sp1][st], stat_vectors[sp2][st]).statistic
      res = pd.Series([sp1,sp2,ks])
      ks_similarity.append(res)
      res_rev = pd.Series([sp2,sp1,ks])
      ks_similarity.append(res_rev)
  
    for sp in stat_vectors.keys():
      dist_to_self = pd.Series([sp,sp,0])
      ks_similarity.append(dist_to_self)

    ks_similarity_df = pd.concat(ks_similarity, axis=1).T
    ks_similarity_df.columns = ['species1','species2','KS']
    ks_similarity_df.set_index('species1', inplace=True)
    matrix = pd.pivot_table(ks_similarity_df, values='KS', index='species1', columns='species2')
    matrix.to_csv(out_matrix, sep='\t')
