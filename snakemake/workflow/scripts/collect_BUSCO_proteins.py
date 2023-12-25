"""
Given a list of species, collect
orthologs into FASTA files
(one file per BUSCO)
Removes outliers based on
protein length and total
intron length MADs.
"""

import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation
from Bio import SeqIO

in_species_list = sys.argv[1]
in_dir = sys.argv[2]
out_dir = sys.argv[3]

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

# Read species list
with open(in_species_list) as f:
    species_list = [l.strip() for l in f.readlines()]

# Create DF with BUSCO names and SeqRecords

def binom_species(s):
    genus, species = s.split('_')[:2]
    return f'{genus}_{species}'

busco_prot = []
for sp in species_list:
    sp_binom = binom_species(sp)
    prot_fasta = os.path.join(in_dir, sp, 'prot.canon.fasta')
    prot_records = SeqIO.to_dict(SeqIO.parse(prot_fasta, 'fasta'))
    for rec_id in prot_records:
        prot_records[rec_id].id = f'{sp_binom}__{rec_id}'
        prot_records[rec_id].description = ''
    sp_busco_stats = os.path.join(in_dir, sp, 'BUSCO.stats')
    sp_busco_stats_df = pd.read_csv(sp_busco_stats, sep='\t', usecols=[1,5,8,12])
    sp_busco_stats_df['species'] = sp_binom
    sp_busco_stats_df['species_full'] = sp
    sp_busco_stats_df['prot_record'] = sp_busco_stats_df['canon_mRNA'].apply(lambda id: prot_records[id] if not pd.isnull(id) else np.nan)
    busco_prot.append(sp_busco_stats_df)
busco_prot_df = pd.concat(busco_prot)
# write table to file - for debugging
out_tsv = os.path.join(out_dir, 'BUSCO_stats.tsv')
out_cols = busco_prot_df.columns[:-1]
busco_prot_df[out_cols].to_csv(out_tsv, sep='\t')

# Extract proteins per BUSCO

MIN_SIZE = 5

def mad(s):
    s = s.dropna()
    return median_abs_deviation(s.values)

def mad_outliers(s):
    smed = s.median()
    smad = mad(s)
    smed10 = 0.1*smed
    x = max(smed10, 3*smad)
    lower = smed - x
    upper = smed + x
    return lower, upper

def create_BUSCO_fasta(busco_id):
    busco_id_df = busco_prot_df.query('Busco_id == @busco_id').dropna()
    prot_len_mad_min, prot_len_mad_max = mad_outliers(busco_id_df['Protein_length'])
    intron_len_mad_min, intron_len_mad_max = mad_outliers(busco_id_df['tot_intron_len'])
    busco_id_df = busco_id_df.query('Protein_length > @prot_len_mad_min & Protein_length < @prot_len_mad_max & tot_intron_len > @intron_len_mad_min & tot_intron_len < @intron_len_mad_max')
    busco_id_prot_records_mad = busco_id_df['prot_record'].to_list()
    busco_id_fasta = os.path.join(out_dir, busco_id + '.faa')
    if len(busco_id_prot_records_mad) < MIN_SIZE:
        busco_id_prot_records_mad = []
    SeqIO.write(busco_id_prot_records_mad, busco_id_fasta, 'fasta')

for busco_id in busco_prot_df['Busco_id'].unique():
    create_BUSCO_fasta(busco_id)
