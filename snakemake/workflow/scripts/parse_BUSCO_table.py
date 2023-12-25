import sys
import numpy as np
import pandas as pd
import gffutils

busco_tbl = sys.argv[1]
canon_tbl = sys.argv[2]
in_gff = sys.argv[3]
out_tbl = sys.argv[4]

busco_df = pd.read_csv(busco_tbl, sep='\t', comment='#', usecols=[0,2,3], names=['Busco_id', 'Sequence', 'Score'])

def select_best_score(busco_group):
    # if missing or single-copy BUSCO, return as is
    if busco_group.shape[0] == 1:
      return busco_group.iloc[0,:]
    # if duplicated BUSCO, choose row with max score
    max_score = busco_group['Score'].max()
    return busco_group.query('Score == @max_score').iloc[0,:]

#busco_df = busco_df.groupby('Busco_id').apply(select_best_score)
busco_df = busco_df.groupby('Busco_id').max()
busco_df.reset_index(inplace=True)

canon_df = pd.read_csv(canon_tbl, sep='\t', header=None, names=['mRNA', 'canon_mRNA'])
busco_df = busco_df.merge(canon_df, how='left', left_on='Sequence', right_on='mRNA')

db_path = in_gff + '.db'
gff_db = gffutils.create_db(in_gff, db_path, force=True, merge_strategy="create_unique")
gff = gffutils.FeatureDB(db_path)

def get_stats(row):
    mrna_id = row['canon_mRNA']
    if pd.isnull(mrna_id):
        return [np.nan]*7
    mrna = gff[mrna_id]
    mrna_strand = mrna.strand
    mrna_len = mrna.end - mrna.start + 1
    cds = list(gff.children(mrna_id, featuretype='CDS'))
    tot_cds_len = sum([c.end-c.start+1 for c in cds])
    prot_len = tot_cds_len/3
    introns = list(gff.children(mrna_id, featuretype='intron'))
    intron_len = [i.end-i.start+1 for i in introns]
    if mrna_strand == '+':
        mrna_start = mrna.start
        intron_pos = [i.start-mrna_start for i in introns]
    elif mrna_strand == '-':
        introns.reverse()
        mrna_start = mrna.end
        intron_pos = [mrna_start-i.end for i in introns]
        intron_len.reverse()

    n_introns = len(intron_len)
    tot_intron_len = sum(intron_len)
    intron_frac = tot_intron_len/mrna_len
    intron_len_str = ','.join(map(str,intron_len))
    intron_pos_str = ','.join(map(str,intron_pos))

    return mrna_len, n_introns, tot_intron_len, intron_frac, intron_len_str, intron_pos_str, prot_len

busco_df[['transcript_len', 'n_introns', 'tot_intron_len', 'intron_frac', 'intron_lens', 'intron_pos', 'Protein_length']] = busco_df.apply(get_stats, axis=1, result_type='expand')
busco_df.set_index('Busco_id')

busco_df.to_csv(out_tbl, sep='\t')
