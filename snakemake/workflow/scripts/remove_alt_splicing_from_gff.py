"""
Takes a gff file and for genes with
multiple mRNAs, only keeps the longest
(largest sum of CDSs lengths).
Also removes UTRs and non protein-
coding genes.
Prints out a new gff.
"""

from __future__ import print_function
import gffutils
import sys

in_gff = sys.argv[1]
out_gff = sys.argv[2]
gene_mrna_out = out_gff + '.gene_to_mRNA'
mrna_canon_out = out_gff + '.mRNA_to_canon'

db_path = in_gff + '.db'
gff_db = gffutils.create_db(in_gff, db_path, force=True, merge_strategy="create_unique")
gff = gffutils.FeatureDB(db_path)

with open(out_gff, 'w') as fo, open(gene_mrna_out,'w') as fo2, open(mrna_canon_out,'w') as fo3:
  for feature in gff.all_features():
    if feature.featuretype not in {'gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'}:
      continue
    if feature.featuretype != 'gene':	# mRNA and exon features
      continue
    gene = feature
    # find longest mRNA (max sum of CDSs)
    mrnas = list(gff.children(gene, featuretype='mRNA'))
    if len(mrnas) == 0:
      continue
    longest_transcript = None
    max_len = 0
    for mrna in mrnas:
      cdss = list(gff.children(mrna, featuretype='CDS'))
      if len(cdss) == 0:
        continue
      total_cds_len = sum([cds.end - cds.start + 1 for cds in cdss])
      if total_cds_len > max_len:
        max_len = total_cds_len
        longest_transcript = mrna
    # if no CDSs found in any mRNA - this is not a protein-coding gene
    if not longest_transcript:
      continue
    
    # print gene - canon transcript matching
    print("%s\t%s" %(gene['ID'][0], longest_transcript['ID'][0]),file=fo2)
    # print transcript - canon transcript matching
    for mrna in mrnas:
      print("%s\t%s" %(mrna['ID'][0],longest_transcript['ID'][0]), file=fo3)

    # remove UTRs from gene, mRNA, and exons
    # find start and end coordinates (start of first and end of last CDSs)
    cds_coords = [(cds.start, cds.end) for cds in gff.children(longest_transcript, featuretype='CDS')]
    mrna_start = min([c[0] for c in cds_coords])
    mrna_end = max([c[1] for c in cds_coords])
    # adjust gene and mRNA coords and print
    gene.start = mrna_start
    gene.end = mrna_end
    print(str(gene), file=fo)
    longest_transcript.start = mrna_start
    longest_transcript.end = mrna_end
    print(str(longest_transcript), file=fo)
    # adjust CDS and exon coords, remove UTR features and UTR exons
    for feat in gff.children(longest_transcript):
      if feat.end < mrna_start or feat.start > mrna_end:	# 5'/3' UTR
        continue
      elif feat.start < mrna_start:	# exon partly 5' UTR
        feat.start = mrna_start
      elif feat.end > mrna_end:	# exon partly 3' UTR
        feat.end = mrna_end
      print(str(feat), file=fo)
