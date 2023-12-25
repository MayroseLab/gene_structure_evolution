"""
Given a fasta file and a list of IDs,
only prints out records with rec.id
included in the list.
The list may be given as a file in the
second command line argument, or taken
from STDOIN
If the list contains a second field
(separate by tab), then it is used
for record name replacement.
"""

from Bio import SeqIO
import sys
import re

# if ids file given
if len(sys.argv) == 3:
  inp = open(sys.argv[2]).readlines()
# if should be taken from STDIN
elif len(sys.argv) == 2:
  inp = sys.stdin
id_dict = {}
for l in inp:
  l = l.strip()
  fields = l.split('\t')
  if len(fields) == 0:
    continue
  elif len(fields) == 1:
    id_dict[fields[0]] = fields[0]
  else:
    id_dict[fields[0]] = fields[1]

rec_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))
dot_regex = re.compile('.+\.\d+$')
rec_dict_nodot = {rec_id.split('.')[0] : rec_dict[rec_id] for rec_id in rec_dict}
valid_nodot = (len(rec_dict) == len(rec_dict_nodot))
c = 0
for mrna_id in id_dict:
  if mrna_id in rec_dict:
    rec = rec_dict[mrna_id]
  elif mrna_id.split('.')[0] in rec_dict:
    rec = rec_dict[mrna_id.split('.')[0]]
  elif not dot_regex.match(mrna_id) and valid_nodot and mrna_id in rec_dict_nodot:
    rec = rec_dict_nodot[mrna_id]
    rec.id = mrna_id
  else:
    continue
  rec.name = ''
  rec.description = ''
  print(rec.format('fasta').strip())
  c += 1

if c != len(id_dict):
  sys.exit(1)
  print(f'Something went wrong: {len(id_dict)} transcript IDs provided, but {c} records written.')
