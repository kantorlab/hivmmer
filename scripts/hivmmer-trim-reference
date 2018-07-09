#!/usr/bin/env python
import sys

from Bio import AlignIO
from Bio import Seq

align = AlignIO.read(sys.argv[1], "fasta")

# find hxb2 alignment
for a in align:
  if a.id == "B.FR.83.HXB2_LAI_IIIB_BRU.K03455":
    hxb2 = a
    break

# find gaps in hxb2 sequence
coords = [i for i, nt in enumerate(hxb2.seq) if nt != '-']

# trim to hxb2 coordinates, and translate sequence
for a in align:
  if a.id == "B.FR.83.HXB2_LAI_IIIB_BRU.K03455": continue
  aa = Seq.translate(''.join(a.seq[i] for i in coords).replace('-', 'N'))
  if not "*" in aa:
    print(">%s" % a.description)
    print(aa)
