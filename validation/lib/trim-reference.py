#!/usr/bin/env python
import sys
from Bio import AlignIO, SeqIO, Seq
from io import StringIO

aafile, start, end, gapped = sys.argv[1:]

# fix short sequences
records = list(SeqIO.parse(aafile, "fasta"))
maxlen = max(len(record.seq) for record in records)

# pad sequences so that they all have the same length
for record in records:
  if len(record.seq) != maxlen:
    record.seq = Seq.Seq(str(record.seq).ljust(maxlen, "-"))
assert all(len(record.seq) == maxlen for record in records)

# write to temp string
tmp = StringIO()
SeqIO.write(records, tmp, "fasta")
tmp.seek(0)

# open temp string
align = AlignIO.read(tmp, "fasta")

# adjust 1-index include to 0-index exclusive range
start = int(start) - 1
end = int(end)

# find hxb2 alignment
for a in align:
  if a.id == "B.FR.83.HXB2_LAI_IIIB_BRU.K03455":
    hxb2 = a
    break

# remove gaps in hxb2 sequence
coords = [i for i, aa in enumerate(hxb2.seq) if aa != "-" and aa != "#"]

# trim to start/end coordinates
coords = coords[start:end]

if gapped == "gapped":
  for a in align:
    aa = str(a.seq[coords[0]:coords[-1]]).replace("#", "-")
    if not "*" in aa:
      print(">%s" % a.description)
      print(aa)
else:
  for a in align:
    aa = "".join(a.seq[i] for i in coords).replace("#", "-")
    if not "*" in aa:
      print(">%s" % a.description)
      print(aa)

# vim: expandtab sw=2 ts=2
