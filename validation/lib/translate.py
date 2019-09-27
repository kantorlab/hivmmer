#!/usr/bin/env python
import argparse
import sys

from Bio import Seq
from Bio import SeqIO

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument("FASTA",
                      nargs=1,
                      help="path to nucleotide FASTA file, or '-' to read from stdin")
  parser.add_argument("--allow-stop-codons",
                      action="store_true",
                      default=False,
                      help="allow stop codons [default: False]")
  return parser.parse_args()

def translate(fasta, allow_stop_codons):

  nskipped = {"N": 0, "*": 0}

  for n, record in enumerate(SeqIO.parse(fasta, "fasta")):

    nstop = 0
    seq = str(record.seq)

    if 'N' in seq:
      nskipped["N"] += 1
      continue

    for i in range(3):
      aa = Seq.translate(seq[i:])
      if aa.count("*") > 0:
        nstop += 1
      else:
        print(">%s-%d" % (record.id, i))
        print(aa)

    seq = str(record.seq.reverse_complement())

    for i in range(3):
      aa = Seq.translate(seq[i:])
      if (not allow_stop_codons) and aa.count("*") > 0:
        nstop += 1
      else:
        print(">%s-%d'" % (record.id, i))
        print(aa)

    if nstop == 6:
      nskipped["*"] += 1

  print("nreads", n, file=sys.stderr)
  print("nskipped N", nskipped["N"], file=sys.stderr)
  print("nskipped *", nskipped["*"], file=sys.stderr)

if __name__ == "__main__":
  args = parse_args()
  fasta = args.FASTA[0]
  if fasta == '-':
    translate(sys.stdin, args.allow_stop_codons)
  else:
    with open(fasta) as f:
      translate(f, args.allow_stop_codons)

# vim: expandtab sw=2 ts=2
