import pandas as pd
import sys
from Bio import SeqIO

threshold, csv_file, fa_file = sys.argv[1:]

ids = pd.read_csv(csv_file, usecols=["query_name", "evalue"])
ids = frozenset(ids.loc[ids.evalue <= float(threshold), "query_name"])

for seq in SeqIO.parse(open(fa_file), "fasta"):
    if seq.id in ids:
        print(">{}".format(seq.id))
        print(seq.seq)
