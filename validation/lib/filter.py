import numpy as np
import pandas as pd
import sys
from Bio import SeqIO

threshold, csv_file, fa_file = sys.argv[1:]

dblength = {
    "env": 854*5226,
    "gag": 429*5183,
    "nef": 207*5188,
    "pol": 983*2956,
    "tat": 71*3361,
    "vif": 171*3937,
    "vpr": 89*3931,
    "vpu": 53*4573
}

ids = pd.read_csv(csv_file, usecols=["target_name", "query_name", "score", "target_length"])
ids["dblength"] = ids["query_name"].apply(dblength.get)
ids["evalue"] = np.log(ids["dblength"] * ids["target_length"] / np.exp2(ids["score"]))
ids = frozenset(ids.loc[ids["evalue"] <= float(threshold), "target_name"])

for seq in SeqIO.parse(open(fa_file), "fasta"):
    if seq.id in ids:
        print(">{}".format(seq.id))
        print(seq.seq)
