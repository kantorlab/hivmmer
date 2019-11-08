"""
"""

import numpy as np
import os
import pandas as pd
import sys
from Bio import Seq

aa_header = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*", "X", "del", "ins"]

ranges = {
    "PR": range(2253, 2550, 3),
    "RT": range(2550, 3870, 3),
    "IN": range(4230, 5037, 3)
}

def aa_table(codonfile, outfile):
    """
    """
    codons = pd.read_csv(codonfile, sep="\t", index_col="hxb2").fillna("")
    subtables = []

    for region in ranges:
        subtable = pd.DataFrame(columns=["region", "position", "coverage"]+aa_header,
                                index=ranges[region]).fillna(0)
        subtable["region"] = region
        subtable["position"] = np.arange(1, len(subtable)+1)
        for hxb2 in subtable.index:
            if hxb2 in codons.index:
                rows = codons.loc[[hxb2]]
                for codon, count in zip(rows["codon"], rows["count"]):
                    if codon == "":
                        subtable.loc[hxb2, "del"] += count
                    elif len(codon) > 3:
                        subtable.loc[hxb2, "ins"] += count
                        aa = str(Seq.translate(codon[:3]))
                        subtable.loc[hxb2, aa] += count
                    else:
                        aa = str(Seq.translate(codon))
                        subtable.loc[hxb2, aa] += count
        subtable["coverage"] = subtable[aa_header].sum(axis=1)
        subtables.append(subtable)

    pd.concat(subtables).to_excel(outfile, index_label="hxb2")

# vim: expandtab sw=4 ts=4
