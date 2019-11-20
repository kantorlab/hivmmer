"""
"""

import pandas as pd
from Bio import Seq

_ambiguous = dict(("".join(sorted(b)), a) for a, b in Seq.IUPAC.IUPACData.ambiguous_dna_values.items())
_frequencies = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.4]

def consensus(codonfile, outfile, min_coverage=10):
    """
    Write FASTA `outfile` with consensus sequences at varying
    thresholds, using only the variants above `min_coverage`.
    """
    codons = pd.read_csv(codonfile, sep="\t", index_col="hxb2").fillna("")

    # Filter out low-coverage variants
    codons = codons[codons["count"] >= min_coverage]

    # Normalize frequencies
    norm = 1.0 / codons.groupby(level=0)["count"].sum()
    norm = norm.to_frame(name="norm")
    codons = codons.join(norm, how="left")
    codons["frequency"] = codons["count"] * codons["norm"]
    
    index = codons.index.tolist()
    with open(outfile, "w") as f:
        for freq in _frequencies:
            codons = codons[codons["frequency"] >= freq]
            print(">{}".format(freq), file=f)
            for i in index:
                if not i in codons.index:
                    f.write("---")
                else:
                    for nt in zip(*codons.loc[[i], "codon"].tolist()):
                        f.write(_ambiguous["".join(sorted(set(nt)))])
            f.write("\n")


# vim: expandtab sw=4 ts=4
