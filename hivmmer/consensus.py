"""
"""

import argparse
import hivmmer
import pandas as pd
from Bio import Seq

_ambiguous = dict(("".join(sorted(b)), a) for a, b in Seq.IUPAC.IUPACData.ambiguous_dna_values.items())
_frequencies = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.4]

def consensus(codonfile, outfile, min_coverage=1000):
    """
    Write FASTA `outfile` with consensus sequences at varying
    thresholds, using only the variants above `min_coverage`.
    """
    codons = pd.read_csv(codonfile, sep="\t", index_col="hxb2").fillna("")

    # Sum counts by site
    sums = codons.groupby(level=0)["count"].sum()
    sums = sums.to_frame(name="sum")
    codons = codons.join(sums, how="left")

    # Filter out low-coverage sites
    codons = codons[codons["sum"] >= min_coverage]

    # Normalize counts
    codons["frequency"] = codons["count"] * (1.0 / codons["sum"])

    hxb2 = hivmmer.list_hxb2()
    with open(outfile, "w") as f:
        for freq in _frequencies:
            codons = codons[codons["frequency"] >= freq]
            print(">{}".format(freq), file=f)
            seq = []
            for i in hxb2:
                if not i in codons.index:
                    seq.append("---")
                else:
                    for nt in zip(*codons.loc[[i], "codon"].tolist()):
                        seq.append(_ambiguous["".join(sorted(set(nt)))])
            print("".join(seq).strip("-"), file=f)

def _run():

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--min-coverage",
                        default=1000,
                        metavar="N",
                        type=int,
                        help="minimum coverage for sites included in the consensus (default: 1000)")
    parser.add_argument("CODONS",
                        help="input tab-separated file with codon frequencies")
    parser.add_argument("FASTA",
                        help="output FASTA file with consensus sequences")
    args = parser.parse_args()

    consensus(args.CODONS, args.FASTA, args.min_coverage)

# vim: expandtab sw=4 ts=4
