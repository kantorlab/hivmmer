"""
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

_genes = {
    "gag": 790,
    "pol": 2085,
    "vif": 5041,
    "vpr": 5559,
    "tat": 5831,
    "vpu": 6062,
    "env": 6225,
    "nef": 8797
}
_end = 9417

def plot_coverage(codonfile, outfile, min_coverage=10):
    """
    """
    codons = pd.read_csv(codonfile, sep="\t", index_col="hxb2").fillna("")
    coverage = codons.groupby(level=0)["count"].sum()

    # Initialize plot
    sns.set_style("ticks")
    fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    ax.set_title("Coverage", loc="left", fontsize=16)

    # Draw coverage
    xlim = (min(_genes.values()), _end)
    ax.fill_between(xlim, (min_coverage, min_coverage), color="r")
    ax.bar(coverage.index, coverage.values, width=1.0, color="k", edgecolor="none")

    # Setup x axis
    ax.set_xlabel("Gene and HXB2 Coordinates")
    ax.set_xlim(xlim)
    ax.set_xticks(list(_genes.values()))
    ax.set_xticklabels(["{} {}".format(gene[0], str(gene[1]).rjust(4)) for gene in _genes.items()], rotation=90, fontfamily="monospace")
    ax.grid(axis="x")

    # Setup y axis
    ax.set_ylabel("# of Reads")
    yticks = [0, min_coverage]
    if coverage.max() > 1.1*min_coverage:
        yticks.append(coverage.max())
    ax.set_ylim(yticks[0], yticks[-1])
    ax.set_yticks(yticks)

    # Finalize and write to file
    sns.despine(trim=True)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

