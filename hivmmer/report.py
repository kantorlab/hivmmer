"""
"""

import matplotlib
matplotlib.rcParams["font.sans-serif"] = ["Nimbus Sans L", "Helvetica", "Arial"]

import hivmmer
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import time
from collections import defaultdict
from datetime import datetime
from matplotlib.ticker import FixedLocator
from subprocess import run

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


def plot_coverage(codonfile, outfile, min_coverage=1000):
    """
    Write a PDF plot to `outfile` showing the coverage at each HXB2 position
    based on the codon counts in `codonfile`.
    """
    codons = pd.read_csv(codonfile, sep="\t", index_col="hxb2").fillna("")
    coverage = codons.groupby(level=0)["count"].sum()
    max_coverage = coverage.max()

    # Initialize plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 2))

    # Draw coverage
    xlim = (min(_genes.values()), _end)
    ax.fill_between(xlim, (min_coverage, min_coverage), color="gray")
    ax.fill_between(coverage.index, coverage.values, step="pre", color="k")

    # Setup x axis
    ax.set_xlim(xlim)
    ax.set_xticks(list(_genes.values()))
    ax.set_xticklabels(list(_genes.values()), rotation=90)
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(list(_genes.values()))
    ax2.set_xticklabels(list(_genes.keys()), rotation=90)
    ax.grid(axis="x")

    # Setup y axis
    ax.set_ylim(0, max(min_coverage, max_coverage))
    ax.set_yticks([min_coverage, max_coverage])

    # Finalize and write to file
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


def plot_coverage_prrt(aafile, outfile, min_coverage=1000):
    """
    Write a PDF plot to `outfile` showing the coverage in the PR/RT/IN
    region of the pol gene.
    """
    aa = pd.read_excel(aafile)

    # Initialize plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 1.5))

    # Serialize the X coordinates for PR, RT, IN
    aa["x"] = aa["position"]
    aa.loc[aa["region"] == "RT", "x"] = aa.loc[aa["region"] == "RT", "position"] + 100
    aa.loc[aa["region"] == "IN", "x"] = aa.loc[aa["region"] == "IN", "position"] + 541
    del aa["region"]
    del aa["position"]

    # Calculate total coverage
    coverage = aa.groupby("x")["coverage"].sum()
    max_coverage = coverage.max()

    # Plot coverage
    xlim = (1, 811)
    ax.fill_between(xlim, (min_coverage, min_coverage), color="gray")
    ax.fill_between(coverage.index, coverage.values, step="pre", color="k")

    # Setup x axis
    xticks = [1, 101, 542]
    xlabels = ["PR", "RT", "IN"]
    ax.set_xlim(*xlim)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=12)
    ax.grid(axis="x")

    # Setup y axis
    ax.set_ylim(0, max(min_coverage, max_coverage))
    ax.set_yticks([min_coverage, max_coverage])

    # Finalize and write to file
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


def plot_drms(aafile, drmfile, column, outfile):
    """
    Write a PDF plot to `outfile` showing DRMs identified by `column`.
    Return a dictionary of DRMs by PI/NRTI/NNRTI/INSTI.
    """
    aa = pd.read_excel(aafile)
    drm = pd.read_csv(drmfile)
    drm = drm[drm[column] == 1]

    # Initialize plot
    fig, ax = plt.subplots(1, 1, figsize=(15, 3))

    # Serialize the X coordinates for PR, RT, IN
    aa["x"] = aa["position"]
    aa.loc[aa["region"] == "RT", "x"] = aa.loc[aa["region"] == "RT", "position"] + 100
    aa.loc[aa["region"] == "IN", "x"] = aa.loc[aa["region"] == "IN", "position"] + 541
    del aa["region"]
    del aa["position"]

    drm["x"] = drm["position"]
    drm.loc[drm["region"] == "RT", "x"] = drm.loc[drm["region"] == "RT", "position"] + 100
    drm.loc[drm["region"] == "IN", "x"] = drm.loc[drm["region"] == "IN", "position"] + 541

    # Calculate total coverage
    coverage = aa.groupby("x")["coverage"].sum()

    # Unpivot and remove zero- or low-coverage calls
    del aa["hxb2"]
    aa = aa.melt(id_vars=["x", "coverage"], var_name="variant", value_name="count")
    aa["frequency"] = aa["count"] / aa["coverage"]

    # X axis
    margin = 1.5
    xticks = [1] + list(range(10, 100, 10)) + [101] + list(range(110, 540, 10)) + [542] + list(range(551, 811, 10))
    xlabels = [1] + list(range(10, 100, 10)) + [1] + list(range(10, 440, 10)) + [1] + list(range(10, 269, 10))
    ax.set_xlim(1-margin, 811+margin)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=9, rotation=90)
    ax.set_xlabel("AA Position", size=14)
    ax.xaxis.set_minor_locator(FixedLocator(list(range(1, 811))))
    ax.xaxis.set_tick_params(direction="out")

    # Y axis
    yticks = [0.1, 1, 5, 10, 20, 50, 100]
    ax.set_ylabel("AA Frequency", size=14)
    ax.set_yscale("log")
    ax.set_ylim(0.09, 110)
    ax.set_yticks(yticks)
    ax.set_yticklabels(['%g%%' % i for i in yticks], size=14)
    ax.get_yaxis().set_tick_params(direction="out")

    # Grid
    ax.grid(True, which="major")
    ax.axvline(x=100, lw=1.0, color="k")
    ax.axvline(x=541, lw=1.0, color="k")

    # Scatterplot
    ax.scatter(aa["x"], 100*aa["frequency"], ec="k", alpha=0.5, marker=".", fc="none")

    # Annotate DRMs
    if len(drm) > 0:
        ax.scatter(drm["x"], 100*drm["frequency"], fc="r", marker=".", ec="none")
        ax.scatter(drm["x"], 100*drm["frequency"], ec="r", marker="o", fc="none")

    # Finalize and write to file
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

    drms = defaultdict(list)
    for row in drm.itertuples():
        drms[row.drug].append("".join(map(str, [row.consensus, row.position, row.variant])))
    return dict((drug, ", ".join(drms[drug])) for drug in drms)


def compile(fastq, coveragefile, coverageprrtfile,
            drmsfile, drmifile, sdrmfile,
            drms, drmi, sdrm,
            workdir, outfile):
    """
    Write a PDF report to `outfile` containing coverage, DRM, and SDRM plots.
    """

    name = os.path.basename(fastq[0]).partition("_")[0]

    # Summarize FASTQ input files
    fastq_summary = []
    for f in fastq:
        st = os.stat(f)
        mb = st.st_size / 1048576
        ts = time.strftime("%-d %b %Y", time.localtime(st.st_mtime))
        fastq_summary.append(r"%s (%.1f MB, %s) \\" % (f, mb, ts))
    fastq_summary[-1] = fastq_summary[-1] + "[3pt]"
    fastq_summary = "\n".join(fastq_summary)

    tex = r"""\documentclass{article}
\usepackage[letterpaper, margin=0.5in]{geometry}
\usepackage{graphicx}
\usepackage{fontspec}
\setmainfont{Nimbus Sans L}
\usepackage[small]{titlesec}
\titlespacing*{\section}{0pt}{0pt}{0pt}
\titlespacing*{\subsection}{0pt}{0pt}{0pt}
\setcounter{secnumdepth}{0}

\begin{document}
\pagenumbering{gobble}
\footnotesize

\section{Dataset: %s}
{\em Report generated on %s by hivmmer %s from input:} \\
%s

\subsection{Coverage: Whole Genome}
\includegraphics[width=7.5in]{%s}

\subsection{Coverage: PR/RT/IN}
\includegraphics[width=7.5in]{%s}

\footnotesize
\subsection{DRMs: Stanford}
\textbf{PI}: %s $\;$ \textbf{NRTI}: %s $\;$ \textbf{NNRTI}: %s $\;$ \textbf{INSTI}: %s \\
\includegraphics[width=7.5in]{%s}

\subsection{DRMs: IAS}
\textbf{PI}: %s $\;$ \textbf{NRTI}: %s $\;$ \textbf{NNRTI}: %s $\;$ \textbf{INSTI}: %s \\
\includegraphics[width=7.5in]{%s}

\subsection{SDRMs: Stanford}
\textbf{PI}: %s $\;$ \textbf{NRTI}: %s $\;$ \textbf{NNRTI}: %s $\;$ \textbf{INSTI}: %s \\
\includegraphics[width=7.5in]{%s}
\end{document}
""" % (name, datetime.now().ctime(), hivmmer.__version__, fastq_summary, coveragefile, coverageprrtfile,
       drms.get("PI", "n/a"), drms.get("NRTI", "n/a"), drms.get("NNRTI", "n/a"), drms.get("INSTI", "n/a"),
       drmsfile,
       drmi.get("PI", "n/a"), drmi.get("NRTI", "n/a"), drmi.get("NNRTI", "n/a"), drmi.get("INSTI", "n/a"),
       drmifile,
       sdrm.get("PI", "n/a"), sdrm.get("NRTI", "n/a"), sdrm.get("NNRTI", "n/a"), sdrm.get("INSTI", "n/a"),
       sdrmfile)

    # Write tex to flie
    texfile = os.path.join(workdir, "report.tex")
    with open(texfile, "w") as f:
        f.write(tex.replace("_", "\_"))

    # Run tectonic
    with open(os.path.join(workdir, "tectonic.log"), "w") as f:
        run(["tectonic", texfile], stdout=f, stderr=f, check=True)

    outfile = os.path.abspath(outfile)
    os.rename(os.path.join(workdir, "report.pdf"), outfile)
    return outfile


# vim: expandtab sw=4 ts=4
