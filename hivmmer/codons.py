"""
"""
import math
import pandas as pd
import sys
from Bio import SearchIO
from Bio import Seq
from Bio import SeqIO
from importlib import resources
from itertools import chain

# For more information on how thresholds were estimated,
# see `validation/README.md` in the hivmmer git repo.
thresholds = {
    "env": 18.0,
    "gag": 6.5,
    "nef": 4.4,
    "pol": 6.6,
    "tat": 0.8,
    "vif": 3.1,
    "vpr": 2.7,
    "vpu": 2.3
}

dblengths = {
    "env": 854*5226,
    "gag": 429*5183,
    "nef": 207*5188,
    "pol": 983*2956,
    "tat": 71*3361,
    "vif": 171*3937,
    "vpr": 89*3931,
    "vpu": 53*4573
}

def _load_hxb2(gene):
    """
    Load a pre-computed index that converts HMM position to
    HXB2 coordinates for `gene`.
    """
    with resources.open_binary("hivmmer", "{}.hxb2.tsv".format(gene)) as f:
        return pd.read_csv(f, sep="\t", usecols=["hmm", "hxb2", "ins", "del"], index_col="hmm")

def codons(readfile, hmmerfile, gene):
    """
    """

    dblength  = dblengths[gene]
    threshold = thresholds[gene]
    coords    = _load_hxb2(gene)

    reads = SeqIO.index(readfile, "fasta")
    hmmer = SearchIO.read(hmmerfile, "hmmer3-text")

    counts = dict((hxb2, {}) for hxb2 in range(coords.hxb2.iloc[0], coords.hxb2.iloc[-1] + 1, 3))

    for hit in hmmer.hits:

        # Skip hits that contain stop codons (indicates wrong frame)
        if "*" in chain(hsp.aln[1].seq for hsp in hit.hsps): continue

        id, _, frame = hit.id.rpartition("-")
        count = int(id.partition("-")[2])

        if frame.endswith("'"):
            seq = str(reads[id].seq.reverse_complement())
            offset = int(frame[:-1])
        else:
            seq = str(reads[id].seq)
            offset = int(frame)

        for hsp in hit.hsps:

            if math.log((dblength * hsp.hit_span) / 2**hsp.bitscore) >= threshold: continue

            i    = 0                         # tracks position in the alignment (0-indexed)
            hmm  = hsp.query_start + 1       # tracks position in the HMM reference sequence (1-indexed)
            read = 3*hsp.hit_start + offset  # tracks position in the read sequence (0-indexed)

            # read/reference sequences should have same length in alignment
            n = len(hsp.aln[0].seq)
            assert len(hsp.aln[1].seq) == n

            # alignment should not start or end with an insertion
            assert hsp.aln[0].seq[0] != "."
            assert hsp.aln[0].seq[n-1] != "."

            while i < n:

                aa_frame = []
                codon_frame = []
                hxb2 = coords.loc[hmm, "hxb2"]

                aa = hsp.aln[1].seq[i]
                if aa != "-":
                    aa_frame.append(aa)
                    codon_frame.append(seq[read:read+3])
                    read += 3

                # Extend frame with insertions relative to HMM
                while i < (n-1) and hsp.aln[0].seq[i+1] == ".":
                    aa_frame.append(hsp.aln[1].seq[i+1])
                    codon_frame.append(seq[read:read+3])
                    assert aa_frame[-1] != "-"
                    read += 3
                    i += 1

                # Extend frame with deletions relative to HXB2
                for _ in range(coords.loc[hmm, "del"]):

                    if i < (n-1):
                        aa = hsp.aln[1].seq[i+1]
                        if aa != "-":
                            aa_frame.append(aa)
                            codon_frame.append(seq[read:read+3])
                            read += 3
                            i += 1
                            hmm += 1

                    while i < (n-1) and hsp.aln[0].seq[i+1] == ".":
                        aa_frame.append(hsp.aln[1].seq[i+1])
                        codon_frame.append(seq[read:read+3])
                        assert aa_frame[-1] != "-"
                        read += 3
                        i += 1

                # Assign codons, adjusting coordinates relative to HXB2 insertions
                assert len(aa_frame) == len(codon_frame)
                for _ in range(coords.loc[hmm, "ins"] + 1):
                    # Iterate through the codon frame
                    if codon_frame:
                        codon = codon_frame.pop(0)
                        aa = aa_frame.pop(0)
                        if aa != 'X' and 'N' not in codon:
                            assert str(Seq.translate(codon)) == aa.upper()
                            counts[hxb2][codon] = counts[hxb2].get(codon, 0) + count
                    # Deletions occur when the codon frame is empty
                    else:
                        counts[hxb2][""] = counts[hxb2].get("", 0) + count
                    hxb2 += 3

                i += 1
                hmm += 1

    # output
    lines = []
    for hxb2 in sorted(counts):
        for codon in sorted(counts[hxb2], key=counts[hxb2].get, reverse=True):
            lines.append("\t".join([str(hxb2), codon, str(counts[hxb2][codon])]))
    return lines

# vim: expandtab sw=4 ts=4
