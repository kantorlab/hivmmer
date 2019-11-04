"""
"""
import math
import pandas as pd
import sys
from Bio import SearchIO
from Bio import Seq
from Bio import SeqIO
from importlib import resources

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

def codons(readfile, hmmerfile, gene, out=sys.stdout):
    """
    """

    dblength  = dblengths[gene]
    threshold = thresholds[gene]
    hxb2      = _load_hxb2(gene)

    reads = SeqIO.index(readfile, "fasta")
    hmmer = SearchIO.read(hmmerfile, "hmmer3-text")

    counts = [{} for _ in range(hxb2.hxb2.min(), hxb2.hxb2.max()+1, 3)]

    for hit in hmmer.hits:

        id, _, frame = hit.id.rpartition("-")
        count = int(id.partition("-")[2])

        if frame.endswith("'"):
            seq = reads[id].seq.reverse_complement()
            offset = int(frame[:-1])
        else:
            seq = str(reads[id].seq)
            offset = int(frame)

        for hsp in hit.hsps:

            if math.log((dblength * hsp.hit_span) / 2**hsp.bitscore) >= threshold: continue

            i    = 0                         # tracks position in the alignment
            hmm  = hsp.query_start           # tracks position in the HMM reference sequence
            read = 3*hsp.hit_start + offset  # tracks position in the read sequence

            n = len(hsp.aln[0].seq)

            # read/reference sequences should have same length in alignment
            assert len(hsp.aln[1].seq) == n

            # alignment should not start with an insertion or deletion
            assert hsp.aln[0].seq[0] != "."
            assert hsp.aln[1].seq[0] != "-"

            while i < n:

                hmm_aa_frame  = [hsp.aln[1].seq[i]]
                read_aa_frame = [hsp.aln[0].seq[i]]
                read_nt_frame = [seq[read:read+3]]

                # Find insertions relative to HMM
                while i < (n-1) and hsp.aln[0].seq[i+1] == ".":
                    read_aa_frame.append(hsp.aln[1].seq[i+1])
                    read += 3
                    read_nt_frame.append(seq[read:read+3])
                    i += 1

                # Find deletions relative to HMM
                while i < (n-1) and hsp.aln[1].seq[i+1] == "-":
                    hmm_frame.append(hsp.aln[0].seq[i+1])
                    i += 1

                # Adjust for HXB2 insertions
                

                # Adjust for HXB2 deletions

                # Test frame agreement
                    assert read_aa.islower() or read_aa == "*"
                    read_codon = seq[read:read+3]
                    assert read_aa == Seq.translate(read_codon)

                # Assign codons

    # output
    print("hxb2", "codon", "count", sep="\t", file=out)
    for i, count in enumerate(counts):
        for codon in sorted(count, key=count.get, reverse=True):
            print(i, codon, count[codon], sep="\t", file=out)


# vim: expandtab sw=4 ts=4
