"""
"""

import math
import sys
from Bio import SearchIO
from Bio import Seq
from Bio import SeqIO

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

def codons(readfile, hmmerfile, gene, outfile, log=sys.stderr):
    """
    """

    dblength  = dblengths[gene]
    threshold = thresholds[gene]

    print("Indexing reads in:", readfile, file=log)
    reads = SeqIO.index(readfile, "fasta")

    print("Parsing hmmsearch output in:", hmmerfile, file=log)
    hmmer = SearchIO.read(hmmerfile, "hmmer3-text")

    print("Counting codons", file=sys.stderr)

    # All genes are <1000 AA
    counts = [{} for _ in range(1000)]

    for hit in hmmer.hits:

        id, _, frame = hit.id.rpartition("-")
        count = int(id.partition("-")[2])

        if frame.endswith("'"):
            seq = reads[id].seq.reverse_complement()
            offset = int(frame[:-1])
        else:
            seq = reads[id].seq
            offset = int(frame)

        for hsp in hit.hsps:

            if math.log((dblength * hsp.hit_span) / 2**hsp.bitscore) >= threshold: continue

            i = 3*hsp.hit_start + offset # tracks position in the read sequence
            j = 0                        # tracks position in the reference sequence
            k = 0                        # tracks position in the hmmer alignment

            n = len(hsp.aln[0].seq)

            # read/reference sequences should have same length in alignment
            assert len(hsp.aln[1].seq) == n

            # alignment should not start with an insertion
            assert hsp.aln[0].seq[0] != "."

            while k < n:
                aa = hsp.aln[1].seq[k]
                if aa == "-":
                    # deletion
                    counts[hsp.query_start+j][""] = counts[hsp.query_start+j].get("", 0) + count
                    j += 1 # increment the reference, but not the read
                    k += 1 # increment the alignment
                else:
                    # substitution
                    codon = [str(seq[i:i+3])]
                    assert aa == Seq.translate(codon[0])
                    i += 3 # increment the read
                    ii = 0 # tracks the position in an insertion
                    while (k+ii+1) < n and hsp.aln[0].seq[k+ii+1] == ".":
                        # insertion
                        aa = hsp.aln[1].seq[k+ii+1]
                        assert aa.islower() or aa == "*"
                        codon.append(str(seq[i:i+3]))
                        assert aa.upper() == Seq.translate(codon[-1])
                        i += 3 # increment the read
                        ii += 1
                    codon = "".join(codon)
                    counts[hsp.query_start+j][codon] = counts[hsp.query_start+j].get(codon, 0) + count
                    j += 1 # increment the reference
                    k += ii + 1 # increment the alignment

    # output
    with open(outfile, "w") as f:
        print("pos", "codon", "count", sep="\t", file=f)
        for i, count in enumerate(counts):
            for codon in sorted(count, key=count.get, reverse=True):
                print(i, codon, count[codon], sep="\t", file=f)

    return True

# vim: expandtab sw=4 ts=4
