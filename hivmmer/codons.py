import argparse
import numpy as np
import os
import sys
from Bio import SearchIO
from Bio import Seq
from Bio import SeqIO
from datetime import datetime
from itertools import product

nt = ("A", "C", "G", "T")

aa_lookup = {
    "A": ("GCT", "GCC", "GCA", "GCG"),
    "C": ("TGT", "TGC"),
    "D": ("GAT", "GAC"),
    "E": ("GAA", "GAG"),
    "F": ("TTT", "TTC"),
    "G": ("GGT", "GGC", "GGA", "GGG"),
    "H": ("CAT", "CAC"),
    "I": ("ATT", "ATC", "ATA"),
    "K": ("AAA", "AAG"),
    "L": ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    "M": ("ATG",),
    "N": ("AAT", "AAC"),
    "P": ("CCT", "CCC", "CCA", "CCG"),
    "Q": ("CAA", "CAG"),
    "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    "S": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    "T": ("ACT", "ACC", "ACA", "ACG"),
    "V": ("GTT", "GTC", "GTA", "GTG"),
    "W": ("TGG",),
    "Y": ("TAT", "TAC"),
    "*": ("TAA", "TGA", "TAG")
}

codons = list(map("".join, product(nt, repeat=3))) + ["ins", "del"]
codon_index = dict((codon, i) for (i, codon) in enumerate(codons))


def extract_codons(reads, hmmer):

    nref = max(f.query_end for f in hmmer.fragments)

    table = np.zeros((nref, len(codons)), dtype=np.uint32)
    txt = [{} for _ in range(nref)]

    n_hsps = {}

    for hit in hmmer.hits:

        id, _, frame = hit.id.rpartition("-")
        count = int(id.partition("-")[2])

        if frame.endswith("'"):
            seq = reads[id].seq.reverse_complement()
            offset = int(frame[:-1])
        else:
            seq = reads[id].seq
            offset = int(frame)

        # Track number of hsps in each hit
        n_hsps[len(hit.hsps)] = n_hsps.get(len(hit.hsps), 0) + 1

        for hsp in hit.hsps:

            if (hsp.bitscore / hsp.hit_span) < min_score: continue

            i = 3*hsp.hit_start + offset # tracks position in the read sequence
            j = 0                                                # tracks position in the reference sequence
            k = 0                                                # tracks position in the hmmer alignment

            n = len(hsp.aln[0].seq)
            assert n == len(hsp.aln[1].seq)

            # alignment should not start with an insertion
            assert hsp.aln[0].seq[0] != "."

            while k < n:
                aa = hsp.aln[1].seq[k]
                if aa == "-":
                    # deletion
                    table[hsp.query_start+j, codon_index["del"]] += count
                    txt[hsp.query_start+j][""] = txt[hsp.query_start+j].get("", 0) + count
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
                    if len(codon) == 1:
                        table[hsp.query_start+j, codon_index[codon[0]]] += count
                    else:
                        table[hsp.query_start+j, codon_index["ins"]] += count
                    codon = "".join(codon)
                    txt[hsp.query_start+j][codon] = txt[hsp.query_start+j].get(codon, 0) + count
                    j += 1 # increment the reference
                    k += ii + 1 # increment the alignment

    print("hsps:", n_hsps, file=sys.stderr)

    return table, txt

def to_txt(codons, f):


def to_aavf(table, ref, out, region="pol"):
    """
    """
    refseq = next(SeqIO.parse(ref, "fasta"))
    chrom = refseq.id
    refseq = str(refseq.seq)

    # Metadata
    print("""\
##fileformat=AAVFv1.0
##fileDate={date}
##source={source}
##reference={ref}
##INFO=<ID=AC,Number=.,Type=String,Description="Alternate Codon">
##INFO=<ID=ACC,Number=.,Type=Float,Description="Alternate Codon Count, for each Alternate Codon, in the same order as listed.">
##INFO=<ID=ACF,Number=.,Type=Float,Description="Alternate Codon Frequency, for each Alternate Codon, in the same order as listed.">
##FILTER=<ID=cov10,Description="Alternate Count < 10">
##FILTER=<ID=freq0.01,Description="Alternate Frequency < 0.01">""".format(
            date=datetime.now().strftime("%Y%m%d %H:%M:%S"),
            source="hivmmer-{}".format(version),
            ref=ref),
        file=out)

    # Header
    print("#CHROM", "GENE", "POS", "REF", "ALT",
                "FILTER", "ALT_FREQ", "COVERAGE", "INFO",
                sep="\t", file=out)

    # Rows
    for i, row in enumerate(table):
        cov = row.sum()
        norm = 1.0 / cov
        if region == "pol":
            if i >= 99:
                gene = "RT"
                pos = i - 98
            else:
                gene = "PR"
                pos = i + 1
        else:
            gene = region[:2].upper()
            pos = i + 1
        for aa, codonlist in aa_lookup.items():
            ac = []
            acc = []
            acf = []
            for codon in codonlist:
                count = table[i, codon_index[codon]]
                if (count > 0):
                    ac.append(codon.lower())
                    acc.append(count)
                    acf.append(count * norm)
            if len(ac) > 0:
                count = sum(acc)
                freq = count * norm
                if count < 10:
                        filt = "cov10"
                elif freq < 0.01:
                        filt = "freq0.01"
                else:
                        filt = "PASS"
                print(chrom, gene, pos, refseq[i], aa, filt,
                            freq, cov, "AC={};ACC={};ACF={}".format(",".join(ac),
                                                                                                            ",".join(map(str, acc)),
                                                                                                            ",".join(map(str, acf))),
                            sep="\t", file=out)


def _run():

    parser = argparse.ArgumentParser()
    parser.add_argument("--hmmer",
                        required=True,
                        help="path to HMMER text output")
    parser.add_argument("--reads",
                        required=True,
                        help="path to amino acid FASTA file aligned with HMMER")
    parser.add_argument("--ref",
                                            required=True,
                                            help="path to the reference NT sequence (e.g. hxb2)")
    parser.add_argument("--region",
                                            choices=("pol", "int", "env"),
                                            default="pol",
                                            help="region of the HIV genome (for indexing)")
    args = parser.parse_args()

    reads = SeqIO.index(args.reads, "fasta")
    print("indexed reads", file=sys.stderr)

    hmmer = SearchIO.read(args.hmmer, "hmmer3-text")
    print("parsed hmmer", file=sys.stderr)

    table, txt = codon_table(reads, hmmer)
    print("constructed codon table", file=sys.stderr)

    basename = os.path.splitext(args.hmmer)[0]

    # CSV output
    np.savetxt(basename + ".codons.csv", table, fmt="%g", comments="",
                         delimiter=",", header=",".join(codons))

    # TXT output
    with open(basename + ".codons.txt", "w") as f:
        print("pos", "codon", "freq", sep="\t", file=f)
        for i, freq in enumerate(txt):
            for codon in sorted(freq, key=freq.get, reverse=True):
                print(i, codon, freq[codon], sep="\t", file=f)

    # AAVF output
    with open(basename + ".aavf", "w") as f:
        print_aavf(table, args.ref, f, args.region)

# vim: expandtab sw=4 ts=4
