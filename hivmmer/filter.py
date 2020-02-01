import argparse
import sys
from Bio import SeqIO
from operator import itemgetter


def add(filename, min_length, min_quality, keep={}):
    """
    Filters out Illumina reads in input FASTA file `filename` with
    length < `min_length` and mean quality score < `min_quality`.

    Deduplicates identical sequences.

    Adds distinct sequences to dictionary `keep`, with the sequence
    count as the value.
    """
    for seq in SeqIO.parse(filename, "fastq"):
        if len(seq) > min_length and sum(seq.letter_annotations["phred_quality"]) / len(seq) > min_quality:
            seq = str(seq.seq)
            keep[seq] = keep.get(seq, 0) + 1


def tofasta(keep, f):
    """
    Writes distinct sequences/counts in dictionary `keep` to FASTA file `f`.
    """
    for i, (seq, n) in enumerate(sorted(keep.items(), key=itemgetter(1), reverse=True)):
        print(">{}-{}".format(i, n), file=f)
        print(seq, file=f)


def _run():

    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--min-length",
                        default=75,
                        metavar="L",
                        type=int,
                        help="minimum read length to retain")
    parser.add_argument("-q", "--min-quality",
                        default=25,
                        metavar="Q",
                        type=float,
                        help="minimum mean quality score to retain")
    parser.add_argument("FASTQ",
                        nargs="+",
                        help="list of FASTQ input files to filter")
    args= parser.parse_args()

    keep = {}

    for fastq in args.FASTQ:
        add(fastq, args.min_length, args.min_quality, keep)

    tofasta(keep, sys.stdout)


# vim: expandtab sw=4 ts=4
