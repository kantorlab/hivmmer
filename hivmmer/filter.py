import argparse
import sys
from Bio import SeqIO
from operator import itemgetter


def mean_filter(filename, min_length, min_quality, keep={}):
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


def split_filter(filename, min_length, min_quality, keep={}):
    """
    Split reads in input FASTA file `filename` into subsequences with
    quality score > `min_quality` and length > `min_length`.

    Deduplicates identical subsequences.

    Adds distinct subsequences to dictionary `keep`, with the subsequence
    count as the value.
    """
    for seq in SeqIO.parse(filename, "fastq"):
        quality = seq.letter_annotations["phred_quality"]
        seq = str(seq.seq)
        start = 0
        i = 0
        for i, (nt, q) in enumerate(zip(seq, quality)):
            if nt == 'N' or q < min_quality:
                subseq = seq[start:i+1]
                if len(subseq) > min_length:
                    keep[subseq] = keep.get(subseq, 0) + 1
                start = i


def mask_filter(filename, min_quality, keep={}):
    """
    Mask bases as N in sequences from `filename` if quality score
    is < `min_quality`.

    Deduplicates identical subsequences.

    Adds distinct subsequences to dictionary `keep`, with the subsequence
    count as the value.
    """
    for seq in SeqIO.parse(filename, "fastq"):
        quality = seq.letter_annotations["phred_quality"]
        masked = []
        for nt, q in zip(str(seq.seq), quality):
            if q >= min_quality:
                masked.append(nt)
            else:
                masked.append('N')
        masked = "".join(masked)
        keep[masked] = keep.get(masked, 0) + 1


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
                        help="minimum quality score to retain")
    parser.add_argument("-m", "--mode",
                        choices=["mean", "split"],
                        default="mean",
                        help="filter on mean quality score vs. split into high-quality subsequences")
    parser.add_argument("FASTQ",
                        nargs="+",
                        help="list of FASTQ input files to filter")
    args= parser.parse_args()

    keep = {}
    modes = {"mean": mean_filter, "split": split_filter}

    for fastq in args.FASTQ:
        modes[args.mode](fastq, args.min_length, args.min_quality, keep)

    tofasta(keep, sys.stdout)


# vim: expandtab sw=4 ts=4
