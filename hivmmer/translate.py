import argparse
import sys
from Bio import Seq
from Bio import SeqIO


def translate(filename, allow_stop_codons, out=sys.stdout, log=sys.stderr):
    """
    Translate nucleotide sequences in FASTA file `filename` to all six possible
    frames, excluding any that have > `allow_stop_codons` stop codons.

    Write amino acid sequences to FASTA file `out`, with the frame number
    appended to the sequence header.

    Log summary statistics to file `log`.
    """

    nskipped = {"N": 0, "*": 0}

    for n, record in enumerate(SeqIO.parse(fasta, "fasta")):

        nstop = 0
        seq = str(record.seq)

        if 'N' in seq:
            nskipped["N"] += 1
            continue

        for i in range(3):
            aa = Seq.translate(seq[i:])
            if aa.count("*") > 0:
                nstop += 1
            else:
                print(">%s-%d" % (record.id, i), file=out)
                print(aa, file=out)

        seq = str(record.seq.reverse_complement())

        for i in range(3):
            aa = Seq.translate(seq[i:])
            if aa.count("*") > allow_stop_codons:
                nstop += 1
            else:
                print(">%s-%d'" % (record.id, i), file=out)
                print(aa, file=out)

        if nstop == 6:
            nskipped["*"] += 1

    print("nreads", n, file=log)
    print("nskipped N", nskipped["N"], file=log)
    print("nskipped *", nskipped["*"], file=log)


def _run():

    parser = argparse.ArgumentParser()
    parser.add_argument("FASTA",
                        nargs=1,
                        help="path to nucleotide FASTA file, or '-' to read from stdin")
    parser.add_argument("--allow-stop-codons",
                        default=0,
                        metavar="N",
                        type=int,
                        help="# of stop codons to allow [default: 0]")
    args = parser.parse_args()

    fasta = args.FASTA[0]
    if fasta == '-':
        translate(sys.stdin, args.allow_stop_codons)
    else:
        with open(fasta) as f:
            translate(f, args.allow_stop_codons)


# vim: expandtab sw=4 ts=4
