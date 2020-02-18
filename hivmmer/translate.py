import argparse
import sys
from Bio import Seq
from Bio import SeqIO


def translate(filename, out=sys.stdout, log=sys.stderr):
    """
    Translate nucleotide sequences in FASTA file `filename` to all six possible
    frames.

    Write amino acid sequences to FASTA file `out`, with the frame number
    appended to the sequence header.

    Log summary statistics to file `log`.
    """

    nskipped = 0

    for n, record in enumerate(SeqIO.parse(filename, "fasta")):

        seq = str(record.seq)

        if 'N' in seq:
            nskipped += 1
            continue

        for i in range(3):
            j = 3 * ((len(seq) - i) // 3) + i
            print(">%s-%d" % (record.id, i), file=out)
            print(Seq.translate(seq[i:j]), file=out)

        seq = str(record.seq.reverse_complement())

        for i in range(3):
            j = 3 * ((len(seq) - i) // 3) + i
            print(">%s-%d'" % (record.id, i), file=out)
            print(Seq.translate(seq[i:j]), file=out)

    print("nreads", n, file=log)
    print("nskipped (N)", nskipped, file=log)


def translate_unambiguous(filename, out=sys.stdout, log=sys.stderr, fraction=0.5):
    """
    Translate nucleotide sequences in FASTA file `filename` to all six possible
    frames.

    Only keep translated sequences with > `fraction` of unambigous amino acids.

    Write amino acid sequences to FASTA file `out`, with the frame number
    appended to the sequence header.

    Log summary statistics to file `log`.
    """

    nskipped = 0

    for n, record in enumerate(SeqIO.parse(filename, "fasta")):

        seq = str(record.seq)

        for i in range(3):
            j = 3 * ((len(seq) - i) // 3) + i
            tseq = Seq.translate(seq[i:j])
            if sum(aa != 'X' for aa in tseq) / len(tseq) > fraction:
                print(">%s-%d" % (record.id, i), file=out)
                print(tseq, file=out)

        seq = str(record.seq.reverse_complement())

        for i in range(3):
            j = 3 * ((len(seq) - i) // 3) + i
            tseq = Seq.translate(seq[i:j])
            if sum(aa != 'X' for aa in tseq) / len(tseq) > fraction:
                print(">%s-%d'" % (record.id, i), file=out)
                print(tseq, file=out)

    print("nreads", n, file=log)


def _run():

    parser = argparse.ArgumentParser()
    parser.add_argument("FASTA",
                        nargs=1,
                        help="path to nucleotide FASTA file, or '-' to read from stdin")
    args = parser.parse_args()

    fasta = args.FASTA[0]
    if fasta == '-':
        translate(sys.stdin, args.allow_stop_codons)
    else:
        with open(fasta) as f:
            translate(f)


# vim: expandtab sw=4 ts=4
