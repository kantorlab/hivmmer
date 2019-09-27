import sys
from Bio import SeqIO
from numpy.random import choice

nts = ["A", "C", "G", "T"]

for seq in SeqIO.parse(sys.argv[1], "fasta"):
    print(seq.id)
    print("".join(choice(nts, size=len(seq.seq)))

