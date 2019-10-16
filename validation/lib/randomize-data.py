import sys
from Bio import SeqIO
from numpy.random import choice, seed

seed(87197155)

nts = ["A", "C", "G", "T"]

for seq in SeqIO.parse(sys.argv[1], "fasta"):
    print(">{}".format(seq.id))
    print("".join(choice(nts, size=len(seq.seq))))

