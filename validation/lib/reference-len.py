import sys
from Bio import SeqIO
n = 0
for record in SeqIO.parse(sys.argv[1], "fasta"):
    n += len(str(record.seq).replace("-", ""))
print(n)
