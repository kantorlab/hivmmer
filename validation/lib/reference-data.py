import os
import sys
from Bio import SeqIO
from numpy.random import choice, randint, seed

seed(81543728)

reads_file = sys.argv[1]
ref_files = sys.argv[2:]

refs = {}
for ref_file in ref_files:
    gene = os.path.basename(ref_file).partition(".")[0]
    refs[gene] = []
    for seq in SeqIO.parse(open(ref_file), "fasta"):
        refs[gene].append(str(seq.seq).replace("-", ""))
    print(len(refs[gene]), gene, "sequences", file=sys.stderr)

genes = list(refs.keys())
lengths = [sum(map(len, refs[gene])) / len(refs[gene]) for gene in genes]
for gene, avglen in zip(genes, lengths):
    print(gene, 3*avglen, file=sys.stderr)
sumlen = sum(lengths)
print("all", 3*sumlen, file=sys.stderr)
p = [x/sumlen for x in lengths]

for i, seq in enumerate(SeqIO.parse(reads_file, "fasta")):
    seqlen = len(seq.seq) // 3
    gene = choice(genes, p=p)
    refnum = randint(0, len(refs[gene]))
    start = len(refs[gene][refnum]) - seqlen
    if start > 0:
        start = randint(0, start)
    print(">{}-{}".format(gene, i))
    print(refs[gene][refnum][start:start+seqlen])

