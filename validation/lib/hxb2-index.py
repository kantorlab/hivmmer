import sys
from Bio import SearchIO
from io import StringIO
from subprocess import run, PIPE

hxb2file, hmmfile, offset, outfile = sys.argv[1:]

# run hmmscan
hmmer = run(["hmmscan", hmmfile, hxb2file], stdout=PIPE, check=True, encoding="utf8")

# parse hmmscan output
print(hmmer.stdout)
hmmer = SearchIO.read(StringIO(hmmer.stdout), "hmmer3-text")
assert len(hmmer.hits) == 1
assert len(hmmer.hits[0].hsps) == 1
hsp = hmmer.hits[0].hsps[0]
scores = hsp.aln_annotation["PP"]

with open(outfile, "w") as f:

  # Header
  print("hmm", "hxb2", "ins", "del", "hmm_aa", "hxb2_aa", "score", sep="\t", file=f)

  # Coordinate starts
  hmm  = 1
  hxb2 = int(offset)

  i = 0
  while i < len(scores):

    hmm_aa  = hsp.aln[1].seq[i]
    hxb2_aa = hsp.aln[0].seq[i]
    score   = scores[i]

    # Find insertions relative to HMM coordinates
    ins = 0
    while i < len(scores) - 1 and hsp.aln[1].seq[i+1] == ".":
      assert hsp.aln[0].seq[i+1].islower()
      hxb2_aa += hsp.aln[0].seq[i+1]
      score += scores[i]
      ins += 1
      i += 1

    print(hmm,
          hxb2,
          ins,
          Del,
          hmm_aa,
          hxb2_aa,
          score,
          sep="\t", file=f)

    hmm += 1
    hxb2 += 3*(1 + ins)

# vim: expandtab sw=2 ts=2
