import sys
from Bio import SearchIO
from io import StringIO
from subprocess import run, PIPE

hxb2file, hmmfile, gene, offset, outfile = sys.argv[1:]

# setup the DRM lookup table for pol
drm = {}
drmgene = {}
if gene == "pol":
  for i, hxb2 in enumerate(range(2253, 2550, 3)):
    drm[hxb2] = i+1
    drmgene[hxb2] = "PR"
  for i, hxb2 in enumerate(range(2550, 3870, 3)):
    drm[hxb2] = i+1
    drmgene[hxb2] = "RT"
  for i, hxb2 in enumerate(range(4230, 5096, 3)):
    drm[hxb2] = i+1
    drmgene[hxb2] = "IN"

# run hmmscan
hmmer = run(["hmmscan", "--max", hmmfile, hxb2file], stdout=PIPE, check=True, encoding="utf8")

# parse hmmscan output
print(hmmer.stdout)
hmmer = SearchIO.read(StringIO(hmmer.stdout), "hmmer3-text")
hsp = hmmer.hits[0].hsps[0]
scores = hsp.aln_annotation["PP"]

with open(outfile, "w") as f:

  # Header
  print("hmm", "hxb2", "drm", "drmgene", "ins", "del", "hmmaa", "hxb2aa", "score", sep="\t", file=f)

  # Coordinate starts
  hmm  = hsp.hit_start + 1
  hxb2 = int(offset)

  # Can't start with an insertion or deletion
  assert hsp.aln[1].seq[0] != "."
  assert hsp.aln[0].seq[0] != "-"

  i = 0
  while i < len(scores):

    hmm_aa  = hsp.aln[1].seq[i]
    hxb2_aa = hsp.aln[0].seq[i]
    score   = scores[i]

    # Find insertions relative to HMM coordinates
    ins = 0
    while i < (len(scores) - 1) and hsp.aln[1].seq[i+1] == ".":
      assert hsp.aln[0].seq[i+1].islower()
      hxb2_aa += hsp.aln[0].seq[i+1]
      score   += scores[i]
      ins     += 1
      i       += 1

    # Find deletions relative to HMM coordinates
    de = 0
    while i < (len(scores) - 1) and hsp.aln[0].seq[i+1] == "-":
      hmm_aa += hsp.aln[1].seq[i+1]
      score  += scores[i]
      de     += 1
      i      += 1

    print(hmm,
          hxb2,
          drm.get(hxb2, ""),
          drmgene.get(hxb2, ""),
          ins,
          de,
          hmm_aa,
          hxb2_aa,
          score,
          sep="\t", file=f)

    # Reprint at each hmm coordinate in a deletion
    for j in range(1, de+1):
        print(hmm+j,
              hxb2,
              drm.get(hxb2, ""),
              drmgene.get(hxb2, ""),
              ins,
              de-j,
              hmm_aa,
              hxb2_aa,
              score,
              sep="\t", file=f)

    i    += 1
    hmm  += 1 + de
    hxb2 += 3*(1 + ins)

# vim: expandtab sw=2 ts=2
