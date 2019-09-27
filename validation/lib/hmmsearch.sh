#!/bin/bash
set -e
ID=$1
GENE=$2
REF=$3
PFA=$4
CPU=$5

PREFIX=scratch/hmmsearch1/$ID/$GENE

hmmsearch --cpu $CPU -E 1000000 --domE 1000000 --max --notextw -A ${PREFIX}.sto -o ${PREFIX}.txt "$REF" "$PFA"
hmmemit -c "$REF" > ${PREFIX}.ref.consensus.fa
hmmbuild --cpu $CPU ${PREFIX}.hmm ${PREFIX}.sto
