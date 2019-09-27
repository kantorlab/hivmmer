#!/bin/bash
set -e
FQ1=$1
FQ2=$2
DIR=$3
pear -y 3G -f $FQ1 -r $FQ2 -o $DIR/pear -k
cat $DIR/pear.*.fastq | fastx_collapser -Q33 -o $DIR/collapsed.fa
