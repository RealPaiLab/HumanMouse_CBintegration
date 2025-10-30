#!/bin/bash

inFile=/u/spai/data/CBL-dev/exprMatrix.tsv.gz
ncols=`zcat $inFile | awk '{print NF}' | sort -nu | tail -n 1`
echo "Num cols: $ncols"
nrows=`zcat $inFile | wc -l`
echo "Num rows: $nrows"
