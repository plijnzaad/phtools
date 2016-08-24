#!/bin/sh
# -*- mode: sh; -*-
## This script is just so that the queueing system can supply the $TMPDIR
## Typically used as e.g. 
## 
##  ... | samtoolssort.sh -@ 4  -n -m 8G  -Obam -o output/x.bam
## 
tmpdir=$TMPDIR/sort$RAND
mkdir -p $tmpdir || exit 28
samtools sort -T$tmpdir "$@"
