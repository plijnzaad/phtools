#!/bin/sh
# -*- mode: sh; -*-
##
## drop-in replacement of samtools sort that makes use of the proper
## local TMPDIR on the compute node. See also the Makefile in this
## directory (and/or sambambasort.sh)

## This script is just so that the queueing system can supply the $TMPDIR
## Typically used as e.g. 
## 
##  ... | samtoolssort.sh -@ 4  -n -m 8G  -Obam -o output/x.bam
## 
tmpdir=$TMPDIR/sort$RAND
mkdir -p $tmpdir || exit 28
samtools sort -T$tmpdir "$@"
