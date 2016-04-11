#!/bin/sh
# -*- mode: sh; -*-
## This script is just so that the queueing system can supply the $TMPDIR
## Typically used as  
## 
##  ... | samtools sort -@ 4  -Obam > x.bam
## 
samtools sort -T$TMPDIR/sort "$@"
