#!/bin/sh -x
# -*- mode: sh; -*-

# silly script to split bam files into separate {R1,R2}{mapped,unmapped} bam file

bam="$@"

name=$(basename $bam .bam)

## R1: -f 64
## R2: -f 128
## unmapped: -f 4
## mapped: -F 4

samtools view -h  -f 64 -F 4   $bam -b > $name-R1mapped.bam
samtools view -h  -f 64 -f 4   $bam -b > $name-R1unmapped.bam
samtools view -h  -f 128 -F 4   $bam -b > $name-R2mapped.bam
samtools view -h  -f 128 -f 4   $bam -b > $name-R2unmapped.bam
