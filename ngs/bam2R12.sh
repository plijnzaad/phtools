#!/bin/sh -x
# -*- mode: sh; -*-

# silly script to split bam files into separate {R1,R2}{mapped,unmapped} bam file

bam="$@"

base=$(basename $bam .bam)

## R1: -f 64
## R2: -f 128
## unmapped: -f 4
## mapped: -F 4

resources="-l h_rt=0:30:00,h_vmem=1G"

name="$base-R1mapped";   qrun.sh -N $name  $resources "samtools view -h  -f 64 -F 4   $bam -b > $name.bam"
name="$base-R1unmapped"; qrun.sh -N $name  $resources "samtools view -h  -f 64 -f 4   $bam -b > $name.bam"
name="$base-R2mapped";   qrun.sh -N $name  $resources "samtools view -h  -f 128 -F 4   $bam -b > $name.bam"
name="$base-R2unmapped"; qrun.sh -N $name  $resources "samtools view -h  -f 128 -f 4   $bam -b > $name.bam"
