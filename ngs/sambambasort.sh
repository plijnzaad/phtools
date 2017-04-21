#!/bin/sh
# -*- mode: sh; -*-
## drop-in replacement of samtools sort (but using sambamba for speed)
## that makes use of the proper local TMPDIR on the compute node. See also 
## the Makefile in this directory.
##
## Usage in pipeline e.g.: 
##    somepipeline |  sambambasort.sh --nthreads=4 --memory-limit=2GB /dev/stdin -o stuff.bam
##
## (NOTE: needs a bam file, not a sam file!)
##
sambamba sort --tmpdir=$TMPDIR "$@"
