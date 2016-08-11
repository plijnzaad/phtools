#!/bin/sh
# -*- mode: sh; -*-
## drop-in replacement of sambamba sort that makes use of the proper local TMPDIR on the compute node
##
## Usage in pipeline e.g.: 
##    somepipeline |  sambambasort.sh --nthreads=4 --memory-limit=2GB /dev/stdin -o stuff.bam
##
## (NOTE: needs a bam file, not a sam file!)
##
sambamba sort --tmpdir=$TMPDIR "$@"
