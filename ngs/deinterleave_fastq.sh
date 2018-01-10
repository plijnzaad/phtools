#!/bin/bash
# Usage: see below
# 
# Deinterleaves a FASTQ file of paired reads into two FASTQ
# files specified on the command line. Optionally GZip compresses the output
# FASTQ files using pigz if the 3rd command line argument is the word "compress"
# 
# Can deinterleave 100 million paired reads (200 million total
# reads; a 43Gbyte file), in memory (/dev/shm), in 4m15s (255s)
# 
# Adapted from https://gist.github.com/3521724
# also compress the 4th lines to just '+'
# (Also see his interleaving script: https://gist.github.com/4544979)
# 
# Inspired by Torsten Seemann's blog post:
# http://thegenomefactory.blogspot.com.au/2012/05/cool-use-of-unix-paste-with-ngs.html

if [ $# -ne 1 ] ;then
  echo "Usage: $0 OUTNAME < INPUT.fastq   # output goes to OUTNAME_R1.fastq.gz and OUTNAME_R2.fastq.gz" 1>&2
  exit 2
fi

paste - - - - - - - -  | tee >(
  cut -f 1-4 | tr "\t" "\n" | sed 's/^+.*/+/' | gzip > $1_R1.fastq.gz) |\
  cut -f 5-8 | tr "\t" "\n" | sed 's/^+.*/+/' | gzip > $1_R2.fastq.gz

