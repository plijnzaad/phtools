#!/bin/sh
# -*- mode: sh; -*-
# stolen from dariober, https://www.biostars.org/p/15011/#103041

if [ $# -ne 1 ] ;then
  echo "Usage: $0 FILE.gz # FILE.gz is renamed FILE-unsorted.gz (and left there), sorted output to FILE.gz" >&2 
  exit 2
fi

buffersize=12G # or something

file=$1
name=${file/.gz}
unsorted="$name-unsorted.gz"

if [  -f $unsorted ] ; then
  echo "Not overwriting existing file $unsorted; nothing sorted." >&2 
  exit 5
fi 

if mv $file $unsorted ; then
: #  echo "Renamed $file to $unsorted" >&2   
else
  echo "Could not mv $file $unsorted" >&2 
  exit 3
fi

nice zcat $unsorted \
  | paste - - - - \
  | nice sort -k1,1 -S $buffersize \
  | tr '\t' '\n' \
  | awk '{ if(NR%4==1)sub(/^@[^ ][^ ]* */, "@"); if(NR%4==3)gsub(/.*/, "+"); print}' \
  | nice gzip > $file

## to do some on the fly editing on the fastq file (e.g. to get rid of useless SRR-numbers), consider
## sticking in an awk line just before the final 'nice gzip', e.g.
##  | awk '{ if(NR%4==1)sub(/^@[^ ][^ ]* */, "@"); if(NR%4==3)gsub(/.*/, "+"); print}' \


### To silently removed the original files:
## if rm $unsorted; then
## :
## else
##   echo "Could not remove unsorted file" >&2 
##   exit 4
## fi 

