#!/bin/sh
# -*- mode: sh; -*-

for bed in "$@"; do
  name=$(basename $bed .bed)
  (echo "track type=bedGraph name='$name'";  
   awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$5}' < $bed |\
    sort -k1,1 -k2,2n ) > "$name.bedGraph"
done
