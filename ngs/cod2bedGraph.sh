#!/bin/sh -x 

# Note: you may want to filter on FDR; this is column $7 in cod file

for i in "$@"; do 
  name=$(echo $i | sed 's/_peak.cod//')
  ( echo "track type=bedGraph name='$name'"; grep -v '^#rank' $i \
    | awk -F "\t" -v OFS="\t" '{print $2,$3,$4, $15;}' ) \
    | sort -k1,1 -k2n > $(basename $i .cod).bedGraph
done

