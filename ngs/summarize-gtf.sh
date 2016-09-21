#!/bin/sh
# -*- mode: sh; -*-

# silly thing to count stuff in gtf file. Usage: summarize.gtf.sh somefile.gtf

echo -n "chromosomes: "
cut -f 1  $1 | sort -u  | wc -l 
echo -n "genes: "
cat $1 | sed 's/.*gene_id/gene_id/g; s/;.*//' | sort -u | wc -l
echo -n "transcripts: "
cat $1 | sed 's/.*transcript_id/transcript_id/g; s/;.*//' |sort -u | wc -l
echo -n "exons: "
grep "	exon	" $1 |wc -l 
