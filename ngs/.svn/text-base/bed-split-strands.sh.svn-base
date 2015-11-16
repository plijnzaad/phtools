#!/bin/sh

if [ $# -lt 1 ]; then
    echo "Usage: $0 *.bed" >&2
    exit 3
fi

if [ ! -f "$chromsizes" ]; then
    echo "Could not find file  '$chromsizes'" >&2
    exit 3
fi

scriptname=`basename $0`

declare -a strand_symbol strand_name
strand_symbols=('+' '-')
strand_names=('fwd' 'rev')

for bed in "$@"; do 
    if [[ ! $bed =~ \.bed$ ]]; then
        echo "$i: Expected a .bed suffix" >&2
        exit 4
    fi
    name=$(basename $bed .bed)
    for i in 0 1; do 
        sym=${strand_symbols[$i]}
        strand=${strand_names[$i]}
        splitbed="$name-strand.bed"
        nice awk -F "\t" -v OFS="\t" '$6 == "'$sym'"' $bed > $splitbed
        echo "Created $splitbed"
    done
done

