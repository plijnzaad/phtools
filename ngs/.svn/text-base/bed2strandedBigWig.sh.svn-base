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

tawk () {
    awk -F "\t" -v OFS="\t" "$@"
}

declare -a strand_symbol strand_name
strand_symbols=('+' '-')
strand_names=('fwd' 'rev')

for bed in "$@"; do 
    if [[ ! $bed =~ \.bed$ ]]; then
        echo "$i: Expected a .bed suffix" >&2
        exit 4
    fi
    base=`basename $bed .bed`

    for i in 0 1; do 
        sym=${strand_symbols[$i]}
        strand=${strand_names[$i]}
        bedgraph=$strand,$base.bedGraph
        bigwig=$strand,$base.bw
        tawk '$6 == "'$sym'"' $bed \
            | sort -k1,1 \
            | bedItemOverlapCount dummydatabase -chromSize=$chromsizes stdin \
            | sort -k1,1 -k2,2n > $bedgraph
        bedGraphToBigWig $bedgraph $chromsizes $bigwig
        rm $bedgraph
        echo "Created $bigwig"
    done
done

