#!/bin/sh -x
# -*- mode: sh; -*-

## Add 'chr' in front of seqnames (ensembl genome doesn't use them but allows
## them , UCSC genome requires them, and so do we. )

addchr=sam-addchr.pl

if [ $# -ne 2 ] ; then
    echo "Usage: $0 file.bame directory" >&2
    exit 2
fi

bamfile="$1"
outdir=`pwd`/"$2"

if [ ! -d $outdir ]; then
    if mkdir -p $outdir; then
        echo "Created dir $outdir ..."
    else
        echo "Failed creating dir $outdir" >&2
        exit 3
    fi
fi

nice samtools view -h $bamfile \
    | nice $addchr \
    | nice samtools view -h -S -b - \
    > $outdir/$bamfile

cd $outdir
nice samtools index $bamfile
exit 0

