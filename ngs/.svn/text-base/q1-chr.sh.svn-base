#!/bin/bash

### takes directory and NAME.bam , and creates NAME-q1chr.bam
### containing unique sorted reads named as chrI, chrII etc.

echo Untested
exit 10

if [ $# -ne 2 ]; then
    echo "Usage: $0 directory bamfile"
    exit 1;
fi

dir="$1"
if ! cd "$dir"; then
    echo "Directory $dir not an existing accessible directory" >&2
    exit 2
fi

if ! [ -f "$bamfile" ]; then
    echo "File $bamfile not an existing BAM file" >&2
    exit 3
fi


name=$(basename $bamfile .bam)

new="$name-q1chr"

exec > $new.log 2>&1

addchr=sam-addchr.pl

nice samtools view -q1 -h $bamfile \
    | nice $addchr \
    | nice samtools view -h -S -b  \
    | nice samtools sort - $new.sorted

nice samtools index $new.sorted.bam

echo "Produced  $new.sorted.bam and $new.sorted.bam.bai"
exit 0

### mv $new.sorted.bam* quality-filtering/
