#!/bin/bash -x

## Script to get rid of non-unique reads. Assuming bwa-mapped files, you'd typically
## keep anything with a mapping quality of 1 or higher, using samtools view -q 1

if [ $# -ne 3 ]; then
    echo "Usage: $0 directory bamfile qual"
    exit 1;
fi

dir="$1"
if ! cd "$dir"; then
    echo "Directory $dir not an existing accessible directory" >&2
    exit 2
fi

bamfile="$2"
if ! [ -f "$bamfile" ]; then
    echo "File $bamfile not an existing BAM file" >&2
    exit 3
fi

if ! samtools view -H "$bamfile" | grep -l '^@HD.*coordinate'; then
    echo "bamfile $bamfile is not sorted"
    exit 4
fi

qual=$3
if [[ ! $qual =~ ^[0-9]+$ ]]; then
    echo "$qual not a numeric quality value" >&2
    exit 5
fi


name=$(basename $bamfile .bam)
new="$name-q$qual"

exec > $new.log 2>&1

nice samtools view -q $qual -h $bamfile -b > $new.bam
nice samtools index $new.bam

echo "Produced  $new.bam and $new.bam.bai"
exit 0

