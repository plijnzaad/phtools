#!/bin/sh -x
# -*- mode: sh; -*-

echo "Use bam2fastq, much faster and less memory hungry. Exiting now." >&2
exit 99

if [ $# -eq 0 ]; then
    echo "Usage: $0 [ --paired_end ] input.[sb]am [ output.fastq ]" >&2
    exit 2
fi

if [ $1 = '--paired_end' ]; then
    paired_end=yes
    shift
fi

input="$1"

if [ $# -eq 2 ]; then
    output="$2" # uh ... what about paired end?
else
    output=$(echo "$input" | sed 's/\.[bs]am$//g')
    if [ x$paired_end = xyes ]; then
        output1="${output}_R1.fastq"
        output2="${output}_R2.fastq"
    else
        output="$output.fastq"
    fi
fi

opts="VALIDATION_STRINGENCY=SILENT"

if [ x$paired_end = xyes ]; then
    ptrun.sh SamToFastq $opts INPUT="$input" FASTQ="$output1" SECOND_END_FASTQ="$output2"
else
    ptrun.sh SamToFastq $opts INPUT="$input" FASTQ="$output"
fi

### NOTE: also check bedtools bamtofastq
