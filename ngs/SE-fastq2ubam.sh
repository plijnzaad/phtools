#!/bin/bash
# -*- mode: sh; -*-
#
# Convert one fastq.gz file to ubam (written to stdout)
# 
# Based on ~/git/sharq2/miscUtilities/fastq2ubam.sh, but for only one mate
#

module load picardtools/2.20.1

## valid_vars="LIBRARY_NAME PLATFORM SEQUENCING_CENTER READ_GROUP_NAME"
valid_vars="
COMMENT
DESCRIPTION
LIBRARY_NAME
PLATFORM
PLATFORM_MODEL
READ_GROUP_NAME
RUN_DATE
SEQUENCING_CENTER
"
## aka CO DS LB PL PM PU PI PG RG DT CN
## (Kemmeren group uses only library_name, platform_unit, platform, seqcenter and readgroup)

if [ $# -ne 3 ]; then
    usage="\n\
Usage: [ var=val ... ] fastq2ubam.sh SAMPLE LANE FILE-R1.fastq.gz   | ... \n\
\n\
var=var clauses are passed to picard's FastqToSam. Valid variables are to to be included
in the header:\n\
$valid_vars\n\
\n\
SAMPLE_NAME and PLATFORM_UNIT are determined from the input filename. OUTPUT is always stdout.\n
"
    echo -e $usage >&2
    exit 1
fi

sample="$1"
lane="$2"
read1="$3"

TMPDIR=${TMPDIR-/var/tmp}

set +eu
other_vars=""

for var in $valid_vars; do
    if [ ${!var-undefined} != "undefined" ]; then
        other_vars="$other_vars -$var ${!var}"
    fi
done

## validation_stringency="VALIDATION_STRINGENCY=SILENT" # STRICT, LENIENT or SILENT
validation_stringency="-VALIDATION_STRINGENCY STRICT"
        
set -eu

nice java -Dpicard.useLegacyParser=false \
  -Djava.io.tmpdir=$TMPDIR -Xmx8G -jar $PICARD FastqToSam \
  $validation_stringency \
  -FASTQ $read1 \
  -OUTPUT /dev/stdout \
  -SAMPLE_NAME $sample \
  -READ_GROUP_NAME $lane \
  -PLATFORM_UNIT $lane \
  $other_vars

exit $?
### to turn into SAM format, pipe through ' | samtools view -h '
