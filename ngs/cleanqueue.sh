#!/bin/sh

echo "Untested!" >&2
exit 7

hosts=`qhost |cut -f1 -d " "|grep -P "n\d{4}"`

for host in $hosts; do 
    qsub -shell no -b yes -l hostname=$host sh -c 'find /tmp/udcCache -depth -user philip -print -exec rm -fr {} \;'
    fi
done

