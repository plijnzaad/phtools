#!/bin/sh

hosts=`seq -f 'n%04g' 1 200`
for host in $hosts; do 
    if ping -c 1 $host > /dev/null 2>&1 ;  then
        ssh $host find /tmp/udcCache -depth -user philip -print -exec rm -fr {} \\\;
    fi
done

