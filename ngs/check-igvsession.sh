#!/bin/bash
### Check validity of URLs in side an igv_session.xml file.
### note: this is based on the xml file as produced by IGV 2.3.72
### Written by plijnzaad@gmail.com

file="$@"

urlnotfound() {
    local url=$1
    curl -s "$url" | strings | head -40  | egrep -iq 'error|not.*found'
    return $?
}

### NOTE: we canot keep track of totals, since the while-construct below
### runs in a sub-shell ...

### Resources first:
echo "==== Resources ===="
sed -n '/Resource.*http/{s|<.*http|http|;s|".*$||;p;}' $file |\
    while read url; do
        echo -n $url
        if urlnotfound $url; then
            echo "**** NOT FOUND ****"
        else
            echo " OK"
        fi
    done

echo "==== Tracks ===="
sed -n '/Track.*http/{s/^.*http/http/;s/" name.*$//;p;}' $file |\
    while read url; do
        echo -n $url
        if urlnotfound $url; then
            echo "**** NOT FOUND ****"
        else
            echo " OK"
        fi
    done
