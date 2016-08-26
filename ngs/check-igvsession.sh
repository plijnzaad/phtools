#!/bin/sh

#### note: this is based on the xml file as produced by IGV 2.3.72
#### also: we rely on the string '/yeast' being present in the URL,
### otherwise we have to check too much ...
### Written by plijnzaad@gmail.com

file="$@"

urlnotfound() {
    local url=$1
    curl -s "$url" | strings | head -40  | egrep -iq 'error|not.*found'
    return $?
}

### Resources first:
sed -n '/Resource.*\/yeast/{s|<.*http|http|;s|".*$||;p;}' $file |\
    while read url; do
        echo $url
        if urlnotfound $url; then
            echo "**** Resource $url not found ****"
        fi
    done

sed -n '/Track.*\/yeast/{s/^.*http/http/;s/" name.*$//;p;}' $file |\
    while read url; do
        echo $url
        if urlnotfound $url; then
            echo "**** Track contains URL $url which is not found ****"
        fi
    done
