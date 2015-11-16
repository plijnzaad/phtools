#!/bin/bash -x

### recombine chromosome-specific bigwig (.bw) files back into complete genomes
### (simple converts to bedGraph, concatenates, then converts back to bigwig)
### This is done for MNase data where peaks have been shifted to get a better nucleosome
### binding peaks


if [ $# -lt 2 ]; then
    echo "Usage: $0 *.bw" >&2
    exit 3
fi

if [ ! -f "$chromsizes" ]; then
    echo "Could not find file  '$chromsizes'" >&2
    exit 3
fi

scriptname=`basename $0`

### tmpdir=` mktemp -d --suffix=$scriptname `
tmpdir=` mktemp -d -t `

trap "rm -fr $tmpdir; exit 2"  SIGINT SIGQUIT SIGTERM SIGHUP SIGKILL SIGABRT SIGCHLD SIGSEGV SIGSYS SIGILL SIGIO SIGALRM


for i in "$@"; do 
    if [[ ! $i =~ ^chr[IXVm][IXVt]*,.*bw$ ]]; then
        echo "$i: Expected a chromosome prefix (e.g. chrXI,) and .bw suffix" >&2
        exit 4
    fi
    base=`basename $i .bw`
    nice bigWigToBedGraph $i $tmpdir/$base.bedGraph
done

# invent new name:

allname=$(echo "$1" | sed 's/^chr[IXVm][IXVt]*,/ALLchr/' | sed 's/\.bw$//')

cp /dev/null $tmpdir/$allname.bedGraph
for i in "$@"; do 
  base=`basename $i .bw`
  cat $tmpdir/$base.bedGraph >> $tmpdir/$allname.bedGraph
done

if nice bedGraphToBigWig -- $tmpdir/$allname.bedGraph $chromsizes $allname.bw; then
  echo "Succesfully created $allname.bw"
fi

rm -fr $tmpdir

## rm "$@"  # throw away originals?
