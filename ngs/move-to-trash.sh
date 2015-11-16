#!/bin/sh
# -*- mode: sh; -*-

### thing to be used as a replacement for rm, e.g. in Makefiles

trash=./trash

if [ ! -d  $trash ]; then 
    if ! mkdir $trash; then 
      echo "$0: could not mkdir $trash"
      exit 2
    fi
fi

if ! \mv -f "$@" $trash/ ; then 
    echo "$0: could not move $@ to to $trash" 1>&2
    exit 3
fi
