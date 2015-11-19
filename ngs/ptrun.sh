#!/bin/bash
# -*- mode: sh; -*-

### little wrapper script for the picard-tools (see
### http://sourceforge.net/projects/picard/)

### To get help on options specific to the subcommand, just 
### do ptrun.sh subcommand. For help on options common to 
### all subcommands, do ptrun subcommand --stdhelp. 

### Note that you can often use INPUT=/dev/stdin for reading from std input, and
### OUPUT=/devstdout QUIET=true for writing to stdout.

require_var() { 
    for var in "$@"; do
      eval "val=\$$var"
      if [ -z "$val"  ]; then
          echo "variable $var is not defined; exiting" >&2
          return 3
      fi
    done
}

require_var pthome
## e.g. /usr/local/picard-tools


pthelp() {
  pushd $pthome
  echo "
  Usage: ptrun.sh jarfile, where FILE is one of:
"
  ls -x *.jar | sed 's/\.jar//g'
  exit 1
}

if [ $# -lt 1 ] ; then
        pthelp
fi

jar="$1"
shift

java -d64 -Xmx512M -jar $pthome/$jar.jar "$@" VALIDATION_STRINGENCY=LENIENT

