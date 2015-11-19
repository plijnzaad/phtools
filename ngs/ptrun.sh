#!/bin/bash
# -*- mode: sh; -*-

### little wrapper script for the picard-tools (see
### http://sourceforge.net/projects/picard/)

### To get help on options specific to the subcommand, just 
### do ptrun.sh subcommand. For help on options common to 
### all subcommands, do ptrun subcommand --stdhelp. 

### Note that you can often use INPUT=/dev/stdin for reading from std input, and
### OUPUT=/devstdout QUIET=true for writing to stdout.

### VALIDATION_STRINGENCY=LENIENT is always used.

validation='VALIDATION_STRINGENCY=LENIENT'

## avoid "ERROR: Option 'VALIDATION_STRINGENCY' cannot be specified more than once."
for i in "$@"; do 
    if [ "$i" == "$validation" ]; then
        validation=""
    fi
done

pthelp() {
  pushd $pthome
  echo "
  Usage: ptrun.sh jarfile, where FILE is one of:
"
  ls -x *.jar | sed 's/\.jar//g'
  exit 1
}

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

if [ $# -lt 1 ] ; then
        pthelp
fi

jar="$1"
shift

java -d64 -Xmx512M -jar $pthome/$jar.jar "$@" $validation

