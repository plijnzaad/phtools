#!/bin/bash
# -*- mode: sh; -*-
## little wrapper script for  the picard-tools, to be used in conjunction with the alias in module.lua file
## (see also ptrun-lmod.lua.eg)
## which takes care of using the right jar, memory and tmpdir

if [ x$_ptjar = x -o ! -f $_ptjar ] ; then # defined in the 
  echo "picard tools jar file not defined or not found" >&2
  exit 8
fi

java $_pt_javaopts -jar $_ptjar "$@"  $_ptdefaults

