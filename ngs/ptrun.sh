#!/bin/bash
# -*- mode: sh; -*-
## little wrapper script for  the picard-tools, to be used in conjunction with the alias in module.lua file
## (see also ptrun-lmod.lua.eg)
## which takes care of using the right jar, memory and tmpdir

if [ x$_pt_jar = x]; then
  echo "_pt_jar variable not defined" >&2
  exit 8
fi

if [  ! -f "$_pt_jar" ]; then
  echo "picard tools jar not found" >&2
  exit 9
fi

java $_pt_javaopts -jar $_pt_jar "$@"  $_pt_defaults

