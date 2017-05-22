#!/bin/sh
# -*- mode: sh; -*-
# show statistics of a past job. Note that this may take fairly long

if [ -z "$1" ]; then
    echo "Usage: $0 jobid" >&2
    exit 2
fi

qacct -j  $@  | egrep '^(jobname|jobnumber|ru_wallclock|ru_utime|ru_stime|ru_maxrss|mem|maxvmem|io|jobargs|submission_time)'
## See http://gridscheduler.sourceforge.net/htmlman/htmlman5/accounting.html for interpretation of 
## these fields
