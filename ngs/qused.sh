#!/bin/sh
# -*- mode: sh; -*-
# show statistics of a past job. Note that this may take fairly long

if [ -z "$1" ]; then
    echo "Usage: $0 jobid" >&2
    exit 2
fi

qacct -j  $@  | egrep '^(jobname|jobnumber|ru_wallclock|mem|maxvmem|io|jobargs|submission_time)'
