#!/bin/sh
# -*- mode: sh; -*-
# show statistics of a past job
qacct -j  $@  | egrep '^(jobname|jobnumber|ru_wallclock|mem|maxvmem|io|jobargs|submission_time)'
