#!/bin/sh
# -*- mode: sh; -*-
## drop-in replacement of sambamba sort that makes use of the proper local TMPDIR on the compute node
sambamba sort --tmpdir=$TMPDIR "$@"
