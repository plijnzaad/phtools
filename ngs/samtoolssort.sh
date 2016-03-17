#!/bin/sh
# -*- mode: sh; -*-

samtools sort -T$TMPDIR/sort "$@"
