#!/bin/sh
# -*- mode: sh; -*-

## extract chromosome sizes from bam file

samtools view -H "$1" | sed -n '/^@SQ/{s/@SQ.*SN://;s/LN://;p}'
