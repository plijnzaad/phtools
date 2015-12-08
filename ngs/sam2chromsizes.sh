#!/bin/sh
# -*- mode: sh; -*-

## extract chromosome sizes from sam/bam file (use '-' to use stdin)

samtools view -S -H "$1" | sed -n '/^@SQ/{s/@SQ.*SN://;s/LN://;p}'
