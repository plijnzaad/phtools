#!/bin/env perl
## Extract insert lengths from SAM file
## Reads from stdin, writes lengths to stdout
## 
## Usage e.g. 
##
##  samtools view foo.bam | sam-insertlengths.pl > foo.insertlen
##
## To select reads based on insert length, best use something like 
## 
##   sambamba view -h --filter 'template_length > 0 and template_length < 100'
## 
## To process this file in to a histrogram of insertsizes, see phtools/ngs/R/insertlen-distro.R
## reused parts from summarize-sam.pl commit 329eab4f0349ed2ab (2015-12-21 14:44:01)

## The script is not exactly fast (better use Bio::DB::Sam library?)

## Written by plijnzaad@gmail.com

use strict;

use Number::Format;
my $fmt=new Number::Format(-thousands_sep => ',');
sub commafy {   $fmt->format_number($_[0]); }

my $flags = ['PE',                      # 0
            'properly aligned',         # 1
            'unmapped',                 # 2
            'mate unmapped',            # 3
            'reverse strand',           # 4
            'mate on reverse strand',   # 5
            'read1', # i.e. first of two halfpairs # 6
            'read2', # i.e. second of two halfpairs # 7
            'part of a secondary alignment', # 8
            'not passing quality controls',  # 9
            'PCR or optical duplicate',      # 10
            'chimeric/supplementary alignment']; # 11

## make code slightly more readable:
my($FPE, $Fproper, $Funmapped, $Fmateunmapped, $Frev, $Fmaterev, $Fread1, $Fread2);
($FPE, $Fproper, $Funmapped, $Fmateunmapped, $Frev, $Fmaterev, $Fread1, $Fread2)=
    map { 1<<$_ } 0..7;

while(<>) { 
  s/[\n\r]*$//;
  next if /^@/;
  my($qname,$flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen,
     $seq, $qual, @optionals)=split("\t", $_);
  print "$tlen\n" unless $tlen <= 0;                # rest is 'the other mate', or unmapped or discordant or
}


