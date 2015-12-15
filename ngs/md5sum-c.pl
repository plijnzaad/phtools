#!/bin/env perl

## Written by plijnzaad@gmail.com

use Getopt::Long;
use strict;


my $usage="md5sum-c.pl [-b] [ -e 's/pattern/replacement/' ] [ somefile.md5sum]

Same as md5sum -c, but allows substitutions on the filename inside the
file with md5sum.

Options:

 -b	Use the file's basename (i.e., filename sans directories)
 -e CMD	Edit filename contained in \$_ using perl expression, typically involving s///;

";

my $help=undef;
my $opt_eval=undef;
my $opt_basename=undef;

die $usage if GetOptions(
  "h|help"=> \$help,
  "e=s" => \$opt_eval,
  "b" => \$opt_basename,
    )==0 || $help;

$opt_eval="s|.*/||" if $opt_basename;

while(<>) { 
  s/[\n\r]*$//;
  my($sum, $file)=split(/ {2}/, $_);    # yes, two spaces
  $_=$file;
  eval $opt_eval if $opt_eval;
  system("echo '$sum  $_' | md5sum -c - \n"); # note: two spaces
}
