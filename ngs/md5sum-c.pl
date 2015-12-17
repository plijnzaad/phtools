#!/bin/env perl

## Written by plijnzaad@gmail.com

use Getopt::Long;
use strict;
use File::Basename;

my $usage="md5sum-c.pl [-f | -b | -p path | -e 's/pattern/replacement/' ] [ somefile.md5sum]

Same as md5sum -c, but allows substitutions on the filename inside the
file with md5sum.

Options:


 -b	Use the file's basename (i.e., filename sans directories)
 -f     Use the md5sum file's dirname as the path for the checksummed file(s)
 -p PATH Prepend files's basename with different path
 -e CMD	Edit filename contained in \$_ using perl expression, typically involving s///;

To recursively check a directory tree containing md5sum files alongside
the checksummed files, do something like

  find . -name '*.md5sum' -exec md5sum-c.pl -f {} \; 

";

my $help=undef;
my $opt_eval=undef;
my $opt_basename=undef;
my $opt_path=undef;
my $opt_findmode=undef;

die $usage if GetOptions(
  "h|help"=> \$help,
  "e=s" => \$opt_eval,
  "p=s" => \$opt_path,
  "b" => \$opt_basename,
  "f" => \$opt_findmode,
    )==0 || $help || (defined($opt_basename) + defined($opt_path) + defined($opt_findmode) > 1);

$opt_eval="s|.*/||" if $opt_basename || $opt_path || $opt_findmode;

my($ignore, $dir);
if ($opt_findmode) { 
  die "need named md5sum file when using -f" unless @ARGV==1;
  ($ignore, $dir)=fileparse($ARGV[0]);
}

$dir=$opt_path if $opt_path;
$dir="." unless $dir;                # otherwise "/localfile"

while(<>) { 
  s/[\n\r]*$//;
  my($sum, $file)=split(/ {2}/, $_);    # yes, two spaces
  $_=$file;
  eval $opt_eval if $opt_eval;
  $file=$_;
  $file = "$dir/$file" if $opt_path || $opt_findmode ;
  system("echo '$sum  $file' | md5sum -c - \n"); # note: two spaces
}
