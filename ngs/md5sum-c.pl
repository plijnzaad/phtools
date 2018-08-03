#!/bin/env perl

## Written by plijnzaad@gmail.com

use Getopt::Long;
use strict;
use File::Basename;

my $usage="md5sum-c.pl [ -q | -f | -b | -p path | -e 's/pattern/replacement/' ] [ somefile.md5sum ]

Same as md5sum -c, but allows more flexibility with the files mentioned in the
2nd column of *.md5sum file.

Options:

 -b Use the file's basename (i.e., filename without directory-part). 
    I.o.w.: assume the file to be checked is in the current directory

 -q quiet mode: don't complain about files file that cannot be
    found ( typically when using -b )

 -f find(1) mode: use the md5sum file's dirname as the path for the
    checksummed file(s)

 -p PATH Prepend files's basename with PATH

 -e CMD Transform full path name of file (contained in \$_) using a perl
    expression (typically involving s///)

To recursively check a directory tree containing md5sum files alongside
the checksummed files, do something like

  find . -name '*.md5sum' -exec md5sum-c.pl -f {} \; 

For *.gzipped files, an extra check is done to make sure that the .gz file
does not include a time-stamp.  (In that case the gzip was run in standard
mode, which renders the .gz useless for checksumming. Detecting this
before checksumming is of course faster)";

#'"; # quote to fool emacs

my $help=undef;
my $opt_eval=undef;
my $opt_basename=undef;
my $opt_path=undef;
my $opt_findmode=undef;
my $opt_quiet=undef;

die $usage if GetOptions(
  "h|help"=> \$help,
  "q|quiet" => \$opt_quiet,
  "e=s" => \$opt_eval,
  "p=s" => \$opt_path,
  "b" => \$opt_basename,
  "f" => \$opt_findmode) ==0 
    || $help || (defined($opt_basename) + defined($opt_path) + defined($opt_findmode) > 1);

$opt_eval="s|.*/||" if $opt_basename || $opt_path || $opt_findmode;

my($ignore, $dir);
if ($opt_findmode) { 
  die "need named md5sum file when using -f" unless @ARGV==1;
  ($ignore, $dir)=fileparse($ARGV[0]);
}

$dir=$opt_path if $opt_path;
$dir="." unless $dir;                # otherwise "/localfile"

LINE:
while(<>) { 
  s/[\n\r]*$//;
  my($sum, $file)=split(/ {2}/, $_);    # yes, two spaces
  $_=$file;
  eval $opt_eval if $opt_eval;
  $file=$_;
  $file = "$dir/$file" if $opt_path || $opt_findmode ;

  if (! -r $file )  {
    warn "File $file not found or readable" unless $opt_quiet;
    next LINE;
  }

  if ($file =~ /\.gz$/ ) { 
    my $cmd="dd bs=4 skip=1 count=1 < $file 2>/dev/null | od -t d4 | sed -n '1s/0*[ \t]*//;p;q;'";
    ## partly nicked from Gilles, http://unix.stackexchange.com/a/79546
    open(FILE, "$cmd |") or die "$cmd: $!";
    my $line=<FILE>;chomp($line);
    close(FILE) or die "$cmd: $!";     # should always work

    if ($line ne '0' ) { 
      warn "File $file was not gzipped properly. It includes a time-stamp ($line) which renders it useless
for checksumming (consider checking such gzipped files using gzip -t)";
      next LINE;
    }
  }
  system("echo '$sum  $file' | nice md5sum -c - \n"); # note: two spaces
}
