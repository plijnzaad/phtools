#!/bin/env perl
my $version = 'unknown';
## (leave as is, id is automagically adjusted by svn)

## The DanPos nucleosome caller requires paired mate ID's to end in 1 and
## 2. This does that (needs name-sorted sam file, and outputs the same

die "Better to use 'bamTobed -i file.bam > file.bed' and run it on that ...";

my $pg_printed=0;

my $cmdline= $0; # . join(" ", @ARGV);

my $second=0;

LINE:
while(<>) { 
  if (/^@/) { 
    print;
    next LINE;
  } else  {
    print "\@PG\tID:pairnumber.pl\tPN:pairnumber.pl\tVN:$version CL:\"$cmdline\"\n"
        unless $pg_printed++;
  }
  my @a=split("\t", $_);
  $a[0] = $second? "$a[0]:2" : "$a[0]:1";
  print join("\t", @a);
  $second = ! $second;
}
