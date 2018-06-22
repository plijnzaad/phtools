#!/usr/bin/env perl
use strict;

## Exceedingly stupid script to clip content of bedGraph files that lie beyond the specified sizes.
## Clipping is a bit of a misnomer: any line with start or end coordinates beyond the chromosome size
## is fully obliterated.

my $usage = "Usage: gunzip< sample.bedgraph.gz | bedGraph-clip.pl chromsizes.txt | gzip > sample-clipped.bedgraph.gz\n";

die $usage unless @ARGV==1;

my $table=shift;

open(TAB, $table) || die "$table: $!";

my $chrlens={};

CHR:
while(<TAB>) { 
  s/#.*//;
  next CHR unless /(\S+)\s+(\d+)/;
  my ($chr, $len)= ($1,$2);
  $chrlens->{$chr}=$len;
}
close(TAB);

my $nlines=0;
my $nclipped=0;

LINE:
while(<>) {
  if (/^track/) {
    print;
    next LINE;
  }
  $nlines++;
  my ($chr,$start, $end,$value)=split("\t", $_);
  my $len = $chrlens->{$chr};
  die "Unknown chromosome '$chr'" unless defined($len);
  if ($start > $len || $end > $len) {
    $nclipped++;
    next LINE;
  }
  print;
}
warn sprintf("Clipped %d out %d lines (%.1f%%)\n", $nclipped, $nlines, 100*$nclipped/$nlines);
