#!/usr/bin/env perl

# Usage: bam2wig.pl foo.bam > foo.wig

use strict;

my $bam=shift;
## my $cmd="samtools depth $bam";
my $window=1000;
my $cmd="sambamba depth window -w $window $bam ";

$bam =~ s/\.bam/wig/;

print "track type=wiggle_0 name=$bam description=$bam\n";

my $prevchrom = 'none';

open(BAM, "$cmd | ") or die "$cmd: $!";

while(<BAM>) {
  my ($chrom, $start, $end, $count, $mean, $sample) = split("\t");
  next if /^#/;
  next unless $count > 0; 
  print "variableStep chrom=$chrom span=$window\n"  if $chrom ne $prevchrom;
  $prevchrom=$chrom;
  ### next unless $. % $window ==0;
  print "$start\t$count\n"; #  unless $depth<3;
}
close(BAM); #  or die "$cmd: $!";
