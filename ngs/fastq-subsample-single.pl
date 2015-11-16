#!/bin/env perl

### NOTE: also check samtools view -s 322.15 -b file.bam > random_15%_of_file.bam

require 'getopts.pl';

use strict;
use Number::Format;
## use Math::Random;                       # broken! rand() is much better

use vars qw($opt_h $opt_p $opt_s);

my $Usage = q{
Usage:  

    fastq-subsample-single.pl  -p percentage  < all.fastq  > subset.fastq

Does approximate subsampling in one pass for single reads. Ordering will
be the same as that of the input file. For paired end reads use
fastq-subsample-paired.pl


Options: 
  -s  N   Seed the random number generator with N,N (for reproduceability purposes)
};


&Getopts('hp:s:') || die $Usage;
die $Usage if ($opt_h or !$opt_p);


if ($opt_s) {
##   random_set_seed(($seed,$seed));
  srand($seed);
}
## warn "random seed was ", join(" ", random_get_seed()), "\n";

my $frac= 1 - ($opt_p/100);

my $fmt=new Number::Format(-thousands_sep => ',');
 
sub commafy {
  $fmt->format_number($_[0]);
}

my ($id_line, $seq, $plus, $qual, $nseqs, $nselected);

while(1 &&  ! eof(STDIN) ) { 

  $id_line= <>;  
  die "incomplete block," if eof(STDIN);
  $seq= <>;
  $plus = <>;
  die "incomplete block," if eof(STDIN);
  $qual = <>;

  $nseqs++;

##  if( random_uniform() > $frac) { 
  if( rand(1) > $frac) { 
    $nselected++;
    print "$id_line$seq$plus$qual";
  }
}
warn(sprintf("Selected %s out of %s reads (%4.2f%%)\n", 
             commafy($nselected), commafy($nseqs), 100*$nselected/$nseqs));
