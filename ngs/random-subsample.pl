#!/bin/env perl

### NOTE: also check samtools view -s 322.15 -b file.bam > random_15%_of_file.bam, use seed 322.15
###
use strict;

use Getopt::Long;
use FileHandle;
use Number::Format;

use vars qw($help $seed $perc $nlines $nselected);

my $Usage = q{
Usage: 

   ... | random-subsample.pl  --perc percentage | ...

Does _approximate_ subsampling in one pass. Ordering will be the same as that of the 
input files. If the latter is undesirable, consider using the Unix command shuf (1) 
and/or reservoir sampling, see http://alexpreynolds.github.io/sample/

Options: 
  --perc  N   Seed the random number generator with N,N (for reproduceability purposes)
  --seed  N   Seed the random number generator with N,N (for reproduceability purposes)
};


die $Usage unless &GetOptions('help' => \$help,
                              'seed=i'=> \$seed,
                              'perc=f'=> \$perc);

die $Usage if ($help or !$perc);

if ($seed) {
#   random_set_seed(($seed,$seed));
  srand($seed);
}
## warn "random seed was ", join(" ", random_get_seed()), "\n";

my $frac= 1 - ($perc/100);

my $fmt=new Number::Format(-thousands_sep => ',');
 
sub commafy {   $fmt->format_number($_[0]); }

while(<>) {
  $nlines++;
  if( rand(1) > $frac) { 
    $nselected++;
    print;
  }
}                                       # while

warn(sprintf("Selected %s out of %s lines (%4.2f%%)\n", 
             commafy($nselected), commafy($nlines), 100*$nselected/$nlines));
