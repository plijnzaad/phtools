#!/usr/bin/env perl

# Print complete table of edit distances

use strict;
use Getopt::Std;
use Text::Levenshtein qw(distance);

use vars qw($opt_h);

### for list of strings (id \t barcode) on standard in, print 
### all pairwise edit distances;

my $Usage="$0 [-h] < file > output";

if ( !getopts("h") || $opt_h ) {
    die $Usage; 
}

my  $names=[];
my  $strings={};
my  $codes={};

### read codes
LINE:
while(<>) { 
  s/[ \t\n\r]*$//g;
  next LINE unless $_;
  my ($name, $code);
  if (/^#/) {
    ($name, $code)=($_, '');
  } else { 
    ($name, $code)=split(' ');            # e.g. 'G7 \t CCAACAAT'
    if (! $code) {                        # just one column of strings
      $code=$name;
      $name= "(line $.)";
    }
  }
  push @$names, $name;
  $codes->{$name}=$code;
  if ($strings->{$code}) { 
    die "duplicate string: '$code' (named '$name', line $.,)" 
        if $code ne '';
  }
  $strings->{$code}=$name if $code;
}                                       # while

my $n=int(@$names);
 I:
for(my $i=0; $i<$n; $i++) {
 J:
  my ($a, $x)=($names->[$i], $codes->{$names->[$i]});
  if ($x eq '') {
    print "$a\n";
    next;
  }
  for(my $j=0; $j<$n; $j++) { 
    my ($b, $y)=($names->[$j], $codes->{$names->[$j]});
    if ($y eq '') {
      print "\t\t";
      next;
    }
    print distance($x,$y), "\t";
  }
  print "\n";
}

