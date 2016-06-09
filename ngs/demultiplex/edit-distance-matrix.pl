#!/usr/bin/env perl
#
use strict;
use Getopt::Std;
use Text::Levenshtein qw(distance);

use vars qw($opt_h);

### for list of strings (id \t barcode) on standard in, print 
### all pairwise edit distances;

my $Usage="Usage: 

  $0 [-h] < file.txt > output.txt

Prints complete table of edit distances. File must be tab-delimited, two-columns (name TAB barcode)
";

if ( !getopts("h") || $opt_h ) {
    die $Usage; 
}

my  $names=[];
my  $strings={};
my  $codes={};

### read codes
LINE:
while(<>) { 
  s/#.*//;
  s/[ \t\n\r]*$//g;
  next LINE unless $_;
  my ($name, $code);
  ($name, $code)=split(' ');            # e.g. 'G7 \t CCAACAAT'
  if (! $code) {                        # just one column of strings, must be code
    $code=$name;
    $name= "(line $.)";
  }
  push @$names, $name;
  $codes->{$name}=$code;
  warn "duplicate string: '$code' (named '$name', line $.,)" 
      if $strings->{$code};
  $strings->{$code}=$name if $code;
}                                       # while

my $n=int(@$names);
print join("\t", ('', @$names)). "\n";
for(my $i=0; $i<$n; $i++) {
  my ($a, $x)=($names->[$i], $codes->{$names->[$i]});
  print "$a\t";
  for(my $j=0; $j<$n; $j++) { 
    my ($b, $y)=($names->[$j], $codes->{$names->[$j]});
    print distance($x,$y), "\t";
  }                                     # for $j
  print "\n";
}                                       # for $i

