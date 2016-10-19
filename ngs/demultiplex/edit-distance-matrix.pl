#!/usr/bin/env perl
#
use strict;
use Getopt::Std;
## use Text::Levenshtein qw(distance);

use vars qw($opt_h);

my $Usage="Usage: 

  $0 [-h] < file.txt > output.txt

Prints complete table of edit distances. File must be tab-delimited, two-columns (name TAB barcode).
Output (to stdout) is likewise tab-delimited.
";

if ( !getopts("h") || $opt_h ) {
    die $Usage; 
}

my  $names=[];
my  $strings={};
my  $codes={};

sub hammingdist {                       
## honestly stolen from http://www.perlmonks.org/?node_id=500244
  length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] );
}

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
    print distance($x,$y), "\t"; ### Levenshtein edit distance, but we don't allow indels
###    print hammingdistance($x,$y), "\t";
  }                                     # for $j
  print "\n";
}                                       # for $i

