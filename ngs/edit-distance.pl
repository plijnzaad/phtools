#!/usr/bin/env perl

# Code to find near-duplicates among barcodes

use strict;
use Getopt::Std;
use Text::Levenshtein qw(distance);

use vars qw($opt_h $opt_l);

### for list of strings on standard in, print the closest distance (and
### corresponding string)

my $Usage="$0 [ -l number ] < file > output";

if ( !getopts("hl:") || $opt_h ) {
    die $Usage; 
}

my $limit= $opt_l || 1;

my  $strings={};

### read codes
LINE:
while(<>) { 
  s/#.*//;
  s/[ \t\n\r]*$//g;
  next LINE unless $_;
  my ($name, $code)=split(' ');            # e.g. 'G7 \t CCAACAAT'
  if (! $code) {                           # just one column of strings
    $code=$name;
    $name= "(line $.)";
  }
  die "duplicate string: '$code' (named '$name', line $.,)" if $strings->{$code};
  $strings->{$code}=$name;
}                                       # while

my @strings = sort keys %$strings;

print "# distance with edit distance <= $limit:\n";
SEQ:
for(my $i=0; $i<@strings; $i++) { 
  my @s=@strings;
  my $s=$s[$i];
  splice(@s, $i, 1);
  my @d = distance($s , @s);
  my @hits = grep( $d[$_] <= $limit, 0..$#d);
  next SEQ unless @hits;
  print $strings->{$s}. "  $s: "; 
  for(my $h=0; $h<@hits; $h++) {
    my $hs=$s[ $hits[$h] ];
    print $strings->{$hs} . " ($hs, d=$d[ $hits[$h] ]); ";
  }
  print "\n";
}

#  print $strings->{$s}. "\t$s:" . join(", ", @d). "\n";
#  my @o= sort { $d[$a] <=>  $d[$b] } a0..int(@d);      # mimic R's order() function
