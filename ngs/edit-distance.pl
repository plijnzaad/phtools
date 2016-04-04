#!/usr/bin/env perl

use strict;
use Getopt::Std;

use vars qw($opt_h $opt_l);

### for list of strings on standard in, print the closest distanance (and correpsonding string)
my $Usage="$0 [ -l number ] < file > output";

if ( !getopts("hl:") || $opt_h ) {
    die $Usage; 
}

my $limit= $opt_l || 1;

use Text::Levenshtein qw(distance);

my  $strings={};

LINE:
while(<>) { 
  s/[\n\r]*$//g;
  s/#.*//;
  next LINE unless $_;
  my ($name, $s)=split(' ');            # e.g. 'G7 \t CCAACAAT'
  $s=$name unless $name;                # just one column of strings
  $name= "(line $.)" unless $name;
  die "duplicate string: '$s' (named '$name', line $.,)" if $strings->{$_};
  $strings->{$s}=$name;
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
