#!/usr/bin/env perl
use strict;
use Getopt::Std;
## use Text::Levenshtein qw(distance);     

use vars qw($opt_h $opt_l);

my $Usage="Find near-duplicates among barcodes. 

For list of strings on stdin (format: id \\t barcode), print all the
sets of strings with distance less than -l to each other.

Usage: 
  $0 [ -l number ] < file > output
";

if ( !getopts("hl:") || $opt_h ) {
    die $Usage; 
}

my $limit= $opt_l || 1;


sub hammingdist {                       
## honestly stolen from http://www.perlmonks.org/?node_id=500244
  length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] );
}

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
  ## my @d = distance($s , @s); ### Levenshtein edit distance, but we don't allow indels
  my @d = map { hammingdist($_, $s); } @s;
  my @hits = grep( $d[$_] <= $limit, 0..$#d);
  next SEQ unless @hits;
  print $strings->{$s}. "\t$s\n"; 
  for(my $h=0; $h<@hits; $h++) {
    my $hs=$s[ $hits[$h] ];
    print $strings->{$hs} . "\t$hs (d=$d[ $hits[$h] ])\n";
  }
  print "\n";
}
