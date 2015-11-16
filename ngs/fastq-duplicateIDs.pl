#!/bin/env perl


require 'getopts.pl';

use strict;
use Number::Format;
use List::Util qw(sum max);             # pairmap not available in older versions ...

use vars qw($opt_h $opt_s $opt_r);

my $Usage = q{
Usage:  fastq-duplicateIDs.pl [-r]  < input.fastq  > output.txt 2> output.log

  Find duplicate IDs (NOT read sequences) in FASTQ files, and prints out
  the the multiplicities (in terms of different sequences per ID).  If
  requested, also does the reverse of that: prints out the non-unique
  sequences (far more!)  and prints out the multiplicities in terms of
  different IDs per sequence). Use fastq-uniquify.pl to get rid of duplicate IDs.

  To get rid of real duplicate IDs (i.e. those also having the same sequences), use
  fastq-uniquify.pl.

Options: 
  -s      Strip off the machine identifier to save memory
  -r      Also do the reverse (finding duplicate IDs per sequence);



};


&Getopts('hrs') || die $Usage;
die $Usage if $opt_h;

my $fmt=new Number::Format(-thousands_sep => ',');
 
sub commafy {
  $fmt->format_number($_[0]);
}


my($id_line, $seq, $plus, $qual, $id, $rest);

my $idmap={};
my $seqmap={};

my $nseqs=0;

while(1 &&  ! eof(STDIN) ) { 

  $id_line= <>;  
  my($id, $rest)=split(' ', $id_line);
  
  $id =~ s/^@[^:]*:// if $opt_s;

  die "incomplete block," if eof(STDIN);

  $seq= <>; $seq =~ s/[\n\r]*$//;

## use tie'd BerkeleyDB hash if it gets too big
  push( @{$idmap->{$id}}, $seq);

  push( @{$seqmap->{$seq}}, $id) if $opt_r;

  $plus = <>;
  die "incomplete block," if eof(STDIN);
  $qual = <>;

  $nseqs++;
  warn "Read ".commafy($nseqs)." seqs\n" if ($nseqs % 10000) == 0;
}

warn "Done reading, found ".commafy(int(keys %$idmap))." non-unique ids\n";

sub multiplicities {
  my($lst) = @_;

  my $occ={};

  for my $elt  (@$lst) { $occ->{$elt} ++; }
  map { "$_: $occ->{$_}" } keys %$occ
}

print "#ID\tsequence1_multiplicity1\tsequence2_multiplicity2\t etc.\n";

my ($n_uniq_id, $n_nonuniq_id, $sum_nonuniq_id);

for(keys %$idmap) { 
  if ( @{$idmap->{$_}} > 1 ) {
    print join("\t", $_, multiplicities($idmap ->{$_})). "\n";
    $n_nonuniq_id++;
    $sum_nonuniq_id += int(@{$idmap->{$_}});
  } else {
    $n_uniq_id++;
  }
}

warn sprintf("Total seqs: %s, unique ids: %s (%4.2f%%), non-unique ids: %s (%4.2f%%), seqs with non-unique id: %s (%4.2f%%)\n",
             commafy($nseqs),
             commafy($n_uniq_id), 100*($n_uniq_id/($n_uniq_id + $n_nonuniq_id)),
             commafy($n_nonuniq_id), 100*($n_nonuniq_id/($n_uniq_id + $n_nonuniq_id)), 
             commafy($sum_nonuniq_id),100*$sum_nonuniq_id/$nseqs);

if($opt_r) { 
  print "########################################################################\n";
  print "########################################################################\n";
  print "########################################################################\n";
  print "#seq\tID1_multiplicity1\tID2_multiplicity2\t etc.\n";
  
  for(keys %$seqmap) { 
    if ( @{$seqmap->{$_}} > 1 ) {
      print join("\t", $_, multiplicities($seqmap ->{$_})). "\n";
    } 
  }
}

