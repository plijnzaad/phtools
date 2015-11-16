#!/bin/env perl


require 'getopts.pl';

use strict;
use Number::Format;

use vars qw($opt_h);

my $Usage = q{
Usage:  fastq-uniquify.pl < withduplicates.fastq  > unique.fastq 2> output.log

  Get rid of spurious duplicates (i.e. same ID and same sequence+score) in fastq
  file.  Genuine duplicates (different ID but same sequence) are kept. An
  error is signalled if the same ID occurs with different sequences. The output is
  sorted by ID. 

  See also fastq-duplicateIDs.pl

};


&Getopts('hs') || die $Usage;
die $Usage if $opt_h;

my $fmt=new Number::Format(-thousands_sep => ',');
 
sub commafy {
  $fmt->format_number($_[0]);
}


my($id_line, $seq, $plus, $qual, $id, $rest);

my $idmap={};
my ($nseqs, $ndups);

my $machine;

while(1 &&  ! eof(STDIN) ) { 

  $id_line= <>;  
  my($id, $rest)=split(' ', $id_line, 2);
  
  $id =~ s/^@([^:]*)://;                # save space in hash

  if(! $machine ) { 
    $machine = $1;
  } else {
    die "Machine field changes from $machine to $1, line $.," unless $machine eq $1;
  }

  die "incomplete block," if eof(STDIN);

  $seq= <>;
  $plus = <>;
  die "incomplete block," if eof(STDIN);
  $qual = <>;

  my $value="$rest$seq$plus$qual";

## use tie'd BerkeleyDB hash if it gets too big
  if ( defined ($idmap->{$id}) ) { 
    die "id \@$machine:$id\n has more than one sequence/score: \n$idmap->{$id}\n --vs-- \n$value\n, line $.,"
        if $idmap->{$id} ne $value;
    $ndups++;
  }else {
    $idmap->{$id}=$value;
  }

  $nseqs++;
  warn "Read ".commafy($nseqs)." seqs\n" if ($nseqs % 10000) == 0;
}

warn "Done reading ".commafy($nseqs)." reads, found ".commafy(int(keys %$idmap))." unique ids and ".commafy($ndups)." superfluous sequences\n";

foreach my $id (sort keys %$idmap) { 
  print "\@$machine:$id $idmap->{$id}";
}

