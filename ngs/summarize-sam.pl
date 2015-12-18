#!/bin/env perl

## simple tool to inspect a SAM file and make the flags human-readable.
## Reads from stdin, writes to stdout

## Bowtie outputs insertlen ==0 for anything that is not 'properly aligned'.
## 'unmapped' means: this read, or its mate, or both are unmapped.
## If one of the mates was mapped, chr, pos and matepos are that of the
## mapped mate (weird).

## The only combinations of flags for properly aligned read pairs are:
## PE; properly aligned; mate is on reverse strand; read1
## PE; properly aligned; mate is on reverse strand; read2
## PE; properly aligned; reverse strand; read1
## PE; properly aligned; reverse strand; read2

use strict;
use Digest::MD5 qw(md5_base64);

my $flags = ['PE',                      # 0
            'properly aligned',         # 1
            'unmapped',                 # 2
            'mate unmapped',            # 3
            'reverse strand',           # 4
            'mate on reverse strand',   # 5
            'read1', # i.e. first of two halfpairs # 6
            'read2', # i.e. second of two halfpairs # 7
            'part of a secondary alignment', # 8
            'not passing quality controls',  # 9
            'PCR or optical duplicate',      # 10
            'chimeric/supplementary alignment']; # 11

sub idhash { 
  my ($id)=@_;
  substr(md5_base64($id), 0, 4); # 16 million possibilities
}

sub explain_flags { 
  my($flag)=@_;

  my $s;
  for(my $i=0; $i<int(@$flags); $i++) { 
    my $bit= 1<<$i;
    if ($flag & $bit) { 
      $s = "$s; " if $s;
      $s = $s . $flags->[$i];
    }
  }
  $s;
}                                       # explain_flags

sub seqsummary { 
  my($seq)=@_;

  my $occ={};
  for my $l (split('', $seq)) { 
    $occ->{$l}++;
  }

  my $sum="";
  for my $k (sort keys %$occ) { 
    $sum = $sum .  $occ->{$k}. $k;
  }
  $sum =~ s/(.*)/\L$1/;                 # for readability
  $sum;
}                                       # seqsummary


print "#idhash	chr	pos	matepos	insertlen	seqsummary	flags\n";
while(<>) { 
  s/[\n\r]*$//;
  next if /^@/;

  my($qname,$flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen,
     $seq, $qual, @optionals)=split("\t", $_);
  my $extra="";
  $extra = '** mate on other chromosome!' if $rnext ne '=' && $rname ne $rnext;
  print join("\t", (idhash($qname), $rname, $pos, $pnext, $tlen, seqsummary($seq),explain_flags($flag)), $extra) . "\n";
}
