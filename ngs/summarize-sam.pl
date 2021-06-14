#!/usr/bin/env perl

## simple tool to inspect a SAM file and make the flags human-readable.
## Reads from stdin, writes to stdout, stats to stderr. Aim is primarily
## to check paired-end stuff.
## 
## Usage e.g. 
##  samtools view foo.bam | summarize-sam.pl 2> foo.sumstat | gzip > foo.sum.gz

## It outputs a very short hash of the id (for fast/visual searching),
## chromosome, position, mateposition, insertlength, a summary of the
## sesquence (numbers of A,C,T,G,N), detail of paired-end status, a human
## readable version of all the flags, and lastly the original seqid.

## The thing has been written with bowtie2 output/terminology in mind.
## For one-off inspection of flags, check https://broadinstitute.github.io/picard/explain-flags.html

## Bowtie outputs insertlen ==0 for anything that is not 'properly aligned'.
## 'unmapped' means: this halfpair, or its mate, or both are unmapped.
## If one of the mates was mapped, chr, pos and matepos are that of the
## mapped mate (weird).

## The only combinations of flags for properly aligned read pairs are:
## PE; properly aligned; mate is on reverse strand; read1
## PE; properly aligned; mate is on reverse strand; read2
## PE; properly aligned; reverse strand; read1
## PE; properly aligned; reverse strand; read2

## The script is not exactly fast (better use Bio::DB::Sam library?). If
## you just need the percentage mapped etc., better use samtools flagstats in
## conjunction with pe-stats.pl (it outputs nearly the same numbers, but
## cannot distinguish between read1 unmapped vs. read2 unmapped).

## Written by plijnzaad@gmail.com

use strict;
use Digest::MD5 qw(md5_base64);
use FileHandle;

use Number::Format;
my $fmt=new Number::Format(-thousands_sep => ',');
sub commafy {   $fmt->format_number($_[0]); }

my $samview="samtools view";            # alternatively, specify full path, or e.g. sambamba

my $version="summarize-sam.pl v0.4";

my $flags = ['PE',                      # 0
            'proper pair',              # 1
            'unmapped',                 # 2
            'mate unmapped',            # 3
            'reverse strand',           # 4
            'mate on reverse strand',   # 5
            'read1', # i.e. first of two halfpairs # 6
            'read2', # i.e. second of two halfpairs # 7
            'part of a secondary alignment', # 8
            'not passing QC',                # 9
            'PCR or optical duplicate',      # 10
            'chimeric/supplementary alignment']; # 11

## make code slightly more readable:
my($FPE, $Fproper, $Funmapped, $Fmateunmapped, $Frev, $Fmaterev, $Fread1, $Fread2);
($FPE, $Fproper, $Funmapped, $Fmateunmapped, $Frev, $Fmaterev, $Fread1, $Fread2)=
    map { 1<<$_ } 0..7;

sub idhash {                            # make easier to spot what goes with what
  my ($id)=@_;
  substr(md5_base64($id), 0, 4); # 16 million possibilities
### for speed, consider using a checksum made of unpack('%...', $id)
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

sub seq_composition { 
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
}                                       # seq_composition

my $stats={};

sub print_stats { 
    my($stats)=@_;

    my $order=['total reads',
               'good',
               'bad',
               'one of mates unmapped',
               "  of which read1 unmapped",
               "  of which read2 unmapped",
               'both mates unmapped',
               'mate on other chromosome',
        ];
    my $ntot=$stats->{'total reads'};
    print STDERR "# $version. All numbers are per mate, not per readpair!\n";
    for my $key (@$order) { 
      print STDERR join("\t", ($key, commafy($stats->{$key}), sprintf("%.1f%%", 100*$stats->{$key}/$ntot)))."\n";
    }
}

sub pval { 
  my($q)=@_;
  sprintf("%4.1g", 10**(($q/-10)));
}

my $fh;
if (@ARGV) { 
  die "Only single file allowed " unless @ARGV==1;
  my $f=$ARGV[0];
  if ($f =~ /\.u?(b|cr)am$/) {
    my $cmd="$samview $f |";
    $fh = FileHandle->new($cmd) or die "$cmd: $!";
  } else {
    $fh = FileHandle->new("< $f") or die "$f: $!";
  }
} else {                                # stdin
    $fh = FileHandle->new_from_fd(0, "<") or die "stdin: $!";
}

sub header{
  my @fields=qw(idhash chr pos matepos inslen seq_composition hits score mapq pval PEsummary flags orig_id);
  print "#" . join("\t", @fields) . "\n";
}

header();

while(<$fh>) { 
  s/[\n\r]*$//;
  next if /^@/;

  my($qname,$flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen,
     $seq, $qual, @tags)=split("\t", $_);
  my $tags = join(" ", @tags);
  my($NH, $AS)=('?', '?');
  my $PEsummary="";
  $stats->{'total reads'}++;
  if($flag & $FPE) { 
    if ( ($flag & $Funmapped) && ($flag & $Fmateunmapped) ) {
      $PEsummary .= "; doubly unmapped"; 
      $stats->{'both mates unmapped'}++;
    }
    if (!!($flag & $Funmapped) != !!($flag & $Fmateunmapped)) { 
      $PEsummary .= "; singleton";
      $stats->{'one of mates unmapped'}++;
      $stats->{"  of which read1 unmapped"}++ if ($flag & $Fread1 && ($flag &$Funmapped));
      $stats->{"  of which read1 unmapped"}++ if ($flag & $Fread2 && ($flag &$Fmateunmapped));
      $stats->{"  of which read2 unmapped"}++ if ($flag & $Fread1 && ($flag &$Fmateunmapped));
      $stats->{"  of which read2 unmapped"}++ if ($flag & $Fread2 && ($flag &$Funmapped));
    }
    if ( !($flag & $Fproper) && !($flag & $Funmapped) && !($flag & $Fmateunmapped)) { 
      ## proper meaning: according to aligner
      $PEsummary .= "; bad";
      $stats->{'bad'}++;
    }
    if ( $rnext ne '=' && $rname ne $rnext) { 
      $PEsummary .= '; mate on other chromosome';
      $stats->{'mate on other chromosome'}++;
    };
    $PEsummary =~ s/^; //;
    if ( !$PEsummary) { 
      $PEsummary = "good";
      $stats->{good}++;
    }
  }
  $tags =~ /NH:i:(\d+)/; $NH=$1;
  $tags =~ /AS:i:(\d+)/; $AS=$1;
  $pos = '?' if $flag & $Funmapped;
  $pnext = '?' if $flag & $Fmateunmapped;
  $tlen= '?' if $flag & ($Funmapped | $Fmateunmapped);
  
  print join("\t", (idhash($qname), $rname, $pos, $pnext, $tlen, seq_composition($seq),
                    $NH, $AS,
                    $mapq, pval($mapq), $PEsummary, explain_flags($flag), $qname)) . "\n";
}                                       # while
$fh->close();
print_stats($stats);
