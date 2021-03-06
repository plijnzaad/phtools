#!/bin/env perl
## Output all metrics in bam file (bwa output for now) from single-ended sequencing
## started with /home/gen/philip/git/phtools/ngs/summarize-sam.pl, 
## 9 may 2017

## Written by plijnzaad@gmail.com

use strict;
use Digest::MD5 qw(md5_base64);
use FileHandle;

use Number::Format;
my $fmt=new Number::Format(-thousands_sep => ',');
sub commafy {   $fmt->format_number($_[0]); }

my $samview="samtools view";            # alternatively, specify full path, or e.g. sambamba

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
my($FPE, $Fproper, $Funmapped, $Fmateunmapped, $Frev, $Fmaterev, $Fread1, $Fread2, 
 $Fsecondary, $FfailedQC, $Foptdup, $Fchimeric)=
    map { 1<<$_ } 0..11;

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

print join("\t", qw{mapq sense rank editdist score score2 seq pos read}) . "\n";

LINE:
while(<$fh>) { 
  s/[\n\r]*$//;
  next if /^@/;
  my($qname,$flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, 
     $seq, $qual, @tags)=split("\t", $_);
  die "paired-end only for now ... " if($flag & $FPE);
  next LINE if $flag & $Funmapped;
##  next LINE if $rname =~ /^ERCC-/;
  my $sense = ($flag & $Frev)? 'as' : 'sense';
  my $rank = ($flag & $Fsecondary)? 'secondary' : 'primary';
  my $h={};
  foreach my $el (@tags) {
    my ($name,$type,$val) = split(":",$el);
    ## $val = !!$val if $name eq 'SA';     # secondary alignment
    $h->{$name}=$val;
  }
  print join("\t", ($mapq, $sense, $rank, (map { my $v=$h->{$_}; defined($v)?$v:'NA'; } qw(NM AS XS)) , $rname, $pos,$qname)) . "\n";
}                                       # while
$fh->close();
