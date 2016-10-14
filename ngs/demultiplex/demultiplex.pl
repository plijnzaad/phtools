#!/usr/bin/env perl
## written <plijnzaad@gmail.com>
## to test, do e.g. 
##   ./demultiplex.pl -m 1 < testdata/one-mismatch.fastq  -b testdata/testbarcodes.txt -p DEMUL

use strict;
use Getopt::Std;
use FileHandle;
use Math::Combinatorics;
use Regexp::Optimizer;

use vars qw($opt_h $opt_b $opt_m $opt_p $opt_o);

my $Usage="Usage:

   ... | $0 -b barcodes.txt [-m mismatches] [ -p outputprefix ] [ -o outputdir ] 

NOTE: the script does *not* check if mismatched barcodes are unambiguous!
Use edit-distance.pl and/or edit-distance-matrix.pl for that.

";

if ( !getopts("b:p:o:m:h") || ! $opt_b ||  $opt_h ) {
    die $Usage; 
}

my  $mismatches_allowed = 1;
$mismatches_allowed = $opt_m if defined($opt_m);  # 0 also possible

my $o=Regexp::Optimizer->new;


sub _getmismatch_REs {
  ## set up regular expressions to allow mismatches
  my($code, $max_mm)=@_;
  
  return () if ! $max_mm;

  my @mmcodes=();
  my(@code)=split(//, $code);

  ## set up array of arrays with '.' where to do the replacements:
  for(my $i=0; $i<$max_mm; $i++) { 
    my @combs = combine(($i+1), 0..$#code);
    foreach my $comb ( @combs ) { 
      my @mm=@code;
      @mm[ @$comb ] = split(//, '.' x int(@$comb) ); # yay, splicing
      push(@mmcodes, join("", @mm));
    }
  }
  @mmcodes;
}                                       # getmismatch_REs

sub readbarcodes {
  ### die "@@@@ : split reading and converting into two subs";

  ## returns list  ($barcodes, $mm_REs);
  ## $barcodes maps bardcodes to IDs; $mm_REs maps barcodes to mismatch regexps
  ## eg. $barcodes->{'AGCGTT') => 'M3'                  }
  ## eg. $mismatch_REs->{'AGCGTT') =>  REGEXP(0x25a7788)

  my $args = ref $_[0] eq 'HASH' ? shift : {@_};
  my ($file, $allowed_mms)=($args->{file}, $args->{allowed_mms});
  my $libs={};
  my $barcodes = {};
  my $mm_REs = {};

  open(FILE, "$file") or die "Barcode '$file': $!";
LINE:
  while(<FILE>) {
    s/[\n\r]*$//g;
    s/#.*//;
    next LINE unless $_;
    my ($lib, $code)=split(' ');            # e.g. 'G7 \t CCAACAAT'
    die "Library '$lib' not unique" if $libs->{$lib}++;
    die "Barcode '$code' not unique" if $barcodes->{$code};
    $barcodes->{$code}=$lib;

    my @res= _getmismatch_REs($code, $mismatches_allowed); # empty if $mismatches_allowed==0

    my $r='^'.join("|", @res).'$';
    $r=$o->optimize(qr/$r/);
    $mm_REs->{$code}= $r;         # just one big regexp!
  }                                     # while LINE
  close(FILE);
  ( $barcodes, $mm_REs);
}                                       # readbarcodes

sub rescue { 
  my($foundcode, $barcodes, $mm_REs)=@_;

  foreach my $code (keys %$barcodes) {
    my $re=$mm_REs->{$code};
    return  $barcodes->{$code} if $foundcode =~ $re;
  }
  return undef;
}                                       # rescue

sub open_infile {
  die "not used nor tested";
  my($file)=@_;
  my $fh=FileHandle->new();
  if ($file =~ /\.gz/) { 
    $fh->open("zcat $file | ", "r")  or die "'$file': $!";
  } else { 
    $fh->open("< $file")  or die "'$file': $!";
  }
  $fh;
}

sub open_outfiles { 
  my(@libs)=@_;
  my $fhs={};

  for my $lib (@libs) { 
    my $name=sprintf("%s.fastq.gz", $lib);
##    my $name=sprintf("%s.fastq", $lib);
    $name="$opt_p$name" if $opt_p;
    $name="$opt_o/$name" if $opt_o;
    my $fh = FileHandle->new("| gzip > $name") or die "library $lib, file $name: $!";
##    my $fh = FileHandle->new(" > $name") or die "library $lib, file $name: $! (did you create the output directory?)";
    warn "Creating/overwriting file $name ...\n";
    $fhs->{$lib}=$fh;
  }
  $fhs;
}                                       # open_outfiles

sub close_outfiles {
  my($fhs)=@_;
  for my $lib (keys %$fhs) {
    $fhs->{$lib}->close() or die "could not close demultiplexed file for library $lib; investigate";
  }
}

my ($barcodes, $mismatch_REs) = readbarcodes(file=>$opt_b, allowed_=>$mismatches_allowed);
## eg. $barcodes->{'AGCGTT') => 'M3'
## eg. $mismatch_REs->{'AGCGTT') =>  REGEXP(0x25a7788)

my @files=(values %$barcodes, 'UNKNOWN');
my $filehandles=open_outfiles(@files);      # opens M3.fastq.gz, ambiguous.fastq.gz etc.

my $nexact=0;
my $nmismatched=0;                         # having at most $mismatch mismatches
my $nunknown=0;

RECORD:
while(1) { 
  my $record=<>;
  last RECORD if (eof(STDIN) || !$record);
  ### e.g.:  ^@NS500413:172:HVFHWBGXX:1:11101:4639:1062 1:N:0:CCGTCCAT$
  my ($foundcode)=(split(':', $record))[-1];
  $foundcode =~ s/[\n\r]*$//;
  $record .= <>; # sequence line
  $record .= <>; # '+'
  $record .= <>; # quality line
  
  my $lib;
 CASE:
  while(1) {
    $lib=$barcodes->{$foundcode};            # majority of cases
    if ($lib) {
      $nexact++;
      last CASE;
    }
    if (! $mismatches_allowed) {
      $nunknown++;
      $lib='UNKNOWN';
      last CASE;
    }
    $lib=rescue($foundcode, $barcodes, $mismatch_REs); # takes longish (regexp matching)
    if($lib) {
      $nmismatched++;
      last CASE;
    } else { 
      $nunknown++;
      $lib='UNKNOWN';
      last CASE;
    }
    die "should not reach this point";
  }                                     # CASE
  $filehandles->{$lib}->print($record);
}                                       # RECORD
close_outfiles($filehandles);

sub commafy {
  # insert comma's to separate powers of 1000
  my($i)=@_;
  my $r = join('',reverse(split('',$i)));
  $r =~ s/(\d{3})/$1,/g;
  $r =~ s/,$//;
  join('',reverse(split('',$r)));
}

warn sprintf("exact: %s\nmismatched: %s\nunknown: %s\n", 
             map { commafy $_ } ($nexact, $nmismatched, $nunknown ));
