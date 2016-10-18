#!/usr/bin/env perl
## written by <plijnzaad@gmail.com>
##
## To test, do e.g. 
##   ./demultiplex.pl -m 1 < testdata/one-mismatch.fastq  -b testdata/testbarcodes.txt -p DEMUL
## 
## Given a barcode file, demultiplexes a FASTQ file (on stdin) while potentially allowing for mismatches.
## For format of barcode file, see testdata/testbarcodes.txt.
##

use strict;
use Getopt::Std;
use FileHandle;

use mismatch;

use vars qw($opt_h $opt_b $opt_m $opt_p $opt_o);

my $Usage="Usage:

   ... | $0 -b barcodes.txt [-m mismatches] [ -p outputprefix ] [ -o outputdir ] 

NOTE: the script does *not* check if mismatched barcodes are unambiguous!
Use edit-distance.pl and/or edit-distance-matrix.pl for that.

";

if ( !getopts("b:p:o:m:h") || ! $opt_b ||  $opt_h ) {
    die $Usage; 
}

my  $allowed_mismatches = 1;
$allowed_mismatches = $opt_m if defined($opt_m);  # 0 also possible, meaning no mismatched allowed

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

my $barcodes = mismatch::readbarcodes($opt_b); ## eg. $h->{'AGCGTT') => 'M3'
my $mismatch_REs = mismatch::convert2mismatchREs(barcodes=>$barcodes, allowed_mismatches =>$allowed_mismatches);# eg. $h->{'AGCGTT') =>  REGEXP(0x25a7788)

my @files=map { "\U$_" } (values %$barcodes, 'UNKNOWN');
my $filehandles=open_outfiles(@files);      # opens M3.fastq.gz, ambiguous.fastq.gz etc.

my $nexact=0;
my $nmismatched=0;                         # having at most $mismatch mismatches
my $nunknown=0;

## lastly, process the actual input:
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
    if (! $allowed_mismatches) {
      $nunknown++;
      $lib='UNKNOWN';
      last CASE;
    }
    $lib=mismatch::rescue($foundcode, $barcodes, $mismatch_REs); # takes longish (regexp matching)
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
