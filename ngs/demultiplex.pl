#!/usr/bin/env perl

use strict;
use Getopt::Std;
use Text::Levenshtein qw(distance);
use FileHandle;

use vars qw($opt_h $opt_b $opt_m $opt_p $opt_o);

my $Usage="Usage:

   ... | $0 -b barcodes.txt [-m mismatches] [ -p outputprefix] [ -o outputdir ] 
";

if ( !getopts("b:p:o:m:h") || ! $opt_b ||  $opt_h ) {
    die $Usage; 
}

my  $mismatches_allowed = 1;
$mismatches_allowed = $opt_m if defined($opt_m);  # 0 also possible

die "only zero-mismatches allowed for now ..." if $mismatches_allowed >0;

my $special=0;

sub readbarcodes {
  my ($file)=@_;
  my $bc={};
  my $libs={};

  open(FILE, "$file") or die "Barcode '$file': $!";
  while(<FILE>) {
    s/[\n\r]*$//g;
    s/#.*//;
    next LINE unless $_;
    my ($lib, $code)=split(' ');            # e.g. 'G7 \t CCAACAAT'
    die "Library '$lib' not unique" if $libs->{$lib}++;
    die "Barcode '$code' not unique" if $bc->{$code};
    $bc->{$code}=$lib;
  }
  close(FILE);
  $bc;
}                                       # readbarcodes

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
    $name="$opt_p$name" if $opt_p;
    $name="$opt_o/$name" if $opt_o;
    my $fh = FileHandle->new("| gzip > $name") or die "library $lib, file $name: $!";
    warn "Creating file $name ...\n";
    $fhs->{$lib}=$fh;
  }
  $fhs;
}

sub close_outfiles {
  my($fhs)=@_;
  for my $lib (keys %$fhs) {
    $fhs->{$lib}->close() or die "could not close demultiplexed file for library $lib; investigate";
  }
}

sub rescue {
  die "to be written";
}

sub ambiguous {
  die "to be written";
}

my $codes = readbarcodes($opt_b);       # eg. $code->{'AGCGTT') => 'M3'
my @files=(values %$codes, 'AMBIGUOUS', 'UNKNOWN');
my $filehandles=open_outfiles(@files);      # opens M3.fastq.gz, ambiguous.fastq.gz etc.

my $nexact=0;
my $nrescued=0;                         # having at most $mismatch mismatches
my $nunknown=0;
my $nambiguous=0;

RECORD:
while(1) { 
  my $record=<>;
  last RECORD if (eof(STDIN) || !$record);
  ### e.g.:  ^@NS500413:172:HVFHWBGXX:1:11101:4639:1062 1:N:0:CCGTCCAT$
  my ($code)=(split(':', $record))[-1];
  $code =~ s/[\n\r]*$//;
  $record .= <>; # sequence line
  $record .= <>; # '+'
  $record .= <>; # quality line
  
  my $lib;
 CASE:
  while(1) {
    $lib=$codes->{$code};
    if ($lib) {
      $nexact++;
      last CASE;
    }
    if ($mismatches_allowed == 0) {
      $nunknown++;
      $lib='UNKNOWN';
      last CASE;
    }
    $lib=rescue($code, $codes, $mismatches_allowed);
    if(!$lib) {
      $nunknown++;
      $lib='UNKNOWN';
      last CASE;
    }
    if($special) {
      # check if mismatch is in 7th bp; if so, call it ambiguous
      if (ambiguous() ) {
        $lib='AMBIGUOUS';
        $nambiguous++;
        last CASE;
      } else {
        $nrescued++;
        last CASE;
      }
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

warn sprintf("exact: %s\nrescued:%s\nambiguous:%s\nunknown: %s\n",
             map { commafy $_ } ($nexact, $nrescued, $nambiguous, $nunknown ));
