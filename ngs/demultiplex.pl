#!/usr/bin/env perl

use strict;
use Getopt::Std;
use FileHandle;
use Math::Combinatorics;

use vars qw($opt_h $opt_b $opt_m $opt_p $opt_o);

my $Usage="Usage:

   ... | $0 -b barcodes.txt [-m mismatches] [ -p outputprefix] [ -o outputdir ] 
";

if ( !getopts("b:p:o:m:h") || ! $opt_b ||  $opt_h ) {
    die $Usage; 
}

my  $mismatches_allowed = 1;
$mismatches_allowed = $opt_m if defined($opt_m);  # 0 also possible

warn "allowing mismatches, untested yet ..." if $mismatches_allowed >0;

our $barcodes = {};        # eg. $barcodes->{'AGCGTT') => 'M3'
our $mismatch_REs = {};    # eg. $mismatch_REs->{'AG.GTT') => 'M3'
our $compiled_REs = {};    # eg. $compiled_REs->{'AG.GTT') => REGEXP(0x25a7788)
our @all_REs;

sub getmismatch_REs {
  my($code, $max_mm)=@_;

  my @mmcodes=();
  my(@code)=split(//, $code);

  ## set up array of arrays with '.' where to do the replacements:
  for(my $i=0; $i<$max_mm; $i++) { 
    my @combs = combine(($i+1), 0..$#code);
    foreach my $comb ( @combs ) { 
      my @mm=@code;
      @mm[ @$comb ] = split(//, '.' x int(@$comb) ); # yay, splicing
      push(@mmcodes, [@mm]);
    }
  }
  map { join("", @$_);} @mmcodes;
}                                       # getmismatch_REs

sub readbarcodes {
  my ($file)=@_;
  my $libs={};

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
    next LINE if $mismatches_allowed==0;

    my @res=getmismatch_REs($code, $mismatches_allowed);

    for my $re (@res) { 
      die "Mismatch RE '$re' for code $code with $mismatches_allowed mismatches ( library $lib) is not unique" if $mismatch_REs->{$re};
      $mismatch_REs->{$re}=$lib;
      my $r="^$re\$";
      $compiled_REs->{$re}=  qr/$r/;
    }
  }
  close(FILE);
  undef;
}                                       # readbarcodes

sub rescue { 
  my($code)=@_;
 
  foreach my $re (@all_REs) { 
    return $mismatch_REs->{$re} if  $code =~ $compiled_REs->{$re};
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
    $name="$opt_p$name" if $opt_p;
    $name="$opt_o/$name" if $opt_o;
    my $fh = FileHandle->new("| gzip > $name") or die "library $lib, file $name: $!";
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

readbarcodes($opt_b);

@all_REs = keys %$mismatch_REs;

my @files=(values %$barcodes, 'AMBIGUOUS', 'UNKNOWN');
my $filehandles=open_outfiles(@files);      # opens M3.fastq.gz, ambiguous.fastq.gz etc.

my $nexact=0;
my $nmismatched=0;                         # having at most $mismatch mismatches
my $nunknown=0;

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
    $lib=$barcodes->{$code};
    if ($lib) {
      $nexact++;
      last CASE;
    }
    if ($mismatches_allowed == 0) {
      $nunknown++;
      $lib='UNKNOWN';
      last CASE;
    }
    $lib=rescue($code);
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

warn sprintf("exact: %s\nrescued:%s\nambiguous:%s\nunknown: %s\n", map { commafy $_ } ($nexact, $nmismatched, $nunknown ));
