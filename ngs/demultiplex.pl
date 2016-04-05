#!/usr/bin/env perl
die "not yet ready";
use strict;
use Getopt::Std;
use Text::Levenshtein qw(distance);
use FileHandle;

use vars qw($opt_h $opt_b $opt_p $opt_o);

my $Usage="... | $0 -b barcodes.txt -p outputprefix -o outputdir";

if ( !getopts("b:p:o:h") || $opt_h ) {
    die $Usage; 
}

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
    my $file=sprintf("%s.fastq", $lib);
    $file="$opt_p$file" if $opt_p;
    $file="$opt_o/$file" if $opt_o;
    $file = FileHandle->new("> $file") or die "library $lib, file $file: $!";
    warn "creating file $file ...\n";
    $fhs->{$lib}=$file;
  }
  $fhs;
}

sub close_outfiles {
  my($fhs)=@_;
  for my $lib (key %$fhs) {
    close($file->{$lib}) or die "could not close demultiplexed file for library $lib; investigate";
  }
}

my $codes = readbarcodes($opt_b);

my @reads;
## $read[0]=open_infile($opt_1);
## $read[1]=open_infile($opt_2);
my @codes=(values %$codes, 'unknown');
my $filehandles=open_files(@codes);      # code maps barcode to lib name. key is e.g. AGTCAA, value is e.g. M12, output is M12.fastq

my $nexact=0;
my $onemismatch=0;
my $nunknown=0;

RECORD:
while(1) { 
  my $record=<>;
  ### e.g.:  ^@NS500413:172:HVFHWBGXX:1:11101:4639:1062 1:N:0:CCGTCCAT$
  my ($code)=(split(':', $record))[-1];
  $record .= <>; # sequence line
  $record .= <>; # '+'
  $record . = <>; # quality line
  
  my $lib=$codes->{$code};
  
  if ( $lib) { 
    $nexact++;
  } elsif() { 
    $lib=rescue($code);
  } else {
    $lib='unknown';
  }
  print $filehandles->{$lib} $record;
}                                       # RECORD
close_files($filehandles);
