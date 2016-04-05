#!/usr/bin/env perl

use strict;
use Getopt::Std;
use Text::Levenshtein qw(distance);
use FileHandle;

use vars qw($opt_h $opt_b $opt_1 $opt_2 $opt_p $opt_o);

my $Usage="$0 -b barcodes.txt -1 all_R1.fastq[.gz] -2 all_R2.fastq[.gz] -p outputprefix -o outputdir";

if ( !getopts("b:1:2:h") || $opt_h ) {
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

sub close_fies {
  my($fhs)=@_;
  for my $lib (key %$fhs) {
    close($file->{$lib}) or die "could not close demultiplexed file for library $lib; investigate";
  }
}

my $codes = readbarcodes($opt_b);

my @reads;
$read[0]=open_infile($opt_1);
$read[1]=open_infile($opt_2);

my $filehandles=open_files(values %$codes);      # code maps barcode to lib name. key is e.g. AGTCAA, value is e.g. M12, output is M12.fastq

my $nexact=0;
my $onemismatch=0;
my $nlost=0;

RECORD:
while(1) { 
  for() ...

  my $idline=<>;
  my $seqline=<>;
  my $dummy=<>;
  my $qualline=<>;

  my $code= ???;
  my $lib=$codes->{$code};

  if ( $lib) { 
    $nexact++;
    print $filehandles->{$lib} $idline . $seqline . "+\n" . $qualline;
    next RECORD;
    DO WE relate this TO GET THIS FROM THE OTHER READ???
  } elsif { 
    $code=rescue();
}

}
close_files($filehandles);

OR USE sabre:     http://www.vcru.wisc.edu/simonlab/bioinformatics/programs/sabre/README.md.txt

(https://github.com/najoshi/sabre ? ) 

OR Casava ? (check e.g. http://seqanswers.com/forums/showthread.php?t=18736 )

OR fastx_barcode_splitter.pl (no, not paired-end)

