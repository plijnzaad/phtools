package mismatch;

### Usage: see demultiplex.pl

use strict;

use Math::Combinatorics;
use Regexp::Optimizer;


sub readbarcodes {
  ## utility function to read barcodes
  ## returns hash with $barcodes->{'AGCGTT') => 'M3' }
  ## note that lower case letters (used for disallowing specific mismatches) are
  ## still there, and won't match actual barcodes!
  my ($file)=@_;
  my $barcodeids={};
  my $barcodes = {};
  my $uppercase_codes={};

  open(FILE, "$file") or die "Barcode '$file': $!";
LINE:
  while(<FILE>) {
    s/[\n\r]*$//g;
    s/#.*//;
    next LINE unless $_;
    my ($barcodeid, $code)=split(' ');            # e.g. 'G7 \t CCAACAAT'
    if( $code =~ /[a-z]/) { 
      warn "barcode $barcodeid contains lower case letters, these will be uppercased and will not be allowed to mismatch";
    }
    die "Barcode id '$barcodeid' not unique" if $barcodeids->{$barcodeid}++;
    die "Barcode '$code' not unique" if $uppercase_codes->{"\U$code"}++;
    $barcodes->{$code}=$barcodeid;
  }                                     # while LINE
  close(FILE);
  $barcodes;
}                                       # readbarcodes

sub convert2mismatchREs {
## takes hash with barcodes (e.g. $h->{'AGCGtT') => 'M3' )  and returns e.g. $h->{'AGCGtT') =>  REGEXP(0x25a7788)
## The hash returned contains, per barcode, one regexp representing all possible mismatches of that barcode.
## Lowercase letters are uppercased and the regexp does not allow these letters to mismatch.
  my $args = ref $_[0] eq 'HASH' ? shift : {@_}; # args: barcodes, allowed
  my $o=Regexp::Optimizer->new;

  my $mm_REs={};
  for my $code (keys %{$args->{barcodes}}) {
    my @res= _getmismatch_REs($code, $args->{allowed_mismatches}); # empty if allowed_mismatches==0
    my $r='^'.join("|", @res).'$';
    $r=$o->optimize(qr/$r/);
    $mm_REs->{$code}= $r;         # just one big regexp!
  }                               # for $code
  $mm_REs;  
}                                       # convert2mismatchREs

sub _getmismatch_REs {
  ## for one barcode, set up the regular expressions that allows mismatches
  my($code, $max_mm)=@_;

  return () if ! $max_mm;

  my @fixed = ();
  if ($code =~ /[a-z]/)  {
    my $fixed= $code;
    $fixed =~ s/[a-z]/!/g;
    @fixed = split(//, $fixed);
    $code = "\U$code";
  }

  my @mmcodes=();
  my(@code)=split(//, $code);

  ## set up array of arrays with '.' where to do the replacements:
  for(my $i=0; $i<$max_mm; $i++) { 
    my @combs = combine(($i+1), 0..$#code);
  COMB:
    foreach my $comb ( @combs ) { 
      my @mm=@code;
      @mm[ @$comb ] = split(//, '.' x int(@$comb) ); # yay, splicing
      my $mm=join("", @mm);
      for my $i (0 .. $#fixed) { 
        if ($fixed[$i] eq '!' && $mm[$i] eq '.') { 
          warn "regexp $mm conflicts with a fixed position, skipped\n";
          next COMB;
        }
      }
      push(@mmcodes, $mm);
    }
  }
  @mmcodes;
}                                       # _getmismatch_REs

sub rescue { 
  my($foundcode, $barcodes, $mm_REs)=@_;

  foreach my $code (keys %$barcodes) {
    my $re=$mm_REs->{$code};
    return  $barcodes->{$code} if $foundcode =~ $re;
  }
  return undef;
}                                       # rescue


1;
