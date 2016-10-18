package mismatch;

### Usage: see demultiplex.pl

use strict;

use Math::Combinatorics;
use Regexp::Optimizer;

sub readbarcodes_mixedcase {
  ## utility function to read barcodes, returns hash with eg. $barcodes->{'AGCGtT') => 'M3' }
  ## Note that lower case letters (used for disallowing specific mismatches) are
  ## still there (and won't match actual barcodes).
  my ($file)=@_;
  my $barcodeids={};
  my $barcodes_mixedcase = {};
  my $uppercase_codes={};

  open(FILE, "$file") or die "Barcode '$file': $!";
LINE:
  while(<FILE>) {
    s/[\n\r]*$//g;
    s/#.*//;
    next LINE unless $_;
    my ($barcodeid, $code)=split(' ');  # e.g. 'G7 \t CCAACAAT'
    if( $code =~ /[a-z]/) { 
      warn "barcode $barcodeid contains lower case letters, these will be uppercased and will not be allowed to mismatch";
    }
    die "Barcode id '$barcodeid' not unique" if $barcodeids->{$barcodeid}++;
    die "Barcode '$code' not unique" if $uppercase_codes->{"\U$code"}++;
    $barcodes_mixedcase->{$code}=$barcodeid;
  }                                     # while LINE
  close(FILE);
  $barcodes_mixedcase;
}                                       # readbarcodes_mixedcase

sub mixedcase2upper { 
  ## utility function to convert the mixed case hash (which is used for the mismatch regular expressions) to an uppercased
  ## hash
  my ($mixed) = @_;
  my $barcodes={};
  for my $code ( keys %$mixed ) { $barcodes->{"\U$code"}=$mixed->{$code}}
  $barcodes;
}

sub convert2mismatchREs {
## takes hash with barcodes (e.g. $h->{'AGCGtT') => 'M3' ) and returns
## e.g. $h->{'AGCGtT') => REGEXP(0x25a7788) The hash returned contains,
## per barcode, one regexp representing all possible mismatches of that
## barcode.  In the values (i.e. regexps), lowercase letters (if any) are
## uppercased and the regexp does not allow these letters to mismatch.

  my $args = ref $_[0] eq 'HASH' ? shift : {@_}; # args: barcodes, allowed
  my $o=Regexp::Optimizer->new;

  my $mm_REs={};
  for my $code (keys %{$args->{barcodes}}) {
    my @res= _getmismatch_REs($code, $args->{allowed_mismatches}); # empty if allowed_mismatches==0
    my $r='^'.join("|", @res).'$';
    $r=$o->optimize(qr/$r/);
    $mm_REs->{$code}= $r;               # just one big regexp!
  }                                     # for $code
  $mm_REs;  
}                                       # convert2mismatchREs

sub rescue { 
  ### return the barcode without mismatches (not its id!)
  my($foundcode, $mm_REs)=@_;

  foreach my $code (keys %$mm_REs) {
    my $re=$mm_REs->{$code};
    return  $code if $foundcode =~ $re;
  }
  return undef;
}                                       # rescue

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



1;
