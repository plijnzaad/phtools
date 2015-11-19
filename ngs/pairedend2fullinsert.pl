#!/usr/bin/perl 
my $version = '0.0';

my $scriptname="pairedend2fullinsert.pl";

use Getopt::Long;
use Number::Format;
my $fmt=new Number::Format(-thousands_sep => ',');
sub commafy {  $fmt->format_number($_[0]); }

my $usage="

Usage: Convert paired-end reads to full-length inserts for the purpose of
  visualizing the coverage across their length (typically ChIP-Seq or
  MNaseSEQ data). Input and output must be/is in SAM format (with header),
  and is read from stdin/written to stdout.  The resulting CIGAR string is
  <insertlength>M; the template length is zero, and the reesulting SEQs
  are all N's.

  Be aware of the way the mapping was done. If e.g. the mapper was told to
  discard the first 15 bp from each read (to improve the mapping), you
  prolly want to adjust them back to their original length (this is what
  the --untrim option is for).

  The --minlen and --maxlen arguments apply to the length after any untrimming.

Options:
  --untrim <number>  Extend all reads on their 5'-side (needed if they were trimmed during mapping)
  --minlen <number>  Skip reads where the fragment length (after any untrimming) is less than this (default: 0)
  --maxlen <number>  Skip reads where the fragment length (after any untrimming) is greater than this (default: 1e6)
  --chrom_sizes <file> tab-delimited file with ^chromosome_name\\tchromosome_length$ (needed if not in header of SAM file, or if those are wrong)

Written by <plijnzaad\@gmail.com>

";

use strict;

my $help=0;
my $chrom_sizes=undef;
my $minlen=0;
my $maxlen=1000000;

my $untrim=0;
my $ninserts=0;

my $pg_printed=0;

my $chromos=undef;

my ($nfirst, $nsecond, $too_short, $too_long, $unmapped,
    $skipped_left, $skipped_right)=(0,0,0,0,0,0,0);

my @argv_copy=@ARGV;                    # eaten by GetOptions
die $usage if  GetOptions('help'=> \$help,
                          'chrom_sizes|c=s' => \$chrom_sizes,
                          'untrim|u=i' => \$untrim,
                          'minlen|G=i' => \$minlen,
                          'maxlen|L=i' => \$maxlen
    ) ==0 || $help;

my $cmdline= "$0 " . join(" ", @argv_copy);

if ($chrom_sizes) {
  $chromos = read_chromo_sizes( $chrom_sizes );
} else {    
  $chromos = {};                        # read during parsing
}


my $single_read_mask= 0x10 | 0xf00;
### (any bits not in this refer to paired end reads, so must be unset)

LINE:
    while(<>) { 
      if (/^@/) { 
        print;
        if (!$chrom_sizes && /^\@SQ\s+SN:(\S+)\s+LN:(\d+)/ )  { # length record
          $chromos->{$1}=$2;
        }
        next LINE;
      } else {                          # append our own line
        print "\@PG\tID:$scriptname\tPN:$scriptname\tVN:$version CL:\"$cmdline\"\n"
            unless $pg_printed++;
      }
      s/[\n\r]*$//;
      
      my($qname,$flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen,
         $seq, $qual, @optionals)=split("\t", $_);
      next if $rname eq  '*';
      my $chr_length=$chromos->{$rname};
      die "Unknown chromosome or chromosome length: chr='$rname', input line $.
 (SAM input file must contain the sequence lengths, e.g. from samtools view -h; otherwise, supply using --chrom_sizes option)" 
      unless $chr_length;

      if(!$tlen) { 
        die "$0: SAM file contains no template lengths (is this paired-end data?)\n";
      }

      my $readlen=length($seq);
      if ($flag & 0x4) {
        $unmapped++;
        next LINE;
      }

      $pos = $pos - $untrim;
      if ($pos <= 0) {
        $skipped_left++;
        next LINE;
        
      }
      my $newlen=abs($tlen) + 2*$untrim;
      if ($pos + $newlen -1 > $chr_length) { 
        $skipped_right++;
        next LINE;
      }
      
      if ( $newlen < $minlen ) {
        $too_short++;
        next LINE;
      }
      if ($newlen > $maxlen ) {
        $too_long++;
        next LINE;
      }

      if($tlen < 0) {                   # second of a mate
        $nsecond++;
        next LINE;
      }
      
      $nfirst++;                        # first of a mate
      $rnext='*';                       # i.e. not available
      $pnext=0;
      $flag = $flag & $single_read_mask;
      $cigar=sprintf('%dM', $newlen);
      $seq=  'N' x $newlen;               # whole sequence is just N's
      $qual='*';
      $ninserts++;
      my @fields=($qname,$flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext,
                  0, $seq, $qual, @optionals);
      print join("\t", @fields) . "\n";
}                                       # LINE

if ($nfirst ne $nsecond) { 
  warn "Expected equal number of first and second mates, instead found "
      . commafy($nfirst) . " first mates but ".commafy($nsecond) ." second mates\n";
}

die "No inserts where written" unless $ninserts;
warn commafy($unmapped) . " reads were unmapped and skipped\n";
warn commafy($skipped_left) . " reads dropped off the left side, ". commafy($skipped_right) . " off the right side of the chromosome\n";
## (non of the above should happen, really)
warn "Wrote ". commafy($ninserts) . " inserts\n";
warn "Dropped ". commafy($too_short) . " fragments because too short, ". commafy($too_long)  ." because too long\n";

sub read_chromo_sizes {
    my($file)=@_;

    my $table={};
    
    open(FILE, $file) or die "Could not read file with chromosome sizes: '$file'\n";
    while(<FILE>) { 
        s/[\n\r]//g;
        my ($chr, $len)=split("\t");
        if (!$chr || $len !~ /^\d+$/) { 
            die "Wrong format for chromosome sizes: chr='$chr', length='$len'\n";
        }
        $table->{$chr}=$len;
    }
    close(FILE);
    $table;
}                                       # read_chromo_sizes
