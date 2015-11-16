#!/usr/bin/perl 
my $version = 'unknown';
## (leave as is, id is automagically adjusted by svn)

## This script is based on center+smooth.pl revision -r1019 (24 aug 2015)

my $scriptname="untrim.pl";
use Getopt::Long;

my $usage=q{
Usage: $scriptname [ --untrim5 NUMBER  ]  [ --untrim3 NUMBER ]

  Reads a SAM file from stdin, and writes (to stdout) a SAM file, all
  'untrimmed' to undo the effect of 5'- or 3'-trimming done by e.g.
  bowtie. After untrimming the read coordinates are changed and consistent
  (although not sorted anymore). The resulting read sequence consists of
  only 'N's, the CIGAR string of only 'M's.

  The program is typically part of a
  pipeline that converts from BAM and back again, as

      prefix=/dev/shm/$USER.$$.center
      samtools view -h reads.bam | $scriptname --untrim5 15 |\\
          samtools sort -T$prefix -Obam > centered_reads.bam

Options:

  --untrim5 <number>  Extend all reads on their 5'-side. Can be negative (in which case the reads become shorter at their 5'-side)
  --untrim3 <number>  Extend all reads on their 3'-side. Can be negative (in which case the reads become shorter at their 3'-side)
  --chrom_sizes <file> tab-delimited file with ^chromosome_name\\tchromosome_length$ (needed if not in header of SAM file, or if those are wrong)
  --strict            Skip shifted reads (enlarged or not) that don't fall entirely within the chromosome. The default is to shorten such reads so that they do again.

For more speed, see bbcfutils::bam2wig.
};

my $help=0;
my $untrim5=0;
my $untrim3=0;
my $chrom_sizes=undef;
my $strict=undef;
my @argv_copy=@ARGV;                    # eaten by GetOptions
die $usage if  GetOptions('help'=> \$help,
                          'chrom_sizes|c=s' => \$chrom_sizes,
                          'untrim5|5=i' => \$untrim5,
                          'untrim3|3=i' => \$untrim3,
                          'strict' => \$strict,
    ) ==0 || $help;

my $cmdline= "$0 " . join(" ", @argv_copy);

use strict;

my $chromos=undef;

if ($chrom_sizes) {
  $chromos = read_chromo_sizes( $chrom_sizes );
} else {    
  $chromos = {};                        # read during parsing
}

my $nreads=0;
my ($trimmed_left, $trimmed_right, $skipped_left, $skipped_right)=(0,0,0,0);
my $unmapped=0;

my $pg_printed=0;

LINE:
    while(<>) { 

      if (/^@/) { 
        print;
        if (! $chrom_sizes && /^\@SQ\s+SN:(\S+)\s+LN:(\d+)/ )  { # length record
          $chromos->{$1}=$2;
        }
        next LINE;
      } else { 
        print "\@PG\tID:$scriptname\tPN:$scriptname\tVN:$version CL:\"$cmdline\"\n"
            unless $pg_printed++;
      }

      s/[\n\r]*$//;

      my($qname,$flag, $rname, $start, $mapq, $cigar, $rnext, $pnext, $tlen,
         $seq, $qual, @optionals)=split("\t", $_);
      next if $rname eq  '*';
      my $chr_length=$chromos->{$rname};
      die "Unknown chromosome or chromosome length: chr='$rname', input line $. \
 (SAM input file must contain the sequence lengths; otherwise, supply using --chrom_sizes option)"  unless $chr_length;

      if ($flag & 0x4) {
        $unmapped++;
        next LINE;
      }

      my $newlen= length($seq) + $untrim5 + $untrim3;
      die "New read-length of '$qname' is less than 1: $newlen\ninput line $."
          if $newlen < 0;

      my $reverse_strand = ($flag & 0x10);
      
      if ($reverse_strand) { 
        $start = $start - $untrim3;
      } else {
        $start = $start - $untrim5;
      }
      my $end=$start+$newlen;

      ## following cannot be salvaged by trimming, kill them:
      if ( $end < 1)  { 
        $skipped_left++;
        next LINE;
      }
      if ( $start > $chr_length ) {     # should not happen ...
        $skipped_right++;
        next LINE;
      }

      # next ones straddle start or end, so trimming will salvage them:
      if ($start < 1 ) { 
        $trimmed_left++;
        next LINE if $strict;
        $newlen= $newlen + $start;
        $start=1;
      }
      if ( $end > $chr_length) {
        $trimmed_right++;
        next LINE if $strict;
        $newlen= $newlen + ($chr_length - $end);
      }

      $nreads++;
      $cigar=sprintf('%dM', $newlen);
      $seq=  'N' x $newlen;               # whole sequence is just N's
      $tlen=$newlen; # note: converting everything to single-end sequences here!
      $qual='*';

      my @fields=($qname,$flag, $rname, $start, $mapq, $cigar, $rnext, $pnext,
                  $tlen, $seq, $qual, @optionals);

      print join("\t", @fields) . "\n";
}                                       # LINE

warn "Shifted $nreads reads, skipped $unmapped unmapped reads\n";
warn "$skipped_left reads skipped on the left side, $skipped_right on the right side of the chromosome\n";
warn "$trimmed_left reads trimmed on the left side, $trimmed_right on the right side of the chromosome\n";

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
