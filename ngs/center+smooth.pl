#!/usr/bin/perl 
my $version = '0.0';

use Getopt::Long;
use Number::Format;
my $fmt=new Number::Format(-thousands_sep => ',');
sub commafy {
  $fmt->format_number($_[0]);
}

my $usage='
Usage: center+smooth.pl --type [paired|single] [ --shift NUMBER ] [ --smooth NUMBER ] [ options ] 

  The script is aimed at MNase-seq data (used to map nucleosome positions)
  and is meant to get sharper peaks at the likely position of the
  nucleosome dyad. It reads a SAM file from stdin, and writes (to stdout)
  a SAM file, all shifted into their 3\'-direction. For paired-end (PE) reads,
  the amount of shifting is determined by the template length (field
  9) found in the SAM file.  It can also be specified manually using the
  --shift paramter, typically when dealing with single-end data. For PE reads,
  only the first of the two reads is kept, and the output data is *single ended*.

  Make sure that the fragments in the SAM file correspond to single
  nucleosomes, otherwise shifting by half the template length does not
  make sense (see also the --minlen and --maxlen options).

  The reads are made into "single-basepair reads". The reason for this is
  that basepair 2, 3, etc. are not additional independent evidence of the
  nucleosome dyad being located at the 5\'-end of the current read, and
  including additional basepairs only smooths the data, thus decreasing
  the resolution. If this smoothing is desired, use the --smooth option
  (31 is a a reasonable value). After centering and/or smoothing, the read
  coordinates are changed and consistent but not sorted
  anymore. The resulting read "sequence" consists of only \'N\'s, the
  CIGAR string of only \'M\'s.

  The program is typically part of a pipeline that converts from BAM and
  back again, as

      prefix=/dev/shm/$USER.$$.center
      samtools view -h reads.bam | center+smooth.pl --smooth 31 |\
          samtools sort -T$prefix -Obam > centered_reads.bam

  Be aware of the way the mapping was done. If e.g. the first 15 bp from
  the reads were discarded to improve the mapping, the shifting should not
  be by 73 bp, but by 73-15 = 58 bp (this point is of course moot for
  paired-end data, where the shift is best inferred from the template length).
  Similarly, if we would be interested in nucleosome boundaries
  (separately for the left and right side), one should in this case shift
  by -15, to pretend that the trimmed reads _did_ actually begin at the
  boundary.

  This caution also applies to the --min and --max length: the defaults
  correspond to single-nucleosome fragments provided the mapping was done
  without any trimming.

  The script refuses to run on data that was previously output by it, as it
  leads to intractable shifting errors (the reason is that the script
  operates essentially on the first basepair of each read, not on its
  middle).

Options:

  --type  paired|single Indicate if the SAM file contains single-end or 
                        paired-end reads
  --shift <number>    Shift all reads by this amount in their own
                      3\'-direction.  Can be negative. If not specified,
                      shifting is by half the fragmentlength found in
                      field 9 of the SAM file. This only makes sense for
                      paired-end reads, so the program warns about reads
                      having length zero, and skips them.

  --minlen <number> Skip reads where the fragment length is less than this
                    (default: 100)

  --maxlen <number> Skip reads where the fragment length is greater than
                    this (default: 200)

  --smooth <number> Number of basepairs to smooth over. Must be uneven
                    (since "middle basepair" +- n gives 2n+1)

  --chrom_sizes <file> tab-delimited file containing
    ^chromosome_name\\tchromosome_length$ (needed if not in header of SAM
                       file, or if those are wrong)

  --strict Skip shifted reads (enlarged or not) that don\'t fall entirely
           within the chromosome. The default is to shorten such reads so
           that they do again.

For more speed, see bbcfutils::bam2wig.

Note that as a result of the change in read length, the coverage will also
change. I.e., if single-end reads of 55 bp are centered to 1 bp, then
smoothed to 31, the coverage per position will be 56% of the original.

';

my $help=0;
my $shift=undef;                        # meaning: automatic
my $smooth=1;                           # i.e. none
my $chrom_sizes=undef;
my $minlen=100;
my $maxlen= 200;                    # i.e. at most one nucleosome!
my $strict=undef;
my $auto=1;
my $seqtype=undef;
## my $nodrop=undef;

my @argv_copy=@ARGV;                    # eaten by GetOptions
die $usage if  GetOptions('help'=> \$help,
                          'type=s' => \$seqtype,
##                          'nodrop' => \$nodrop,
                          'shift=i' => \$shift,
                          'chrom_sizes=s' => \$chrom_sizes,
                          'minlen=i' => \$minlen,
                          'maxlen=i' => \$maxlen,
                          'shift=i' => \$shift,
                          'smooth=i' => \$smooth,
                          'strict' => \$strict,
    ) ==0 || $help;

my $cmdline= "$0 " . join(" ", @argv_copy);

use strict;

die "--type argument required, must be 'paired' or 'single' " 
    unless ($seqtype =~ /single/i || $seqtype =~ /paired/i ) ;

my $single= ($seqtype =~ /single/i);

if ($single) {
  die "--shift argument is required for single-end reads" unless  $shift;
###  die "--nodrop option only valid when using --type paired" if $nodrop;
}

my $pg_printed=0;

my $chromos=undef;

if ($chrom_sizes) {
  $chromos = read_chromo_sizes( $chrom_sizes );
} else {    
  $chromos = {};                        # read during parsing
}

my $nreads=0;
my ($trimmed_left, $trimmed_right, $skipped_left, $skipped_right)=(0,0,0,0);
my ($unmapped, $too_short, $too_long, $no_length, $mate2dropped)=(0,0,0,0,0);

die "smoothing window must be uneven" unless ($smooth % 2);

my $halfsmooth= int( ($smooth -1) /2);

LINE:
    while(<>) {
      if (/^@/) { 
        if ( /PG.*center\+smooth.pl/ ) {
          die "Refusing to run center+smooth.pl on data that was previously centered and/or smoothed\n";
        }
        print;
        if (!$chrom_sizes && /^\@SQ\s+SN:(\S+)\s+LN:(\d+)/ )  { # length record
          $chromos->{$1}=$2;
        }
        next LINE;
      } else { 
        print "\@PG\tID:center+smooth.pl\tPN:center+smooth.pl\tVN:$version CL:\"$cmdline\"\n"
            unless $pg_printed++;
      }
      s/[\n\r]*$//;
      
      my($qname,$flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen,
         $seq, $qual, @optionals)=split("\t", $_);
      if ($rname eq  '*') { 
        $unmapped++;
        next LINE;
      }
      my $chr_length=$chromos->{$rname};
      die "Unknown chromosome or chromosome length: chr='$rname', input line $.
 (SAM input file must contain the sequence lengths; otherwise, supply using --chrom_sizes option)" 
      unless $chr_length;

      if(!$tlen && !$shift ) { 
        $no_length++;
        next LINE;
      }
      my $tlen_sign = ($tlen <=> 0);
      $tlen=abs($tlen);

      my $readlen=length($seq); ### cannot trust this if there are indels!
      if ($flag & (0x4 | 0x8) ) {        # note: they may have the $rname of their mate, so have yet been skipped
        $unmapped++;
        next LINE;
      }
      
      my $s;
      if ($single) {
        die "read is paired, but you specified --type single" if ($flag & 0x1);
        $s=$shift;
      } else {
        if ($flag & 0x80) {             # looking at R2
          $mate2dropped++;
          ## die "not implemented properly: coordinates will be wrong for 2nd mate if there are indels ..." if $nodrop;
          next LINE;
          ## note: R1 will be turned into a SE read!
        }
        if ($tlen < $minlen ) {
          $too_short++;
          next LINE;
        }
        if ( $tlen > $maxlen ) {
          $too_long++;
          next LINE;
        }
        $s=int( ($tlen-1)/2);     # exact for uneven seqs
        $s = $s + (unpack("%B*",$seq)%2) unless ($tlen % 2); # quasi-randomly  add 1
        # (the unpack expression is quasi-uniform hash based on sequence content)
        $s = $shift if $shift;         # allow override
      }

      $rnext='*';
      $pnext=0;
      $flag = ($flag & ~(1<<0));        # make R1 into a single-end read

      my $reverse_strand = ($flag & 0x10);

      if ($reverse_strand) { 
        die "Not yet properly implemented (overlooking indels)" if $single;
        $pos = $pos + ($readlen -1) - $s - $halfsmooth;
      } else {
        $pos = $pos                 + $s - $halfsmooth;
      }
      my $end=$pos + $smooth -1;
      
      ## following cannot be salvaged by trimming, kill them:
      if ( $end < 1)  { 
        $skipped_left++;
        next LINE;
      }
      if ( $pos > $chr_length ) { 
        $skipped_right++;
        next LINE;
      }
      # next ones straddle start or end, so trimming will salvage them:
      my $newlen=$smooth;
      if ($pos < 1 ) { 
        if ($strict) { 
          $skipped_left++;
          next LINE;
        }
        $trimmed_left++;
        $newlen=$smooth + ($pos -1);
        $pos=1;
      }
      if ( $end > $chr_length) {
        if ($strict) { 
          $skipped_right++;
          next LINE;
        }
        $trimmed_right++;
        $newlen=$smooth - ($end - $chr_length);
      }
      $cigar=sprintf('%dM', $newlen);
      $seq=  'N' x $newlen;               # whole sequence is just N's
      $tlen=0; # since everything now a single-end read!
      $qual='*';
      
      my @fields=($qname,$flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext,
                  $tlen_sign*$tlen, $seq, $qual, @optionals);
      
      print join("\t", @fields) . "\n";
      $nreads++;
}                                       # LINE

warn "Shifted and output ". commafy($nreads) . " reads, skipped ". commafy($unmapped). " unmapped reads\n";
warn commafy($skipped_left) . " reads skipped on the left side, ". commafy($skipped_right) . " on the right side of the chromosome\n";
warn commafy($trimmed_left) . " reads trimmed on the left side, " . commafy($trimmed_right) . " on the right side of the chromosome\n";
if(!$single) { ## paired-end only:
  warn "Dropped ". commafy($mate2dropped) . " mate2 lines because uninformative\n";
  warn "Dropped ". commafy($too_short) . " fragments because too short, ". commafy($too_long)  ." because to long\n";
  warn "Found " . commafy($no_length) . " unpaired reads, skipped because no --shift was specified\n";
}
die "No reads were output!" unless $nreads > 0;

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
