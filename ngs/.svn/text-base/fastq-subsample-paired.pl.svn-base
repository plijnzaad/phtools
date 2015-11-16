#!/bin/env perl

### NOTE: also check samtools view -s 322.15 -b file.bam > random_15%_of_file.bam

use strict;

use Getopt::Long;
use FileHandle;
use Number::Format;
## use Math::Random;                       # broken! rand() is much better

use vars qw($help $in1 $in2 $out1 $out2 $seed $perc);

my $Usage = q{
Usage: 

    fastq-subsample-paired.pl --perc percentage --in1 input_R1.fastq --in2 input_R2.fastq --out1 subset_R1.fastq --out2 subset_R2.fastq

Does approximate subsampling in one pass for paired reads. Ordering will
be the same as that of the input files. Note that the ordering in both
input files must be identical; an error is thrown if this is not the case.

Input files can be gzip-compressed (judged by their names ending in '.gz'); likewise,
output files ending in '.gz' will be automatigically gzip-compressed on the fly.

For single-end reads use fastq-subsample-single.pl

Options: 
  --seed  N   Seed the random number generator with N,N (for reproduceability purposes)
};


die $Usage unless &GetOptions('help' => \$help,
                              'in1=s'=> \$in1,
                              'in2=s'=> \$in2,
                              'out1=s'=> \$out1,
                              'out2=s'=> \$out2,
                              'seed=i'=> \$seed,
                              'perc=f'=> \$perc);

die $Usage if ($help or !$perc);

die "\n *** Need total of four file names ***\n$Usage" 
    if (!! $in1 + !! $in2  + !!$out1 + !!$out2) < 4;

if ($seed) {
#   random_set_seed(($seed,$seed));
  srand($seed);
}
## warn "random seed was ", join(" ", random_get_seed()), "\n";

my $frac= 1 - ($perc/100);

my $fmt=new Number::Format(-thousands_sep => ',');
 
sub commafy {   $fmt->format_number($_[0]); }

my (@ins, @outs, @innames, @outnames, $id1, $id2, $rest);

$innames[0] = $in1;
$innames[1] = $in2;
$outnames[0] = $out1;
$outnames[1] = $out2;

for my $i (0 .. 1) { 
  if( $innames[$i] =~ /\.gz$/)  { 
      $ins[$i] = FileHandle->new("zcat $innames[$i] | ") or die "$innames[$i]: $!";
  } else { 
      $ins[$i] = FileHandle->new("< $innames[$i] ") or die "$innames[$i]: $!";
  }

  if( $outnames[$i] =~ /\.gz$/)  { 
      $outs[$i] = FileHandle->new(" | gzip > $outnames[$i]") or die "$outnames[$i]: $!";
  } else { 
      $outs[$i] = FileHandle->new("> $outnames[$i] ") or die "$outnames[$i]: $!";
  }

}

my (@id_lines, @seqs, @pluses, @quals, $nseqs, $nselected);

while(! eof($ins[0]) && !eof($ins[1]) ) { 

  foreach my $i (0..1) { 
    $id_lines[$i] = $ins[$i]->getline();
    die "incomplete block for file $innames[$i]" if eof($ins[$i]);
    $seqs[$i]= $ins[$i]->getline();
    die "incomplete block for file $innames[$i]" if eof($ins[$i]);
    $pluses[$i] = $ins[$i]->getline();
    die "incomplete block for file $innames[$i]" if eof($ins[$i]);
    $quals[$i] = $ins[$i]->getline();
  }
  ($id1, $rest) = split(' ', $id_lines[0], 2);
  ($id2, $rest) = split(' ', $id_lines[1], 2);

  die "read IDs not identical: $id1 vs $id2 at $in1 line "
      .$ins[0]->input_line_number(). " / $in2 line ".$ins[1]->input_line_number()."\n"
      unless $id1 eq $id2;

  $nseqs++;

##  if( random_uniform() > $frac) { 
  if( rand(1) > $frac) { 
    $nselected++;
    foreach my $i (0 .. 1) { 
      $outs[$i]->print($id_lines[$i] . $seqs[$i] . $pluses[$i] . $quals[$i]);
    }
  }
}                                       # while

foreach my $i (0 .. 1) { 
  $outs[$i]->close() or die "$outnames[$i]: $!";
}

warn(sprintf("Selected %s out of %s reads (%4.2f%%)\n", 
             commafy($nselected), commafy($nseqs), 100*$nselected/$nseqs));
