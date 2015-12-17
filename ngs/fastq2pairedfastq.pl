#!/bin/env perl

use strict;
use Getopt::Long;
use Number::Format;
use FileHandle;

use vars qw($help @outfiles @outs);

my $dflt_unpaired = "unpaired.fastq.gz";

my $Usage = q{
  Usage:  fastq2pairedfastq.pl  < pairedendreads.fastq  --out1 out_R1.fastq[.gz] --out2 out_R2.fastsq[.gz]

This script is meant as a postprocessing step when converting a paired-end
read BAM file back to FASTQ (e.g. if you lost it). In the first step, the
BAM file must be sorted by name (using e.g. samtools sort -n). This
name-sorted BAM file is converted to FASTQ using bamToFastq, but
bamToFastq gets confused by multimapping reads. Therefore, use bamToFastq
in single-end mode, and recover the pairs from the resulting FASTQ using
the current script. Duplicates due to multi-mapping reads will be
discarded.  The matching is done on the basis of the read IDs. Only IDs
that have exactly two different sequences (although their numbers do not
matter) are output as matching mate pairs. Unpaired reads are warned
about, and output to a third file.

If the  output filenames end in '.gz', the result will be gzipped on the fly. 

Options:

  --unpaired FILE   Write unpaired reads here (default is $dflt_unpaired). 

Maybe also have a look at 
https://github.com/enormandeau/Scripts/blob/master/fastqCombinePairedEnd.py

};

$outfiles[2]=$dflt_unpaired;

die $Usage unless &GetOptions(
  'help' => \$help,
  'out1=s' => \$outfiles[0],
  'out2=s' => \$outfiles[1],
  'unpaired=s' => \$outfiles[2],
    );

die $Usage if ($help || !$outfiles[0] || !$outfiles[1]);

my $fmt=new Number::Format(-thousands_sep => ',');
 
sub commafy {
  $fmt->format_number($_[0]);
}

sub do_output { 
  ## also return 1 for unpaired, 2 for paired
  my($id, $values)=@_;
  
  my @vals=keys %$values;

  if ( int(@vals)  == 1 ) { 
    chomp($id);
    warn "Unpaired: $id\n";
    $outs[2]->print($id."\n" . $vals[0]);
    return 1;
  }

  if ( int(@vals)  > 2 ) { 
    chomp($id);
    warn "*** found more than 2 reads for '$id', not supposed to happen!!!\n**** ----\n";
    for my $v (@vals) { 
      print STDERR "    $v\n--------\n";
    }
    print STDERR "**** --------\n";
    return 3;
  }
  
  for my $i (0 .. 1) { 
    $outs[$i]->print($id, $vals[$i]);
  }
  return 2;
}                                       # sub do_output

for my $i ( 0 .. 2) {                   # open output files
  if( $outfiles[$i] =~ /\.gz$/)  { 
      $outs[$i] = FileHandle->new(" | gzip > $outfiles[$i]") or die "$outfiles[$i]: $!";
  } else { 
      $outs[$i] = FileHandle->new("> $outfiles[$i] ") or die "$outfiles[$i]: $!";
  }
}

my ($id, $values, $previd, $seq, $plus, $qual);
my ($nfrags, $ninput, $nunpaired, $nerror);
$values={};

while(1 &&  ! eof(STDIN) ) { 
  $id= <>;

  if ($previd && ($id ne $previd)) {    # found new id, output previous one
##    die "input not sorted by name, at line $.: $previd should come later than $id\n"
##        unless (( $previd cmp $id ) == -1);
## samtools sort -n uses a weird sorting rule that breaks this check. Ignore.
    my $res=do_output($previd, $values);

    $nunpaired += ($res == 1);
    $nfrags += ($res == 2);
    $nerror += ($res == 3);

    $values = {};
  }
  $previd=$id;

  $seq= <>; 
  $plus = <>; warn "incomplete record at end," if eof(STDIN);
  $qual = <>;

  my $value="$seq$plus$qual";
  $ninput++;
  warn "Read ".commafy($ninput)." reads\n" if ($ninput % 100_000) == 0;

  $values->{$value}++;
}

my $res=do_output($id, $values);                # every read is sacred (cue Meaning of Life)
$nunpaired += ($res == 1);
$nfrags += ($res == 2);
$nerror += ($res == 3);

warn "Done reading " . commafy($ninput) . " reads\n";
warn "Wrote " . commafy($nunpaired) . " unpaired reads, ignored $nerror ids having errors (see log)\n";
warn "Wrote " . commafy($nfrags*2) . " paired reads, i.e. ".commafy($nfrags)." fragments\n";

foreach my $i (0 .. 2) {                # closing needed for gzip pipes
  $outs[$i]->close() or die "$outfiles[$i]: $!";
}
