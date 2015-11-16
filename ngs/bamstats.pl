#!/usr/bin/perl

use strict;
use Number::Format;

my $usage="
Usage: $0 somefile.bam

Prints total reads, mapped reads, unmapped reads, percentage mapped
(tab-delimited).  It uses samtools idxstats, which requires the bam file
to have been indexed (i.e.  there should be a somefile.bam.bai in the same
directory).

See also bbcfutils::bamstat, which gives a bit more info

";

die $usage unless @ARGV ==1;


my $fmt=new Number::Format(-thousands_sep => ',');
 
sub commafy { $fmt->format_number($_[0]); }

sub read_idxstats {
    my($file)=@_;
    my $tot_mapped=0;
    my $tot_unmapped;

    my $cmd = "samtools idxstats  $file "; 
    open(BAM,  "$cmd | ") or die "$cmd: $!";

    while(<BAM>) { 
        my ($chrom, $reflen, $mapped, $unmapped)=split("\t");
        $tot_mapped += $mapped;
        $tot_unmapped += $unmapped;
    }
    close(BAM) or die "Could not close $file";
    my $tot=$tot_mapped+$tot_unmapped;
    ( map {commafy $_} ($tot, $tot_mapped, $tot_unmapped), 
      sprintf("%.1f",100*$tot_mapped/$tot));
}

print "#total\tmapped\tunmapped\t%mapped\n";
print join("\t", read_idxstats($ARGV[0])). "\n";
