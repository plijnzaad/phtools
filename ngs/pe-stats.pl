#!/bin/env perl

## reads samtools flagstat output and reformat it to something more
## readable, including all percentages. only tested on output from
## 'samtools flagstat', run on bowtie2 results.

## Written by plijnzaad@gmail.com

use strict;
use Number::Format;

my $version='v0.1';

sub get_line { 
  my($re)=@_;
  
  $_=<>;
  chomp;
  die "expected /$re/, got: $_ (line $.)\n" unless /$re/;
  ($1,$2);
}

###### input looks like:
## 36640494 + 0 in total (QC-passed reads + QC-failed reads) 
## 0 + 0 secondary                                           
## 0 + 0 supplementary                                       
## 0 + 0 duplicates                                          
## 33902083 + 0 mapped (92.53%:-nan%)                        
## 36640494 + 0 paired in sequencing                         
## 18320247 + 0 read1                                        
## 18320247 + 0 read2                                        
## 17801358 + 0 properly paired (48.58%:-nan%)               
## 33083344 + 0 with itself and mate mapped                  
## 818739 + 0 singletons (2.23%:-nan%)                       
## 1697054 + 0 with mate mapped to a different chr           
## 1022341 + 0 with mate mapped to a different chr (mapQ>=5) 
######

my ($total,)=get_line('^(\d+) \+ (\d+) in total');
my ($secondary,)=get_line('^(\d+) \+ (\d+) secondary'); # ignored
my ($supplementary,)=get_line('^(\d+) \+ (\d+) supplementary'); # ignored
my ($duplicates,)=get_line('^(\d+) \+ (\d+) duplicates');       # ignored
my ($mapped,)=get_line('^(\d+) \+ (\d+) mapped');
my ($pe,)=get_line('^(\d+) \+ (\d+) paired'); # ignored
my ($read1,)=get_line('^(\d+) \+ (\d+) read1'); # ignored
my ($read2,)=get_line('^(\d+) \+ (\d+) read2');
my ($properly,)=get_line('^(\d+) \+ (\d+) properly paired');
my ($paired,)=get_line('^(\d+) \+ (\d+) with itself and mate mapped');
my ($singly_unm,)=get_line('^(\d+) \+ (\d+) singletons');
my ($other_chromo,)=get_line('^(\d+) \+ (\d+) with mate mapped to a different chr');
my ($other_chromo_Q5,)=get_line('^(\d+) \+ (\d+) with mate mapped to a different chr \(mapQ>=5\)');

## Note: samtools flagstats counts most thing per mate, not per pair.
## We'll do the same here and divide later on. For singly_unmp and
## other_chromo, flagstats gives stats per single mate, so we have to
## double those (can read as: mate is part of a pair that ...)

$other_chromo *=2;
$other_chromo_Q5 *=2;
$singly_unm *= 2;

sub perc { sprintf("%6.1f%%", 100*$_[0]); }
sub out { print join("\t", @_) . "\n"; }
my $fmt=new Number::Format(-thousands_sep => ',');
sub commafy {   sprintf("%12s", $fmt->format_number($_[0])); }

my $discordant=$paired-$properly;
my $unmapped=$total - $paired;
my $doubly_unm=$unmapped - $singly_unm;

my $factor=2;
# use 1 for numbers based on readpair mates, 2 for numbers based on readpairs

print "# pe-stats.pl $version\n";
out("total pairs", commafy($total/$factor), perc(1) );
out("  mapped", commafy($paired/$factor), perc($paired/$total) );
out("    concordant", commafy($properly/$factor), perc($properly/$total) );
out("    discordant", commafy($discordant/$factor), perc($discordant/$total) );
out("  unmapped", commafy($unmapped/$factor), perc($unmapped/$total) );
out("    one mate unmapped", commafy($singly_unm/$factor), perc($singly_unm/$total) );
out("    both unmapped", commafy($doubly_unm/$factor), perc($doubly_unm/$total) );
out("mate on diff. chromosome (Q>=5)", commafy($other_chromo/$factor), perc($other_chromo_Q5/$total) );


