#!/usr/bin/perl 

### silly script to quickly produce a GFF file. Currently only knows about
### perfect single-segment matches

use strict;
use warnings;
no warnings qw(uninitialized);

use Carp;

use Getopt::Long;

my $source = 'unknown';
my $type = 'exon';
my $help=0;

my $usage = "$0 [-h] [ -s  source ] [ -t featuretype ] < input.sam  > output.gff\n";

die $usage if ( !GetOptions('h' => \$help,
                            's=s' => \$source,
                            't=s' => \$type)
                or $help
                or @ARGV);



my $flags= { 0x1 => 'multiple segments',
             0x2 => 'properly aligned',
             0x4 => 'unmapped',
             0x8 => 'mate unmapped',
             0x10 => 'revcomped',
             0x20 => 'mate revcomped',
             0x40 => 'first segment',
             0x80 => 'last segment',
             0x100 => 'secondary algnmt',
             0x200 => 'bad quality',
             0x400 => 'PCR or optical duplicate'
            };

for my $key ( keys %$flags) { $flags->{ $flags->{$key} } = $key; } # also provide reverse

print "##gff-version 3\n";

my $seen={};

LINE:
while(<>) { 
    chomp;
    my ($qname, $flag, $refname, $pos, 
        $qual, $cigar, $matename, $matepos, 
        $length, $seq, $phred)=split("\t");
    next LINE if $flag & $flags->{unmapped};

    my $score= '.';
    my $end=$pos + length($seq);
    my $strand =  ($flag & $flags->{revcomped})  ? '-' : '+';
    my $frame = '.';
    my $id=$qname;
    if ($seen->{$qname} ++) { 
##         $id = $qname . chr( ord('a')-1 + $seen->{$qname} );
        $id = "$qname-$seen->{$qname}";
    }
    print join("\t", $refname, $source, $type, $pos, $end, $score, $strand, $frame, "ID=$id"), "\n";
        
}
