#!/usr/bin/perl

# Add 'chr' in front of seqnames (ensembl genome doesn't use but allows
# them , UCSC genome requires them ...)

# Usage: $0 < file.bedGraph > file-withChr.bedGraph
#  (also works for wig flies)


use strict;


while(<>) { 
    if ( /^#/ ) {
         print; next;
    }
    if (/^[XIV]+\t/ ) { 
        print "chr$_";
    } elsif ( /^Mito/ ) { 
        s/^Mito/chrmt/;
        print;
    } else { 
        die "don't recognize this line: '$_'\n (line $.)";
    }
}
