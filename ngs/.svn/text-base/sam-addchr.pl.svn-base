#!/usr/bin/perl

# Add 'chr' in front of seqnames (ensembl genome doesn't use but allows
# them , UCSC genome requires them ...)

# Usage: $0 < file.sam > file-withChr.sam
# Or:  samtools view -h file.bam | $0  | samtools view -h -S -b > file-chr.bam


use strict;


while(<>) { 
    if ( /^\@/ ) { 
        if ( /SN:Mito/ ) { 
            s/SN:Mito/SN:chrmt/;
        } else {
            s/SN:/SN:chr/;
        }
        print;
        next;
    }
    
    my @f=split("\t");
    if($f[2] eq 'Mito') { 
        $f[2]='chrmt'; 
    }  elsif ($f[2] =~ /^[XIV]+$/ ) { 
        $f[2]="chr$f[2]";
    }
    
    print join("\t", @f);
}
