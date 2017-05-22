#!/usr/bin/perl

# Usage: split-bychr.pl *.bed
#
# This script 'sorts' each chromosome into the right file as it is read, which takes just
# one pass over the data. A simple-minded grep would be much slower (unless you want to 
# select just one chromosome).

use strict;
use warnings;
no warnings 'uninitialized'; 

use Carp;

use File::Basename;
use FileHandle;

## our @chroms=qw(chrI    chrII   chrIII  chrIV   chrV    chrVI   chrVII 
##                chrVIII chrIX   chrX    chrXI   chrXII  chrXIII chrXIV 
##                chrXV   chrXVI);##  chrmt  
## 

confess "Usage: $0 *.{bed,bedGraph}" if(@ARGV == 0); 

my @known_exts=qw(.bed);##  .bedGraph);

my $dir="byChr";
opendir(DIR, "$dir") or confess "Need directory $dir to write per-chromosome files to";
close(DIR);


for my $file (@ARGV) { 
    my ($name, $path, $ext) = fileparse($file, @known_exts);
    confess "Unknown suffix $ext. Known (so far) are ". join(' ', @known_exts)
      unless $ext;
    my $newdir="$path$dir";
    opendir(DIR, $newdir) or confess "No directory $newdir to write per-chromosome files to,";
    close(DIR);
    
    if ($ext eq '.bed' ) { 
        split_bed($newdir, $file, $name, $ext);
    } else {
      confess "Oops, not implemented";
  }
}

sub split_bed { 
    my($newdir, $file, $name, $ext)=@_;
    my $fileHandles={};
    my $fileNames={};
    my $lines={};

    open(READ, $file) or confess "$file: $!,";
    while(<READ>) {
        ### my($chr, @rest)=split("\t", $_);
        my ($chr)=substr($_, 0, index($_, "\t")); # should be faster
        my ($fh)=$fileHandles->{$chr};
        if(!$fh) { 
            my $newfile="$newdir/$chr,${name}$ext";
            $fh = FileHandle->new("> $newfile") or confess "$newfile: $!";
            $fileHandles->{$chr}=$fh;
            $fileNames->{$chr}=$newfile; 
        }
        $lines->{$chr}++;
        print $fh $_;
    }

    close(READ) or warn "Something unusual when closing file $file: $!";
    foreach my $chr (sort keys %$fileHandles) { 
        my $name=$fileNames->{$chr};
        close($fileHandles->{$chr})
          or warn "Something usual when closing file $name:  $!";
        warn("Created $name (". $lines->{$chr}. " lines)\n");
    }
}

