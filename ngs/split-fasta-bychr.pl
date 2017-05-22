#!/usr/bin/perl
# From one big fasta file containing many sequences, extract separate *fa.gz files, each containing one
# sequence. The new files are named by the first word on the '>'-line, overwriting anything that may
# already be present.
# 
# Usage:   zcat bigfasta-with-many-seqs.fa.gz | split-fasta-bychr.pl
#
use strict;
use FileHandle;

my $fh= undef;
my $n=0;

while (<>) {
  if ( /^>\s*(\S+)/  ) { 
    my $name="$1.fa.gz";
    $fh->close() if defined($fh);
    warn "(over)writing $name ...\n";
    $fh = FileHandle->new("| gzip -n > $name");
    $n++;
  }
  print $fh $_;
}
$fh->close();
warn "Created $n *.fa.gz files\n";
