#!/usr/bin/env perl

my $stats={};
my $txpt_lens={};

LINE:
while(<>) { 
  s/[\n\r]*$//;

  if (/^\@/) { 
    next LINE unless /^\@SQ/;
    
    my($tag, $name, $len)=split(' ');
    $name =~ s/^SN://;
    $len =~ s/^LN://;
    $txpt_lens->{$name}=$len;
  }

  my($qname,$flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen,
     $seq, $qual, @rest)=split("\t", $_);
  next LINE if $rname eq '*';

  ## from scseq/process_sam_cel384v2.pl:
  my $X0 = 0;
  my $dum = 'NA';
  foreach my $el (@rest){
    # ($dum,$dum,$NM) = split(":",$el) if ($el =~ /^NM\:/); # NM: number of mismatches
    # ($dum,$dum,$XA) = split(":",$el) if ($el =~ /^XA\:/); # XA: number of alternative hits (chr,pos,CIGAR,NM;)+
    ($dum,$dum,$X0) = split(":",$el) if ($el =~ /^X0\:/); # X0: number of best hits (bwa-specific!)
  }
  next LINE if ($X0 != 1) || ($flag & 16); # get rid of multimappers and antisense mappers (for now)
  my $len=$txpt_lens->{$rname};
  print "$pos\t$len\t$rname\n";
}                                       # while
