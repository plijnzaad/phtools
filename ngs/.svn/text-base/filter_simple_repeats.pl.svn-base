#!/usr/bin/env perl

############     filter_simple_repeats.pl     ############
#
#	This Perl script reads sequence files in FASTA format
#	and filters out homo- and/or dipolymeric stretches. This
#	can be useful prior to mapping reads to reference genomes.
#
#	The script is part of the "NGS tools for the novice"
#	authored by David Rosenkranz, Institue of Anthropology
#	Johannes Gutenberg University Mainz, Germany
#
#	Contact: rosenkrd@uni-mainz.de

############     HOW TO USE     ############
#
#	You can pass file names of the files to be processed and
#	options for the filtering process via arguments to the script.
#
#	Input files have to be passed to the scripts via -i:
#	-i file1.fas -i file2.fas
#
#	A bare digit from 1-3 (1 or 2 or 3) will tell the script which
#	reads to filter out. By default it is set to 3.
#	-1 = sequence reads with homopolymeric stretches.
#	-2 = sequence reads with dipolymeric stretches.
#	-3 = sequence reads with homo- or dipolymeric stretches.
#
#	You can set a threshold for the filtering process that
#	determines the minimum size or relative amount of a
#	homo- or dipolymeric stretch via -pt [%] (for a percentual
#	thresholt) or -at [nt] (for an absolute thresholt). I you
#	use both arguments, only one parameter has to comply
#	with the requirements to cause the rejection of the read:
#	-rt 75 -at 20
#
#	By default the script runs with -rt 100 and -at 999999 which
#	in practice means that the whole read has to be a simple repeat to
#	be filtered out.
#
#	You can choose to allow one mismatch including insertion/
#	deletion within the homo- or dipolymeric stretch when searching
#	for simple repeats via -mm. This will make the filtering more
#	strict. The default is 0 (= no mismatch allowed):
#	-m 1
#
#	The script will create the following selfexplanatory output
#	files:
#	- simple_repeats.fas
#	- no_simple_repeats.fas
#
#	For example you can type the following command:
#	perl filter_simple_repeats.pl -i input_file.fas -rt 75 -at 20
#
#	If you do not want to enter each file name seperately, you
#	can provide a file that contains a list all file names (one
#	file name per line). Pass the file name to the script via -I:
#	-I list_of_files.txt
#
#	Multiple files and combinations of all kinds of arguments are
#	allowed:
#	perl filter_simple_repeats.pl -i input_file.fas -I list_of_files.txt -rt 75 -at 20



@input_files=();
$how_to_filter=3;
$mismatch=0;
$relative_threshold=100;
$absolute_threshold=999999;
$|=1;

###   CHECK COMMAND LINE ARGUMENTS   ###
if(@ARGV==0)
	{
	print"No arguments passed to the script!\nIf you entered arguments try the following command:\nperl filter_simple_repeats.pl -argument1 -argument2 ...\n\n";
	exit;
	}

$argv="";
foreach(@ARGV)
	{
	$argv.=$_;
	}
@arguments=split('-',$argv);

foreach(@arguments)
	{
	if($_=~/^ *i/)
		{
		$_=~s/^ *i//;
		$_=~s/ //g;
		push(@input_files,$_);
		}
	elsif($_=~/^ *I/)
		{
		$_=~s/^ *I//;
		$_=~s/ //g;
		open(FILES_IN,"$_");
		while(<FILES_IN>)
			{
			unless($_=~/^\s*$/)
				{
				$_=~s/\s//sg;
				push(@input_files,$_);
				}
			}
		}
	elsif($_=~/^ *[123] *$/)
		{
		$_=~s/ *//g;
		$_=~s/ //g;
		$how_to_filter=$_;
		}
	elsif($_=~/^ *m *[01]/)
		{
		$_=~s/^ *m//;
		$_=~s/ //g;
		$mismatch=$_;
		}
	elsif($_=~/^ *at/)
		{
		$_=~s/^ *at//;
		$_=~s/ //g;
		$absolute_threshold=$_;
		}	
	elsif($_=~/^ *rt/)
		{
		$_=~s/^ *rt//;
		$_=~s/ //g;
		$relative_threshold=$_;
		}
	elsif($_!~/^\s*$/)
		{
		print"Don't know how to treat argument $_!\nIt will be ignored.\n\n";
		}
	}
if(@input_files==0)
	{
	print"No input file specified!\n";
	exit;
	}
unless($absolute_threshold=~/^\d+$/)
	{
	print"Absolute threshold has to be numerical!\n";
	exit;
	}
unless($relative_threshold=~/^\d+$/)
	{
	print"Relative threshold has to be numerical!\n";
	exit;
	}

###   PRINT ARGUMENTS   ###
print"The following files will be processed:\n";
foreach(@input_files)
	{
	if(-e $_)
		{
		print"$_\n";
		push(@input_files_ok,$_);
		}
	else
		{
		print"could not find file: $_. It will be ignored.\n";
		}
	}
if($how_to_filter==1)
	{
	print"\nFiltering homopolymeric stretches.\n";
	}
elsif($how_to_filter==2)
	{
	print"\nFiltering dipolymeric stretches.\n";
	}
elsif($how_to_filter==3)
	{
	print"\nFiltering homo- and dipolymeric stretches.\n";
	}
print"Relative threshold [%]:  $relative_threshold\nAbsolute threshold [nt]: $absolute_threshold\n\n";

###   START   ###
@homopolymers_0mm=('A+','T+','G+','C+');
@homopolymers_1mm=('A+.?A*','T+.?T*','G+.?G*','C+.?C*');
@dipolymers_0mm=('T?(AT)+A?','G?(AG)+A?','C?(AC)+A?','G?(TG)+T?','C?(TC)+T?','C?(GC)+G?');
@dipolymers_1mm=('T?(AT)+A?.?T?(AT)+A?','G?(AG)+A?.?G?(AG)+A?','C?(AC)+A?.?C?(AC)+A?','G?(TG)+T?.?G?(TG)+T?','C?(TC)+T?.?C?(TC)+T?','C?(GC)+G?.?C?(GC)+G?');
if($mismatch==0)
	{
	@homopolymers=@homopolymers_0mm;
	@dipolymers=@dipolymers_0mm;
	}
elsif($mismatch==1)
	{
	@homopolymers=@homopolymers_1mm;
	@dipolymers=@dipolymers_1mm;
	}

open(OUT_SIMPLE,">simple_repeats.fas");
open(OUT_OK,">no_simple_repeats.fas");
$ok=0;
$rejected=0;
foreach$file(@input_files_ok)
	{
	print"processing $file";
	open(IN,$file);
	while(<IN>)
		{
		if($_=~/^>/)
			{
			$title=$_;
			}
		else
			{
			$seq=$_;
			chomp$seq;
			$seq_length=length$seq;
			$seq_for_print=$seq;
			$reject_sequence=0;
			if($how_to_filter==1||$how_to_filter==3)
				{
				foreach$homopolymer(@homopolymers)
					{
					last if $reject_sequence==1;
					while($seq=~s/$homopolymer//i)
						{
						last if $reject_sequence==1;
						if(length$&>=$absolute_threshold)
							{
							$reject_sequence=1;
							}
						elsif(((length$&)*100)/$seq_length>=$relative_threshold)
							{
							$reject_sequence=1;
							$test1=length$&;
							$test2=length$seq;							
							}
						}
					}
				}
			$seq=$seq_for_print;
			if($how_to_filter==2||$how_to_filter==3)
				{
				foreach$dipolymer(@dipolymers)
					{
					last if $reject_sequence==1;
					while($seq=~s/$dipolymer//i)
						{
						last if $reject_sequence==1;
						if(length$&>=$absolute_threshold)
							{
							$reject_sequence=1;
							}
						elsif(((length$&)*100)/$seq_length>=$relative_threshold)
							{
							$reject_sequence=1;
							}
						}
					}
				}
			if($reject_sequence==0)
				{
				$ok++;
				print OUT_OK "$title$seq_for_print\n";
				}
			elsif($reject_sequence==1)
				{
				$rejected++;
				print OUT_SIMPLE "$title$seq_for_print\n";
				}
			}
		}
	close IN;
	print" done.\n";
	}
close OUT_OK;
close OUT_SIMPLE;
print"Passed sequences:\t$ok\nRejected sequences:\t$rejected\n\n";
exit;
