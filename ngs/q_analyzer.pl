#!/usr/bin/perl

### taken from http://sourceforge.net/projects/ngs-toolbox/ (as per 4sep2013)

############     q_analyzer.pl     ############
#
#	This Perl script reads sequence files in FASTQ format
#	and outputs some statistical quality related information
#	based on Phred quality scores.
#
#	The script is part of the "NGS tools for the novice"
#	authored by David Rosenkranz, Institue of Anthropology
#	Johannes Gutenberg University Mainz, Germany
#
#	Contact: rosenkrd@uni-mainz.de


############     HOW TO USE     ############
#
#	You can pass file names of the files to be processed, the
#	output file name and options via arguments to the script.
#
#	Input files have to be passed to the scripts via -i:
#	-i file1.fas -i file2.fas
#
#	If you do not want to enter each file name seperately, you
#	can provide a file that contains a list all file names (one
#	file name per line). Pass the file name to the script via -I:
#	-I list_of_files.txt
#
#	By default, the script assumes your FASTQ file to be in
#	Sanger format (scores from 0-93 using ASCII 33-126).
#	You can change this setting via -f:
#	-f Sanger	| if your FASTQ file is in Sanger format.
#	-f Illumina	| if your FASTQ file is in Illumina format.
#
#	Selecting the correct format is crucial for correct computation!
#	Info:
#	Illumina 1.0+ = Illumina format
#	Illumina 1.8+ = Sanger format
#
#	By default, the results will be displayed in the command
#	promt / terminal. Alternatively, you can pass an output file
#	name to the script via -o:
#	-o output_file.txt
#
#	Multiple files and combinations of all kinds of arguments
#	are allowed:
#	perl q_analyzer.pl -i input_file.fas -I list_of_files.txt -f Illumina -o output_file.txt




@input_files=();
$format="Sanger";
@Phreds=();
$|=1;

###   CHECK COMMAND LINE ARGUMENTS   ###
if(@ARGV==0)
	{
	print"No arguments passed to the script!\nIf you entered arguments try the following command:\nperl q_analyzer.pl -argument1 -argument2 ...\n\n";
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
	elsif($_=~/^ *f/)
		{
		$_=~s/^ *f//;
		$_=~s/ //g;
		$format=$_;
		}
	elsif($_=~/^ *o/)
		{
		$_=~s/^ *o//;
		$_=~s/ //g;
		$output=$_;
		}
	elsif($_!~/^\s*$/)
		{
		print"Don't know how to treat argument $_!\nIt will be ignored.\n\n";
		}
	}
if($format!~/^ *Illumina *$/i&&$format!~/^ *Sanger *$/i)
	{
	print"Input format has to be 'Sanger' or 'Illumina'\n";
	exit;
	}
if(@input_files==0)
	{
	print"No input file specified!\n";
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
if($format=~/Illumina/)
	{
	foreach(0..62)
		{
		$Phreds[$_]=0;
		}
	print"\nInput files are assumed to be Illumina format.\n";
	}
elsif($format=~/Sanger/)
	{
	foreach(0..93)
		{
		$Phreds[$_]=0;
		}
	print"\nInput files are assumed to be Sanger format.\n";
	}

###   START   ###
$bases_total=0;
$terminal_Bs=0;
$seqs_with_terminal_Bs=0;
$terminal_rhombs=0;
$seqs_with_terminal_rhombs=0;
@positional_Phred=();
@n_bases_at_position=();
$total_reads=0;
$average_P0errors=0;
foreach$file(@input_files_ok)
	{
	print"processing $file";
	open(IN,$file);
	while(<IN>)
		{
		if($_=~/^@/&&$title_index==0)
			{
			$title=$_;
			$title=~s/^@/>/;
			$title_index=1;
			}
		elsif($_=~/^[ATGCN]+\n?$/&&$seq_index==0)
			{
			$sequence=$_;
			$seq_index=1;
			}
		elsif($_=~/^\+/&&$quality_prefix_index==0)
			{
			$quality_header=$_;
			$quality_prefix_index=1;
			}
		elsif($quality_index==0)
			{
			$quality=$_;
			$quality_index=1;
			}
		if($title_index+$seq_index+$quality_prefix_index+$quality_index==4)
			{
			$total_reads++;
			$title_index=0;
			$seq_index=0;
			$quality_prefix_index=0;
			$quality_index=0;
			chomp$quality;
			if($quality=~/B+$/)
				{
				$terminal_Bs+=length$&;
				$seqs_with_terminal_Bs++;
				}
			if($quality=~/#+$/)
				{
				$terminal_rhombs+=length$&;
				$seqs_with_terminal_rhombs++;
				}
			$pos=0;
			$prob_0_errors=1;
			@quality=split('',$quality);
			if($format=~/Illumina/)
				{
				foreach(@quality)
					{
					$pos++;
					$Phred=ord($_)-64;
					$Phreds[$Phred]++;
					$prob=1-(10**(($Phred*-1)/10));
					$prob_0_errors=$prob_0_errors*$prob;
					$positional_Phred[$pos]+=$Phred;
					$n_bases_at_position[$pos]++;
					$average_score+=$Phred;
					$bases_total++;
					}
				}
			else
				{
				foreach(@quality)
					{
					$pos++;
					$Phred=ord($_)-33;
					$Phreds[$Phred]++;
					$prob=1-(10**(($Phred*-1)/10));
					$prob_0_errors=$prob_0_errors*$prob;
					$positional_Phred[$pos]+=$Phred;
					$n_bases_at_position[$pos]++;
					$average_score+=$Phred;
					$bases_total++;
					}
				}
			$average_P0errors+=$prob_0_errors;
			}
		}
	close IN;
	print" done.\n";
	}
$average_P0errors=$average_P0errors/$total_reads;
$average_score=$average_score/$bases_total;
$longest_read=@positional_Phred;
foreach(1..$longest_read)
	{
	if($n_bases_at_position[$_]>0)
		{
		$positional_Phred[$_]=$positional_Phred[$_]/$n_bases_at_position[$_];
		}
	}

$pop=1;
while($pop==1)
	{
	if($Phreds[-1]==0)
		{
		pop@Phreds;
		}
	else
		{
		$pop=0;
		}
	}
unless($output)
	{
	print"\n\nRESULTS:\nAverage probability for overall accuracy of a read: $average_P0errors";
	print"\nAverage Phred score for base calling: $average_score\n";
	print"\nPhred score\tNumber of base calls\n";
	$element=-1;
	foreach(@Phreds)
		{
		$element++;
		if($element<10)
			{
			$element=" ".$element;
			}
		print"$element\t\t$_\n";
		}
	print"\n";
	print"Sequences ending with 'B'-scored bases*: $seqs_with_terminal_Bs\n";
	print"Sequences ending with '#'-scored bases*: $seqs_with_terminal_rhombs\n";
	print"Total terminal 'B'-scored bases*: $terminal_Bs\n";
	print"Total terminal '#'-scored bases*: $terminal_rhombs\n\n";
	print"*\nIn Illumina 1.5+ stretches of B corresponding to a Phred score of 2\n";
	print"are used to indicate, that a specific final proportion of the read\n";
	print"should not be used in further analysis. In Illumina 1.8+ (Sanger format)\n";
	print"low quality ends are indicated by strectches of #\n\nSeq. position\tAverage score for base calling:\n";
	$element=0;
	foreach(@positional_Phred)
		{
		$element++;
		print"$element\t\t$_\n";
		}
	print"\n";
	}
else
	{
	open(OUT,">$output");
	print OUT"RESULTS:\nAverage probability for overall accuracy of a read: $average_P0errors";
	print OUT"\nAverage Phred score for base calling: $average_score\n";
	print OUT"\nPhred score\tNumber of base calls\n";
	$element=-1;
	foreach(@Phreds)
		{
		$element++;
		if($element<10)
			{
			$element=" ".$element;
			}
		print OUT"$element\t\t$_\n";
		}
	print OUT"\n";
	print OUT"Sequences ending with 'B'-scored bases*: $seqs_with_terminal_Bs\n";
	print OUT"Sequences ending with '#'-scored bases*: $seqs_with_terminal_rhombs\n";
	print OUT"Total terminal 'B'-scored bases*: $terminal_Bs\n";
	print OUT"Total terminal '#'-scored bases*: $terminal_rhombs\n\n";
	print OUT"*In Illumina 1.5+ stretches of B corresponding to a Phred score of 2\n";
	print OUT" are used to indicate, that a specific final proportion of the read\n";
	print OUT" should not be used in further analysis. In Illumina 1.8+ (Sanger format)\n";
	print OUT" low quality ends are indicated by strectches of #\n\nSeq. position\tAverage score for base calling:\n";
	$element=0;
	foreach(@positional_Phred)
		{
		$element++;
		print OUT"$element\t\t$_\n";
		}
	print OUT"\n";
	close OUT;
	}
exit;
