# -*- mode: Makefile; -*-

### !! Always run as make --warn-undefined-variables (best redefine make as:
### !! alias make='make --warn-undefined-variables' )

## The following environment variables can be used:
##
## extra_mapping_opts # not used by any rule, but available to e.g. ./do_single_read_mapping
## sort_opts 	      # passed to samtools sort. E.g, '-n' for
## filter_opts 	      # passed to samtools view. e.g. '-f 2' to discard non-paired reads
## chromsizes	      # mandatory, needed by bedtools and bedGraphToBigWig
##
## To use them, do something like
##
##   env variable=value make -f Thismakefile file.bam

## NOTE! make (or: this Makefile) has difficulty with filenames containing '='. Dont use them!

### -------- general stuff --------

makefile_path := $(shell unset CDPATH && cd -P $(dir $(word 1, $(MAKEFILE_LIST))) && pwd -P)
thismakefile := $(makefile_path)/$(word 1, $(MAKEFILE_LIST))

tmpdir := $(TMPDIR)
# (TMPDIR is created and provided by the queueing environment at time of job startup)

## Functions

## check existence of a variable, fail if absent:
check_var = if [ -z "$(strip $($(strip $(1))))" ]; then echo "Missing variable '$(1)', exiting" >&2;  exit 7; fi
## Use as e.g. $(call check_var, chromsizes)

## just warn about missing variable:
warn_var = if [ -z "$(strip $($(strip $(1))))" ]; then echo "Missing variable '$(1)'" >&2;  fi

## check existence of file, fail if absent:
check_exists = if [ ! -s $(1) ]; then echo "File $(1) does not exist or is empty" >&2;  exit 8; fi

## check for created empty target (fail, but don't remove the target)
check_notempty = if [ ! -s $(1) ]; then echo "Empty file '$(1)' created, not removing it, may be latency." >&2; exit 6; fi
## Use exactly (one line, don't forget \t and semi) as
##	COMPLICATED_RECIPE ; $(call check_notempty, $@)

## as check_notempty, but also check exit status:
check_result = if [ $(1) -ne 0 -o ! -s $(2) ]; then echo "Non-zero exit status or empty result for file '$(2)'. This may be latency, not removing it." >&2; exit 6; fi
## Use exactly (one line, don't forget \t and semi) as
## 	COMPLICATED_RECIPE ; $(call check_result, $$?, $@)

## (re)definitions of commands

## nice = /bin/nice
nice = 
movetotrash = move-to-trash.sh		#script that moves, rather than removes
CO = 					#skip version control
## rm = rm -f				# production
rm = $(movetotrash)			# for debugging

bamstat = /usr/bin/true			# bamstat is from bbcfutils. too much work 
                                        # to recompile this, replace with dummy 

## misc. directives
.SUFFIXES: 				# this deletes the traditional .c .h etc. extensions
.SUFFIXES: .fastq .fastq.gz -dusted.fastq .sam .bam .bed .bedGraph .aln .bw
.SECONDARY:				# means: don't remove any intermediate files
.DELETE_ON_ERROR:			# means: any file created by a failing recipe is deleted
.EXPORT_ALL_VARIABLES:
## .ONESHELL:				# pass all lines of a recipe form to one shell invocation


.PHONY: help usage debug smooth extra_mapping_opts sort_opts filter_opts chromsizes 

### -------- mapping-related stuff (formerly in SRA.makefile) -------- 

### General direction is sra -> fastq -> sam -> bam -> bedGraph -> bw

### 

%.fastq: %.sra
	@echo Creating $@ 1>&2
	date
	fastq-dump $*.sra;  $(call check_result, $$?, $@)
	date
	@echo Done creating $@ 1>&2

### removed [bs]am -> fastq, since it goes the wrong way, but leave in 
### for documentation purposes:
## %.fastq: %.bam
## 	@echo Creating $@ 1>&2
## 	bam2fastq --output $*"_R#.fastq" $*.bam # untested
## 	@echo Done creating $@ 1>&2
## 
## ### make sure there's no circular dependency here ...

## %-dusted.fastq: %.fastq
## 	@echo Creating $@ 1>&2
## 	printseq-lite.pl -lc_method dust -lc_method 7  \
##             -fastq $*.fastq -out_format 3 -out_good $*-dusted.fastq -out_bad null;  $(call check_result, $$?, $@)
## 	@echo Done creating $@ 1>&2

%_fastqc.zip: %.bam
	date
	@echo Creating $@ 1>&2
	fastqc --threads 1 --quiet --noextract --format bam $*.bam
##	--thread 1 is useless: limits the number of files, not actual threads!!!
	date
	@echo Done creating $@ 1>&2

%_fastqc.zip: %.fastq.gz
	date
	@echo Creating $@ 1>&2
	fastqc --quiet --noextract --format fastq $*.fastq.gz
	date
	@echo Done creating $@ 1>&2


## single read mapping (this should be done in parallel ...)
## NOTE: code duplication here, the rules are identical apart from prereq suffix
%.sam: %.fastq
	@echo Creating $@ 1>&2
	date
	$(call warn_var, extra_mapping_opts)
	@echo -n 'Script is: '; ls -l `pwd`/do_single_read_mapping
	@echo '============================================================'
	cat do_single_read_mapping
	@echo '============================================================'
	./do_single_read_mapping $< $@; $(call check_result, $$?, $@)
	date
	@echo Done creating $@ 1>&2
## NOTE: code duplication here, same rule as %.sam : %.fastq
%.sam: %.fastq.gz 
	@echo Creating $@ 1>&2
	date
	$(call warn_var, extra_mapping_opts)
	@echo -n 'Script is: '; ls -l `pwd`/do_single_read_mapping
	@echo '============================================================'
	cat do_single_read_mapping
	@echo '============================================================'
	./do_single_read_mapping $< $@; $(call check_result, $$?, $@)
	date
	@echo Done creating $@ 1>&2
## note: extra_mapping_opts is not passed as argument (but through env),
## but it is still useful to check it.

## %.bam: %.fastq.gz extra_mapping_opts
##	@echo  'does not work with large files, fifos choke on the sort step ...'
## 	exit 3
##	date
##	$(call warn_var, extra_mapping_opts)
## 	@echo Creating $@ 1>&2
## 	@echo -n 'Script is: '; ls -l `pwd`/do_single_read_mapping
## 	@echo '============================================================'
## 	cat do_single_read_mapping
## 	@echo '============================================================'
## 	./do_single_read_mapping $*.fastq.gz $*.bam; $(call check_result, $$?, $@)
## 	@echo Done creating $@ 1>&2
## ## note: extra_mapping_opts is not passed as argument (but through env), but it is still
## useful to check it.

## paired-end mapping (this should be done in parallel ...)
%.sam: %_R1.fastq %_R2.fastq 
	@echo Creating $@ 1>&2
	$(call warn_var, extra_mapping_opts)
	@echo -n 'Script is: '; ls -l `pwd`/do_pairedend_read_mapping
	@echo '============================================================'
	cat do_pairedend_read_mapping
	@echo '============================================================'
	./do_pairedend_read_mapping $^ $*.sam; $(call check_result, $$?, $@)
	date
	@echo Done creating $@ 1>&2

### bowtie can also read fastq.gz files, for now just copy the previous rule :-(
%.sam: %_R1.fastq.gz %_R2.fastq.gz
	@echo Creating $@ 1>&2
	$(call warn_var, extra_mapping_opts)
	@echo -n 'Script is: '; ls -l `pwd`/do_pairedend_read_mapping
	@echo '============================================================'
	cat do_pairedend_read_mapping
	@echo '============================================================'
	./do_pairedend_read_mapping $^ $*.sam; $(call check_result, $$?, $@)
	date
	@echo Done creating $@ 1>&2


## simple conversion + sorting + indexing and read-shifting
## (The sorting can also be done in parallel)
%.bam: %.sam 
	@echo Creating $@ 1>&2
	$(call warn_var, filter_opts)
	$(call warn_var, sort_opts)
	samtools view -h $(filter_opts) -b -u $*.sam | samtools sort -@ 4 $(sort_opts) -T$(tmpdir)/sort -Obam > $@
	date
	samtools index $@
	(bamstats.pl $@ ; echo '----'; $(bamstat) $@) > $*.bamstat
	date
	@echo Done creating $@ 1>&2

%Q2.bam: %.bam 
	@echo Creating $@ 1>&2
	samtools view $(filter_opts) -q2 -h -b $< > $@
	samtools index $@
	(bamstats.pl $@ ; echo '----'; $(bamstat) $@) > $*Q2.bamstat
	date
	@echo Done creating $@ 1>&2

%-fwd.bam: %.bam			#assume it's sorted!
	@echo Creating $@ 1>&2
	date
	samtools view -b $(filter_opts) -F20 -h $*.bam -o $@
	date
	samtools index $@
	(bamstats.pl $@ ; echo '----'; $(bamstat) $@) > $*-fwd.bamstat
	date
	@echo Done creating $@ 1>&2

%-rev.bam: %.bam			#assume it's sorted
	@echo Creating $@ 1>&2
	samtools view -b -f18 -F4 -h $*.bam -o $@
	samtools index $@
	(bamstats.pl $@ ; echo '----'; $(bamstat) $@) > $*-rev.bamstat
	@echo Done creating $@ 1>&2


### for completeness:
%.bam.bai: %.bam
	@echo Creating $@ 1>&2
	date
	samtools index $*.bam
	(bamstats.pl $< ; echo '----'; $(bamstat) $< ) > $*.bamstat
	date
	@echo Done creating $@ 1>&2

### stats:

%.stats: %.bam
	samtools stats $*.bam > $*.stats

%.flagstats: %.bam
	samtools flagstat $*.bam > $*.flagstats

%.idxstats: %.bam
	samtools idxstats $*.bam > $*.idxstats

%.pestats: %.flagstats
	pe-stats.pl < $*.flagstats > $*.pestats

%.bamstats: %.bam %.bam.bai %.stats %.idxstats
	@echo Creating $@ 1>&2
	(bamstats.pl $< ; echo '----'; $(bamstat) $< ) > $*.bamstats
	@echo Done creating $@ 1>&2

#### -------- coverage related things (formerly in coverage.makefile) --------

%.bedGraph: %.bam
	@echo Creating $@ 1>&2
	date
	$(call check_var, chromsizes)
	$(call check_exists, $(chromsizes))
	bedtools genomecov -bg -split -ibam $*.bam -g $(chromsizes) -bg | env LC_COLLATE=C sort -k1,1 -k2,2n  > $*.bedGraph
	date
	@echo Done creating $@ 1>&2

%.bw: %.bedGraph
	@echo Creating $@ 1>&2
	date
	$(call check_var, chromsizes)
	$(call check_exists, $(chromsizes))
	bedGraphToBigWig $< $(chromsizes) $@; $(call check_result, $$?, $@)
	date
	@echo Done creating $@ 1>&2

### bed files not needed anymore, so much, out-comment them

## %-fwd.bed: %.bed
## 	@echo Creating $@ 1>&2
## 	grep '+$$' $*.bed > $@
## 	@echo Done creating $@ 1>&2
## 
## %-rev.bed: %.bed
## 	@echo Creating $@ 1>&2
## 	grep '.-$$' $*.bed > $@
## 	@echo Done creating $@ 1>&2

### Note: .bed files are rather big, try to clean up now and again.
### We can get rid of all of the .bed and .bedGraph rules (with fifos), but
### is bit of work
%.bed: %.bam
	@echo Creating $@ 1>&2
	date
	bamToBed -i $< -split > $@ ; ## $(call check_result, $$?, $@)
	date
	@echo Done creating $@ 1>&2

## %.bedGraph: %.bed
## 	@echo Creating $@ 1>&2
## 	date
## 	$(call check_var, chromsizes)
## 	$(call check_exists, $(chromsizes))
## 	bedItemOverlapCount < $< none -chromSize=$(chromsizes) stdin > $@ ; $(call check_result, $$?, $@)
## 	date
## 	@echo Done creating $@ 1>&2

#### -------- Nucleosome data ------------------------------------------------

## Shifted Smoothed Singlebasepair reads
%.insertlen: %.bam
	@echo Creating $@ 1>&2
	samtools view $(filter_opts) $< | sam-insertlengths.pl > $@
	@echo Done creating $@ 1>&2

%-cen.bam: %.bam
	@echo Creating $@ 1>&2
	@echo -n 'Script is: '; ls -l `pwd`/do_center_reads
	@echo '============================================================'
	cat do_center_reads
	@echo '============================================================'
	samtools view $(filter_opts) -h $< | ./do_center_reads > $@; $(call check_result, $$?, $@)
	samtools index $@
	@echo Done creating $@ 1>&2

## Creating stranded bigwigs (shifted and/or 1 bp) 
## The rev-strand bw files have NEGATIVE
## coverages for aesthetic reasons. (that these cannot be log()-ed anymore ! )
## TOO DANGEROUS: if centered reads are again centered, they shift by half-length!!
## Note that IGB and UCSC allow negating values, so there is less need for this ...
### %-fwd.bw: %.bam
### 	@echo Creating $@ ...  1>&2
### 	$(call check_var, chromsizes)
### 	$(call check_exists, $(chromsizes))
### 	$(call check_var, smooth)
### 	date
### 	samtools view  $< -F20 -h |\
### 	  center+smooth.pl --shift 0 --smooth $(smooth) |\
### 	  samtools view -b -u - | bamToBed -split |\
### 	  bedItemOverlapCount none -chromSize=$(chromsizes) stdin >\
###              $*-fwd.bedGraph; $(call check_result, $$?, $*-fwd.bedGraph)
### 	bedGraphToBigWig $*-fwd.bedGraph $(chromsizes) $@; $(call check_result, $$?, $@)
### #	$(rm) $*-fwd.bedGraph
### 	@echo Done creating $@ 1>&2
### 
### %-rev.bw: %.bam
### 	@echo Creating $@ ...  1>&2
### 	$(call check_var, chromsizes)
### 	$(call check_exists, $(chromsizes))
### 	$(call check_var, smooth)
### 	date
### 	samtools view  $< -F4 -f16 -h |\
### 	  center+smooth.pl --shift 0 --smooth $(smooth) |\
### 	  samtools view -b -u -  | bamToBed -split |\
### 	  bedItemOverlapCount none -chromSize=$(chromsizes) stdin |\
### 	  awk -F"\\t" -v OFS="\\t" '{print $$1, $$2, $$3, -$$4}' > \
###              $*-rev.bedGraph; $(call check_result, $$?, $*-rev.bedGraph)
### 	bedGraphToBigWig $*-rev.bedGraph $(chromsizes) $@; $(call check_result, $$?, $@)
### #	$(rm) $*-rev.bedGraph
### 	date
### 	@echo Done creating $@ 1>&2

## -------- Misc. --------

%.aln: %.bed				#for CisGenome
	@echo Creating $@ 1>&2
	file_bed2aln -i $< -o $@ ; $(call check_result, $$?, $@)
	tablesorter_str $@ $*.sort; $(call check_result, $$?, $*.sort)
	mv -f $*.sort $@
	@echo Done creating sorted $@ 1>&2

## ------------------------------------------------------------------------

usage help:
	@echo 
	@echo Usage: make -f $(thismakefile) [ chromsizes=file.txt ] FILE.ext 
	@echo  where ext is one of .fastq, sam, bam, sorted.bam, sorted.bam.bai, .bed
	@echo 

debug: 
	@echo bowtie_genome: $(bowtie_genome) 
	@echo chromsizes: $(chromsizes)
	@echo > /dev/null


