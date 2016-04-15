#!/bin/sh
#
# This lesspipe.sh is adapted from the /usr/bin/lesspipe.sh shipping
# with CentOS7. It is extended to be able to read .bam and .cram files
# (genome alignment files used in next generation sequencing) This
# requires the installation of samtools; if not available, less will
# revert to the old behaviour (saying '%s may be a binary file.  See it
# anyway?')
#
# To use this filter with less, define LESSOPEN as
#
# export LESSOPEN="| INSTALLATIONDIRECTORY/lesspipe.sh %s"
#
# (where INSTALLATIONDIRECTORY is the location of the file you're now reading)
#
# (If you use this as part of a lmod module, do something like
#
#     local base = "/usr/local/src/git/phtools/ngs"
#     prepend_path("PATH", base)
#     load("samtools")  # makes sure that samtools view is on the path
#     setenv("LESSOPEN", "| " .. pathJoin(base, "lesspipe.sh") .. " %s"))
#
# Extension by plijnzaad@gmail.com; see https://github.com/plijnzaad/phtools/tree/master/ngs


# Note: the script should return zero if the output was valid and non-zero
# otherwise, so less could detect even a valid empty output
# (for example while uncompressing gzipped empty file).
# For backward-compatibility, this is not required by default. To turn
# this functionality there should be another vertical bar (|) straight
# after the first one in the LESSOPEN environment variable:
# export LESSOPEN="||/usr/bin/lesspipe.sh %s"

if [ ! -e "$1" ] ; then
	exit 1
fi

if [ -d "$1" ] ; then
	ls -alF -- "$1"
	exit $?
fi

# exec 2>/dev/null

# Allow for user defined filters
if [ -x ~/.lessfilter ]; then
	~/.lessfilter "$1"
	if [ $? -eq 0 ]; then
		exit 0
	fi
fi

case "$1" in
*.[1-9n].bz2|*.[1-9]x.bz2|*.man.bz2|*.[1-9n].[gx]z|*.[1-9]x.[gx]z|*.man.[gx]z|*.[1-9n].lzma|*.[1-9]x.lzma|*.man.lzma)
	case "$1" in
	*.gz)		DECOMPRESSOR="gzip -dc" ;;
	*.bz2)		DECOMPRESSOR="bzip2 -dc" ;;
	*.xz|*.lzma)	DECOMPRESSOR="xz -dc" ;;
	esac
	if [ -n "$DECOMPRESSOR" ] && $DECOMPRESSOR -- "$1" | file - | grep -q troff; then
		$DECOMPRESSOR -- "$1" | groff -Tascii -mandoc -
		exit $?
	fi ;;&
*.[1-9n]|*.[1-9]x|*.man)
	if file "$1" | grep -q troff; then
		man -l "$1" | cat -s
		exit $?
	fi ;;&
*.tar) tar tvvf "$1" ;;
*.tgz|*.tar.gz|*.tar.[zZ]) tar tzvvf "$1" ;;
*.tar.xz) tar Jtvvf "$1" ;;
*.xz|*.lzma) xz -dc -- "$1" ;;
*.tar.bz2|*.tbz2) bzip2 -dc -- "$1" | tar tvvf - ;;
*.[zZ]|*.gz) gzip -dc -- "$1" ;;
*.bz2) bzip2 -dc -- "$1" ;;
*.zip|*.jar|*.nbm) zipinfo -- "$1" ;;
*.rpm) rpm -qpivl --changelog -- "$1" ;;
*.cpi|*.cpio) cpio -itv < "$1" ;;
*.gpg) gpg -d "$1" ;;
*.gif|*.jpeg|*.jpg|*.pcd|*.png|*.tga|*.tiff|*.tif)
	if [ -x /usr/bin/identify ]; then
		identify "$1"
	elif [ -x /usr/bin/gm ]; then
		gm identify "$1"
	else
		echo "No identify available"
		echo "Install ImageMagick or GraphicsMagick to browse images"
		exit 1
	fi ;;
### bioinformatics:
*.bam|*.cram) 
        ## Use Samtools to view a next generation sequencing genome alignment file (.bam or .cram)
        if [ -x "`which samtools 2>/dev/null`" ]; then
            samtools view -h "$1"
        else 
	    echo -e "$0: cannot find samtools" >&2
            exit 1
        fi ;;
*.bw|*.bigwig)
        ## Use bigWigInfo
        if [ -x "`which bigWigInfo 2>/dev/null`" ]; then
            bigWigInfo "$1"
        else 
	    echo -e "$0: cannot find bigWigInfo (part of the UCSC suite) " >&2
            exit 1
        fi ;;

### last resort:
*)
	if [ -x /usr/bin/file ] && [ -x /usr/bin/iconv ] && [ -x /usr/bin/cut ]; then
		case `file -b "$1"` in
		*UTF-16*) conv='UTF-16' ;;
		*UTF-32*) conv='UTF-32' ;;
		esac
		if [ -n "$conv" ]; then
			env=`echo $LANG | cut -d. -f2`
			if [ -n "$env" -a "$conv" != "$env" ]; then
				iconv -f $conv -t $env "$1"
				exit $?
			fi
		fi
	fi
	exit 1
esac
exit $?

