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

warn() { 
    echo "*** $0:  $@ *** " 1>&2 
}
die() { 
    warn "$@"
    exit 28 
}

havebinary() {
    if [ -x "`which $1 2>/dev/null`" ]; then
        return 0
    else
        warn "Could not find $1"
        return 1
    fi
}

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
	fi ;;
*.[1-9n]|*.[1-9]x|*.man)
	if file "$1" | grep -q troff; then
		man -l "$1" | cat -s
		exit $?
	fi ;;
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
		warn "No identify available"
		die "Install ImageMagick or GraphicsMagick to browse images"
	fi ;;

*.json)                                 # pretty-print the json
    if havebinary jq; then
        jq < "$1"
    else
        cat "$1"
    fi ;;

*.gitbundle|*.bundle)                   # list heads and pre-requisite commits
    if havebinary git; then
        if file "$1" | grep -q 'Git bundle'; then
            git bundle verify "$1"
        else
            warn 'I thought this was a git bundle ... '
            cat "$1"
        fi
    else
        cat "$1"
    fi;;

### bioinformatics:
*.ubam|*.bam|*.cram) 
        ## Use Samtools to view a next generation sequencing genome alignment file (.[u]bam or .cram)
        if havebinary samtools; then
            (echo '-- First 20 header lines: --'
             samtools view -H "$1" | head -20
             echo '-- First 10 readgroups (if any): --'
             samtools view -H "$1" | grep '^@RG' | head -10
             echo '-- Program headers: --'
             samtools view -H "$1" | egrep '^@(PG|CO)'
             echo '-- Skipping to first read: --'
             samtools view "$1")
        else
            warn "Need samtools for this"
	    cat "$1"
        fi ;;

*.bw|*.bigwig)
        if havebinary bigWigInfo; then
            bigWigInfo "$1"
        else 
	    die "need bigWigInfo for this (part of the UCSC suite)"
            cat "$1"
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
