#!/bin/sh
# -*- mode: sh; -*-

queue=veryshort                        # default queue
logdirname=logs                        # always relative to current dir!

usage_msg="
\n
\n  Usage: $(basename $0) [qsub-options] command arg1 arg2 ...
\n
\n  qrun.sh is a qsub replacement for simple job submissions.  It frees
\n  you from having to write a wrapper scripts around the simple command
\n  line that you use interactively: simply put 'qrun.sh' in front of it,
\n  and it will run on the cluster, rather than locally. The qrun.sh
\n  script takes care of all the directories, environment variables
\n  (including the various PATHs and LIBs). Stdout and stderr go to subdir
\n  ./$logdirname, which is created if not present in the current directory. Default
\n  queue is $queue, but you can override this with the -q QUEUENAME
\n  option. A job name, if not specified using -N JOBNAME, is 'invented'
\n  based on the name of submitted command + arguments.
\n
\n  For simplicity, there is only a limited number of qsub options
\n  available, most importantly -q (to specify which queue) and -N (to
\n  specify the jobname). Other qsub options currently recognized are
\n  -j, -M -m, -p and -l. More may be added in the future (feel free to
\n  request them from plijnzaad@gmail.com). For the -p option
\n  (which corresponds to qsub's -pe) you can specify the environment
\n  and the number separated by '='. I.e. both -p 'threaded 4' and 
\n  -p threaded=4 are fine (tha latter may save you quoting trouble).
\n
\n  Unix pipe lines can be specified as a job, but the pipeline must be
\n  quoted as a whole, which would look like
\n
\n  $ qrun.sh \"cmd1 | cmd2 > out.txt\"
\n
\n  Be careful in this case with the dollar signs of variable references (including
\n  things like \$1, \$2 in a one-liner awk or perl-script). Escape them to
\n  avoid premature substitution by the invoking shell.
\n
\n  Due to some brain damage on the part of both CentOS and SGE, you may
\n  see errors like 'sh: which: line 1: syntax error: unexpected end of
\n  file', which are harmless.
\n"

##
## TODO: make sure any whitepace inside directory and filenames remains
##   intact; 
##
## Written by Philip Lijnzaad <plijnzaad@gmail.com>, November 2014

## Avoid SGE/OGE, and/or tcsh and/or Redhat's Modules brain damage:
unset module which BASH_FUNC_which >/dev/null 2>&1
unalias module which BASH_FUNC_which >/dev/null 2>&1
## Without this, you get the following in all stderr of all scripts submitted
## with the current wrapper script:
##
##   /bin/sh: module: line 1: syntax error: unexpected end of file
##   /bin/sh: error importing function definition for `BASH_FUNC_module'
##

error(){ 
    echo -e "$*" >&2; exit 6 
}

usage(){ 
    error $usage_msg
}

if [ $# -eq 0 ]; then usage; fi

# default mail adres:
user=$(id -un)
mail_addres=$(ldapsearch -x "(& (objectClass=person)(uid=$user))" mail |  awk '/^mail/{print $2}')
opt_M="-M $mail_addres"

# default: sent mail when jobs is aborted
opt_m="-m a"

jobname=""
while getopts "N:q:o:e:j:M:m:p:l:" opt; do
    case $opt in
        N)
            jobname=$OPTARG
            ;;
        q)
            queue=$OPTARG
            ;;
        o|e)
            error "You cannot use the -o and -e options; output always goes to ./$logdirname/ (relative to current diretory)"
            ;;
        j)
            opt_j="-j $OPTARG"
            ;;
        m)
            opt_m="-m $OPTARG"
            ;;
        l)
            opt_m="-l $OPTARG"
            ;;
        M)
            opt_M="-M $OPTARG"
            ;;
        p)
            opt_pe="-pe $OPTARG"        # sorry, -p has to be one-letter options
            opt_pe=$(echo $opt_pe | tr '=' ' ')
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "$jobname" ]; then
    jobname=$(echo $@)
fi
jobname=$(echo $jobname | sed 's/[ /<>|][ /<>|]*/_/g;s/^[._]*//;')

### find/set up the directory for stdout and stderr logs:
here="$(pwd -P)"
if [ $(basename "$here")  == $logdirname ]; then
    error "You're in a $logdirname/ directory ( $here ), you prolly meant to invoke it from one level up. Job not submitted"
fi

logdir="$here/$logdirname"
if [  -d "$logdir" ]; then
:
else
    if mkdir "$logdir" ; then
        echo "Created directory $logdir for stdout and stderr" >&2
    else
        eror "Could not create directory $logdir for stdout and stderr. Job not submitted" >&2
    fi
fi
echo "Using directory $logdir for stdout and stderr" >&2

other_opts="$opt_j $opt_m $opt_M $opt_pe"
qsub_opts="-shell no -b yes -cwd -V -o $logdir -e $logdir -q $queue -N $jobname $other_opts"

## finally, submit:
cmd=$(echo "$@")                        # remove one level of quoting ...
qsub $qsub_opts sh -c "$cmd"
