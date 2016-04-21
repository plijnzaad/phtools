#!/bin/sh
# -*- mode: sh; -*-

if lsb_release -r | grep -q '[	 ]7'; then 
  queue=all.q
else
  queue=veryshort                        # default queue
fi

logdirname=logs                        # always relative to current dir!

usage_msg="
\n
\n  Usage: $(basename $0) [ qsub-options ] command arg1 arg2 ...
\n
\n  qrun.sh is a qsub replacement for simple job submissions.  It frees
\n  you from having to write wrapper scripts around the simple command
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
\n  available, all of them single-letter options (SGE's qsub are not
\n  always!).  The most important ones are: -q (specifies queue) and -N
\n  (specifies jobname). Other qsub options currently recognized are -j
\n  (join stdout and stderr), -M (mail address to use) -m (when to mail),
\n  -p (priority) and -l (resources). For the -p option (corresponds
\n  to qsub's \"-pe\") you can specify the environment and the number
\n  separated by '='. I.e. both -p 'threaded 4' and -p threaded=4 are
\n  fine (the latter may save you quoting trouble). Option -h
\n  correpsonds to qsub's option \"-hold_jid\"; option -H to qsub's
\n  \"-hold_jid_ad\"
\n
\n  For the -l option, you have to specify all resources in one string,
\n  separated by commas. E.g.: -l \"h_rt=2:30:00,h_vmem=32G,tmpspace=100G\".
\n
\n  Unix pipe lines can be specified as a job, but the pipeline must be
\n  quoted as a whole, which would look like
\n
\n  $ qrun.sh \"cmd1 | cmd2 > out.txt\"
\n
\n  Be careful in this case with the dollar signs of variable references
\n  (including things like \$1, \$2 in a one-liner awk or
\n  perl-script). Escape them to avoid premature substitution by the
\n  invoking shell.
\n
\n  If your system has no ldapsearch for looking up the e-mail address, 
\n  consider creating a file $HOME/.forward that contains the right
\n  e-mail address (file must be unreadable to group and world).
\n
\n  Due to some brain damage on the part of both CentOS and SGE, you may
\n  see errors like 'sh: which: line 1: syntax error: unexpected end of
\n  file'. These errors are harmless.
\n
\n Written by Philip Lijnzaad <plijnzaad@gmail.com>
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
    echo -e $usage_msg
    exit 0
}

if [ $# -eq 0 ]; then usage; fi

# default mail adres:
user=$(id -un)

if [ -f $HOME/.forward ]; then
  mail_address=$user@localhost
else
    mail_address=$(ldapsearch -x "(& (objectClass=person)(uid=$user))" mail |  awk '/^mail/{print $2}')
    if [ -z "$mail_address" ]; then
        echo "$0: Could not find mail address for user $user (consider creating a $HOME/.forward). Exiting." >&2
        exit 5
    fi
fi

opt_M="-M $mail_address"

# default: sent mail when jobs is aborted
opt_m="-m a"

jobname=""
while getopts "N:q:o:e:j:M:m:p:l:h:H:" opt; do
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
            opt_l="-l $OPTARG"          # note: for now, specify as one comma-separed set of resources
            ;;
        M)
            opt_M="-M $OPTARG"
            ;;

### multi-letter qsub options:
        p)
            opt_pe="-pe $OPTARG"        # sorry, -p has to be one-letter options
            opt_pe=$(echo $opt_pe | tr '=' ' ')
            ;;

        h)
            opt_hold_jid="-hold_jid $OPTARG"
            opt_hold_jid=$(echo $opt_hold_jid | tr '=' ' ')
            ;;

        H)
            opt_hold_jid_ad="-hold_jid_ad $OPTARG"
            opt_hold_jid_ad=$(echo $opt_hold_jid_ad | tr '=' ' ')
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
jobname=$(echo $jobname | sed 's/[ /<>|:][ /<>|:]*/_/g;s/^[._]*//;')

## check the length
maxlen=245 # 512 for SGE, but stdout and stderr go to files whose names have length limit of 255 chars ...
if [[ $(echo $jobname  | wc -c) -gt $maxlen ]] ; then
   error "Job name is longer than $maxlen. Use an explicit -N jobname argument
The jobname (which may have been constructed from input arguments) is:
$jobname" 
fi

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
        error "Could not create directory $logdir for stdout and stderr. Job not submitted" >&2
    fi
fi
echo "Using directory $logdir for stdout and stderr" >&2

other_opts="$opt_j $opt_m $opt_l $opt_M $opt_pe $opt_hold_jid $opt_hold_jid_ad"
qsub_opts="-shell no -b yes -cwd -V -o $logdir -e $logdir -q $queue -N $jobname $other_opts"

## finally, submit:
cmd=$(echo "$@")                        # remove one level of quoting ...
qsub $qsub_opts sh -c "(echo -n 'Started ';date)>&2 ; $cmd;  (echo -n 'Ended ';date)>&2"
