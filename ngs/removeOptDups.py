#!/bin/env python3
### #!/usr/local/bin/python

## Remove optical duplicates from a SAM file
## Original by Anna Salzberg at https://gist.github.com/annasa/eef7c30152ac296bb49b
## See also http://sourceforge.net/p/samtools/mailman/message/32915073/
##  (http://samtools-help.narkive.com/iU2fiVxO/reporting-bug-optical-duplicates-of-picard-markduplicates)
##
## Additions by plijnzaad@gmail.com

## (it assumes a the reads to be on a single tile?)

import sys
import string
import os
import time
import math
import re
import networkx as nx
import argparse

# debugging:
import pdb
import pprint
pp=pprint.PrettyPrinter(indent=2).pprint # call as pp(some.object)

global ARGS

def convertStr(s):
    """Convert string to either int or float."""
    try:
        ret = int(s)
    except ValueError:
        try :
            ret = float(s)
        except ValueError :
            ret = -1
                      
    return ret

usage="""Usage: python RemoveOpticalDuplicates.py  ARGS

Given a SAM file (sorted by position), removes all but one of the optical duplicates.

"""
epilog="Written by plijnzaad@gmail.com, based on original by Anna Salzberg"

def parseArgs():
    parser=argparse.ArgumentParser(usage=usage, epilog=epilog,
                                   add_help=False)
    parser.add_argument("--dist", dest='dist', type=int, default=10,
                        ## Use around 100 pixels for later versions of the Illumina software
                        help='Reads closer than this number of pixels are considered optical duplicates')
    parser.add_argument("--input", dest="inFile", nargs='?', type=argparse.FileType('r'), default=sys.stdin, help="name of input file in coordinate-sorted SAM format. (default: stdin)")
    parser.add_argument("--out", dest="outFile", nargs='?', type=argparse.FileType('w'), default=sys.stdout, help="name of output file (will be in SAM format, including the header (default: stdout)")
    parser.add_argument("--log", dest="logFile", nargs='?', type=argparse.FileType('w'), default=sys.stderr, help="name of log file. Contains the replicates (default: stderr). Use /dev/null to supress")
    parser.add_argument("--justdist", action="store_true", help="only output the distances within each position")
    parser.add_argument("-h", "--help", action="store_true", help="obvious no?")
    args=parser.parse_args()
    args.dist2 = args.dist*args.dist
    if args.help:
        parser.print_help()
        sys.exit(1)
    args.argv=sys.argv
    return args
## parseArgs

def output_best(dups):
    # x,y,mapq, line
    d= sorted(dups, key=lambda x:x[2], reverse=True)
    ARGS.outFile.write(d[0][3])
    ARGS.logFile.write("# selected: "  + d[0][3]) # use /dev/null if unwanted
    for i in range(1,len(d)-1):
        ARGS.logFile.write("# duplicate: "  + d[i][3])
## output_best

def output_uniq(dupCands):
    # dict contains tile+cigar combinations having reads with same start
    # position. dist2 is the squared distance (for faster comparison)
    noutput=0
    nuniq=0
    ndups=0
    for tile in dupCands.keys() :
        reads = dupCands[tile]          # x, y, mapq, line
        G=nx.Graph()
        G.add_nodes_from(range(len(reads)))
        for i in range(0, len(reads)) :  # len(reads) can just be 1
            xi = reads[i][0]
            yi = reads[i][1]
            for j in range(i+1, len(reads)) :
                xj = reads[j][0]
                yj = reads[j][1]
                d2 = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi)
                if(ARGS.justdist):
                    ARGS.outFile.write("%.2f\n" % math.sqrt(d2))
                if d2 < ARGS.dist2:
                    G.add_edge(i,j)
        if not ARGS.justdist:
            for comp in nx.connected_components(G):
                noutput +=1
                if len(comp) ==  1:
                    ARGS.outFile.write(reads[0][3])
                    nuniq += 1
                else:
                    ndups += len(comp)
                output_best(  [ reads[i] for i in comp ] )
    return (noutput,nuniq,ndups)
## output_uniq


def main():
    pdb.set_trace()
    nout=0
    nuniq=0
    ndup=0
    dupCands_perPos = {}            # reinitialized every pos
    curChr = ""
    curStartPos = ""
    PG_output = False
    line = ARGS.inFile.readline()
    while (1) :
        if not line: break
        nextLine = ARGS.inFile.readline()
        if line.startswith("@") :
            ARGS.outFile.write(line)
            line = nextLine
            continue
        if not PG_output:
            PG_output=True
            ARGS.outFile.write("@PG\tID:removeOptDups.py\tPN:removeOptDups.py\tVN:0 CL:"+ " ".join(ARGS.argv)+"\n")
        tokens = line.split()
        chr = tokens[2]
        startPos = tokens[3]
        mapq = tokens[4] 
        cigar = tokens[5]
        subtokens = tokens[0].split(":")
        # machine_name = subtokens[0]; run_name = subtokens[1]; flowcell_id = subtokens[2]; lane = subtokens[3]
        tile = subtokens[4]
        x = subtokens[5]
        y = subtokens[6]
        tile_cigar = tile + "_" + cigar
        read = [ convertStr(x), convertStr(y), mapq, line]
        if chr == "*" :
            ARGS.outFile.write(line)
            curChr = chr
            curStartPos = startPos
            continue
        if (not nextLine or chr != curChr or startPos != curStartPos ) :
            ## EOF, or found new non-duplicate, output old ones and start over
            (n, u, d) = output_uniq( dupCands_perPos)
            nout+= n; nuniq += u; ndup += d
            dupCands_perPos = {} 
            dupCands_perPos[tile_cigar] = [ read ]
        else:
            dupCands_perPos[tile_cigar] =  dupCands_perPos.get(tile_cigar, []) # make sure it exists
            dupCands_perPos[tile_cigar].append(read)
        curChr = chr
        curStartPos = startPos
        line = nextLine
    ARGS.logFile.write("Wrote %d reads, found %d unique, and %d replicates\n" % (nout, nuniq, ndup))
    ARGS.inFile.close()
    ARGS.outFile.close()
    ARGS.logFile.close()
## main
    
# Execute as application
if __name__ == '__main__' :
    ARGS=parseArgs()
    pp(ARGS)
    main() 

