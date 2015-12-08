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

import pprint                           # debugging
pp=pprint.PrettyPrinter(indent=2)

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

usage="""Usage: python RemoveOpticalDuplicates.py  ARGS"""

def parseArgs():
    parser=argparse.ArgumentParser(usage=usage, epilog="note that the input file must be a coordinate-sorted SAM file")
    parser.add_argument(flag="--dist", type=int, default=10,
                        ## Use around 100 pixels for later versions of the Illumina software
                        help='Reads closer than this number of pixels are considered optical duplicates')
    parser.add_argument(flag="--input", nargs='?', type=argparse.FileType('r'), default=sys.stdin, help="name of input file in coordinate-sorted SAM format. (default: stdout)")
    parser.add_argument(flag="--uniq", nargs='?', type=argparse.FileType('w'), default=sys.stdout, help="name of output file (will be in SAM format, including the header (default: stdout)")
    parser.add_argument(flag="--repl", nargs='?', type=argparse.FileType('w'), default=sys.stdout, help="name of output file containing the replicates (default: stderr)")
    parser.print_help()
    args=parser.parse_args()
    pp(args)
    return  args
    
def output_best(dups):
    # x,y,mapq, line
    d= sorted(dups, key=lambda x:x[2], revers=True)
    outFile.write(d[0][3])
    if outputDups:
        args.logFile.write("# selected: "  + d[0][3])
        for i in range(1,len(d)-1):
            args.logFile.write("# dup: "  + d[i][3])
        
def output_uniq(dupCands) :
    # dict contains tile+cigar combinations having reads with same start
    # position
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
                if d2 < optDist2:
                    G.add_edge(i,j)
        for comp in nx.connected_components(G):
            if len(comp) ==  1:
                outFile.write(reads[0][3])
                nuniq += 1
            else:
                ndups += len(comp)
                output_best(  [ reads[i] for i in comp ] )
    return (nuniq,ndups)

def removeOpticalDuplicates(param) :
    inFN = param.getInFN()
    outFile = sys.stdout
    optDist = param.getOpticalDuplicatePixelDistance()
    optDist2= optDist*optDist # squared Euclidean distance criterion
    nuniq=0
    ndup=0

    inFile = args.input
    outFile = args.uniq
    logFile = args.repl
    
    dupCands_perPos = {}            # reinitialized every pos
    curChr = ""
    curStartPos = ""

    PG_output = False
    line = args.inFile.readline()
    while (1) :
        if not line: break
        nextLine = args.inFile.readline()
        if line.startswith("@") :
            args.outFile.write(line)
            line = nextLine
            continue
        if not PG_output:
            PG_output=True
            args.outFile.write("@PG\tID:removeOptDups.py\tPN:removeOptDups.py\tVN:0 CL:"+" ".join(param.argv)+"\n")

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
            args.outFile.write(line)
            curChr = chr
            curStartPos = startPos
            continue

        if (not nextLine or chr != curChr or startPos != curStartPos ) :
            ## EOF, or found new non-duplicate, output old ones and start over
            (uniq, dup) = output_unique( dupCands_perPos )
            nuniq += uniq; ndup += dup
            dupCands_perPos = {} 
            dupCands_perPos[tile_cigar] = [ read ]
        else:
            dupCands_perPos[tile_cigar] =  dupCands_perPos.get(tile_cigar, []) # make sure it exists
            dupCands_perPos[tile_cigar].append(read)

        curChr = chr
        curStartPos = startPos
        line = nextLine

    args.inFile.close()
    args.outFile.close()


# Execute as application
if __name__ == '__main__' :
    args=parseArgs()
    removeOpticalDuplicates(args) 

