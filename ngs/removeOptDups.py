#!/usr/local/bin/python

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

print "Untested, and prolly broken"
sys.exit(-1)

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

class Parameter :
    def __init__(self, argv) :
        self.argv = argv
        self.inFN = ""
        self.opticalDuplicatePixelDistance = 10

    def getOpticalDuplicatePixelDistance(self) :
        return self.opticalDuplicatePixelDistance 

    def getInFN(self) :
        return self.inFN

    def getOutFN(self) :
        return self.outFN

    # Prints usage
    def usage(self):
        print "Usage: python RemoveOpticalDuplicates.py  [-d pixels]  [ <input file> ] " 
        print "where the input is a sam file sorted by chr position.  Default optical_duplicate_pixel_distance: 10"
        print "    (Note: Use around 100 pixels for later versions of the Illumina software)"
        print "Example: python RemoveOpticalDuplicates.py  sampleA.bwa.sorted.sam  >  sampleA.bwa.rmoptdup.sam"
        print "The -d argument specifies the minimum x or y-difference (in pixels). Any pair of read closer than this is considered duplicates"

        sys.exit(1)
    
    def checkArgs(self) :
        i = 1
        argc = len(sys.argv)
        errorMsg = ""
        while i < argc and errorMsg == "":
            if self.argv[i] == '-d' :
                i += 1
                if i < argc and not self.argv[i].startswith('-') :
                    self.opticalDuplicatePixelDistance = convertStr(self.argv[i])
                    if self.opticalDuplicatePixelDistance < 1 :
                        errorMsg = 'Error: optical duplicate pixel distance must be an integer >= 1'
                else :
                    errorMsg = 'Error: optical duplicate pixel distance expected after -d flag'
            elif i == argc-1 :
                self.inFN = self.argv[i]
            else :
                errorMsg = 'Error: unknown flag: ' + self.argv[i]
                
            i += 1

        if errorMsg == "" and self.inFN != "" and \
                                         not re.match(r'.*\.sam$', self.inFN):
            errorMsg = 'Need a SAM file'
        
        if errorMsg <> "" :
            print errorMsg
            self.usage()


def output_best_replicates(dupCands) :
    # dict contains tile+cigar combinations having reads with same start position
    ndups=0
    for tile in dupCands.keys() :
        reads = dupCands[tile]
        best = find_best_replicate(reads)
        outFile.write(dups[best][3])
        # ndups += ???
        ## find connected components among set of reads; for each of them, output the best one ...
        for i in range(0, len(dups)) :  # len(dups) can just be 1
            xi = convertStr(dups[i][0])
            yi = convertStr(dups[i][1])
            mapqi = convertStr(dups[i][2])
            
            best = i
            bestMapq = mapqi
            for j in range(i+1, len(dups)) :
                xj = convertStr(dups[j][0])
                yj = convertStr(dups[j][1])
                d2 = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi)
                if d2 > optDist2:
                    continue
                mapqj = convertStr(dups[j][2])
                if mapqj > bestMapq :
                    if outputDups: # @@add this option
                        dupFile.write( dups[best][3] )
                        best = j
                        bestMapq = mapqj
            output best j for i?
         output best pair of all i,j?
    ndups=0
    return ndups

def removeOpticalDuplicates(param) :
    inFN = param.getInFN()
    outFile = sys.stdout
    optDist = param.getOpticalDuplicatePixelDistance()
    optDist2= optDist*optDist # squared Euclidean distance criterion

    if inFN == "":
        inFile = sys.stdin
    else:
        try:
            inFile = open(inFN, 'r')
        except:
            print "Unable to open file " + inFN
            sys.exit(-1)

    dupCands_perTile = {}            # reinitialized every line
    curChr = ""
    curStartPos = ""

    PG_output = False
    line = inFile.readline()
    while (1) :
        if not line: break
        nextLine = inFile.readline()
        if line.startswith("@") :
            outFile.write(line)
            line = nextLine
            continue
        if not PG_output:
            PG_output=True
            outFile.write("@PG\tID:removeOptDups.py\tPN:removeOptDups.py\tVN:0 CL:"+" ".join(param.argv)+"\n")

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
        read = [x, y, mapq, line]

        if chr == "*" :
            outFile.write(line)
            curChr = chr
            curStartPos = startPos
            continue

        if (not nextLine or chr <> curChr or startPos <> curStartPos ) :
            ## EOF, or found new non-duplicate, output potential old ones and start over
            ndups = ndups + output_best_replicates( dupCands_perTile )
            dupCands_perTile = {} 
            dupCands_perTile[tile_cigar] = [ read ]
        else:
            dupCands_perTile[tile_cigar] =  dupCands_perTile.get(tile_cigar, []) # make sure it exists
            dupCands_perTile[tile_cigar].append(read)

        curChr = chr
        curStartPos = startPos
        line = nextLine

    inFile.close()
    outFile.close()


# Execute as application
if __name__ == '__main__' :
    param = Parameter(sys.argv)
    param.checkArgs()

    removeOpticalDuplicates(param) 

