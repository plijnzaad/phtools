#!/usr/local/bin/python

"Remove optical duplicates from a SAM file"
"Original by Anna Salzberg at https://gist.github.com/annasa/eef7c30152ac296bb49b"
"Additions by plijnzaad@gmail.com"

import sys
import string
import os
import time
import math
import re

def convertStr(s):
    """Convert string to either int or float."""
    try:
        ret = int(s)
    except ValueError:
        try :
            ret = float(s);
        except ValueError :
            ret = -1;
                      
    return ret

class Parameter :
    def __init__(self, argv) :
        self.argv = argv;
        self.inFN = "";
        self.opticalDuplicatePixelDistance = 10;

    def getOpticalDuplicatePixelDistance(self) :
        return self.opticalDuplicatePixelDistance 

    def getInFN(self) :
        return self.inFN

    def getOutFN(self) :
        return self.outFN;    

    # Prints usage
    def usage(self):
        print "Usage: python RemoveOpticalDuplicates.py  [-d pixels]  [ <input file> ] " 
        print "where the input is a sam file sorted by chr position.  Default optical_duplciate_pixel_distance: 10"
        print "    (Note: Use around 100 pixels for later versions of the Illumina software)"
        print "Example: python RemoveOpticalDuplicates.py  sampleA.bwa.sorted.sam  >  sampleA.bwa.rmoptdup.sam"
        print "The -d argument specifies the minimum x or y-difference (in pixels). Any pair of read closer than this is considered duplicates"

        sys.exit(1)
    
    def checkArgs(self) :
        i = 1;
        argc = len(sys.argv);
        errorMsg = "";
        while i < argc and errorMsg == "":
            if self.argv[i] == '-d' :
                i += 1;
                if i < argc and not self.argv[i].startswith('-') :
                    self.opticalDuplicatePixelDistance = convertStr(self.argv[i]);
                    if self.opticalDuplicatePixelDistance < 1 :
                        errorMsg = 'Error: optical duplicate pixel distance must be an integer >= 1'
                else :
                    errorMsg = 'Error: optical duplicate pixel distance expected after -d flag'
            elif i == argc-1 :
                self.inFN = self.argv[i];
            else :
                errorMsg = 'Error: unknown flag: ' + self.argv[i];
                
            i += 1;

        if errorMsg == "" and self.inFN != "" and \
                                         not re.match(r'.*\.sam$', self.inFN):
            errorMsg = 'Need a SAM file'
        
        if errorMsg <> "" :
            print errorMsg;
            self.usage();


def removeOpticalDuplicates(param) :
    inFN = param.getInFN()
    outFile = sys.stdout
    optDist = param.getOpticalDuplicatePixelDistance()

    if inFN == "":
        inFile = sys.stdin
    else:
        try:
            inFile = open(inFN, 'r')
        except:
            print "Unable to open file " + inFN
            sys.exit(-1)

    tile_cigarToDupDict = {} 
    curChr = ""
    curStartPos = ""

    first = True
    PG_output = False
    line = inFile.readline()
    while (1) :
        if not line: break
        nextLine = inFile.readline()
        if line.startswith("@") :
            outFile.write(line)
        else:
            if not PG_output:
                PG_output=True
                outFile.write("@PG\tID:removeOptDups.py\tPN:removeOptDups.py\tVN:0 CL:"+" ".join(param.argv)+"\n")
            tokens = line.split()
            chr = tokens[2]
            startPos = tokens[3]
            mapq = tokens[4] 
            cigar = tokens[5]
            subtokens = tokens[0].split(":")
            machine_name = subtokens[0]
            run_name = subtokens[1]
            flowcell_id = subtokens[2]
            flowcell_lane = subtokens[3]
            tile = subtokens[4]
            x = subtokens[5]
            y = subtokens[6]
            tile_cigar = tile + "_" + cigar; 
            if first :
                curChr = chr
                curStartPos = startPos
                first = False;

            if chr <> "*" and (chr == curChr and startPos == curStartPos) :
                    dups = tile_cigarToDupDict.get(tile_cigar, []);
                    dups.append([x, y, mapq, line]);
                    tile_cigarToDupDict[tile_cigar] = dups;

            if chr <> "*" and (chr <> curChr or startPos <> curStartPos or not nextLine) : 
                for key in tile_cigarToDupDict.keys() :
                   dups = tile_cigarToDupDict[key];
                   dups.sort();
                   if len(dups) == 1 :
                      nondupline = dups[0][3]
                      outFile.write(nondupline);
                   elif len(dups) > 1 :
                       processed = []
                       for k in range(0, len(dups)) :
                          processed.append(False)

                       for i in range(0, len(dups)) :
                          if not processed[i] : 
                              xi = convertStr(dups[i][0])
                              yi = convertStr(dups[i][1])
                              mapqi = convertStr(dups[i][2])
                              processed[i] = True;
                        
                              best = i;
                              bestMapq = mapqi;
                              for j in range(i+1, len(dups)) :
                                  if not processed[j] : 
                                      xj = convertStr(dups[j][0])
                                      yj = convertStr(dups[j][1])
                                      mapqj = convertStr(dups[j][2])
                       
                                      if abs(xj-xi) > optDist : break; # no dup
                                      if abs(yj-yi) <= optDist : # dup, keep best one
                                          processed[j] = True 
                                          if mapqj > bestMapq :
                                              # (add output of the 'loosing' line here?)
                                              best = j;
                                              bestMapq = mapqj;
                              bestline = dups[best][3]
                              outFile.write(bestline);
                tile_cigarToDupDict = {} 
                tile_cigarToDupDict[tile_cigar] = [[x, y, mapq, line]]

            if chr <> "*" and not nextLine and (chr <> curChr or startPos <> curStartPos) : 
                outFile.write(line);
            if chr == "*" :
                outFile.write(line);
            curChr = chr
            curStartPos = startPos;

        line = nextLine;
                
    inFile.close();
    outFile.close();


# Execute as application
if __name__ == '__main__' :
    param = Parameter(sys.argv)
    param.checkArgs();

    removeOpticalDuplicates(param) 

