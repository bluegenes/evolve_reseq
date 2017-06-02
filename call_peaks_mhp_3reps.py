#!/usr/bin/env python
#######################################
### Tessa Pierce and Thiago Lima ###
### 06/02/2016 ###
#######################################

import sys
import math
import argparse


def call_peaks(inFile, outFile):
    out = open(outFile, 'w')
    out.write("Scaffold"+'\t'+"Pos"+'\t'+"SD_rep1_count"+'\t'+"SC_rep1_count"+'\t'+"SD_rep2_count"+'\t'+"SC_rep2_count"+'\t'+"SD_rep3_count"+'\t'+"SC_rep3_count"+'\t'+"p"+'\t'+"logP"+'\t'+'Pop'+ '\t' +'Peak_ID' + '\n')
    minSNPs = 10
    maxDist = 1000
    logPthreshTop = 5
    logPthreshBase = 2
    with open(inFile, 'r') as f:
        next(f)
        previousLine,previousSC,previousSD,sdPeak,scnPeak  = [],[],[],[],[]
        sdPeakCount,scnPeakCount = 0,0
        for line in f:
            splitLine = line.strip().split('\t')
            Scaffold = splitLine[0]
            Position = int(splitLine[1])
            rep1_SD,rep1_SC = float(splitLine[2]),float(splitLine[3])
            rep2_SD,rep2_SC  = float(splitLine[6]),float(splitLine[7])
            rep3_SD,rep3_SC = float(splitLine[10]),float(splitLine[11])
            pval = float(splitLine[14])
            logP = -(math.log10(pval))
            currentLine = [Scaffold,Position,rep1_SD,rep1_SC,rep2_SD,rep2_SC,rep3_SD,rep3_SC,pval,logP]
            # if we're on a new scaffold, print any peaks that we found but haven't been written out yet
            if  len(previousLine) > 1 and Scaffold != previousLine[0]:
                #any sd Peak
                if len(sdPeak) >= minSNPs and any([True for x in sdPeak if x[-1] >=logPthreshTop]):
                    sdPeakCount+=1
                    info_to_add = ['SD','SD_' + str(sdPeakCount)]
                else:
                    info_to_add = ['NS']
                for entry in sdPeak:
                    out.write('\t'.join(map(str,entry + info_to_add)) + '\n')
                #scn peaks
                if len(scnPeak) >= minSNPs and any([True for x in scnPeak if x[-1] >=logPthreshTop]):
                    scnPeakCount+=1
                    info_to_add = ['SC','SC_' + str(scnPeakCount)]
                else:
                    info_to_add = ['NS']
                for entry in scnPeak:
                    out.write('\t'.join(map(str,entry + info_to_add)) + '\n')
                sdPeak,scnPeak,previousSD,previousSC = [],[],[],[]
            # find SD peaks
            if ((rep1_SD > rep1_SC) + (rep2_SD > rep2_SC) + (rep3_SD > rep3_SC) >=2)  and logP > logPthreshBase:
                if len(previousSD) > 1:
                    if Position - previousSD[1] >= maxDist:
                        if len(sdPeak) >= minSNPs and any([True for x in sdPeak if x[-1] >=logPthreshTop]):
                            sdPeakCount+=1
                            info_to_add = ['SD','SD_' + str(sdPeakCount)]
                        else:
                            info_to_add = ['NS']
                        for entry in sdPeak:
                            out.write('\t'.join(map(str,entry + info_to_add)) + '\n')
                        previousSD,sdPeak = [],[]
                sdPeak.append(currentLine)
                previousSD = currentLine
            #find SC peaks
            elif ((rep1_SD < rep1_SC) + (rep2_SD < rep2_SC) + (rep3_SD < rep3_SC) >=2) and logP > logPthreshBase:
                if len(previousSC) > 1: #not first count --> we can check for distance to last position
                    if Position - previousSC[1] >= maxDist:
                        if len(scnPeak) >= minSNPs and any([True for x in scnPeak if x[-1] >=logPthreshTop]):
                            scnPeakCount+=1
                            info_to_add = ['SC','SC_' + str(scnPeakCount)]
                        else:
                            info_to_add = ['NS']
                        for entry in scnPeak:
                            out.write('\t'.join(map(str,entry+info_to_add)) + '\n')
                        previousSC,scnPeak = [],[] #reset
                scnPeak.append(currentLine)
                previousSC = currentLine
            else:
                info_to_add = ['NS']
                out.write('\t'.join(map(str,currentLine+info_to_add)) + '\n')
            previousLine = currentLine

        if len(sdPeak) >= minSNPs and any([True for x in sdPeak if x[-1] >=logPthreshTop]):
            sdPeakCount+=1
            info_to_add = ['SD','SD_' + str(sdPeakCount)]
        else:
            info_to_add = ['NS']
        for entry in sdPeak:
            out.write('\t'.join(map(str,entry + info_to_add)) + '\n')
        #scn peaks
        if len(scnPeak) >= minSNPs and any([True for x in scnPeak if x[-1] >=logPthreshTop]):
            scnPeakCount+=1
            info_to_add = ['SC','SC_' + str(scnPeakCount)]
        else:
            info_to_add = ['NS']
        for entry in scnPeak:
            out.write('\t'.join(map(str,entry + info_to_add)) + '\n')
#        sdPeak,scnPeak,previousSD,previousSC = [],[],[],[]
    out.close()


if(__name__=='__main__'):
    parser = argparse.ArgumentParser(description="Call Peaks")
    parser.add_argument('-f', '--inFile', help='input File')
    parser.add_argument('-o', '--outFile', help='output File')
    args = parser.parse_args()
    call_peaks(args.inFile, args.outFile)
