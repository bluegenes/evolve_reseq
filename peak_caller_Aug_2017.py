#!/usr/bin/env python
#######################################
### Tessa Pierce (code) and Thiago Lima (concept)###
### 06/02/2016 ###
#######################################
#try this: https://scipher.wordpress.com/2010/12/02/simple-sliding-window-iterator-in-python/
import sys
import math
import argparse

def call_peaks(inFile,outFile,logPthreshBase,logPthreshTop, minSNPs):
    out = open(outFile, 'w')
    out.write("Scaffold"+'\t'+"Pos"+'\t'+"SD_rep1_count"+'\t'+"SC_rep1_count"+'\t'+"SD_rep2_count"+'\t'+"SC_rep2_count"+'\t'+"SD_rep3_count"+'\t'+"SC_rep3_count"+'\t'+"p"+'\t'+"logP"+'\t'+'Pop'+ '\t' +'Peak_ID' + '\n')
    with open(inFile, 'r') as f:
        next(f)
        previousLine,sdPeak,scnPeak  = [],[],[]
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
                sdPeakCount, scnPeakCount = writePeak(sdPeak, scnPeak, out, minSNPs, logPthreshTop,sdPeakCount,scnPeakCount)
                sdPeak,scnPeak = [],[]
            # find SD peaks
            if ((rep1_SD > rep1_SC) + (rep2_SD > rep2_SC) + (rep3_SD > rep3_SC) >=2)  and logP > logPthreshBase:
                if len(scnPeak) >= 1:
                # if we were on an SCN peak, write it out and reset scnPeak list to []
                    sdPeakCount, scnPeakCount = writePeak(sdPeak, scnPeak, out, minSNPs, logPthreshTop,sdPeakCount,scnPeakCount)
                    scnPeak = []
                # ok, now add this line to SD peak list
                sdPeak.append(currentLine)

            elif ((rep1_SD < rep1_SC) + (rep2_SD < rep2_SC) + (rep3_SD < rep3_SC) >=2) and logP > logPthreshBase:
                if len(sdPeak) >= 1:
                # if we were on an SD peak, write it out and reset sdPeak to []
                    sdPeakCount, scnPeakCount = writePeak(sdPeak, scnPeak, out, minSNPs, logPthreshTop,sdPeakCount,scnPeakCount)
                    sdPeak =[]
                scnPeak.append(currentLine) # append current line to scn peak
            else: # current line is not significant
                info_to_add = ['NS', 'NA']
                # write out any previous peaks we had, and reset peak lists to []
                sdPeakCount, scnPeakCount = writePeak(sdPeak, scnPeak, out, minSNPs, logPthreshTop,sdPeakCount,scnPeakCount)
                scnPeak, sdPeak = [],[]
                # write out current NS line
                out.write('\t'.join(map(str,currentLine+info_to_add)) + '\n')
            previousLine = currentLine
        # catch the last peak
        sdPeakCount, scnPeakCount = writePeak(sdPeak, scnPeak, out, minSNPs, logPthreshTop,sdPeakCount,scnPeakCount)
    out.close()

def writePeak(sdPeakSNPs, scnPeakSNPs, outF, minSNPs, logPthreshTop,sdPeakCount,scnPeakCount):
    if (len(sdPeakSNPs) >0 and len(scnPeakSNPs) > 0):
        sys.stdout.write("improperly keeping track of peaks")
    elif len(sdPeakSNPs) > 0:
        if len(sdPeakSNPs) >= minSNPs and any([True for x in sdPeakSNPs if x[-1] >=logPthreshTop]):
            sdPeakCount+=1
            info_to_add = ['SD','SD_' + str(sdPeakCount)]
        else:
            info_to_add = ['NS', 'NA']
        for entry in sdPeakSNPs:
            outF.write('\t'.join(map(str,entry+info_to_add)) + '\n')

    elif len(scnPeakSNPs) > 0:
        if len(scnPeakSNPs) >= minSNPs and any([True for x in scnPeakSNPs if x[-1] >=logPthreshTop]):
            scnPeakCount+=1
            info_to_add = ['SC','SC_' + str(scnPeakCount)]
        else:
            info_to_add = ['NS', 'NA']
        for entry in scnPeakSNPs:
            outF.write('\t'.join(map(str,entry+info_to_add)) + '\n')
    return sdPeakCount,scnPeakCount


if(__name__=='__main__'):
    parser = argparse.ArgumentParser(description="Call Peaks post smoothing with sliding window")
    parser.add_argument('-f', '--inFile', help='input File')
    parser.add_argument('-o', '--outFile', help='output File')
    parser.add_argument('-n', '--min', help='minimum # of SNPs in a peak', default=1)
    parser.add_argument('-t', '--logPthreshTop', help='log pval threshold of top of peak', default=2)
    parser.add_argument('-b', '--logPthreshBase', help='log pval threshold of base of peak', default=1)
    args = parser.parse_args()
    call_peaks(args.inFile, args.outFile, args.logPthreshBase, args.logPthreshTop, args.min)
