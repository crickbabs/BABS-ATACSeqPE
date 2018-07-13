#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on June 29th 2018 to annotate merged peaks
#######################################################################
#######################################################################

import os
import argparse

import funcs

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Add sample boolean files and aggregate columns from merged MACS narrow or broad peak file.'
Epilog = """Example usage: python expand_merged_macs.py <MERGED_INTERVAL_FILE> <SAMPLE_NAME_LIST> <OUTFILE> --is_narrow_peak --min_samples 0"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('MERGED_INTERVAL_FILE', help="Merged MACS2 interval file created using linux sort and mergeBed.")
argParser.add_argument('SAMPLE_NAME_LIST', help="Comma-separated list of sample names as named in individual MACS2 broadPeak/narrowPeak output file e.g. SAMPLE_1 for SAMPLE_1_peak_1.")
argParser.add_argument('OUTFILE', help="Full path to output directory.")

## OPTIONAL PARAMETERS
argParser.add_argument('-in', '--is_narrow_peak', dest="IS_NARROW_PEAK", help="Whether merged interval file was generated from narrow or broad peak files (default: False).",action='store_true')
argParser.add_argument('-ms', '--min_samples', type=int, dest="MIN_SAMPLES", default=0, help="Minumum number of samples required per merged peak (default: 0).")
args = argParser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

## MergedIntervalTxtFile is file created using commands below:
## 1) broadPeak
## sort -k1,1 -k2,2n <MACS_BROADPEAK_FILES_LIST> | mergeBed -c 2,3,4,5,6,7,8,9 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > merged_peaks.txt
## 2) narrowPeak
## sort -k1,1 -k2,2n <MACS_NARROWPEAK_FILE_LIST> | mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > merged_peaks.txt

def expand_merged_macs(MergedIntervalTxtFile,SampleNameList,OutFile,isNarrow=False,minSamples=0):

    funcs.makedir(os.path.dirname(OutFile))

    SampleNameList = sorted(SampleNameList)
    totalInIntervals = 0; totalOutIntervals = 0
    sampleFreqDict = {}; combFreqDict = {}
    fin = open(MergedIntervalTxtFile,'r')
    fout = open(OutFile,'w')
    oFields = ['chr','start','end','interval_id','num_peaks','num_samples'] + [x+'.bool' for x in SampleNameList] + [x+'.fc' for x in SampleNameList] + [x+'.qval' for x in SampleNameList] + [x+'.pval' for x in SampleNameList] + [x+'.start' for x in SampleNameList] + [x+'.end' for x in SampleNameList]
    if isNarrow:
        oFields += [x+'.summit' for x in SampleNameList]
    fout.write('\t'.join(oFields) + '\n')
    while True:
        line = fin.readline()
        if line:
            lspl = line.strip().split('\t')

            chromID = lspl[0]; mstart = int(lspl[1]); mend = int(lspl[2]);
            starts = [int(x) for x in lspl[3].split(',')]; ends = [int(x) for x in lspl[4].split(',')]
            names = lspl[5].split(','); fcs = [float(x) for x in lspl[8].split(',')]
            pvals = [float(x) for x in lspl[9].split(',')]; qvals = [float(x) for x in lspl[10].split(',')]
            summits = []
            if isNarrow:
                summits = [int(x) for x in lspl[11].split(',')]

            fcDict = {}; qvalDict = {}; pvalDict = {}; startDict = {}; endDict = {}; summitDict = {}
            for idx in range(len(names)):
                sample = '_'.join(names[idx].split('_')[:-2])
                if not fcDict.has_key(sample):
                    fcDict[sample] = []
                fcDict[sample].append(str(fcs[idx]))
                if not qvalDict.has_key(sample):
                    qvalDict[sample] = []
                qvalDict[sample].append(str(qvals[idx]))
                if not pvalDict.has_key(sample):
                    pvalDict[sample] = []
                pvalDict[sample].append(str(pvals[idx]))
                if not startDict.has_key(sample):
                    startDict[sample] = []
                startDict[sample].append(str(starts[idx]))
                if not endDict.has_key(sample):
                    endDict[sample] = []
                endDict[sample].append(str(ends[idx]))
                if isNarrow:
                    if not summitDict.has_key(sample):
                        summitDict[sample] = []
                    summitDict[sample].append(str(summits[idx]))

            samples = sorted(fcDict.keys())
            numSamples = len(samples)
            if numSamples >= minSamples:
                boolList  = ['TRUE' if x in samples else 'FALSE' for x in SampleNameList]
                fcList = [';'.join(fcDict[x]) if x in samples else 'NA' for x in SampleNameList]
                qvalList = [';'.join(qvalDict[x]) if x in samples else 'NA' for x in SampleNameList]
                pvalList = [';'.join(pvalDict[x]) if x in samples else 'NA' for x in SampleNameList]
                startList = [';'.join(startDict[x]) if x in samples else 'NA' for x in SampleNameList]
                endList = [';'.join(endDict[x]) if x in samples else 'NA' for x in SampleNameList]
                oList = [str(x) for x in [chromID,mstart,mend,'Interval_'+str(totalOutIntervals+1),len(names),numSamples]+boolList+fcList+qvalList+pvalList+startList+endList]
                if isNarrow:
                    oList += [';'.join(summitDict[x]) if x in samples else 'NA' for x in SampleNameList]
                fout.write('\t'.join(oList) + '\n')

                if not sampleFreqDict.has_key(numSamples):
                    sampleFreqDict[numSamples] = 0
                sampleFreqDict[numSamples] += 1
                totalOutIntervals += 1

            tsamples = tuple(sorted(samples))
            if not combFreqDict.has_key(tsamples):
                combFreqDict[tsamples] = 0
            combFreqDict[tsamples] += 1
            totalInIntervals += 1
        else:
            fin.close()
            fout.close()
            break

    fout = open(OutFile[:-4]+'.sampleCounts.log','w')
    for k,v in sorted(sampleFreqDict.items(),reverse=True):
        fout.write('%s\t%s\n' % (k,v))
    fout.write('\n\n')
    combFreqItems = sorted([(combFreqDict[x],x) for x in combFreqDict.keys()],reverse=True)
    for k,v in combFreqItems:
        fout.write('%s\t%s\n' % (k,','.join(v)))
    fout.close()

    print 'Total input intervals = ' + str(totalInIntervals)
    print 'Total output intervals = ' + str(totalOutIntervals)

############################################
############################################
## RUN FUNCTION
############################################
############################################

expand_merged_macs(MergedIntervalTxtFile=args.MERGED_INTERVAL_FILE,SampleNameList=args.SAMPLE_NAME_LIST.split(','),OutFile=args.OUTFILE,isNarrow=args.IS_NARROW_PEAK,minSamples=args.MIN_SAMPLES)

############################################
############################################
############################################
############################################
