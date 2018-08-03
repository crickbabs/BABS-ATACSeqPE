#!/usr/bin/env python

###############################################################################
###############################################################################
## Created on February 1st 2017 to remove singletons from paired-end BAM file
###############################################################################
###############################################################################

import os
import pysam
import argparse

import funcs

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Remove singleton reads from paired-end BAM file i.e if read1 is present in BAM file without read 2 and vice versa.'
Epilog = """Example usage: bampe_rm_orphan.py <BAM_INPUT_FILE> <BAM_OUTPUT_FILE>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('BAM_INPUT_FILE', help="Input BAM file sorted by name.")
argParser.add_argument('BAM_OUTPUT_FILE', help="Output BAM file sorted by name.")

argParser.add_argument('-edc', '--excl_diff_chrom', dest="EXCL_DIFF_CHROM", help="Remove pairs that map to different chromosomes.",action='store_true')
args = argParser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

def bampe_rm_orphan(BAMIn,BAMOut,ExclDiffChrom=False):

    ## SETUP DIRECTORY/FILE STRUCTURE
    OutDir = os.path.dirname(BAMOut)
    funcs.makedir(OutDir)

    ## COUNT VARIABLES
    TotalReads = 1; TotalOutputPairs = 0; TotalSingletons = 0; TotalDiffChromPairs = 0

    ## ITERATE THROUGH SAM FILE
    EOF = 0
    SAMFin = pysam.Samfile(BAMIn,"rb")
    SAMFout = pysam.Samfile(BAMOut, "wb",header=SAMFin.header)
    currRead = SAMFin.next()
    for read in SAMFin:
        TotalReads += 1
        if currRead.qname == read.qname:
            pair1 = currRead; pair2 = read
            if ExclDiffChrom:
                if pair1.tid == pair2.tid:
                    TotalOutputPairs += 1
                    SAMFout.write(pair1)
                    SAMFout.write(pair2)
                else:
                    TotalDiffChromPairs += 1
            else:
                TotalOutputPairs += 1
                SAMFout.write(pair1)
                SAMFout.write(pair2)

            ## RESET COUNTER
            try:
                currRead = SAMFin.next()
            except:
                StopIteration
                EOF = 1

        ## READS WHERE ONLY ONE OF A PAIR IS IN FILE
        else:
            TotalSingletons += 1
            pair1 = currRead
            currRead = read

    if not EOF:
        TotalSingletons += 1
        pair1 = currRead

    ## CLOSE ALL FILE HANDLES
    SAMFin.close()
    SAMFout.close()

    LogFile = os.path.join(OutDir,'%s_filterBAMPE_ExclSingletons.log' % (os.path.basename(BAMOut[:-4])))
    SamLogFile = open(LogFile,'w')
    SamLogFile.write('\n##############################\n')
    SamLogFile.write('FILES/DIRECTORIES')
    SamLogFile.write('\n##############################\n\n')
    SamLogFile.write('Input File: ' + BAMIn + '\n')
    SamLogFile.write('Output File: ' + BAMOut + '\n')
    SamLogFile.write('\n##############################\n')
    SamLogFile.write('OVERALL COUNTS')
    SamLogFile.write('\n##############################\n\n')
    SamLogFile.write('Total Input Reads = ' + str(TotalReads) + '\n')
    SamLogFile.write('Total Output Pairs = ' + str(TotalOutputPairs) + '\n')
    SamLogFile.write('Total Singletons Excluded = ' + str(TotalSingletons) + '\n')
    SamLogFile.write('Total Pairs Mapped To Different Chromosomes Excluded = ' + str(TotalDiffChromPairs) + '\n')
    SamLogFile.write('\n##############################\n')
    SamLogFile.close()

############################################
############################################
## RUN FUNCTION
############################################
############################################

bampe_rm_orphan(BAMIn=args.BAM_INPUT_FILE,BAMOut=args.BAM_OUTPUT_FILE,ExclDiffChrom=args.EXCL_DIFF_CHROM)

############################################
############################################
############################################
############################################
