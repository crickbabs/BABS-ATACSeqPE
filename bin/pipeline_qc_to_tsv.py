#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on August 6th 2018 to create QC metrics file
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

Description = 'Create a tab-delimited file with various QC metrics generated by BABS-ATACSeqPE pipeline. This will be specific to directory structure of BABS-ATACSeqPE nextflow pipeline.'
Epilog = """Example usage: python pipeline_qc_to_tsv.py <RESULTS_DIR> <MITO_NAME> <OUT_FILE>"""
argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('RESULTS_DIR', help="Results directory. The directory structure used to find files will be specific to BABS-ATACSeqPE nextflow pipeline.")
argParser.add_argument('MITO_NAME', help="Name of Mitochondrial chomosome in genome fasta (e.g. chrM).")
argParser.add_argument('OUT_FILE', help="Path to output file.")

args = argParser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

def pipeline_qc_to_tsv(ResultsDir,MitoName,OutFile):

    funcs.makedir(os.path.dirname(OutFile))

    fileInfoList = [('RUN-LEVEL', 'cutadapt', '', ResultsDir, '.cutadapt.log'),
                    ('RUN-LEVEL', 'flagstat', 'unfiltered', ResultsDir, '.mkD.sorted.bam.flagstat'),
                    ('RUN-LEVEL', 'idxstats', 'unfiltered', ResultsDir, '.mkD.sorted.bam.idxstats'),
                    ('RUN-LEVEL', 'picard_insert_metrics', 'unfiltered', ResultsDir, '.mkD.CollectMultipleMetrics.insert_size_metrics'),
                    ('RUN-LEVEL', 'flagstat', 'filter', ResultsDir, '.mkD.flT.sorted.bam.flagstat'),
                    ('RUN-LEVEL', 'idxstats', 'filter', ResultsDir, '.mkD.flT.sorted.bam.idxstats'),

                    ('REPLICATE-LEVEL', 'flagstat', '', os.path.join(ResultsDir,'align/mergeReplicate/'), '.mRp.rmD.sorted.bam.flagstat'),
                    ('REPLICATE-LEVEL', 'macs2', '', os.path.join(ResultsDir,'align/mergeReplicate/'), '_peaks.broadPeak'),
                    ('REPLICATE-LEVEL', 'frip', '', os.path.join(ResultsDir,'align/mergeReplicate/'), '_peaks.frip.txt'),

                    ('SAMPLE-LEVEL', 'flagstat', '', os.path.join(ResultsDir,'align/mergeSample/'), '.mSm.rmD.sorted.bam.flagstat'),
                    ('SAMPLE-LEVEL', 'macs2', '', os.path.join(ResultsDir,'align/mergeSample/'), '_peaks.broadPeak'),
                    ('SAMPLE-LEVEL', 'frip', '', os.path.join(ResultsDir,'align/mergeSample/'), '_peaks.frip.txt')]

    headerDict = {}
    qcDict = {}
    for section,tool,header_prefix,search_dir,extension in fileInfoList:
        fileList = funcs.recursive_glob(search_dir, '*%s' % (extension))

        if not qcDict.has_key(section):
            qcDict[section] = {}
        if not headerDict.has_key(section):
            headerDict[section] = []

        for idx in range(len(fileList)):
            sample = os.path.basename(fileList[idx])[:-len(extension)]
            if not qcDict[section].has_key(sample):
                qcDict[section][sample] = []

            if tool == 'cutadapt':
                fields = ['totalPairs','passTrimmedPairs','passTrimmedBases']
                ofields = fields
                cutadaptDict = funcs.cutadaptPELogToDict(fileList[idx])
                qcDict[section][sample] += [cutadaptDict[x] for x in fields]
                if idx == 0:
                    headerDict[section] += [header_prefix+' '+x for x in ofields]

            elif tool == 'flagstat':
                fields = ['mapped','properly paired','duplicates','read1','read2']
                ofields = ['mapped','properlyPaired','duplicates','read1','read2']
                flagstatDict = funcs.flagstatToDict(fileList[idx])
                qcDict[section][sample] += [str(flagstatDict[x][0]) for x in fields]
                if idx == 0:
                    headerDict[section] += [header_prefix+' '+x for x in ofields]

            elif tool == 'idxstats':
                fields = [MitoName]
                ofields = fields
                idxstatsDict = funcs.idxstatsToDict(fileList[idx])
                qcDict[section][sample] += [str(idxstatsDict[x][1]) for x in fields]
                if idx == 0:
                    headerDict[section] += [header_prefix+' '+x for x in ofields]

            elif tool == 'picard_insert_metrics':
                fields = ['MEAN_INSERT_SIZE', 'STANDARD_DEVIATION', 'MAX_INSERT_SIZE']
                ofields = ['insertMean', 'insertStdDev', 'insertMax']
                metricsDict = funcs.picardInsertMetricsToDict(fileList[idx])
                qcDict[section][sample] += [metricsDict[x] for x in fields]
                if idx == 0:
                    headerDict[section] += [header_prefix+' '+x for x in ofields]

            elif tool == 'macs2':
                fields = ['numPeaks']
                ofields = fields
                numLines = str(funcs.numLinesInFile(fileList[idx]))
                qcDict[section][sample] += [numLines]
                if idx == 0:
                    headerDict[section] += [header_prefix+' '+x for x in ofields]

            elif tool == 'frip':
                fields = ['fripScore']
                ofields = fields
                fin = open(fileList[idx],'r')
                frip = fin.readline().strip()
                fin.close()
                qcDict[section][sample] += [frip]
                if idx == 0:
                    headerDict[section] += [header_prefix+' '+x for x in ofields]

    sectionOrder = ['RUN-LEVEL', 'REPLICATE-LEVEL', 'SAMPLE-LEVEL']
    fout = open(OutFile,'w')
    for section in sectionOrder:
        fout.write('## %s\n' % (section))
        fout.write('\t'.join(['sample'] + headerDict[section]) + '\n')
        for sample in sorted(qcDict[section].keys()):
            fout.write('\t'.join([sample] + qcDict[section][sample]) + '\n')
        fout.write('\n')
    fout.close()

############################################
############################################
## RUN FUNCTION
############################################
############################################

pipeline_qc_to_tsv(ResultsDir=args.RESULTS_DIR,MitoName=args.MITO_NAME,OutFile=args.OUT_FILE)

############################################
############################################
############################################
############################################