#######################################################################
#######################################################################
## Created on July 6th 2018 to store utility functions for .py exe
#######################################################################
#######################################################################

import os
import fnmatch
# from decimal import *

############################################
############################################
## UTILITY FUNCTIONS
############################################
############################################

def makedir(path):

    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

############################################

def recursive_glob(treeroot, pattern):

    results = []
    for base, dirs, files in os.walk(os.path.abspath(treeroot)):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)

    return results

############################################

# def floatToRoundDPStr(Float,DP=2):
#
#     return str(Decimal(str(Float)).quantize(Decimal(str(1/(float(10**DP)))),rounding='ROUND_HALF_EVEN'))

############################################
############################################
############################################
############################################
