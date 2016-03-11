#!/usr/bin/env python

import argparse
import sys
from snppipeline import __version__
from snppipeline import snppipeline
from snppipeline import utils

#==============================================================================
# Command line driver
#==============================================================================
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Combine the SNP positions across all samples into a single unified SNP
                                                    list file identifing the postions and sample names where SNPs were called.""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(                   dest='sampleDirsFile', type=str,                        help='Relative or absolute path to file containing a list of directories -- one per sample')
    parser.add_argument('-f', '--force',   dest='forceFlag',      action='store_true',             help='Force processing even when result file already exists and is newer than inputs')
    parser.add_argument('-n', '--vcfname', dest='vcfFileName',    type=str, default='var.flt.vcf', metavar='NAME', help='File name of the VCF files which must exist in each of the sample directories')
    parser.add_argument(      '--maxsnps', dest='maxSnps',        type=int, default=-1,            metavar='INT',  help='Exclude samples having more than this maximum allowed number of SNPs. Set to -1 to disable this function.')
    parser.add_argument('-o', '--output',  dest='snpListFile',    type=str, default='snplist.txt', metavar='FILE', help='Output file.  Relative or absolute path to the SNP list file')
    parser.add_argument('-v', '--verbose', dest='verbose',        type=int, default=1,             metavar='0..5', help='Verbose message level (0=no info, 5=lots)')
    parser.add_argument('--version', action='version', version='%(prog)s version ' + __version__)
    args_dict = vars(parser.parse_args())

    sys.excepthook = utils.handle_global_exception
    utils.set_logging_verbosity(args_dict)
    snppipeline.set_logging_verbosity(args_dict)
    snppipeline.create_snp_list(args_dict)
