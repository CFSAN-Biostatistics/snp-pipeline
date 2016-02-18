#!/usr/bin/env python

import argparse
from snppipeline import __version__
from snppipeline import snppipeline
from snppipeline import utils

#==============================================================================
# Command line driver
#==============================================================================
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Create the SNP pileup file for a sample -- the pileup file
                                                    at the positions where SNPs were called in any of the samples.""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--force',         dest='forceFlag',     action='store_true',                                  help='Force processing even when result file already exists and is newer than inputs')
    parser.add_argument('-l', '--snpListFile',   dest='snpListFile',   type=str, default='snplist.txt',      metavar='FILE', help='Relative or absolute path to the SNP list file across all samples')
    parser.add_argument('-a', '--allPileupFile', dest='allPileupFile', type=str, default='reads.all.pileup', metavar='FILE', help='Relative or absolute path to the genome-wide pileup file for this sample')
    parser.add_argument('-o', '--output',        dest='snpPileupFile', type=str, default='reads.snp.pileup', metavar='FILE', help='Output file. Relative or absolute path to the sample SNP pileup file')
    parser.add_argument('-v', '--verbose',       dest='verbose',       type=int, default=1,                  metavar='0..5', help='Verbose message level (0=no info, 5=lots)')
    parser.add_argument('--version', action='version', version='%(prog)s version ' + __version__)
    args_dict = vars(parser.parse_args())

    utils.set_logging_verbosity(args_dict)
    snppipeline.set_logging_verbosity(args_dict)
    snppipeline.create_snp_pileup(args_dict)
