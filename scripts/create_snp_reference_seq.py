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

    parser = argparse.ArgumentParser(description='Write reference sequence bases at SNP locations to a fasta file.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(                          dest='referenceFile', type=str,                                               help='Relative or absolute path to the reference bases file in fasta format')
    parser.add_argument('-f', '--force',          dest='forceFlag',     action='store_true',                                    help='Force processing even when result file already exists and is newer than inputs')
    parser.add_argument('-l', '--snpListFile',    dest='snpListFile',   type=str, default='snplist.txt',        metavar='FILE', help='Relative or absolute path to the SNP list file')
    parser.add_argument('-o', '--output',         dest='snpRefFile',    type=str, default='referenceSNP.fasta', metavar='FILE', help='Output file.  Relative or absolute path to the SNP reference sequence file')
    parser.add_argument('-v', '--verbose',        dest='verbose',       type=int, default=1,                    metavar='0..5', help='Verbose message level (0=no info, 5=lots)')
    parser.add_argument('--version', action='version', version='%(prog)s version ' + __version__)
    args_dict = vars(parser.parse_args())

    sys.excepthook = utils.handle_global_exception
    utils.set_logging_verbosity(args_dict)
    snppipeline.set_logging_verbosity(args_dict)
    snppipeline.create_snp_reference_seq(args_dict)
