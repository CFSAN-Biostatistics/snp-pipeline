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

    parser = argparse.ArgumentParser(description="""Calculate pairwise SNP distances from the multi-fasta SNP matrix. Generates a file of pairwise
                                                    distances and a file containing a matrix of distances.""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(                  dest='inputFile',    type=str,                           metavar='snpMatrixFile', help='Relative or absolute path to the input multi-fasta SNP matrix file.')
    parser.add_argument('-f', '--force',  dest='forceFlag',    action='store_true',                                         help='Force processing even when result file already exists and is newer than inputs')
    parser.add_argument('-p', '--pairs',  dest='pairwiseFile', type=str, default=None,             metavar='FILE',          help='Relative or absolute path to the pairwise distance output file.')
    parser.add_argument('-m', '--matrix', dest='matrixFile',   type=str, default=None,             metavar='FILE',          help='Relative or absolute path to the distance matrix output file.')
    parser.add_argument('-v', '--verbose', dest='verbose',     type=int, default=1,                metavar='0..5',          help='Verbose message level (0=no info, 5=lots)')
    parser.add_argument('--version', action='version', version='%(prog)s version ' + __version__)
    args_dict = vars(parser.parse_args())

    sys.excepthook = utils.handle_global_exception
    utils.set_logging_verbosity(args_dict)
    snppipeline.set_logging_verbosity(args_dict)
    snppipeline.calculate_snp_distances(args_dict)
