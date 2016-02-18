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

    parser = argparse.ArgumentParser(description="""Create the SNP matrix containing the consensus base for each of the samples
                                                    at the positions where SNPs were called in any of the samples.  The matrix
                                                    contains one row per sample and one column per SNP position.  Non-SNP
                                                    positions are not included in the matrix.  The matrix is formatted as a fasta
                                                    file, with each sequence (all of identical length) corresponding to the SNPs
                                                    in the correspondingly named sequence.""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(                          dest='sampleDirsFile',     type=str,                                               help='Relative or absolute path to file containing a list of directories -- one per sample')
    parser.add_argument('-f', '--force',          dest='forceFlag',          action='store_true',                                    help='Force processing even when result file already exists and is newer than inputs')
    parser.add_argument('-c', '--consFileName',   dest='consFileName',       type=str,   default='consensus.fasta',  metavar='NAME', help='File name of the previously created consensus SNP call file which must exist in each of the sample directories')
    parser.add_argument('-o', '--output',         dest='snpmaFile',          type=str,   default='snpma.fasta',      metavar='FILE', help='Output file.  Relative or absolute path to the SNP matrix file')
    parser.add_argument('-v', '--verbose',        dest='verbose',            type=int,   default=1,                  metavar='0..5', help='Verbose message level (0=no info, 5=lots)')
    parser.add_argument('--version', action='version', version='%(prog)s version ' + __version__)
    args_dict = vars(parser.parse_args())

    sys.excepthook = utils.handle_global_exception
    utils.set_logging_verbosity(args_dict)
    snppipeline.set_logging_verbosity(args_dict)
    snppipeline.create_snp_matrix(args_dict)
