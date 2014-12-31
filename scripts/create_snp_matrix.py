#!/usr/bin/env python2.7

import argparse
from snppipeline import __version__
from snppipeline import snppipeline
from snppipeline import utils

#==============================================================================
# Command line driver
#==============================================================================
if __name__ == '__main__':

    def minConsFreq(value):
        fvalue = float(value)
        if fvalue <= 0.5 or fvalue > 1:
            raise argparse.ArgumentTypeError("Consensus threshold must be > 0.5 and <= 1.0")
        return fvalue

    parser = argparse.ArgumentParser(description="""Create the SNP matrix containing the consensus base for each of the samples 
                                                    at the positions where SNPs were called in any of the samples.  The matrix 
                                                    contains one row per sample and one column per SNP position.  Non-SNP 
                                                    positions are not included in the matrix.  The matrix is formatted as a fasta
                                                    file, with each sequence (all of identical length) corresponding to the SNPs 
                                                    in the correspondingly named sequence.""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(                          dest='sampleDirsFile',     type=str,                                               help='Relative or absolute path to file containing a list of directories -- one per sample')
    parser.add_argument('-f', '--force',          dest='forceFlag',          action='store_true',                                    help='Force processing even when result file already exists and is newer than inputs')
    parser.add_argument('-l', '--snpListFile',    dest='snpListFile',        type=str,   default='snplist.txt',      metavar='FILE', help='Relative or absolute path to the SNP list file')
    parser.add_argument('-p', '--pileupFileName', dest='pileupFileName',     type=str,   default='reads.snp.pileup', metavar='NAME', help='File name of the SNP pileup files which must exist in each of the sample directories')
    parser.add_argument('-o', '--output',         dest='snpmaFile',          type=str,   default='snpma.fasta',      metavar='FILE', help='Output file.  Relative or absolute path to the SNP matrix file')
    parser.add_argument('-c', '--minConsFreq',    dest='minConsFreq',  type=minConsFreq, default=0.60,               metavar='FREQ', help='Mimimum fraction of reads that must agree to make a consensus call')
    parser.add_argument('-v', '--verbose',        dest='verbose',            type=int,   default=1,                  metavar='0..5', help='Verbose message level (0=no info, 5=lots)')
    parser.add_argument('--version', action='version', version='%(prog)s version ' + __version__)
    args_dict = vars(parser.parse_args())

    utils.set_logging_verbosity(args_dict)
    snppipeline.set_logging_verbosity(args_dict)
    snppipeline.create_snp_matrix(args_dict)
