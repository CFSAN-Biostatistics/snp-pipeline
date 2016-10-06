#!/usr/bin/env python

import argparse
import sys
from snppipeline import __version__
from snppipeline import snppipeline
from snppipeline import utils

#==============================================================================
# Command line driver
#==============================================================================


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Remove abnormally dense SNPs from the input VCF file, save the reserved SNPs into a new VCF file,
                                                    and save the removed SNPs into another VCF file.""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(                   dest='sampleDirsFile', type=str,                        help='Relative or absolute path to file containing a list of directories -- one per sample')
    parser.add_argument(                   dest='refFastaFile', type=str,                        help='Relative or absolute path to the reference fasta file')
    parser.add_argument('-f', '--force',   dest='forceFlag',      action='store_true',             help='Force processing even when result files already exist and are newer than inputs')
    parser.add_argument('-n', '--vcfname', dest='vcfFileName',    type=str, default='var.flt.vcf', metavar='NAME', help='File name of the input VCF files which must exist in each of the sample directories')
    parser.add_argument('-l', '--edge_length',  dest='edgeLength',    type=int, default=500, metavar='EDGE_LENGTH', help='The length of the edge regions in a contig, in which all SNPs will be removed.')
    parser.add_argument('-w', '--window_size',  dest='windowSize',  type=int, default=1000, metavar='WINDOW_SIZE', help='the length of the window in which the number of SNPs should be no more than max_num_snp.')
    parser.add_argument('-m', '--max_snp',  dest='maxSNP',  type=int, default=3, metavar='MAX_NUM_SNPs', help='The maximum number of SNPs allowed in a window.')
    parser.add_argument('-g', '--out_group',  dest='outGroupFile',  type=str, default=None, metavar='OUT_GROUP', help='Relative or absolute path to the file indicating outgroup samples, one sample ID per line.')
    parser.add_argument('-v', '--verbose', dest='verbose',        type=int, default=1,             metavar='0..5', help='Verbose message level (0=no info, 5=lots)')
    parser.add_argument('--version', action='version', version='%(prog)s version ' + __version__)

    args_dict=vars(parser.parse_args())

    sys.excepthook = utils.handle_global_exception
    utils.set_logging_verbosity(args_dict)
    snppipeline.set_logging_verbosity(args_dict)


    #==========================================================================
    # Validate arguments
    #==========================================================================
    if (args_dict["edgeLength"] <1):
        utils.global_error("Error: the length of the edge regions must be a positive integer, and the input is %d." % args_dict["edgeLength"])
        sys.exit("Error: the length of the edge regions must be a positive integer, and the input is %d." % args_dict["edgeLength"])

    if (args_dict["windowSize"] <1):
        utils.global_error("Error: the length of the window must be a positive integer, and the input is %d." % args_dict["windowSize"])
        sys.exit("Error: the length of the window must be a positive integer, and the input is %d." % args_dict["windowSize"])

    if (args_dict["maxSNP"] <1):
        utils.global_error("Error: the maximum number of SNPs allowed must be a positive integer, and the input is %d." % args_dict["maxSNP"])
        sys.exit("Error: the maximum number of SNPs allowed must be a positive integer, and the input is %d." % args_dict["maxSNP"])

    snppipeline.remove_bad_snp(args_dict)
