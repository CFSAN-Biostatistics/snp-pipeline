#!/usr/bin/env python2.7

import argparse
from snppipeline import snppipeline

#==============================================================================
# Command line driver
#==============================================================================
if __name__ == '__main__':

#TODO Change -i and -o aruguments so that they use action='store_true' and act as just flags.

    parser = argparse.ArgumentParser(description='Run SNP pipeline.')
    parser.add_argument('-n', '--n-processes',      dest='maxThread',                        type=int,   default=4,                 help='Max number of concurrent processes launched')
    parser.add_argument('-d', '--mainPath',         dest='mainPath',                         type=str,   default='',                help='Directory containing the reference subdirectory, path file name, snp list, snp matrix')
    parser.add_argument('-r', '--Reference',        dest='Reference',                        type=str,   default='reference.fasta', help='Reference file for mapping')
    parser.add_argument('-f', '--pathFileName',     dest='pathFileName',                     type=str,   default='path.txt',        help='File containing a list of directories -- one per sample')
    parser.add_argument('-l', '--snplistFileName',  dest='snplistFileName',                  type=str,   default='snplist.txt',     help='Output Snplist file name')
    parser.add_argument('-a', '--snpmaFileName',    dest='snpmaFileName',                    type=str,   default='snpma.fa',        help='Output fasta file name')
    parser.add_argument('-b', '--bamFileName',      dest='bamFileName',                      type=str,   default='reads.bam',       help='bam file name')
    parser.add_argument('-p', '--pileupFileName',   dest='pileupFileName',                   type=str,   default='reads.pileup',    help='pileup file name')
    parser.add_argument('-v', '--verbose',          dest='verbose',                          type=int,   default=1,                 help='Verbose flag (0=no info, 5=lots)')
    parser.add_argument('-i', '--includeReference', dest='includeReference',                 type=bool,  default=False,             help='Write reference sequence bases at SNP positions in fasta format.')
    parser.add_argument('-o', '--useOldPileups',    dest='useOldPileups',                    type=bool,  default=False,             help='Use available pileup files.')
    args_dict = vars(parser.parse_args())

    snppipeline.run_snp_pipeline(args_dict)
