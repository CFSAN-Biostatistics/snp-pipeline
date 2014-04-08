#!/usr/bin/env python2.7

import argparse
import imp
#TODO fix these two lines to make paths relative
snppipeline = imp.load_source('snppipeline', '/home/hugh.rand/projects/snppipeline/snppipeline/snppipeline.py')

#==============================================================================
# Command line driver
#==============================================================================
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run SNP pipeline.')
    parser.add_argument('-n', '--n-processes',      dest='maxThread',        type=int,  default=4,                 help='Max number of concurrent jobs.')
    parser.add_argument('-d', '--mainPath',         dest='mainPath',         type=str,  default='', help='Path for all files')
    parser.add_argument('-r', '--Reference',        dest='Reference',        type=str,  default='reference.fasta', help='reference for mapping')
    parser.add_argument('-f', '--pathFileName',     dest='pathFileName',     type=str,  default='path.txt',        help='Path file name')
    parser.add_argument('-l', '--snplistFileName',  dest='snplistFileName',  type=str,  default='snplist.txt',     help='Snplist file name')
    parser.add_argument('-a', '--snpmaFileName',    dest='snpmaFileName',    type=str,  default='snpma.fa',        help='fasta file name')
    parser.add_argument('-b', '--bamFileName',      dest='bamFileName',      type=str,  default='reads.bam',       help='bam file name')
    parser.add_argument('-p', '--pileupFileName',   dest='pileupFileName',   type=str,  default='reads.pileup',    help='pileup file name')
    parser.add_argument('-v', '--verbose',          dest='verbose',          type=int,  default=1,                 help='Verbose flag (0=no info, 5=lots')
    parser.add_argument('-i', '--includeReference', dest='includeReference', type=bool, default=False,             help='Write reference sequence bases at SNP positions in fasta format.')
    parser.add_argument('-o', '--useOldPileups',    dest='useOldPileups',    type=bool, default=False,             help='Use available pileup files.')
    parser.add_argument(      '--DP',               dest='combinedDepthAcrossSamples',       type=int,   default=10,  help='Combined depth across samples.')
    parser.add_argument(      '--AF1',              dest='alleleFrequencyForFirstALTAllele', type=float, default=1.0, help='Allele frequency for first allele.')
    parser.add_argument(      '--AR',               dest='arFlagValue',                      type=float, default=1.0, help='AR flag value.')
    args_dict = vars(parser.parse_args())

    snppipeline.run_snp_pipeline(args_dict)

