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

    def minConsFreq(value):
        fvalue = float(value)
        if fvalue <= 0.5 or fvalue > 1:
            raise argparse.ArgumentTypeError("Minimum consensus frequency must be > 0.5 and <= 1.0")
        return fvalue

    def minConsStrdBias(value):
        fvalue = float(value)
        if fvalue < 0.0 or fvalue > 0.5:
            raise argparse.ArgumentTypeError("Minimum consensus strand bias must be >= 0.0 and <= 0.5")
        return fvalue

    parser = argparse.ArgumentParser(description="""Call the consensus base for a sample at the specified positions
                                                    where SNPs were previously called in any of the samples.  Generates
                                                    a single-sequence fasta file with one base per specified position.""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    help = dict()
    help['allPileupFile']  = """Relative or absolute path to the genome-wide pileup file for this sample."""
    help['force']          = """Force processing even when result file already exists and is newer than inputs."""
    help['snpListFile']    = """Relative or absolute path to the SNP list file across all samples."""
    help['output']         = """Output file. Relative or absolute path to the consensus fasta file for this sample."""
    help['minBaseQual']    = """Mimimum base quality score to count a read. All other snp filters take effect after the low-quality reads are discarded."""
    help['minConsFreq']    = """Consensus frequency. Mimimum fraction of high-quality reads supporting the consensus to make a call."""
    help['minConsStrdDpth']= """Consensus strand depth. Minimum number of high-quality reads supporting the consensus which must be present on both the
                                forward and reverse strands to make a call."""
    help['minConsStrdBias']= """Strand bias. Minimum fraction of the high-quality consensus-supporting reads which must be present on both the
                                forward and reverse strands to make a call. The numerator of this fraction is the number of high-quality
                                consensus-supporting reads on one strand.  The denominator is the total number of high-quality consensus-supporting
                                reads on both strands combined."""
    help['vcfFileName']    = """VCF Output file name. If specified, a VCF file with this file name will be created in the same directory as the
                                consensus fasta file for this sample."""
    help['vcfRefName']     = """Name of the reference file.  This is only used in the generated VCF file header."""
    help['vcfAllPos']      = """Flag to cause VCF file generation at all positions, not just the snp positions.  This has no effect on
                                the consensus fasta file, it only affects the VCF file.  This capability is intended primarily as a diagnostic tool and
                                enabling this flag will greatly increase execution time."""
    help['vcfPreserveRefCase'] = """Flag to cause the VCF file generator to emit each reference base in uppercase/lowercase as it appears in the reference
                                    sequence file.  If not specified, the reference base is emitted in uppercase."""
    help['vcfFailedSnpGt'] = """Controls the VCF file GT data element when a snp fails filters.  Possible values:
                                .) The GT element will be a dot, indicating unable to make a call.
                                0) The GT element will be 0, indicating the reference base.
                                1) The GT element will be the ALT index of the most commonly occuring base, usually 1."""

    help['verbose']        = """Verbose message level (0=no info, 5=lots)"""

    parser.add_argument(                              dest='allPileupFile',      type=str,                                                        help=help['allPileupFile'])
    parser.add_argument('-f', '--force',              dest='forceFlag',          action='store_true',                                             help=help['force'])
    parser.add_argument('-l', '--snpListFile',        dest='snpListFile',        type=str,            default='snplist.txt',      metavar='FILE', help=help['snpListFile'])
    parser.add_argument('-o', '--output',             dest='consensusFile',      type=str,            default='consensus.fasta',  metavar='FILE', help=help['output'])
    parser.add_argument('-q', '--minBaseQual',        dest='minBaseQual',        type=int,            default=0,                  metavar='INT',  help=help['minBaseQual'])
    parser.add_argument('-c', '--minConsFreq',        dest='minConsFreq',        type=minConsFreq,    default=0.60,               metavar='FREQ', help=help['minConsFreq'])
    parser.add_argument('-d', '--minConsStrdDpth',    dest='minConsStrdDpth',    type=int,            default=0,                  metavar='INT',  help=help['minConsStrdDpth'])
    parser.add_argument('-b', '--minConsStrdBias',    dest='minConsStrdBias',    type=minConsStrdBias,default=0,                  metavar='FREQ', help=help['minConsStrdBias'])
    parser.add_argument(      '--vcfFileName',        dest='vcfFileName',        type=str,            default=None,               metavar='NAME', help=help['vcfFileName'])
    parser.add_argument(      '--vcfRefName',         dest='vcfRefName',         type=str,            default='Unknown reference',metavar='NAME', help=help['vcfRefName'])
    parser.add_argument(      '--vcfAllPos',          dest='vcfAllPos',          action='store_true',                                             help=help['vcfAllPos'])
    parser.add_argument(      '--vcfPreserveRefCase', dest='vcfPreserveRefCase', action='store_true',                                             help=help['vcfPreserveRefCase'])
    parser.add_argument(      '--vcfFailedSnpGt',     dest='vcfFailedSnpGt',     type=str,            default=".",                choices=['.','0','1'], help=help['vcfFailedSnpGt'], )
    parser.add_argument('-v', '--verbose',            dest='verbose',            type=int,            default=1,                  metavar='0..5', help=help['verbose'])
    parser.add_argument('--version', action='version', version='%(prog)s version ' + __version__)
    args_dict = vars(parser.parse_args())

    sys.excepthook = utils.handle_sample_exception
    utils.set_logging_verbosity(args_dict)
    snppipeline.set_logging_verbosity(args_dict)
    snppipeline.call_consensus(args_dict)
