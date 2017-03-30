#!/usr/bin/env python

"""
This is the top-level script for the collection of cfsan_snp_pipeline tools.
"""

from __future__ import absolute_import
import argparse
import sys
from snppipeline import __version__
from snppipeline import snppipeline
from snppipeline import utils
from snppipeline.utils import verbose_print

from snppipeline import map_reads
from snppipeline import merge_vcfs

#==============================================================================
# Command line driver
#==============================================================================

def not_implemented(args):
    print("The %s command will be added in a future release" % args.subparser_name)


def parse_arguments(system_args):
    """
    Parse command line arguments.

    Parameters
    ----------
    system_args : list
        List of command line arguments, usually sys.argv[1:].

    Returns
    -------
    Namespace
        Command line arguments are stored as attributes of a Namespace.
    """
    # Create the top-level parser
    description = """The CFSAN SNP Pipeline is a collection of tools using reference-based
                     alignments to call SNPs for a set of samples."""
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--version", action="version", version="%(prog)s version " + __version__)
    subparsers = parser.add_subparsers(dest="subparser_name", help="Command help:")

    # -------------------------------------------------------------------------
    # Create the parser for the "map_reads" command
    # -------------------------------------------------------------------------
    description = """Align the sequence reads for a specified sample to a specified reference genome.
                     The output is written to the file "reads.sam" in the sample directory."""
    subparser = subparsers.add_parser("map_reads", help="Align reads to the reference", description=description, formatter_class=formatter_class)
    subparser.add_argument(dest="referenceFile",    type=str, help="Relative or absolute path to the reference fasta file")
    subparser.add_argument(dest="sampleFastqFile1", type=str, help="Relative or absolute path to the fastq file")
    subparser.add_argument(dest="sampleFastqFile2", type=str, help="Optional relative or absolute path to the mate fastq file, if paired", nargs="?")
    subparser.add_argument("-f", "--force",   dest="forceFlag", action="store_true", help="Force processing even when result files already exist and are newer than inputs")
    subparser.add_argument("-v", "--verbose", dest="verbose",   type=int, default=1, metavar="0..5", help="Verbose message level (0=no info, 5=lots)")
    subparser.add_argument("--version", action="version", version="%(prog)s version " + __version__)
    subparser.set_defaults(func=map_reads.map_reads)
    subparser.set_defaults(excepthook=utils.handle_sample_exception)

    # -------------------------------------------------------------------------
    # Create the parser for the "merge_vcfs" command
    # -------------------------------------------------------------------------
    description = """Merge the consensus vcf files from all samples into a single multi-vcf file for all samples."""
    subparser = subparsers.add_parser("merge_vcfs", help="Merge the per-sample VCF files", description=description, formatter_class=formatter_class)
    subparser.add_argument(dest="sampleDirsFile", type=str, help="Relative or absolute path to file containing a list of directories -- one per sample")
    subparser.add_argument("-f", "--force",   dest="forceFlag", action="store_true", help="Force processing even when result files already exist and are newer than inputs")
    subparser.add_argument("-n", "--vcfname", dest="vcfFileName",   type=str, default="consensus.vcf", metavar="NAME", help="File name of the vcf files which must exist in each of the sample directories")
    subparser.add_argument("-o", "--output",  dest="mergedVcfFile", type=str, default="snpma.vcf",     metavar="FILE", help="Output file.  Relative or absolute path to the merged multi-vcf file")
    subparser.add_argument("-v", "--verbose", dest="verbose",   type=int, default=1, metavar="0..5", help="Verbose message level (0=no info, 5=lots)")
    subparser.add_argument("--version", action="version", version="%(prog)s version " + __version__)
    subparser.set_defaults(func=merge_vcfs.merge_vcfs)
    subparser.set_defaults(excepthook=utils.handle_global_exception)

    # -------------------------------------------------------------------------
    # parse the args
    # -------------------------------------------------------------------------
    args = parser.parse_args()

    # Special validation
    if args.subparser_name == "my-subcommand":
        pass  # special validation here

    return args


def main():
    """This is the main function which is magically turned into an executable
    script by the setuptools entry_points.  See setup.py.

    To run this function as a script, first install the package:
        $ python setup.py develop
        or
        $ pip install --user snp-pipeline


    Parameters
    ----------
    This function must not take any parameters

    Returns
    -------
    The return value is passed to sys.exit()
    """
    args = parse_arguments(sys.argv[1:])

    # Call the sub-command function
    sys.excepthook = args.excepthook
    utils.set_logging_verbosity(args)
    args.func(args)  # this executes the function previously associated with the subparser with set_defaults

    verbose_print("")
    verbose_print("# %s %s %s finished" % (utils.timestamp(), utils.program_name(), args.subparser_name))
    return 0
