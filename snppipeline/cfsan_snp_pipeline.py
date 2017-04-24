#!/usr/bin/env python

"""
This is the top-level script for the collection of cfsan_snp_pipeline tools.
"""

from __future__ import absolute_import
import argparse
import sys
from snppipeline import __version__
from snppipeline import utils
from snppipeline.utils import verbose_print

from snppipeline import data
from snppipeline import index_ref
from snppipeline import map_reads
from snppipeline import call_sites
from snppipeline import merge_vcfs
from snppipeline import collect_metrics
from snppipeline import combine_metrics

#==============================================================================
# Command line driver
#==============================================================================

def not_implemented(args):
    print("The %s command will be added in a future release" % args.subparser_name)


def parse_arguments(argv):
    """
    Parse command line arguments.

    Parameters
    ----------
    argv : list
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
    # Create the parser for the "data" command
    # -------------------------------------------------------------------------
    description = """Copy data included with the CFSAN SNP Pipeline to a specified directory."""
    subparser = subparsers.add_parser("data", help="Copy included data to a specified directory", description=description, formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Example:
# create a new directory "testLambdaVirus" and copy the Lambda virus input data there
$ cfsan_snp_pipeline data lambdaVirusInputs testLambdaVirus
"""
    )

    subparser.add_argument("whichData",
    metavar="whichData",
    choices=["lambdaVirusInputs", "lambdaVirusExpectedResults",
             "agonaInputs", "agonaExpectedResults",
             "listeriaInputs", "listeriaExpectedResults",
             "configurationFile"],
    help="""    Which of the supplied data sets to copy.  The choices are:
        lambdaVirusInputs          : Input reference and fastq files
        lambdaVirusExpectedResults : Expected results files
        agonaInputs                : Input reference file
        agonaExpectedResults       : Expected results files
        listeriaInputs             : Input reference file
        listeriaExpectedResults    : Expected results files
        configurationFile          : File of parameters to customize the
                                     SNP pipeline

    Note: the lambda virus data set is complete with input data and expected
    results.  The agona and listeria data sets have the reference genome and
    the expected results, but not the input fastq files, because the files are
    too large to include with the package.
    """)

    subparser.add_argument("destDirectory",
    nargs="?",
    type=str,
    default=".",
    help="""    Destination directory into which the SNP pipeline data files will be copied.
    The data files are copied into the destination directory if the directory
    already exists.  Otherwise the destination directory is created and the
    data files are copied there.  (default: current directory)""")
    subparser.add_argument("--version", action="version", version="%(prog)s version " + __version__)
    subparser.set_defaults(func=data.copy_data)
    subparser.set_defaults(excepthook=None) # keep default exception handler
    subparser.set_defaults(verbose=0)

    # -------------------------------------------------------------------------
    # Create the parser for the "index_ref" command
    # -------------------------------------------------------------------------
    description = """Index the reference genome for subsequent read mapping, and create
                     the faidx index file for subsequent pileups. The output is written
                     to the reference directory."""
    subparser = subparsers.add_parser("index_ref", help="Index the reference", description=description, formatter_class=formatter_class)
    subparser.add_argument(dest="referenceFile",    type=str, help="Relative or absolute path to the reference fasta file")
    subparser.add_argument("-f", "--force",   dest="forceFlag", action="store_true", help="Force processing even when result files already exist and are newer than inputs")
    subparser.add_argument("-v", "--verbose", dest="verbose",   type=int, default=1, metavar="0..5", help="Verbose message level (0=no info, 5=lots)")
    subparser.add_argument("--version", action="version", version="%(prog)s version " + __version__)
    subparser.set_defaults(func=index_ref.index_ref)
    subparser.set_defaults(excepthook=utils.handle_global_exception)

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
    # Create the parser for the "call_sites" command
    # -------------------------------------------------------------------------
    description = """Find the sites with SNPs in a sample."""
    subparser = subparsers.add_parser("call_sites", help="Find the sites with SNPs in a sample", description=description, formatter_class=formatter_class)
    subparser.add_argument(dest="referenceFile",    type=str, help="Relative or absolute path to the reference fasta file")
    subparser.add_argument(dest="sampleDir", type=str, help="Relative or absolute directory of the sample")
    subparser.add_argument("-f", "--force",   dest="forceFlag", action="store_true", help="Force processing even when result files already exist and are newer than inputs")
    subparser.add_argument("-v", "--verbose", dest="verbose",   type=int, default=1, metavar="0..5", help="Verbose message level (0=no info, 5=lots)")
    subparser.add_argument("--version", action="version", version="%(prog)s version " + __version__)
    subparser.set_defaults(func=call_sites.call_sites)
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
    # Create the parser for the "collect_metrics" command
    # -------------------------------------------------------------------------
    description = """Collect alignment, coverage, and variant metrics for a single specified sample."""
    subparser = subparsers.add_parser("collect_metrics", help="Collect quality and SNP metrics for a sample", description=description, formatter_class=formatter_class)
    subparser.add_argument(dest="sampleDir", type=str, help="Relative or absolute directory of the sample")
    subparser.add_argument(dest="referenceFile",    type=str, help="Relative or absolute path to the reference fasta file")
    subparser.add_argument("-f", "--force",   dest="forceFlag", action="store_true", help="Force processing even when result files already exist and are newer than inputs")
    subparser.add_argument("-o", "--output",  dest="metricsFile",            type=str, default="metrics",                   metavar="FILE", help="Output file.  Relative or absolute path to the metrics file")
    subparser.add_argument("-m", "--maxsnps", dest="maxSnps",                type=int, default=-1,                          metavar="INT",  help="Maximum allowed number of SNPs per sample")
    subparser.add_argument("-c", dest="consensusFastaFileName",              type=str, default="consensus.fasta",           metavar="NAME", help="File name of the consensus fasta file which must exist in the sample directory")
    subparser.add_argument("-C", dest="consensusPreservedFastaFileName",     type=str, default="consensus_preserved.fasta", metavar="NAME", help="File name of the consensus preserved fasta file which must exist in the sample directory")
    subparser.add_argument("-v", dest="consensusVcfFileName",                type=str, default="consensus.vcf",             metavar="NAME", help="File name of the consensus vcf file which must exist in the sample directory")
    subparser.add_argument("-V", dest="consensusPreservedVcfFileName",       type=str, default="consensus_preserved.vcf",   metavar="NAME", help="File name of the consensus preserved vcf file which must exist in the sample directory")
    subparser.add_argument("--verbose", dest="verbose",                      type=int, default=1,                           metavar="0..5", help="Verbose message level (0=no info, 5=lots)")
    subparser.add_argument("--version", action="version", version="%(prog)s version " + __version__)
    subparser.set_defaults(func=collect_metrics.collect_metrics)
    subparser.set_defaults(excepthook=utils.handle_sample_exception)

    # -------------------------------------------------------------------------
    # Create the parser for the "combine_metrics" command
    # -------------------------------------------------------------------------
    description = """Combine the metrics from all samples into a single table of metrics for all samples.
                     The output is a tab-separated-values file with a row for each sample and a column
                     for each metric.

                     Before running this command, the metrics for each sample must be created by the
                     collectSampleMetrics.sh script."""

    subparser = subparsers.add_parser("combine_metrics", help="Merge the per-sample metrics", description=description, formatter_class=formatter_class)
    subparser.add_argument(dest="sampleDirsFile", type=str, help="Relative or absolute path to file containing a list of directories -- one per sample")
    subparser.add_argument("-f", "--force",   dest="forceFlag", action="store_true", help="Force processing even when result files already exist and are newer than inputs")
    subparser.add_argument("-n", "--metrics", dest="metricsFileName", type=str, default="metrics", metavar="NAME", help="File name of the metrics files which must exist in each of the sample directories.")
    subparser.add_argument("-o", "--output",  dest="mergedMetricsFile", type=str, default="metrics.tsv", metavar="FILE", help="Output file. Relative or absolute path to the combined metrics file.")
    subparser.add_argument("-s", "--spaces",  dest="spaceHeadings", action="store_true", help="Emit column headings with spaces instead of underscores")
    subparser.add_argument("-v", "--verbose", dest="verbose",   type=int, default=1, metavar="0..5", help="Verbose message level (0=no info, 5=lots)")
    subparser.add_argument("--version", action="version", version="%(prog)s version " + __version__)
    subparser.set_defaults(func=combine_metrics.combine_metrics)
    subparser.set_defaults(excepthook=utils.handle_global_exception)

    # -------------------------------------------------------------------------
    # parse the args
    # -------------------------------------------------------------------------
    args = parser.parse_args(argv)

    # Special validation
    if args.subparser_name == "my-subcommand":
        pass  # special validation here

    return args


def run_command(argv):
    """Run a subcommand with arguments in the argv list.

    Parameters
    ----------
    argv : list of str
        List of command line arguments.  Usually sys.argv[1:].
        The first item in the list should be the subcommand name.

    Returns
    -------
    Returns 0 on success if it completes with no exceptions.
    """
    args = parse_arguments(argv)

    # Call the sub-command function
    if args.excepthook:
        sys.excepthook = args.excepthook
    utils.set_logging_verbosity(args)
    args.func(args)  # this executes the function previously associated with the subparser with set_defaults

    verbose_print("")
    verbose_print("# %s %s %s finished" % (utils.timestamp(), utils.program_name(), args.subparser_name))
    return 0


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
    return run_command(sys.argv[1:])
