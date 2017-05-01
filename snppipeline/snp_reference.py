"""This module is part of the CFSAN SNP Pipeline. It contains the code to
write reference sequence bases at SNP locations to a fasta file.
"""

from __future__ import print_function
from __future__ import absolute_import

from snppipeline import utils
from snppipeline.utils import verbose_print


def create_snp_reference_seq(args):
    """Write reference sequence bases at SNP locations to a fasta file.

    Write reference sequence bases at SNP locations to a fasta file.
    This function expects, or creates '(*)', the following files:
            reference.fasta
            snplist.txt
            referenceSNP.fasta (*)

    The files are used as follows:
        1. The reference.fasta input file contains the whole-genome reference
           bases.
        2. The snplist.txt input file contains the list of SNP positions across
           all the samples.
        2. The referenceSNP.fasta output file contains the reference bases at
           the identified SNP locations.

    The snplist.txt file is created outside of this function.  The package
        documentation provides an example of creating this file based on the
        lambda_virus sequence that is used as one test for this package.

    Parameters
    ----------
    args : Namespace
        referenceFile: File path (not just file name) for reference sequence in fasta format
        snpListFile: File path (not just file name) of text format list of SNP positions
        snpRefFile: File path (not just file name) for the SNP reference sequence file.

    Raises:

    Examples:
    args = argparse.Namespace
    args.referenceFile = 'reference.fasta'
    args.snpListFile = 'snplist.txt'
    args.snpRefFile = 'referenceSNP.fasta'
    create_snp_reference_seq(args)
    """
    utils.print_log_header()
    utils.print_arguments(args)

    #==========================================================================
    #    Write reference sequence bases at SNP locations to a fasta file.
    #==========================================================================
    reference_file = args.referenceFile
    snp_list_file_path = args.snpListFile
    snp_ref_seq_path = args.snpRefFile

    #==========================================================================
    # Verify input files exist
    #==========================================================================
    bad_file_count = utils.verify_existing_input_files("Snplist file", [snp_list_file_path])
    if bad_file_count > 0:
        utils.global_error("Error: cannot create the snp reference sequence without the snplist file.")

    bad_file_count = utils.verify_non_empty_input_files("Reference file", [reference_file])
    if bad_file_count > 0:
        utils.global_error("Error: cannot create the snp reference sequence without the reference fasta file.")

    #==========================================================================
    # Find the reference bases at the snp positions
    #==========================================================================
    source_files = [reference_file, snp_list_file_path]
    if args.forceFlag or utils.target_needs_rebuild(source_files, snp_ref_seq_path):
        utils.write_reference_snp_file(reference_file, snp_list_file_path, snp_ref_seq_path)
    else:
        verbose_print("SNP reference sequence %s has already been freshly built.  Use the -f option to force a rebuild." % snp_ref_seq_path)
