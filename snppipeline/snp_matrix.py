"""This module is part of the CFSAN SNP Pipeline. It contains the code to
create the matrix of snps.
"""

from __future__ import print_function
from __future__ import absolute_import

import os

from snppipeline import utils


def create_snp_matrix(args):
    """Create SNP matrix.

    Create the SNP matrix containing the consensus base for each of the samples
    at the positions where SNPs were found in any of the samples.  The matrix
    contains one row per sample and one column per SNP position.  Non-SNP
    positions are not included in the matrix.
    This function expects, or creates '(*)', the following
        files arranged in the following way:
            sampleDirectories.txt
            samples
                sample_name_one/consensus.fasta
                ...
            snpma.fasta (*)

    The files are used as follows:
        1. The sampleDirectories.txt input file contains a list of the paths to
           the sample directories.
        2. The consensus.fasta input files are previously called consensus
           for each sample to construct the SNP matrix fasta file.
        3. The snpma.fasta output file contains the SNP calls for each
           sequence, arranged as a multi-fasta file with one sequence per
           sample.

    The sampleDirectories.txt, and consensus.fasta are created outside of this
        function. The package documentation provides an example of creating
        these files based on the lambda_virus sequence that is used as one
        test for this package.

    Parameters
    ----------
    args : Namespace
        sampleDirsFile : str
            File path (not just file name) of file containing paths
            to directories containing consensus.fasta file for each sequence.
        snpListFile : str
            File path (not just file name) of text format list of SNP positions
        consFileName : str
            File name of the previously called consensus fasta files which must
            exist in each of the sample directories
        snpmaFile : str
            File path (not just file name) of the output snp matrix, formatted
            as a fasta file, with each sequence (all of identical length)
            corresponding to the SNPs in the correspondingly named sequence.

    Raises:

    Examples:
    args = argparse.Namespace
    args.sampleDirsFile = 'sampleDirectories.txt'
    args.consFileName = 'consensus.fasta'
    args.snpmaFile = 'snpma.fasta'
    args.minConsFreq = 0.6
    create_snp_matrix(args)
    """
    utils.print_log_header()
    utils.print_arguments(args)

    #==========================================================================
    # Prep work
    #==========================================================================
    sample_directories_list_filename = args.sampleDirsFile
    bad_file_count = utils.verify_non_empty_input_files("File of sample directories", [sample_directories_list_filename])
    if bad_file_count > 0:
        utils.global_error(None)

    with open(sample_directories_list_filename, "r") as sample_directories_list_file:
        list_of_sample_directories = [line.rstrip() for line in sample_directories_list_file]
    list_of_sample_directories = sorted([d for d in list_of_sample_directories if d])

    #==========================================================================
    # Verify input consensus.fasta files exist
    #==========================================================================
    consensus_files = []
    bad_file_count = 0
    for sample_directory in list_of_sample_directories:
        consensus_file_path = os.path.join(sample_directory, args.consFileName)
        bad_count = utils.verify_non_empty_input_files("Consensus fasta file", [consensus_file_path])
        if bad_count == 1:
            bad_file_count += 1
        else:
            consensus_files.append(consensus_file_path)  # keep the list of good files

    if bad_file_count == len(list_of_sample_directories):
        utils.global_error("Error: all %d consensus fasta files were missing or empty." % bad_file_count)
    elif bad_file_count > 0:
        utils.sample_error("Error: %d consensus fasta files were missing or empty." % bad_file_count, continue_possible=True)

    #==========================================================================
    # Check if the result is already fresh
    #==========================================================================
    snpma_file_path = args.snpmaFile
    source_files = consensus_files
    if not args.forceFlag:
        if not utils.target_needs_rebuild(source_files, snpma_file_path):
            utils.verbose_print("SNP matrix %s has already been freshly built.  Use the -f option to force a rebuild." % snpma_file_path)
            return

    #==========================================================================
    #   Create snp matrix. Write results to file.
    #==========================================================================
    with open(snpma_file_path, "w") as output_file:
        for consensus_file_path in consensus_files:
            utils.verbose_print("Merging " + consensus_file_path)
            with open(consensus_file_path, "r") as input_file:
                for line in input_file:
                    output_file.write(line)
