"""This module is part of the CFSAN SNP Pipeline. It contains the code to
prepare the snp list file with the sites having snps in any of the samples.
"""

from __future__ import print_function
from __future__ import absolute_import

import os
from snppipeline import utils


def merge_sites(args):
    """Create SNP list file

    Description:
    Create the SNP list -- the list of positions where variants were found
    and the corresponding list of samples having a variant at each position.
    This function expects, or creates '(*)', the following files arranged
    in the following way:
            sampleDirectories.txt
            samples
                sample_name_one/var.flt.vcf
                ...
            snplist.txt (*)

    The files are used as follows:
        1. The sampleDirectories.txt input file contains a list of the paths to
           the sample directories.
        2. The var.flt.vcf variant input files are used to construct the
           SNP position list.
        3. The snplist.txt output file contains the union of the SNP positions
           and sample names extracted from all the var.flt.vcf files.

    The sampleDirectories.txt and var.flt.vcf files are created outside of
    this function. The package documentation provides an example of creating
    these files based on the lambda_virus sequence that is used as one test
    for this package.

    Parameters
    ----------
    args : Namespace
        sampleDirsFile: File path (not just file name) of file containing paths
            to directories containing var.flt.vcf file for each sequence.
        vcfFileName: File name of the VCF files which must exist in each of the
            sample directories
        snpListFile: File path (not just file name) of text format list
            of SNP positions

    Raises:

    Examples:
    args = argparse.Namespace
    args.sampleDirsFile = 'sampleDirectories.txt'
    args.vcfFileName = 'var.flt.vcf'
    args.snpListFile = 'snplist.txt'
    merge_sites(args)
    """
    utils.print_log_header()
    utils.print_arguments(args)

    #==========================================================================
    # Prep work
    #==========================================================================
    sample_directories_list_path = args.sampleDirsFile
    bad_file_count = utils.verify_non_empty_input_files("File of sample directories", [sample_directories_list_path])
    if bad_file_count > 0:
        utils.global_error(None)

    with open(sample_directories_list_path, "r") as sample_directories_list_file:
        unsorted_list_of_sample_directories = [line.rstrip() for line in sample_directories_list_file]
    unsorted_list_of_sample_directories = [d for d in unsorted_list_of_sample_directories if d]
    sorted_list_of_sample_directories = sorted(unsorted_list_of_sample_directories)

    #==========================================================================
    # Validate inputs
    #==========================================================================
    snp_list_file_path = args.snpListFile
    vcf_file_name = args.vcfFileName
    list_of_vcf_files = [os.path.join(dir, vcf_file_name) for dir in sorted_list_of_sample_directories]

    bad_file_count = utils.verify_non_empty_input_files("VCF file", list_of_vcf_files)
    if bad_file_count == len(list_of_vcf_files):
        utils.global_error("Error: all %d VCF files were missing or empty." % bad_file_count)
    elif bad_file_count > 0:
        utils.sample_error("Error: %d VCF files were missing or empty." % bad_file_count, continue_possible=True)

    #==========================================================================
    # Read in all vcf files and process into dict of SNPs passing various
    # criteria. Do this for each sample. Write to file.
    #==========================================================================
    if args.forceFlag or utils.target_needs_rebuild(list_of_vcf_files, snp_list_file_path):
        snp_dict = dict()
        excluded_sample_directories = set()
        for sample_dir, vcf_file_path in zip(sorted_list_of_sample_directories, list_of_vcf_files):

            if not os.path.isfile(vcf_file_path):
                continue
            if os.path.getsize(vcf_file_path) == 0:
                continue

            utils.verbose_print("Processing VCF file %s" % vcf_file_path)
            sample_name = os.path.basename(os.path.dirname(vcf_file_path))
            snp_set = utils.convert_vcf_file_to_snp_set(vcf_file_path)
            max_snps = args.maxSnps
            if max_snps >= 0 and len(snp_set) > max_snps:
                utils.verbose_print("Excluding sample %s having %d snps." % (sample_name, len(snp_set)))
                excluded_sample_directories.add(sample_dir)
                continue

            for key in snp_set:
                if key not in snp_dict:
                    sample_list = [sample_name]
                    snp_dict[key] = sample_list
                else:
                    sample_list = snp_dict[key]
                    sample_list.append(sample_name)

        utils.verbose_print('Found %d snp positions across %d sample vcf files.' % (len(snp_dict), len(list_of_vcf_files)))
        utils.write_list_of_snps(snp_list_file_path, snp_dict)

        #==========================================================================
        # Write the filtered list of sample directories
        #==========================================================================
        #sample_directories_list_path = sample_directories_list_path + ".filtered"
        filtered_sample_directories_list_path = args.filteredSampleDirsFile
        with open(filtered_sample_directories_list_path, "w") as filtered_samples_file_object:
            # Loop over the unsorted list to keep the order of samples the same as the original.
            # This will keep the same HPC log file suffix number.
            for sample_dir in unsorted_list_of_sample_directories:
                if sample_dir not in excluded_sample_directories:
                    filtered_samples_file_object.write("%s\n" % sample_dir)
    else:
        utils.verbose_print("SNP list %s has already been freshly built.  Use the -f option to force a rebuild." % snp_list_file_path)
