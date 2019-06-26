"""This module is part of the CFSAN SNP Pipeline. It contains the code to
remove the abnormal snps.
"""

from __future__ import print_function
from __future__ import absolute_import

from Bio import SeqIO
from collections import defaultdict
import os
import sys
import shutil
import vcf
from snppipeline import utils


def find_dense_regions(max_allowed_snps, window_size, snps):
    """Scan a list of snp positions to find regions where the snp density exceeds the
    allowed thershold.

    Parameters
    ----------
    max_allowed_snps : int
        Maximum allowed number of snps in a given rolling window.
    window_size : int
        Size of rolling window along the length of a genome which
        is scanned for excessive snps.
    snps : list of int
        Sorted list of snp positions

    Returns
    -------
    dense_region_list : list of tuples
        List of (start_position, end_position) tuples identifying the list of
        dense snp regions.

    Examples
    --------
    # Empty list
    >>> find_dense_regions(3, 1000, [])
    []

    # Not dense window
    >>> find_dense_regions(3, 1000, [1, 2, 3, 1001])
    []

    # One more than max_allowed_snps at window boundaries
    >>> find_dense_regions(3, 1000, [1, 20, 30, 1000])
    [(1, 1000)]

    # Two more than max_allowed_snps at window boundaries
    >>> find_dense_regions(3, 1000, [1, 20, 30, 40, 1000])
    [(1, 1000)]

    # Overlapping dense regions with combined size greater than window_size
    >>> find_dense_regions(3, 1000, [1, 20, 30, 40, 501, 600, 1000, 1500])
    [(1, 1500)]

    # Multiple dense regions
    >>> find_dense_regions(3, 1000, [1, 2, 3, 1000, 1500, 3001, 3002, 3003, 4000])
    [(1, 1000), (3001, 4000)]
    """
    snp_count = len(snps)
    dense_region_list = []
    for idx, pos_start in enumerate(snps):
        if (idx + max_allowed_snps) < snp_count:
            pos_end = snps[idx + max_allowed_snps]
            if (pos_start + window_size - 1) >= pos_end:
                dense_region_list.append((pos_start, pos_end))
    dense_region_list = utils.merge_regions(dense_region_list)
    return dense_region_list


def filter_regions(args):
    """Remove bad SNPs from original vcf files

    Remove bad SNPs -- this function finds bad regions, including the edges
    and probable prophage regions; then remove SNPs in these regions in
    original vcf files of all samples.

    This function expects, or creates '(*)', the following files arranged
    in the following way:
            sampleDirectories.txt
            samples
                sample_name_one/var.flt.vcf
                sample_name_one/var.flt_removed.vcf (*)
                sample_name_one/var.flt_preserved.vcf (*)
                ...

    The files are used as follows:
        1. The sampleDirectories.txt input file contains a list of the paths to
           the sample directories.
        2. The var.flt.vcf variant input files (i.e., the original vcf file).
        3. The var.flt_removed.vcf and var.flt_preserved.vcf output files contain the removed SNPs and
           preserved SNPs.

    The sampleDirectories.txt and var.flt.vcf files are created outside of
    this function. The package documentation provides an example of creating
    these files based on the lambda_virus sequence that is used as one test
    for this package.

    Parameters
    ----------
    Args:
        sampleDirsFile: File path (not just file name) of file containing paths
            to directories containing var.flt.vcf file for each sequence.
        vcfFileName: File name of the VCF files which must exist in each of the
            sample directories
        refFastaFile: File path (not just file name) of reference fasta file
        edgeLength: the length of edge of a contig in which SNPs will be removed.
            Default is 500.
        windowSize: the size of the window in which max number of SNPs are allowed.
            Default is 1000.
        maxSNP: the maximum number of SNPs allowed in a window of a size defined in
            windowSize. Default is 3.
        mode:
            all = Dense regions found in any sample are filtered from all samples.
            each = Dense regions found in any sample are filtered independently from samples.

    Raises:

    Examples:
    args = argparse.Namespace
    args.sampleDirsFile = 'sampleDirectories.txt'
    args.vcfFileName = 'var.flt.vcf'
    args.refFastaFile = 'snplist.txt'
    remove_bad_snp(args)
    """
    utils.print_log_header()
    utils.print_arguments(args)

    #==========================================================================
    # Get arguments from Argparse namespace
    #==========================================================================
    sample_directories_list_path = args.sampleDirsFile
    ref_fasta_path = args.refFastaFile
    force_flag = args.forceFlag
    vcf_file_name = args.vcfFileName
    edge_length = args.edgeLength
    window_size_list = args.windowSizeList
    max_num_snps_list = args.maxSnpsList
    out_group_list_path = args.outGroupFile
    filter_across_samples = args.mode == "all"

    #==========================================================================
    # Validate inputs
    #==========================================================================
    bad_file_count = utils.verify_non_empty_input_files("File of sample directories", [sample_directories_list_path])
    if bad_file_count > 0:
        utils.global_error(None)

    with open(sample_directories_list_path, "r") as sample_directories_list_file:
        unsorted_list_of_sample_directories = [line.rstrip() for line in sample_directories_list_file]
    unsorted_list_of_sample_directories = [d for d in unsorted_list_of_sample_directories if d]
    sorted_list_of_sample_directories = sorted(unsorted_list_of_sample_directories)

    list_of_vcf_files = [os.path.join(dir, vcf_file_name) for dir in sorted_list_of_sample_directories]
    bad_file_count = utils.verify_non_empty_input_files("VCF file", list_of_vcf_files)
    if bad_file_count == len(list_of_vcf_files):
        utils.global_error("Error: all %d VCF files were missing or empty." % bad_file_count)
    elif bad_file_count > 0:
        utils.sample_error("Error: %d VCF files were missing or empty." % bad_file_count, continue_possible=True)

    bad_file_count = utils.verify_non_empty_input_files("Reference file", [ref_fasta_path])
    if bad_file_count > 0:
        utils.global_error(None)

    sorted_list_of_outgroup_samples = list()
    if out_group_list_path is not None:
        bad_file_count = utils.verify_non_empty_input_files("File of outgroup samples", [out_group_list_path])
        if bad_file_count > 0:
            utils.global_error(None)
        try:
            # There are outgroup samples
            with open(out_group_list_path, "r") as out_group_list_file:
                unsorted_list_of_outgroup_samples = [line.rstrip() for line in out_group_list_file]
            sorted_list_of_outgroup_samples = sorted(unsorted_list_of_outgroup_samples)
        except:
            utils.global_error("Error: Cannot open the file containing the list of outgroup samples!")

    #==========================================================================
    # Get contigs' length from the reference fasta file
    #==========================================================================
    try:
        handle = open(ref_fasta_path, "r")
        contig_length_dict = dict()
        for record in SeqIO.parse(handle, "fasta"):
            # build contig_length_dict
            contig_length_dict[record.id] = len(record.seq)
    except:
        utils.global_error("Error: cannot open the reference fastq file, or fail to read the contigs in the reference fastq file.")
    else:
        if handle:
            handle.close()

    #==========================================================================
    # Filter regions
    #==========================================================================
    if filter_across_samples:
        filter_regions_across_samples(list_of_vcf_files, contig_length_dict, sorted_list_of_outgroup_samples, force_flag, edge_length, window_size_list, max_num_snps_list, ref_fasta_path, out_group_list_path)
    else:
        filter_regions_per_sample(list_of_vcf_files, contig_length_dict, sorted_list_of_outgroup_samples, force_flag, edge_length, window_size_list, max_num_snps_list, ref_fasta_path, out_group_list_path)


def filter_regions_across_samples(list_of_vcf_files, contig_length_dict, sorted_list_of_outgroup_samples, force_flag, edge_length, window_size_list, max_num_snps_list, ref_fasta_path, out_group_list_path):
    """Detect abnormal regions in each sample and filter those regions from all samples.

    Parameters
    ----------
    list_of_vcf_files : list of str
        List of input VCF file paths -- one per sample.
    contig_length_dict : dict, str --> int
        Mapping of contig id to int length of contig.
    sorted_list_of_outgroup_samples : list of str
        List of sample IDs for samples that are outgroup samples.
    force_flag : bool
        Force processing even when result files already exist and are newer than inputs.
    edge_length : int
        The length of the edge regions in a contig, in which all SNPs will be removed.
    window_size_list : list of int
        The length of the window in which the number of SNPs should be no more than max_num_snp.
    max_num_snps_list : list of int
        The maximum number of SNPs allowed in a window.  This list has the same size as window_size_list
        and the entries correspond to one another.
    ref_fasta_path : str
        Path to the reference fasta file.
    out_group_list_path : str
        Path to the file indicating outgroup samples, one sample ID per line.
    """
    #==========================================================================
    # Prep work
    #==========================================================================
    input_file_list = list()
    input_file_list.append(ref_fasta_path)
    if out_group_list_path:
        input_file_list.append(out_group_list_path)
    input_file_list.extend(list_of_vcf_files)

    #==========================================================================
    # Which samples need rebuild?
    #
    # Any changed or new input file will trigger rebuild for all samples because
    # the bad regions are combined across all samples.  However, a missing
    # output file will only cause rebuild of the missing file.
    #==========================================================================
    need_rebuild_dict = dict()
    for vcf_file_path in list_of_vcf_files:
        preserved_vcf_file_path = vcf_file_path[:-4] + "_preserved.vcf"
        removed_vcf_file_path = vcf_file_path[:-4] + "_removed.vcf"
        preserved_needs_rebuild = utils.target_needs_rebuild(input_file_list, preserved_vcf_file_path)
        removed_needs_rebuild = utils.target_needs_rebuild(input_file_list, removed_vcf_file_path)
        need_rebuild_dict[vcf_file_path] = force_flag or preserved_needs_rebuild or removed_needs_rebuild

    if not any(need_rebuild_dict.values()):
        utils.verbose_print("All preserved and removed vcf files are already freshly built.  Use the -f option to force a rebuild.")
        return

    #==========================================================================
    # Find all bad regions.
    #==========================================================================
    # The bad_regions_dict holds all the bad regions across all samples mixed together.
    # Key is the contig ID, and the value is a list of bad region tuples (start_position, end_position).
    bad_regions_dict = dict()
    for vcf_file_path in list_of_vcf_files:
        try:
            vcf_reader_handle = open(vcf_file_path, 'r')
            vcf_reader = vcf.Reader(vcf_reader_handle)
        except:
            utils.sample_error("Error: Cannot open the input vcf file: %s." % vcf_file_path, continue_possible=True)
            continue

        sample_ID = utils.sample_id_from_file(vcf_file_path)
        utils.verbose_print("Processing sample %s" % sample_ID)
        if sample_ID in sorted_list_of_outgroup_samples:
            write_outgroup_preserved_and_removed_vcf_files(vcf_file_path, vcf_reader)
        else:
            collect_dense_regions(vcf_reader, bad_regions_dict, contig_length_dict, edge_length, max_num_snps_list, window_size_list)
        vcf_reader_handle.close()

    # Combine all bad regions for each contig
    for contig, regions in bad_regions_dict.items():
        combined_regions = utils.merge_regions(regions)
        bad_regions_dict[contig] = combined_regions

    #==========================================================================
    # Write the output files
    #==========================================================================
    # Scan vcf files to remove SNPs
    for vcf_file_path in list_of_vcf_files:
        if not need_rebuild_dict[vcf_file_path]:
            continue
        sample_ID = utils.sample_id_from_file(vcf_file_path)

        if sample_ID in sorted_list_of_outgroup_samples:
            continue

        write_preserved_and_removed_vcf_files(vcf_file_path, bad_regions_dict)


def filter_regions_per_sample(list_of_vcf_files, contig_length_dict, sorted_list_of_outgroup_samples, force_flag, edge_length, window_size_list, max_num_snps_list, ref_fasta_path, out_group_list_path):
    """Detect abnormal regions in each sample and filter those regions from all samples.

    Parameters
    ----------
    list_of_vcf_files : list of str
        List of input VCF file paths -- one per sample.
    contig_length_dict : dict, str --> int
        Mapping of contig id to int length of contig.
    sorted_list_of_outgroup_samples : list of str
        List of sample IDs for samples that are outgroup samples.
    force_flag : bool
        Force processing even when result files already exist and are newer than inputs.
    edge_length : int
        The length of the edge regions in a contig, in which all SNPs will be removed.
    window_size_list : list of int
        The length of the window in which the number of SNPs should be no more than max_num_snp.
    max_num_snps_list : list of int
        The maximum number of SNPs allowed in a window.  This list has the same size as window_size_list
        and the entries correspond to one another.
    ref_fasta_path : str
        Path to the reference fasta file.
    out_group_list_path : str
        Path to the file indicating outgroup samples, one sample ID per line.
    """
    #==========================================================================
    # Prep work
    #==========================================================================
    input_file_list = list()
    input_file_list.append(ref_fasta_path)
    if out_group_list_path:
        input_file_list.append(out_group_list_path)

    #==========================================================================
    # Which samples need rebuild?
    #
    # Any changed or new input file will trigger rebuild for all samples because
    # the bad regions are combined across all samples.  However, a missing
    # output file will only cause rebuild of the missing file.
    #==========================================================================
    need_rebuild_dict = dict()
    for vcf_file_path in list_of_vcf_files:
        preserved_vcf_file_path = vcf_file_path[:-4] + "_preserved.vcf"
        removed_vcf_file_path = vcf_file_path[:-4] + "_removed.vcf"
        input_files = input_file_list + [vcf_file_path]
        preserved_needs_rebuild = utils.target_needs_rebuild(input_files, preserved_vcf_file_path)
        removed_needs_rebuild = utils.target_needs_rebuild(input_files, removed_vcf_file_path)
        need_rebuild_dict[vcf_file_path] = force_flag or preserved_needs_rebuild or removed_needs_rebuild

    if not any(need_rebuild_dict.values()):
        utils.verbose_print("All preserved and removed vcf files are already freshly built.  Use the -f option to force a rebuild.")
        return

    #==========================================================================
    # Find all bad regions in one sample at a time
    #==========================================================================
    for vcf_file_path in list_of_vcf_files:
        if not need_rebuild_dict[vcf_file_path]:
            continue
        try:
            vcf_reader_handle = open(vcf_file_path, 'r')
            vcf_reader = vcf.Reader(vcf_reader_handle)
        except:
            utils.sample_error("Error: Cannot open the input vcf file: %s." % vcf_file_path, continue_possible=True)
            continue

        sample_ID = utils.sample_id_from_file(vcf_file_path)
        utils.verbose_print("Processing sample %s" % sample_ID)
        if sample_ID in sorted_list_of_outgroup_samples:
            write_outgroup_preserved_and_removed_vcf_files(vcf_file_path, vcf_reader)
        else:
            # The bad_regions_dict holds the bad regions for this sample
            # Key is the contig ID, and the value is a list of bad region tuples (start_position, end_position).
            bad_regions_dict = dict()
            collect_dense_regions(vcf_reader, bad_regions_dict, contig_length_dict, edge_length, max_num_snps_list, window_size_list)

            # Combine all bad regions for each contig
            for contig, regions in bad_regions_dict.items():
                combined_regions = utils.merge_regions(regions)
                bad_regions_dict[contig] = combined_regions

            # Write the output files
            write_preserved_and_removed_vcf_files(vcf_file_path, bad_regions_dict)

        vcf_reader_handle.close()


def collect_dense_regions(vcf_reader, bad_regions_dict, contig_length_dict, edge_length, max_num_snps_list, window_size_list):
    """Collect the abnormal regions in a VCF file and store the results in the bad_regions_dict.

    Parameters
    ----------
    vcf_reader : PyVcf vcf.Reader
        Previously opened VCF reader object.
    bad_regions_dict : dict
        Key is the contig ID, and the value is a list of bad region tuples (start_position, end_position).
        This dictionary is modified by this function.
    contig_length_dict : dict, str --> int
        Mapping of contig id to int length of contig.
    edge_length : int
        The length of the edge regions in a contig, in which all SNPs will be removed.
    window_size_list : list of int
        The length of the window in which the number of SNPs should be no more than max_num_snp.
    max_num_snps_list : list of int
        The maximum number of SNPs allowed in a window.  This list has the same size as window_size_list
        and the entries correspond to one another.
    """
    # snp_dict key is contig name, value is list of positions
    # The CHROM should be a contig name in the format of Velvet/SPAdes output.
    snp_dict = defaultdict(list)
    for vcf_data_line in vcf_reader:
        snp_dict[vcf_data_line.CHROM].append(vcf_data_line.POS)

    # Find bad regions and add them into bad_region_dict
    for contig, snp_list in snp_dict.items():

        # First collect the ends of each contig
        if contig not in bad_regions_dict:
            contig_length = contig_length_dict.get(contig, sys.maxsize)

            if (contig_length <= (edge_length * 2)):
                bad_regions_dict[contig] = [(0, contig_length)]
            else:
                bad_regions_dict[contig] = [(0, edge_length), (contig_length - edge_length, contig_length)]

        # Process dense snp regions
        sorted_snps = sorted(snp_list)
        for max_allowed_snps, window_size in zip(max_num_snps_list, window_size_list):
            dense_regions = find_dense_regions(max_allowed_snps, window_size, sorted_snps)
            bad_regions_dict[contig].extend(dense_regions)


def write_outgroup_preserved_and_removed_vcf_files(vcf_file_path, vcf_reader):
    """The dense snps are not filtered from outgroup samples.  Instead, we
    copy the original vcf file to _preserved.vcf, and create an empty _removed.vcf.

    Parameters
    ----------
    vcf_file_path : str
        Path to a sample VCF file.
    vcf_reader : PyVcf vcf.Reader
        Previously opened VCF reader object.
    """
    preserved_vcf_file_path = vcf_file_path[:-4] + "_preserved.vcf"
    removed_vcf_file_path = vcf_file_path[:-4] + "_removed.vcf"

    try:
        vcf_writer_removed = None
        vcf_writer_removed = vcf.Writer(open(removed_vcf_file_path, 'w'), vcf_reader)
    except:
        # close vcf_writer_reserved and remove the file reserved_vcf_file_path
        if vcf_writer_removed is not None:
            vcf_writer_removed.close()
        os.remove(removed_vcf_file_path)
        utils.sample_error("Error: Cannot create the file for removed SNPs: %s." % removed_vcf_file_path, continue_possible=True)
        return

    vcf_writer_removed.close()
    shutil.copyfile(vcf_file_path, preserved_vcf_file_path)


def write_preserved_and_removed_vcf_files(vcf_file_path, bad_regions_dict):
    """Given a VCF file and a collection of abnormal regions, scan the snps in
    the VCF file and write each snps to either the preserved or removed output VCF file.

    Parameters
    ----------
    vcf_file_path : str
        Path to a sample VCF file.
    bad_regions_dict : dict
        Key is the contig ID, and the value is a list of bad region tuples (start_position, end_position).
    """
    try:
        vcf_reader_handle = open(vcf_file_path, 'r')
        vcf_reader = vcf.Reader(vcf_reader_handle)
    except:
        utils.sample_error("Error: Cannot open the input vcf file: %s." % vcf_file_path, continue_possible=True)
        return

    # SNP list, saved as (Contig_Name, [(SNP_Position, SNP_Record),]), where SNP_Record is a line in VCF.

    preserved_vcf_file_path = vcf_file_path[:-4] + "_preserved.vcf"
    removed_vcf_file_path = vcf_file_path[:-4] + "_removed.vcf"

    try:
        vcf_writer_preserved = None
        vcf_writer_preserved = vcf.Writer(open(preserved_vcf_file_path, 'w'), vcf_reader)
    except:
        if vcf_writer_preserved is not None:
            vcf_writer_preserved.close()
        os.remove(preserved_vcf_file_path)
        vcf_reader_handle.close()
        utils.sample_error("Error: Cannot create the file for preserved SNPs: %s." % preserved_vcf_file_path, continue_possible=True)
        return

    try:
        vcf_writer_removed = None
        vcf_writer_removed = vcf.Writer(open(removed_vcf_file_path, 'w'), vcf_reader)
    except:
        # close vcf_writer_reserved and remove the file reserved_vcf_file_path
        if vcf_writer_removed is not None:
            vcf_writer_removed.close()
        os.remove(removed_vcf_file_path)
        vcf_writer_preserved.close()
        vcf_reader_handle.close()
        utils.sample_error("Error: Cannot create the file for removed SNPs: %s." % removed_vcf_file_path, continue_possible=True)
        return

    for vcf_data_line in vcf_reader:
        # Create a dict to store all SNPs in this sample
        # get contig length from contig name.The CHROM should be a contig name in the format of Velvet/SPAdes output.
        contig = vcf_data_line.CHROM
        if utils.in_region(vcf_data_line.POS, bad_regions_dict[contig]):
            # Remove this SNP
            vcf_writer_removed.write_record(vcf_data_line)
        else:
            # Preserve this SNP
            vcf_writer_preserved.write_record(vcf_data_line)

    vcf_writer_preserved.close()
    vcf_writer_removed.close()
    vcf_reader_handle.close()
