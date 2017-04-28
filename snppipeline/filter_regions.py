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
    # Validate some parameters
    #==========================================================================
    edge_length = args.edgeLength
    window_size = args.windowSize
    max_num_snp = args.maxSNP

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

    input_file_list = list()
    out_group_list_path = args.outGroupFile
    sorted_list_of_outgroup_samples = list()
    if out_group_list_path is not None:
        bad_file_count = utils.verify_non_empty_input_files("File of outgroup samples", [out_group_list_path])
        if bad_file_count > 0:
            utils.global_error(None)
        try:
            #There are outgroup samples
            input_file_list.append(out_group_list_path)
            with open(out_group_list_path, "r") as out_group_list_file:
                unsorted_list_of_outgroup_samples = [line.rstrip() for line in out_group_list_file]
            sorted_list_of_outgroup_samples = sorted(unsorted_list_of_outgroup_samples)
        except:
            utils.global_error("Error: Cannot open the file containing the list of outgroup samples!")

    #==========================================================================
    # Validate inputs
    #==========================================================================
    vcf_file_name = args.vcfFileName
    list_of_vcf_files = [os.path.join(dir, vcf_file_name) for dir in sorted_list_of_sample_directories]
    input_file_list.extend(list_of_vcf_files)

    bad_file_count = utils.verify_non_empty_input_files("VCF file", list_of_vcf_files)
    if bad_file_count == len(list_of_vcf_files):
        utils.global_error("Error: all %d VCF files were missing or empty." % bad_file_count)
    elif bad_file_count > 0:
        utils.sample_error("Error: %d VCF files were missing or empty." % bad_file_count, continue_possible=True)

    bad_file_count = utils.verify_non_empty_input_files("Reference file", [args.refFastaFile])
    if bad_file_count > 0:
        utils.global_error(None)

    #==========================================================================
    # Get contigs' length from the reference fasta file
    #==========================================================================
    try:
        handle = open(args.refFastaFile, "r")
        contig_length_dict = dict()
        for record in SeqIO.parse(handle, "fasta"):
            #build contig_length_dict
            contig_length_dict[record.id] = len(record.seq)
        input_file_list.append(args.refFastaFile)
    except:
        utils.global_error("Error: cannot open the reference fastq file, or fail to read the contigs in the reference fastq file.")
    else:
        if handle:
            handle.close()

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
        need_rebuild_dict[vcf_file_path] = args.forceFlag or preserved_needs_rebuild or removed_needs_rebuild

    if not any(need_rebuild_dict.values()):
        utils.verbose_print("All preserved and removed vcf files are already freshly built.  Use the -f option to force a rebuild.")
        return

    #==========================================================================
    # Find all bad regions.
    #==========================================================================
    bad_regions_dict = dict() # Key is the contig ID, and the value is a list of bad regions.
    for vcf_file_path in list_of_vcf_files:
        try:
            vcf_reader_handle = open(vcf_file_path, 'r')
            vcf_reader = vcf.Reader(vcf_reader_handle)
        except:
            utils.sample_error("Error: Cannot open the input vcf file: %s." % vcf_file_path, continue_possible=True)
            continue

        #Get sample ID
        ss = vcf_file_path.split('/')
        sample_ID = ss[-2]

        if sample_ID in sorted_list_of_outgroup_samples:
            if not need_rebuild_dict[vcf_file_path]:
                vcf_reader_handle.close()
                continue
            #Copy original vcf file to _preserved.vcf, and created an empty _removed.vcf

            #SNP list, saved as (Contig_Name, [(SNP_Position, SNP_Record),]), where SNP_Record is a line in VCF.

            preserved_vcf_file_path = vcf_file_path[:-4] + "_preserved.vcf"
            removed_vcf_file_path = vcf_file_path[:-4] + "_removed.vcf"

            try:
                vcf_writer_removed = None
                vcf_writer_removed = vcf.Writer(open(removed_vcf_file_path, 'w'), vcf_reader)
            except:
                #print "Cannot create the file for removed SNPs: %d." % removed_vcf_file_path
                #close vcf_writer_reserved and remove the file reserved_vcf_file_path
                if vcf_writer_removed is not None:
                    vcf_writer_removed.close()
                os.remove(removed_vcf_file_path)
                vcf_reader_handle.close()
                utils.sample_error("Error: Cannot create the file for removed SNPs: %s." % removed_vcf_file_path, continue_possible=True)
                continue

            vcf_writer_removed.close()
            vcf_reader_handle.close()
            shutil.copyfile(vcf_file_path, preserved_vcf_file_path)
        else:
            #SNP list, saved as (Contig_Name, [(SNP_Position, SNP_Record),]), where SNP_Record is a line in VCF.
            snp_dict = defaultdict(list)
            for vcf_data_line in vcf_reader:
                #Create a dict to store all SNPs in this sample
                #get contig length from contig name.The CHROM should be a contig name in the format of Velvet/SPAdes output.
                record = (vcf_data_line.POS, vcf_data_line)
                snp_dict[vcf_data_line.CHROM].append(record)

            #Find bad regions and add them into bad_region
            for contig, snp_list in snp_dict.items():

                #sort all SNPs in this contig by position
                sorted_list = sorted(snp_list, key=lambda SNPs: SNPs[0])

                #total number of SNPs
                num_of_snp = len(sorted_list)

                if contig not in bad_regions_dict:
                    #New contig
                    try:
                        contig_length = contig_length_dict[contig]
                    except:
                        #cannot find contig length. Use the sys.maxsize.
                        contig_length = sys.maxsize

                    if (contig_length <= (edge_length * 2)):
                        bad_regions_dict[contig] = [(0, contig_length)]
                    else:
                        region = [(0, edge_length), (contig_length - edge_length, contig_length)]
                        bad_regions_dict[contig] = region

                #Process SNPs
                for idx, snp in enumerate(sorted_list):
                    if (idx + max_num_snp) < num_of_snp:
                        pos_start = snp[0]
                        pos_end = sorted_list[idx + max_num_snp][0]
                        if (pos_start + window_size) >= pos_end:
                            #Add bad region
                            regions = bad_regions_dict[contig]
                            temp_region = (pos_start, pos_end)
                            regions.append(temp_region)
        vcf_reader_handle.close()

    #Combine all bad regions for each contig
    for contig, regions in bad_regions_dict.items():
        sorted_regions = utils.sort_coord(regions)
        combined_regions = utils.consensus(sorted_regions)
        bad_regions_dict[contig] = combined_regions

    #Scan vcf files to remove SNPs
    for vcf_file_path in list_of_vcf_files:
        if not need_rebuild_dict[vcf_file_path]:
            continue
        #Get sample ID
        ss = vcf_file_path.split('/')
        sample_ID = ss[-2]

        if sample_ID not in sorted_list_of_outgroup_samples:
            try:
                vcf_reader_handle = open(vcf_file_path, 'r')
                vcf_reader = vcf.Reader(vcf_reader_handle)
            except:
                utils.sample_error("Error: Cannot open the input vcf file: %s." % vcf_file_path, continue_possible=True)
                continue

            #SNP list, saved as (Contig_Name, [(SNP_Position, SNP_Record),]), where SNP_Record is a line in VCF.

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
                continue

            try:
                vcf_writer_removed = None
                vcf_writer_removed = vcf.Writer(open(removed_vcf_file_path, 'w'), vcf_reader)
            except:
                #close vcf_writer_reserved and remove the file reserved_vcf_file_path
                if vcf_writer_removed is not None:
                    vcf_writer_removed.close()
                os.remove(removed_vcf_file_path)
                vcf_writer_preserved.close()
                vcf_reader_handle.close()
                utils.sample_error("Error: Cannot create the file for removed SNPs: %s." % removed_vcf_file_path, continue_possible=True)
                continue

            for vcf_data_line in vcf_reader:
                #Create a dict to store all SNPs in this sample
                #get contig length from contig name.The CHROM should be a contig name in the format of Velvet/SPAdes output.
                contig = vcf_data_line.CHROM
                if utils.in_region(vcf_data_line.POS, bad_regions_dict[contig]):
                    #Remove this SNP
                    vcf_writer_removed.write_record(vcf_data_line)
                else:
                    #Preserve this SNP
                    vcf_writer_preserved.write_record(vcf_data_line)

            vcf_writer_preserved.close()
            vcf_writer_removed.close()
            vcf_reader_handle.close()
