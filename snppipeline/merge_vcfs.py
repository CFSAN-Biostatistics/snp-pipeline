"""This module is part of the CFSAN SNP Pipeline. It contains the code to create
a multi-sample VCF file from the per-sample consensus.vcf files.
"""

from __future__ import print_function
from __future__ import absolute_import

import os
import shutil
import sys
import tempfile

from snppipeline import command
from snppipeline import utils
from snppipeline.utils import verbose_print

def merge_vcfs(args):
    """Merge the per-sample VCF files.

    Execute an external program (bcftools merge)) to merge the VCF files.

    This function expects, or creates '(*)', the following files arranged
    in the following way:
            samples
                sample_name_one/consensus.vcf
            snpma.vcf*

    All the input files are created outside of this function.  Before
    running this command, the vcf file for each sample must be created by the
    call_consensus.py script.

    The package documentation provides an example of preparing these files based
    on the lambda_virus sequence that is used as one test for this package.

    Parameters
    ----------
    args : argparse.Namespace
        sampleDirsFile : Path to file containing a list of directories -- one per sample
        vcfFileName : File name of the vcf files which must exist in each of the sample directories
        mergedVcfFile : Path to the output merged multi-vcf file
    """
    utils.print_log_header()
    utils.print_arguments(args)

    #==========================================================================
    # Validate inputs
    #==========================================================================

    sample_directories_list_path = args.sampleDirsFile
    vcf_file_name = args.vcfFileName
    merged_vcf_file = args.mergedVcfFile

    utils.verify_non_empty_input_files("File of sample directories", [sample_directories_list_path], error_handler="global")

    with open(sample_directories_list_path, "r") as f:
        sample_directories = [line.rstrip() for line in f]
    sample_directories = [d for d in sample_directories if d]
    vcf_files = [os.path.join(d, vcf_file_name) for d in sample_directories]

    good_vcf_files = []
    for vcf_file in vcf_files:
        bad = utils.verify_non_empty_input_files("Sample vcf file", [vcf_file], error_handler="sample", continue_possible=True)
        if not bad:
            good_vcf_files.append(vcf_file)

    if len(good_vcf_files) == 0:
        utils.global_error("There are no vcf files to merge.")

    #==========================================================================
    # Check if merge has already been done
    #==========================================================================
    needs_rebuild = utils.target_needs_rebuild(vcf_files, merged_vcf_file)
    if not args.forceFlag and not needs_rebuild:
        verbose_print("# Multi-VCF file is already freshly created.  Use the -f option to force a rebuild.")
        return

    #==========================================================================
    # Copy, Compress, Index, Merge
    #==========================================================================

    # If there is only one good sample, just copy the consensus VCF file to the snpma.vcf file
    if len(good_vcf_files) ==  1:
        shutil.copy(good_vcf_files[0], merged_vcf_file)
        return

    # Copy single VCF files to a common directory where the files will be edited
    verbose_print("# %s Copying VCF files to temp directory" % utils.timestamp())
    parent_of_temp_dir = os.path.dirname(merged_vcf_file)
    temp_dir = tempfile.mkdtemp(prefix="tmp.vcf.", dir=parent_of_temp_dir)
    file_copies = []
    for d in sample_directories:
        src_file = os.path.join(d, vcf_file_name)
        if src_file in good_vcf_files:
            dst_file = os.path.join(temp_dir, os.path.basename(d) + ".vcf")
            file_copies.append(dst_file)
            verbose_print("copy %s %s" % (src_file, dst_file))
            #if not os.path.isfile(dst_file) or os.stat(src_file).st_mtime > os.stat(dst_file).st_mtime:
            shutil.copy2(src_file, dst_file)

    # bgzip all the sample vcf files
    verbose_print("# %s Compressing VCF files" % utils.timestamp())
    for file in file_copies:
        verbose_print("bgzip -c %s > %s" % (file, file + ".gz"))
        command.run("bgzip -c " + file, file + ".gz")

    # Index all the zipped sample vcf file
    verbose_print("# %s Indexing VCF files" % utils.timestamp())
    for file in file_copies:
        file += ".gz"
        verbose_print("tabix -f -p vcf " + file)
        command.run("tabix -f -p vcf " + file, sys.stdout)

    # Substitute the default parameters if the user did not specify bcftools parameters
    default_params = "--merge all --info-rules NS:sum"
    bcf_tools_extra_params = os.environ.get("BcftoolsMerge_ExtraParams") or default_params

    # Merge the VCFs
    verbose_print("# %s Merging VCF files" % utils.timestamp())
    command_line = "bcftools merge -o " + merged_vcf_file + ' ' + bcf_tools_extra_params + ' ' + temp_dir + "/*.gz"
    verbose_print(command_line)
    command.run(command_line, sys.stdout)

    # Clean up
    shutil.rmtree(temp_dir)
