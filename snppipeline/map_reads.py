"""This module is part of the CFSAN SNP Pipeline. It contains the code to
align reads to a reference genome using an external mapper -- bowtie2 or smalt.
Read group tags are also assigned during this step.
"""

from __future__ import print_function
from __future__ import absolute_import

import os
import psutil
import re
import shutil
import sys

from snppipeline import command
from snppipeline import fastq
from snppipeline import utils
from snppipeline.utils import verbose_print


def map_reads(args):
    """Align reads to the reference.

    Execute an external program (bowtie2 or smalt) to map the fastq reads
    to a reference file.

    The environment variable SnpPipeline_Aligner selects between bowtie2 and smalt.

    This function expects, or creates '(*)', the following files arranged
    in the following way:
            reference
                referenceFile.fasta
            samples
                sample_name_one/sampleFastqFile_1.fastq
                sample_name_one/sampleFastqFile_2.fastq
                sample_name_one/reads.sam*

    The reverse fastq file is optional.
    The fastq files may be either compressed with gzip or uncompressed.

    The reverse fastq file is optional.

    All the input files are created outside of this function. The package
    documentation provides an example of preparing these files based on the
    lambda_virus sequence that is used as one test for this package.

    Parameters
    ----------
    args : argparse.Namespace
        referenceFile : File path of the reference fasta file
        sampleFastqFile1 : File path of the forward fastq file
        sampleFastqFile2 : Optional file path of the reverse fastq file
    """
    utils.print_log_header()
    utils.print_arguments(args)

    #==========================================================================
    # Validate inputs
    #==========================================================================

    # Verify reference fasta file exists and is not empty
    reference_file_path = args.referenceFile
    utils.verify_non_empty_input_files("Reference file", [reference_file_path], error_handler="global")

    # Verify fastq files exist and are not empty
    sample_fastq_file1 = args.sampleFastqFile1
    sample_fastq_file2 = args.sampleFastqFile2
    fastq_files = [sample_fastq_file1]
    if sample_fastq_file2:
        fastq_files.append(sample_fastq_file2)

    utils.verify_non_empty_input_files("Sample file", fastq_files, error_handler="sample")

    # The environment variable SnpPipeline_Aligner selects between bowtie2 and smalt
    snp_pipeline_aligner = os.environ.get("SnpPipeline_Aligner") or "bowtie2"
    snp_pipeline_aligner = snp_pipeline_aligner.lower()
    if snp_pipeline_aligner not in ["bowtie2", "smalt"]:
        utils.global_error("Error: only bowtie2 and smalt aligners are supported.")

    sample_dir = os.path.dirname(sample_fastq_file1)
    sample_id = utils.sample_id_from_file(sample_fastq_file1)
    reference_base_path = os.path.splitext(reference_file_path)[0] # strip the file extension
    reference_id = os.path.basename(reference_base_path)

    #==========================================================================
    # Check if alignment to reference has already been done
    #==========================================================================
    sam_file = os.path.join(sample_dir, "reads.sam")
    source_files = [sample_fastq_file1]
    if sample_fastq_file2:
        source_files.append(sample_fastq_file2)
    if snp_pipeline_aligner == "bowtie2":
        source_files.append(reference_base_path + ".rev.1.bt2")
    elif snp_pipeline_aligner == "smalt":
        source_files.append(reference_base_path + ".smi")
    needs_rebuild = utils.target_needs_rebuild(source_files, sam_file)

    if not args.forceFlag and not needs_rebuild:
        verbose_print("# %s has already been aligned to %s.  Use the -f option to force a rebuild." % (sample_id, reference_id))
        return

    #==========================================================================
    # Construct the command line to execute bowtie2 or smalt
    #==========================================================================

    # The read group identifies reads from a single run and lane
    read_group_tags = fastq.construct_read_group_tags(sample_fastq_file1, sample_id)

    # Default to 8 cores on HPC or all cpu cores on workstation
    if os.environ.get("JOB_ID") or os.environ.get("PBS_JOBID"):
        num_cores = 8
    else:
        num_cores = psutil.cpu_count()

    num_cores_param = ""

    if snp_pipeline_aligner == "bowtie2":
        version_str = utils.extract_version_str("bowtie2", "bowtie2 --version")

        # Parse the user-specified bowtie parameters to determine if the user specified the number of CPU cores
        bowtie2_align_extra_params = os.environ.get("Bowtie2Align_ExtraParams") or ""
        if not utils.detect_numeric_option_in_parameters_str(bowtie2_align_extra_params, "-p"):
            num_cores_param = "-p " + str(num_cores)

        # Specify the read group and sample tags here, --rg tags cannot be specified without ID.
        # The read group tags are used by some downstream tools, like Picard and GATK.
        read_group_params = ""
        if read_group_tags:
            read_group_params += " --rg-id " + read_group_tags.ID
            read_group_params += " --rg SM:" + read_group_tags.SM
            read_group_params += " --rg LB:" + read_group_tags.LB
            read_group_params += " --rg PL:" + read_group_tags.PL
            read_group_params += " --rg PU:" + read_group_tags.PU

        # Substitute the default parameters if the user did not specify bowtie parameters
        bowtie2_align_params = bowtie2_align_extra_params or "--reorder -q"

        # Build the command with options depending on whether the fastq files are paired
        command_line = "bowtie2 " + num_cores_param + " " + read_group_params + " " + bowtie2_align_params + " -x " + reference_base_path
        if sample_fastq_file2:
            command_line += " -1 " + sample_fastq_file1 + " -2 " + sample_fastq_file2
        else:
            command_line += " -U " + sample_fastq_file1

    elif snp_pipeline_aligner == "smalt":
        version_str = utils.extract_version_str("smalt", "smalt version")

        # Parse the user-specified smalt parameters to determine if the user specified the number of CPU cores
        smalt_align_extra_params = os.environ.get("SmaltAlign_ExtraParams") or ""
        if not utils.detect_numeric_option_in_parameters_str(smalt_align_extra_params, "-n"):
            num_cores_param = "-n " + str(num_cores)

        # Substitute the default parameters if the user did not specify smalt parameters
        smalt_align_params = smalt_align_extra_params or "-O"

        # Don't use the -i 1000 option if the fastq file is unpaired
        if not sample_fastq_file2:
            smalt_align_params = re.sub("-i[ ]+[0-9]+", '', smalt_align_extra_params) # regex substitute

        command_line = "smalt map " + num_cores_param + " " + smalt_align_params + " " + reference_base_path + " " + sample_fastq_file1 + " " + (sample_fastq_file2 or "")

    #==========================================================================
    # Run the command to execute bowtie2 or smalt
    #==========================================================================
    verbose_print("# Align sequence %s to reference %s" % (sample_id, reference_id))
    verbose_print("# %s %s" % (utils.timestamp(), command_line))
    verbose_print("# %s" % version_str)
    command.run(command_line, sam_file)

    #==========================================================================
    # When using smalt, assign read groups in a separate step.
    # This is already done when using bowtie2.
    #==========================================================================
    if snp_pipeline_aligner == "smalt" and read_group_tags:
        smalt_sam_file = os.path.join(sample_dir, "reads.smalt.sam")
        shutil.move(sam_file, smalt_sam_file)
        version_str = utils.extract_version_str("Picard", "java  picard.cmdline.PicardCommandLine AddOrReplaceReadGroups --version 2>&1")
        jvm_params = os.environ.get("PicardJvm_ExtraParams") or ""
        command_line = "java " + jvm_params + " picard.cmdline.PicardCommandLine AddOrReplaceReadGroups"
        command_line += " I=" + smalt_sam_file
        command_line += " O=" + sam_file
        command_line += " RGID=" + read_group_tags.ID
        command_line += " RGSM=" + read_group_tags.SM
        command_line += " RGLB=" + read_group_tags.LB
        command_line += " RGPL=" + read_group_tags.PL
        command_line += " RGPU=" + read_group_tags.PU
        verbose_print("")
        verbose_print("# Assign read group id %s" % (read_group_tags.ID))
        verbose_print("# %s %s" % (utils.timestamp(), command_line))
        verbose_print("# %s" % version_str)
        command.run(command_line, sys.stdout)
