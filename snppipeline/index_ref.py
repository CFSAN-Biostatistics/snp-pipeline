"""This module is part of the CFSAN SNP Pipeline. It contains the code to
Index the reference genome for subsequent read mapping, and create the faidx
index file for subsequent pileups.
"""

from __future__ import print_function
from __future__ import absolute_import

import os
import sys

from snppipeline import command
from snppipeline import utils
from snppipeline.utils import verbose_print


def index_ref(args):
    """Index the reference genome.

    Execute an external program (bowtie2 or smalt) to create an index for the
    reference genome to be used during subsequent alignment.  Execute samtools
    to create the faidx index file to be used during subsequent pileups.

    The environment variable SnpPipeline_Aligner selects between bowtie2 and smalt.

    This function expects, or creates '(*)', the following files arranged
    in the following way:
            reference
                referenceFile.fasta         # input fasta
                referenceFile.#.bt2*        # bowtie2 output
                referenceFile.rev.#.bt2*    # bowtie2 output
                referenceFile.sma*          # smalt output
                referenceFile.smi*          # smalt output
                referenceFile.fasta.fai*    # samtools faidx output

    The input fasta file is created outside of this function. The package
    documentation provides an example of preparing these files based on the
    lambda_virus sequence that is used as one test for this package.

    Parameters
    ----------
    args : argparse.Namespace
        referenceFile : File path of the reference fasta file
    """
    utils.print_log_header()
    utils.print_arguments(args)

    #==========================================================================
    # Validate inputs
    #==========================================================================

    # Verify reference fasta file exists and is not empty
    reference_file_path = args.referenceFile
    utils.verify_non_empty_input_files("Reference file", [reference_file_path], error_handler="global")

    reference_base_path = os.path.splitext(reference_file_path)[0] # strip the file extension

    # The environment variable SnpPipeline_Aligner selects between bowtie2 and smalt
    snp_pipeline_aligner = os.environ.get("SnpPipeline_Aligner") or "bowtie2"
    snp_pipeline_aligner = snp_pipeline_aligner.lower()
    if snp_pipeline_aligner not in ["bowtie2", "smalt"]:
        utils.global_error("Error: only bowtie2 and smalt aligners are supported.")

    # Create index file for reference
    if snp_pipeline_aligner == "bowtie2":
        target_file = reference_base_path + ".rev.1.bt2"
        needs_rebuild = utils.target_needs_rebuild([reference_file_path], target_file)
        if not args.forceFlag and not needs_rebuild:
            verbose_print("# Bowtie index %s is already freshly built.  Use the -f option to force a rebuild." % target_file)
        else:
            version_str = utils.extract_version_str("bowtie2", "bowtie2 --version")
            bowtie2_build_extra_params = os.environ.get("Bowtie2Build_ExtraParams") or ""
            command_line = "bowtie2-build " + bowtie2_build_extra_params + ' ' + reference_file_path + ' ' + reference_base_path
            verbose_print("# %s %s" % (utils.timestamp(), command_line))
            verbose_print("# %s" % version_str)
            command.run(command_line, sys.stdout)

    elif snp_pipeline_aligner == "smalt":
        target_file = reference_base_path + ".smi"
        needs_rebuild = utils.target_needs_rebuild([reference_file_path], target_file)
        if not args.forceFlag and not needs_rebuild:
            verbose_print("# Smalt index %s is already freshly built.  Use the -f option to force a rebuild." % target_file)
        else:
            version_str = utils.extract_version_str("smalt", "smalt version")
            smalt_index_extra_params = os.environ.get("SmaltIndex_ExtraParams") or ""
            command_line = "smalt index " + smalt_index_extra_params + ' ' + reference_base_path + ' ' + reference_file_path
            verbose_print("# %s %s" % (utils.timestamp(), command_line))
            verbose_print("# %s" % version_str)
            command.run(command_line, sys.stdout)

    # Create the samtools fai index
    verbose_print("")
    target_file = reference_file_path + ".fai"
    needs_rebuild = utils.target_needs_rebuild([reference_file_path], target_file)
    if not args.forceFlag and not needs_rebuild:
        verbose_print("# SAMtools fai index %s is already freshly built.  Use the -f option to force a rebuild." % target_file)
    else:
        version_str = utils.extract_version_str("samtools", "samtools 2>&1 > /dev/null")
        samtools_faidx_extra_params = os.environ.get("SamtoolsFaidx_ExtraParams") or ""
        command_line = "samtools faidx " + samtools_faidx_extra_params + ' ' + reference_file_path
        verbose_print("# %s %s" % (utils.timestamp(), command_line))
        verbose_print("# %s" % version_str)
        command.run(command_line, sys.stdout)
        utils.global_error_on_missing_file(target_file, "samtools faidx")
