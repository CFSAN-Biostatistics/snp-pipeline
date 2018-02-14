"""This module is part of the CFSAN SNP Pipeline. It contains the code to
create a pileup file and find the sites with SNPs for a single sample.
"""

from __future__ import print_function
from __future__ import absolute_import

import os
import sys

from snppipeline import command
from snppipeline import utils
from snppipeline.utils import verbose_print

def call_sites(args):
    """Find the sites with SNPs in a sample.

    The sample alignment is sorted, duplicate reads are removed, a pileup is generated, and
    snps are called.

    This function expects, or creates '(*)', the following files arranged
    in the following way:
            reference
                referenceFile.fasta
            samples
                sample_name_one/reads.sorted.deduped.indelrealigned.bam
                sample_name_one/reads.all.pileup*
                sample_name_one/var.flt.vcf*

    The input files are created outside of this function. The package
    documentation provides an example of preparing these files based on the
    lambda_virus sequence that is used as one test for this package.

    Parameters
    ----------
    args : argparse.Namespace
        referenceFile : File path of the reference fasta file
        sampleDir : Relative or absolute directory of the sample
    """
    utils.print_log_header(classpath=True)
    utils.print_arguments(args)

    #==========================================================================
    # Validate inputs
    #==========================================================================

    # Verify reference fasta file exists and is not empty
    reference_file_path = args.referenceFile
    utils.verify_non_empty_input_files("Reference file", [reference_file_path], error_handler="global")

    sample_dir = args.sampleDir

    remove_duplicate_reads = os.environ.get("RemoveDuplicateReads") or "true"
    remove_duplicate_reads = remove_duplicate_reads.lower()
    if remove_duplicate_reads == "true":
        input_bam_file = os.path.join(sample_dir, "reads.sorted.deduped.indelrealigned.bam")
    else:
        input_bam_file = os.path.join(sample_dir, "reads.sorted.indelrealigned.bam")

    utils.verify_non_empty_input_files("Sample BAM file", [input_bam_file], error_handler="sample")

    sample_id = utils.sample_id_from_dir(sample_dir)

    #==========================================================================
    # Create the pileup file
    #==========================================================================

    # Check for fresh pileup; if not, create it
    pileup_file = os.path.join(sample_dir, "reads.all.pileup")
    needs_rebuild = utils.target_needs_rebuild([input_bam_file, reference_file_path], pileup_file)
    if not args.forceFlag and not needs_rebuild:
        verbose_print("# Pileup file is already freshly created for %s.  Use the -f option to force a rebuild." % sample_id)
    else:
        version_str = utils.extract_version_str("SAMtools", "samtools 2>&1 > /dev/null")
        samtools_mpileup_extra_params = os.environ.get("SamtoolsMpileup_ExtraParams") or ""
        command_line = "samtools mpileup " + samtools_mpileup_extra_params + " -f " + reference_file_path + ' ' + input_bam_file
        verbose_print("# Create pileup from bam file.")
        verbose_print("# %s %s" % (utils.timestamp(), command_line))
        verbose_print("# %s" % version_str)
        command.run(command_line, pileup_file)
        utils.sample_error_on_missing_file(pileup_file, "samtools mpileup")
        verbose_print("")

    #==========================================================================
    # Find the sites with SNPs
    #==========================================================================

    # Check for fresh unfiltered vcf; if not, create it
    vcf_file = os.path.join(sample_dir, "var.flt.vcf")
    needs_rebuild = utils.target_needs_rebuild([pileup_file], vcf_file)
    if not args.forceFlag and not needs_rebuild:
        verbose_print("# VCF file is already freshly created for %s.  Use the -f option to force a rebuild." % sample_id)
    else:
        classpath = os.environ.get("CLASSPATH")
        if not classpath or "varscan" not in classpath.lower():
            utils.global_error("Error: cannot execute VarScan. Define the path to VarScan in the CLASSPATH environment variable.")
        else:
            version_str = utils.extract_version_str("VarScan", "java net.sf.varscan.VarScan 2>&1 > /dev/null | head -n 1 | cut -d ' ' -f 2")
            varscan_jvm_extra_params = os.environ.get("VarscanJvm_ExtraParams") or ""
            varscan_mpileup2snp_extra_params = os.environ.get("VarscanMpileup2snp_ExtraParams") or ""
            command_line = "java " + varscan_jvm_extra_params + " net.sf.varscan.VarScan mpileup2snp " + pileup_file + " --output-vcf 1 " + varscan_mpileup2snp_extra_params
            verbose_print("# Create vcf file")
            verbose_print("# %s %s" % (utils.timestamp(), command_line))
            verbose_print("# %s" % version_str)
            command.run(command_line, vcf_file)
            utils.sample_error_on_missing_file(vcf_file, "VarScan")
            utils.sample_error_on_file_contains(vcf_file, "OutOfMemoryError", "VarScan")
            utils.sample_error_on_file_contains(vcf_file, "Insufficient", "VarScan")
