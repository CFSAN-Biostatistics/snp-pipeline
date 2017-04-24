"""This module is part of the CFSAN SNP Pipeline. It contains the code to
combine the metrics from all samples into a single table of metrics for all samples.
"""

from __future__ import print_function
from __future__ import absolute_import

import os

from snppipeline import utils
from snppipeline.utils import verbose_print


def combine_metrics(args):
    """Combine the per-sample metrics files into a single table of metrics for all samples.

    This function expects, or creates '(*)', the following files arranged
    in the following way:
            samples
                sample_name_one/metrics
            metrics.tsv

    All the input files are created outside of this function.  Before
    running this command, the metrics file for each sample must be created by the
    collectSampleMetrics script.

    The package documentation provides an example of preparing these files based
    on the lambda_virus sequence that is used as one test for this package.

    Parameters
    ----------
    args : argparse.Namespace
        sampleDirsFile : Path to file containing a list of directories -- one per sample
        metricsFileName : File name of the metrics files which must exist in each of the sample directories
        mergedMetricsFile : Path to the output merged metrics file
    """
    utils.print_log_header()
    utils.print_arguments(args)

    #==========================================================================
    # Validate inputs
    #==========================================================================

    sample_directories_list_path = args.sampleDirsFile
    utils.verify_non_empty_input_files("File of sample directories", [sample_directories_list_path], error_handler="global")

    metrics_file_name = args.metricsFileName
    merged_metrics_path = args.mergedMetricsFile

    with open(sample_directories_list_path, "r") as f:
        sample_directories = [line.rstrip() for line in f]
    sample_directories = [d for d in sample_directories if d]
    metrics_files = [os.path.join(d, metrics_file_name) for d in sample_directories]

    #==========================================================================
    # Check if merge has already been done
    #==========================================================================
    needs_rebuild = utils.target_needs_rebuild(metrics_files, merged_metrics_path)
    if not args.forceFlag and not needs_rebuild:
        verbose_print("# The merged metrics file is already freshly created.  Use the -f option to force a rebuild.")
        return

    #==========================================================================
    # Parse the metrics files and print the tabular results
    #==========================================================================

    with open(merged_metrics_path, 'w') as f:
        # Emit the column headings
        column_headings = ["Sample", "Fastq Files", "Fastq File Size", "Machine", "Flowcell", "Number of Reads", "Duplicate Reads", "Percent of Reads Mapped",
                           "Average Insert Size", "Average Pileup Depth", "Phase1 SNPs", "Phase1 Preserved SNPs", "Phase2 SNPs", "Phase2 Preserved SNPs",
                           "Missing SNP Matrix Positions", "Missing Preserved SNP Matrix Positions", "Excluded Sample", "Excluded Preserved Sample", "Warnings and Errors"]
        if not args.spaceHeadings:
            column_headings = [heading.replace(' ', '_') for heading in column_headings]

        tabbed_headings = '\t'.join(column_headings)
        f.write(tabbed_headings + '\n')

        # Reads the metrics from each sample, and emit the values
        for metrics_file in metrics_files:
            verbose_print("Processing " + metrics_file)
            message = None
            if not os.path.isfile(metrics_file):
                message = "Sample metrics file %s does not exist." % metrics_file
            elif os.path.getsize(metrics_file) == 0:
                message = "Sample metrics file %s is empty." % metrics_file
            if message:
                f.write(message + '\n')
                utils.sample_warning(message)
                continue

            metrics = utils.read_properties(metrics_file)

            f.write(quoted(metrics.get("sample", "")) + '\t')
            f.write(quoted(metrics.get("fastqFileList", "")) + '\t')
            f.write(metrics.get("fastqFileSize", "") + '\t')
            f.write(metrics.get("machine", "") + '\t')
            f.write(metrics.get("flowcell", "") + '\t')
            f.write(metrics.get("numberReads", "") + '\t')
            f.write(metrics.get("numberDupReads", "") + '\t')
            f.write(metrics.get("percentReadsMapped", "") + '\t')
            f.write(metrics.get("aveInsertSize", "") + '\t')
            f.write(metrics.get("avePileupDepth", "") + '\t')
            f.write(metrics.get("phase1Snps", "") + '\t')
            f.write(metrics.get("phase1SnpsPreserved", "") + '\t')
            f.write(metrics.get("snps", "") + '\t')
            f.write(metrics.get("snpsPreserved", "") + '\t')
            f.write(metrics.get("missingPos", "") + '\t')
            f.write(metrics.get("missingPosPreserved", "") + '\t')
            f.write(metrics.get("excludedSample", "") + '\t')
            f.write(metrics.get("excludedSamplePreserved", "") + '\t')
            f.write(quoted(metrics.get("errorList", "")) + '\n')


def quoted(text):
    """Return a string surrounded by quotes if the string is not empty.

    Parameters
    ----------
    text : str
        A text string.

    Returns
    -------
    quoted_text : str
        Text string surrounded by quotes if the string is not empty.
    """
    if not text:
        return text
    return '"' + text + '"'
