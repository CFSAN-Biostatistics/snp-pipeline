"""This module is part of the CFSAN SNP Pipeline. It contains the code to
purge the intermediate output files in the "samples" directory.
"""

from __future__ import print_function
from __future__ import absolute_import

import os
import shutil

from snppipeline import utils
from snppipeline.utils import verbose_print


def purge(args):
    """Purge the intermediate output files in the "samples" directory if no errors are detected in the snp-pipeline run.

    The working directory containing the "samples" directory to recursively delete is specified by a command line argument.

    Parameters
    ----------
    args : Namespace
        dir: directory to recursively delete.

    Examples:
    args = argparse.Namespace
    args.purge_dir = "workdir"
    purge(args)
    """
    utils.print_log_header()
    utils.print_arguments(args)

    work_dir = args.work_dir

    purge_dir = os.path.join(work_dir, "samples")

    if not os.path.exists(purge_dir):
        utils.report_error(purge_dir + " does not exist.")
        return
    if not os.path.isdir(purge_dir):
        utils.report_error(purge_dir + " is not a directory.")
        return
    if not utils.is_directory_writeable(purge_dir):
        utils.report_error(purge_dir + " is not writable.")
        return

    # Don't purge if there were errors during the run
    error_output_file = os.path.join(work_dir, "error.log")
    errors_detected = os.path.isfile(error_output_file)
    if errors_detected:
        utils.report_error("Detected errors during the SNP Pipeline run. Skipping the intermediate output file purge.")
        return

    verbose_print("Purging directory %s" % purge_dir)
    shutil.rmtree(purge_dir, ignore_errors=True)
    verbose_print("Purging completed.")
