"""This module is part of the CFSAN SNP Pipeline. It contains the code to
copy included sample data to a specified directory.
"""

from __future__ import print_function
from __future__ import absolute_import

import os
import subprocess
from pkg_resources import resource_filename


def copy_data(args):
    """Copy data included with the CFSAN SNP Pipeline to a specified directory.

    This function uses the pkg_resources package to find the data in the installation.

    The data is organized in the following way:
        snppipeline/data/agonaExpectedResults
        snppipeline/data/agonaInputs
        snppipeline/data/configuration
        snppipeline/data/lambdaVirusExpectedResults
        snppipeline/data/lambdaVirusInputs
        snppipeline/data/listeriaExpectedResults
        snppipeline/data/listeriaExpectedResults/samples
        snppipeline/data/listeriaInputs

    Parameters
    ----------
    args : argparse.Namespace
        destDirectory : Destination directory into which the SNP pipeline data files will be copied.
        whichData : Which of the supplied data sets to copy.  The choices are:
            lambdaVirusInputs          : Input reference and fastq files
            lambdaVirusExpectedResults : Expected results files
            agonaInputs                : Input reference file
            agonaExpectedResults       : Expected results files
            listeriaInputs             : Input reference file
            listeriaExpectedResults    : Expected results files
            configurationFile          : File of parameters to customize the
                                         SNP pipeline
    """
    dest_directory = args.destDirectory

    command = 'mkdir -p ' + dest_directory
    return_code = subprocess.call(command, shell=True)
    if return_code != 0:
        exit(return_code)

    data_directory = resource_filename(__name__, 'data')

    if args.whichData == 'configurationFile':
        source_file = os.path.join(data_directory, 'configuration', 'snppipeline.conf')
        command = 'cp -p ' + source_file + ' ' + dest_directory
        return_code = subprocess.call(command, shell=True)
    else:
        source_directory = os.path.join(data_directory, args.whichData)
        command = 'cp -r -p ' + source_directory + '/* ' + dest_directory
        return_code = subprocess.call(command, shell=True)

    exit(return_code)
