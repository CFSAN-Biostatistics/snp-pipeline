#!/usr/bin/env python

import os
import subprocess
import argparse
from pkg_resources import resource_filename
from snppipeline import snppipeline



#==============================================================================
# Command line driver
#==============================================================================
if __name__ == '__main__':

    data_directory = resource_filename(snppipeline.__name__, 'data')

    parser = argparse.ArgumentParser(
    description='Copy SNP Pipeline data to a specified directory.',
    formatter_class=argparse.RawTextHelpFormatter,
    epilog="""
Example:
# create a new directory "testLambdaVirus" and copy the input data there
$ copy_snppipeline_data.py lambdaVirusInputs testLambdaVirus

"""
    )

    parser.add_argument('whichData',
    metavar='whichData',
    choices=['lambdaVirusInputs', 'lambdaVirusExpectedResults',
             'agonaInputs', 'agonaExpectedResults',
             'listeriaInputs', 'listeriaExpectedResults',
             'configurationFile'],
    help="""    Which of the supplied data sets to copy.  The choices are:
        lambdaVirusInputs          : Input reference and fastq files
        lambdaVirusExpectedResults : Expected results files
        agonaInputs                : Input reference file
        agonaExpectedResults       : Expected results files
        listeriaInputs             : Input reference file
        listeriaExpectedResults    : Expected results files
        configurationFile          : File of parameters to customize the
                                     SNP pipeline

    Note: the lambda virus data set is complete with input data and expected
    results.  The agona and listeria data sets have the reference genome and
    the expected results, but not the input fastq files, because the files are
    too large to include with the package.  (default: None)
    """)

    parser.add_argument('destDirectory',
    nargs='?',
    type=str,
    default='.',
    help="""    Destination directory into which the SNP pipeline data files will be copied.
    The data files are copied into the destination directory if the directory
    already exists.  Otherwise the destination directory is created and the
    data files are copied there.  (default: current directory)
    """)

    args_dict = vars(parser.parse_args())

    dest_directory = args_dict['destDirectory']

    #print 'whichData=' + args_dict['whichData']
    #print 'dest_directory=' + args_dict['destDirectory']

    command = 'mkdir -p ' + dest_directory
    return_code = subprocess.call(command, shell=True)
    if return_code != 0:
        exit(return_code)

    if args_dict['whichData'] == 'configurationFile':
        source_file = os.path.join(data_directory, 'configuration', 'snppipeline.conf')
        command = 'cp -p ' + source_file + ' ' + dest_directory
        return_code = subprocess.call(command, shell=True)
    else:
        source_directory = os.path.join(data_directory, args_dict['whichData'])
        command = 'cp -r -p ' + source_directory + '/* ' + dest_directory
        return_code = subprocess.call(command, shell=True)

    exit(return_code)