#!/usr/bin/env python2.7

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
    description='Copy SNP Pipeline example data to a specified directory.',
    formatter_class=argparse.RawTextHelpFormatter,
    epilog="""
Example:
# create a new directory "testLambdaVirus" and copy the input data there
$ copy_snppipeline_data.py lambdaVirusInputs testLambdaVirus
"""
    )
    
    parser.add_argument('whichData',
    metavar='whichData',
    choices=['lambdaVirusInputs', 'lambdaVirusExpectedResults', 'agonaInputs', 'agonaExpectedResults'],
    help="""    Which of the supplied data sets to copy.  The choices are:
        lambdaVirusInputs          : Input reference and fastq files
        lambdaVirusExpectedResults : Expected results files
        agonaInputs                : Input reference file
        agonaExpectedResults       : Expected results files
    Note: the lambda virus data set is complete with input data and expected results. 
    The agona data set has the reference genome and the expected results, but not 
    the input fastq files, because the files are too large to include with the 
    package.  (default: None)
    """)
    
    parser.add_argument('destDirectory', 
    nargs='?', 
    type=str,   
    default='.',  
    help="""    Destination directory into which the SNP pipeline example data will be copied. 
    (default: current directory)""")

    args_dict = vars(parser.parse_args())

    #print 'whichData=' + args_dict['whichData']
    #print 'destDirectory=' + args_dict['destDirectory']

    command = 'cp -r ' + data_directory + '/' + args_dict['whichData'] + ' ' + args_dict['destDirectory']

    #print command
    return_code = subprocess.call(command, shell=True)
    exit(return_code)