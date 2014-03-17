#!/usr/bin/env python2.7

from __future__ import print_function
import re
import operator
import os
#import subprocess


def pileup_wrapper(args):
    "Wraps pileup to use multiple arguments with multiprocessing package."
    return pileup(*args)

def pileup(filePath,options_dict):
    """Run samtools to generate pileup.
    
    Description:
    Generate pileup files, using snplist file and the reference fasta file.
    
    Args:
        filePath: Path   #TODO - finish
        options_dict: Specified command-line options #TODO - finish
    """
    verbose = False
    verbose_print = print if verbose else lambda *a, **k: None
    
    os.chdir(filePath)
    verbose_print('Generating pileup file '+options_dict['pileupFileName']+ ' in '+filePath)
    pileup_file  = filePath + "/reads.pileup"
    snplist_file = options_dict['mainPath'] + options_dict['snplistFileName']
    
    #TODO - allow for reading of already done pileup?
    if os.path.isfile(pileup_file):
        verbose_print('Removing old pileup file '+pileup_file)
        os.remove(pileup_file)
    
    command_line = (
        'samtools mpileup ' +
            '-l ' + snplist_file +
            ' -f ' + options_dict['mainPath'] +
            options_dict['Reference'] + ' ' +
            options_dict['bamFileName'] +
            ' > ' + options_dict['pileupFileName']
    )
    verbose_print('Executing: '+command_line)
    
    return_code = os.system(command_line)  #TODO - replace with subprocess call at some point
    #subprocess.call(command_line,cwd=filePath)  #TODO - need to make this work

    if not os.path.isfile(pileup_file):
        print('Pileup file not created: '+pileup_file)
    
    verbose_print('pileup function exit')
    
    return(pileup_file)

def get_consensus_base_from_pileup(base,length,data):
    """Call the base for each SNP position
    
    Description:
    Calls the base based on the pipelup data for a SNP position with a given
        length cutoff and with a given reference base.    
    
    Args:
        base: Reference base.
        length: length cutoff.
        data: information from alignment in pileup format.
        
    Returns:
        consensus_base: Consensus base from alignment.
    
    Raises:


    Examples:
    
    >>> print(get_consensus_base_from_pileup('T',10,',.....,,.,.,...,,,.,..A'))
    T
    >>> print(get_consensus_base_from_pileup('A',10,',.$.....,,.,.,...,,,.,..^+.'))
    A
    >>> print(get_consensus_base_from_pileup('C',10,',.....,,.,.,...,,,.,..A'))
    C
    >>> print(get_consensus_base_from_pileup('T',10,',.$....,,.,.,...,,,.,...'))
    T
    >>> print(get_consensus_base_from_pileup('G',10,',$....,,.,.,...,,,.,...^l.'))
    G
    >>> print(get_consensus_base_from_pileup('G',10,'...T,,.,.,...,,,.,....'))
    G
    >>> print(get_consensus_base_from_pileup('T',10,'....,,.,.,.C.,,,.,..G.'))
    T
    >>> print(get_consensus_base_from_pileup('T',10,'....,,.,.,...,,,.,....^k.'))
    T
    >>> print(get_consensus_base_from_pileup('A',10,'A..T,,.,.,...,,,.,.....'))
    A
    """

    consensus_base = ""
    base_count_dict = {'.,': 0, 'A': 0, 'C': 0, 'G': 0, 'N': 0, 'T': 0}
    
    i = 0
    while i < len(data):
        char = data[i]
        if char == '.' or char == ',':
            base_count_dict[".,"] += 1
        elif char == 'A' or char == 'a':
            base_count_dict["A"] += 1
        elif char == 'C' or char == 'c':
            base_count_dict["C"] += 1
        elif char == 'T' or char == 't':
            base_count_dict["T"] += 1
        elif char == 'G' or char == 'g':
            base_count_dict["G"] += 1
        elif char == 'N' or char == 'n':
            base_count_dict["N"] += 1
        elif char == '+' or char == '-':
            countStr = ""
            count = 1
            while re.match("\d" ,data[i + 1]):
                countStr += data[i + 1]
                i += 1
            if countStr != "":
                count = int(countStr)
            i += count
        elif char == '^':
            if data[i+1] != "." and data[i+1] != ',':
                i +=1
        i += 1
            
    consensus_base = max(base_count_dict.iteritems(), key=operator.itemgetter(1))[0]
    if base_count_dict[consensus_base] <= (int(length)/2):
         consensus_base = "-" 
    elif consensus_base ==".,":
         consensus_base = base

    return consensus_base


def create_consensus_dict(pileup_file_path):
    """Create a dict based on the information in a pileup file.
    
    Given a path to a pileup file, create a dict based on the pileup
        data in that file. The dict contains the consensus
        base calls for positions in the reference sequence.
    
    Args:
        pileup_file_path: full path to a pileup file.
        
    Returns:
        A dict mapping a key formed from a sequence and position to a 
        consensus sequence base. For example:

        {'gi|9626243|ref|NC_001416.1|:46842': 'G',
         'gi|9626243|ref|NC_001416.1|:47425': 'T',
         'gi|9626243|ref|NC_001416.1|:47893': 'A'}

    Raises:

    Examples:
    >>> consensus_dict = create_consensus_dict("/home/hugh.rand/projects/snppipeline/test/testForOriginalCode/sample2/reads.pileup")
    >>> consensus_dict['gi|9626243|ref|NC_001416.1|:3678']
    'T'
    >>> consensus_dict['gi|9626243|ref|NC_001416.1|:40984']
    'A'
    """
    position_value_dict = dict()
    pileup_file_object  = open(pileup_file_path, "r")
    for pileup_line in pileup_file_object:
        current_line_data = pileup_line.rstrip().split()
        if len(current_line_data) >5:   #don't process lines without 5 pieces of information or more
            position_value_dict[current_line_data[0] + ":" + current_line_data[1]] = get_consensus_base_from_pileup(current_line_data[2],current_line_data[3],current_line_data[4])
    pileup_file_object.close()
    return position_value_dict
