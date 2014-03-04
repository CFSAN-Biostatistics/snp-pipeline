# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 08:35:42 2014

@author: hugh.rand
"""

###from the pileup file, call the base for each SNP position.

import re
import operator

#Used in pileup:
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#import subprocess
#import os



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


###store each pileup information to a Hash.
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
    """
    position_value_dict = dict()
    pileup_file_object  = open(pileup_file_path, "r")
    for pileup_line in pileup_file_object:
        current_line_data = pileup_line.split()
        if len(current_line_data) >5:   #don't process lines without 5 pieces of information or more
            position_value_dict[current_line_data[0] + ":" + current_line_data[1]] = get_consensus_base_from_pileup(current_line_data[2],current_line_data[3],current_line_data[4])
    pileup_file_object.close()
    return position_value_dict
