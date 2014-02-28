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
    """Call the base for eash SNP position
    
    Calls the base based on the pipelup data for a SNP position with a given
        length cutoff and with a given reference base.    
    
    Args:
        base: Reference base.
        length: length cutoff.
        data: information from alignment in pileup format.
        
    Returns:
    
    
    Raises:


    Examples:
    >>> get_consensus_base_from_pileup('T',10,',.....,,.,.,...,,,.,..A')
    'T'
    """

    ret = ""
    charHash = dict()
    charHash[".,"] = 0
    charHash["A"] = 0
    charHash["C"] = 0
    charHash["T"] = 0
    charHash["G"] = 0
    charHash["N"] = 0

    i = 0
    while i < len(data):
        char = data[i]
        if char == '.' or char == ',':
            charHash[".,"] += 1
        elif char == 'A' or char == 'a':
            charHash["A"] += 1
        elif char == 'C' or char == 'c':
            charHash["C"] += 1
        elif char == 'T' or char == 't':
            charHash["T"] += 1
        elif char == 'G' or char == 'g':
            charHash["G"] += 1
        elif char == 'N' or char == 'n':
            charHash["N"] += 1
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
            
    ret = max(charHash.iteritems(), key=operator.itemgetter(1))[0]
    if charHash[ret] <= (int(length)/2):
         ret = "-" 
    elif ret ==".,":
         ret = base
    return ret


###store each pileup information to a Hash.
def createPositionValueHash(pileup_file_path):
    """Create a dict based on the information in a pileup file.
    
    Given a path to a pileup file, create a dictionary based on the pileup
        data in that file.    
    
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
        curLineData = pileup_line.split()
        if len(curLineData) >5:   #don't process lines without 5 pieces of information or more
            position_value_dict[curLineData[0] + ":" + curLineData[1]] = get_consensus_base_from_pileup(curLineData[2],curLineData[3],curLineData[4])
    pileup_file_object.close()
    return position_value_dict
