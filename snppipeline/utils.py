#!/usr/local/bin/python
# -*- coding: utf-8 -*-

"""
Utility classes and functions for snp pipeline code.

@author: hugh.rand
"""

import re
import operator

class FuncThread(threading.Thread):
    def __init__(self,target,*args):
        self._target = target
        self._args = args
        threading.Thread.__init__(self)

    def run(self):
        self._target(*self._args)

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


def create_position_value_hash(pileupFilePath):
    """Store each pileup information to a Hash.
    """"
    positionValueHash = dict()
    pileupFile = open(pileupFilePath, "r")
    while 1:
        curpileupFileLine = pileupFile.readline()
        if not curpileupFileLine:
            break
        curLineData = curpileupFileLine.split()
        if len(curLineData) <5:
            continue
        positionValueHash[curLineData[0] + ":" + curLineData[1]] = get_consensus_base_from_pileup(curLineData[2],curLineData[3],curLineData[4])
    pileupFile.close()
    return positionValueHash


def pileup(filePath,snplistFilePath,dirName):
    """
    
    Retrieves rows pertaining to the given keys from the Table instance
    represented by big_table.  Silly things may happen if
    other_silly_variable is not None.

    Args:
        file_path: Directory in which to save pileup files.
        snp_list_file_path: A sequence of strings representing the key of each table row
            to fetch.
        dir_name: 

    Returns:
 
    Raises:

    """
    seqString = ""
    os.chdir(filePath)

    ####generate pileup files, using snplist file and the reference fasta file.
    subprocess.call("samtools mpileup -l " + opts.mainPath + opts.snplistFileName +
                    " -f " + opts.mainPath + opts.Reference +
                    " reads.bam > reads.pileup",shell=True )

    ####read in pileup file and store information to a hash
    positionValueHash = create_position_value_hash(filePath + "/reads.pileup")

    ####append the nucleotide to the record
    snplistFile_r = open(snplistFilePath, "r")
    snplistFile_r.seek(0)
    i = 0
    while 1:
        curSnplistLine = snplistFile_r.readline()
        if not curSnplistLine:
            break
        i = i+1
        curSnplistData = curSnplistLine.split()
        if  len(curSnplistData) <2:
            print "snplistfile: bad line# "+i+" line="+curSnplistLine
            continue
        chrom = curSnplistData[0]
        pos = curSnplistData[1]

        if positionValueHash.has_key(chrom + ":" + pos):
            seqString += positionValueHash[chrom + ":" + pos]
        else:
            seqString += "-"
    print "length of seqRecordString="+str(len(seqString))
    seq = Seq(seqString)
    seqRecord = SeqRecord(seq,id=dirName)
    records.append(seqRecord)
    snplistFile_r.close()
