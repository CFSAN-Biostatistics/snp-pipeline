#!/usr/bin/env python2.7

from Bio import SeqIO
import argparse
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
import utilsnew

#==============================================================================
#     """Create SNP matrix
#     
#     Description:
#     Create a SNP matrix. Use variant file var.flt.vcf to construct SNP
#         position list. Use reads.pileup to extract the nucleotide base at
#         each SNP position for each sample to construct the SNP fasta file.
#         Pileups are run in parallel to speed the whole thing up.
#     
#     Args:
#         maxThread: Max number of cocurrent threads (default=15)
#         mainPath:  Directory containing all input files (no default). Output
#             files will also be written here.
#         Reference: File name for reference sequence (in fasta format) for
#             mapping (no default)
#         pathFileName: Name of file containing full paths to directories
#             containing information for each sequence (default="path.txt").
#         snplistFileName: Snplist file name (default="snplist.txt") #TODO - finish this one up
#         snpmaFileName: Name of file containing snp matrix in fasta format
#             (default="snpma.fa"). Written to mainPath directory
#        
#         SNP matrix #TODO - finish this one
#     
#     Side effects:
#
#     Raises:
# 
#     Note:
#         (1)Each directory for each sequence is expected to have the following
#            in it:  #TODO - finish it
#
#     Examples:
#         python 4snplist_matrix_P_01022014.py -n 10 -d ~/projects/snppipeline/test/testForOriginalCode/ -f path.txt -r lambda_virus.fa -l snplist.txt -a snpma.fasta
# args_dict = {'maxThread':2,
#             'mainPath':'',
#             'Reference':'',
#             'pathFileName':'path.txt',
#             'snplistFileName':'snplist.txt',
#             'snpmaFileName':'snpma.fa',
#             'bamFileName':'reads.bam',
#             'pileupFileName':'reads.pileup'
#            }
# run_snp_pipeline(options_dict) 
#    """
#==============================================================================

def run_snp_pipeline(options_dict):
    #==========================================================================
    #Prep work     
    #==========================================================================

    list_of_sample_directories = [line.rstrip() for line in
                                      open(options_dict['mainPath'] + options_dict['pathFileName'], "r")]
    #remove any blank rows that were read in
    list_of_sample_directories = filter(None,list_of_sample_directories)

    #==========================================================================
    #read in all vcf files and process into list of SNPs passing various
    #  criteria. Do this for each sample.
    #==========================================================================
    snplistHash = dict()
    
    print(options_dict['mainPath'] + options_dict['pathFileName'])
    for filePath in list_of_sample_directories:
        dirName  = filePath.split(os.sep)[-1]
    
        #TODO - look at use of PyVCF to process vcf file
        for line in open(filePath+ "/var.flt.vcf","r"):
            curVcfFileLine=line.strip()
            if curVcfFileLine.startswith("#"):
                continue
            curLineData = curVcfFileLine.split()
            chrom = curLineData[0]
            pos   = curLineData[1]
            info  = curLineData[7]
            if str("INDEL") in info:
                continue
            infoFields = info.split(";")
            dpFlag = False
            af1Flag = False
            for infoField in infoFields:
                infoPair = infoField.split("=")
                if infoPair[0] == "DP" and int(infoPair[1]) >= 10:
                    dpFlag = True
                elif (infoPair[0] == "AF1" and infoPair[1] == "1" ) or (infoPair[0] == "AR" and infoPair[1] == "1.00"):
                    af1Flag = True
            # find a good record fo SNP position, save data to hash
            if dpFlag and af1Flag:
                if not snplistHash.has_key(chrom + "\t" + pos):
                    record = [1]
                    record.append(dirName)
                    snplistHash[chrom + "\t" + pos] = record
                else:
                    record = snplistHash[chrom + "\t" + pos]
                    record[0] += 1
                    record.append(dirName)
    #    vcfFile.close()    #TODO - did this get closed or not?
        
    #write out list of snps for all samples to a single file        
    snplistFile = open(options_dict['mainPath'] + options_dict['snplistFileName'], "w")
    for key in sorted(snplistHash.iterkeys()):
        snplistFile.write(key)
        values = snplistHash[key]
        for value in values:
            snplistFile.write("\t" + str(value))
        snplistFile.write("\n")
    snplistFile.close()
    
    #==========================================================================
    #   Generate Pileups of samples (in parallel)
    #==========================================================================
    
    #create a list of tuples containing values need for pileup code (as passed
    #  via pileup code wrapper)
    parameter_list = zip(list_of_sample_directories,
                         len(list_of_sample_directories)*[options_dict])
    
    #the parallel bit. Note that we use map and not map_async so that we block
    #  until all the pileups are done (or bad things will happen in subsequent
    #  parts of the code).
    pool        = Pool(processes=options_dict['maxThread']) # start pool
    result_many = pool.map(utilsnew.pileup_wrapper, parameter_list) #parallel
    #print result_many.get()
    
    print "all pileups are finished"
    
    #==========================================================================
    #   Create snp matrix
    #==========================================================================
    
    snplistFilePath = options_dict['mainPath'] + options_dict['snplistFileName'] 
    records = []

    for filePath in list_of_sample_directories:
        dirName    = filePath.split(os.sep)[-1]
        pileupFile = filePath + "/reads.pileup"

        ###read in pileup file and store information to a dict
        positionValueHash = utilsnew.create_consensus_dict(pileupFile)

        ####append the nucleotide to the record
        seqString = ""
        with open(snplistFilePath,'r') as snplist_file_object:
            for curSnplistLine in snplist_file_object:
                curSnplistData = curSnplistLine.split()
                if  len(curSnplistData) <2:
                    print('snplistfile: bad line: '+curSnplistLine)
                    continue
                chrom = curSnplistData[0]
                pos   = curSnplistData[1]
        
                if positionValueHash.has_key(chrom + ":" + pos):
                    seqString += positionValueHash[chrom + ":" + pos]
                else:
                    seqString += "-"
        seq = Seq(seqString)
        seqRecord = SeqRecord(seq,id=dirName)
        records.append(seqRecord)
    
    ####write the records to fasta file           
    fastaFile = open(options_dict['mainPath'] + options_dict['snpmaFileName'], "w") 
    SeqIO.write(records, fastaFile, "fasta")
    fastaFile.close()

#==============================================================================
# Command line driver
#==============================================================================
if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Run SNP pipeline.')
    parser.add_argument('-n','--n-processes',dest='maxThread',type=int,default=4,help='Max number of concurrent jobs.')
    parser.add_argument('-d','--mainPath',dest='mainPath',type=str,default='/home/yan.luo/Desktop/analysis/Montevideo/XL-C2/bowtie/Matrices/',help='Path for all files')
    parser.add_argument('-r','--Reference',dest='Reference',type=str,default='CFSAN001339_pacbio.fasta',help='reference for mapping')
    parser.add_argument('-f','--pathFileName',dest='pathFileName',type=str,default='path.txt',help='Path file name')
    parser.add_argument('-l','--snplistFileName',dest='snplistFileName',type=str,default='snplist.txt',help='Snplist file name')
    parser.add_argument('-a','--snpmaFileName',dest='snpmaFileName',type=str,default='snpma.fa',help='fasta file name')
    parser.add_argument('-b','--bamFileName',dest='bamFileName',type=str,default='reads.bam',help='bam file name')
    parser.add_argument('-p','--pileupFileName',dest='pileupFileName',type=str,default='reads.pileup',help='pileup file name')
    args = parser.parse_args()
    argsdict = vars(args)

    run_snp_pipeline(argsdict)
    
