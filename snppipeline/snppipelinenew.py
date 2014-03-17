#!/usr/bin/env python2.7

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import os
from multiprocessing import Pool
import utilsnew

def run_snp_pipeline(options_dict):
    """Create SNP matrix
     
    Description:
    Create a SNP matrix. This function expects or creates '(*)' the following
        files arranged in the following way:
            mainPath            
                snplist.txt
                reference.fasta
                path.txt
                snplist.txt (*)
                snpma.fasta (*)
            samples
                sample_name_one/reads.pileup (*)
                sample_name_one/var.flt.vcf
  
   The files are used as follows:
        1. Use variant file var.flt.vcf to construct SNP position list. 
        2. Use reads.pileup to extract the nucleotide base at each SNP position
            for each sample to construct the SNP fasta file.
    Note that pileups are run in parallel to speed the whole thing up.
     
    Args:
        maxThread: (15) Max number of cocurrent threads (default=15)
        mainPath:  (no default) Directory containing all input files. Output
            files will also be written here.
        Reference: (no default) File name for reference sequence (in fasta
            format) for mapping.
        pathFileName: ("path.txt") Name of file containing full paths to
            directories containing information for each sequence.
        snplistFileName: Text format list of SNP positions in samples.
        snpmaFileName: File name for snp matrix, formatted as a fasta file,
            with each sequence (all of identical length) corresponding to the
            SNPs in the correspondingly named sequence.
        bamFileNAme: #TODO - do we actually ever use this?
        pileupFileName: Name for pileup files. One is generated for each
            sample, and placed in the corresponding directory for each sample.
     
    Raises:
 
    Examples:
    args_dict = {'maxThread':2,
                 'mainPath':'',
                 'Reference':'',
                 'pathFileName':'path.txt',
                 'snplistFileName':'snplist.txt',
                 'snpmaFileName':'snpma.fa',
                 'bamFileName':'reads.bam',
                 'pileupFileName':'reads.pileup'
                }
    run_snp_pipeline(options_dict) 
    """

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
#----
#The header line names the 8 fixed, mandatory columns. These columns are as follows:
#
#    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,extra = string.split(current_line_data,maxsplit=9)
#
#If genotype data is present in the file, these are followed by a FORMAT column header, then an arbitrary number of sample IDs. The header line is tab-delimited.
#data_line = 'gi|9626243|ref|NC_001416.1|	43852	.	ACCT	A	214	.	INDEL;DP=40;VDB=0.0160;AF1=1;AC1=2;DP4=0,0,15,16;MQ=42;FQ=-128	GT:PL:GQ	1/1:255,93,0:99'
#gi|9626243|ref|NC_001416.1|	44530	.	GCAACA	GCA	214	.	INDEL;DP=48;VDB=0.0193;AF1=1;AC1=2;DP4=0,0,16,15;MQ=42;FQ=-128	GT:PL:GQ	1/1:255,93,0:99
#gi|9626243|ref|NC_001416.1|	46842	.	G	C	207	.	DP=41;VDB=0.0147;AF1=1;AC1=2;DP4=0,0,14,8;MQ=42;FQ=-90	GT:PL:GQ	1/1:240,63,0:99
#----



    snp_list_dict = dict()
    
    for sample_directory in list_of_sample_directories:
        sample_name  = sample_directory.split(os.sep)[-1]
    
        #TODO - look at use of PyVCF to process vcf file
        for line in open(sample_directory+ "/var.flt.vcf","r"):
            curVcfFileLine=line.strip()
            if curVcfFileLine.startswith("#"):
                continue
            print(curVcfFileLine)
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
                if not snp_list_dict.has_key(chrom + "\t" + pos):
                    record = [1]
                    record.append(sample_name)
                    snp_list_dict[chrom + "\t" + pos] = record
                else:
                    record = snp_list_dict[chrom + "\t" + pos]
                    record[0] += 1
                    record.append(sample_name)
        
    #==========================================================================
    #     write out list of snps for all samples to a single file.
    #==========================================================================
    snplistFile = open(options_dict['mainPath'] + options_dict['snplistFileName'], "w")
    for key in sorted(snp_list_dict.iterkeys()):
        snplistFile.write(key)
        values = snp_list_dict[key]
        for value in values:
            snplistFile.write("\t" + str(value))
        snplistFile.write("\n")
    snplistFile.close()
    
    #==========================================================================
    #   Generate Pileups of samples (in parallel)
    #==========================================================================
    
    #create a list of tuples containing values need for pileup code (as passed
    #  via pileup code wrapper)
    #TODO - maybe allow for reading of already done pileup?
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
    
    records = []

    for sample_directory in list_of_sample_directories:
        sample_name       = sample_directory.split(os.sep)[-1]
        pileup_file_name  = sample_directory + "/reads.pileup"
        positionValueHash = utilsnew.create_consensus_dict(pileup_file_name)

        seqString = ""
        for key in sorted(snp_list_dict.iterkeys()):  #TODO - Why is sorting what we want to do?
            chrom,pos   = key.split()
            if positionValueHash.has_key(chrom + ":" + pos):
                seqString += positionValueHash[chrom + ":" + pos]
            else:
                seqString += "-"

        seq = Seq(seqString)
        seqRecord = SeqRecord(seq,id=sample_name)
        records.append(seqRecord)
    
    ####write bases for snps for each sequence to a fasta file           
    fastaFile = open(options_dict['mainPath'] + options_dict['snpmaFileName'], "w") 
    SeqIO.write(records, fastaFile, "fasta")
    fastaFile.close()

#==============================================================================
# Command line driver
#==============================================================================
if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Run SNP pipeline.')
    parser.add_argument('-n','--n-processes',      dest='maxThread',type=int,default=4,help='Max number of concurrent jobs.')
    parser.add_argument('-d','--mainPath',         dest='mainPath', type=str,default='/home/yan.luo/Desktop/analysis/Montevideo/XL-C2/bowtie/Matrices/',help='Path for all files')
    parser.add_argument('-r','--Reference',        dest='Reference',type=str,default='reference.fasta',help='reference for mapping')
    parser.add_argument('-f','--pathFileName',     dest='pathFileName',type=str,default='path.txt',help='Path file name')
    parser.add_argument('-l','--snplistFileName',  dest='snplistFileName',type=str,default='snplist.txt',help='Snplist file name')
    parser.add_argument('-a','--snpmaFileName',    dest='snpmaFileName',type=str,default='snpma.fa',help='fasta file name')
    parser.add_argument('-b','--bamFileName',      dest='bamFileName',type=str,default='reads.bam',help='bam file name')
    parser.add_argument('-p','--pileupFileName',   dest='pileupFileName',type=str,default='reads.pileup',help='pileup file name')
    parser.add_argument('-v','--verbose',          dest='verbose',       type=int,default=1,help='Verbose flag (0=no info, 5=lots')
    parser.add_argument('-i','--includeReference', dest='includeReference',type=bool,default=False,help='include reference sequence in SNP matrix.')
    parser.add_argument('-o','--useOldPileups',    dest='useOldPileups',type=bool,default=False,help='Use available pileup files.')
    args = parser.parse_args()
    argsdict = vars(args)

    run_snp_pipeline(argsdict)
    

#if verbose:
#    def verboseprint(*args):
#        # Print each argument separately so caller doesn't need to
#        # stuff everything to be printed into a single string
#        for arg in args:
#           print arg,
#        print
#else:   
#    verboseprint = lambda *a: None      # do-nothing function
#
#
#If you're using Python 3, where print is already a function (or if you're willing to use print as a function in 2.x using from __future__ import print_function) it's even simpler:
#
#verboseprint = print if verbose else lambda *a, **k: None