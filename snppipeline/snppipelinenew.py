#!/usr/bin/env python2.7
#from var.flt.vcf to construct SNP position list; from reads.pileup to extract the nucleotide base at each SNP position for each sample to construct the SNP fasta file. Multiple threads.

from Bio import SeqIO
from optparse import OptionParser
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import imp
from multiprocessing import Pool
utilsnew = imp.load_source('utilsnew', '/home/hugh.rand/projects/snppipeline/snppipeline/utilsnew.py')



#==============================================================================
#     """Create SNP matrix
#     
#     Description:
#     Create a SNP matrix based on #TODO - finish this    
#     
#     Args:
#         maxThread: Max number of cocurrent threads (default=15)
#         mainPath:  Directory containing all input files (no default). Output
#             files will also be written here.
#         Reference: File name for reference sequence (in fasta format) for
#             mapping (no default)
#         pathFileName: Name of file containing full paths to directories
#             containing information for each sequence (default="path.txt").
#         snplistFileName: Snplist file name (default="snplist.txt") #TODO: finish this one up
#         snpmaFileName: Name of file containing snp matrix in fasta format
#             (default="snpma.fa"). Written to mainPath directory
#        
#         SNP matrix
#     
#     Side effects:
#
#     Raises:
# 
#     Note:
#         (1)Each directory for each sequence is expected to have the following
#            in it:
#
#     Examples:
#         python 4snplist_matrix_P_01022014.py -n 10 -d ~/projects/snppipeline/test/testForOriginalCode/ -f path.txt -r lambda_virus.fa -l snplist.txt -a snpma.fasta
#     """
#==============================================================================

#### Command line usage
usage = "usage: %prog -n 15 -d /home/yan.luo/Desktop/ -f path.txt -r reference -l snplist.txt -a snpma.fasta -b reads.bam -p reads.pileup"

p = OptionParser(usage)
p.add_option ("-n","--cpu",dest="maxThread",type="int",default=15,help="Max count of cocurrent thread (default=15)")
p.add_option ("-d","--mainPath",dest="mainPath",default="/home/yan.luo/Desktop/analysis/Montevideo/XL-C2/bowtie/Matrices/",help="Path for all files")
p.add_option ("-r","--Reference",dest="Reference",default="CFSAN001339_pacbio.fasta",help="reference for mapping")
p.add_option ("-f","--pathFileName",dest="pathFileName",default="path.txt",help="Path file name")
p.add_option ("-l","--snplistFileName",dest="snplistFileName",default="snplist.txt",help="Snplist file name")
p.add_option ("-a","--snpmaFileName",dest="snpmaFileName",default="snpma.fa",help="fasta file name")
p.add_option ("-b","--bamFileName",dest="bamFileName",default="reads.bam",help="bam file name")
p.add_option ("-p","--pileupFileName",dest="pileupFileName",default="reads.pileup",help="pileup file name")
(opts,args)=p.parse_args()


pathFile_file_object = open(opts.mainPath + opts.pathFileName, "r")
snplistFile = open(opts.mainPath + opts.snplistFileName, "w")
snplistHash = dict()

#read in all vcf files and process into list of SNPs passing various criteria.
#  Do this for each sample. 
for pathFile_file_line in pathFile_file_object:
    filePath = pathFile_file_line[:-1]
    dirName  = filePath.split(os.sep)[-1]

    #TODO - look at use of PyVCF to process vcf file
    for line in open(filePath+ "/var.flt.vcf","r"):
        curVcfFileLine=line.strip()
        if curVcfFileLine.startswith("#"):
            continue
        print(curVcfFileLine)
        curLineData = curVcfFileLine.split()
        chrom = curLineData[0]
        pos = curLineData[1]
        info = curLineData[7]
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
#    vcfFile.close()
pathFile_file_object.close()
    
#write out list of snps for all samples to a single file        
for key in sorted(snplistHash.iterkeys()):
    snplistFile.write(key)
    values = snplistHash[key]
    for value in values:
        snplistFile.write("\t" + str(value))
    snplistFile.write("\n")
snplistFile.close()

#Generate Pileups of samples (in parallel)

#get list of sample directories to run samtools pileup in
list_of_sample_directories = [line.rstrip() for line in
                                  open(opts.mainPath + opts.pathFileName, "r")]
    
#create a list of tuples containing values need for pileup code (as passed
#  via pileup code wrapper)
parameter_list = zip(list_of_sample_directories,
                     len(list_of_sample_directories)*[opts])

#the parallel bit. Note that we use map and not map_async so that we block
#  until all the pileups are done (or bad things will happen in subsequent
#  parts of the code).
pool        = Pool(processes=opts.maxThread) # start pool
result_many = pool.map(utilsnew.pileup_wrapper, parameter_list) #parallel
#print result_many.get()

print "all commands are finished"

#Next stuff

snplistFilePath = opts.mainPath + opts.snplistFileName 
records = []
pathFile_file_object = open(opts.mainPath + opts.pathFileName, "r")
while 1:
    filePath = pathFile_file_object.readline()[:-1]
    dirName = filePath.split(os.sep)[-1]
    if not filePath:
        break

    pileupFile = filePath + "/reads.pileup"
    ###read in pileup file and store information to a dict
    positionValueHash = utilsnew.create_consensus_dict(pileupFile)

    ####append the nucleotide to the record
    snplistFile_r = open(snplistFilePath, "r")
    snplistFile_r.seek(0)
    i = 0
    seqString = ""
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
    #print "length of seqRecordString="+str(len(seqString))
    seq = Seq(seqString)
    seqRecord = SeqRecord(seq,id=dirName)
    records.append(seqRecord)
    snplistFile_r.close()

####write the records to fasta file           
fastaFile = open(opts.mainPath + opts.snpmaFileName, "w") 
SeqIO.write(records, fastaFile, "fasta")
fastaFile.close()

#==============================================================================
# 
#==============================================================================
#if __name__=='__main__':
#    parser = argparse.ArgumentParser(description='Run SNP pipeline.')
#    parser.add_argument('-n', '--n-processes',type=int,
#                        help='Max number of concurrent jobs.',default=1)
#
#    p = OptionParser(usage)
#    p.add_argument ('-n','--cpu',dest='maxThread',type='int',default=15,help='Max count of cocurrent thread (default=15)')
#    p.add_argument ('-d','--mainPath',dest='mainPath',default='/home/yan.luo/Desktop/analysis/Montevideo/XL-C2/bowtie/Matrices/',help='Path for all files')
#    p.add_argument ('-r','--Reference',dest='Reference',default='CFSAN001339_pacbio.fasta',help='reference for mapping')
#    p.add_argument ('-f','--pathFileName',dest='pathFileName',default='path.txt',help='Path file name')
#    p.add_argument ('-l','--snplistFileName',dest='snplistFileName',default='snplist.txt',help='Snplist file name')
#    p.add_argument ('-a','--snpmaFileName',dest='snpmaFileName',default='snpma.fa',help='fasta file name')
#    p.add_argument ('-b','--bamFileName',dest='bamFileName',default='reads.bam',help='bam file name')
#    p.add_argument ('-p','--pileupFileName',dest='pileupFileName',default='reads.pileup',help='pileup file name')
#    (opts,args)=p.parse_args()
#    args = parser.parse_args()
#    argsdict = vars(args)
#    run_snp_pipeline(argsdict)
#    
#def run_snp_pipeline(argsdict):
#    ...
        
