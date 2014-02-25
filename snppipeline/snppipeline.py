#!/usr/local/bin/python
#from var.flt.vcf to construct SNP position list; from reads.pileup to extract the nucleotide base at each SNP position for each sample to construct the SNP fasta file. Multiple threads.

from Bio import SeqIO
from optparse import OptionParser #TODO: Replace deprecated optparse with argparse
import sys,string,os,shutil
import re
import operator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os.path import join
from operator import itemgetter
import subprocess
from subprocess import call
from datetime import datetime
import threading
import time

#### Command line usage
usage = "usage: %prog -n 10 -d /home/yan.luo/Desktop/ -f path.txt -r reference -l snplist.txt -a snpma.fasta"

p = OptionParser(usage)
p.add_option ("-n","--cpu",dest="maxThread",type="int",default=15,help="Max count of cocurrent thread (default=15)")
p.add_option ("-d","--mainPath",dest="mainPath",default="/home/yan.luo/Desktop/analysis/Montevideo/XL-C2/bowtie/Matrices/",help="Path for all files")
p.add_option ("-r","--Reference",dest="Reference",default="CFSAN001339_pacbio.fasta",help="reference for mapping")
p.add_option ("-f","--pathFileName",dest="pathFileName",default="path.txt",help="Path file name")
p.add_option ("-l","--snplistFileName",dest="snplistFileName",default="snplist.txt",help="Snplist file name")
p.add_option ("-a","--snpmaFileName",dest="snpmaFileName",default="snpma.fa",help="fasta file name")
(opts,args)=p.parse_args()

pathFile = open(opts.mainPath + opts.pathFileName, "r")
snplistFile = open(opts.mainPath + opts.snplistFileName, "w")
snplistHash = dict()

###read all *vcf file for SNP list

while 1:
    filePath = pathFile.readline()[:-1]
    dirName = filePath.split(os.sep)[-1]
    if not filePath:
        break
    print filePath
    print dirName
    vcfFile = open(filePath + "/var.flt.vcf","r") 
    while 1:
        curVcfFileLine = vcfFile.readline()
        if not curVcfFileLine:
            break
        if curVcfFileLine.startswith("#"):
            continue

        curLineData = curVcfFileLine.split()
        chrom = curLineData[0]
        pos = curLineData[1]
        info = curLineData[7]
        if str("INDEL") in curLineData[7]:
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
    vcfFile.close()
        
for key in sorted(snplistHash.iterkeys()):
    snplistFile.write(key)
    values = snplistHash[key]
    for value in values:
        snplistFile.write("\t" + str(value))
    snplistFile.write("\n")
snplistFile.close()


pathFile.seek(0)
snplistFilePath = opts.mainPath + opts.snplistFileName 
fastaFile = open(opts.mainPath + opts.snpmaFileName, "w") 

records = []
threads = []

while 1:
    filePath = pathFile.readline()[:-1]
    dirName = filePath.split(os.sep)[-1]
    if not filePath:
        break

    t1 = FuncThread(pileup,filePath,snplistFilePath,dirName)
    threads.append(t1)
    while threading.activeCount() > opts.maxThread:
        time.sleep(15)
    t1.start()

for thread in threads:
    thread.join()

####write the records to fasta file           
SeqIO.write(records, fastaFile, "fasta")
fastaFile.close()

