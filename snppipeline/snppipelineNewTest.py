#!/usr/local/bin/python
# -*- coding: utf-8 -*-

import sys, argparse

def run_snp_pipeline(argv):
    """
    run_snp_pipeline(-n 10 -d /home/yan.luo/Desktop/ -f path.txt -r reference -l snplist.txt -a snpma.fasta')
    
    """
    #### Command line usage
    usage = "usage: %prog -n 10 -d /home/yan.luo/Desktop/ -f path.txt -r reference -l snplist.txt -a snpma.fasta"
    
    p = argparse.ArgumentParser(description='Run snp pipeline.')
    p.add_argument("-n","--cpu",dest="maxThread",type="int",default=15,help="Max count of cocurrent thread (default=15)")
    p.add_argument("-d","--mainPath",dest="mainPath",default="/home/yan.luo/Desktop/analysis/Montevideo/XL-C2/bowtie/Matrices/",help="Path for all files")
    p.add_argument("-r","--Reference",dest="Reference",default="CFSAN001339_pacbio.fasta",help="reference for mapping")
    p.add_argument("-f","--pathFileName",dest="pathFileName",default="path.txt",help="Path file name")
    p.add_argument("-l","--snplistFileName",dest="snplistFileName",default="snplist.txt",help="Snplist file name")
    p.add_argument("-a","--snpmaFileName",dest="snpmaFileName",default="snpma.fa",help="fasta file name")
    (opts,args)=p.parse_args()
    
    #pathFile = open(opts.mainPath + opts.pathFileName, "r") #TODO get working with error checking
    #snplistFile = open(opts.mainPath + opts.snplistFileName, "w") #TODO get working with error checking
    snplistHash = dict()
    
    ###read all *vcf file for SNP list
    
    ##Do Stuff
    
    ####write the records to fasta file           

#if __name__ == "__main__":
#    run_snp_pipeline(sys.argv)


print('test')