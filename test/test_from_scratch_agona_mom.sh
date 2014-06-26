#!/bin/bash
#
#Author: Steve Davis (scd)
#Purpose: 
#    Run set up and test of snppipline code on the lamdba virus test data set.
#Input:
#    This script takes no parameters, but it does expect the following directories 
#    and files:
#    ./testAgonaClean/reference/NC_011149.fasta
#    ./testAgonaClean/samples/CFSAN000448/G0H235M04.RL10.fastq
#    ./testAgonaClean/samples/CFSAN000449/G00JH2D03.RL11.fastq
#    ./testAgonaClean/samples/CFSAN000450/HB4DJL101.RL1.fastq
#    ./testAgonaClean/samples/ERR178930/ERR178930.fastq
#    ./testAgonaClean/samples/ERR178931/ERR178931.fastq
#	 ./codeComparisonFiles/testAgonaMOM/snplist.txt
#	 ./codeComparisonFiles/testAgonaMOM/snpma.fasta
#    ./codeComparisonFiles/testAgonaMOM/samples/CFSAN000448/reads.pileup
#    ./codeComparisonFiles/testAgonaMOM/samples/CFSAN000449/reads.pileup
#    ./codeComparisonFiles/testAgonaMOM/samples/CFSAN000450/reads.pileup
#    ./codeComparisonFiles/testAgonaMOM/samples/ERR178930/reads.pileup
#    ./codeComparisonFiles/testAgonaMOM/samples/ERR178931/reads.pileup

#Output:
#	 This script clones the testAgonaClean directory into a new testAgonaMOM
#    directory under the current working directory.  Within the testAgonaMOM directory, the input 
#    samples are copied and the outputs are generated.  Many intermediate files are
#    generated, but the most important results are:
#    are:
#        ./testAgonaMOM/snplist.txt
#            a SNP list identifying the SNPs found across all samples
#        ./testAgonaMOM/snpma.fasta
#            a SNP matrix with one row per sample and one column per SNP
#        ./testAgonaMOM/samples/CFSAN000448/reads.pileup
#        ./testAgonaMOM/samples/CFSAN000449/reads.pileup
#        ./testAgonaMOM/samples/CFSAN000450/reads.pileup
#        ./testAgonaMOM/samples/ERR178930/reads.pileup
#        ./testAgonaMOM/samples/ERR178931/reads.pileup
#            one pileup file per sample
#    The script compares the generated results to the expected results ands emits
#    differences.
#Use example:
#    cd test
#    test_from_scratch_agona_mom.sh
#History:
#   20140624-scd: Started.
#Notes:
#Bugs:
#

sh test_from_scratch.sh  testAgonaClean  testAgonaMOM  codeComparisonFiles/testAgonaMOM  1
