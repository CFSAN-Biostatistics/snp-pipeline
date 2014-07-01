#!/bin/bash
#
#Author: Steve Davis (scd)
#Purpose: 
#    Run set up and test of snppipline code on the agona test data set.
#Input:
#    This script takes no parameters, but it does expect the following directories 
#    and files:
#    ../snppipeline/data/agonaInputs/reference/NC_011149.fasta
#    ../snppipeline/data/agonaInputs/samples/CFSAN000448/G0H235M04.RL10.fastq
#    ../snppipeline/data/agonaInputs/samples/CFSAN000449/G00JH2D03.RL11.fastq
#    ../snppipeline/data/agonaInputs/samples/CFSAN000450/HB4DJL101.RL1.fastq
#    ../snppipeline/data/agonaInputs/samples/ERR178930/ERR178930.fastq
#    ../snppipeline/data/agonaInputs/samples/ERR178931/ERR178931.fastq
#    ../snppipeline/data/agonaExpectedResults/snplist.txt
#    ../snppipeline/data/agonaExpectedResults/snpma.fasta
#    ../snppipeline/data/agonaExpectedResults/samples/CFSAN000448/reads.pileup
#    ../snppipeline/data/agonaExpectedResults/samples/CFSAN000449/reads.pileup
#    ../snppipeline/data/agonaExpectedResults/samples/CFSAN000450/reads.pileup
#    ../snppipeline/data/agonaExpectedResults/samples/ERR178930/reads.pileup
#    ../snppipeline/data/agonaExpectedResults/samples/ERR178931/reads.pileup

#Output:
#	 This script clones the agonaInputs directory into a new testAgona
#    directory under the current working directory.  Within the testAgona directory, the input 
#    samples are copied and the outputs are generated.  Many intermediate files are
#    generated, but the most important results are:
#    are:
#        ./testAgona/snplist.txt
#            a SNP list identifying the SNPs found across all samples
#        ./testAgona/snpma.fasta
#            a SNP matrix with one row per sample and one column per SNP
#        ./testAgona/samples/CFSAN000448/reads.pileup
#        ./testAgona/samples/CFSAN000449/reads.pileup
#        ./testAgona/samples/CFSAN000450/reads.pileup
#        ./testAgona/samples/ERR178930/reads.pileup
#        ./testAgona/samples/ERR178931/reads.pileup
#            one pileup file per sample
#    The script compares the generated results to the expected results ands emits
#    differences.
#Use example:
#    cd test
#    test_from_scratch_agona_mom.sh
#History:
#   20140624-scd: Started.
#   20140701-scd: move data under the snppipeline package
#Notes:
#Bugs:
#

sh test_from_scratch.sh  ../snppipeline/data/agonaInputs  testAgona  ../snppipeline/data/agonaExpectedResults  1
