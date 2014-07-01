#!/bin/bash
#
#Author: Steve Davis (scd)
#Purpose: 
#    Run set up and test of snppipline code on the lamdba virus test data set.
#Input:
#    This script takes no parameters, but it does expect the following directories 
#    and files:
#    ../snppipeline/data/lambdaVirusInputs/reference/lambda_virus.fasta
#    ../snppipeline/data/lambdaVirusInputs/samples/sample1/sample1_1.fastq
#    ../snppipeline/data/lambdaVirusInputs/samples/sample1/sample1_2.fastq
#    ../snppipeline/data/lambdaVirusInputs/samples/sample2/sample2_1.fastq
#    ../snppipeline/data/lambdaVirusInputs/samples/sample2/sample2_2.fastq
#    ../snppipeline/data/lambdaVirusInputs/samples/sample3/sample3_1.fastq
#    ../snppipeline/data/lambdaVirusInputs/samples/sample3/sample3_2.fastq
#    ../snppipeline/data/lambdaVirusInputs/samples/sample4/sample4_1.fastq
#    ../snppipeline/data/lambdaVirusInputs/samples/sample4/sample4_2.fastq
#    ../snppipeline/data/lambdaVirusExpectedResults/snplist.txt
#    ../snppipeline/data/lambdaVirusExpectedResults/snpma.fasta
#    ../snppipeline/data/lambdaVirusExpectedResults/samples/sample1/reads.pileup
#    ../snppipeline/data/lambdaVirusExpectedResults/samples/sample2/reads.pileup
#    ../snppipeline/data/lambdaVirusExpectedResults/samples/sample3/reads.pileup
#    ../snppipeline/data/lambdaVirusExpectedResults/samples/sample4/reads.pileup
#Output:
#	 This script clones the lambdaVirusInputs directory into a new testLambdaVirus
#    directory under the current working directory.  Within the testLambdaVirus directory, the input 
#    samples are copied and the outputs are generated.  Many intermediate files are
#    generated, but the most important results are:
#    are:
#        ./testLambdaVirus/snplist.txt
#            a SNP list identifying the SNPs found across all samples
#        ./testLambdaVirus/snpma.fasta
#            a SNP matrix with one row per sample and one column per SNP
#        ./testLambdaVirus/samples/sample1/reads.pileup
#        ./testLambdaVirus/samples/sample2/reads.pileup
#        ./testLambdaVirus/samples/sample3/reads.pileup
#        ./testLambdaVirus/samples/sample4/reads.pileup
#            one pileup file per sample
#    The script compares the generated results to the expected results ands emits
#    differences.
#Use example:
#    cd test
#    test_from_scratch_lambda_virus.sh
#History:
#   20140624-scd: Started.
#   20140625-scd: changed to a one-liner calling test_from_scratch.sh
#   20140701-scd: move data under the snppipelien package
#Notes:
#Bugs:
#

sh test_from_scratch.sh  ../snppipeline/data/lambdaVirusInputs  testLambdaVirus  ../snppipeline/data/lambdaVirusExpectedResults  2
