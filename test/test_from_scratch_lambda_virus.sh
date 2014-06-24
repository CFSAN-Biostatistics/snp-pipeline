#!/bin/bash
#
#Author: Steve Davis (scd)
#Purpose: 
#    Run set up and test of snppipline code on the lamdba virus test data set.
#Input:
#    This script takes no parameters, but it does expect the following directories 
#    and files:
#    ./testLambdaVirusClean/reference/lambda_virus.fasta
#    ./testLambdaVirusClean/samples/sample1/sample1_1.fastq
#    ./testLambdaVirusClean/samples/sample1/sample1_2.fastq
#    ./testLambdaVirusClean/samples/sample2/sample2_1.fastq
#    ./testLambdaVirusClean/samples/sample2/sample2_2.fastq
#    ./testLambdaVirusClean/samples/sample1/sample3_1.fastq
#    ./testLambdaVirusClean/samples/sample1/sample3_2.fastq
#    ./testLambdaVirusClean/samples/sample1/sample4_1.fastq
#    ./testLambdaVirusClean/samples/sample1/sample4_2.fastq
#	 ./codeComparisonFiles/testLambdaVirus/snplist.txt
#	 ./codeComparisonFiles/testLambdaVirus/snpma.fasta
#    ./codeComparisonFiles/testLambdaVirus/samples/sample1/reads.pileup
#    ./codeComparisonFiles/testLambdaVirus/samples/sample2/reads.pileup
#    ./codeComparisonFiles/testLambdaVirus/samples/sample3/reads.pileup
#    ./codeComparisonFiles/testLambdaVirus/samples/sample4/reads.pileup
#Output:
#	 This script clones the testLambdaVirusClean directory into a new testLambdaVirus
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
#Notes:
#Bugs:
#


echo "\nStep 1 - Prep work"
rm -rf testLambdaVirus
cp -r testLambdaVirusClean testLambdaVirus
cd testLambdaVirus
ls -d --color=never samples/* > sampleDirectoryNames.txt

echo "\nStep 2 - Prep the reference"
prepReference.sh reference/lambda_virus

echo "\nStep 3 - Prep the samples"
#  Note: This could be run in parallel using gnu parallel (workstation) or qsub (PBS on HPC)
cat sampleDirectoryNames.txt | xargs -n 1 prepSamples.sh reference/lambda_virus
        
echo "\nStep 4 - Run snp pipeline (samtools pileup in parallel and combine alignment and pileup to generate snp matrix)"
runsnppipeline.py -n 10 -d ./ -f sampleDirectoryNames.txt -r reference/lambda_virus.fasta -l snplist.txt -a snpma.fasta -i True

echo "\nStep 5 - compare results"
echo "\ndiff pileup1"
diff  ../codeComparisonFiles/testLambdaVirus/samples/sample1/reads.pileup    samples/sample1/reads.pileup
if [ $? -eq 0 ]; then echo OK; fi
echo "\ndiff pileup2"
diff  ../codeComparisonFiles/testLambdaVirus/samples/sample2/reads.pileup    samples/sample2/reads.pileup
if [ $? -eq 0 ]; then echo OK; fi
echo "\ndiff pileup3"
diff  ../codeComparisonFiles/testLambdaVirus/samples/sample3/reads.pileup    samples/sample3/reads.pileup
if [ $? -eq 0 ]; then echo OK; fi
echo "\ndiff pileup4"
diff  ../codeComparisonFiles/testLambdaVirus/samples/sample4/reads.pileup    samples/sample4/reads.pileup
if [ $? -eq 0 ]; then echo OK; fi
echo "\ndiff snplist.txt"
diff  ../codeComparisonFiles/testLambdaVirus/snplist.txt   snplist.txt
if [ $? -eq 0 ]; then echo OK; fi
echo "\ndiff snpma.fasta"
diff  ../codeComparisonFiles/testLambdaVirus/snpma.fasta   snpma.fasta
if [ $? -eq 0 ]; then echo OK; fi

cd ..
