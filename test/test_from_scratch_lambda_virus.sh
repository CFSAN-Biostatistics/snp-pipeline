#!/bin/bash

echo "\nStep 1 - Prep work"
rm -rf testLambdaVirus
cp -r testLambdaVirusClean testLambdaVirus
cd testLambdaVirus
ls -1 --color=never $PWD/samples > sampleDirectoryNames.txt

echo "\nStep 2 - Prep the reference"
prepReference.sh lambda_virus

echo "\nStep 3 - Prep the samples"
#  Note: This could be run in parallel using gnu parallel (workstation) or qsub (PBS on HPC)
cat sampleDirectoryNames.txt | xargs -n 1 prepSamples.sh lambda_virus
        
echo "\nStep 4 - Run snp pipeline (samtools pileup in parallel and combine alignment and pileup to generate snp matrix)"
#TODO replace 'path.txt' with better name
ls -d -1 --color=never $PWD/samples/* > path.txt
runsnppipeline.py -n 10 -d ./ -f path.txt -r reference/lambda_virus.fasta -l snplist.txt -a snpma.fasta -i True

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
