#!/bin/bash
#
#Author: Hugh A. Rand (har)
#Purpose: Run set up and test of snppipline code.
#Input:
#    baseDirectory
#Output:
#    various
#Use example:
#    runLamdaVirusTest.sh /home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline
#History:
#   20140610-har: Started.
#Notes:
#Bugs:
#

#Process arguments
if [ -z "$1" ]; then
    echo usage: $0 baseDirectory
    exit
fi
BASEDIRECTORY=$1

#Prep work
#  Copy the sequence data to work on to a test directory.
cd $BASEDIRECTORY'/test'
cp -r testLambdaVirusClean testLambdaVirus
cd testLambdaVirus
ls -1 --color=never $PWD/samples > sampleDirectoryNames.txt

#Set up reference sequence
$BASEDIRECTORY'/scripts/prepReference.sh' lambda_virus

#Set up sample sequence
#  Note: This could be run in parallel using gnu parallel (workstation) or
#    qsub (PBS on HPC)
cat sampleDirectoryNames.txt | xargs -n 1 $BASEDIRECTORY'/scripts/prepSamples.sh' lambda_virus
        
#Run snp pipeline (samtools pileup in parallel and combine alignment and pileup to
#   generate snp matrix)
#TODO replace 'path.txt' with better name
ls -d -1 --color=never $PWD/samples/* > path.txt
#$BASEDIRECTORY'/scripts/runsnppipeline.py' -n 10 -d ~/mnt/biob/svn/Biostats/rand/snppipeline/test/testLambdaVirus/ -f path.txt -r reference/lambda_virus.fasta -l snplist.txt -a snpma.fasta -i True
$BASEDIRECTORY/test/test_snppipeline.py -v




