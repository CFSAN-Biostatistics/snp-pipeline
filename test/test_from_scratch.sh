#!/bin/bash
#
#Author: Steve Davis (scd)
#Purpose: 
#    Run set up and test of snppipline code on a specified test data set.
#Input:
#    <cleanDirectory> = directory of reference and sample fastq data (lambdaVirusInputs or agonaInputs)
#    <workDirectory> = directory where the clean data is copied and results will be generated
#    <compareDirectory> = directory containing correct result files (sample pileups, snplist, and snpma.fasta)
#    <numFastq> = 1 or 2.  1=unpaired or 2=paired
#
#    This script expects the following directories and files:
#    <cleanDirectory>/reference/<reference name>.fasta
#    <cleanDirectory>/samples/<multiple sample subdirectories>/*.fastq
#    <compareDirectory>/snplist.txt
#    <compareDirectory>/snpma.fasta
#    <compareDirectory>/samples/<multiple sample subdirectories>/reads.pileup
#Output:
#	 This script clones the <cleanDirectory> into a new <workDirectory>
#    under the current working directory.  Within the <workDirectory>, the input 
#    samples are copied and the outputs are generated.  Many intermediate files are
#    generated, but the most important results are:
#    are:
#        <workDirectory>/snplist.txt
#            a SNP list identifying the SNPs found across all samples
#        <workDirectory>/snpma.fasta
#            a SNP matrix with one row per sample and one column per SNP
#        <workDirectory>/samples/<multiple sample subdirectories>/reads.pileup
#            one pileup file per sample
#    The script compares the generated results to the expected results ands emits differences.
#Use example:
#    cd test
#    sh test_from_scratch.sh  ../snppipeline/data/lambdaVirusInputs  testLambdaVirus  ../snppipeline/data/lambdaVirusExpectedResults  2
#    sh test_from_scratch.sh  ../snppipeline/data/agonaInputs  testAgona  ../snppipeline/data/agonaExpectedResults  1
#History:
#   20140626-scd: Started.
#Notes:
#Bugs:
#

if [ $# -ne 4 ]; then
    echo usage: $0  cleanDirectory workDirectory compareDirectory numFastQ
    echo
    echo 'cleanDirectory = directory of reference and sample fastq data (lambdaVirusInputs or agonaInputs)'
    echo 'workDirectory = directory where the clean data is copied and results will be generated'
    echo 'compareDirectory = directory containing correct result files (sample pileups, snplist, and snpma.fasta)'
    echo 'numFastq = 1 or 2.  1=unpaired or 2=paired'
    exit 1
fi

CLEANDIR=$1
WORKDIR=$2
COMPAREDIR=$3
NUMFASTQ=$4

echo "\nStep 1 - Prep work"
rm -rf $WORKDIR
cp -r $CLEANDIR $WORKDIR
cd $WORKDIR
ls -d --color=never samples/* > sampleDirectoryNames.txt
find samples -type f | grep fastq | sort -u > sampleFullPathNames.txt


echo "\nStep 2 - Prep the reference"
referencePath=$(ls reference/*.fasta)
referenceBasePath=${referencePath%.fasta} # strip the file extension
prepReference.sh $referenceBasePath

echo "\nStep 3 - Prep the samples"
#  Note: This could be run in parallel using gnu parallel (workstation) or qsub (PBS on HPC)
cat sampleFullPathNames.txt | xargs -n $NUMFASTQ prepSamples.sh $referenceBasePath
        
echo "\nStep 4 - Run snp pipeline (samtools pileup in parallel and combine alignment and pileup to generate snp matrix)"
create_snp_matrix.py -n 10 -d ./ -f sampleDirectoryNames.txt -r $referencePath -l snplist.txt -a snpma.fasta -i True

echo "\nStep 5 - compare results"
cat sampleDirectoryNames.txt | while read sampleDir
do
    echo "\ndiff $sampleDir"
    diff  ../$COMPAREDIR/$sampleDir/reads.pileup $sampleDir/reads.pileup
    if [ $? -eq 0 ]; then echo OK; fi
done
echo "\ndiff snplist.txt"
diff  ../$COMPAREDIR/snplist.txt   snplist.txt
if [ $? -eq 0 ]; then echo OK; fi
echo "\ndiff snpma.fasta"
diff  ../$COMPAREDIR/snpma.fasta   snpma.fasta
if [ $? -eq 0 ]; then echo OK; fi
echo "\ndiff referenceSNP.fasta"
diff  ../$COMPAREDIR/referenceSNP.fasta   referenceSNP.fasta
if [ $? -eq 0 ]; then echo OK; fi

cd ..
