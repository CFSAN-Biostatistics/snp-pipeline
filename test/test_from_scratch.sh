#!/bin/bash
#
#Author: Steve Davis (scd)
#Purpose: 
#    Run set up and test of snppipline code on a specified test data set.
#Input:
#    <cleanDirectory>   : directory of reference and sample fastq data (lambdaVirusInputs or agonaInputs)
#    <workDirectory>    : directory where the clean data is copied and results will be generated
#    <compareDirectory> : directory containing correct result files (sample pileups, snplist, and snpma.fasta)
#    [platform]         : "torque" or can be omitted for non-HPC
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
#    ./test_from_scratch.sh  ../snppipeline/data/lambdaVirusInputs  testLambdaVirus  ../snppipeline/data/lambdaVirusExpectedResults
#    ./test_from_scratch.sh  ../snppipeline/data/agonaInputs  testAgona  ../snppipeline/data/agonaExpectedResults  torque
#History:
#   20140626-scd: Started.
#   20140710-scd: Enhanced to conditionally support HPC (Torque).
#   20140721-scd: Enhanced to support mix of paired and unpaired samples.
#Notes:
#Bugs:
#

if (($# < 3)) || (($# > 4)); then
    echo usage: $0  cleanDirectory workDirectory compareDirectory platform
    echo
    echo 'cleanDirectory   : directory of reference and sample fastq data (lambdaVirusInputs or agonaInputs)'
    echo 'workDirectory    : directory where the clean data is copied and results will be generated'
    echo 'compareDirectory : directory containing correct result files (sample pileups, snplist, and snpma.fasta)'
    echo 'platform         : "torque" or can be omitted for regular non-HPC environments'
    exit 1
fi

CLEANDIR=$1
WORKDIR=$2
COMPAREDIR=$3
PLATFORM=$4

NUMCORES=$(grep -c ^processor /proc/cpuinfo)

if (($# < 4)); then
    PLATFORM=Regular
fi

echo -e "\nStep 1 - Prep work"
mkdir -p $WORKDIR
cp -v -u -r $CLEANDIR/* $WORKDIR
find $WORKDIR -name "*.bt2" -type f -delete
find $WORKDIR -name "*.fai" -type f -delete
find $WORKDIR -name "*.pileup" -type f -delete
find $WORKDIR -name "*.sam" -type f -delete
find $WORKDIR -name "*.bam" -type f -delete
find $WORKDIR -name "*.vcf" -type f -delete
rm $WORKDIR/*
cd $WORKDIR
ls -d --color=never samples/* > sampleDirectoryNames.txt
rm sampleFullPathNames.txt
cat sampleDirectoryNames.txt | while read dir; do echo $dir/*.fastq >> sampleFullPathNames.txt; done
sampleCount=$(cat sampleDirectoryNames.txt | wc -l)

echo -e "\nStep 2 - Prep the reference"
referencePath=$(ls reference/*.fasta)
referenceBasePath=${referencePath%.fasta} # strip the file extension
if [[ $PLATFORM == torque ]]; then
    prepReferenceJobId=$(echo prepReference.sh $referenceBasePath | qsub -d $WORKDIR -N job.prepReference -j oe)
else
    prepReference.sh $referenceBasePath
fi

echo -e "\nStep 3 - Align the samples to the reference"
if [[ $PLATFORM == torque ]]; then
    numAlignThreads=8
    alignSamplesJobId=$(echo | qsub -t 1-$sampleCount -d $WORKDIR -N job.alignSamples -j oe -W depend=afterok:$prepReferenceJobId -l nodes=1:ppn=$numAlignThreads << _EOF_
    alignSamplesParameters=\$(cat sampleFullPathNames.txt | head -n \$PBS_ARRAYID | tail -n 1)
    alignSampleToReference.sh $numAlignThreads $referenceBasePath \$alignSamplesParameters
_EOF_
)
else
    cat sampleFullPathNames.txt | xargs --max-args=2 --max-lines=1 alignSampleToReference.sh $NUMCORES $referenceBasePath
fi

echo -e "\nStep 4 - Prep the samples"
if [[ $PLATFORM == torque ]]; then
    alignSamplesJobArray=${alignSamplesJobId%%.*}
    prepSamplesJobId=$(echo | qsub -t 1-$sampleCount -d $WORKDIR -N job.prepSamples -j oe -W depend=afterokarray:$alignSamplesJobArray << _EOF_
    prepSamplesParameters=\$(cat sampleDirectoryNames.txt | head -n \$PBS_ARRAYID | tail -n 1)
    prepSamples.sh $referenceBasePath \$prepSamplesParameters
_EOF_
)
else
    cat sampleDirectoryNames.txt | xargs -n 1 prepSamples.sh $referenceBasePath
fi
        
echo -e "\nStep 5 - Run snp pipeline (samtools pileup in parallel and combine alignment and pileup to generate snp matrix)"
if [[ $PLATFORM == torque ]]; then
    cmd="create_snp_matrix.py -d ./ -f sampleDirectoryNames.txt -r $referencePath -l snplist.txt -a snpma.fasta -i True"
    prepSamplesJobArray=${prepSamplesJobId%%.*}
    createSnpMatrixJobId=$(echo $cmd | qsub -d $WORKDIR -N job.createSnpMatrix -j oe -W depend=afterokarray:$prepSamplesJobArray)
else
    create_snp_matrix.py -d ./ -f sampleDirectoryNames.txt -r $referencePath -l snplist.txt -a snpma.fasta -i True
fi    

echo -e "\nStep 6 - compare results"
if [[ $PLATFORM == torque ]]; then
    snpmaFileCount=0
    while ((snpmaFileCount == 0)); do
        qstat -a
        echo
        sleep 10
        snpmaFileCount=$(ls snpma.fasta 2>/dev/null | wc -l) # wait for the SNP matrix to exist
    done
fi
if [[ $COMPAREDIR != /* ]] && [[ $COMPAREDIR != ~* ]]; then
    COMPAREDIR=../$COMPAREDIR # handle relative path
fi
cat sampleDirectoryNames.txt | while read sampleDir
do
    echo -e "\ndiff $sampleDir/reads/pileup"
    diff  $COMPAREDIR/$sampleDir/reads.pileup $sampleDir/reads.pileup
    if [ $? -eq 0 ]; then echo OK; fi
done
echo -e "\ndiff snplist.txt"
diff  $COMPAREDIR/snplist.txt   snplist.txt
if [ $? -eq 0 ]; then echo OK; fi
echo -e "\ndiff snpma.fasta"
diff  $COMPAREDIR/snpma.fasta   snpma.fasta
if [ $? -eq 0 ]; then echo OK; fi
echo -e "\ndiff referenceSNP.fasta"
diff  $COMPAREDIR/referenceSNP.fasta   referenceSNP.fasta
if [ $? -eq 0 ]; then echo OK; fi

cd - > /dev/null
exit 0
