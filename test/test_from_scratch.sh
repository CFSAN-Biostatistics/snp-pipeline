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
#    [platform] = HPC or can be omitted for non-HPC
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
#   20140710-scd: Enhanced to conditionally support HPC (Torque).
#Notes:
#Bugs:
#

if (($# < 4)) || (($# > 5)); then
    echo usage: $0  cleanDirectory workDirectory compareDirectory numFastQ platform
    echo
    echo 'cleanDirectory = directory of reference and sample fastq data (lambdaVirusInputs or agonaInputs)'
    echo 'workDirectory = directory where the clean data is copied and results will be generated'
    echo 'compareDirectory = directory containing correct result files (sample pileups, snplist, and snpma.fasta)'
    echo 'numFastq = 1 or 2.  1=unpaired or 2=paired'
    echo 'platform = HPC or can be omitted for non-HPC'
    exit 1
fi

CLEANDIR=$1
WORKDIR=$2
COMPAREDIR=$3
NUMFASTQ=$4
PLATFORM=$5

if (($# < 5)); then
    PLATFORM=Regular
fi

echo -e "\nStep 1 - Prep work"
rm -rf $WORKDIR
cp -r $CLEANDIR $WORKDIR
cd $WORKDIR
ls -d --color=never samples/* > sampleDirectoryNames.txt
find samples -type f | grep fastq | sort -u > sampleFullPathNames.txt
sampleCount=$(cat sampleDirectoryNames.txt | wc -l)

echo -e "\nStep 2 - Prep the reference"
referencePath=$(ls reference/*.fasta)
referenceBasePath=${referencePath%.fasta} # strip the file extension
if [[ $PLATFORM == HPC ]]; then
    prepReferenceJobId=$(echo prepReference.sh $referenceBasePath | qsub -d $WORKDIR -N job.prepReference -j oe)
else
    prepReference.sh $referenceBasePath
fi

echo -e "\nStep 3 - Prep the samples"
if [[ $PLATFORM == HPC ]]; then
    cat sampleFullPathNames.txt | xargs -n $NUMFASTQ  > temp.prepSamples.input
    prepSamplesJobId=$(echo | qsub -t 1-$sampleCount -d $WORKDIR -N job.prepSamples -j oe -W depend=afterok:$prepReferenceJobId << _EOF_
    prepSamplesParameters=\$(cat temp.prepSamples.input | head -n \$PBS_ARRAYID | tail -n 1)
    prepSamples.sh $referenceBasePath \$prepSamplesParameters
_EOF_
)
else
    cat sampleFullPathNames.txt | xargs -n $NUMFASTQ prepSamples.sh $referenceBasePath
fi
        
echo -e "\nStep 4 - Run snp pipeline (samtools pileup in parallel and combine alignment and pileup to generate snp matrix)"
if [[ $PLATFORM == HPC ]]; then
    cmd="create_snp_matrix.py -d ./ -f sampleDirectoryNames.txt -r $referencePath -l snplist.txt -a snpma.fasta -i True"
    prepSamplesJobArray=${prepSamplesJobId%%.*}
    createSnpMatrixJobId=$(echo $cmd | qsub -d $WORKDIR -N job.createSnpMatrix -j oe -W depend=afterokarray:$prepSamplesJobArray)
else
    create_snp_matrix.py -d ./ -f sampleDirectoryNames.txt -r $referencePath -l snplist.txt -a snpma.fasta -i True
fi    

echo -e "\nStep 5 - compare results"
if [[ $PLATFORM == HPC ]]; then
    snpmaFileCount=0
    while ((snpmaFileCount == 0)); do
        qstat -a
        echo
        sleep 10
        snpmaFileCount=$(ls snpma.fasta 2>/dev/null | wc -l) # wait for the SNP matrix to exist
    done
fi
#if ! [[ $COMPAREDIR =~ ^[/~].*$ ]]; then
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
