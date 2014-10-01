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
#    This script clones the <cleanDirectory> into a new <workDirectory>
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
#   20140728-scd: Changed the comparison logic include the upstream files to detect reproducibility problems sooner.
#   20140728-scd: Changed the comparison logic to easily allow running a test repeatedly overnight to verify repeatable results.
#   20140818-scd: Split the create_snp_matrix script into 4 smaller scripts.
#   20140919-scd: Changed to use the new script, run_snp_pipeline.sh
#Notes:
#Bugs:
#

if (($# < 4)) || (($# > 5)); then
    echo usage: $0  referenceFile samplesDirectory workDirectory compareDirectory platform
    echo
    echo 'referenceFile    : Relative or absolute path to the reference fasta file.'
    echo 'samplesDirectory : Relative or absolute path to the parent directory of all the sample directories.'
    echo 'workDirectory    : directory where the clean data is copied and results will be generated'
    echo 'compareDirectory : directory containing correct result files (sample pileups, snplist, and snpma.fasta)'
    echo 'platform         : "torque" or can be omitted for regular non-HPC environments'
    exit 1
fi

referenceFile="$1"
samplesDirectory="$2"
WORKDIR="$3"
COMPAREDIR="$4"
PLATFORM="$5"

echo referenceFile=$referenceFile
echo samplesDirectory=$samplesDirectory
echo WORKDIR=$WORKDIR
echo COMPAREDIR=$COMPAREDIR

run_snp_pipeline.sh -m -o "$WORKDIR" -s "$samplesDirectory" "$referenceFile"

#cd $WORKDIR


echo -e "\nStep 9 - compare results"
if [[ "$PLATFORM" == "torque" ]]; then
    snpmaFileCount=0
    while ((snpmaFileCount == 0)); do
        qstat -a
        echo
        sleep 10
        snpmaFileCount=$(ls snpma.fasta 2>/dev/null | wc -l) # wait for the SNP matrix to exist
    done
fi

echo -e "\ndiff reads.sam"
echo -n "reads.sam:" >> overnight_run.txt
cat $WORKDIR/sampleDirectories.txt | while read sampleDir
do
    partialDir=${sampleDir##*/samples/} # strip directories
    diff -q --ignore-matching-lines=bowtie $COMPAREDIR/samples/$partialDir/reads.sam $sampleDir/reads.sam
    stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
    echo -n $stat" " >> overnight_run.txt
done

echo -e "\ndiff reads.all.pileup"
echo -n "  reads.all.pileup:" >> overnight_run.txt
cat $WORKDIR/sampleDirectories.txt | while read sampleDir
do
    partialDir=${sampleDir##*/samples/} # strip directories
    diff -q $COMPAREDIR/samples/$partialDir/reads.all.pileup $sampleDir/reads.all.pileup
    stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
    echo -n $stat" " >> overnight_run.txt
done

echo -e "\ndiff var.flt.vcf"
echo -n "  var.flt.vcf:" >> overnight_run.txt
cat $WORKDIR/sampleDirectories.txt | while read sampleDir
do
    partialDir=${sampleDir##*/samples/} # strip directories
    diff -q $COMPAREDIR/samples/$partialDir/var.flt.vcf $sampleDir/var.flt.vcf
    stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
    echo -n $stat" " >> overnight_run.txt
done

echo -e "\ndiff snplist.txt"
echo -n "  snplist.txt:" >> overnight_run.txt
diff -q $COMPAREDIR/snplist.txt   $WORKDIR/snplist.txt
stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
echo -n $stat" " >> overnight_run.txt

echo -e "\ndiff reads.snp.pileup"
echo -n "  reads.snp.pileup:" >> overnight_run.txt
cat $WORKDIR/sampleDirectories.txt | while read sampleDir
do
    partialDir=${sampleDir##*/samples/} # strip directories
    diff -q $COMPAREDIR/samples/$partialDir/reads.snp.pileup $sampleDir/reads.snp.pileup
    stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
    echo -n $stat" " >> overnight_run.txt
done

echo -e "\ndiff snpma.fasta"
echo -n "  snpma.fasta:" >> overnight_run.txt
diff -q $COMPAREDIR/snpma.fasta   $WORKDIR/snpma.fasta
stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
echo -n $stat" " >> overnight_run.txt

echo -e "\ndiff referenceSNP.fasta"
echo -n "  referenceSNP.fasta:" >> overnight_run.txt
diff -q $COMPAREDIR/referenceSNP.fasta   $WORKDIR/referenceSNP.fasta
stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
echo -n $stat" " >> overnight_run.txt

echo >> overnight_run.txt

exit 0
