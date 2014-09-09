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
#   20140728-scd: Changed the comparison logic include the upstream files to detect reproducibility problems sooner.
#   20140728-scd: Changed the comparison logic to easily allow running a test repeatedly overnight to verify repeatable results.
#   20140818-scd: Split the create_snp_matrix script into 4 smaller scripts.
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
rm -rf $WORKDIR
cp -v -r -s $CLEANDIR $WORKDIR
cd $WORKDIR
# get sample directories sorted by size, largest first
ls -d samples/* | xargs ls -L -s -m | grep -E "(samples|total)" | sed 'N;s/\n//;s/:total//' | sort -k 2 -n -r | cut -f 1 -d " " > sampleDirectories.txt
sampleCount=$(cat sampleDirectories.txt | wc -l)
# get the *.fastq or *.fq files in each sample directory, possibly compresessed, on one line per sample, ready to feed to bowtie
TMPFILE1=$(mktemp tmp.fastqs.XXXXXXXX)
cat sampleDirectories.txt | while read dir; do echo $dir/*.fastq* >> $TMPFILE1; echo $dir/*.fq* >> $TMPFILE1; done
grep -v '*.fq*' $TMPFILE1 | grep -v '*.fastq*' > sampleFullPathNames.txt
rm $TMPFILE1

echo -e "\nStep 2 - Prep the reference"
referenceFilePath=$(ls reference/*.fasta)
if [[ $PLATFORM == torque ]]; then
    prepReferenceJobId=$(echo prepReference.sh $referenceFilePath | qsub -d $WORKDIR -N job.prepReference -j oe)
else
    prepReference.sh $referenceFilePath
fi

echo -e "\nStep 3 - Align the samples to the reference"
if [[ $PLATFORM == torque ]]; then
    numAlignThreads=8
    alignSamplesJobId=$(echo | qsub -t 1-$sampleCount -d $WORKDIR -N job.alignSamples -j oe -W depend=afterok:$prepReferenceJobId -l nodes=1:ppn=$numAlignThreads << _EOF_
    alignSamplesParameters=\$(cat sampleFullPathNames.txt | head -n \$PBS_ARRAYID | tail -n 1)
    alignSampleToReference.sh -p $numAlignThreads $referenceFilePath \$alignSamplesParameters
_EOF_
)
else
    cat sampleFullPathNames.txt | xargs --max-args=2 --max-lines=1 alignSampleToReference.sh -p $NUMCORES $referenceFilePath
fi

echo -e "\nStep 4 - Prep the samples"
if [[ $PLATFORM == torque ]]; then
    alignSamplesJobArray=${alignSamplesJobId%%.*}
    prepSamplesJobId=$(echo | qsub -t 1-$sampleCount -d $WORKDIR -N job.prepSamples -j oe -W depend=afterokarray:$alignSamplesJobArray << _EOF_
    sampleDir=\$(cat sampleDirectories.txt | head -n \$PBS_ARRAYID | tail -n 1)
    prepSamples.sh $referenceFilePath \$sampleDir
_EOF_
)
else
    cat sampleDirectories.txt | xargs -n 1 -P $NUMCORES prepSamples.sh $referenceFilePath
fi

echo -e "\nStep 5 - Combine the SNP positions across all samples into the SNP list file"
if [[ $PLATFORM == torque ]]; then
    prepSamplesJobArray=${prepSamplesJobId%%.*}
    snpListJobId=$(echo create_snp_list.py -n var.flt.vcf -o snplist.txt sampleDirectories.txt | qsub -d $WORKDIR -N job.snpList -j oe -W depend=afterokarray:$prepSamplesJobArray)
else
    create_snp_list.py -n var.flt.vcf -o snplist.txt sampleDirectories.txt
fi

echo -e "\nStep 6 - Create pileups at SNP positions for each sample"
if [[ $PLATFORM == torque ]]; then
    snpPileupJobId=$(echo | qsub -t 1-$sampleCount -d $WORKDIR -N job.snpPileup -j oe -W depend=afterok:$snpListJobId << _EOF_
    sampleDir=\$(cat sampleDirectories.txt | head -n \$PBS_ARRAYID | tail -n 1)
    create_snp_pileup.py -l snplist.txt -a \$sampleDir/reads.all.pileup -o \$sampleDir/reads.snp.pileup
_EOF_
)
else
    cat sampleDirectories.txt | xargs -n 1 -P $NUMCORES -I XX create_snp_pileup.py -l snplist.txt -a XX/reads.all.pileup -o XX/reads.snp.pileup
fi

echo -e "\nStep 7 - Create the SNP matrix"
if [[ $PLATFORM == torque ]]; then
    snpPileupJobArray=${snpPileupJobId%%.*}
    snpMatrixJobId=$(echo | qsub -d $WORKDIR -N job.snpMatrix -j oe -W depend=afterokarray:$snpPileupJobArray << _EOF_
    create_snp_matrix.py -l snplist.txt -p reads.snp.pileup -o snpma.fasta sampleDirectories.txt
_EOF_
)
else
    create_snp_matrix.py -l snplist.txt -p reads.snp.pileup -o snpma.fasta sampleDirectories.txt
fi    

echo -e "\nStep 8 - Create the reference base sequence"
if [[ $PLATFORM == torque ]]; then
    echo create_snp_reference_seq.py -l snplist.txt -o referenceSNP.fasta $referenceFilePath | qsub -d $WORKDIR -N job.snpReference -j oe -W depend=afterokarray:$snpPileupJobArray
else
    create_snp_reference_seq.py -l snplist.txt -o referenceSNP.fasta $referenceFilePath
fi

echo -e "\nStep 9 - compare results"
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

echo -e "\ndiff reads.sam"
echo -n "reads.sam:" >> ../overnight_run.txt
cat sampleDirectories.txt | while read sampleDir
do
    diff -q --ignore-matching-lines=bowtie $COMPAREDIR/$sampleDir/reads.sam $sampleDir/reads.sam
    stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
    echo -n $stat" " >> ../overnight_run.txt
done

echo -e "\ndiff reads.all.pileup"
echo -n "  reads.all.pileup:" >> ../overnight_run.txt
cat sampleDirectories.txt | while read sampleDir
do
    diff -q $COMPAREDIR/$sampleDir/reads.all.pileup $sampleDir/reads.all.pileup
    stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
    echo -n $stat" " >> ../overnight_run.txt
done

echo -e "\ndiff var.flt.vcf"
echo -n "  var.flt.vcf:" >> ../overnight_run.txt
cat sampleDirectories.txt | while read sampleDir
do
    diff -q $COMPAREDIR/$sampleDir/var.flt.vcf $sampleDir/var.flt.vcf
    stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
    echo -n $stat" " >> ../overnight_run.txt
done

echo -e "\ndiff snplist.txt"
echo -n "  snplist.txt:" >> ../overnight_run.txt
diff -q $COMPAREDIR/snplist.txt   snplist.txt
stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
echo -n $stat" " >> ../overnight_run.txt

echo -e "\ndiff reads.snp.pileup"
echo -n "  reads.snp.pileup:" >> ../overnight_run.txt
cat sampleDirectories.txt | while read sampleDir
do
    diff -q $COMPAREDIR/$sampleDir/reads.snp.pileup $sampleDir/reads.snp.pileup
    stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
    echo -n $stat" " >> ../overnight_run.txt
done

echo -e "\ndiff snpma.fasta"
echo -n "  snpma.fasta:" >> ../overnight_run.txt
diff -q $COMPAREDIR/snpma.fasta   snpma.fasta
stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
echo -n $stat" " >> ../overnight_run.txt

echo -e "\ndiff referenceSNP.fasta"
echo -n "  referenceSNP.fasta:" >> ../overnight_run.txt
diff -q $COMPAREDIR/referenceSNP.fasta   referenceSNP.fasta
stat=$([ "$?" -eq 0 ] && echo OK || echo xx)
echo -n $stat" " >> ../overnight_run.txt

echo >> ../overnight_run.txt

cd - > /dev/null
exit 0
