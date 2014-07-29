#!/bin/bash
#
#Directives for Portable Batch System (PBS) if HPC with Torque or equivalent is installed.
#PBS -N job.prepSamples
#PBS -m be
#PBS -j oe
#PBS -M user.name@fda.hhs.gov    #TODO Set this to be your email address
#
#Author: Hugh A. Rand (har)
#        Steven C. Davis (scd)
#Purpose: Preps sample sequence data for snppipline code.
#Input:
#    referenceDir/referenceName (without the .fasta extension)
#    sampleDir
#Output:
#    various files too tedious to explain
#    written into the same directories containing the input sample fastq files
#Use example:
#   On workstation with one sample, unpaired
#       prepSamples.sh Users/NC_011149 Users/ERR178926
#   On workstation with multiple samples
#       ls -d --color=never samples/* > sampleDirectoryNames.txt
#       cat sampleDirectoryNames.txt | xargs -n 1 prepSamples.sh reference/NC_011149
#   On a workstation with gnu parallel:
#       cat sampleDirectoryNames.txt | parallel prepSamples.sh reference/NC_011149
#   With PBS
#       qsub -d $PWD temp.sh ERR178926 NC_011149
#       qsub -d $PWD temp1.sh
#History:
#   20140512-har: Started.
#   20140520-har: Download of sequence moved to different script.
#   20140623-scd: Changes for varscan.
#   20140715-scd: Moved the bowtie align to the alignSampleToReference.sh script.
#   20140728-scd: Log the executed commands with all options to stdout.
#Notes:
#
#Bugs:
#   1. Should add prints to stdout to show progress to user
#

#Process arguments

if (($# < 2)) || (($# > 3)); then
    echo usage: $0 referencePath sampleDir
    echo '      referencePath : relative or absolute path to the reference, without the .fasta extension'
    echo '      sampleDir     : directory containing the sample'
    exit 1
fi

REFERENCEPATH=$1
SAMPLEDIR=$2
SAMPLEID=${SAMPLEDIR##*/} # strip the parent directories

#Check if bam file exists; if not convert to bam file with only mapped positions
if [ -s $SAMPLEDIR/reads.bam ]; then
    echo '**Bam file already exists for '$SAMPLEID
else
    echo '**Convert sam file to bam file with only mapped positions.'
    echo -e samtools view -bS -F 4 -o $SAMPLEDIR/'reads.unsorted.bam' $SAMPLEDIR/'reads.sam'
    samtools 2>&1 > /dev/null | grep Version | sed 's/^/# SAMtools /'
    samtools view -bS -F 4 -o $SAMPLEDIR/'reads.unsorted.bam' $SAMPLEDIR/'reads.sam'
    #Convert to a sorted bam
    echo '**Convert bam to sorted bam file.'
    echo -e samtools sort $SAMPLEDIR/'reads.unsorted.bam' $SAMPLEDIR/'reads'
    samtools 2>&1 > /dev/null | grep Version | sed 's/^/# SAMtools /'
    samtools sort $SAMPLEDIR/'reads.unsorted.bam' $SAMPLEDIR/'reads'
fi

#Check if pileup present; if not create it 
if [ -s $SAMPLEDIR/'reads.all.pileup' ]; then
    echo '**'$SAMPLEID'.pileup already exists'
else
    echo '**Produce bcf file from pileup and bam file.'
    echo -e samtools mpileup -f $REFERENCEPATH'.fasta' $SAMPLEDIR/'reads.bam'
    samtools 2>&1 > /dev/null | grep Version | sed 's/^/# SAMtools /'
    samtools mpileup -f $REFERENCEPATH'.fasta' $SAMPLEDIR/'reads.bam' > $SAMPLEDIR/'reads.all.pileup'
fi

#Check if unfiltered vcf exists; if not create it
if [ -s $SAMPLEDIR/'var.flt.vcf' ]; then
    echo '**vcf file already exists for '$SAMPLEID
else
    echo '**Creating vcf file'
    if [ ! -z "$CLASSPATH" ]; then
        echo -e java net.sf.varscan.VarScan mpileup2snp $SAMPLEDIR/'reads.all.pileup' --min-var-freq 0.90 --output-vcf 1
        java net.sf.varscan.VarScan 2>&1 > /dev/null | head -n 1 | sed 's/^/# /'
        java net.sf.varscan.VarScan mpileup2snp $SAMPLEDIR/'reads.all.pileup' --min-var-freq 0.90 --output-vcf 1 > $SAMPLEDIR/'var.flt.vcf'
    else
        echo '*** Error: cannot execute VarScan. Define the path to VarScan in the CLASSPATH environment variable.'
        exit 2
    fi
fi
