#!/bin/bash
#
#Directives for Portable Batch System (PBS) if HPC with Torque or equivalent is installed.
#PBS -N prepSamples
#PBS -m be
#PBS -j oe
#PBS -M user.name@fda.hhs.gov    #TODO Set this to be your email address
#
#Author: Hugh A. Rand (har)
#        Steven C. Davis (scd)
#Purpose: Preps sample sequence data for snppipline code.
#Input:
#    referenceDir/referenceName
#    samplePath to fastq
#    [optional samplePath to mate fastq if paired]
#Output:
#    various files too tedious to explain
#    written into the same directories containing the input sample fastq files
#Use example:
#   On workstation with one sample, unpaired
#       prepSamples.sh Users/NC_011149 Users/ERR178926.fastq
#   On workstation with one sample, paired
#       prepSamples.sh Users/NC_011149 Users/ERR178926_1.fastq Users/ERR178926_2.fastq
#   On workstation with multiple samples
#       find samples -type f | grep fastq > prepInput
#       cat prepInput | xargs -n 1 prepSamples.sh reference/NC_011149
#	On a workstation with gnu parallel:
#       cat prepInput | parallel prepSamples.sh reference/NC_011149
#   With PBS
#       qsub -d $PWD temp.sh ERR178926 NC_011149
#       qsub -d $PWD temp1.sh
#History:
#   20140512-har: Started.
#   20140520-har: Download of sequence moved to different script.
#   20140623-scd: Changes for varscan.
#Notes:
#   1. Assumes file named 'referenceName.fasta' in the referenceDir directory
#   2. Assumes sequence file(s) are paired end and names '*_1.fastq' and '*_2.fastq'
#Bugs:
#   1. Should add prints to stdout to show progress to user
#

#Process arguments

if (($# < 2)) || (($# > 3)); then
    echo usage: $0 referencePath sampleFastqPath1 [sampleFastqPath2]
    exit 1
fi

REFERENCEPATH=$1
SAMPLEPATH1=$2
if (($# == 3)); then
	SAMPLEPATH2=$3
fi
SAMPLEDIR=${SAMPLEPATH1%/*}
SAMPLEID=${SAMPLEPATH1##*/} # strip the directory
SAMPLEID=${SAMPLEID%_1.fastq} # strip the file extension and leading _1 if any
SAMPLEID=${SAMPLEID%.fastq} # strip the file extension regardless of leading _1
REFERENCEID=${REFERENCEPATH##*/}

#Check if alignment to reference has been done; if not align sequences to reference
if [ -s $SAMPLEDIR/reads.sam ]; then
    echo '**'$SAMPLEID' has already been aligned to '$REFERENCEID
else
    echo '**Align sequence '$SAMPLEID' to reference '$REFERENCEID
    if [ $SAMPLEPATH2 ]; then
        bowtie2 -p 11 -q -x $REFERENCEPATH -1 $SAMPLEPATH1 -2 $SAMPLEPATH2 > $SAMPLEDIR/'reads.sam'
	else
        bowtie2 -p 11 -q -x $REFERENCEPATH $SAMPLEPATH1 > $SAMPLEDIR/'reads.sam'
	fi
fi

#Check if bam file exists; if not convert to bam file with only mapped positions
if [ -s $SAMPLEDIR/reads.bam ]; then
	echo '**Bam file already exists for '$SAMPLEID
else
	echo '**Convert sam file to bam file with only mapped positions.'
	samtools view -bS -F 4 -o $SAMPLEDIR/'reads.unsorted.bam' $SAMPLEDIR/'reads.sam'
	#Convert to a sorted bam
	echo '**Convert bam to sorted bam file.'
	samtools sort $SAMPLEDIR/'reads.unsorted.bam' $SAMPLEDIR/'reads'
fi

#Check if pileup present; if not create it 
if [ -s $SAMPLEDIR/'reads.all.pileup' ]; then
	echo '**'$SAMPLEID'.pileup already exists'
else
	echo '**Produce bcf file from pileup and bam file.'
	samtools mpileup -f $REFERENCEPATH'.fasta' $SAMPLEDIR/'reads.bam' > $SAMPLEDIR/'reads.all.pileup'
fi

#Check if unfiltered vcf exists; if not create it
if [ -s $SAMPLEDIR/'var.flt.vcf' ]; then
	echo '**vcf file already exists for '$SAMPLEID
else
	echo '**Creating vcf file'
	java -jar /usr/bin/VarScan.jar mpileup2snp $SAMPLEDIR/'reads.all.pileup' --min-var-freq 0.90 --output-vcf 1 > $SAMPLEDIR/'var.flt.vcf'
fi
