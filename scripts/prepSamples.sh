#!/bin/bash
#
#Directives for Portable Batch System (PBS) if HPC with Torque or equivalent is installed.
#PBS -N prepSamples
#PBS -m be
#PBS -j oe
#PBS -M hugh.rand@fda.hhs.gov    #TODO Set this to be your email address
#
#Author: Hugh A. Rand (har)
#Purpose: Preps sample sequence data for snppipline code.
#Input:
#    referenceDir/referenceName
#    samplePath
#Output:
#    various files too tedious to explain
#Use example:
#   On workstation with one sample
#       prepSamples.sh Users/NC_011149 Users/ERR178926
#   On workstation with multiple samples
#       ls -d samples/* > prepInput
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
if [ -z "$1" ]; then
    echo usage: $0 referencePath samplePath
    exit
fi
REFERENCEPATH=$1
SAMPLEPATH=$2
SAMPLEID=${SAMPLEPATH##*/}
REFERENCEID=${REFERENCEPATH##*/}

#Check if alignment to reference has been done; if not align sequences to reference
if [ -s $SAMPLEPATH/reads.sam ]; then
	echo '**'$SAMPLEID' has already been aligned to '$REFERENCEID
else
	echo '**Align sequence '$SAMPLEID' to reference '$REFERENCEID
	bowtie2 -p 11 -q -x $REFERENCEPATH -1 $SAMPLEPATH/$SAMPLEID'_1.fastq' -2 $SAMPLEPATH/$SAMPLEID'_2.fastq' > $SAMPLEPATH/'reads.sam'
fi

#Check if bam file exists; if not convert to bam file with only mapped positions
if [ -s $SAMPLEPATH/reads.bam ]; then
	echo '**Bam file already exists for '$SAMPLEID
else
	echo '**Convert sam file to bam file with only mapped positions.'
	samtools view -bS -F 4 -o $SAMPLEPATH/'reads.unsorted.bam' $SAMPLEPATH/'reads.sam'
	#Convert to a sorted bam
	echo '**Convert bam to sorted bam file.'
	samtools sort $SAMPLEPATH/'reads.unsorted.bam' $SAMPLEPATH/'reads'
fi

#Check if pileup present; if not create it 
if [ -s $SAMPLEPATH/'reads.all.pileup' ]; then
	echo '**'$SAMPLEID'.pileup already exists'
else
	echo '**Produce bcf file from pileup and bam file.'
	samtools mpileup -f $REFERENCEPATH'.fasta' $SAMPLEPATH/'reads.bam' > $SAMPLEPATH/'reads.all.pileup'
fi

#Check if unfiltered vcf exists; if not create it
if [ -s $SAMPLEPATH/'var.flt.vcf' ]; then
	echo '**vcf file already exists for '$SAMPLEID
else
	echo '**Creating vcf file'
	java -jar /usr/bin/VarScan.jar mpileup2snp $SAMPLEPATH/'reads.all.pileup' --min-var-freq 0.90 --output-vcf 1 > $SAMPLEPATH/'var.flt.vcf'
fi
