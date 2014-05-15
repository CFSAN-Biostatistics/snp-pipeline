#!/bin/bash
#
#Directives for Portable Batch System (PBS) if HPC with Torque or equivalent is installed.
#PBS -N prepSamples
#PBS -m be
#PBS -j oe
#PBS -M hugh.rand@fda.hhs.gov    #TODO Set this to be your email address
#
#Author: Hugh A. Rand (har)
#Purpose: Set up sample input for snppipline code.
#Input:
#    referenceName
#    sampleName
#    various files too tedious to explain
#Output:
#    various files too tedious to explain
#Use example:
#   On workstation with one sample
#       prepSamples.sh ERR178926 NC_011149
#   On workstation with multiple samples
#       echo -e "ERR178926\nERR178927\nERR178928\nERR178929\nERR178930\n" > prepInput
#       cat prepInput | xargs -n 1 prepSamples.sh NC_011149
#   With PBS
#       qsub -d $PWD temp.sh ERR178926 NC_011149
#       qsub -d $PWD temp1.sh
#History:
#  20140512-har: Started.
#Notes:
#Bugs:
#
#

#Process arguments
if [ -z "$1" ]; then
    echo usage: $0 referenceName sampleName
    exit
fi
REFERENCENAME=$1
SAMPLENAME=$2

#Set up directories
mkdir -p samples

#Get the sample sequences
fastq-dump --outdir samples/$SAMPLENAME --split-files $SAMPLENAME

#Align sequences to reference
~/software/bowtie2-2.2.2/bowtie2 -p 11 -q -x reference/$REFERENCENAME -1 samples/$SAMPLENAME/$SAMPLENAME'_1.fastq' -2 samples/$SAMPLENAME/$SAMPLENAME'_2.fastq' > samples/$SAMPLENAME/'reads.sam'

#Convert to bam file with only mapped positions
samtools view -bS -F 4 -o samples/$SAMPLENAME/'reads.bam' samples/$SAMPLENAME/'reads.sam'

#Convert to a sorted bam 
samtools sort samples/$SAMPLENAME/'reads.bam' samples/$SAMPLENAME/'reads.sorted'

#Get a bcf file from the pileup and bam file
samtools mpileup -uf reference/NC_011149.fasta samples/$SAMPLENAME/'reads.sorted.bam' | bcftools view -bvcg - > samples/$SAMPLENAME/'reads.bcf'

#Convert bcf to vcf
bcftools view samples/$SAMPLENAME/'reads.bcf' | vcfutils.pl varFilter -D1000 > samples/$SAMPLENAME/'var.flt.vcf'
