#!/bin/bash
#Set up input for snppipline code.
#Some prep work needs to be done before running this, and is as follows:
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
  ~/software/bowtie2-2.2.2/bowtie2 -p 11 -q -x reference/$REFERENCENAME -1 samples/$SAMPLENAME/$SAMPLENAME'_1.fastq' -2 samples/$SAMPLENAME/$SAMPLENAME'_2.fastq' > samples/$SAMPLENAME/$SAMPLENAME'.sam'

#Convert to bam file with only mapped positions
  samtools view -bS -F 4 -o samples/$SAMPLENAME/$SAMPLENAME'.bam' samples/$SAMPLENAME/$SAMPLENAME'.sam'

#Convert to a sorted bam 
  samtools sort samples/$SAMPLENAME/$SAMPLENAME'.bam' samples/$SAMPLENAME/$SAMPLENAME'.sorted.bam'

#Get a bcf file from the pileup and bam file
  samtools mpileup -uf reference/NC_011149.fasta samples/$SAMPLENAME/$SAMPLENAME'.sorted.bam.bam' | bcftools view -bvcg - > samples/$SAMPLENAME/$SAMPLENAME'.bcf'

#Convert bcf to vcf
  bcftools view samples/$SAMPLENAME/$SAMPLENAME'.bcf' | vcfutils.pl varFilter -D1000 > samples/$SAMPLENAME/$SAMPLENAME'.vcf'
