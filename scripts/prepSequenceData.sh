#!/bin/bash
#
#Author: Hugh A. Rand (har)
#Purpose: Set up reference and sample sequence for snppipline code.
#Input:
#    referenceName
#    sampleName
#Output:
#    reference and sample sequence files
#Use example:
#   With one sample
#       prepSequenceData.sh NC_011149 ERR178926
#   With multiple samples
#       echo -e "ERR178926\nERR178927\nERR178928\nERR178929\nERR178930\n" > prepInput
#       cat prepInput | xargs -n 1 prepSequenceData.sh NC_011149
#History:
#   20140520-har: Started.
#Notes:
#Bugs:
#

#Process arguments
if [ -z "$1" ]; then
    echo usage: $0 referenceName sampleName
    exit
fi
REFERENCENAME=$1
SAMPLENAME=$2

#Set up directories
mkdir -p reference samples

#Get the Reference sequence
~/mnt/biob/svn/Biostats/rand/cfsanutils/scripts/fetch.py $REFERENCENAME -e hugh.rand@fda.hhs.gov > reference/$REFERENCENAME'.fasta'

#Get the sample sequences
fastq-dump --outdir samples/$SAMPLENAME --split-files $SAMPLENAME


