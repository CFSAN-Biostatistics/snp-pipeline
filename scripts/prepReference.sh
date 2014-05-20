#!/bin/bash
#
#Author: Hugh A. Rand (har)
#Purpose: Prep the reference sequence for snppipline code.
#Input:
#    referenceName
#Output:
#    bowtie index files from reference sequence
#Use example:
#   prepReference.sh ERR178926
#History:
#   20140512-har: Started.
#   20140520-har: Download of sequence moved to different script.
#Notes:
#   1.Assumes file named 'referenceName.fasta' is in a directory 'reference' 
#Bugs:
#

#Process arguments
if [ -z "$1" ]; then
    echo usage: $0 referenceName
    exit
fi
REFERENCENAME=$1

#Create index file for reference
~/software/bowtie2-2.2.2/bowtie2-build reference/$REFERENCENAME'.fasta' reference/$REFERENCENAME

