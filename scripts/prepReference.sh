#!/bin/bash
#
#Author: Hugh A. Rand (har)
#Purpose: Prep the reference sequence for snppipline code.
#Input:
#    referenceName
#Output:
#    bowtie index files from reference sequence written to the reference subdirectory
#Use example:
#   prepReference.sh ERR178926
#History:
#   20140512-har: Started.
#   20140520-har: Download of sequence moved to different script.
#   20140612-scd: Removed the hardcoded path to bowtie2.  It must be on the $PATH now.
#Notes:
#   1. Assumes a subdirectory named 'reference' below the current working directory
#   2. Assumes a file named 'referenceName.fasta' is in the 'reference' directory
#Bugs:
#

#Process arguments
if [ -z "$1" ]; then
    echo usage: $0 referenceName
    exit
fi
REFERENCENAME=$1

#Create index file for reference
bowtie2-build reference/$REFERENCENAME'.fasta' reference/$REFERENCENAME

