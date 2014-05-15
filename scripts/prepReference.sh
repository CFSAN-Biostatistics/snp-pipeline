#!/bin/bash
#
#Author: Hugh A. Rand (har)
#Purpose: Set up reference sequence input for snppipline code.
#Input:
#    referenceName
#Output:
#    various files too tedious to explain
#Use example:
#   prepReference.sh ERR178926
#History:
#   20140512-har: Started.
#Notes:
#Bugs:
#

#Process arguments
if [ -z "$1" ]; then
    echo usage: $0 referenceName
    exit
fi
REFERENCENAME=$1

#Set up directories
mkdir -p reference

#Get the Reference sequence
~/mnt/biob/svn/Biostats/rand/cfsanutils/scripts/fetch.py $REFERENCENAME -e hugh.rand@fda.hhs.gov > reference/$REFERENCENAME'.fasta'

#Create index file for reference
~/software/bowtie2-2.2.2/bowtie2-build reference/$REFERENCENAME'.fasta' reference/$REFERENCENAME

