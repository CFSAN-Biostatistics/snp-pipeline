#!/bin/bash
#
#Directives for Portable Batch System (PBS) if HPC with Torque or equivalent is installed.
#PBS -N job.prepReference
#PBS -m be
#PBS -j oe
#PBS -M user.name@fda.hhs.gov    #TODO Set this to be your email address
#
#Author: Hugh A. Rand (har)
#        Steven C. Davis (scd)
#Purpose: Prep the reference sequence for snppipline code.
#Input:
#    referenceDir/referenceName (without the fasta extension)
#Output:
#    bowtie index files from reference sequence written to the reference subdirectory
#Use example:
#   prepReference.sh reference/ERR178926
#History:
#   20140512-har: Started.
#   20140520-har: Download of sequence moved to different script.
#   20140612-scd: Removed the hardcoded path to bowtie2.  It must be on the $PATH now.
#   20140623-scd: Changed calling convention to match prepSamples.sh -- referenceDir is expected in the command parameter
#   20140721-scd: Print the bowtie version and command line to facilitate troubleshooting
#Notes:
#   1. Assumes a file named 'referenceName.fasta' is in the referenceDir directory
#Bugs:
#
#References:
#   http://stackoverflow.com/questions/14008125/shell-script-common-template
#

#Setup-------------------------------------------------------------

USAGE="-h referencePath"

#Options processing------------------------------------------------

if [ $# == 0 ] ; then
    echo usage: $0 referencePath
    exit 1;
fi

while getopts "h" optname
  do
    case "$optname" in
      "h")
	echo usage: $0 referencePath
        exit 0;
        ;;
      "?")
        echo "Unknown option $OPTARG"
        exit 0;
        ;;
      ":")
        echo "No argument value for option $OPTARG"
        exit 0;
        ;;
      *)
        echo "Unknown error while processing options"
        exit 0;
        ;;
    esac
  done

shift $(($OPTIND - 1))

REFERENCEPATH=$1

#Body--------------------------------------------------------------

#Create index file for reference
echo -e "\n"bowtie2-build $REFERENCEPATH'.fasta' $REFERENCEPATH"\n"
bowtie2-build --version | sed 's/^/# /'
bowtie2-build $REFERENCEPATH'.fasta' $REFERENCEPATH

#Create fai index
echo -e "\n"samtools faidx $REFERENCEPATH'.fasta'"\n"
samtools 2>&1 > /dev/null | grep Version | sed 's/^/# SAMtools /'
samtools faidx $REFERENCEPATH'.fasta'

#Wrap Up-----------------------------------------------------------
