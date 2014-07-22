#!/bin/bash
#
#Directives for Portable Batch System (PBS) if HPC with Torque or equivalent is installed.
#PBS -N job.alignSamples
#PBS -m be
#PBS -j oe
#PBS -M user.name@fda.hhs.gov    #TODO Set this to be your email address
#
#Author: Hugh A. Rand (har)
#        Steven C. Davis (scd)
#Purpose: Aligns sample sequence data to reference for snppipline code.
#Input:
#    numBowtieThreads
#    referenceDir/referenceName (without the .fasta extension)
#    samplePath to fastq file
#    [optional samplePath to mate fastq file if paired]
#Output:
#    SAM file
#    written into the same directories containing the input sample fastq files
#Use example:
#   On workstation with one sample, unpaired
#       alignSampleToReference.sh Users/NC_011149 Users/ERR178926.fastq
#   On workstation with one sample, paired
#       alignSampleToReference.sh Users/NC_011149 Users/ERR178926_1.fastq Users/ERR178926_2.fastq
#   On workstation with multiple samples
#       NUMCORES=$(grep -c ^processor /proc/cpuinfo)
#       ls -d --color=never samples/* > sampleDirectoryNames.txt
#       rm sampleFullPathNames.txt
#       cat sampleDirectoryNames.txt | while read dir; do echo $dir/*.fastq >> sampleFullPathNames.txt; done
#       cat sampleFullPathNames.txt | xargs --max-args=2 --max-lines=1 alignSampleToReference.sh $NUMCORES reference/NC_011149
#   On a workstation with gnu parallel:
#       cat sampleFullPathNames.txt | parallel alignSampleToReference.sh $NUMCORES reference/NC_011149
#   With PBS
#       qsub -d $PWD temp.sh ERR178926 NC_011149
#       qsub -d $PWD temp1.sh
#History:
#   20140715-scd: Started - based on previous code copied from prepSamples.sh
#   20140721-scd: Print the bowtie version and command line to facilitate troubleshooting
#Notes:
#   1. Assumes file named 'referenceName.fasta' in the referenceDir directory
#   2. Assumes the sequence file(s) have names '*_1.fastq' and '*_2.fastq', if paired.
#   2. Assumes the sequence file(s) have names '*.fastq', if unpaired.
#Bugs:
#

#Process arguments

if (($# < 3)) || (($# > 4)); then
    echo usage: $0 numBowtieThreads referencePath sampleFastqPath1 [sampleFastqPath2]
    echo '      numBowtieThreads : number of CPU cores to use concurrently during bowtie alignment'
    echo '      referencePath    : relative or absolute path to the reference, without the .fasta extension'
    echo '      sampleFastqPath1 : relative or absolute path to the to fastq file'
    echo '      sampleFastqPath2 : optional relative or absolute path to the mate fastq file, if paired'
    exit 1
fi

NUMCORES=$1
REFERENCEPATH=$2
SAMPLEPATH1=$3
if (($# == 4)); then
    SAMPLEPATH2=$4
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
    bowtie2 --version | sed 's/^/# /'
    if [ $SAMPLEPATH2 ]; then
        echo bowtie2 -p $NUMCORES -q -x $REFERENCEPATH -1 $SAMPLEPATH1 -2 $SAMPLEPATH2
        bowtie2 -p $NUMCORES -q -x $REFERENCEPATH -1 $SAMPLEPATH1 -2 $SAMPLEPATH2 > $SAMPLEDIR/'reads.sam'
    else
        echo bowtie2 -p $NUMCORES -q -x $REFERENCEPATH $SAMPLEPATH1
        bowtie2 -p $NUMCORES -q -x $REFERENCEPATH $SAMPLEPATH1 > $SAMPLEDIR/'reads.sam'
    fi
fi
