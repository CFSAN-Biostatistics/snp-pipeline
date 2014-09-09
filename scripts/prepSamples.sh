#!/bin/bash
#
#Directives for Portable Batch System (PBS) if HPC with Torque or equivalent is installed.
#PBS -N job.prepSamples
#PBS -m be
#PBS -j oe
#PBS -M user.name@fda.hhs.gov    #TODO Set this to be your email address
#
#Author: Hugh A. Rand (har)
#        Steven C. Davis (scd)
#Purpose: Preps sample sequence data for snppipline code.
#Input:
#    referenceDir/referenceName
#    sampleDir
#Output:
#    various files too tedious to explain
#    written into the same directories containing the input sample fastq files
#Use example:
#   On workstation with one sample, unpaired
#       prepSamples.sh Users/NC_011149.fasta Users/ERR178926
#   On workstation with multiple samples
#       ls -d --color=never samples/* > sampleDirectoryNames.txt
#       cat sampleDirectoryNames.txt | xargs -n 1 prepSamples.sh reference/NC_011149.fasta
#   On a workstation with gnu parallel:
#       cat sampleDirectoryNames.txt | parallel prepSamples.sh reference/NC_011149.fasta
#   With PBS
#       qsub -d $PWD temp.sh ERR178926 NC_011149
#       qsub -d $PWD temp1.sh
#History:
#   20140512-har: Started.
#   20140520-har: Download of sequence moved to different script.
#   20140623-scd: Changes for varscan.
#   20140715-scd: Moved the bowtie align to the alignSampleToReference.sh script.
#   20140728-scd: Log the executed commands with all options to stdout.
#   20140905-scd: Use getopts to parse the command arguments.  Improved online help.
#   20140905-scd: Expects the full file name of the reference on the command line.
#Notes:
#
#Bugs:
#   1. Should add prints to stdout to show progress to user
#

usage()
{
    echo usage: $0 [-h] referenceFile sampleDir
    echo
    echo 'Find variants in a specified sample.'
    echo 'The output files are written to the sample directory.'
    echo
    echo 'positional arguments:'
    echo '  referenceFile    : Relative or absolute path to the reference fasta file'
    echo '  sampleDir        : Relative or absolute directory of the sample'
    echo
    echo 'options:'
    echo '  -h               : Show this help message and exit'
    echo
}

# --------------------------------------------------------
# getopts command line option handler: 

# For each valid option, 
#   If it is given, create a var dynamically to
#   indicate it is set: $opt_name_set = 1

#   If var gets an arg, create another var to
#   hold its value: $opt_name_arg = some value

# For invalid options given, 
#   Invoke Usage routine

# precede option list with a colon
# option list is a list of allowed option characters
# options that require an arg are followed by a colon

# example: ":abc:d"
# -abc 14 -d

while getopts ":h" option; do
  if [ "$option" = "h" ]; then
    usage
    exit 0
  elif [ "$option" = "?" ]; then
    echo
    echo "Invalid option -- '$OPTARG'"
    usage
    exit 1
  elif [ "$option" = ":" ]; then
    echo
    echo "Missing argument for option -- '$OPTARG'"
    usage
    exit 2
  else
    declare opt_"$option"_set="1"
    if [ "$OPTARG" != "" ]; then
      declare opt_"$option"_arg="$OPTARG"
    fi
  fi
done

# --------------------------------------------------------
# get the arguments

shift $((OPTIND-1))
referenceFilePath="$1"
if [ "$referenceFilePath" = "" ]; then
  echo "Missing reference file"
  echo
  usage
  exit 3
fi

sampleDir="$2"
if [ "$sampleDir" = "" ]; then
  echo "Missing sample directory"
  echo
  usage
  exit 4
fi

sampleId=${sampleDir##*/} # strip the parent directories

#Check if bam file exists; if not convert to bam file with only mapped positions
if [ -s $sampleDir/reads.bam ]; then
    echo '**Bam file already exists for '$sampleId
else
    echo '**Convert sam file to bam file with only mapped positions.'
    echo "# "$(date +"%Y-%m-%d %T") samtools view -bS -F 4 -o $sampleDir/'reads.unsorted.bam' $sampleDir/'reads.sam'
    echo "# SAMtools "$(samtools 2>&1 > /dev/null | grep Version)
    samtools view -bS -F 4 -o $sampleDir/'reads.unsorted.bam' $sampleDir/'reads.sam'
    #Convert to a sorted bam
    echo '**Convert bam to sorted bam file.'
    echo "# "$(date +"%Y-%m-%d %T") samtools sort $sampleDir/'reads.unsorted.bam' $sampleDir/'reads'
    echo "# SAMtools "$(samtools 2>&1 > /dev/null | grep Version)
    samtools sort $sampleDir/'reads.unsorted.bam' $sampleDir/'reads'
fi

#Check if pileup present; if not create it 
if [ -s $sampleDir/'reads.all.pileup' ]; then
    echo '**'$sampleId'.pileup already exists'
else
    echo '**Create pileup from bam file.'
    echo "# "$(date +"%Y-%m-%d %T") samtools mpileup -f $referenceFilePath $sampleDir/'reads.bam'
    echo "# SAMtools "$(samtools 2>&1 > /dev/null | grep Version)
    samtools mpileup -f $referenceFilePath $sampleDir/'reads.bam' > $sampleDir/'reads.all.pileup'
fi

#Check if unfiltered vcf exists; if not create it
if [ -s $sampleDir/'var.flt.vcf' ]; then
    echo '**vcf file already exists for '$sampleId
else
    echo '**Create vcf file'
    if [ ! -z "$CLASSPATH" ]; then
        echo "# "$(date +"%Y-%m-%d %T") java net.sf.varscan.VarScan mpileup2snp $sampleDir/'reads.all.pileup' --min-var-freq 0.90 --output-vcf 1
        echo "# "$(java net.sf.varscan.VarScan 2>&1 > /dev/null | head -n 1)
        java net.sf.varscan.VarScan mpileup2snp $sampleDir/'reads.all.pileup' --min-var-freq 0.90 --output-vcf 1 > $sampleDir/'var.flt.vcf'
    else
        echo '*** Error: cannot execute VarScan. Define the path to VarScan in the CLASSPATH environment variable.'
        exit 2
    fi
fi

echo "# "$(date +"%Y-%m-%d %T") prepSamples.sh finished
