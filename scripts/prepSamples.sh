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
#       ls -d samples/* > sampleDirectoryNames.txt
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
#   20140910-scd: Outputs are not rebuilt when already fresh, unless the -f (force) option is specified.
#   20141003-scd: Enhance log output
#   20141019-scd: Use the configuration parameter environment variables.
#   20141020-scd: Check for empty target files.
#   20141030-scd: Fix Python 2.6 compatibility issue when logging RAM size.
#   20150109-scd: Log the Grid Engine job ID.
#   20151214-scd: Detect errors and prevent execution of unwanted processing when earlier processing steps fail.
#   20160414-scd: Support SAMtools 1.3.
#Notes:
#
#Bugs:
#   1. Should add prints to stdout to show progress to user
#

# source the utility functions
. snp_pipeline_inc.sh

usage()
{
    echo usage: $0 [-h] [-f] referenceFile sampleDir
    echo
    echo 'Find variants in a specified sample.'
    echo 'The output files are written to the sample directory.'
    echo
    echo 'Positional arguments:'
    echo '  referenceFile    : Relative or absolute path to the reference fasta file'
    echo '  sampleDir        : Relative or absolute directory of the sample'
    echo
    echo 'Options:'
    echo '  -h               : Show this help message and exit'
    echo '  -f               : Force processing even when result files already exist and '
    echo '                     are newer than inputs'
    echo
}

# --------------------------------------------------------
# Log the starting conditions
logSysEnvironment()
{
    echo "# Command           : $0 $@"
    echo "# Working Directory : $(pwd)"
    if [[ "$PBS_JOBID" != "" ]]; then
    echo "# Job ID            : $PBS_JOBID"
    elif [[ "$JOB_ID" != "" ]]; then
    echo "# Job ID            : $JOB_ID[$SGE_TASK_ID]"
    fi
    echo "# Hostname          :" $(hostname)
    echo "# RAM               :" $(python -c 'from __future__ import print_function; import psutil; import locale; locale.setlocale(locale.LC_ALL, ""); print("%s MB" % locale.format("%d", psutil.virtual_memory().total / 1024 / 1024, grouping=True))')
    echo "# CLASSPATH         :" $CLASSPATH
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

while getopts ":hf" option; do
  if [ "$option" = "h" ]; then
    usage
    exit 0
  elif [ "$option" = "?" ]; then
    echo "Invalid option -- '$OPTARG'"
    usage
    exit 1
  elif [ "$option" = ":" ]; then
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

logSysEnvironment $@
setupSampleErrorHandler

# Verify reference fasta file exists and is not empty
if [[ ! -e "$referenceFilePath" ]]; then globalError "Reference file $referenceFilePath does not exist."; fi
if [[ ! -f "$referenceFilePath" ]]; then globalError "Reference file $referenceFilePath is not a file."; fi
if [[ ! -s "$referenceFilePath" ]]; then globalError "Reference file $referenceFilePath is empty."; fi

# Verify SAM file exists and is not empty
samFile="$sampleDir/reads.sam"
if [[ ! -e "$samFile" ]]; then sampleError "Sample SAM file $samFile does not exist." false; fi
if [[ ! -f "$samFile" ]]; then sampleError "Sample SAM file $samFile is not a file." false; fi
if [[ ! -s "$samFile" ]]; then sampleError "Sample SAM file $samFile is empty." false; fi

sampleId=${sampleDir##*/} # strip the parent directories

# Substitute the default parameters if the user did not specify samtools view parameters
defaultParams="-F 4"
SamtoolsSamFilter_Params=${SamtoolsSamFilter_ExtraParams:-$defaultParams}

# Check for fresh bam file; if not, convert to bam file with only mapped positions
if [[ "$opt_f_set" != "1" && -s "$sampleDir/reads.unsorted.bam" && "$sampleDir/reads.unsorted.bam" -nt "$sampleDir/reads.sam" ]]; then
    echo "# Unsorted bam file is already freshly created for $sampleId.  Use the -f option to force a rebuild."
else
    echo "# Convert sam file to bam file with only mapped positions."
    echo "# "$(date +"%Y-%m-%d %T") samtools view -S -b $SamtoolsSamFilter_Params -o "$sampleDir/reads.unsorted.bam" "$sampleDir/reads.sam"
    echo "# SAMtools "$(samtools 2>&1 > /dev/null | grep Version)
    samtools view -S -b $SamtoolsSamFilter_Params -o "$sampleDir/reads.unsorted.bam" "$sampleDir/reads.sam"
    sampleErrorOnMissingFile "$sampleDir/reads.unsorted.bam" "samtools view"
    echo
fi


# Check for fresh sorted bam file; if not, sort it
if [[ "$opt_f_set" != "1" && -s "$sampleDir/reads.sorted.bam" && "$sampleDir/reads.sorted.bam" -nt "$sampleDir/reads.unsorted.bam" ]]; then
    echo "# Sorted bam file is already freshly created for $sampleId.  Use the -f option to force a rebuild."
else
    #Convert to a sorted bam
    echo "# Convert bam to sorted bam file."
    # Inspect the samtools version to determine how to execute samtools
    # Use the -o FILE command line option with SAMtools 1.3 and higher
    samtoolsVersion=$(samtools 2>&1 > /dev/null | grep Version | cut -d " " -f 2)
    if [[ $samtoolsVersion < '1.3' ]]; then
        echo "# "$(date +"%Y-%m-%d %T") samtools sort $SamtoolsSort_ExtraParams "$sampleDir/reads.unsorted.bam" "$sampleDir/reads.sorted"
        echo "# SAMtools "$(samtools 2>&1 > /dev/null | grep Version)
        samtools sort $SamtoolsSort_ExtraParams "$sampleDir/reads.unsorted.bam" "$sampleDir/reads.sorted"
    else
        echo "# "$(date +"%Y-%m-%d %T") samtools sort $SamtoolsSort_ExtraParams -o "$sampleDir/reads.sorted.bam" "$sampleDir/reads.unsorted.bam"
        echo "# SAMtools "$(samtools 2>&1 > /dev/null | grep Version)
        samtools sort $SamtoolsSort_ExtraParams -o "$sampleDir/reads.sorted.bam" "$sampleDir/reads.unsorted.bam"
    fi
    sampleErrorOnMissingFile "$sampleDir/reads.sorted.bam" "samtools sort"
    echo
fi

#Check for fresh pileup; if not, create it 
if [[ "$opt_f_set" != "1" && -s "$sampleDir/reads.all.pileup" && "$sampleDir/reads.all.pileup" -nt "$sampleDir/reads.sorted.bam" && "$sampleDir/reads.all.pileup" -nt "$referenceFilePath" ]]; then
    echo "# Pileup file is already freshly created for $sampleId.  Use the -f option to force a rebuild."
else
    echo "# Create pileup from bam file."
    echo "# "$(date +"%Y-%m-%d %T") samtools mpileup $SamtoolsMpileup_ExtraParams -f "$referenceFilePath" "$sampleDir/reads.sorted.bam"
    echo "# SAMtools "$(samtools 2>&1 > /dev/null | grep Version)
    samtools mpileup $SamtoolsMpileup_ExtraParams -f "$referenceFilePath" "$sampleDir/reads.sorted.bam" > "$sampleDir/reads.all.pileup"
    sampleErrorOnMissingFile "$sampleDir/reads.all.pileup" "samtools mpileup"
    echo
fi

#Check for fresh unfiltered vcf; if not, create it
if [[ "$opt_f_set" != "1" && -s "$sampleDir/var.flt.vcf" && "$sampleDir/var.flt.vcf" -nt "$sampleDir/reads.all.pileup" ]]; then
    echo "# Vcf file is already freshly created for $sampleId.  Use the -f option to force a rebuild."
else
    echo "# Create vcf file"
    if [ ! -z "$CLASSPATH" ]; then
        echo "# "$(date +"%Y-%m-%d %T") java $VarscanJvm_ExtraParams net.sf.varscan.VarScan mpileup2snp "$sampleDir/reads.all.pileup"  --output-vcf 1 $VarscanMpileup2snp_ExtraParams
        echo "# "$(java net.sf.varscan.VarScan 2>&1 > /dev/null | head -n 1)
        java $VarscanJvm_ExtraParams net.sf.varscan.VarScan mpileup2snp "$sampleDir/reads.all.pileup" --output-vcf 1 $VarscanMpileup2snp_ExtraParams > "$sampleDir/var.flt.vcf"
        sampleErrorOnMissingFile "$sampleDir/var.flt.vcf" "VarScan"
        sampleErrorOnFileContains "$sampleDir/var.flt.vcf" "OutOfMemoryError" "VarScan"
        sampleErrorOnFileContains "$sampleDir/var.flt.vcf" "Insufficient" "VarScan"
        echo
    else
        globalError "Error: cannot execute VarScan. Define the path to VarScan in the CLASSPATH environment variable."
    fi
fi

echo "# "$(date +"%Y-%m-%d %T") prepSamples.sh finished
