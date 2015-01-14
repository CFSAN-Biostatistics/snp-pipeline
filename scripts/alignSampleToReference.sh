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
#    referenceDir/referenceName
#    samplePath to fastq file
#    [optional samplePath to mate fastq file if paired]
#Output:
#    SAM file
#    written into the same directories containing the input sample fastq files
#Use example:
#   On workstation with one sample, unpaired
#       alignSampleToReference.sh Users/NC_011149.fasta Users/ERR178926.fastq
#   On workstation with one sample, paired
#       alignSampleToReference.sh Users/NC_011149.fasta Users/ERR178926_1.fastq Users/ERR178926_2.fastq
#   On workstation with multiple samples
#       numCores=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)
#       ls -d samples/* | xargs ls -L -s -m | grep -E "(samples|total)" | sed 'N;s/\n//;s/:total//' | sort -k 2 -n -r | cut -f 1 -d " " > sampleDirectories.txt
#       TMPFILE1=$(mktemp tmp.fastqs.XXXXXXXX)
#       cat sampleDirectories.txt | while read dir; do echo $dir/*.fastq* >> $TMPFILE1; echo $dir/*.fq* >> $TMPFILE1; done
#       grep -v '*.fq*' $TMPFILE1 | grep -v '*.fastq*' > sampleFullPathNames.txt
#       rm $TMPFILE1
#       cat sampleFullPathNames.txt | xargs -n 2 -L 1 alignSampleToReference.sh -p $NUMCORES Users/NC_011149.fasta
#   On a workstation with gnu parallel:
#       cat sampleFullPathNames.txt | parallel alignSampleToReference.sh -p 1 reference/NC_011149.fasta
#   With PBS
#       qsub -d $PWD temp.sh ERR178926 NC_011149
#       qsub -d $PWD temp1.sh
#History:
#   20140715-scd: Started - based on previous code copied from prepSamples.sh
#   20140721-scd: Print the bowtie version and command line to facilitate troubleshooting
#   20140903-scd: Use getopts to parse the command arguments.  Changed the calling convention.
#   20140905-scd: Expects the full file name of the reference on the command line.
#   20140910-scd: Outputs are not rebuilt when already fresh, unless the -f (force) option is specified.
#   20140919-scd: Handle spaces in file names.
#   20141003-scd: Enhance log output
#   20141009-scd: Remove the -p option to specify the number of CPU cores.
#   20141020-scd: Substitute the default parameters if the user did not specify bowtie parameters.
#   20141030-scd: Fix Python 2.6 compatibility issue when logging RAM size.
#   20150109-scd: Log the Grid Engine job ID.
#Notes:
#
#Bugs:
#


usage()
{
    echo usage: $0 [-h] [-f] referenceFile sampleFastqFile1 [sampleFastqFile2]
    echo
    echo 'Align the sequence reads for a specified sample to a specified reference genome.'
    echo 'The output is written to the file "reads.sam" in the sample directory.'
    echo
    echo 'Positional arguments:'
    echo '  referenceFile    : Relative or absolute path to the reference fasta file'
    echo '  sampleFastqFile1 : Relative or absolute path to the fastq file'
    echo '  sampleFastqFile2 : Optional relative or absolute path to the mate fastq file, if paired'
    echo
    echo 'Options:'
    echo '  -h               : Show this help message and exit'
    echo '  -f               : Force processing even when result files already exist and '
    echo '                     are newer than inputs'
    echo
}

# --------------------------------------------------------
# Log the starting conditions
echo "# Command           : $0 $@"
echo "# Working Directory : $(pwd)"
if [[ "$PBS_JOBID" != "" ]]; then
echo "# Job ID            : $PBS_JOBID"
elif [[ "$JOB_ID" != "" ]]; then
echo "# Job ID            : $JOB_ID[$SGE_TASK_ID]"
fi
echo "# Hostname          :" $(hostname)
echo "# RAM               :" $(python -c 'from __future__ import print_function; import psutil; import locale; locale.setlocale(locale.LC_ALL, ""); print("%s MB" % locale.format("%d", psutil.virtual_memory().total / 1024 / 1024, grouping=True))')
echo

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
referenceBasePath=${referenceFilePath%.fasta} # strip the file extension

sampleFilePath1="$2"
if [ "$sampleFilePath1" = "" ]; then
  echo "Missing sample fastq file 1"
  echo
  usage
  exit 4
fi

sampleFilePath2="$3"

sampleDir=${sampleFilePath1%/*}
sampleId=${sampleFilePath1##*/} # strip the directory
sampleId=${sampleId%.gz*} # strip the .gz file extension if any
sampleId=${sampleId%_1.fastq} # strip the _1.fastq file extension if any
sampleId=${sampleId%.fastq} # strip the .fastq file extension regardless of leading _1
sampleId=${sampleId%_1.fq} # strip the _1.fq file extension if any
sampleId=${sampleId%.fq} # strip the .fq file extension regardless of leading _1
referenceId=${referenceBasePath##*/} # strip the directory

# Parse the user-specified bowtie parameters to determine if the user specified the number of CPU cores
regex="(-p[[:blank:]]*)([[:digit:]]+)"
if [[ "$Bowtie2Align_ExtraParams" =~ $regex ]]; then
    numCoresParam=""
else
    # if not user-specified, default to 8 on HPC or all cpu cores on workstation
    if [[ "$PBS_JOBID" != "" ]]; then
        numCores=8
    else
        numCores=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)
    fi
    numCoresParam="-p $numCores"
fi


# echo ==============================
# echo numCores=$numCores
# echo referenceFilePath=$referenceFilePath
# echo sampleFilePath1=$sampleFilePath1
# echo sampleFilePath2=$sampleFilePath2
# if [ $sampleFilePath2 ]; then
#     echo sampleFilePath2 exists
# fi
# echo sampleDir=$sampleDir
# echo sampleId=$sampleId
# echo referenceId=$referenceId
# echo ==============================
# exit 0 


# Substitute the default parameters if the user did not specify bowtie parameters
defaultParams="--reorder -q"
Bowtie2Align_Params=${Bowtie2Align_ExtraParams:-$defaultParams}

#Check if alignment to reference has been done; if not align sequences to reference
if [ $sampleFilePath2 ]; then
    if [[ "$opt_f_set" != "1" && -s "$sampleDir/reads.sam" && "$sampleDir/reads.sam" -nt "$referenceBasePath.rev.1.bt2" && "$sampleDir/reads.sam" -nt "$sampleFilePath1" && "$sampleDir/reads.sam" -nt "$sampleFilePath2" ]]; then
        echo "# $sampleId has already been aligned to $referenceId.  Use the -f option to force a rebuild."
    else
        echo "# Align sequence $sampleId to reference $referenceId"
        echo "# "$(date +"%Y-%m-%d %T") bowtie2 $numCoresParam $Bowtie2Align_Params -x "$referenceBasePath" -1 "$sampleFilePath1" -2 "$sampleFilePath2"
        echo "# "$(bowtie2 --version | grep -i -E "bowtie.*version")
        bowtie2 $numCoresParam $Bowtie2Align_Params -x "$referenceBasePath" -1 "$sampleFilePath1" -2 "$sampleFilePath2" > "$sampleDir/reads.sam"
        echo
    fi
else
    if [[ "$opt_f_set" != "1" && -s "$sampleDir/reads.sam" && "$sampleDir/reads.sam" -nt "$referenceBasePath.rev.1.bt2" && "$sampleDir/reads.sam" -nt "$sampleFilePath1" ]]; then
        echo "# $sampleId has already been aligned to $referenceId.  Use the -f option to force a rebuild."
    else
        echo "# Align sequence $sampleId to reference $referenceId"
        echo "# "$(date +"%Y-%m-%d %T") bowtie2 $numCoresParam $Bowtie2Align_Params -x "$referenceBasePath" -U "$sampleFilePath1"
        echo "# "$(bowtie2 --version | grep -i -E "bowtie.*version")
        bowtie2 $numCoresParam $Bowtie2Align_Params -x "$referenceBasePath" -U "$sampleFilePath1" > "$sampleDir/reads.sam"
        echo
    fi
fi
echo $(date +"# %Y-%m-%d %T") alignSampleToReference.sh finished
