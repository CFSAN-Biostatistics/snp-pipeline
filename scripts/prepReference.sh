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
#    referenceDir/referenceFilePath (with the fasta extension)
#Output:
#    bowtie index files from reference sequence written to the reference subdirectory
#Use example:
#   prepReference.sh reference/ERR178926.fasta
#History:
#   20140512-har: Started.
#   20140520-har: Download of sequence moved to different script.
#   20140612-scd: Removed the hardcoded path to bowtie2.  It must be on the $PATH now.
#   20140623-scd: Changed calling convention to match prepSamples.sh -- referenceDir is expected in the command parameter
#   20140721-scd: Print the bowtie version and command line to facilitate troubleshooting
#   20140905-scd: Use getopts to parse the command arguments.  Improved online help.
#   20140905-scd: Expects the full file name of the reference on the command line.
#   20140910-scd: Outputs are not rebuilt when already fresh, unless the -f (force) option is specified.
#   20141003-scd: Enhance log output
#   20141019-scd: Use the configuration parameter environment variable.
#   20141020-scd: Check for empty target files.
#   20141030-scd: Fix Python 2.6 compatibility issue when logging RAM size.
#   20150109-scd: Log the Grid Engine job ID.
#Notes:
#
#Bugs:
#
#References:
#   http://stackoverflow.com/questions/14008125/shell-script-common-template
#

usage()
{
    echo usage: $0 [-h] [-f] referenceFile
    echo
    echo 'Index the reference genome for subsequent alignment, and create'
    echo 'the faidx index file for subsequent pileups. The output is written'
    echo 'to the reference directory.'
    echo
    echo 'Positional arguments:'
    echo '  referenceFile    : Relative or absolute path to the reference fasta file'
    echo
    echo 'Options:'
    echo '  -h               : Show this help message and exit'
    echo '  -f               : Force processing even when result files already exist and'
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
echo "# Job ID            : $JOB_ID"
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

#Create index file for reference
if [[ $opt_f_set != "1" && -s "$referenceBasePath.rev.1.bt2" && "$referenceBasePath.rev.1.bt2" -nt "$referenceFilePath" ]]; then
    echo "# Bowtie index $referenceBasePath.rev.1.bt2 is already freshly built.  Use the -f option to force a rebuild."
else
    echo "# "$(date +"%Y-%m-%d %T") bowtie2-build $Bowtie2Build_ExtraParams "$referenceFilePath" "$referenceBasePath"
    echo "# "$(bowtie2-build --version | grep -i -E "bowtie.*version")
    bowtie2-build $Bowtie2Build_ExtraParams "$referenceFilePath" "$referenceBasePath"
    echo
fi

#Create fai index
if [[ "$opt_f_set" != "1" && -s "$referenceFilePath.fai" && "$referenceFilePath.fai" -nt "$referenceFilePath" ]]; then
    echo "# SAMtools fai index $referenceFilePath.fai is already freshly built.  Use the -f option to force a rebuild."
else
    echo "# "$(date +"%Y-%m-%d %T") samtools faidx $SamtoolsFaidx_ExtraParams "$referenceFilePath"
    echo "# SAMtools "$(samtools 2>&1 > /dev/null | grep Version)
    samtools faidx $SamtoolsFaidx_ExtraParams "$referenceFilePath"
    echo
fi

echo "# "$(date +"%Y-%m-%d %T") prepReference.sh finished
