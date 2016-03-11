#!/bin/bash
#
#Authors: Steve Davis (scd)
#
#Purpose: Combine the various alignment and snp call metrics for all samples into one file.
#
#Input:
#    file containing a list of all the sample directories
#Output:
#    tab separated value file with one line per sample, and columns for each metric
#Use example:
#   combineSampleMetrics.sh sampleDirsFile.txt
#History:
#   20150319-scd: Started
#   20150324-scd: Integrated into the snp pipeline.
#   20150413-scd: Fix the sun grid engine "undefined" task id.
#   20151230-scd: Detect errors and prevent execution of unwanted processing when earlier processing steps fail.
#   20160226-scd: Add the average insert size metric
#   20160301-scd: Emit column headings with underscores.
#   20160309-scd: Add the Excluded Sample metric.
#Notes:
#
#Bugs:
#

# source the utility functions
. snp_pipeline_inc.sh


#------
# Usage
#------

usage()
{
  echo usage: $0 [-h] [-n NAME] [-o FILE] sampleDirsFile
  echo
  echo 'Combine the metrics from all samples into a single table of metrics for all samples.'
  echo 'The output is a tab-separated-values file with a row for each sample and a column'
  echo 'for each metric.'
  echo
  echo 'Before running this command, the metrics for each sample must be created by the'
  echo 'collectSampleMetrics.sh script.'
  echo
  echo 'Positional arguments:'
  echo '  sampleDirsFile   : Relative or absolute path to file containing a list of'
  echo '                     directories -- one per sample'
  echo
  echo 'Options:'
  echo '  -h               : Show this help message and exit'
  echo '  -n NAME          : File name of the metrics files which must exist in each of'
  echo '                     the sample directories. (default: metrics)'
  echo '  -o FILE          : Output file. Relative or absolute path to the combined metrics'
  echo '                     file. (default: stdout)'
  echo '  -s               : Emit column headings with spaces instead of underscores'
}

# --------------------------------------------------------
# Log the starting conditions
# --------------------------------------------------------
logSysEnvironment()
{
  echo "# Command           : $0 $@" 1>&2
  echo "# Working Directory : $(pwd)" 1>&2
  if [[ "$PBS_JOBID" != "" ]]; then
  echo "# Job ID            : $PBS_JOBID" 1>&2
  elif [[ "$JOB_ID" != "" ]]; then
  echo "# Job ID            : $JOB_ID" 1>&2
  fi
  echo "# Hostname          :" $(hostname) 1>&2
  echo "# RAM               :" $(python -c 'from __future__ import print_function; import psutil; import locale; locale.setlocale(locale.LC_ALL, ""); print("%s MB" % locale.format("%d", psutil.virtual_memory().total / 1024 / 1024, grouping=True))') 1>&2
  echo 1>&2
}

#--------
# Options
#--------

while getopts ":hn:o:s" option; do
  if [ "$option" = "h" ]; then
    usage
    exit 0
  elif [ "$option" = "?" ]; then
    echo "Invalid option -- '$OPTARG'" 1>&2
    usage
    exit 1
  elif [ "$option" = ":" ]; then
    echo "Missing argument for option -- '$OPTARG'" 1>&2
    usage
    exit 2
  else
    declare opt_"$option"_set="1"
    if [ "$OPTARG" != "" ]; then
      declare opt_"$option"_arg="$OPTARG"
    fi
  fi
done

# Determine the name of the metrics file
if [ "$opt_n_set" = "1" ]; then
  metricsFileName="$opt_n_arg"
else
  metricsFileName="metrics"
fi

# Create file descriptor 3 for output
if [ "$opt_o_set" = "1" ]; then
  # 3 points to a file
  exec 3>"$opt_o_arg"
else
  # 3 points to stdout
  exec 3>&1
fi

#----------
# Arguments
#----------
shift $((OPTIND-1))

# Get the sample directories file
sampleDirsFile="$1"
if [ "$sampleDirsFile" = "" ]; then
  echo "Missing sample directories file." 1>&2
  echo 1>&2
  usage
  exit 10
fi

# Extra arguments not allowed
if [[ "$2" != "" ]]; then
  echo "Unexpected argument \"$2\" specified after the sample directory" 1>&2
  echo 1>&2
  usage
  exit 20
fi

logSysEnvironment $@
setupGlobalErrorHandler

if [[ ! -e "$sampleDirsFile" ]]; then globalError "Sample directories file $sampleDirsFile does not exist."; fi
if [[ ! -f "$sampleDirsFile" ]]; then globalError "Sample directories file $sampleDirsFile is not a file."; fi
if [[ ! -s "$sampleDirsFile" ]]; then globalError "Sample directories file $sampleDirsFile is empty."; fi


#-------------------------------------------------------
# Parse the metrics files and print the tabular results
#-------------------------------------------------------

if [ "$opt_s_set" = "1" ]; then
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'  "Sample" "Fastq Files" "Fastq File Size" "Machine" "Flowcell" "Number of Reads" "Percent of Reads Mapped" "Average Insert Size" "Average Pileup Depth" "Phase1 SNPs" "Phase2 SNPs" "Missing SNP Matrix Positions" "Excluded Sample" "Warnings and Errors" >&3
else
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'  "Sample" "Fastq_Files" "Fastq_File_Size" "Machine" "Flowcell" "Number_of_Reads" "Percent_of_Reads_Mapped" "Average_Insert_Size" "Average_Pileup_Depth" "Phase1_SNPs" "Phase2_SNPs" "Missing_SNP_Matrix_Positions" "Excluded_Sample" "Warnings_and_Errors" >&3
fi

cat "$sampleDirsFile" | while IFS='' read -r dir || [[ -n "$dir" ]]
do
  metricsFilePath="$dir/$metricsFileName"
  echo Processing $metricsFilePath 1>&2
  if [[ ! -e "$metricsFilePath" ]]; then echo "Sample metrics file $metricsFilePath does not exist." >&3; sampleWarning "Sample metrics file $metricsFilePath does not exist."; continue; fi
  if [[ ! -f "$metricsFilePath" ]]; then echo "Sample metrics file $metricsFilePath is not a file."  >&3; sampleWarning "Sample metrics file $metricsFilePath is not a file." ; continue; fi
  if [[ ! -s "$metricsFilePath" ]]; then echo "Sample metrics file $metricsFilePath is empty."       >&3; sampleWarning "Sample metrics file $metricsFilePath is empty."      ; continue; fi
  while IFS='' read -r line || [[ -n "$line" ]]
  do
    if [[ "$line" =~ "=" ]]; then
      param=${line%%=*}  # strip everything after =
      value=${line##*=}  # strip everything before =
      declare "$param"="$value"
    fi
  done  < "$metricsFilePath" # This syntax without piping is needed to retain the values of variables declared in the loop
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'  "$sample" "$fastqFileList" "$fastqFileSize" "$machine" "$flowcell" "$numberReads" "$percentReadsMapped" "$aveInsertSize" "$avePileupDepth" "$phase1Snps" "$snps" "$missingPos" "$excludedSample" "$errorList" >&3
done

echo "# "$(date +"%Y-%m-%d %T") combineSampleMetrics.sh finished 1>&2
