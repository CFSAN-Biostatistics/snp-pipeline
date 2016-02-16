#!/bin/bash
#
#Author: Steve Davis (scd)
#
#Purpose: Create a multi-sample VCF file from per-sample VCF files.
#
#Input:
#    sampleDir
#
#Output:
#    snpma.vcf
#
#Use example:
#   ls -d samples/* > sampleDirectoryNames.txt
#   cat sampleDirectoryNames.txt | xargs -n 1 mergeVcf.sh sampleDirectoryNames.txt
#
#History:
#   20150915-scd: Started.
#   20151217-scd: Detect errors and prevent execution of unwanted processing when earlier processing steps fail.
#
#Notes:
#
#Bugs:
#

# source the utility functions
. snp_pipeline_inc.sh

usage()
{
  echo usage: $0 [-h] [-f] [-n NAME] [-o FILE] sampleDirsFile
  echo
  echo 'Merge the vcf files from all samples into a single multi-vcf file for all samples.'
  echo
  echo 'Before running this command, the vcf file for each sample must be created by the'
  echo 'call_consensus.py script.'
  echo
  echo 'Positional arguments:'
  echo '  sampleDirsFile   : Relative or absolute path to file containing a list of'
  echo '                     directories -- one per sample'
  echo
  echo 'Options:'
  echo '  -h               : Show this help message and exit'
  echo '  -f               : Force processing even when result files already exist and '
  echo '                     are newer than inputs'
  echo '  -n NAME          : File name of the vcf files which must exist in each of'
  echo '                     the sample directories. (default: consensus.vcf)'
  echo '  -o FILE          : Output file. Relative or absolute path to the merged'
  echo '                     multi-vcf file. (default: snpma.vcf)'
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
    echo "# Job ID            : $JOB_ID"
    fi
    echo "# Hostname          :" $(hostname)
    echo "# RAM               :" $(python -c 'from __future__ import print_function; import psutil; import locale; locale.setlocale(locale.LC_ALL, ""); print("%s MB" % locale.format("%d", psutil.virtual_memory().total / 1024 / 1024, grouping=True))')
    echo
}

#--------
# Options
#--------

while getopts ":hfn:o:" option; do
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

# Determine the name of the per-sample vcf files
if [ "$opt_n_set" = "1" ]; then
  inFileName="$opt_n_arg"
else
  inFileName="consensus.vcf"
fi

# Determine the output file
if [ "$opt_o_set" = "1" ]; then
  outFilePath="$opt_o_arg"
else
  outFilePath="snpma.vcf"
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
  echo "Unexpected argument \"$2\" specified after the sample directories file" 1>&2
  echo 1>&2
  usage
  exit 20
fi

logSysEnvironment $@
setupGlobalErrorHandler

if [[ ! -e "$sampleDirsFile" ]]; then globalError "Sample directories file $sampleDirsFile does not exist."; fi
if [[ ! -f "$sampleDirsFile" ]]; then globalError "Sample directories file $sampleDirsFile is not a file."; fi
if [[ ! -s "$sampleDirsFile" ]]; then globalError "Sample directories file $sampleDirsFile is empty."; fi

# Validate the single sample VCF files
(( numGoodVcfFiles = 0 )) || true
(( needRebuild = 0 )) || true
while IFS='' read -r dir || [[ -n "$dir" ]]
do
  inFilePath="$dir/$inFileName"
  if   [[ ! -e "$inFilePath" ]]; then sampleError "Sample vcf file $inFilePath does not exist." true;
  elif [[ ! -f "$inFilePath" ]]; then sampleError "Sample vcf file $inFilePath is not a file." true;
  elif [[ ! -s "$inFilePath" ]]; then sampleError "Sample vcf file $inFilePath is empty." true;
  else
    (( numGoodVcfFiles++ )) || true;
    lastGoodVcfFile="$inFilePath"
  fi
  # Check if rebuild is needed
  if [[ "$opt_f_set" == "1" || ! -s "$outFilePath" || "$inFilePath" -nt "$outFilePath" ]]; then
    (( needRebuild++ )) || true
  fi
done < "$sampleDirsFile"

if [ $numGoodVcfFiles = 0 ]; then
  globalError "There are no vcf files to merge."
fi

if [ $needRebuild -eq 0 ]; then
  echo "# Multi-VCF file is already freshly created.  Use the -f option to force a rebuild."
  echo "# "$(date +"%Y-%m-%d %T") mergeVcf.sh finished
  exit 0
fi

# If there is only one good sample, just copy the consensus VCF file to the snpma.vcf file
if [ $numGoodVcfFiles = 1 ]; then
  cp "$lastGoodVcfFile" "$outFilePath"
  echo "# "$(date +"%Y-%m-%d %T") mergeVcf.sh finished
  exit 0
fi

# Copy single VCF files to a common directory where the files will be edited
echo "# "$(date +"%Y-%m-%d %T") Starting merge operations
tempDir=$(mktemp -d tmp.vcf.XXXXXXXX)
cat "$sampleDirsFile" | while IFS='' read -r dir || [[ -n "$dir" ]]
do
  inFilePath="$dir/$inFileName"
  echo Copying $inFilePath to $tempDir/$(basename $dir).vcf
  cp -p -u "$inFilePath" $tempDir/$(basename $dir).vcf || true
done

# Zip and create index for all sample files
for filestozip in $tempDir/*.vcf; do
  echo Compressing $filestozip
  bgzip -c $filestozip > $filestozip.gz
done
for filestoidx in $tempDir/*.gz; do
  echo Indexing $filestoidx
  tabix -f -p vcf $filestoidx
done

# Merge the VCFs
bcftools merge --info-rules NS:sum -o "$outFilePath" $tempDir/*.gz

# Clean up
if [[ -e "$tempDir" ]]; then 
    rm -rf "$tempDir"
fi

echo "# "$(date +"%Y-%m-%d %T") mergeVcf.sh finished
