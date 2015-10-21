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
#
#Notes:
#
#Bugs:
#


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
  echo '                     the sample directories. Default: consensus.vcf)'
  echo '  -o FILE          : Output file. Relative or absolute path to the merged'
  echo '                     multi-vcf file. Default: snpma.vcf)'
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

if [[ ! -e "$sampleDirsFile" ]]; then echo "Sample directories file $sampleDirsFile does not exist." 1>&2; exit 30; fi
if [[ ! -f "$sampleDirsFile" ]]; then echo "Sample directories file $sampleDirsFile is not a file." 1>&2; exit 30; fi
if [[ ! -s "$sampleDirsFile" ]]; then echo "Sample directories file $sampleDirsFile is empty." 1>&2; exit 30; fi


# Validate the single sample VCF files
(( needRebuild = 0 ))
while IFS='' read -r dir || [[ -n "$dir" ]]
do
  inFilePath="$dir/$inFileName"
  if [[ ! -e "$inFilePath" ]]; then echo "Sample vcf file $inFilePath does not exist." 1>&2; exit 40; fi
  if [[ ! -f "$inFilePath" ]]; then echo "Sample vcf file $inFilePath is not a file." 1>&2; exit 40; fi
  if [[ ! -s "$inFilePath" ]]; then echo "Sample vcf file $inFilePath is empty." 1>&2; exit 40; fi
  # Check if rebuild is needed
  if [[ "$opt_f_set" == "1" || ! -s "$outFilePath" || "$inFilePath" -nt "$outFilePath" ]]; then
    (( needRebuild++ ))
  fi
done < "$sampleDirsFile" 

if [ $needRebuild -eq 0 ]; then
  echo "# Multi-VCF file is already freshly created.  Use the -f option to force a rebuild."
  echo "# "$(date +"%Y-%m-%d %T") mergeVcf.sh finished
  exit 0
fi

# Copy single VCF files to a common directory where the files will be edited
echo "# "$(date +"%Y-%m-%d %T") Starting merge operations
tempDir=$(mktemp -d tmp.vcf.XXXXXXXX)
cat "$sampleDirsFile" | while IFS='' read -r dir || [[ -n "$dir" ]]
do
  inFilePath="$dir/$inFileName"
  echo Processing $inFilePath 1>&2
  if [[ ! -e "$inFilePath" ]]; then echo "Sample vcf file $inFilePath does not exist." 1>&2; exit 40; fi
  if [[ ! -f "$inFilePath" ]]; then echo "Sample vcf file $inFilePath is not a file." 1>&2; exit 40; fi
  if [[ ! -s "$inFilePath" ]]; then echo "Sample vcf file $inFilePath is empty." 1>&2; exit 40; fi
  cp -p -u "$inFilePath" $tempDir/$(basename $dir).vcf
done

# Replace "Sample1" with sample names in all single sample files
cd $tempDir
for filename in *.vcf; do
  samplename="${filename%.*}"
  sed -i s/Sample1/$samplename/ $filename
done

# Zip and create index for all sample files
for filestozip in *.vcf; do bgzip -c $filestozip > $filestozip.gz; done
for filestoidx in *.gz; do tabix -f -p vcf $filestoidx; done
cd - > /dev/null

# Merge the VCFs
bcftools merge --info-rules NS:sum -o "$outFilePath" $tempDir/*.gz

# Clean up
if [[ -e "$tempDir" ]]; then 
    rm -rf "$tempDir"
fi

echo "# "$(date +"%Y-%m-%d %T") mergeVcf.sh finished

