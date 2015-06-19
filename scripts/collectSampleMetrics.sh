#!/bin/bash
#
#Authors: Joseph D. Baugher (jdb)
#         Steve Davis (scd)
#
#Purpose: Collect various alignment and snp call metrics for a single specified sample.
#
#Input:
#    sampleDir
#Output:
#    file of name=value metrics
#Use example:
#   collectSampleMetrics.sh sampleDir
#History:
#   20150311-jdb: Initial version authored by Joseph D. Baugher.
#   20150324-scd: Modified to work as an integrated step within the CFSAN SNP Pipeline.
#   20150413-scd: Calculate average depth using all reference positions, including those with zero coverage.
#Notes:
#
#Bugs:
#


#------
# Usage
#------

usage()
{
  echo usage: $0 [-h] [-f] [-m FILE] [-o FILE] sampleDir referenceFile
  echo
  echo 'Collect alignment, coverage, and variant metrics for a single specified sample.'
  echo
  echo 'Positional arguments:'
  echo '  sampleDir        : Relative or absolute directory of the sample'
  echo '  referenceFile    : Relative or absolute path to the reference fasta file'
  echo
  echo 'Options:'
  echo '  -h               : Show this help message and exit'
  echo '  -f               : Force processing even when result files already exist and'
  echo '                     are newer than inputs'
  echo '  -m FILE          : Relative or absolute path to the SNP matrix file'
  echo '                     (default: snpma.fasta)'
  echo '  -o FILE          : Output file. Relative or absolute path to the metrics file'
  echo '                     (Default: stdout)'
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
  echo "# Job ID            : $JOB_ID[$SGE_TASK_ID]" 1>&2
  fi
  echo "# Hostname          :" $(hostname) 1>&2
  echo "# RAM               :" $(python -c 'from __future__ import print_function; import psutil; import locale; locale.setlocale(locale.LC_ALL, ""); print("%s MB" % locale.format("%d", psutil.virtual_memory().total / 1024 / 1024, grouping=True))') 1>&2
  echo 1>&2
}

#--------
# Options
#--------

while getopts ":hfm:o:" option; do
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

logSysEnvironment $@

# SNP Matrix file
if [ "$opt_m_set" = "1" ]; then
  snpmaFile="$opt_m_arg"
else
  snpmaFile="snpma.fasta"
fi
if [[ ! -e "$snpmaFile" ]]; then echo "SNP matrix file $snpmaFile does not exist." 1>&2; exit 3; fi
if [[ ! -f "$snpmaFile" ]]; then echo "SNP matrix file $snpmaFile is not a file." 1>&2; exit 3; fi
if [[ ! -s "$snpmaFile" ]]; then echo "SNP matrix file $snpmaFile is empty." 1>&2; exit 3; fi


#----------
# Arguments
#----------
shift $((OPTIND-1))

# Get the sample directory
sampleDir="$1"
sampleDir=${1%/}
if [ "$sampleDir" = "" ]; then
  echo
  echo "Missing sample directory." 1>&2
  usage
  exit 10
fi
if [[   -f "$sampleDir" ]]; then echo "$sampleDir is a file, expecting a directory."; exit 10 1>&2; fi
if [[ ! -d "$sampleDir" ]]; then echo "Sample directory $sampleDir does not exist."; exit 10 1>&2; fi
sampleDirFiles=$(ls -A "$sampleDir")
if [[ -z "$sampleDirFiles" ]]; then echo "Sample directory $sampleDir is empty."; exit 10 1>&2; fi

# Get the reference file
referenceFile="$2"
if [ "$referenceFile" = "" ]; then
  echo
  echo "Missing reference file." 1>&2
  usage
  exit 10
fi
if [[ ! -e "$referenceFile" ]]; then echo "Reference file $referenceFile does not exist." 1>&2; exit 10; fi
if [[ ! -f "$referenceFile" ]]; then echo "Reference file $referenceFile is not a file." 1>&2; exit 10; fi
if [[ ! -s "$referenceFile" ]]; then echo "Reference file $referenceFile is empty." 1>&2; exit 10; fi


# Extra arguments not allowed
if [[ "$3" != "" ]]; then
  echo "Unexpected argument \"$2\" specified after the sample directory" 1>&2
  echo
  usage
  exit 20
fi


errorList=""
sampleDirBasename=${sampleDir##*/}

#------------------------------------
# Metrics already freshly collected?
#------------------------------------
if [[ "$opt_f_set" != "1" && "$opt_o_set" = "1" && -s "$opt_o_arg" && "$opt_o_arg" -nt "$snpmaFile" ]]; then
    echo "Metrics are already freshly created for $sampleDirBasename.  Use the -f option to force a rebuild." 1>&2
    exit 0
fi


#-----------------
# Redirect output
#-----------------

# Create file descriptor 3 for output
if [ "$opt_o_set" = "1" ]; then
  # 3 points to a file
  exec 3>"$opt_o_arg"
else
  # 3 points to stdout
  exec 3>&1
fi


#-------------------------
echo "# "$(date +"%Y-%m-%d %T") Get machine and flowcell from fastq header 1>&2
#-------------------------

# Check if fastq.gz exist and are populated
fastqFileArray=(${sampleDir}/*.fastq.gz)
# if no fastq.gz files, look for fq.gz files
if [ ! -s ${fastqFileArray[0]} ]; then
  fastqFileArray=(${sampleDir}/*.fq.gz)
fi

if [ -s ${fastqFileArray[0]} ]; then
  fq_head=$(gzip -cd ${fastqFileArray[0]} | head -1)
  # @M01402:25:000000000-A6P1T:1:1101:11489:1986 1:N:0:2
  machine=${fq_head#@}
  machine=${machine%%:*}
  fields=(${fq_head//:/ })
  flowcell=${fields[2]}
else
  error="No compressed fastq.gz or fq.gz files were found."
  echo "$error" 1>&2
  errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
fi


#-------------------------
echo "# "$(date +"%Y-%m-%d %T") Sum file sizes of paired fastq files 1>&2
#-------------------------

# if no compressed fastq files, look for uncompressed fastq files
if [ ! -s ${fastqFileArray[0]} ]; then
  fastqFileArray=(${sampleDir}/*.fastq)
fi
# if no fastq files, look for fq files
if [ ! -s ${fastqFileArray[0]} ]; then
  fastqFileArray=(${sampleDir}/*.fq)
fi
if [ -s ${fastqFileArray[0]} ]; then
  fastqFileList=""
  file_size=0
  for file in ${fastqFileArray[@]}; do
    size=$(du --dereference --bytes ${file} | cut -f1)
    file_size=$((file_size+size))
    file=${file##*/}  # strip any leading directories, leaving the file name
    fastqFileList=${fastqFileList}${fastqFileList:+", "}$"$file"  # Insert commas between file names
  done
else
  error="No fastq or fq files were found."
  echo "$error" 1>&2
  errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
fi

#-------------------------
echo "# "$(date +"%Y-%m-%d %T") Calculate number of reads and %mapped from sam file 1>&2
#-------------------------
samFile=$sampleDir/reads.sam
if [ -s "$samFile" ]; then
  nreads=$(samtools view -S -c ${sampleDir}/reads.sam)
  mapped=$(samtools view -S -c -F 4 ${sampleDir}/reads.sam)
  perc_mapped=$(bc <<< "scale=4; ($mapped/$nreads)*100")
  perc_mapped=$(printf '%.2f' $perc_mapped)
else
  error="SAM file was not found."
  echo "$error" 1>&2
  errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
fi

#-------------------------
echo "# "$(date +"%Y-%m-%d %T") Calculate mean depth from pileup file 1>&2
#-------------------------
pileupFile=${sampleDir}/reads.all.pileup
if [ -s "$pileupFile" ]; then
  sumDepth=$(awk '{sum+=$4}END{print sum}' "$pileupFile")
  refLength=$(grep -v '>' "$referenceFile" | wc -m)
  depth=$(printf %.2f $(echo "$sumDepth / $refLength" | bc -l))
else
  error="Pileup file was not found."
  echo "$error" 1>&2
  errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
fi

#-------------------------
echo "# "$(date +"%Y-%m-%d %T") Count number of SNPs from vcf file 1>&2
#-------------------------
vcfFile=${sampleDir}/var.flt.vcf
if [ -s "$samFile" ]; then
  snps=$(grep -c -v '^#' "$vcfFile")
else
  error="VCF file was not found."
  echo "$error" 1>&2
  errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
fi

#------------------------------------------
echo "# "$(date +"%Y-%m-%d %T") Count missing positions in the snp matrix 1>&2
#------------------------------------------

missingPos=$(python << END
from Bio import SeqIO
handle = open("$snpmaFile", "rU")
for record in SeqIO.parse(handle, "fasta"):
    if record.id == "$sampleDirBasename":
        missing = record.seq.count('-')
        print missing
        break
handle.close()
END
)


#-------------------------
echo "# "$(date +"%Y-%m-%d %T") Print Results 1>&2
#-------------------------
echo "sample=\"$sampleDirBasename\"" >&3
echo "fastqFileList=\"$fastqFileList\"" >&3
echo "fastqFileSize=$file_size" >&3
echo "machine=$machine" >&3
echo "flowcell=$flowcell" >&3
echo "numberReads=$nreads" >&3
echo "percentReadsMapped=$perc_mapped" >&3
echo "avePileupDepth=$depth" >&3
echo "snps=$snps" >&3
echo "missingPos=$missingPos" >&3
echo "errorList=\"$errorList"\" >&3

#printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'  "$sampleDirBasename" "$fastqFileList" "$file_size" "$machine" "$flowcell" "$nreads" "$perc_mapped" "$depth" "$snps" "$errorList"

echo "# "$(date +"%Y-%m-%d %T") collectSampleMetrics.sh finished 1>&2
