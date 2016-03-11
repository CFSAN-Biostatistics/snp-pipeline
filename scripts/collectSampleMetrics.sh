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
#   20151229-scd: Detect errors and prevent execution of unwanted processing when earlier processing steps fail.
#   20160224-scd: Capture separate metrics counting phase1 and phase2 snps.
#   20160224-scd: Reuse previously computed sam file and pileup metrics when re-running the pipeline.
#   20160225-scd: Default output to the metrics file, not stdout.
#   20160226-scd: Add the average insert size metric
#   20160304-scd: Count missing positions in the consensus.fasta file instead of the snp matrix.
#   20160309-scd: Add the -m option to specify the maximum allowable snp threshold.  Add the excludedSample flag.
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
  echo usage: $0 [-h] [-f] [-c FILE] [-m INT ] [-o FILE] [-v FILE] sampleDir referenceFile
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
  echo '  -c FILE          : Relative or absolute path to the consensus fasta file'
  echo '                     (default: consensus.fasta in the sampleDir)'
  echo '  -m INT           : Maximum allowed number of SNPs per sample. (default: -1)'
  echo '  -o FILE          : Output file. Relative or absolute path to the metrics file'
  echo '                     (default: metrics in the sampleDir)'
  echo '  -v FILE          : Relative or absolute path to the consensus vcf file'
  echo '                     (default: consensus.vcf in the sampleDir)'
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

while getopts ":hfc:m:o:v:" option; do
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

let maxSnps=-1
if [ "$opt_m_set" = 1 ]; then
  if [[ "$opt_m_arg" =~ ^-?[0-9]+$ ]]; then
    let maxSnps=$opt_m_arg
  else
    echo "Invalid argument for option -m, expected an integer"
    usage
    exit 3
  fi
fi

logSysEnvironment $@
setupSampleErrorHandler


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
if [[   -f "$sampleDir" ]]; then sampleError "$sampleDir is a file, expecting a directory."; false; fi
if [[ ! -d "$sampleDir" ]]; then sampleError "Sample directory $sampleDir does not exist."; false; fi
sampleDirFiles=$(ls -A "$sampleDir")
if [[ -z "$sampleDirFiles" ]]; then sampleError "Sample directory $sampleDir is empty."; false; fi

# Get the reference file
referenceFile="$2"
if [ "$referenceFile" = "" ]; then
  echo
  echo "Missing reference file." 1>&2
  usage
  exit 10
fi
if [[ ! -e "$referenceFile" ]]; then globalError "Reference file $referenceFile does not exist."; fi
if [[ ! -f "$referenceFile" ]]; then globalError "Reference file $referenceFile is not a file."; fi
if [[ ! -s "$referenceFile" ]]; then globalError "Reference file $referenceFile is empty."; fi


# Extra arguments not allowed
if [[ "$3" != "" ]]; then
  echo "Unexpected argument \"$2\" specified after the sample directory" 1>&2
  echo
  usage
  exit 20
fi


errorList=""
sampleDirBasename=${sampleDir##*/}


#-----------------
# Redirect output
#-----------------

if [ "$opt_o_set" = "1" ]; then
  outfile="$opt_o_arg"
else
  outfile=${sampleDir}/metrics
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
# Calculate number of reads and %mapped from sam file
#-------------------------
samFile=$sampleDir/reads.sam
if [ -s "$samFile" ]; then
  # Metrics already freshly collected?
  unset nreads
  unset perc_mapped
  if [[ "$opt_f_set" != "1" && -s "$outfile" && "$outfile" -nt "$samFile" ]]; then
    readParameter "$outfile" "numberReads" # reuse already fresh metrics
    readParameter "$outfile" "percentReadsMapped" # reuse already fresh metrics
    nreads=$numberReads
    perc_mapped=$percentReadsMapped
  fi
  if [[ -n $nreads && -n $perc_mapped ]]; then
    echo "# "$(date +"%Y-%m-%d %T") Reusing previously calculated number of reads and %mapped 1>&2
  else
    echo "# "$(date +"%Y-%m-%d %T") Calculate number of reads and %mapped from sam file 1>&2
    nreads=$(samtools view -S -c $samFile)
    mapped=$(samtools view -S -c -F 4 $samFile)
    perc_mapped=$(bc <<< "scale=4; ($mapped/$nreads)*100")
    perc_mapped=$(printf '%.2f' $perc_mapped)
  fi
else
  error="SAM file was not found."
  echo "$error" 1>&2
  errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
fi


#-------------------------
# Calculate mean insert size
#-------------------------
bamFile=${sampleDir}/reads.sorted.bam
if [ -s "$bamFile" ]; then
  # Metrics already freshly collected?
  unset insertSize
  if [[ "$opt_f_set" != "1" && -s "$outfile" && "$outfile" -nt "$bamFile" ]]; then
    readParameter "$outfile" "aveInsertSize" # reuse already fresh metrics
    insertSize=$aveInsertSize
  fi
  if [[ -n $insertSize ]]; then
    echo "# "$(date +"%Y-%m-%d %T") Reusing previously calculated mean insert size 1>&2
  else
    echo "# "$(date +"%Y-%m-%d %T") Calculate mean insert size from bam file 1>&2
    # Extract inferred insert sizes (TLEN, column 9 of BAM file) for reads "mapped in proper pair" (2) and "first in pair" (64) = 66
    tmpFile=$(mktemp -p "$sampleDir" tmp.inserts.XXXXXXXX)
    samtools view -f 66 "$bamFile" | cut -f 9 | sed 's/^-//' > $tmpFile
    insertSum=$(awk '{s+=$1} END {print s}' $tmpFile)
    insertCnt=$(wc -l $tmpFile | cut -d ' ' -f 1)
    rm $tmpFile
    if [[ -z $insertSum || -z $insertCnt || $insertCnt == 0 ]]; then
      insertSize=""
      error="Cannot calculate mean insert size."
      echo "$error" 1>&2
      errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
    else
      insertSize=$(printf %.2f $(echo "$insertSum / $insertCnt" | bc -l))
    fi
  fi
else
  error="$bamFile file was not found."
  echo "$error" 1>&2
  errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
fi


#-------------------------
# Calculate mean depth from pileup file
#-------------------------
pileupFile=${sampleDir}/reads.all.pileup
if [ -s "$pileupFile" ]; then
  # Metrics already freshly collected?
  unset depth
  if [[ "$opt_f_set" != "1" && -s "$outfile" && "$outfile" -nt "$pileupFile" ]]; then
    readParameter "$outfile" "avePileupDepth" # reuse already fresh metrics
    depth=$avePileupDepth
  fi
  if [[ -n $depth ]]; then
    echo "# "$(date +"%Y-%m-%d %T") Reusing previously calculated mean pileup depth 1>&2
  else
    echo "# "$(date +"%Y-%m-%d %T") Calculate mean depth from pileup file 1>&2
    sumDepth=$(awk '{sum+=$4}END{print sum}' "$pileupFile")
    refLength=$(grep -v '>' "$referenceFile" | wc -m)
    depth=$(printf %.2f $(echo "$sumDepth / $refLength" | bc -l))
  fi
else
  error="Pileup file was not found."
  echo "$error" 1>&2
  errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
fi

#-------------------------
echo "# "$(date +"%Y-%m-%d %T") Count number of high confidence SNP positions from phase 1 vcf file 1>&2
#-------------------------
excludedSample=""
vcfFile=${sampleDir}/var.flt.vcf
if [ -s "$vcfFile" ]; then
  phase1Snps=$(grep -c -v '^#' "$vcfFile") || true
  # Flag excessive snps
  if [[ $maxSnps -ge 0 && $phase1Snps -gt $maxSnps ]]; then
    excludedSample="Excluded"
    error="Excluded: exceeded $maxSnps maxsnps."
    echo "$error" 1>&2
    errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
  fi
else
  error="VCF file $vcfFile was not found."
  echo "$error" 1>&2
  errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
fi

#-------------------------
echo "# "$(date +"%Y-%m-%d %T") Count number of consensus snps from consensus vcf file 1>&2
#-------------------------
# Get the consensus VCF file
if [ "$opt_v_set" = "1" ]; then
  consensusVcfFile="$opt_v_arg"
else
  consensusVcfFile=${sampleDir}/consensus.vcf
fi
if [ "$excludedSample" = "Excluded" ]; then
  phase2Snps="" # Omit the phase2 snp count. It will be meaningless since this sample's phase1 snps are excluded from the snplist.
elif [ -s "$consensusVcfFile" ]; then
  phase2Snps=$(python << END
from __future__ import print_function
import vcf
num_snps = 0
with open("$consensusVcfFile") as inp:
    reader = vcf.VCFReader(inp)
    for record in reader:
        if not record.is_snp: # is ALT not in [A,C,G,T,N,*] ?
            continue
        if record.ALT == record.REF:
            continue
        for sample in record.samples:
            if not sample.is_variant: # is GT == REF ?
                continue
            if sample.data.FT != "PASS":
                continue
            num_snps += 1
print(num_snps)
END
)
else
  error="Consensus VCF file $consensusVcfFile was not found."
  echo "$error" 1>&2
  errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
fi

#------------------------------------------
echo "# "$(date +"%Y-%m-%d %T") Count missing positions in the snp matrix 1>&2
#------------------------------------------

# Get the consensus fasta file
if [ "$opt_c_set" = "1" ]; then
  consensusFastaFile="$opt_c_arg"
else
  consensusFastaFile=${sampleDir}/consensus.fasta
fi
if [ "$excludedSample" = "Excluded" ]; then
  missingPos="" # Omit the gap count. It will be meaningless since this sample's phase1 snps are excluded from the snplist.
elif [ -s "$consensusFastaFile" ]; then
  missingPos=$(python << END
from __future__ import print_function
from Bio import SeqIO
handle = open("$consensusFastaFile", "rU")
for record in SeqIO.parse(handle, "fasta"):
    if record.id == "$sampleDirBasename":
        missing = record.seq.count('-')
        print(missing)
        break
handle.close()
END
)
else
  error="Consensus fasta file $consensusFastaFile was not found."
  echo "$error" 1>&2
  errorList=${errorList}${errorList:+" "}$"$error"  # Insert spaces between errors
fi


#-------------------------
echo "# "$(date +"%Y-%m-%d %T") Print Results 1>&2
#-------------------------
echo "sample=\"$sampleDirBasename\"" > "$outfile"
echo "fastqFileList=\"$fastqFileList\"" >> "$outfile"
echo "fastqFileSize=$file_size" >> "$outfile"
echo "machine=$machine" >> "$outfile"
echo "flowcell=$flowcell" >> "$outfile"
echo "numberReads=$nreads" >> "$outfile"
echo "percentReadsMapped=$perc_mapped" >> "$outfile"
echo "aveInsertSize=$insertSize" >> "$outfile"
echo "avePileupDepth=$depth" >> "$outfile"
echo "phase1Snps=$phase1Snps" >> "$outfile"
echo "snps=$phase2Snps" >> "$outfile"
echo "missingPos=$missingPos" >> "$outfile"
echo "excludedSample=$excludedSample" >> "$outfile"
echo "errorList=\"$errorList"\" >> "$outfile"

echo "# "$(date +"%Y-%m-%d %T") collectSampleMetrics.sh finished 1>&2
