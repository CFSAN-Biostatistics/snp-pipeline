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
#       numCores=$(grep -c ^processor /proc/cpuinfo)
#       ls -d samples/* | xargs ls -L -s -m | grep -E "(samples|total)" | sed 'N;s/\n//;s/:total//' | sort -k 2 -n -r | cut -f 1 -d " " > sampleDirectories.txt
#       TMPFILE1=$(mktemp tmp.fastqs.XXXXXXXX)
#       cat sampleDirectories.txt | while read dir; do echo $dir/*.fastq* >> $TMPFILE1; echo $dir/*.fq* >> $TMPFILE1; done
#       grep -v '*.fq*' $TMPFILE1 | grep -v '*.fastq*' > sampleFullPathNames.txt
#       rm $TMPFILE1
#       cat sampleFullPathNames.txt | xargs --max-args=2 --max-lines=1 alignSampleToReference.sh -p $NUMCORES Users/NC_011149.fasta
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
#Notes:
#
#Bugs:
#


usage()
{
    echo usage: $0 [-h] [-p INT] referenceFile sampleFastqFile1 [sampleFastqFile2]
    echo
    echo 'Align the sequence reads for a specified sample to a specified reference genome.'
    echo 'The output is written to the file "reads.sam" in the sample directory.'
    echo
    echo 'positional arguments:'
    echo '  referenceFile    : Relative or absolute path to the reference fasta file'
    echo '  sampleFastqFile1 : Relative or absolute path to the fastq file'
    echo '  sampleFastqFile2 : Optional relative or absolute path to the mate fastq file, if paired'
    echo
    echo 'options:'
    echo '  -h               : Show this help message and exit'
    echo '  -p INT           : Number of CPU cores to use concurrently during bowtie alignment (default: all)'
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

while getopts ":p:h" option; do
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

if [ "$opt_p_set" = "1" ]; then
    numCores="$opt_p_arg"
else
    numCores=$(grep -c ^processor /proc/cpuinfo)
fi

# echo
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
# echo
# exit 0 

#Check if alignment to reference has been done; if not align sequences to reference
if [ -s $sampleDir/reads.sam ]; then
    echo '**'$sampleId' has already been aligned to '$referenceId
else
    echo '**Align sequence '$sampleId' to reference '$referenceId
    if [ $sampleFilePath2 ]; then
        echo "# "$(date +"%Y-%m-%d %T") bowtie2 -p $numCores --reorder -q -x $referenceBasePath -1 $sampleFilePath1 -2 $sampleFilePath2
        echo "# "$(bowtie2 --version | grep -i -E "bowtie.*version")
        bowtie2 -p $numCores --reorder -q -x $referenceBasePath -1 $sampleFilePath1 -2 $sampleFilePath2 > $sampleDir/'reads.sam'
    else
        echo "# "$(date +"%Y-%m-%d %T") bowtie2 -p $numCores --reorder -q -x $referenceBasePath $sampleFilePath1
        echo "# "$(bowtie2 --version | grep -i -E "bowtie.*version")
        bowtie2 -p $numCores --reorder -q -x $referenceBasePath $sampleFilePath1 > $sampleDir/'reads.sam'
    fi
    echo $(date +"# %Y-%m-%d %T") alignSampleToReference.sh finished
fi
