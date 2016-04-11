#!/bin/bash
#
#Author: Steve Davis (scd)
#  alterations for grid engine by Al Shpuntoff (afs)
#Purpose: 
#    Run the SNP Pipeline on a specified data set.
#Input:
#    reference          : fasta genome reference file
#    samples            : collection of samples, each in a separate directory
#    workDirectory      : directory where the data is copied and results will be generated
#
#    This script expects the following directories and files:
#    <reference name>.fasta
#    <multiple sample subdirectories>/*.fastq
#Output:
#	 If requested, this script mirrors the reference and samples into a new 
#    <workDirectory>.  Within the <workDirectory>, the input files are
#    linked and the outputs are generated.  Many files are generated, but the 
#    most important results are:
#        <workDirectory>/snplist.txt
#            a SNP list identifying the SNPs found across all samples
#        <workDirectory>/snpma.fasta
#            a SNP matrix with one row per sample and one column per SNP
#        <workDirectory>/samples/<multiple sample subdirectories>/reads.snp.pileup
#            one pileup file per sample
#        <workDirectory>/referenceSNP.fasta
#            a fasta file containing the reference sequence bases at all the SNP locations
#Use example:
#
#History:
#   20140910-scd: Started.
#   20141008-scd: Added support for configuration files.
#   20141015-scd: Support soft, hard, copy mirror modes.
#   20141124-afs: Support added for grid engine in addition to torque
#   20150327-scd: Added sample metrics collection and tabulation.
#   20150330-scd: Process sample directories in size order, considering only by the size of fastq files and ignoring all other files.
#   20150331-scd: Don't skip the last sample when run with -S option and the file of directories is not terminated with a newline.
#   20150407-scd: Exit with an error message when run with -s option and the samples directory is empty or contains no subdirectories with fastq files
#   20150408-scd: Don't exit upon du command error when there are no fastq files or no fq files.
#   20150410-scd: Verify each of the sample directories in the sample directory file is not empty and contains fastq's.
#   20150618-scd: Allow trailing slashes in the file of sample directories.
#   20150618-scd: Allow blank lines in the file of sample directories.
#   20150630-scd: Add support for Smalt.
#   20150824-scd: Replace create_snp_pileup.py with call_consensus.py.
#   20150911-scd: Generate the consensus.vcf file for all samples.
#   20150928-scd: Merge the per-sample VCF files into a multi-vcf file.
#   20151013-scd: Added a configuration parameter to control stripping the array jobid suffix.
#   20151210-scd: Trap errors and shutdown the pipeline if configured to do so.
#   20160301-scd: Allow configurable extra paremeters for collectSampleMetrics and combineSampleMetrics.
#   20160304-scd: CollectSampleMetrics : Count missing positions in the consensus.fasta file instead of the snp matrix.
#   20160308-scd: Exclude samples with excessive number of snps from the snp matrix and snpma.vcf.
#   20160405-scd: Add the calculate_snp_distances.py pipeline stage.
#Notes:
#
#Bugs:
#

# Source the utility functions
. snp_pipeline_inc.sh

# Trap errors
set -o pipefail  # Set the exit status of a unix pipe to the exit code of the last program to exit non-zero

# Log a fatal error to the error summary file and to stderr and then exit.
fatalError()
{
    reportError "$1"

    exit 1
}

# This function is automatically called by bash upon any command error.
# It Logs the error and exits the SNP Pipeline if configured to do so.
handleTrappedErrors()
{
    errorCode=$?
    bashCommand="$BASH_COMMAND"
    bashLineNo=${BASH_LINENO[0]}

    # A subprocess failed and was already trapped.
    # Actually, we cannot be 100% certain the error was trapped if the error code is 123.  This
    # indicates an error in an array of sample jobs launched in parallel by xargs; but since 
    # SnpPipeline_StopOnSampleError is true, we will assume the error was already trapped.
    if [[ $errorCode = 100 || $errorCode = 123 && $SnpPipeline_StopOnSampleError = true ]]; then
        logError "See also the log files in directory $logDir"
        logError "Shutting down the SNP Pipeline."
        logError "================================================================================"
        cat "$errorOutputFile" 1>&2
        exit 1
    fi

    # When configured to ignore single-sample errors, try to continue execution if
    # xargs returns 123 when executing an array of jobs, one per sample
    if [[ $errorCode = 123 && $SnpPipeline_StopOnSampleError != true ]]; then
        return 0
    fi

    # Error code 98 indicates an error affecting a single sample when the pipeline is configured to
    # ignore such errors
    #if [[ $errorCode = 98 ]]; then
    #fi

    # An error occured and was not already trapped.  Error code 98 indicates an error
    # affecting a single sample only.
    if [[ $errorCode != 98 ]]; then
        logError "Error detected while running $(basename $0)."
        logError ""
        logError "The command at line $bashLineNo returned error code $errorCode:"
        logError "    $bashCommand"
        logError ""
        logError "Shutting down the SNP Pipeline."
        logError "================================================================================"
        cat "$errorOutputFile" 1>&2
        exit 1
    fi
}

trap handleTrappedErrors ERR


usageshort()
{
    echo "usage: run_snp_pipeline.sh [-h] [-f] [-m MODE] [-c FILE] [-Q \"torque\"|\"grid\"]  [-o DIR]  (-s DIR | -S FILE)  referenceFile"
    echo '  -h for detailed help message'
    echo 
}

usagelong()
{
    echo "usage: run_snp_pipeline.sh [-h] [-f] [-m MODE] [-c FILE] [-Q torque|grid] [-o DIR] (-s DIR|-S FILE)"
    echo "                           referenceFile"
    echo 
    echo 'Run the SNP Pipeline on a specified data set.'
    echo
    echo 'Positional arguments:'
    echo '  referenceFile  : Relative or absolute path to the reference fasta file.'
    echo
    echo 'Options:'
    echo '  -h             : Show this help message and exit.'
    echo
    echo '  -f             : Force processing even when result files already exist and '
    echo '                   are newer than inputs.'
    echo
    echo '  -m MODE        : Create a mirror copy of the reference directory and all the sample '
    echo '                   directories.  Use this option to avoid polluting the reference directory and '
    echo '                   sample directories with the intermediate files generated by the snp pipeline. '
    echo '                   A "reference" subdirectory and a "samples" subdirectory are created under '
    echo '                   the output directory (see the -o option).  One directory per sample is created '
    echo '                   under the "samples" directory.  Three suboptions allow a choice of how the '
    echo '                   reference and samples are mirrored.'
    echo '                     -m soft : creates soft links to the fasta and fastq files instead of copying'
    echo '                     -m hard : creates hard links to the fasta and fastq files instead of copying'
    echo '                     -m copy : copies the fasta and fastq files' 
    echo
    echo '  -c FILE        : Relative or absolute path to a configuration file for overriding defaults '
    echo '                   and defining extra parameters for the tools and scripts within the pipeline. '
    echo '                   Note: A default parameter configuration file named "snppipeline.conf" is '
    echo '                         used whenever the pipeline is run without the -c option.  The '
    echo '                         configuration file used for each run is copied into the log directory, '
    echo '                         capturing the parameters used during the run.'
    echo 
    echo '  -Q torque|grid : Job queue manager for remote parallel job execution in an HPC environment.'
    echo '                   Currently "torque" and "grid" are supported.  If not specified, the pipeline '
    echo '                   will execute locally.'
    echo
    echo '  -o DIR         : Output directory for the snp list, snp matrix, and reference snp files.'
    echo '                   Additional subdirectories are automatically created under the output '
    echo '                   directory for logs files and the mirrored samples and reference files '
    echo '                   (see the -m option).  The output directory will be created if it does '
    echo '                   not already exist.  If not specified, the output files are written to '
    echo '                   the current working directory.  If you re-run the pipeline on previously'
    echo '                   processed samples, and specify a different output directory, the '
    echo '                   pipeline will not rebuild everything unless you either force a rebuild '
    echo '                   (see the -f option) or you request mirrored inputs (see the -m option).'
    echo
    echo '  -s DIRECTORY   : Relative or absolute path to the parent directory of all the sample '
    echo '                   directories.  The -s option should be used when all the sample directories'
    echo '                   are in subdirectories immediately below a parent directory.'
    echo '                   Note: You must specify either the -s or -S option, but not both.'
    echo '                   Note: The specified directory should contain only a collection of sample'
    echo '                         directories, nothing else.'
    echo '                   Note: Unless you request mirrored inputs, see the -m option, additional files'
    echo '                         will be written to each of the sample directories during the execution '
    echo '                         of the SNP Pipeline'
    echo 
    echo '  -S FILE        : Relative or absolute path to a file listing all of the sample directories.'
    echo '                   The -S option should be used when the samples are not under a common parent '
    echo '                   directory.  '
    echo '                   Note: If you are not mirroring the samples (see the -m option), you can'
    echo '                         improve parallel processing performance by sorting the the list of '
    echo '                         directories descending by size, largest first.  The -m option '
    echo '                         automatically generates a sorted directory list.'
    echo '                   Note: You must specify either the -s or -S option, but not both.'
    echo '                   Note: Unless you request mirrored inputs, see the -m option, additional files'
    echo '                         will be written to each of the sample directories during the execution '
    echo '                         of the SNP Pipeline'
    echo
}

get_abs_filename() 
{
  # $1 : relative filename
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}


stripTorqueJobArraySuffix()
{
  if [[ $Torque_StripJobArraySuffix = true ]]; then
    echo ${1%%.*}
  else
    echo $1
  fi
}

stripGridEngineJobArraySuffix()
{
  if [[ $GridEngine_StripJobArraySuffix = true ]]; then
    echo ${1%%.*}
  else
    echo $1
  fi
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

while getopts ":hfm:c:Q:o:s:S:" option; do
    if [ "$option" = "h" ]; then
        usagelong
        exit 0
    elif [ "$option" = "?" ]; then
        echo "Invalid option -- '$OPTARG'"
        echo
        usageshort
        exit 1
    elif [ "$option" = ":" ]; then
        echo "Missing argument for option -- '$OPTARG'"
        echo
        usageshort
        exit 2
    else
        declare opt_"$option"_set="1"
        if [ "$OPTARG" != "" ]; then
            declare opt_"$option"_arg="$OPTARG"
        fi
    fi
done



# Handle output working directory.  Create the directory if it does not exist.
if [[ "$opt_o_set" = "1" ]]; then
    export workDir="$opt_o_arg"
    if ! mkdir -p "$workDir"; then echo "Could not create the output directory $workDir"; exit 50; fi
    if [[ ! -w "$workDir" ]]; then echo "Output directory $workDir is not writable."; exit 50; fi
else
    export workDir="$(pwd)"
fi

# The error log is in the main workdir
export errorOutputFile="$workDir/error.log"
rm  "$errorOutputFile" 2> /dev/null || true


# --------------------------------------------------------
# get the arguments
shift $((OPTIND-1))

# Reference fasta file
export referenceFilePath="$1"
if [ "$referenceFilePath" = "" ]; then
    echo "Missing reference file."
    echo
    usageshort
    exit 10
fi
if [[ ! -f "$referenceFilePath" ]]; then fatalError "Reference file $referenceFilePath does not exist."; fi
if [[ ! -s "$referenceFilePath" ]]; then fatalError "Reference file $referenceFilePath is empty."; fi
export referenceFileName=${referenceFilePath##*/} # strip directories

# Extra arguments not allowed
if [[ "$2" != "" ]]; then 
    echo "Unexpected argument \"$2\" specified after the reference file."
    echo
    usageshort
    exit 20
fi

# Force rebuild flag
if [[ "$opt_f_set" = "1" ]]; then
    export forceFlag="-f"
else
    unset forceFlag
fi

# Mirror copy input files flag
if [[ "$opt_m_set" = "1" ]]; then
    mirrorMode="$opt_m_arg"
    mirrorMode=$(echo "$mirrorMode" | tr [:upper:] [:lower:])
    if [[ "$mirrorMode" != "soft" && "$mirrorMode" != "hard" && "$mirrorMode" != "copy" ]]; then 
        echo "Invalid mirror mode: $mirrorMode"
        echo
        usageshort
        exit 30
    fi
fi

# Job queue manager for remote parallel job execution
if [[ "$opt_Q_set" = "1" ]]; then
    platform=$(echo "$opt_Q_arg" | tr '[:upper:]' '[:lower:]')
    if [[ "$platform" != "torque"  &&  "$platform" != "grid" ]]; then
        echo "Only the torque and grid job queues are currently supported."
        echo
        usageshort
        exit 40
    fi
fi

# Handle sample directories
if [[ "$opt_s_set" = "1" && "$opt_S_set" = "1" ]]; then
    echo "Options -s and -S are mutually exclusive."
    echo
    usageshort
    exit 60
fi
if [[ "$opt_s_set" != "1" && "$opt_S_set" != "1" ]]; then
    echo "You must specify one of the -s or -S options to identify the samples."
    echo
    usageshort
    exit 60
fi


# Create the logs directory
runTimeStamp=$(date +"%Y%m%d.%H%M%S")
export logDir="$workDir/logs-$runTimeStamp"
mkdir -p "$logDir"

# Handle configuration file, use the specified file, or create a default file
if [[ "$opt_c_set" = "1" ]]; then
    configFilePath="$opt_c_arg"
    if [[ ! -f "$configFilePath" ]]; then fatalError "Configuration file $configFilePath does not exist."; fi
    if [[ ! -s "$configFilePath" ]]; then fatalError "Configuration file $configFilePath is empty."; fi
    cp -p "$configFilePath" "$logDir"
    source "$configFilePath"
else
    # true below is to ignore preserve timestamp error
    copy_snppipeline_data.py configurationFile "$logDir" || true
    source "$logDir/snppipeline.conf"
fi

# Validate the configured aligner choice
if [[ "$SnpPipeline_Aligner" == "" ]]; then
    SnpPipeline_Aligner="bowtie2"
else
    # make lowercase 
    SnpPipeline_Aligner=$(echo "$SnpPipeline_Aligner" | tr '[:upper:]' '[:lower:]')
    if [[ "$SnpPipeline_Aligner" != "bowtie2" && "$SnpPipeline_Aligner" != "smalt" ]]; then
        fatalError "Config file error in SnpPipeline_Aligner parameter: only bowtie2 and smalt aligners are supported."
    fi
fi

# Stop the pipeline by default upon single sample errors if not configured either way
if [ -z $SnpPipeline_StopOnSampleError ]; then
    SnpPipeline_StopOnSampleError=true
fi

export SnpPipeline_StopOnSampleError
export SnpPipeline_Aligner
export Bowtie2Build_ExtraParams
export SmaltIndex_ExtraParams
export SamtoolsFaidx_ExtraParams
export Bowtie2Align_ExtraParams
export SmaltAlign_ExtraParams
export SamtoolsSamFilter_ExtraParams
export SamtoolsSort_ExtraParams
export SamtoolsMpileup_ExtraParams
export VarscanMpileup2snp_ExtraParams
export VarscanJvm_ExtraParams
export CreateSnpList_ExtraParams
export CallConsensus_ExtraParams
export CreateSnpMatrix_ExtraParams
export CreateSnpReferenceSeq_ExtraParams
export CollectSampleMetrics_ExtraParams
export CombineSampleMetrics_ExtraParams
export GridEngine_PEname


# Verify the scripts and needed tools were properly installed on the path
(( dependencyErrors = 0 )) || true
onPath=$(verifyOnPath "prepReference.sh"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "alignSampleToReference.sh"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "prepSamples.sh"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "create_snp_list.py"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "call_consensus.py"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "create_snp_matrix.py"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "create_snp_reference_seq.py"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "collectSampleMetrics.sh"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "combineSampleMetrics.sh"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "calculate_snp_distances.py"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "$SnpPipeline_Aligner"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "samtools"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "java"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "tabix"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "bgzip"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "bcftools"); if [[ $onPath != true ]]; then (( dependencyErrors += 1 )); fi
result=$(java net.sf.varscan.VarScan 2>&1) || true
if [[ $result =~ .*Error.* ]]; then
    reportError "CLASSPATH is not configured with the path to VarScan"
    (( dependencyErrors += 1 ))
fi
if (( $dependencyErrors > 0 )); then
    fatalError 'Check the SNP Pipeline installation instructions here: http://snp-pipeline.readthedocs.org/en/latest/installation.html'
fi


# Rewrite the file of sample directories, removing trailing slashes and blank lines.
rewriteCleansedFileOfSampleDirs()
{
    local inSampleDirsFile=$1
    local outSampleDirsFile=$2

    # Is the file of sample dirs the same as the file to be created?
    if [ "$inSampleDirsFile" -ef "$outSampleDirsFile" ]; then
        # Remove trailing slashes and blank lines
        sed --in-place -e "s,/\+$,," -e '/^[[:space:]]*$/d' "$inSampleDirsFile"
    else
        sed            -e "s,/\+$,," -e '/^[[:space:]]*$/d' "$inSampleDirsFile" > "$outSampleDirsFile"
    fi
}

# Verify each of the sample directories in the sample directory file is not empty and contains fastq's
validateFileOfSampleDirs()
{
    local sampleDirsFile=$1
    local foundError=false

    while IFS=$'\n' read -r dir || [[ -n "$dir" ]]
    do 
        if [[ ! -d "$dir" ]]; then 
            reportError "Sample directory $dir does not exist."
            foundError=true
        elif [[ -z $(ls -A "$dir") ]]; then 
            reportError "Sample directory $dir is empty."
            foundError=true
        else
            fastqFiles=$({ find "$dir" -path "$dir"'/*.fastq*'; find "$dir" -path "$dir"'/*.fq*'; })
            if [[ -z "$fastqFiles" ]]; then 
                reportError "Sample directory $dir does not contain any fastq files."
                foundError=true
            fi
        fi
    done  < "$sampleDirsFile" 
    if [[ "$foundError" == true ]]; then 
        if [[ -z $SnpPipeline_StopOnSampleError || $SnpPipeline_StopOnSampleError == true ]]; then
            exit 1
        else
            logError "================================================================================"
        fi
    fi
}


# --------------------------------------------------------
# get sample directories sorted by size, largest first

persistSortedSampleDirs()
{
    local samplesDir=$1
    local workDir=$2

    tmpFile=$(mktemp -p "$workDir" tmp.sampleDirs.XXXXXXXX)
    du -b -L "$samplesDir"/*/*.fastq* 2> /dev/null > "$tmpFile" || true
    du -b -L "$samplesDir"/*/*.fq* 2> /dev/null >> "$tmpFile" || true
    < "$tmpFile" xargs -n 2 sh -c 'echo $0 $(dirname "$1")' | \
    awk '{sizes[$2]+=$1} END {for (dir in sizes) {printf "%s %.0f\n", dir, sizes[dir]}}' | \
    sort -k 2 -n -r | \
    cut -f 1 -d" " > "$workDir/sampleDirectories.txt"  
    rm "$tmpFile"
}


if [[ "$opt_s_set" = "1" ]]; then
    samplesDir="${opt_s_arg%/}"  # strip trailing slash
    if [[ ! -d "$samplesDir" ]]; then fatalError "Samples directory $samplesDir does not exist."; fi
    if [[   -z $(ls -A "$samplesDir") ]]; then fatalError "Samples directory $samplesDir is empty."; fi
    fastqFiles=$({ find "$samplesDir" -path "$samplesDir"'/*/*.fastq*'; find "$samplesDir" -path "$samplesDir"'/*/*.fq*'; })
    if [[   -z "$fastqFiles" ]]; then fatalError "Samples directory $samplesDir does not contain subdirectories with fastq files."; fi
    persistSortedSampleDirs "$samplesDir" "$workDir"
    sampleDirsFile="$workDir/sampleDirectories.txt"
fi
if [[ "$opt_S_set" = "1" ]]; then
    sampleDirsFile="$opt_S_arg"
    if [[ ! -f "$sampleDirsFile" ]]; then fatalError "The file of samples directories, $sampleDirsFile, does not exist."; fi
    if [[ ! -s "$sampleDirsFile" ]]; then fatalError "The file of samples directories, $sampleDirsFile, is empty."; fi
    rewriteCleansedFileOfSampleDirs "$sampleDirsFile" "$workDir/sampleDirectories.txt"
    sampleDirsFile="$workDir/sampleDirectories.txt"
    validateFileOfSampleDirs "$sampleDirsFile"
fi
sampleCount=$(cat "$sampleDirsFile" | wc -l)

# --------------------------------------------------------
# Mirror the input reference and samples if requested
if [[ "$opt_m_set" = "1" ]]; then
    if [[ "$mirrorMode" == "soft" ]]; then
        # soft link, subsequent freshness checks use the timestamp of original file, not the soft link
        mirrorFlag="-s"
    elif [[ "$mirrorMode" == "hard" ]]; then
        # hard link, automatically preserves attributes of the original file
        mirrorFlag="-l"
    else
        # regular copy, -p explicitly preserves attributes of the original file
        mirrorFlag="-p"
    fi

    # Mirror/link the reference
    mkdir -p "$workDir/reference"
    absoluteReferenceFilePath=$(get_abs_filename "$referenceFilePath")
    cp -v -u -f $mirrorFlag "$absoluteReferenceFilePath" "$workDir/reference"
    # since we mirrored the reference, we need to update our reference location
    referenceFileName=${referenceFilePath##*/} # strip directories
    referenceFilePath="$workDir/reference/$referenceFileName"

    # Mirror/link the samples
    cat "$sampleDirsFile" | while IFS=$'\n' read -r dir || [[ -n "$dir" ]]
    do
        baseDir=${dir##*/} # strip the parent directories
        mkdir -p "$workDir/samples/$baseDir"
        # copy without stderr message and without exit error code because the fastq or fq files might not exist
        absoluteSampleDir=$(get_abs_filename "$dir")
        cp -r -v -u -f $mirrorFlag "$absoluteSampleDir"/*.fastq* "$workDir/samples/$baseDir" 2> /dev/null || true
        cp -r -v -u -f $mirrorFlag "$absoluteSampleDir"/*.fq* "$workDir/samples/$baseDir" 2> /dev/null || true
    done
    # since we mirrored the samples, we need to update our samples location and sorted list of samples
    samplesDir="$workDir/samples"
    persistSortedSampleDirs "$samplesDir" "$workDir"
    sampleDirsFile="$workDir/sampleDirectories.txt"
fi


# --------------------------------------------------------
echo -e "\nStep 1 - Prep work"
export numCores=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)
# get the *.fastq or *.fq files in each sample directory, possibly compresessed, on one line per sample, ready to feed to bowtie
tmpFile=$(mktemp -p "$workDir" tmp.fastqs.XXXXXXXX)
cat "$sampleDirsFile" | while IFS=$'\n' read -r dir || [[ -n "$dir" ]]; do echo $dir/*.fastq* >> "$tmpFile"; echo "$dir"/*.fq* >> "$tmpFile"; done
grep -v '*.fq*' "$tmpFile" | grep -v '*.fastq*' > "$workDir/sampleFullPathNames.txt"
rm "$tmpFile"

echo -e "\nStep 2 - Prep the reference"
if [[ "$platform" == "grid" ]]; then
    prepReferenceJobId=$(echo | qsub -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N prepReference
#$ -V
#$ -j y
#$ -cwd
#$ -o $logDir/prepReference.log
    prepReference.sh $forceFlag "$referenceFilePath" 
_EOF_
)
elif [[ "$platform" == "torque" ]]; then
    prepReferenceJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N prepReference
    #PBS -j oe
    #PBS -d $(pwd)
    #PBS -o $logDir/prepReference.log
    #PBS -V
    prepReference.sh $forceFlag "$referenceFilePath"
_EOF_
)
else
    prepReference.sh $forceFlag "$referenceFilePath" 2>&1 | tee $logDir/prepReference.log
fi

echo -e "\nStep 3 - Align the samples to the reference"
# Parse the user-specified aligner parameters to find the number of CPU cores requested, for example, "-p 16" or "-n 16"
# Set the default number of  CPU cores if the user did not configure a value.
if [[ "$platform" == "grid" || "$platform" == "torque" ]]; then
    if [[ "$SnpPipeline_Aligner" == "smalt" ]]; then
        regex="(-n[[:blank:]]*)([[:digit:]]+)"
        if [[ "$SmaltAlign_ExtraParams" =~ $regex ]]; then
            numAlignThreads=${BASH_REMATCH[2]}
        else
            numAlignThreads=8
            SmaltAlign_ExtraParams="$SmaltAlign_ExtraParams -n $numAlignThreads"
        fi
    else
        regex="(-p[[:blank:]]*)([[:digit:]]+)"
        if [[ "$Bowtie2Align_ExtraParams" =~ $regex ]]; then
            numAlignThreads=${BASH_REMATCH[2]}
        else
            numAlignThreads=8
            Bowtie2Align_ExtraParams="$Bowtie2Align_ExtraParams -p $numAlignThreads"
        fi
    fi
fi

if [[ "$platform" == "grid" ]]; then
    alignSamplesJobId=$(echo | qsub -terse -t 1-$sampleCount $GridEngine_QsubExtraParams << _EOF_
#$   -N alignSamples
#$   -cwd
#$   -V
#$   -j y
#$   -pe $GridEngine_PEname $numAlignThreads
#$   -hold_jid $prepReferenceJobId
#$   -o $logDir/alignSamples.log-\$TASK_ID
    alignSampleToReference.sh $forceFlag "$referenceFilePath" \$(cat "$workDir/sampleFullPathNames.txt" | head -n \$SGE_TASK_ID | tail -n 1)
_EOF_
)
elif [[ "$platform" == "torque" ]]; then
    alignSamplesJobId=$(echo | qsub -t 1-$sampleCount $Torque_QsubExtraParams << _EOF_
    #PBS -N alignSamples
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -l nodes=1:ppn=$numAlignThreads
    #PBS -W depend=afterok:$prepReferenceJobId
    #PBS -o $logDir/alignSamples.log
    #PBS -V
    samplesToAlign=\$(cat "$workDir/sampleFullPathNames.txt" | head -n \$PBS_ARRAYID | tail -n 1)
    alignSampleToReference.sh $forceFlag "$referenceFilePath" \$samplesToAlign
_EOF_
)
else
    nl "$workDir/sampleFullPathNames.txt" | xargs -n 3 -L 1 bash -c 'set -o pipefail; alignSampleToReference.sh $forceFlag "$referenceFilePath" $1 $2 2>&1 | tee $logDir/alignSamples.log-$0'
fi

echo -e "\nStep 4 - Prep the samples"
if [[ "$platform" == "grid" ]]; then
    sleep $((1 + sampleCount / 150)) # workaround potential bug when submitting two large consecutive array jobs
    alignSamplesJobArray=$(stripGridEngineJobArraySuffix $alignSamplesJobId)
    prepSamplesJobId=$(echo | qsub -terse -t 1-$sampleCount $GridEngine_QsubExtraParams << _EOF_
#$   -N prepSamples
#$   -cwd
#$   -V
#$   -j y
#$   -hold_jid_ad $alignSamplesJobArray
#$   -l h_rt=05:00:00
#$   -o $logDir/prepSamples.log-\$TASK_ID
    prepSamples.sh $forceFlag "$referenceFilePath" "\$(cat "$sampleDirsFile" | head -n \$SGE_TASK_ID | tail -n 1)"
_EOF_
)
elif [[ "$platform" == "torque" ]]; then
    sleep $((1 + sampleCount / 150)) # workaround torque bug when submitting two large consecutive array jobs
    alignSamplesJobArray=$(stripTorqueJobArraySuffix $alignSamplesJobId)
    prepSamplesJobId=$(echo | qsub -t 1-$sampleCount $Torque_QsubExtraParams << _EOF_
    #PBS -N prepSamples
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterokarray:$alignSamplesJobArray
    #PBS -l walltime=05:00:00
    #PBS -o $logDir/prepSamples.log
    #PBS -V
    sampleDir=\$(cat "$sampleDirsFile" | head -n \$PBS_ARRAYID | tail -n 1)
    prepSamples.sh $forceFlag "$referenceFilePath" "\$sampleDir"
_EOF_
)
else
    if [[ "$MaxConcurrentPrepSamples" != "" ]]; then
        numPrepSamplesCores=$MaxConcurrentPrepSamples
    else
        numPrepSamplesCores=$numCores
    fi
    nl "$sampleDirsFile" | xargs -n 2 -P $numPrepSamplesCores bash -c 'set -o pipefail; prepSamples.sh $forceFlag "$referenceFilePath" $1 2>&1 | tee $logDir/prepSamples.log-$0'
fi

echo -e "\nStep 5 - Combine the SNP positions across all samples into the SNP list file"
if [[ "$platform" == "grid" ]]; then
    prepSamplesJobArray=$(stripGridEngineJobArraySuffix $prepSamplesJobId)
    snpListJobId=$(echo | qsub  -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N snpList
#$ -cwd
#$ -j y
#$ -V
#$ -hold_jid $prepSamplesJobArray
#$ -o $logDir/snpList.log
    create_snp_list.py $forceFlag -n var.flt.vcf -o "$workDir/snplist.txt" $CreateSnpList_ExtraParams "$sampleDirsFile" 
_EOF_
)
elif [[ "$platform" == "torque" ]]; then
    prepSamplesJobArray=$(stripTorqueJobArraySuffix $prepSamplesJobId)
    snpListJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N snpList
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterokarray:$prepSamplesJobArray
    #PBS -o $logDir/snpList.log
    #PBS -V
    create_snp_list.py $forceFlag -n var.flt.vcf -o "$workDir/snplist.txt" $CreateSnpList_ExtraParams "$sampleDirsFile" 
_EOF_
)
else
    create_snp_list.py $forceFlag -n var.flt.vcf -o "$workDir/snplist.txt" $CreateSnpList_ExtraParams "$sampleDirsFile" 2>&1 | tee $logDir/snpList.log
fi

# The create_snp_list process creates the filtered list of sample directories.  It is the list of samples having removed the samples with excessive snps.
# When running on a workstation, the file exists at this point during the script execution, but on grid or torque, it has not yet been created. However,
# we know the path to the file regardless of whether it exists yet.
filteredSampleDirsFile="${sampleDirsFile}.filtered"


echo -e "\nStep 6 - Call the consensus SNPs for each sample"
if [[ "$platform" == "grid" ]]; then
    callConsensusJobId=$(echo | qsub -terse -t 1-$sampleCount $GridEngine_QsubExtraParams << _EOF_
#$ -N callConsensus
#$ -cwd
#$ -V
#$ -j y
#$ -hold_jid $snpListJobId
#$ -o $logDir/callConsensus.log-\$TASK_ID
    sampleDir=\$(cat "$sampleDirsFile" | head -n \$SGE_TASK_ID | tail -n 1)
    call_consensus.py $forceFlag -l "$workDir/snplist.txt" -o "\$sampleDir/consensus.fasta" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams "\$sampleDir/reads.all.pileup"
_EOF_
)
elif [[ "$platform" == "torque" ]]; then
    callConsensusJobId=$(echo | qsub -t 1-$sampleCount $Torque_QsubExtraParams << _EOF_
    #PBS -N callConsensus
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterok:$snpListJobId
    #PBS -o $logDir/callConsensus.log
    #PBS -V
    sampleDir=\$(cat "$sampleDirsFile" | head -n \$PBS_ARRAYID | tail -n 1)
    call_consensus.py $forceFlag -l "$workDir/snplist.txt" -o "\$sampleDir/consensus.fasta" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams "\$sampleDir/reads.all.pileup"
_EOF_
)
else
    if [[ "$MaxConcurrentCallConsensus" != "" ]]; then
        numCallConsensusCores=$MaxConcurrentCallConsensus
    else
        numCallConsensusCores=$numCores
    fi
    nl "$sampleDirsFile" | xargs -n 2 -P $numCallConsensusCores bash -c 'set -o pipefail; call_consensus.py $forceFlag -l "$workDir/snplist.txt" -o "$1/consensus.fasta" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams "$1/reads.all.pileup" 2>&1 | tee $logDir/callConsensus.log-$0'
fi

echo -e "\nStep 7 - Create the SNP matrix"
if [[ "$platform" == "grid" ]]; then
    callConsensusJobArray=$(stripGridEngineJobArraySuffix $callConsensusJobId)
    snpMatrixJobId=$(echo | qsub -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N snpMatrix
#$ -cwd
#$ -V
#$ -j y
#$ -hold_jid $callConsensusJobArray
#$ -l h_rt=05:00:00
#$ -o $logDir/snpMatrix.log
    create_snp_matrix.py $forceFlag -c consensus.fasta -o "$workDir/snpma.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile"
_EOF_
)
elif [[ "$platform" == "torque" ]]; then
    callConsensusJobArray=$(stripTorqueJobArraySuffix $callConsensusJobId)
    snpMatrixJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N snpMatrix
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterokarray:$callConsensusJobArray
    #PBS -l walltime=05:00:00
    #PBS -o $logDir/snpMatrix.log
    #PBS -V
    create_snp_matrix.py $forceFlag -c consensus.fasta -o "$workDir/snpma.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile"
_EOF_
)
else
    create_snp_matrix.py $forceFlag -c consensus.fasta -o "$workDir/snpma.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile" 2>&1 | tee $logDir/snpMatrix.log
fi    

echo -e "\nStep 8 - Create the reference base sequence"
if [[ "$platform" == "grid" ]]; then
    snpReferenceJobId=$(echo | qsub -terse $GridEngine_QsubExtraParams << _EOF_
#$ -V
#$ -N snpReference 
#$ -cwd
#$ -j y 
#$ -hold_jid $callConsensusJobArray
#$ -o $logDir/snpReference.log
    create_snp_reference_seq.py $forceFlag -l "$workDir/snplist.txt" -o "$workDir/referenceSNP.fasta" $CreateSnpReferenceSeq_ExtraParams "$referenceFilePath"
_EOF_
)
elif [[ "$platform" == "torque" ]]; then
    snpReferenceJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N snpReference 
    #PBS -d $(pwd)
    #PBS -j oe 
    #PBS -W depend=afterokarray:$callConsensusJobArray
    #PBS -o $logDir/snpReference.log
    #PBS -V
    create_snp_reference_seq.py $forceFlag -l "$workDir/snplist.txt" -o "$workDir/referenceSNP.fasta" $CreateSnpReferenceSeq_ExtraParams "$referenceFilePath"
_EOF_
)
else
    create_snp_reference_seq.py $forceFlag -l "$workDir/snplist.txt" -o "$workDir/referenceSNP.fasta" $CreateSnpReferenceSeq_ExtraParams "$referenceFilePath" 2>&1 | tee $logDir/snpReference.log
fi


echo -e "\nStep 9 - Create the Multi-VCF file"
if [[ $CallConsensus_ExtraParams =~ .*vcfFileName.* ]]; then
    if [[ "$platform" == "grid" ]]; then
        mergeVcfJobId=$(echo | qsub  -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N mergeVcf
#$ -cwd
#$ -j y
#$ -V
#$ -hold_jid $callConsensusJobArray
#$ -o $logDir/mergeVcf.log
        mergeVcf.sh $forceFlag -o "$workDir/snpma.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile"
_EOF_
)
    elif [[ "$platform" == "torque" ]]; then
        mergeVcfJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
        #PBS -N mergeVcf
        #PBS -d $(pwd)
        #PBS -j oe
        #PBS -W depend=afterokarray:$callConsensusJobArray
        #PBS -o $logDir/mergeVcf.log
        #PBS -V
        mergeVcf.sh $forceFlag -o "$workDir/snpma.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile"
_EOF_
)
    else
        mergeVcf.sh $forceFlag -o "$workDir/snpma.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile" 2>&1 | tee $logDir/mergeVcf.log
    fi
else
    echo -e "Skipped per CallConsensus_ExtraParams configuration"
fi

echo -e "\nStep 10 - Collect metrics for each sample"
if [[ "$platform" == "grid" ]]; then
    collectSampleMetricsJobId=$(echo | qsub -terse -t 1-$sampleCount $GridEngine_QsubExtraParams << _EOF_
#$ -N collectMetrics
#$ -cwd
#$ -V
#$ -j y
#$ -hold_jid_ad $callConsensusJobArray
#$ -l h_rt=02:00:00
#$ -o $logDir/collectSampleMetrics.log-\$TASK_ID
    sampleDir=\$(cat "$sampleDirsFile" | head -n \$SGE_TASK_ID | tail -n 1)
    collectSampleMetrics.sh -o "\$sampleDir/metrics" $CollectSampleMetrics_ExtraParams "\$sampleDir"  "$referenceFilePath"
_EOF_
)
elif [[ "$platform" == "torque" ]]; then
    collectSampleMetricsJobId=$(echo | qsub -t 1-$sampleCount $Torque_QsubExtraParams << _EOF_
    #PBS -N collectMetrics
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterokarray:$callConsensusJobArray
    #PBS -l walltime=02:00:00
    #PBS -o $logDir/collectSampleMetrics.log
    #PBS -V
    sampleDir=\$(cat "$sampleDirsFile" | head -n \$PBS_ARRAYID | tail -n 1)
    collectSampleMetrics.sh -o "\$sampleDir/metrics" $CollectSampleMetrics_ExtraParams "\$sampleDir"  "$referenceFilePath"
_EOF_
)
else
    if [[ "$MaxConcurrentCollectSampleMetrics" != "" ]]; then
        numCollectSampleMetricsCores=$MaxConcurrentCollectSampleMetrics
    else
        numCollectSampleMetricsCores=$numCores
    fi
    nl "$sampleDirsFile" | xargs -n 2 -P $numCollectSampleMetricsCores bash -c 'set -o pipefail; collectSampleMetrics.sh -o "$1/metrics" $CollectSampleMetrics_ExtraParams "$1" "$referenceFilePath" 2>&1 | tee $logDir/collectSampleMetrics.log-$0'
fi

echo -e "\nStep 11 - Combine the metrics across all samples into the metrics table"
if [[ "$platform" == "grid" ]]; then
    collectSampleMetricsJobArray=$(stripGridEngineJobArraySuffix $collectSampleMetricsJobId)
    combineSampleMetricsJobId=$(echo | qsub  -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N combineMetrics
#$ -cwd
#$ -j y
#$ -V
#$ -hold_jid $collectSampleMetricsJobArray
#$ -o $logDir/combineSampleMetrics.log
    combineSampleMetrics.sh -n metrics -o "$workDir/metrics.tsv" $CombineSampleMetrics_ExtraParams "$sampleDirsFile"
_EOF_
)
elif [[ "$platform" == "torque" ]]; then
    collectSampleMetricsJobArray=$(stripTorqueJobArraySuffix $collectSampleMetricsJobId)
    combineSampleMetricsJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N combineMetrics
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterokarray:$collectSampleMetricsJobArray
    #PBS -o $logDir/combineSampleMetrics.log
    #PBS -V
    combineSampleMetrics.sh -n metrics -o "$workDir/metrics.tsv" $CombineSampleMetrics_ExtraParams "$sampleDirsFile"
_EOF_
)
else
    combineSampleMetrics.sh -n metrics -o "$workDir/metrics.tsv" $CombineSampleMetrics_ExtraParams "$sampleDirsFile" 2>&1 | tee $logDir/combineSampleMetrics.log
fi


echo -e "\nStep 12 - Calculate SNP distance matrix"
if [[ "$platform" == "grid" ]]; then
    calcSnpDistanceJobId=$(echo | qsub  -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N snpDistance
#$ -cwd
#$ -j y
#$ -V
#$ -hold_jid $snpMatrixJobId
#$ -o $logDir/calcSnpDistances.log
    calculate_snp_distances.py $forceFlag -p "$workDir/snp_distance_pairwise.tsv" -m "$workDir/snp_distance_matrix.tsv" "$workDir/snpma.fasta"
_EOF_
)
elif [[ "$platform" == "torque" ]]; then
    calcSnpDistanceJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N snpDistance
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterok:$snpMatrixJobId
    #PBS -o $logDir/calcSnpDistances.log
    #PBS -V
    calculate_snp_distances.py $forceFlag -p "$workDir/snp_distance_pairwise.tsv" -m "$workDir/snp_distance_matrix.tsv" "$workDir/snpma.fasta"
_EOF_
)
else
    calculate_snp_distances.py $forceFlag -p "$workDir/snp_distance_pairwise.tsv" -m "$workDir/snp_distance_matrix.tsv" "$workDir/snpma.fasta" 2>&1 | tee $logDir/calcSnpDistances.log
fi


# Step 13 - Notify user of any non-fatal errors accumulated during processing
if [[ -s "$errorOutputFile" && $SnpPipeline_StopOnSampleError != true ]]; then
    echo "" 1>&2
    echo "There were errors processing some samples." 1>&2
    echo "See the log file $errorOutputFile for a summary of errors." 1>&2
fi

