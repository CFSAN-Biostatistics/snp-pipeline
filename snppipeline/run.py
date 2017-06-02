#!/usr/bin/env python

"""This module is part of the CFSAN SNP Pipeline. It contains the code to
run all the sequential steps of the pipeline in the correct order.
"""

from __future__ import print_function
from __future__ import absolute_import

import os
import shutil
import subprocess
import sys
import traceback

from snppipeline.job_runner import JobRunner
from snppipeline.job_runner import JobRunnerException
from snppipeline import utils
from snppipeline.utils import log_error


# Module globals
log_dir = None


def handle_called_process_exception(exc_type, exc_value, exc_traceback):
    """This function handles exceptions in the child processes executed by
    the snp-pipeline.

    It Logs the error and exits the SNP Pipeline if configured to do so.
    """
    global log_dir

    external_program_command = exc_value.cmd
    error_code = exc_value.returncode
    #print("error_code =", error_code)

    trace_entries = traceback.extract_tb(exc_traceback)
    file_name, line_number, function_name, code_text = trace_entries[-1]
    exc_type_name = exc_type.__name__

    stop_on_error_env = os.environ.get("SnpPipeline_StopOnSampleError")
    stop_on_error = stop_on_error_env is None or stop_on_error_env == "true"

    error_output_file = os.environ.get("errorOutputFile")

    # Error code 98 indicates an error affecting a single sample when the pipeline is configured to ignore such errors.
    # This error has already been reported, but it should not stop execution.
    if error_code == 98 and not stop_on_error:
        return

    # When configured to ignore single-sample errors, try to continue execution if
    # xargs returns 123 when executing an array of jobs, one per sample
    if error_code == 123 and not stop_on_error:
        return # TODO: where does execution continue after returning?

    # A subprocess failed and was already trapped.
    # Actually, we cannot be 100% certain the error was trapped if the error code is 123.  This
    # indicates an error in an array of sample jobs launched in parallel by xargs; but since
    # SnpPipeline_StopOnSampleError is true, we will assume the error was already trapped.
    if stop_on_error and (error_code == 100 or error_code == 123):
        log_error("See also the log files in directory " + log_dir)
        log_error("Shutting down the SNP Pipeline.")
        log_error("================================================================================")

        # Send the error log contents to stderr
        if error_output_file:
            with open(error_output_file, "r") as err_log:
                shutil.copyfileobj(err_log, sys.stderr)
        sys.exit(1)

    # An error occured and was not already trapped.
    # Error code 98 indicates an already handled error affecting a single sample when the pipeline is configured to ignore such errors
    if error_code != 98:
        external_program_command = external_program_command.replace("set -o pipefail; ", "")
        log_error("Error detected while running " + utils.program_name_with_command())
        log_error("")
        log_error("The error occured while running:")
        log_error("    %s" % external_program_command)
        log_error("")
        log_error("Shutting down the SNP Pipeline.")
        log_error("================================================================================")

        # Send the error log contents to stderr
        if error_output_file:
            with open(error_output_file, "r") as err_log:
                shutil.copyfileobj(err_log, sys.stderr)
        sys.exit(1)


def handle_internal_exception(exc_type, exc_value, exc_traceback):
    """This function handles exceptions in the snp-pipeline main driver code -- the code
    that sequentially executes all the sub-tasks. It does not handle exceptions in the
    child processes executed by the snp-pipeline.
    """
    if exc_type == JobRunnerException:
        utils.report_error(str(exc_value))
    else:
        #r = default_exception_handler(exc_type, exc_value, exc_traceback)
        trace_entries = traceback.extract_tb(exc_traceback)
        file_name, line_number, function_name, code_text = trace_entries[-1]
        exc_type_name = exc_type.__name__
        utils.report_error("Exception in function %s at line %d in file %s" % (function_name, line_number, file_name))
        utils.report_error("    %s" % code_text)
        utils.report_error("%s: %s" % (exc_type_name, str(exc_value)))

    utils.report_error("\nShutting down the SNP Pipeline.\n" + 80 * '=')
    sys.exit(1)


def handle_exception(exc_type, exc_value, exc_traceback):
    """This function replaces the default python unhandled exception handler.
    It distinguishes between exceptions in child processes and exceptions
    in the snp-pipeline driver itself.
    """
    if exc_type == subprocess.CalledProcessError:
        handle_called_process_exception(exc_type, exc_value, exc_traceback)
    else:
        handle_internal_exception(exc_type, exc_value, exc_traceback)


def run():
    error_output_file = os.environ.get("errorOutputFile")

    stop_on_error_env = os.environ.get("SnpPipeline_StopOnSampleError")
    stop_on_error = stop_on_error_env is None or stop_on_error_env == "true"


"""
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
  if $Torque_StripJobArraySuffix = true:
    echo ${1%%.*}
  else
    echo $1
  fi
}

stripGridEngineJobArraySuffix()
{
  if $GridEngine_StripJobArraySuffix = true:
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
if "$opt_o_set" = "1":
    export workDir="$opt_o_arg"
    if ! mkdir -p "$workDir"; then echo "Could not create the output directory $workDir"; exit 50; fi
    if ! -w "$workDir": echo "Output directory $workDir is not writable."; exit 50; fi
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
if ! -f "$referenceFilePath": fatalError "Reference file $referenceFilePath does not exist."; fi
if ! -s "$referenceFilePath": fatalError "Reference file $referenceFilePath is empty."; fi
export referenceFileName=${referenceFilePath##*/} # strip directories

# Extra arguments not allowed
if "$2" != "":
    echo "Unexpected argument \"$2\" specified after the reference file."
    echo
    usageshort
    exit 20
fi

# Force rebuild flag
if "$opt_f_set" = "1":
    export forceFlag="-f"
else
    unset forceFlag
fi

# Mirror copy input files flag
if "$opt_m_set" = "1":
    mirrorMode="$opt_m_arg"
    mirrorMode=$(echo "$mirrorMode" | tr [:upper:] [:lower:])
    if "$mirrorMode" != "soft" && "$mirrorMode" != "hard" && "$mirrorMode" != "copy":
        echo "Invalid mirror mode: $mirrorMode"
        echo
        usageshort
        exit 30
    fi
fi

# Job queue manager for remote parallel job execution
if "$opt_Q_set" = "1":
    platform=$(echo "$opt_Q_arg" | tr '[:upper:]' '[:lower:]')
    if platform != "torque"  &&  platform != "grid":
        echo "Only the torque and grid job queues are currently supported."
        echo
        usageshort
        exit 40
    fi
fi

# Handle sample directories
if "$opt_s_set" = "1" && "$opt_S_set" = "1":
    echo "Options -s and -S are mutually exclusive."
    echo
    usageshort
    exit 60
fi
if "$opt_s_set" != "1" && "$opt_S_set" != "1":
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
if "$opt_c_set" = "1":
    configFilePath="$opt_c_arg"
    if ! -f "$configFilePath": fatalError "Configuration file $configFilePath does not exist."; fi
    if ! -s "$configFilePath": fatalError "Configuration file $configFilePath is empty."; fi
    cp -p "$configFilePath" "$logDir"
    source "$configFilePath"
else
    # true below is to ignore preserve timestamp error
    cfsan_snp_pipeline data configurationFile "$logDir" || true
    source "$logDir/snppipeline.conf"
fi

# Validate the configured aligner choice
if "$SnpPipeline_Aligner" == "":
    SnpPipeline_Aligner="bowtie2"
else
    # make lowercase
    SnpPipeline_Aligner=$(echo "$SnpPipeline_Aligner" | tr '[:upper:]' '[:lower:]')
    if "$SnpPipeline_Aligner" != "bowtie2" && "$SnpPipeline_Aligner" != "smalt":
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
export SnpPipeline_RemoveDuplicateReads
export PicardMarkDuplicates_ExtraParams
export PicardJvm_ExtraParams
export SamtoolsMpileup_ExtraParams
export VarscanMpileup2snp_ExtraParams
export VarscanJvm_ExtraParams
export RemoveAbnormalSnp_ExtraParams
export CreateSnpList_ExtraParams
export CallConsensus_ExtraParams
export CreateSnpMatrix_ExtraParams
export BcftoolsMerge_ExtraParams
export CreateSnpReferenceSeq_ExtraParams
export CollectSampleMetrics_ExtraParams
export CombineSampleMetrics_ExtraParams
export GridEngine_PEname


# Verify the scripts and needed tools were properly installed on the path
(( dependencyErrors = 0 )) || true
onPath=$(verifyOnPath "cfsan_snp_pipeline"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "prepReference.sh"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "prepSamples.sh"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "snp_filter.py"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "create_snp_list.py"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "call_consensus.py"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "create_snp_matrix.py"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "create_snp_reference_seq.py"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "collectSampleMetrics.sh"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "combineSampleMetrics.sh"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "calculate_snp_distances.py"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "$SnpPipeline_Aligner"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "samtools"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "java"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "tabix"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "bgzip"); if $onPath != true: (( dependencyErrors += 1 )); fi
onPath=$(verifyOnPath "bcftools"); if $onPath != true: (( dependencyErrors += 1 )); fi
result=$(java net.sf.varscan.VarScan 2>&1) || true
if $result =~ .*Error.*:
    reportError "CLASSPATH is not configured with the path to VarScan"
    (( dependencyErrors += 1 ))
fi
if -z $SnpPipeline_RemoveDuplicateReads || $SnpPipeline_RemoveDuplicateReads = true:
    result=$(java picard.cmdline.PicardCommandLine 2>&1) || true
    if $result =~ .*Error.*:
        reportError "CLASSPATH is not configured with the path to Picard"
        (( dependencyErrors += 1 ))
    fi
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
        if ! -d "$dir":
            reportError "Sample directory $dir does not exist."
            foundError=true
        elif -z $(ls -A "$dir"):
            reportError "Sample directory $dir is empty."
            foundError=true
        else
            fastqFiles=$({ find "$dir" -path "$dir"'/*.fastq*'; find "$dir" -path "$dir"'/*.fq*'; })
            if -z "$fastqFiles":
                reportError "Sample directory $dir does not contain any fastq files."
                foundError=true
            fi
        fi
    done  < sample_dirs_file
    if "$foundError" == true:
        if -z $SnpPipeline_StopOnSampleError || $SnpPipeline_StopOnSampleError == true:
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


# TODO: detect broken fastq symlinks
if "$opt_s_set" = "1":
    samplesDir="${opt_s_arg%/}"  # strip trailing slash
    if ! -d "$samplesDir": fatalError "Samples directory $samplesDir does not exist."; fi
    if   -z $(ls -A "$samplesDir"): fatalError "Samples directory $samplesDir is empty."; fi
    fastqFiles=$({ find "$samplesDir" -path "$samplesDir"'/*/*.fastq*'; find "$samplesDir" -path "$samplesDir"'/*/*.fq*'; })
    if   -z "$fastqFiles": fatalError "Samples directory $samplesDir does not contain subdirectories with fastq files."; fi
    persistSortedSampleDirs "$samplesDir" "$workDir"
    sampleDirsFile="$workDir/sampleDirectories.txt"
fi
if "$opt_S_set" = "1":
    sampleDirsFile="$opt_S_arg"
    if ! -f sample_dirs_file: fatalError "The file of samples directories, $sampleDirsFile, does not exist."; fi
    if ! -s sample_dirs_file: fatalError "The file of samples directories, $sampleDirsFile, is empty."; fi
    rewriteCleansedFileOfSampleDirs sample_dirs_file "$workDir/sampleDirectories.txt"
    sampleDirsFile="$workDir/sampleDirectories.txt"
    validateFileOfSampleDirs sample_dirs_file
fi
sampleCount=$(cat sample_dirs_file | wc -l)

# --------------------------------------------------------
# Mirror the input reference and samples if requested
if "$opt_m_set" = "1":
    if "$mirrorMode" == "soft":
        # soft link, subsequent freshness checks use the timestamp of original file, not the soft link
        mirrorFlag="-s"
    elif "$mirrorMode" == "hard":
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
    cat sample_dirs_file | while IFS=$'\n' read -r dir || [[ -n "$dir" ]]
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
cat sample_dirs_file | while IFS=$'\n' read -r dir || [[ -n "$dir" ]]; do echo $dir/*.fastq* >> "$tmpFile"; echo "$dir"/*.fq* >> "$tmpFile"; done
grep -v '*.fq*' "$tmpFile" | grep -v '*.fastq*' > "$workDir/sampleFullPathNames.txt"
rm "$tmpFile"

echo -e "\nStep 2 - Prep the reference"
if platform == "grid":
    prepReferenceJobId=$(echo | qsub -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N prepReference
#$ -V
#$ -j y
#$ -cwd
#$ -o $logDir/prepReference.log
    cfsan_snp_pipeline index_ref $forceFlag "$referenceFilePath"
_EOF_
)
elif platform == "torque":
    prepReferenceJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N prepReference
    #PBS -j oe
    #PBS -d $(pwd)
    #PBS -o $logDir/prepReference.log
    #PBS -V
    cfsan_snp_pipeline index_ref $forceFlag "$referenceFilePath"
_EOF_
)
else
    cfsan_snp_pipeline index_ref $forceFlag "$referenceFilePath" 2>&1 | tee $logDir/prepReference.log
fi

echo -e "\nStep 3 - Align the samples to the reference"
# Parse the user-specified aligner parameters to find the number of CPU cores requested, for example, "-p 16" or "-n 16"
# Set the default number of  CPU cores if the user did not configure a value.
if platform == "grid" || platform == "torque":
    if "$SnpPipeline_Aligner" == "smalt":
        regex="(-n[[:blank:]]*)([[:digit:]]+)"
        if "$SmaltAlign_ExtraParams" =~ $regex:
            numAlignThreads=${BASH_REMATCH[2]}
        else
            numAlignThreads=8
            SmaltAlign_ExtraParams="$SmaltAlign_ExtraParams -n $numAlignThreads"
        fi
    else
        regex="(-p[[:blank:]]*)([[:digit:]]+)"
        if "$Bowtie2Align_ExtraParams" =~ $regex:
            numAlignThreads=${BASH_REMATCH[2]}
        else
            numAlignThreads=8
            Bowtie2Align_ExtraParams="$Bowtie2Align_ExtraParams -p $numAlignThreads"
        fi
    fi
fi

if platform == "grid":
    alignSamplesJobId=$(echo | qsub -terse -t 1-$sampleCount $GridEngine_QsubExtraParams << _EOF_
#$   -N alignSamples
#$   -cwd
#$   -V
#$   -j y
#$   -pe $GridEngine_PEname $numAlignThreads
#$   -hold_jid $prepReferenceJobId
#$   -o $logDir/alignSamples.log-\$TASK_ID
    cfsan_snp_pipeline map_reads $forceFlag "$referenceFilePath" \$(cat "$workDir/sampleFullPathNames.txt" | head -n \$SGE_TASK_ID | tail -n 1)
_EOF_
)
elif platform == "torque":
    alignSamplesJobId=$(echo | qsub -t 1-$sampleCount $Torque_QsubExtraParams << _EOF_
    #PBS -N alignSamples
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -l nodes=1:ppn=$numAlignThreads
    #PBS -W depend=afterok:$prepReferenceJobId
    #PBS -o $logDir/alignSamples.log
    #PBS -V
    samplesToAlign=\$(cat "$workDir/sampleFullPathNames.txt" | head -n \$PBS_ARRAYID | tail -n 1)
    cfsan_snp_pipeline map_reads $forceFlag "$referenceFilePath" \$samplesToAlign
_EOF_
)
else
    nl "$workDir/sampleFullPathNames.txt" | xargs -n 3 -L 1 bash -c 'set -o pipefail; cfsan_snp_pipeline map_reads $forceFlag "$referenceFilePath" $1 $2 2>&1 | tee $logDir/alignSamples.log-$0'
fi

echo -e "\nStep 4 - Prep the samples"
if platform == "grid":
    sleep $((1 + sampleCount / 150)) # workaround potential bug when submitting two large consecutive array jobs
    alignSamplesJobArray=$(stripGridEngineJobArraySuffix $alignSamplesJobId)
    prepSamplesJobId=$(echo | qsub -terse -t 1-$sampleCount $GridEngine_QsubExtraParams << _EOF_
#$   -N prepSamples
#$   -cwd
#$   -V
#$   -j y
#$   -hold_jid_ad $alignSamplesJobArray
#$   -o $logDir/prepSamples.log-\$TASK_ID
    cfsan_snp_pipeline call_sites $forceFlag "$referenceFilePath" "\$(cat sample_dirs_file | head -n \$SGE_TASK_ID | tail -n 1)"
_EOF_
)
elif platform == "torque":
    sleep $((1 + sampleCount / 150)) # workaround torque bug when submitting two large consecutive array jobs
    alignSamplesJobArray=$(stripTorqueJobArraySuffix $alignSamplesJobId)
    prepSamplesJobId=$(echo | qsub -t 1-$sampleCount $Torque_QsubExtraParams << _EOF_
    #PBS -N prepSamples
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterokarray:$alignSamplesJobArray
    #PBS -o $logDir/prepSamples.log
    #PBS -V
    sampleDir=\$(cat sample_dirs_file | head -n \$PBS_ARRAYID | tail -n 1)
    cfsan_snp_pipeline call_sites $forceFlag "$referenceFilePath" "\$sampleDir"
_EOF_
)
else
    if "$MaxConcurrentPrepSamples" != "":
        numPrepSamplesCores=$MaxConcurrentPrepSamples
    else
        numPrepSamplesCores=$numCores
    fi
    nl sample_dirs_file | xargs -n 2 -P $numPrepSamplesCores bash -c 'set -o pipefail; cfsan_snp_pipeline call_sites $forceFlag "$referenceFilePath" $1 2>&1 | tee $logDir/prepSamples.log-$0'
fi

#Filter abnormal SNPs if needed
echo -e "\nStep 5 - Remove abnormal SNPs"
if platform == "grid":
	prepSamplesJobArray=$(stripGridEngineJobArraySuffix $prepSamplesJobId)
	filterAbnSNPJobId=$(echo | qsub  -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N snpFiltering
#$ -cwd
#$ -j y
#$ -V
#$ -hold_jid $prepSamplesJobArray
#$ -o $logDir/filterAbnormalSNP.log
cfsan_snp_pipeline filter_regions -n var.flt.vcf sample_dirs_file "$referenceFilePath" $RemoveAbnormalSnp_ExtraParams
_EOF_
)
elif platform == "torque":
	prepSamplesJobArray=$(stripTorqueJobArraySuffix $prepSamplesJobId)
	filterAbnSNPJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
	#PBS -N snpFiltering
	#PBS -d $(pwd)
	#PBS -j oe
	#PBS -W depend=afterokarray:$prepSamplesJobArray
	#PBS -o $logDir/filterAbnormalSNP.log
	#PBS -V
	cfsan_snp_pipeline filter_regions -n var.flt.vcf sample_dirs_file "$referenceFilePath" $RemoveAbnormalSnp_ExtraParams
_EOF_
)
else
	cfsan_snp_pipeline filter_regions -n var.flt.vcf sample_dirs_file "$referenceFilePath" $RemoveAbnormalSnp_ExtraParams 2>&1 | tee $logDir/filterAbnormalSNP.log
fi

#Starting from here, there are 2 threads:
#Thread X.1: the thread processing the original VCF files and corresponding downstream results
#Thread X.2: the thread processing the preserved VCF files and corresponding downstream results


echo -e "\nStep 6.1 - Combine the SNP positions across all samples into the SNP list file"
# The create_snp_list process creates the filtered list of sample directories.  It is the list of samples having removed the samples with excessive snps.
# When running on a workstation, the file exists at this point during the script execution, but on grid or torque, it has not yet been created. However,
# we know the path to the file regardless of whether it exists yet.
filteredSampleDirsFile="${sampleDirsFile}.OrigVCF.filtered"
touch $filteredSampleDirsFile

if platform == "grid":
    snpListJobId=$(echo | qsub  -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N snpList
#$ -cwd
#$ -j y
#$ -V
#$ -hold_jid $filterAbnSNPJobId
#$ -o $logDir/snpList.log
cfsan_snp_pipeline merge_sites $forceFlag -n var.flt.vcf -o "$workDir/snplist.txt" $CreateSnpList_ExtraParams sample_dirs_file "$filteredSampleDirsFile"
_EOF_
)
elif platform == "torque":
    snpListJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N snpList
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterok:$filterAbnSNPJobId
    #PBS -o $logDir/snpList.log
    #PBS -V
    cfsan_snp_pipeline merge_sites $forceFlag -n var.flt.vcf -o "$workDir/snplist.txt" $CreateSnpList_ExtraParams sample_dirs_file "$filteredSampleDirsFile"
_EOF_
)
else
    cfsan_snp_pipeline merge_sites $forceFlag -n var.flt.vcf -o "$workDir/snplist.txt" $CreateSnpList_ExtraParams sample_dirs_file "$filteredSampleDirsFile" 2>&1 | tee $logDir/snpList.log
fi

echo -e "\nStep 7.1 - Call the consensus SNPs for each sample"
if platform == "grid":
    callConsensusJobId=$(echo | qsub -terse -t 1-$sampleCount $GridEngine_QsubExtraParams << _EOF_
#$ -N callConsensus
#$ -cwd
#$ -V
#$ -j y
#$ -hold_jid $snpListJobId
#$ -o $logDir/callConsensus.log-\$TASK_ID
    sampleDir=\$(cat sample_dirs_file | head -n \$SGE_TASK_ID | tail -n 1)
    cfsan_snp_pipeline call_consensus $forceFlag -l "$workDir/snplist.txt" -o "\$sampleDir/consensus.fasta" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams  --vcfFileName consensus.vcf "\$sampleDir/reads.all.pileup"
_EOF_
)
elif platform == "torque":
    callConsensusJobId=$(echo | qsub -t 1-$sampleCount $Torque_QsubExtraParams << _EOF_
    #PBS -N callConsensus
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterok:$snpListJobId
    #PBS -o $logDir/callConsensus.log
    #PBS -V
    sampleDir=\$(cat sample_dirs_file | head -n \$PBS_ARRAYID | tail -n 1)
    cfsan_snp_pipeline call_consensus $forceFlag -l "$workDir/snplist.txt" -o "\$sampleDir/consensus.fasta" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams  --vcfFileName consensus.vcf "\$sampleDir/reads.all.pileup"
_EOF_
)
else
    if "$MaxConcurrentCallConsensus" != "":
        numCallConsensusCores=$MaxConcurrentCallConsensus
    else
        numCallConsensusCores=$numCores
    fi
    nl sample_dirs_file | xargs -n 2 -P $numCallConsensusCores bash -c 'set -o pipefail; cfsan_snp_pipeline call_consensus $forceFlag -l "$workDir/snplist.txt" -o "$1/consensus.fasta" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams  --vcfFileName consensus.vcf "$1/reads.all.pileup" 2>&1 | tee $logDir/callConsensus.log-$0'
fi

echo -e "\nStep 8.1 - Create the SNP matrix"
if platform == "grid":
    callConsensusJobArray=$(stripGridEngineJobArraySuffix $callConsensusJobId)
    snpMatrixJobId=$(echo | qsub -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N snpMatrix
#$ -cwd
#$ -V
#$ -j y
#$ -hold_jid $callConsensusJobArray
#$ -o $logDir/snpMatrix.log
    cfsan_snp_pipeline snp_matrix $forceFlag -c consensus.fasta -o "$workDir/snpma.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile"
_EOF_
)
elif platform == "torque":
    callConsensusJobArray=$(stripTorqueJobArraySuffix $callConsensusJobId)
    snpMatrixJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N snpMatrix
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterokarray:$callConsensusJobArray
    #PBS -o $logDir/snpMatrix.log
    #PBS -V
    cfsan_snp_pipeline snp_matrix $forceFlag -c consensus.fasta -o "$workDir/snpma.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile"
_EOF_
)
else
    cfsan_snp_pipeline snp_matrix $forceFlag -c consensus.fasta -o "$workDir/snpma.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile" 2>&1 | tee $logDir/snpMatrix.log
fi

echo -e "\nStep 9.1 - Create the reference base sequence"
if platform == "grid":
    snpReferenceJobId=$(echo | qsub -terse $GridEngine_QsubExtraParams << _EOF_
#$ -V
#$ -N snpReference
#$ -cwd
#$ -j y
#$ -hold_jid $callConsensusJobArray
#$ -o $logDir/snpReference.log
    cfsan_snp_pipeline snp_reference $forceFlag -l "$workDir/snplist.txt" -o "$workDir/referenceSNP.fasta" $CreateSnpReferenceSeq_ExtraParams "$referenceFilePath"
_EOF_
)
elif platform == "torque":
    snpReferenceJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N snpReference
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterokarray:$callConsensusJobArray
    #PBS -o $logDir/snpReference.log
    #PBS -V
    cfsan_snp_pipeline snp_reference $forceFlag -l "$workDir/snplist.txt" -o "$workDir/referenceSNP.fasta" $CreateSnpReferenceSeq_ExtraParams "$referenceFilePath"
_EOF_
)
else
    cfsan_snp_pipeline snp_reference $forceFlag -l "$workDir/snplist.txt" -o "$workDir/referenceSNP.fasta" $CreateSnpReferenceSeq_ExtraParams "$referenceFilePath" 2>&1 | tee $logDir/snpReference.log
fi


echo -e "\nStep 10.1 - Create the Multi-VCF file"
if $CallConsensus_ExtraParams =~ .*vcfFileName.*:
    if platform == "grid":
        mergeVcfJobId=$(echo | qsub  -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N mergeVcf
#$ -cwd
#$ -j y
#$ -V
#$ -hold_jid $callConsensusJobArray
#$ -o $logDir/mergeVcf.log
        cfsan_snp_pipeline merge_vcfs $forceFlag -o "$workDir/snpma.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile"
_EOF_
)
    elif platform == "torque":
        mergeVcfJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
        #PBS -N mergeVcf
        #PBS -d $(pwd)
        #PBS -j oe
        #PBS -W depend=afterokarray:$callConsensusJobArray
        #PBS -o $logDir/mergeVcf.log
        #PBS -V
        cfsan_snp_pipeline merge_vcfs $forceFlag -o "$workDir/snpma.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile"
_EOF_
)
    else
        cfsan_snp_pipeline merge_vcfs $forceFlag -o "$workDir/snpma.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile" 2>&1 | tee $logDir/mergeVcf.log
    fi
else
    echo -e "Skipped per CallConsensus_ExtraParams configuration"
fi

echo -e "\nStep 11.1 - Calculate SNP distance matrix"
if platform == "grid":
    calcSnpDistanceJobId=$(echo | qsub  -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N snpDistance
#$ -cwd
#$ -j y
#$ -V
#$ -hold_jid $snpMatrixJobId
#$ -o $logDir/calcSnpDistances.log
    cfsan_snp_pipeline distance $forceFlag -p "$workDir/snp_distance_pairwise.tsv" -m "$workDir/snp_distance_matrix.tsv" "$workDir/snpma.fasta"
_EOF_
)
elif platform == "torque":
    calcSnpDistanceJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N snpDistance
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterok:$snpMatrixJobId
    #PBS -o $logDir/calcSnpDistances.log
    #PBS -V
    cfsan_snp_pipeline distance $forceFlag -p "$workDir/snp_distance_pairwise.tsv" -m "$workDir/snp_distance_matrix.tsv" "$workDir/snpma.fasta"
_EOF_
)
else
    cfsan_snp_pipeline distance $forceFlag -p "$workDir/snp_distance_pairwise.tsv" -m "$workDir/snp_distance_matrix.tsv" "$workDir/snpma.fasta" 2>&1 | tee $logDir/calcSnpDistances.log
fi

# Step 14.1 - Notify user of any non-fatal errors accumulated during processing
if -s "$errorOutputFile" && $SnpPipeline_StopOnSampleError != true:
    echo "" 1>&2
    echo "There were errors processing some samples." 1>&2
    echo "See the log file $errorOutputFile for a summary of errors." 1>&2
fi

#Starting now are codes processing preserved SNPs after SNP filtering.
echo -e "\nStep 6.2 - Combine the SNP positions across all samples into the SNP list file"
###Create another copy of sample directories file, for the thread processing preserved snp files.
filteredSampleDirsFile2="${sampleDirsFile}.PresVCF.filtered"
touch $filteredSampleDirsFile2

#cp $sampleDirsFile $sampleDirsFile_Preserved
if platform == "grid":
    snpListJobId2=$(echo | qsub  -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N snpList_preserved
#$ -cwd
#$ -j y
#$ -V
#$ -hold_jid $filterAbnSNPJobId
#$ -o $logDir/snpList_preserved.log
    cfsan_snp_pipeline merge_sites $forceFlag -n var.flt_preserved.vcf -o "$workDir/snplist_preserved.txt" $CreateSnpList_ExtraParams sample_dirs_file "$filteredSampleDirsFile2"
_EOF_
)
elif platform == "torque":
    snpListJobId2=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N snpList_preserved
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterok:$filterAbnSNPJobId
    #PBS -o $logDir/snpList_preserved.log
    #PBS -V
    cfsan_snp_pipeline merge_sites $forceFlag -n var.flt_preserved.vcf -o "$workDir/snplist_preserved.txt" $CreateSnpList_ExtraParams sample_dirs_file "$filteredSampleDirsFile2"
_EOF_
)
else
    cfsan_snp_pipeline merge_sites $forceFlag -n var.flt_preserved.vcf -o "$workDir/snplist_preserved.txt" $CreateSnpList_ExtraParams sample_dirs_file "$filteredSampleDirsFile2" 2>&1 | tee $logDir/snpList_preserved.log
fi


echo -e "\nStep 7.2 - Call the consensus SNPs for each sample"
if platform == "grid":
    callConsensusJobId2=$(echo | qsub -terse -t 1-$sampleCount $GridEngine_QsubExtraParams << _EOF_
#$ -N callConsensus_preserved
#$ -cwd
#$ -V
#$ -j y
#$ -hold_jid $snpListJobId2
#$ -o $logDir/callConsensus_preserved.log-\$TASK_ID
    sampleDir=\$(cat sample_dirs_file | head -n \$SGE_TASK_ID | tail -n 1)
    cfsan_snp_pipeline call_consensus $forceFlag -l "$workDir/snplist_preserved.txt" -o "\$sampleDir/consensus_preserved.fasta" -e "\$sampleDir/var.flt_removed.vcf" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams  --vcfFileName consensus_preserved.vcf "\$sampleDir/reads.all.pileup"
_EOF_
)
elif platform == "torque":
    callConsensusJobId2=$(echo | qsub -t 1-$sampleCount $Torque_QsubExtraParams << _EOF_
    #PBS -N callConsensus_preserved
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterok:$snpListJobId2
    #PBS -o $logDir/callConsensus_preserved.log
    #PBS -V
    sampleDir=\$(cat sample_dirs_file | head -n \$PBS_ARRAYID | tail -n 1)
    cfsan_snp_pipeline call_consensus $forceFlag -l "$workDir/snplist_preserved.txt" -o "\$sampleDir/consensus_preserved.fasta" -e "\$sampleDir/var.flt_removed.vcf" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams  --vcfFileName consensus_preserved.vcf "\$sampleDir/reads.all.pileup"
_EOF_
)
else
    if "$MaxConcurrentCallConsensus" != "":
        numCallConsensusCores=$MaxConcurrentCallConsensus
    else
        numCallConsensusCores=$numCores
    fi
    nl sample_dirs_file | xargs -n 2 -P $numCallConsensusCores bash -c 'set -o pipefail; cfsan_snp_pipeline call_consensus $forceFlag -l "$workDir/snplist_preserved.txt" -o "$1/consensus_preserved.fasta" -e "$1/var.flt_removed.vcf" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams  --vcfFileName consensus_preserved.vcf "$1/reads.all.pileup" 2>&1 | tee $logDir/callConsensus_preserved.log-$0'
fi

echo -e "\nStep 8.2 - Create the SNP matrix"
if platform == "grid":
    callConsensusJobArray2=$(stripGridEngineJobArraySuffix $callConsensusJobId2)
    snpMatrixJobId2=$(echo | qsub -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N snpMatrix_preserved
#$ -cwd
#$ -V
#$ -j y
#$ -hold_jid $callConsensusJobArray2
#$ -o $logDir/snpMatrix_preserved.log
    cfsan_snp_pipeline snp_matrix $forceFlag -c consensus_preserved.fasta -o "$workDir/snpma_preserved.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile2"
_EOF_
)
elif platform == "torque":
    callConsensusJobArray2=$(stripTorqueJobArraySuffix $callConsensusJobId2)
    snpMatrixJobId2=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N snpMatrix_preserved
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterokarray:$callConsensusJobArray2
    #PBS -o $logDir/snpMatrix_preserved.log
    #PBS -V
    cfsan_snp_pipeline snp_matrix $forceFlag -c consensus_preserved.fasta -o "$workDir/snpma_preserved.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile2"
_EOF_
)
else
    cfsan_snp_pipeline snp_matrix $forceFlag -c consensus_preserved.fasta -o "$workDir/snpma_preserved.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile2" 2>&1 | tee $logDir/snpMatrix_preserved.log
fi

echo -e "\nStep 9.2 - Create the reference base sequence"
if platform == "grid":
    snpReferenceJobId2=$(echo | qsub -terse $GridEngine_QsubExtraParams << _EOF_
#$ -V
#$ -N snpReference_preserved
#$ -cwd
#$ -j y
#$ -hold_jid $callConsensusJobArray2
#$ -o $logDir/snpReference_preserved.log
    cfsan_snp_pipeline snp_reference $forceFlag -l "$workDir/snplist_preserved.txt" -o "$workDir/referenceSNP_preserved.fasta" $CreateSnpReferenceSeq_ExtraParams "$referenceFilePath"
_EOF_
)
elif platform == "torque":
    snpReferenceJobId2=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N snpReference_preserved
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterokarray:$callConsensusJobArray2
    #PBS -o $logDir/snpReference_preserved.log
    #PBS -V
    cfsan_snp_pipeline snp_reference $forceFlag -l "$workDir/snplist_preserved.txt" -o "$workDir/referenceSNP_preserved.fasta" $CreateSnpReferenceSeq_ExtraParams "$referenceFilePath"
_EOF_
)
else
    cfsan_snp_pipeline snp_reference $forceFlag -l "$workDir/snplist_preserved.txt" -o "$workDir/referenceSNP_preserved.fasta" $CreateSnpReferenceSeq_ExtraParams "$referenceFilePath" 2>&1 | tee $logDir/snpReference_preserved.log
fi


echo -e "\nStep 10.2 - Create the Multi-VCF file"
if $CallConsensus_ExtraParams =~ .*vcfFileName.*:
    if platform == "grid":
        mergeVcfJobId2=$(echo | qsub  -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N mergeVcf_preserved
#$ -cwd
#$ -j y
#$ -V
#$ -hold_jid $callConsensusJobArray2
#$ -o $logDir/mergeVcf_preserved.log
        cfsan_snp_pipeline merge_vcfs $forceFlag -n consensus_preserved.vcf -o "$workDir/snpma_preserved.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile2"
_EOF_
)
elif platform == "torque":
        mergeVcfJobId2=$(echo | qsub $Torque_QsubExtraParams << _EOF_
        #PBS -N mergeVcf_preserved
        #PBS -d $(pwd)
        #PBS -j oe
        #PBS -W depend=afterokarray:$callConsensusJobArray2
        #PBS -o $logDir/mergeVcf_preserved.log
        #PBS -V
        cfsan_snp_pipeline merge_vcfs $forceFlag -n consensus_preserved.vcf -o "$workDir/snpma_preserved.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile2"
_EOF_
)
    else
        cfsan_snp_pipeline merge_vcfs $forceFlag -n consensus_preserved.vcf -o "$workDir/snpma_preserved.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile2" 2>&1 | tee $logDir/mergeVcf_preserved.log
    fi
else
    echo -e "Skipped per CallConsensus_ExtraParams configuration"
fi

echo -e "\nStep 11.2 - Calculate SNP distance matrix"
if platform == "grid":
    calcSnpDistanceJobId2=$(echo | qsub  -terse $GridEngine_QsubExtraParams << _EOF_
#$ -N snpDistance_preserved
#$ -cwd
#$ -j y
#$ -V
#$ -hold_jid $snpMatrixJobId2
#$ -o $logDir/calcSnpDistances_preserved.log
    cfsan_snp_pipeline distance $forceFlag -p "$workDir/snp_distance_pairwise_preserved.tsv" -m "$workDir/snp_distance_matrix_preserved.tsv" "$workDir/snpma_preserved.fasta"
_EOF_
)
elif platform == "torque":
    calcSnpDistanceJobId2=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N snpDistance_preserved
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterok:$snpMatrixJobId2
    #PBS -o $logDir/calcSnpDistances_preserved.log
    #PBS -V
    cfsan_snp_pipeline distance $forceFlag -p "$workDir/snp_distance_pairwise_preserved.tsv" -m "$workDir/snp_distance_matrix_preserved.tsv" "$workDir/snpma_preserved.fasta"
_EOF_
)
else
    cfsan_snp_pipeline distance $forceFlag -p "$workDir/snp_distance_pairwise_preserved.tsv" -m "$workDir/snp_distance_matrix_preserved.tsv" "$workDir/snpma_preserved.fasta" 2>&1 | tee $logDir/calcSnpDistances_preserved.log
fi

echo -e "\nStep 12 - Collect metrics for each sample"
if platform == "grid":
    collectSampleMetricsJobId=$(echo | qsub -terse -t 1-$sampleCount $GridEngine_QsubExtraParams << _EOF_
#$ -N collectMetrics
#$ -cwd
#$ -V
#$ -j y
#$ -hold_jid_ad $callConsensusJobArray,$callConsensusJobArray2
#$ -o $logDir/collectSampleMetrics.log-\$TASK_ID
    sampleDir=\$(cat sample_dirs_file | head -n \$SGE_TASK_ID | tail -n 1)
    cfsan_snp_pipeline collect_metrics -o "\$sampleDir/metrics" $CollectSampleMetrics_ExtraParams "\$sampleDir"  "$referenceFilePath"
_EOF_
)
elif platform == "torque":
    collectSampleMetricsJobId=$(echo | qsub -t 1-$sampleCount $Torque_QsubExtraParams << _EOF_
    #PBS -N collectMetrics
    #PBS -d $(pwd)
    #PBS -j oe
    #PBS -W depend=afterokarray:$callConsensusJobArray:$callConsensusJobArray2
    #PBS -o $logDir/collectSampleMetrics.log
    #PBS -V
    sampleDir=\$(cat sample_dirs_file | head -n \$PBS_ARRAYID | tail -n 1)
    cfsan_snp_pipeline collect_metrics -o "\$sampleDir/metrics" $CollectSampleMetrics_ExtraParams "\$sampleDir"  "$referenceFilePath"
_EOF_
)
else
    if "$MaxConcurrentCollectSampleMetrics" != "":
        numCollectSampleMetricsCores=$MaxConcurrentCollectSampleMetrics
    else
        numCollectSampleMetricsCores=$numCores
    fi
    nl sample_dirs_file | xargs -n 2 -P $numCollectSampleMetricsCores bash -c 'set -o pipefail; cfsan_snp_pipeline collect_metrics -o "$1/metrics" $CollectSampleMetrics_ExtraParams "$1" "$referenceFilePath" 2>&1 | tee $logDir/collectSampleMetrics.log-$0'
fi


print("\nStep 13 - Combine the metrics across all samples into the metrics table")
command_line = 'cfsan_snp_pipeline combine_metrics -n metrics -o "$workDir/metrics.tsv" $CombineSampleMetrics_ExtraParams sample_dirs_file'
log_file = '$logDir/combineSampleMetrics.log'
combine_metrics_job_id = runner.run(command_line, "combineMetrics", log_file, wait_for_array=[collectSampleMetricsJobArray])

# Step 14.2 - Notify user of any non-fatal errors accumulated during processing
if os.path.getsize(error_output_file) > 0 and not stop_on_error:
    print("\nThere were errors processing some samples.\nSee the log file %s for a summary of errors." % error_output_file, file=sys.stderr)

"""
