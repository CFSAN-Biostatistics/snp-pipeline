#! /bin/bash

# Purpose: Common utility functions used by the snppipline code.
# History:
#   20150417-scd: Started.
# Notes:
#
# Bugs:
#
# References:
# http://stackoverflow.com/questions/3966048/access-arguments-to-bash-script-inside-a-function
# http://stackoverflow.com/questions/3915040/bash-fish-command-to-print-absolute-path-to-a-file/23002317#23002317
#

# Print the arguments of the calling script, in order.
# Call like this:
#       scriptArgs=$(getScriptArgs)
# Coded by Dean Hall on StackOverflow
# http://stackoverflow.com/questions/3966048/access-arguments-to-bash-script-inside-a-function
getScriptArgs()
{
    # Get the number of arguments passed to this script
    # (The BASH_ARGV array does not include $0.)
    local n=${#BASH_ARGV[@]}

    if (( $n > 0 ))
    then
        # Get the last index of the args in BASH_ARGV
        local n_index=$(( $n - 1 ))

        # Loop through the indexes from largest to smallest
        for i in $(seq ${n_index} -1 0)
        do
            # Print a space if necessary
            if (( $i < $n_index ))
            then
                echo -n ' '
            fi

            # Print the actual argument
            echo -n "${BASH_ARGV[$i]}"
        done

        # Print a newline
        echo
    fi
}


# Convert a relative path to absolute path.
# This works for both directories and files.
# Call like this:
#       absolutePath=$(getAbsPath "$relativePath")
# Coded by Alexander Klimetschek on StackOverflow
# http://stackoverflow.com/questions/3915040/bash-fish-command-to-print-absolute-path-to-a-file/23002317#23002317
getAbsPath()
{
    # $1     : relative path
    # return : absolute path
    if [ -d "$1" ]; then
        # dir
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        # file
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    fi
}

# Get the last token in a string.
# Can be used to get the last argument on a command line.
# Call like this:
#       token=$(getLastToken "$*")
getLastToken()
{
    local i
    for i in $1; do :; done
    echo "$i"
}


# Given a sample directory or a sample fastq file, return the sample id
# which is just the basename of the sample directory.
# Call like this:
#       sampleId=$(getSampleId "$samplePath")
getSampleId()
{
    # $1     : sample directory or file
    # return : basename of the sample directory
    local sampleId
    sampleId=$(getAbsPath "$1")
    if [ -d "$sampleId" ]; then
        # dir
        sampleId="${sampleId%/}"        # Strip trailing /, if any
        sampleId="${sampleId##*/}"      # Strip all leading directories
        echo "$sampleId"
    elif [ -f "$sampleId" ]; then
        # file
        sampleId="${sampleId%/*}"       # Strip file name
        sampleId="${sampleId##*/}"      # Strip all leading directories
        echo "$sampleId"
    fi
}

# Read a single parameter from a file containing parameters of the form name=value.
# The parameter becomes a READ-ONLY global bash variable.
readParameter()
{
  parameterFile=$1
  parameterName=$2
  while IFS='' read -r line || [[ -n "$line" ]]
  do
    if [[ "$line" =~ "=" ]]; then
      param=${line%%=*}  # strip everything after =
      if [[ $param == $parameterName ]]; then
        value=${line##*=}  # strip everything before =
        readonly "$param"="$value"
      fi
    fi
  done  < "$parameterFile" # This syntax without piping is needed to retain the values of variables declared in the loop
}


# Write an error message to the error log if enabled.
logError()
{
    # Log to error file if it is defined, which happens automatically if running run_snp_pipeline.sh.
    # The errorOutputFile may also be defined manually if running scripts without run_snp_pipeline.sh.
    if [ -n "$errorOutputFile" ]; then
        echo "$@" >> "$errorOutputFile"
    fi
}


# Log an error message to the error summary file and to stderr.
reportError()
{
    errorMsg="$1"

    # Send the error to the error log
    logError "$errorMsg"

    # Also send the error message to stderr
    echo "$errorMsg" 1>&2
}


# Log a fatal error to the error summary file and exit with error code 100 to cause Sun Grid Engine to
# also detect the error.  Global errors prevent subsequent processing for all samples.
globalError()
{
    errorMsg="$1"

    logError "$(basename $0) failed."
    logError "$errorMsg"
    logError "================================================================================"

    # Also send the detail error message to stderr -- this will put the error message in the 
    # process-specific log file.
    echo "$errorMsg" 1>&2

    # Exit 100 does two things:
    # 1. Sun Grid Engine will stop execution of dependent jobs
    # 2. run_snp_pipeline.sh will know this error has already been reported
    exit 100
}


# Log an error to the error summary file and exit with error code 100 to cause Sun Grid Engine to
# also detect the error.  Sample errors only affect one sample and can be conditionally ignored.
sampleError()
{
    errorMsg="$1"
    continuePossible="$2"

    if [[ -z $SnpPipeline_StopOnSampleError || $SnpPipeline_StopOnSampleError = true || "$continuePossible" != true ]]; then
        logError "$(basename $0) failed."
    else
        logError "$(basename $0)"
    fi
    logError "$errorMsg"
    logError "================================================================================"

    # Also send the detail error message to stderr -- this will put the error message in the 
    # process-specific log file.
    echo "$errorMsg" 1>&2

    # Exit 100 tells Sun Grid Engine to stop execution of dependent jobs
    if [[ -z $SnpPipeline_StopOnSampleError || $SnpPipeline_StopOnSampleError = true ]]; then
        exit 100 # run_snp_pipeline.sh will know this error has already been reported
    elif [[ "$continuePossible" != true ]]; then
        exit 98 # run_snp_pipeline.sh will know this error has already been reported, but it should not stop execution of the rest of the pipeline
    fi
}

# Log a warning to the error summary file but do not exit regardless of the SnpPipeline_StopOnSampleError setting
sampleWarning()
{
    errorMsg="$1"

    logError "$(basename $0) warning:"
    logError "$errorMsg"
    logError "================================================================================"

    # Also send the detail error message to stderr -- this will put the error message in the 
    # process-specific log file.
    echo "$errorMsg" 1>&2
}


# Generate a global error if a specified file is missing or empty.
# $1 = file
# $2 = process that should have created the file
globalErrorOnMissingFile()
{
    local file
    local process
    file=$1
    process=$2
    if [[ ! -e "$file" ]]; then globalError "Error: $file does not exist after running $process." false; fi
    if [[ ! -s "$file" ]]; then globalError "Error: $file is empty after running $process." false; fi
}


# Generate a sample error if a specified file is missing or empty.
# $1 = file
# $2 = process that should have created the file
sampleErrorOnMissingFile()
{
    local file
    local process
    file=$1
    process=$2
    if [[ ! -e "$file" ]]; then sampleError "Error: $file does not exist after running $process." false; fi
    if [[ ! -s "$file" ]]; then sampleError "Error: $file is empty after running $process." false; fi
}


# Generate a sample error if a specified file contains a specified string
# $1 = file
# $2 = string
# $3 = process that created the file
sampleErrorOnFileContains()
{
    local file
    local target
    local process
    file=$1
    target=$2
    process=$3
    grep -i "$target" "$file" &> /dev/null
    if [[ $? == 0 ]]; then sampleError "Error: $file contains unexpected text: '$target' after running $process." false; fi
}


# This function is automatically called by bash when error traps are setup.
# It Logs the error and returns error code 100 to cause Sun Grid Engine to
# also detect the error.
handleTrappedGlobalErrors()
{
    errorCode=$?
    bashCommand="$BASH_COMMAND"
    scriptName=$(basename $0)
    scriptArgs=$(getScriptArgs)
    logError "Error detected while running $scriptName."
    logError ""
    logError "The command line was:"
    logError "    $scriptName $scriptArgs"
    logError ""
    logError "The command at line ${BASH_LINENO[0]} returned error code $errorCode:"
    logError "    $bashCommand"
    logError "================================================================================"
    # $BASH_COMMAND contains the command that was being executed at the time of the trap
    # ${BASH_LINENO[0]} contains the line number in the script of that command

    # Sun Grid Engine looks for error code 100
    exit 100
}


# This function is automatically called by bash when error traps are setup.
# It Logs the error and returns error code 100 to cause Sun Grid Engine to
# also detect the error.
handleTrappedSampleErrors()
{
    errorCode=$?
    bashCommand="$BASH_COMMAND"
    scriptName=$(basename $0)
    scriptArgs=$(getScriptArgs)
    logError "Error detected while running $scriptName."
    logError ""
    logError "The command line was:"
    logError "    $scriptName $scriptArgs"
    logError ""
    logError "The command at line ${BASH_LINENO[0]} returned error code $errorCode:"
    logError "    $bashCommand"
    logError "================================================================================"
    # $BASH_COMMAND contains the command that was being executed at the time of the trap
    # ${BASH_LINENO[0]} contains the line number in the script of that command

    # Sun Grid Engine looks for error code 100
    if [[ -z $SnpPipeline_StopOnSampleError || $SnpPipeline_StopOnSampleError = true ]]; then
        exit 100 # run_snp_pipeline.sh will know this error has already been reported
    else
        exit 98 # run_snp_pipeline.sh will know this error has already been reported, but it should not stop execution of the rest of the pipeline
    fi
}


# Register the error handling function for any errors detected during script execition.
setupGlobalErrorHandler()
{
    trap handleTrappedGlobalErrors ERR
}


# Register the error handling function for any errors detected during script execition.
setupSampleErrorHandler()
{
    trap handleTrappedSampleErrors ERR
}


# Log an error to the error log and stderr if a specified program is not on the path.
# Print a true/false status to stdout.
# Implementation note: returning 1 will not cause a trap when error traps are enabled, that is why it prints a status.
# Call like this:
#       onPath=$(verifyOnPath "grep")
verifyOnPath()
{
    result=$(which "$1" || true)
    if [[ -z $result ]]; then
        reportError "$1 is not on the path"
        echo false
    else
        echo true
    fi
}


