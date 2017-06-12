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
import time
import traceback

from snppipeline import command
from snppipeline import fastq
from snppipeline.job_runner import JobRunner
from snppipeline.job_runner import JobRunnerException
from snppipeline import utils
from snppipeline.utils import log_error


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

    elif True: # Send the exception traceback to stderr and to the error.log file

        # Format a stack trace and the exception information into a list of new-line terminated strings
        trace_lines = traceback.format_exception(exc_type, exc_value, exc_traceback, limit=None)

        # Format the exception report
        exc_string = ''.join(trace_lines)
        utils.report_error(exc_string)

    else:  # TODO: this block of code can probably be deleted
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


def run(args):
    """Run all the steps of the snp pipeline in th correct order.

    Parameters
    ----------
    args : Namespace
        referenceFile : str
            Relative or absolute path to the reference fasta file
        forceFlag : bool
            Force processing even when result files already exist and are newer than inputs
        mirror : str
            Mode to create a mirror copy of the reference directory and all the sample directories.
            Possible values: {soft, hard, copy}
        configFile : str
            Relative or absolute path to a configuration file for overriding defaults and defining
            extra parameters for the tools and scripts within the pipeline.
        jobQueueMgr : str
            Job queue manager for remote parallel job execution in an HPC environment.  Currently
            "torque" and "grid" are supported.  If not specified, the pipeline will execute locally.
        workDir : str
            Output directory for the result files.
        samplesDir : str
            Relative or absolute path to the parent directory of all the sample directories.
        samplesFile : str
            Relative or absolute path to a file listing all of the sample directories.
    """

    #stop_on_error_env = os.environ.get("SnpPipeline_StopOnSampleError")
    #stop_on_error = stop_on_error_env is None or stop_on_error_env == "true"

    # Erase any left-over error log environment variable from a previous run
    os.environ.pop("errorOutputFile", None) # the 2nd arg avoids an exception when not in dict

    # Handle output working directory.  Create the directory if it does not exist.
    # Any errors creating the work_dir will not be logged to the error log because
    # the error log belongs in the work_dir.
    work_dir = args.workDir
    try:
        utils.mkdir_p(work_dir)
    except OSError as exc:
        utils.fatal_error("Error: could not create the output directory %s" % work_dir)
    if not utils.is_directory_writeable(work_dir):
        utils.fatal_error("Error: output directory % is not writable." % work_dir)

    # The error log is in the main workdir
    error_output_file = os.path.join(work_dir, "error.log")
    os.environ["errorOutputFile"] = error_output_file
    # TODO: copy old error log to old logs directory, because otherwise it will be removed and lost forever
    if os.path.isfile(error_output_file):
        os.remove(error_output_file)

    # Validate reference fasta file
    reference_file_path = args.referenceFile
    if not os.path.isfile(reference_file_path):
        utils.fatal_error("Error: reference file %s does not exist." % reference_file_path)
    if os.path.getsize(reference_file_path) == 0:
        utils.fatal_error("Error: reference file %s is empty." % reference_file_path)
    reference_file_name = os.path.basename(reference_file_path)

    # Force rebuild flag is passed to all the subtask commands below
    force_flag = " -f " if args.forceFlag else ""

    # Create the logs directory with name like "logs-20170215.144253"
    run_time_stamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
    log_dir = os.path.join(work_dir, "logs-"+run_time_stamp)
    try:
        utils.mkdir_p(log_dir)
    except OSError as exc:
        utils.fatal_error("Error: could not create the logs directory %s" % log_dir)
    if not utils.is_directory_writeable(work_dir):
        utils.fatal_error("Error: logs directory % is not writable." % log_dir)

    # Handle configuration file, use the specified file, or create a default file
    if args.configFile:
        config_file_path = args.configFile
        if not os.path.isfile(reference_file_path):
            utils.fatal_error("Error: configuration file %s does not exist." % config_file_path)
        if os.path.getsize(reference_file_path) == 0:
            utils.fatal_error("Error: configuration file %s is empty." % config_file_path)

        shutil.copy2(config_file_path, log_dir)  # copy2 tries to preserve timestamps
        config_params = utils.read_properties(config_file_path, recognize_vars=True)
    else:
        command.run("cfsan_snp_pipeline data configurationFile " + log_dir, outfile=sys.stdout)
        config_file_path = os.path.join(log_dir, "snppipeline.conf")
        config_params = utils.read_properties(config_file_path, recognize_vars=True)

    # Validate the configured aligner choice
    snp_pipeline_aligner = config_params.get("SnpPipeline_Aligner", "").lower() or "bowtie2"
    if snp_pipeline_aligner not in ["bowtie2", "smalt"]:
        utils.fatal_error("Config file error in SnpPipeline_Aligner parameter: only bowtie2 and smalt aligners are supported.")
    os.environ["SnpPipeline_Aligner"] = snp_pipeline_aligner

    # Stop the pipeline by default upon single sample errors if not configured either way
    stop_on_error = config_params.get("SnpPipeline_StopOnSampleError", "").lower() or "true"
    os.environ["SnpPipeline_StopOnSampleError"] = stop_on_error

    # Put the configuration parameters into the process environment variables
    os.environ["Bowtie2Build_ExtraParams"] = config_params.get("Bowtie2Build_ExtraParams", "")
    os.environ["SmaltIndex_ExtraParams"] = config_params.get("SmaltIndex_ExtraParams", "")
    os.environ["SamtoolsFaidx_ExtraParams"] = config_params.get("SamtoolsFaidx_ExtraParams", "")
    os.environ["Bowtie2Align_ExtraParams"] = config_params.get("Bowtie2Align_ExtraParams", "")
    os.environ["SmaltAlign_ExtraParams"] = config_params.get("SmaltAlign_ExtraParams", "")
    os.environ["SamtoolsSamFilter_ExtraParams"] = config_params.get("SamtoolsSamFilter_ExtraParams", "")
    os.environ["SamtoolsSort_ExtraParams"] = config_params.get("SamtoolsSort_ExtraParams", "")
    os.environ["SnpPipeline_RemoveDuplicateReads"] = config_params.get("SnpPipeline_RemoveDuplicateReads", "").lower() or "true"
    os.environ["PicardMarkDuplicates_ExtraParams"] = config_params.get("PicardMarkDuplicates_ExtraParams", "")
    os.environ["PicardJvm_ExtraParams"] = config_params.get("PicardJvm_ExtraParams", "")
    os.environ["SamtoolsMpileup_ExtraParams"] = config_params.get("SamtoolsMpileup_ExtraParams", "")
    os.environ["VarscanMpileup2snp_ExtraParams"] = config_params.get("VarscanMpileup2snp_ExtraParams", "")
    os.environ["VarscanJvm_ExtraParams"] = config_params.get("VarscanJvm_ExtraParams", "")
    os.environ["RemoveAbnormalSnp_ExtraParams"] = config_params.get("RemoveAbnormalSnp_ExtraParams", "")
    os.environ["CreateSnpList_ExtraParams"] = config_params.get("CreateSnpList_ExtraParams", "")
    os.environ["CallConsensus_ExtraParams"] = config_params.get("CallConsensus_ExtraParams", "")
    os.environ["CreateSnpMatrix_ExtraParams"] = config_params.get("CreateSnpMatrix_ExtraParams", "")
    os.environ["BcftoolsMerge_ExtraParams"] = config_params.get("BcftoolsMerge_ExtraParams", "")
    os.environ["CreateSnpReferenceSeq_ExtraParams"] = config_params.get("CreateSnpReferenceSeq_ExtraParams", "")
    os.environ["CollectSampleMetrics_ExtraParams"] = config_params.get("CollectSampleMetrics_ExtraParams", "")
    os.environ["CombineSampleMetrics_ExtraParams"] = config_params.get("CombineSampleMetrics_ExtraParams", "")
    os.environ["GridEngine_PEname"] = config_params.get("GridEngine_PEname", "")

    # Verify the dependencies are available on the path
    dependencies = ["cfsan_snp_pipeline", snp_pipeline_aligner, "samtools", "java", "tabix", "bgzip", "bcftools"]
    found_all_dependencies = True
    for executable in dependencies:
        if not utils.which(executable):
            utils.report_error(executable + " is not on the path")
            found_all_dependencies = False

    stdout = command.run("java net.sf.varscan.VarScan 2>&1")
    if "Error" in stdout:
        utils.report_error("CLASSPATH is not configured with the path to VarScan")
        found_all_dependencies = False

    if os.environ["SnpPipeline_RemoveDuplicateReads"] == "true":
        stdout = command.run("java picard.cmdline.PicardCommandLine 2>&1")
        if "Error" in stdout:
            utils.report_error("CLASSPATH is not configured with the path to Picard")
            found_all_dependencies = False

    if not found_all_dependencies:
        utils.fatal_error("Check the SNP Pipeline installation instructions here: http://snp-pipeline.readthedocs.org/en/latest/installation.html")


    def rewrite_cleansed_file_of_sample_dirs(inpath, outpath):
        """Rewrite the file of sample directories, removing trailing slashes and blank lines.

        Parameters
        ----------
        inpath : str
            Path to the input file of sample directories
        outpath : str
            Path to the input file of sample directories
        """
        with open(inpath) as f:
            lines = f.read().splitlines()
        lines = [line.rstrip('/') for line in lines] # strip trailing slash
        lines = [line.strip() for line in lines] # strip leading and trailing spaces
        lines = [line for line in lines if line] # discard blank lines

        with open(outpath, 'w') as f:
            for line in lines:
                print(line, file=f)


    def validate_file_of_sample_dirs(sample_dirs_file):
        """Verify each of the sample directories in the sample directory file is not empty and contains fastq's.

        Parameters
        ----------
        sample_dirs_file : str
            Path to the file of sample directories
        """
        found_error = False

        with open(sample_dirs_file) as f:
            for directory in f:
                directory = directory.strip()
                if not utils.verify_non_empty_directory("Sample directory", directory):
                    found_error = True
                else:
                    files = fastq.list_fastq_files(directory)
                    if len(files) == 0:
                        utils.report_error("Sample directory %s does not contain any fastq files." % directory)
                        found_error = True

        if found_error:
            if os.environ.get("SnpPipeline_StopOnSampleError") == "true":
                sys.exit(1)
            else:
                log_error("================================================================================")


    # --------------------------------------------------------
    # get sample directories sorted by size, largest first

    # persistSortedSampleDirs()
    # {
    #     local samplesDir=$1
    #     local workDir=$2

    #     tmpFile=$(mktemp -p "$workDir" tmp.sampleDirs.XXXXXXXX)
    #     du -b -L "$samplesDir"/*/*.fastq* 2> /dev/null > "$tmpFile" || true
    #     du -b -L "$samplesDir"/*/*.fq* 2> /dev/null >> "$tmpFile" || true
    #     < "$tmpFile" xargs -n 2 sh -c 'echo $0 $(dirname "$1")' | \
    #     awk '{sizes[$2]+=$1} END {for (dir in sizes) {printf "%s %.0f\n", dir, sizes[dir]}}' | \
    #     sort -k 2 -n -r | \
    #     cut -f 1 -d" " > "$workDir/sampleDirectories.txt"
    #     rm "$tmpFile"
    # }


    # TODO: detect broken fastq symlinks
    if args.samplesDir:
        pass
    #     samplesDir="${opt_s_arg%/}"  # strip trailing slash
    #     if ! -d "$samplesDir": fatalError "Samples directory $samplesDir does not exist."; fi
    #     if   -z $(ls -A "$samplesDir"): fatalError "Samples directory $samplesDir is empty."; fi
    #     fastqFiles=$({ find "$samplesDir" -path "$samplesDir"'/*/*.fastq*'; find "$samplesDir" -path "$samplesDir"'/*/*.fq*'; })
    #     if   -z "$fastqFiles": fatalError "Samples directory $samplesDir does not contain subdirectories with fastq files."; fi
    #     persistSortedSampleDirs "$samplesDir" "$workDir"
    #     sampleDirsFile="$workDir/sampleDirectories.txt"

    if args.samplesFile:
        sample_dirs_file = args.samplesFile
        if not os.path.isfile(sample_dirs_file):
            utils.fatal_error("Error: the file of samples directories, %s does not exist." % sample_dirs_file)
        if os.path.getsize(sample_dirs_file) == 0:
            utils.fatal_error("Error: the file of samples directories, %s is empty." % sample_dirs_file)
        rewrite_cleansed_file_of_sample_dirs(sample_dirs_file, os.path.join(work_dir, "sampleDirectories.txt"))
        sample_dirs_file = os.path.join(work_dir, "sampleDirectories.txt")
        validate_file_of_sample_dirs(sample_dirs_file)

"""

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
    absoluteReferenceFilePath=$(get_abs_filename reference_file_path)
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
    cfsan_snp_pipeline index_ref + force_flag + reference_file_path
_EOF_
)
elif platform == "torque":
    prepReferenceJobId=$(echo | qsub $Torque_QsubExtraParams << _EOF_
    #PBS -N prepReference
    #PBS -j oe
    #PBS -d $(pwd)
    #PBS -o $logDir/prepReference.log
    #PBS -V
    cfsan_snp_pipeline index_ref + force_flag + reference_file_path
_EOF_
)
else
    cfsan_snp_pipeline index_ref + force_flag + reference_file_path 2>&1 | tee $logDir/prepReference.log
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
    cfsan_snp_pipeline map_reads + force_flag + reference_file_path \$(cat "$workDir/sampleFullPathNames.txt" | head -n \$SGE_TASK_ID | tail -n 1)
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
    cfsan_snp_pipeline map_reads + force_flag + reference_file_path \$samplesToAlign
_EOF_
)
else
    nl "$workDir/sampleFullPathNames.txt" | xargs -n 3 -L 1 bash -c 'set -o pipefail; cfsan_snp_pipeline map_reads + force_flag + reference_file_path $1 $2 2>&1 | tee $logDir/alignSamples.log-$0'
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
    cfsan_snp_pipeline call_sites + force_flag + reference_file_path "\$(cat sample_dirs_file | head -n \$SGE_TASK_ID | tail -n 1)"
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
    cfsan_snp_pipeline call_sites + force_flag + reference_file_path "\$sampleDir"
_EOF_
)
else
    if "$MaxConcurrentPrepSamples" != "":
        numPrepSamplesCores=$MaxConcurrentPrepSamples
    else
        numPrepSamplesCores=$numCores
    fi
    nl sample_dirs_file | xargs -n 2 -P $numPrepSamplesCores bash -c 'set -o pipefail; cfsan_snp_pipeline call_sites + force_flag + reference_file_path $1 2>&1 | tee $logDir/prepSamples.log-$0'
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
cfsan_snp_pipeline filter_regions -n var.flt.vcf sample_dirs_file reference_file_path $RemoveAbnormalSnp_ExtraParams
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
	cfsan_snp_pipeline filter_regions -n var.flt.vcf sample_dirs_file reference_file_path $RemoveAbnormalSnp_ExtraParams
_EOF_
)
else
	cfsan_snp_pipeline filter_regions -n var.flt.vcf sample_dirs_file reference_file_path $RemoveAbnormalSnp_ExtraParams 2>&1 | tee $logDir/filterAbnormalSNP.log
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
cfsan_snp_pipeline merge_sites + force_flag + -n var.flt.vcf -o "$workDir/snplist.txt" $CreateSnpList_ExtraParams sample_dirs_file "$filteredSampleDirsFile"
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
    cfsan_snp_pipeline merge_sites + force_flag + -n var.flt.vcf -o "$workDir/snplist.txt" $CreateSnpList_ExtraParams sample_dirs_file "$filteredSampleDirsFile"
_EOF_
)
else
    cfsan_snp_pipeline merge_sites + force_flag + -n var.flt.vcf -o "$workDir/snplist.txt" $CreateSnpList_ExtraParams sample_dirs_file "$filteredSampleDirsFile" 2>&1 | tee $logDir/snpList.log
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
    cfsan_snp_pipeline call_consensus + force_flag + -l "$workDir/snplist.txt" -o "\$sampleDir/consensus.fasta" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams  --vcfFileName consensus.vcf "\$sampleDir/reads.all.pileup"
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
    cfsan_snp_pipeline call_consensus + force_flag + -l "$workDir/snplist.txt" -o "\$sampleDir/consensus.fasta" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams  --vcfFileName consensus.vcf "\$sampleDir/reads.all.pileup"
_EOF_
)
else
    if "$MaxConcurrentCallConsensus" != "":
        numCallConsensusCores=$MaxConcurrentCallConsensus
    else
        numCallConsensusCores=$numCores
    fi
    nl sample_dirs_file | xargs -n 2 -P $numCallConsensusCores bash -c 'set -o pipefail; cfsan_snp_pipeline call_consensus + force_flag + -l "$workDir/snplist.txt" -o "$1/consensus.fasta" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams  --vcfFileName consensus.vcf "$1/reads.all.pileup" 2>&1 | tee $logDir/callConsensus.log-$0'
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
    cfsan_snp_pipeline snp_matrix + force_flag + -c consensus.fasta -o "$workDir/snpma.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile"
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
    cfsan_snp_pipeline snp_matrix + force_flag + -c consensus.fasta -o "$workDir/snpma.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile"
_EOF_
)
else
    cfsan_snp_pipeline snp_matrix + force_flag + -c consensus.fasta -o "$workDir/snpma.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile" 2>&1 | tee $logDir/snpMatrix.log
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
    cfsan_snp_pipeline snp_reference + force_flag + -l "$workDir/snplist.txt" -o "$workDir/referenceSNP.fasta" $CreateSnpReferenceSeq_ExtraParams reference_file_path
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
    cfsan_snp_pipeline snp_reference + force_flag + -l "$workDir/snplist.txt" -o "$workDir/referenceSNP.fasta" $CreateSnpReferenceSeq_ExtraParams reference_file_path
_EOF_
)
else
    cfsan_snp_pipeline snp_reference + force_flag + -l "$workDir/snplist.txt" -o "$workDir/referenceSNP.fasta" $CreateSnpReferenceSeq_ExtraParams reference_file_path 2>&1 | tee $logDir/snpReference.log
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
        cfsan_snp_pipeline merge_vcfs + force_flag + -o "$workDir/snpma.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile"
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
        cfsan_snp_pipeline merge_vcfs + force_flag + -o "$workDir/snpma.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile"
_EOF_
)
    else
        cfsan_snp_pipeline merge_vcfs + force_flag + -o "$workDir/snpma.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile" 2>&1 | tee $logDir/mergeVcf.log
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
    cfsan_snp_pipeline distance + force_flag + -p "$workDir/snp_distance_pairwise.tsv" -m "$workDir/snp_distance_matrix.tsv" "$workDir/snpma.fasta"
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
    cfsan_snp_pipeline distance + force_flag + -p "$workDir/snp_distance_pairwise.tsv" -m "$workDir/snp_distance_matrix.tsv" "$workDir/snpma.fasta"
_EOF_
)
else
    cfsan_snp_pipeline distance + force_flag + -p "$workDir/snp_distance_pairwise.tsv" -m "$workDir/snp_distance_matrix.tsv" "$workDir/snpma.fasta" 2>&1 | tee $logDir/calcSnpDistances.log
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
    cfsan_snp_pipeline merge_sites + force_flag + -n var.flt_preserved.vcf -o "$workDir/snplist_preserved.txt" $CreateSnpList_ExtraParams sample_dirs_file "$filteredSampleDirsFile2"
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
    cfsan_snp_pipeline merge_sites + force_flag + -n var.flt_preserved.vcf -o "$workDir/snplist_preserved.txt" $CreateSnpList_ExtraParams sample_dirs_file "$filteredSampleDirsFile2"
_EOF_
)
else
    cfsan_snp_pipeline merge_sites + force_flag + -n var.flt_preserved.vcf -o "$workDir/snplist_preserved.txt" $CreateSnpList_ExtraParams sample_dirs_file "$filteredSampleDirsFile2" 2>&1 | tee $logDir/snpList_preserved.log
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
    cfsan_snp_pipeline call_consensus + force_flag + -l "$workDir/snplist_preserved.txt" -o "\$sampleDir/consensus_preserved.fasta" -e "\$sampleDir/var.flt_removed.vcf" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams  --vcfFileName consensus_preserved.vcf "\$sampleDir/reads.all.pileup"
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
    cfsan_snp_pipeline call_consensus + force_flag + -l "$workDir/snplist_preserved.txt" -o "\$sampleDir/consensus_preserved.fasta" -e "\$sampleDir/var.flt_removed.vcf" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams  --vcfFileName consensus_preserved.vcf "\$sampleDir/reads.all.pileup"
_EOF_
)
else
    if "$MaxConcurrentCallConsensus" != "":
        numCallConsensusCores=$MaxConcurrentCallConsensus
    else
        numCallConsensusCores=$numCores
    fi
    nl sample_dirs_file | xargs -n 2 -P $numCallConsensusCores bash -c 'set -o pipefail; cfsan_snp_pipeline call_consensus + force_flag + -l "$workDir/snplist_preserved.txt" -o "$1/consensus_preserved.fasta" -e "$1/var.flt_removed.vcf" --vcfRefName "$referenceFileName" $CallConsensus_ExtraParams  --vcfFileName consensus_preserved.vcf "$1/reads.all.pileup" 2>&1 | tee $logDir/callConsensus_preserved.log-$0'
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
    cfsan_snp_pipeline snp_matrix + force_flag + -c consensus_preserved.fasta -o "$workDir/snpma_preserved.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile2"
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
    cfsan_snp_pipeline snp_matrix + force_flag + -c consensus_preserved.fasta -o "$workDir/snpma_preserved.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile2"
_EOF_
)
else
    cfsan_snp_pipeline snp_matrix + force_flag + -c consensus_preserved.fasta -o "$workDir/snpma_preserved.fasta" $CreateSnpMatrix_ExtraParams "$filteredSampleDirsFile2" 2>&1 | tee $logDir/snpMatrix_preserved.log
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
    cfsan_snp_pipeline snp_reference + force_flag + -l "$workDir/snplist_preserved.txt" -o "$workDir/referenceSNP_preserved.fasta" $CreateSnpReferenceSeq_ExtraParams reference_file_path
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
    cfsan_snp_pipeline snp_reference + force_flag + -l "$workDir/snplist_preserved.txt" -o "$workDir/referenceSNP_preserved.fasta" $CreateSnpReferenceSeq_ExtraParams reference_file_path
_EOF_
)
else
    cfsan_snp_pipeline snp_reference + force_flag + -l "$workDir/snplist_preserved.txt" -o "$workDir/referenceSNP_preserved.fasta" $CreateSnpReferenceSeq_ExtraParams reference_file_path 2>&1 | tee $logDir/snpReference_preserved.log
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
        cfsan_snp_pipeline merge_vcfs + force_flag + -n consensus_preserved.vcf -o "$workDir/snpma_preserved.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile2"
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
        cfsan_snp_pipeline merge_vcfs + force_flag + -n consensus_preserved.vcf -o "$workDir/snpma_preserved.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile2"
_EOF_
)
    else
        cfsan_snp_pipeline merge_vcfs + force_flag + -n consensus_preserved.vcf -o "$workDir/snpma_preserved.vcf" $MergeVcf_ExtraParams "$filteredSampleDirsFile2" 2>&1 | tee $logDir/mergeVcf_preserved.log
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
    cfsan_snp_pipeline distance + force_flag + -p "$workDir/snp_distance_pairwise_preserved.tsv" -m "$workDir/snp_distance_matrix_preserved.tsv" "$workDir/snpma_preserved.fasta"
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
    cfsan_snp_pipeline distance + force_flag + -p "$workDir/snp_distance_pairwise_preserved.tsv" -m "$workDir/snp_distance_matrix_preserved.tsv" "$workDir/snpma_preserved.fasta"
_EOF_
)
else
    cfsan_snp_pipeline distance + force_flag + -p "$workDir/snp_distance_pairwise_preserved.tsv" -m "$workDir/snp_distance_matrix_preserved.tsv" "$workDir/snpma_preserved.fasta" 2>&1 | tee $logDir/calcSnpDistances_preserved.log
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
    cfsan_snp_pipeline collect_metrics -o "\$sampleDir/metrics" $CollectSampleMetrics_ExtraParams "\$sampleDir"  reference_file_path
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
    cfsan_snp_pipeline collect_metrics -o "\$sampleDir/metrics" $CollectSampleMetrics_ExtraParams "\$sampleDir"  reference_file_path
_EOF_
)
else
    if "$MaxConcurrentCollectSampleMetrics" != "":
        numCollectSampleMetricsCores=$MaxConcurrentCollectSampleMetrics
    else
        numCollectSampleMetricsCores=$numCores
    fi
    nl sample_dirs_file | xargs -n 2 -P $numCollectSampleMetricsCores bash -c 'set -o pipefail; cfsan_snp_pipeline collect_metrics -o "$1/metrics" $CollectSampleMetrics_ExtraParams "$1" reference_file_path 2>&1 | tee $logDir/collectSampleMetrics.log-$0'
fi


print("\nStep 13 - Combine the metrics across all samples into the metrics table")
command_line = 'cfsan_snp_pipeline combine_metrics -n metrics -o "$workDir/metrics.tsv" $CombineSampleMetrics_ExtraParams sample_dirs_file'
log_file = '$logDir/combineSampleMetrics.log'
combine_metrics_job_id = runner.run(command_line, "combineMetrics", log_file, wait_for_array=[collectSampleMetricsJobArray])

# Step 14.2 - Notify user of any non-fatal errors accumulated during processing
if os.path.getsize(error_output_file) > 0 and not stop_on_error:
    print("\nThere were errors processing some samples.\nSee the log file %s for a summary of errors." % error_output_file, file=sys.stderr)

"""
