#!/usr/bin/env python

"""This module is part of the CFSAN SNP Pipeline. It contains the code to
run all the sequential steps of the pipeline in the correct order.
"""

from __future__ import print_function
from __future__ import absolute_import

import os
import psutil
import re
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

# Globals
log_dir = ""
job_queue_mgr = None

def progress(message):
    """Print a progress message.

    Parameters
    ----------
    message : str
        The text of the progress message.
    """
    if job_queue_mgr is None: # local
        print("")
        print("*" * 80)
        print(message)
        print("*" * 80)
    else:
        print("Submitting " + message)



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


def get_sorted_sample_dirs_fastq_sizes(samples_parent_dir):
    """Given a parent directory containing multiple sample directories, return a list
    of the sample subdirectories with the size of the fastq files in each.

    Parameters
    ----------
    samples_parent_dir : str
        Path to the parent directory containing multiple sample directories

    Returns
    -------
    dir_sizes : list of tuples
        Sorted list of (size, path) tuples, largest first.
    """
    sub_dirs = [os.path.join(samples_parent_dir, d) for d in os.listdir(samples_parent_dir) if os.path.isdir(os.path.join(samples_parent_dir, d))]
    dir_sizes = []
    for d in sub_dirs:
        size = sum(map(os.path.getsize, fastq.list_fastq_files(d)))
        dir_sizes.append((size, d))
    dir_sizes.sort(reverse=True)
    return dir_sizes


def persist_sorted_sample_dirs_file(samples_parent_dir, sample_dirs_file):
    """Given a parent directory containing multiple sample directories, create
    a file listing the sample directories sorted by fastq file size, largest first.

    Subdirectories without fastq files are ignored.

    Parameters
    ----------
    samples_parent_dir : str
        Path to the parent directory containing multiple sample directories
    sample_dirs_file : str
        Path to the file of sample directories which will be created
    """
    dir_sizes = get_sorted_sample_dirs_fastq_sizes(samples_parent_dir)

    # ignore subdirectories with no fastq files
    dir_sizes = [(size, path) for size, path in dir_sizes if size > 0]

    with open(sample_dirs_file, 'w') as f:
        for size, directory in dir_sizes:
            print(directory, file=f)


def configure_process_threads(extra_params_env_var, threads_option, default_threads_per_process, max_cpu_cores):
    """Detect the user-configured number of allowed threads for a process and compute the
    corresponding number of allowed processes given the maximum allowed number of CPUs.

    The named environment variable is parsed to detect a user-setting for the number of threads per process.
    If not found, the environment variable is modified with a default setting.

    If the user requests a number of threads that is greater than the max_cpu_cores, the number of allowed
    threads will be set to the max_cpu_cores.

    Parameters
    ----------
    extra_params_env_var : str
        Name of an environment variable the user can set with embedded command line options
        to configure the number of threads for a process.
    threads_option : str
        The exact spelling of a command line option to set the number of threads, for example "-n".
    default_threads_per_process : int
        The number of threads to use if the user did not set a preferred value in the named
        environment variable.
    max_cpu_cores : int or None
        The maximum allowed number of CPU cores to consume by all instances of the process.
        If set the None, it implies no limit to the number of CPU cores that can be used.

    Returns
    -------
    max_processes : int or None
        The computed maximum number of allowed concurrent processes, or None to allow unlimited.
    threads_per_process : int
        The number of threads per process instance.  This will either be the value previously
        configured by the user in the named environment variable, or the default value passed
        as an argument.

    Examples
    --------
    # env var not set, max CPU not set
    >>> _ = os.environ.pop("SmaltAlign_ExtraParams", None)
    >>> configure_process_threads("SmaltAlign_ExtraParams", "-n", 8, None)
    (None, 8)
    >>> os.environ["SmaltAlign_ExtraParams"]
    '-n 8'

    # env var not set, max CPU set, not multiple of threads
    >>> _ = os.environ.pop("SmaltAlign_ExtraParams", None)
    >>> configure_process_threads("SmaltAlign_ExtraParams", "-n", 8, 20)
    (2, 8)
    >>> os.environ["SmaltAlign_ExtraParams"]
    '-n 8'

    # env var not set, max CPU set, multiple of threads
    >>> _ = os.environ.pop("SmaltAlign_ExtraParams", None)
    >>> configure_process_threads("SmaltAlign_ExtraParams", "-n", 8, 24)
    (3, 8)
    >>> os.environ["SmaltAlign_ExtraParams"]
    '-n 8'

    # env var exists but option not set, max CPU not set
    >>> os.environ["SmaltAlign_ExtraParams"] = "--version"
    >>> configure_process_threads("SmaltAlign_ExtraParams", "-n", 9, None)
    (None, 9)
    >>> os.environ["SmaltAlign_ExtraParams"]
    '--version -n 9'

    # env var exists, option set, max CPU not set
    >>> os.environ["SmaltAlign_ExtraParams"] = "--version -n 10"
    >>> configure_process_threads("SmaltAlign_ExtraParams", "-n", 8, None)
    (None, 10)
    >>> os.environ["SmaltAlign_ExtraParams"]
    '--version -n 10'

    # env var exists, option set, max CPU set, not multiple of threads
    >>> os.environ["SmaltAlign_ExtraParams"] = "-n 10 --version"
    >>> configure_process_threads("SmaltAlign_ExtraParams", "-n", 8, 22)
    (2, 10)
    >>> os.environ["SmaltAlign_ExtraParams"]
    '-n 10 --version'

    # env var exists, option set, max CPU set, multiple of threads
    >>> os.environ["SmaltAlign_ExtraParams"] = "-n 10 --version"
    >>> configure_process_threads("SmaltAlign_ExtraParams", "-n", 8, 30)
    (3, 10)
    >>> os.environ["SmaltAlign_ExtraParams"]
    '-n 10 --version'

    # env var exists, option set, max CPU set, max CPU == desired threads
    >>> os.environ["SmaltAlign_ExtraParams"] = "-n 4 --version"
    >>> configure_process_threads("SmaltAlign_ExtraParams", "-n", 8, 4)
    (1, 4)
    >>> os.environ["SmaltAlign_ExtraParams"]
    '-n 4 --version'

    # env var exists, option set, max CPU set, max CPU less than desired threads
    >>> os.environ["SmaltAlign_ExtraParams"] = "-n 10 --version"
    >>> configure_process_threads("SmaltAlign_ExtraParams", "-n", 8, 2)
    (1, 2)
    >>> os.environ["SmaltAlign_ExtraParams"]
    '-n 2 --version'

    # env var exists but option not set, max CPU less than default number of threads
    >>> os.environ["SmaltAlign_ExtraParams"] = "--version"
    >>> configure_process_threads("SmaltAlign_ExtraParams", "-n", 9, 7)
    (1, 7)
    >>> os.environ["SmaltAlign_ExtraParams"]
    '--version -n 7'
    """
    regex_str = threads_option + "[ \t]*([0-9]+)"
    extra_params = os.environ.get(extra_params_env_var, "")
    match = re.search(regex_str, extra_params)
    if match:
        configured_threads_per_process = int(match.group(1))
        threads_per_process = configured_threads_per_process
    else:
        threads_per_process = default_threads_per_process

    if max_cpu_cores is None:
        max_processes = None
    elif max_cpu_cores >= threads_per_process:
        max_processes = int(max_cpu_cores / threads_per_process)
    else:
        max_processes = 1
        threads_per_process = max_cpu_cores

    threads_option += ' ' + str(threads_per_process)
    if match and threads_per_process != configured_threads_per_process:
        extra_params = re.sub(regex_str, threads_option, extra_params)
        os.environ[extra_params_env_var] = extra_params
    elif not match:
        if extra_params:
            extra_params += ' '
        os.environ[extra_params_env_var] = extra_params + threads_option

    return max_processes, threads_per_process


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
    global log_dir
    global job_queue_mgr

    # Where are we running: grid, torque, or None (local)
    job_queue_mgr = args.jobQueueMgr

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
    force_flag = " -f " if args.forceFlag else " "

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

    # How many CPU cores can we use?
    max_cpu_cores = config_params.get("MaxCpuCores", None)
    if max_cpu_cores == "":
        max_cpu_cores = None
    if max_cpu_cores:
        try:
            max_cpu_cores = int(max_cpu_cores)
            if max_cpu_cores < 1:
                utils.fatal_error("Config file error in MaxCpuCores parameter: %s is less than one." % max_cpu_cores)
        except ValueError:
            utils.fatal_error("Config file error in MaxCpuCores parameter: %s is not a valid number." % max_cpu_cores)

    if job_queue_mgr is None: # workstation
        num_local_cpu_cores = psutil.cpu_count()
        max_cpu_cores = min(num_local_cpu_cores, max_cpu_cores) if max_cpu_cores else num_local_cpu_cores

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
    os.environ["FilterRegions_ExtraParams"] = config_params.get("FilterRegions_ExtraParams", "")
    os.environ["MergeSites_ExtraParams"] = config_params.get("MergeSites_ExtraParams", "")
    os.environ["CallConsensus_ExtraParams"] = config_params.get("CallConsensus_ExtraParams", "")
    os.environ["CreateSnpMatrix_ExtraParams"] = config_params.get("CreateSnpMatrix_ExtraParams", "")
    os.environ["BcftoolsMerge_ExtraParams"] = config_params.get("BcftoolsMerge_ExtraParams", "")
    os.environ["CreateSnpReferenceSeq_ExtraParams"] = config_params.get("CreateSnpReferenceSeq_ExtraParams", "")
    os.environ["CollectSampleMetrics_ExtraParams"] = config_params.get("CollectSampleMetrics_ExtraParams", "")
    os.environ["CombineSampleMetrics_ExtraParams"] = config_params.get("CombineSampleMetrics_ExtraParams", "")
    #os.environ["GridEngine_PEname"] = config_params.get("GridEngine_PEname", "")

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

    # Process the sample directory command line option
    # TODO: detect broken fastq symlinks
    if args.samplesDir:
        samples_parent_dir = args.samplesDir.rstrip('/') # strip trailing slash
        if not utils.verify_non_empty_directory("Samples directory", samples_parent_dir):
            sys.exit(1)

        # verify at least one of the subdirectories contains fastq files.
        dir_sizes = get_sorted_sample_dirs_fastq_sizes(samples_parent_dir)
        dir_sizes = [(size, path) for size, path in dir_sizes if size > 0]
        if len(dir_sizes) == 0:
            utils.fatal_error("Samples directory %s does not contain subdirectories with fastq files." % samples_parent_dir)

        sample_dirs_file = os.path.join(work_dir, "sampleDirectories.txt")
        persist_sorted_sample_dirs_file(samples_parent_dir, sample_dirs_file)

    # Process the file of sample directories command line option
    # TODO: detect broken fastq symlinks
    if args.samplesFile:
        sample_dirs_file = args.samplesFile
        if not os.path.isfile(sample_dirs_file):
            utils.fatal_error("Error: the file of samples directories, %s does not exist." % sample_dirs_file)
        if os.path.getsize(sample_dirs_file) == 0:
            utils.fatal_error("Error: the file of samples directories, %s is empty." % sample_dirs_file)
        rewrite_cleansed_file_of_sample_dirs(sample_dirs_file, os.path.join(work_dir, "sampleDirectories.txt"))
        sample_dirs_file = os.path.join(work_dir, "sampleDirectories.txt")
        validate_file_of_sample_dirs(sample_dirs_file)

    with open(sample_dirs_file) as f:
       sample_dirs_list = f.read().splitlines()
    sample_count = len(sample_dirs_list)

    # --------------------------------------------------------
    if job_queue_mgr is None:
        progress("Step 1 - Prep work")
    else:
        print("Step 1 - Prep work")

    # --------------------------------------------------------
    # Mirror the input reference and samples if requested
    # TODO: make this a pure python solution
    if args.mirror:
        if args.mirror == "soft":
            # soft link, subsequent freshness checks use the timestamp of original file, not the soft link
            mirror_flag = " -s "
        elif args.mirror == "hard":
            # hard link, automatically preserves attributes of the original file
            mirror_flag = " -l "
        else:
            # regular copy, -p explicitly preserves attributes of the original file
            mirror_flag = " -p "

        # flush stdout to keep the unbuffered stderr in chronological order with stdout
        sys.stdout.flush()

        # Mirror/link the reference
        work_reference_dir = os.path.join(work_dir, "reference")
        utils.mkdir_p(work_reference_dir)
        src_reference_file = os.path.abspath(reference_file_path)
        cmd = "cp -v -u -f" + mirror_flag + src_reference_file + ' ' + work_reference_dir
        subprocess.check_call(cmd, shell=True)

        # since we mirrored the reference, we need to update our reference location
        reference_file_path = os.path.join(work_reference_dir, reference_file_name)

        # Mirror/link the samples
        work_samples_parent_dir = os.path.join(work_dir, "samples")
        for directory in sample_dirs_list:
            basedir = os.path.basename(directory)
            work_sample_dir = os.path.join(work_samples_parent_dir, basedir)
            utils.mkdir_p(work_sample_dir)
            src_sample_dir = os.path.abspath(directory)
            # copy without stderr message and without exit error code because the fastq or fq files might not exist
            cmd = "cp -r -v -u -f" + mirror_flag + src_sample_dir + "/*.fastq* " + work_sample_dir + " 2> /dev/null || true"
            subprocess.check_call(cmd, shell=True)
            cmd = "cp -r -v -u -f" + mirror_flag + src_sample_dir + "/*.fq* " + work_sample_dir + " 2> /dev/null || true"
            subprocess.check_call(cmd, shell=True)

        # since we mirrored the samples, we need to update our sorted list of samples
        sample_dirs_file = os.path.join(work_dir, "sampleDirectories.txt")
        persist_sorted_sample_dirs_file(work_samples_parent_dir, sample_dirs_file)

        # refresh the list of sample dirs -- now in sorted order
        with open(sample_dirs_file) as f:
           sample_dirs_list = f.read().splitlines()

    # get the *.fastq or *.fq files in each sample directory, possibly compresessed, on one line per sample, ready to feed to bowtie
    sample_full_path_names_file = os.path.join(work_dir, "sampleFullPathNames.txt")
    with open(sample_full_path_names_file, 'w') as f:
        for directory in sample_dirs_list:
            file_list = fastq.list_fastq_files(directory)
            print(' '.join(file_list), file=f)

    # Initialize the job runner
    if job_queue_mgr is None:
        runner = JobRunner("local")
    elif job_queue_mgr == "grid":
        strip_job_array_suffix = config_params.get("GridEngine_StripJobArraySuffix", "true").lower()
        runner = JobRunner(job_queue_mgr, strip_job_array_suffix == "true")
    else:
        strip_job_array_suffix = config_params.get("Torque_StripJobArraySuffix", "false").lower()
        runner = JobRunner(job_queue_mgr, strip_job_array_suffix == "true")

    progress("Step 2 - Index the reference")
    log_file = os.path.join(log_dir, "indexRef.log")
    command_line = "cfsan_snp_pipeline index_ref" + force_flag + reference_file_path
    job_id_index_ref = runner.run(command_line, "indexRef", log_file)

    progress("Step 3 - Map the sample reads to the reference")
    # Parse the user-specified aligner parameters to find the number of CPU cores requested, for example, "-p 16" or "-n 16"
    # Set the default number of CPU cores if the user did not configure a value.
    if snp_pipeline_aligner == "smalt":
        extra_params_env_var = "SmaltAlign_ExtraParams"
        threads_option = "-n"
    else:
        extra_params_env_var = "Bowtie2Align_ExtraParams"
        threads_option = "-p"

    max_processes, threads_per_process = configure_process_threads(extra_params_env_var, threads_option, 8, max_cpu_cores)

    parallel_environment = config_params.get("GridEngine_PEname", None)
    log_file = os.path.join(log_dir, "mapReads.log")
    command_line = "cfsan_snp_pipeline map_reads" + force_flag + reference_file_path + " {1} {2}"
    job_id_map_reads = runner.run_array(command_line, "mapReads", log_file, sample_full_path_names_file, max_processes=max_processes, wait_for=[job_id_index_ref], threads=threads_per_process, parallel_environment=parallel_environment)

    progress("Step 4 - Find sites with SNPs in each sample")
    if job_queue_mgr in ["grid", "torque"]:
        time.sleep(1.0 + float(sample_count) / 150) # workaround torque bug when submitting two large consecutive array jobs, potential bug for grid

    log_file = os.path.join(log_dir, "callSites.log")
    command_line = "cfsan_snp_pipeline call_sites" + force_flag + reference_file_path + " {1}"
    job_id_call_sites = runner.run_array(command_line, "callSites", log_file, sample_dirs_file, max_processes=max_cpu_cores, wait_for_array=[job_id_map_reads], slot_dependency=True)

    progress("Step 5 - Filter abnormal SNP regions")
    log_file = os.path.join(log_dir, "filterRegions.log")
    extra_params = os.environ.get("FilterRegions_ExtraParams", "")
    command_line = "cfsan_snp_pipeline filter_regions" + force_flag + "-n var.flt.vcf " + sample_dirs_file + ' ' + reference_file_path + ' ' + extra_params
    job_id_filter_regions = runner.run(command_line, "filterRegions", log_file, wait_for_array=[job_id_call_sites])

    # Starting from here, there are 2 threads:
    # Thread X.1: the thread processing the original VCF files and corresponding downstream results
    # Thread X.2: the thread processing the preserved VCF files and corresponding downstream results

    progress("Step 6.1 - Merge the SNP sites across all samples into the SNP list file")
    # The create_snp_list process creates the filtered list of sample directories.  It is the list of samples not having excessive snps.
    # When running on a workstation, the file exists at this point during the script execution, but on grid or torque, it has not yet been created. However,
    # we know the path to the file regardless of whether it exists yet.
    filtered_sample_dirs_file = sample_dirs_file + ".OrigVCF.filtered"
    # touch $filtered_sample_dirs_file # TODO: why was this touch here in the old run_snp_pipeline.sh script?
    log_file = os.path.join(log_dir, "mergeSites.log")
    output_file = os.path.join(work_dir, "snplist.txt")
    extra_params = os.environ.get("MergeSites_ExtraParams", "")
    command_line = "cfsan_snp_pipeline merge_sites" + force_flag + "-n var.flt.vcf -o " + output_file + ' ' + extra_params + ' ' + sample_dirs_file + ' ' + filtered_sample_dirs_file
    job_id_merge_sites = runner.run(command_line, "mergeSites", log_file, wait_for=[job_id_filter_regions])

    progress("Step 6.2 - Merge the SNP sites across all samples into the SNP list file")
    # Create another copy of sample directories file, for the thread processing preserved snp files.
    filtered_sample_dirs_file2 = sample_dirs_file + ".PresVCF.filtered"
    # touch $filtered_sample_dirs_file2 # TODO: why was this touch here in the old run_snp_pipeline.sh script?
    log_file = os.path.join(log_dir, "mergeSites_preserved.log")
    output_file = os.path.join(work_dir, "snplist_preserved.txt")
    extra_params = os.environ.get("MergeSites_ExtraParams", "")
    command_line = "cfsan_snp_pipeline merge_sites" + force_flag + "-n var.flt_preserved.vcf -o " + output_file + ' ' + extra_params + ' ' + sample_dirs_file + ' ' + filtered_sample_dirs_file2
    job_id_merge_sites2 = runner.run(command_line, "mergeSites_preserved", log_file, wait_for=[job_id_filter_regions])

    progress("Step 7.1 - Call the consensus SNPs for each sample")
    log_file = os.path.join(log_dir, "callConsensus.log")
    list_file = os.path.join(work_dir, "snplist.txt")
    output_file = "{1}/consensus.fasta"
    extra_params = os.environ.get("CallConsensus_ExtraParams", "")
    command_line = "cfsan_snp_pipeline call_consensus" + force_flag + "-l " + list_file + " -o " + output_file + " --vcfRefName " + reference_file_name + ' ' + extra_params + " --vcfFileName consensus.vcf {1}/reads.all.pileup"
    job_id_call_sites = runner.run_array(command_line, "callConsensus", log_file, sample_dirs_file, max_processes=max_cpu_cores, wait_for=[job_id_merge_sites])

    progress("Step 7.2 - Call the consensus SNPs for each sample")
    log_file = os.path.join(log_dir, "callConsensus_preserved.log")
    list_file = os.path.join(work_dir, "snplist_preserved.txt")
    output_file = "{1}/consensus_preserved.fasta"
    extra_params = os.environ.get("CallConsensus_ExtraParams", "")
    command_line = "cfsan_snp_pipeline call_consensus" + force_flag + "-l " + list_file + " -o " + output_file + " -e {1}/var.flt_removed.vcf --vcfRefName " + reference_file_name + ' ' + extra_params + " --vcfFileName consensus_preserved.vcf {1}/reads.all.pileup"
    job_id_call_sites2 = runner.run_array(command_line, "callConsensus", log_file, sample_dirs_file, max_processes=max_cpu_cores, wait_for=[job_id_merge_sites2])

"""

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
    cfsan_snp_pipeline snp_matrix + force_flag + -c consensus.fasta -o "$workDir/snpma.fasta" $CreateSnpMatrix_ExtraParams "$filtered_sample_dirs_file"
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
    cfsan_snp_pipeline snp_matrix + force_flag + -c consensus.fasta -o "$workDir/snpma.fasta" $CreateSnpMatrix_ExtraParams "$filtered_sample_dirs_file"
_EOF_
)
else
    cfsan_snp_pipeline snp_matrix + force_flag + -c consensus.fasta -o "$workDir/snpma.fasta" $CreateSnpMatrix_ExtraParams "$filtered_sample_dirs_file" 2>&1 | tee $logDir/snpMatrix.log
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
        cfsan_snp_pipeline merge_vcfs + force_flag + -o "$workDir/snpma.vcf" $MergeVcf_ExtraParams "$filtered_sample_dirs_file"
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
        cfsan_snp_pipeline merge_vcfs + force_flag + -o "$workDir/snpma.vcf" $MergeVcf_ExtraParams "$filtered_sample_dirs_file"
_EOF_
)
    else
        cfsan_snp_pipeline merge_vcfs + force_flag + -o "$workDir/snpma.vcf" $MergeVcf_ExtraParams "$filtered_sample_dirs_file" 2>&1 | tee $logDir/mergeVcf.log
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
    cfsan_snp_pipeline snp_matrix + force_flag + -c consensus_preserved.fasta -o "$workDir/snpma_preserved.fasta" $CreateSnpMatrix_ExtraParams "$filtered_sample_dirs_file2"
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
    cfsan_snp_pipeline snp_matrix + force_flag + -c consensus_preserved.fasta -o "$workDir/snpma_preserved.fasta" $CreateSnpMatrix_ExtraParams "$filtered_sample_dirs_file2"
_EOF_
)
else
    cfsan_snp_pipeline snp_matrix + force_flag + -c consensus_preserved.fasta -o "$workDir/snpma_preserved.fasta" $CreateSnpMatrix_ExtraParams "$filtered_sample_dirs_file2" 2>&1 | tee $logDir/snpMatrix_preserved.log
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
        cfsan_snp_pipeline merge_vcfs + force_flag + -n consensus_preserved.vcf -o "$workDir/snpma_preserved.vcf" $MergeVcf_ExtraParams "$filtered_sample_dirs_file2"
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
        cfsan_snp_pipeline merge_vcfs + force_flag + -n consensus_preserved.vcf -o "$workDir/snpma_preserved.vcf" $MergeVcf_ExtraParams "$filtered_sample_dirs_file2"
_EOF_
)
    else
        cfsan_snp_pipeline merge_vcfs + force_flag + -n consensus_preserved.vcf -o "$workDir/snpma_preserved.vcf" $MergeVcf_ExtraParams "$filtered_sample_dirs_file2" 2>&1 | tee $logDir/mergeVcf_preserved.log
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
    collectSampleMetricsJobId=$(echo | qsub -terse -t 1-$sample_count $GridEngine_QsubExtraParams << _EOF_
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
    collectSampleMetricsJobId=$(echo | qsub -t 1-$sample_count $Torque_QsubExtraParams << _EOF_
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
