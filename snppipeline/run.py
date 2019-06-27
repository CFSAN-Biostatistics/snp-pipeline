#!/usr/bin/env python

"""This module is part of the CFSAN SNP Pipeline. It contains the code to
run all the sequential steps of the pipeline in the correct order.
"""

from __future__ import print_function
from __future__ import absolute_import

from jobrunner import JobRunner
from jobrunner import JobRunnerException
import os
import psutil
import shutil
import subprocess
import sys
import time
import traceback

from snppipeline import command
from snppipeline import fastq
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

    trace_entries = traceback.extract_tb(exc_traceback)
    file_name, line_number, function_name, code_text = trace_entries[-1]
    exc_type_name = exc_type.__name__

    stop_on_error_env = os.environ.get("StopOnSampleError")
    stop_on_error = stop_on_error_env is None or stop_on_error_env == "true"

    error_output_file = os.environ.get("errorOutputFile")

    # Error code 98 indicates an error affecting a single sample when the pipeline is configured to ignore such errors.
    # This error has already been reported, but it should not stop execution.
    if error_code == 98 and not stop_on_error:
        return

    # When configured to ignore single-sample errors, try to continue execution if
    # xargs returns 123 when executing an array of jobs, one per sample
    if error_code == 123 and not stop_on_error:
        return

    # A subprocess failed and was already trapped.
    # Actually, we cannot be 100% certain the error was trapped if the error code is 123.  This
    # indicates an error in an array of sample jobs launched in parallel by xargs; but since
    # StopOnSampleError is true, we will assume the error was already trapped.
    if error_code == 100 or (stop_on_error and error_code == 123):
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
    in the snp-pipeline driver itself.  This function is installed to
    replace the default system exception handler with sys.excepthook.

    This function is installed as the JobRunner exception handler
    to handle CalledProcessError exceptions.  In this case, the handler can
    return and try to continue execution if the handler decides to do so.

    When this function is executed by the Python as the sys.excepthook,
    there is no hope of continuing execition.  The stack has already
    been unwound and this function is called immediately before exit.
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
        if os.environ.get("StopOnSampleError") == "true":
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


def validate_properties(properties):
    """Verify the previously read properties are from the latest version of the properties file.

    Prints an error message to stderr and exits if the properties file is not the latest.

    Parameters
    ----------
    properties : dict
        Dictionary mapping property names to str values.  This is obtained by calling the utils.read_properties() function.
    """
    message = "Your configuration file does not match this version of the SNP Pipeline.\n" \
              "Please start with a fresh configuration file.\n" \
              "Use this command:\n" \
              "    cfsan_snp_pipeline data configurationFile  # this will overwrite an existing snppipeline.conf file!!!"

    required_parameters = ["StopOnSampleError", "MaxCpuCores", "CpuCoresPerProcessOnHPC", "CpuCoresPerProcessOnWorkstation",
                           "MaxSnps", "RemoveDuplicateReads", "EnableLocalRealignment",
                           ]
    for p in required_parameters:
        if p not in properties:
            utils.fatal_error(message)


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

    start_time = time.time()
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
    log_dir = os.path.join(work_dir, "logs-" + run_time_stamp)
    try:
        utils.mkdir_p(log_dir)
    except OSError as exc:
        utils.fatal_error("Error: could not create the logs directory %s" % log_dir)
    if not utils.is_directory_writeable(work_dir):
        utils.fatal_error("Error: logs directory % is not writable." % log_dir)

    # Handle configuration file, use the specified file, or create a default file
    if args.configFile:
        config_file_path = args.configFile
        if not os.path.isfile(config_file_path):
            utils.fatal_error("Error: configuration file %s does not exist." % config_file_path)
        if os.path.getsize(config_file_path) == 0:
            utils.fatal_error("Error: configuration file %s is empty." % config_file_path)

        shutil.copy2(config_file_path, log_dir)  # copy2 tries to preserve timestamps
        config_params = utils.read_properties(config_file_path, recognize_vars=True)
        validate_properties(config_params)
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
    # The environment variable is used by called processes
    stop_on_error = config_params.get("StopOnSampleError", "").lower() or "true"
    os.environ["StopOnSampleError"] = stop_on_error

    # Convert the stop_on_error flag to boolean for internal use in this function
    stop_on_error = stop_on_error == "true"

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

    # How many CPU cores per process?
    if job_queue_mgr is None: # workstation
        cpu_cores_per_process = config_params.get("CpuCoresPerProcessOnWorkstation", None)
        if cpu_cores_per_process:
            try:
                cpu_cores_per_process = int(cpu_cores_per_process)
                if cpu_cores_per_process < 1:
                    utils.fatal_error("Config file error in CpuCoresPerProcessOnWorkstation parameter: %s is less than one." % cpu_cores_per_process)
            except ValueError:
                utils.fatal_error("Config file error in CpuCoresPerProcessOnWorkstation parameter: %s is not a valid number." % cpu_cores_per_process)
        else:
            cpu_cores_per_process = min(num_local_cpu_cores, max_cpu_cores)
    else: # HPC
        cpu_cores_per_process = config_params.get("CpuCoresPerProcessOnHPC", None)
        if not cpu_cores_per_process:
            utils.fatal_error("Config file error. CpuCoresPerProcessOnHPC parameter must be set to a value.")
        else:
            try:
                cpu_cores_per_process = int(cpu_cores_per_process)
                if cpu_cores_per_process < 1:
                    utils.fatal_error("Config file error in CpuCoresPerProcessOnHPC parameter: %s is less than one." % cpu_cores_per_process)
            except ValueError:
                utils.fatal_error("Config file error in CpuCoresPerProcessOnHPC parameter: %s is not a valid number." % cpu_cores_per_process)

    # Put the configuration parameters into the process environment variables
    os.environ["Bowtie2Build_ExtraParams"] = config_params.get("Bowtie2Build_ExtraParams", "")
    os.environ["SmaltIndex_ExtraParams"] = config_params.get("SmaltIndex_ExtraParams", "")
    os.environ["CreateSequenceDictionary_ExtraParams"] = config_params.get("CreateSequenceDictionary_ExtraParams", "")
    os.environ["SamtoolsFaidx_ExtraParams"] = config_params.get("SamtoolsFaidx_ExtraParams", "")
    os.environ["Bowtie2Align_ExtraParams"] = config_params.get("Bowtie2Align_ExtraParams", "")
    os.environ["SmaltAlign_ExtraParams"] = config_params.get("SmaltAlign_ExtraParams", "")
    os.environ["SamtoolsSamFilter_ExtraParams"] = config_params.get("SamtoolsSamFilter_ExtraParams", "")
    os.environ["SamtoolsSort_ExtraParams"] = config_params.get("SamtoolsSort_ExtraParams", "")
    os.environ["SamtoolsIndex_ExtraParams"] = config_params.get("SamtoolsIndex_ExtraParams", "")
    os.environ["RemoveDuplicateReads"] = config_params.get("RemoveDuplicateReads", "").lower() or "true"
    os.environ["PicardJvm_ExtraParams"] = config_params.get("PicardJvm_ExtraParams", "")
    os.environ["PicardMarkDuplicates_ExtraParams"] = config_params.get("PicardMarkDuplicates_ExtraParams", "")
    os.environ["EnableLocalRealignment"] = config_params.get("EnableLocalRealignment", "").lower() or "true"
    os.environ["GatkJvm_ExtraParams"] = config_params.get("GatkJvm_ExtraParams", "")
    os.environ["RealignerTargetCreator_ExtraParams"] = config_params.get("RealignerTargetCreator_ExtraParams", "")
    os.environ["IndelRealigner_ExtraParams"] = config_params.get("IndelRealigner_ExtraParams", "")
    os.environ["SamtoolsMpileup_ExtraParams"] = config_params.get("SamtoolsMpileup_ExtraParams", "")
    os.environ["VarscanMpileup2snp_ExtraParams"] = config_params.get("VarscanMpileup2snp_ExtraParams", "")
    os.environ["VarscanJvm_ExtraParams"] = config_params.get("VarscanJvm_ExtraParams", "")
    os.environ["FilterRegions_ExtraParams"] = config_params.get("FilterRegions_ExtraParams", "")
    os.environ["MergeSites_ExtraParams"] = config_params.get("MergeSites_ExtraParams", "")
    os.environ["CallConsensus_ExtraParams"] = config_params.get("CallConsensus_ExtraParams", "")
    os.environ["SnpMatrix_ExtraParams"] = config_params.get("SnpMatrix_ExtraParams", "")
    os.environ["BcftoolsMerge_ExtraParams"] = config_params.get("BcftoolsMerge_ExtraParams", "")
    os.environ["SnpReference_ExtraParams"] = config_params.get("SnpReference_ExtraParams", "")
    os.environ["MergeVcfs_ExtraParams"] = config_params.get("MergeVcfs_ExtraParams", "")
    os.environ["CollectMetrics_ExtraParams"] = config_params.get("CollectMetrics_ExtraParams", "")
    os.environ["CombineMetrics_ExtraParams"] = config_params.get("CombineMetrics_ExtraParams", "")

    # Verify the dependencies are available on the path
    print("Checking dependencies...")

    dependencies = ["cfsan_snp_pipeline", snp_pipeline_aligner, "java", "tabix", "bgzip", "bcftools"]
    found_all_dependencies = True
    for executable in dependencies:
        if not utils.which(executable):
            utils.report_error(executable + " is not on the path")
            found_all_dependencies = False

    if not utils.which("samtools"):
        utils.report_error("samtools is not on the path")
        found_all_dependencies = False
    else:
        version_str = utils.extract_version_str("SAMtools", "samtools 2>&1 > /dev/null")
        samtools_version = version_str.split()[-1] # just the number
        if samtools_version < "1.4":
            utils.report_error("The installed %s is not supported.  Version 1.4 or higher is required." % version_str)
            found_all_dependencies = False

    jar_file_path = utils.find_path_in_path_list("VarScan", "CLASSPATH")
    if jar_file_path:
        stdout = command.run("java -jar " + jar_file_path + " 2>&1")
    if not jar_file_path or "error" in stdout.lower():
        utils.report_error("CLASSPATH is not configured with the path to VarScan.jar")
        found_all_dependencies = False

    picard_required = os.environ["RemoveDuplicateReads"] == "true" or os.environ["EnableLocalRealignment"] == "true"
    if picard_required:
        jar_file_path = utils.find_path_in_path_list("picard", "CLASSPATH")
        if not jar_file_path:
            utils.report_error("CLASSPATH is not configured with the path to picard.jar")
            found_all_dependencies = False
        else:
            stdout = command.run("java -jar " + jar_file_path + " 2>&1")
            if stdout.lower().startswith("error"):
                utils.report_error(stdout)
                found_all_dependencies = False

    gatk_required = os.environ["EnableLocalRealignment"] == "true"
    if gatk_required:
        jar_file_path = utils.find_path_in_path_list("GenomeAnalysisTK", "CLASSPATH")
        if not jar_file_path:
            utils.report_error("CLASSPATH is not configured with the path to GenomeAnalysisTK.jar")
            found_all_dependencies = False
        else:
            stdout = command.run("java -jar " + jar_file_path + " --version 2>&1")
            if stdout.lower().startswith("error"):
                utils.report_error(stdout)
                found_all_dependencies = False
            else:
                stdout = command.run("java -jar " + jar_file_path + " -T IndelRealigner --version 2>&1")
                if "not a valid command" in stdout.lower() or "indelrealigner is no longer included" in stdout.lower():
                    utils.report_error("The installed GATK version does not support indel realignment.  Try installing an older release prior to GATK v4.")
                    found_all_dependencies = False
                elif "user error has occurred" in stdout.lower():
                    utils.report_error(stdout)
                    found_all_dependencies = False

    if not found_all_dependencies:
        utils.fatal_error("Check the SNP Pipeline installation instructions here: http://snp-pipeline.readthedocs.org/en/latest/installation.html")
    else:
        print("OK")

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
            utils.fatal_error("Error: the file of samples directories, %s, does not exist." % sample_dirs_file)
        if os.path.getsize(sample_dirs_file) == 0:
            utils.fatal_error("Error: the file of samples directories, %s, is empty." % sample_dirs_file)
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
        runner = JobRunner("local", exception_handler=handle_exception, verbose=args.verbose >= 4)
    elif job_queue_mgr == "grid":
        strip_job_array_suffix = config_params.get("GridEngine_StripJobArraySuffix", "true").lower()
        qsub_extra_params = config_params.get("GridEngine_QsubExtraParams")
        runner = JobRunner(job_queue_mgr, strip_job_array_suffix == "true", qsub_extra_params=qsub_extra_params, verbose=args.verbose >= 4)
    else:
        strip_job_array_suffix = config_params.get("Torque_StripJobArraySuffix", "false").lower()
        qsub_extra_params = config_params.get("Torque_QsubExtraParams")
        runner = JobRunner(job_queue_mgr, strip_job_array_suffix == "true", qsub_extra_params=qsub_extra_params, verbose=args.verbose >= 4)

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

    aligner_max_processes, aligner_threads_per_process = utils.configure_process_threads(extra_params_env_var, threads_option, cpu_cores_per_process, max_cpu_cores)
    samfilter_max_processes, samfilter_threads_per_process = utils.configure_process_threads("SamtoolsSamFilter_ExtraParams", ["-@", "--threads"], cpu_cores_per_process, max_cpu_cores)
    samsort_max_processes, samsort_threads_per_process = utils.configure_process_threads("SamtoolsSort_ExtraParams", ["-@", "--threads"], cpu_cores_per_process, max_cpu_cores)
    samindex_max_processes, samindex_threads_per_process = utils.configure_process_threads("SamtoolsIndex_ExtraParams", ["-@"], cpu_cores_per_process, max_cpu_cores)
    realigner_max_processes, realigner_threads_per_process = utils.configure_process_threads("RealignerTargetCreator_ExtraParams", ["-nt", "--num_threads"], cpu_cores_per_process, max_cpu_cores)

    # There are multiple processes within map_reads, each with multiple threads.
    # The CPU allocation must be enough for the process needing the largest number of threads.
    max_processes_list = [aligner_max_processes, samfilter_max_processes, samsort_max_processes, samindex_max_processes, realigner_max_processes]
    if all([i is None for i in max_processes_list]):
        max_processes = None
    else:
        max_processes = min([i for i in max_processes_list if i is not None])
    threads_per_process = max(aligner_threads_per_process, samfilter_threads_per_process, samsort_threads_per_process, samindex_threads_per_process, realigner_threads_per_process)

    parallel_environment = config_params.get("GridEngine_PEname", None)
    log_file = os.path.join(log_dir, "mapReads.log")
    command_line = "cfsan_snp_pipeline map_reads --threads " + str(threads_per_process) + force_flag + reference_file_path + " {1} {2}"
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
    # The mergeSites process creates the filtered list of sample directories.  It is the list of samples not having excessive snps.
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
    job_id_call_consensus = runner.run_array(command_line, "callConsensus", log_file, sample_dirs_file, max_processes=max_cpu_cores, wait_for=[job_id_merge_sites])

    progress("Step 7.2 - Call the consensus SNPs for each sample")
    log_file = os.path.join(log_dir, "callConsensus_preserved.log")
    list_file = os.path.join(work_dir, "snplist_preserved.txt")
    output_file = "{1}/consensus_preserved.fasta"
    extra_params = os.environ.get("CallConsensus_ExtraParams", "")
    command_line = "cfsan_snp_pipeline call_consensus" + force_flag + "-l " + list_file + " -o " + output_file + " -e {1}/var.flt_removed.vcf --vcfRefName " + reference_file_name + ' ' + extra_params + " --vcfFileName consensus_preserved.vcf {1}/reads.all.pileup"
    job_id_call_consensus2 = runner.run_array(command_line, "callConsensus_preserved", log_file, sample_dirs_file, max_processes=max_cpu_cores, wait_for=[job_id_merge_sites2])

    progress("Step 8.1 - Create the SNP matrix")
    log_file = os.path.join(log_dir, "snpMatrix.log")
    output_file = os.path.join(work_dir, "snpma.fasta")
    extra_params = os.environ.get("SnpMatrix_ExtraParams", "")
    command_line = "cfsan_snp_pipeline snp_matrix" + force_flag + "-c consensus.fasta -o " + output_file + ' ' + extra_params + ' ' + filtered_sample_dirs_file
    job_id_snp_matrix = runner.run(command_line, "snpMatrix", log_file, wait_for_array=[job_id_call_consensus])

    progress("Step 8.2 - Create the SNP matrix")
    log_file = os.path.join(log_dir, "snpMatrix_preserved.log")
    output_file = os.path.join(work_dir, "snpma_preserved.fasta")
    extra_params = os.environ.get("SnpMatrix_ExtraParams", "")
    command_line = "cfsan_snp_pipeline snp_matrix" + force_flag + "-c consensus_preserved.fasta -o " + output_file + ' ' + extra_params + ' ' + filtered_sample_dirs_file2
    job_id_snp_matrix2 = runner.run(command_line, "snpMatrix_preserved", log_file, wait_for_array=[job_id_call_consensus2])

    progress("Step 9.1 - Create the reference sequence at SNP sites")
    log_file = os.path.join(log_dir, "snpReference.log")
    list_file = os.path.join(work_dir, "snplist.txt")
    output_file = os.path.join(work_dir, "referenceSNP.fasta")
    extra_params = os.environ.get("SnpReference_ExtraParams", "")
    command_line = "cfsan_snp_pipeline snp_reference" + force_flag + "-l " + list_file + " -o " + output_file + ' ' + extra_params + ' ' + reference_file_path
    job_id_snp_reference = runner.run(command_line, "snpReference", log_file, wait_for_array=[job_id_call_consensus])

    progress("Step 9.2 - Create the reference sequence at SNP sites")
    log_file = os.path.join(log_dir, "snpReference_preserved.log")
    list_file = os.path.join(work_dir, "snplist_preserved.txt")
    output_file = os.path.join(work_dir, "referenceSNP_preserved.fasta")
    extra_params = os.environ.get("SnpReference_ExtraParams", "")
    command_line = "cfsan_snp_pipeline snp_reference" + force_flag + "-l " + list_file + " -o " + output_file + ' ' + extra_params + ' ' + reference_file_path
    job_id_snp_reference2 = runner.run(command_line, "snpReference_preserved", log_file, wait_for_array=[job_id_call_consensus2])

    progress("Step 10.1 - Merge sample VCFs to create the multi-VCF file")
    if "--vcfFileName" in os.environ.get("CallConsensus_ExtraParams", ""):
        log_file = os.path.join(log_dir, "mergeVcfs.log")
        output_file = os.path.join(work_dir, "snpma.vcf")
        extra_params = os.environ.get("MergeVcfs_ExtraParams", "")
        command_line = "cfsan_snp_pipeline merge_vcfs" + force_flag + "-o " + output_file + ' ' + extra_params + ' ' + filtered_sample_dirs_file
        job_id_merge_vcfs = runner.run(command_line, "mergeVcfs", log_file, wait_for_array=[job_id_call_consensus])
    else:
        print("Skipped per CallConsensus_ExtraParams configuration")

    progress("Step 10.2 - Merge sample VCFs to create the multi-VCF file")
    if "--vcfFileName" in os.environ.get("CallConsensus_ExtraParams", ""):
        log_file = os.path.join(log_dir, "mergeVcfs_preserved.log")
        output_file = os.path.join(work_dir, "snpma_preserved.vcf")
        extra_params = os.environ.get("MergeVcfs_ExtraParams", "")
        command_line = "cfsan_snp_pipeline merge_vcfs" + force_flag + "-n consensus_preserved.vcf -o " + output_file + ' ' + extra_params + ' ' + filtered_sample_dirs_file2
        job_id_merge_vcfs2 = runner.run(command_line, "mergeVcfs_preserved", log_file, wait_for_array=[job_id_call_consensus2])
    else:
        print("Skipped per CallConsensus_ExtraParams configuration")

    progress("Step 11.1 - Calculate SNP distance matrix")
    log_file = os.path.join(log_dir, "distance.log")
    input_file = os.path.join(work_dir, "snpma.fasta")
    pair_output_file = os.path.join(work_dir, "snp_distance_pairwise.tsv")
    matrix_output_file = os.path.join(work_dir, "snp_distance_matrix.tsv")
    command_line = "cfsan_snp_pipeline distance" + force_flag + "-p " + pair_output_file + " -m " + matrix_output_file + ' ' + input_file
    job_id_distance = runner.run(command_line, "distance", log_file, wait_for=[job_id_snp_matrix])

    progress("Step 11.2 - Calculate SNP distance matrix")
    log_file = os.path.join(log_dir, "distance_preserved.log")
    input_file = os.path.join(work_dir, "snpma_preserved.fasta")
    pair_output_file = os.path.join(work_dir, "snp_distance_pairwise_preserved.tsv")
    matrix_output_file = os.path.join(work_dir, "snp_distance_matrix_preserved.tsv")
    command_line = "cfsan_snp_pipeline distance" + force_flag + "-p " + pair_output_file + " -m " + matrix_output_file + ' ' + input_file
    job_id_distance2 = runner.run(command_line, "distance_preserved", log_file, wait_for=[job_id_snp_matrix2])

    progress("Step 12 - Collect metrics for each sample")
    log_file = os.path.join(log_dir, "collectMetrics.log")
    output_file = "{1}/metrics"
    extra_params = os.environ.get("CollectMetrics_ExtraParams", "")
    command_line = "cfsan_snp_pipeline collect_metrics" + force_flag + "-o " + output_file + ' ' + extra_params + " {1} " + reference_file_path
    job_id_collect_metrics = runner.run_array(command_line, "collectMetrics", log_file, sample_dirs_file, max_processes=max_cpu_cores, wait_for_array=[job_id_call_consensus, job_id_call_consensus2], slot_dependency=True)

    progress("Step 13 - Combine the metrics across all samples into the metrics table")
    log_file = os.path.join(log_dir, "combineMetrics.log")
    output_file = os.path.join(work_dir, "metrics.tsv")
    extra_params = os.environ.get("CombineMetrics_ExtraParams", "")
    command_line = "cfsan_snp_pipeline combine_metrics" + force_flag + "-n metrics -o " + output_file + ' ' + extra_params + ' ' + sample_dirs_file
    combine_metrics_job_id = runner.run(command_line, "combineMetrics", log_file, wait_for_array=[job_id_collect_metrics])

    # Step 14 - Notify user of any non-fatal errors accumulated during processing
    if os.path.isfile(error_output_file) and os.path.getsize(error_output_file) > 0 and not stop_on_error:
        print("\nThere were errors processing some samples.\nSee the log file %s for a summary of errors." % error_output_file, file=sys.stderr)

    # Exit here to prevent showing the "cfsan_snp_pipeline run finished" message.  The jobs are queued, not finished yet.
    if job_queue_mgr is not None: # HPC
        sys.exit(0)
    else:
        end_time = time.time()
        elapsed_time = end_time - start_time
        print("Elapsed time =", elapsed_time)
