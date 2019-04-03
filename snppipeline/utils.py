from __future__ import print_function
from __future__ import absolute_import
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

import errno
import locale
import os
import platform
import psutil
import re
import sys
import subprocess
import time
import traceback
import vcf

from snppipeline.__init__ import __version__
from snppipeline import command


#==============================================================================
#Prep work
#==============================================================================

log_verbosity = 0


#==============================================================================
#Define functions
#==============================================================================

def set_logging_verbosity(args):
    """Enable or disable logging.

    Args:
        verbose : Verbosity value, any value greater than 0 enables logging
    """
    global log_verbosity
    log_verbosity = args.verbose


def verbose_print(*args):
    """Print messages conditionally depending on the configured verbosity.
    """
    if log_verbosity > 0:
        print(*args)


def timestamp():
    """Return a timestamp string."""
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())


def program_name():
    """Return the basename of the python script being executed."""
    return os.path.basename(sys.argv[0])


def program_name_with_command():
    """Return the basename of the python script being executed with the
    subcommand if possible."""
    program = os.path.basename(sys.argv[0])
    if program == "cfsan_snp_pipeline":
        program += " " + sys.argv[1]
    return program


def command_line_short():
    """Return the command line string without the full path to the program."""
    return "%s %s" % (program_name(), " ".join(sys.argv[1:]))


def command_line_long():
    """Return the command line string with the full path to the program."""
    return " ".join(sys.argv)


def print_log_header(classpath=False):
    """Print a standardized header for the log with starting conditions.

    Parameters
    ----------
    classpath : bool
        When True, the Java CLASSPATH environment variable is logged.  This should be
        enabled when a script executes a Java program.
    """
    verbose_print("# Command           : %s" % command_line_long())
    verbose_print("# Working Directory : %s" % os.getcwd())
    pbs_jobid = os.environ.get("PBS_JOBID")
    sge_jobid = os.environ.get("JOB_ID")
    sge_task_id = os.environ.get("SGE_TASK_ID")
    if sge_task_id == "undefined":
        sge_task_id = None
    if pbs_jobid:
        verbose_print("# Job ID            : %s" % pbs_jobid)
    elif sge_jobid and sge_task_id:
        verbose_print("# Job ID            : %s[%s]" % (sge_jobid, sge_task_id))
    elif sge_jobid:
        verbose_print("# Job ID            : %s" % sge_jobid)

    verbose_print("# Hostname          : %s" % platform.node())
    locale.setlocale(locale.LC_ALL, '')
    ram_mbytes = psutil.virtual_memory().total / 1024 / 1024
    ram_str = locale.format("%d", ram_mbytes, grouping=True)
    verbose_print("# RAM               : %s MB" % ram_str)
    if classpath:
        verbose_print("# CLASSPATH         : %s" % os.environ.get("CLASSPATH"))
    verbose_print("# Python Version    : %s" % sys.version.replace("\n", " "))
    verbose_print("# Program Version   : %s %s" % (program_name_with_command(), __version__))
    verbose_print("")
    verbose_print("# %s %s" % (timestamp(), command_line_short()))


def print_arguments(args):
    """Print the program options.

    Parameters
    ----------
    args : argparse.Namespace
        Pre-parsed program arguments
    """
    verbose_print("Options:")
    options_dict = vars(args)
    for key in sorted(list(options_dict)):
        if key in ["subparser_name", "func", "excepthook"]:
            continue
        verbose_print("    %s=%s" % (key, options_dict[key]))
    verbose_print("")


def detect_numeric_option_in_parameters_str(parameters, option):
    """Parses a string of options to find a particular option followed by its
    value.

    This function is used to determine whether a user is overriding a default
    parameter.

    Parameters
    ----------
    parameters : str
        String of multiple user-specified options with values
    option : str
        The particular option to look for.  Should have a leading '-' if expected.

    Returns
    -------
    found : bool
        True if the option is detected, false otherwise

    Examples
    --------
    >>> detect_numeric_option_in_parameters_str("-p1", "-p")
    True
    >>> detect_numeric_option_in_parameters_str("-p 1", "-p")
    True
    >>> detect_numeric_option_in_parameters_str("-p 22", "-p")
    True
    >>> detect_numeric_option_in_parameters_str("-a x -p 22 -z 44", "-p")
    True
    >>> detect_numeric_option_in_parameters_str("-p 22", "-n")
    False
    >>> detect_numeric_option_in_parameters_str("-p", "-p")
    False
    >>> detect_numeric_option_in_parameters_str("", "-p")
    False
    """
    regex_str = option + "[ ]*([0-9])+"
    regex = re.compile(regex_str)
    return regex.search(parameters) is not None


def extract_version_str(program_name, command_line):
    """Run a program with options to emit the version and construct
    a string with the program name a version.

    Parameters
    ----------
    program_name : str
        Friendly program name -- this will be returned in the version string
    command_line : str
        Command to be executed to get the version somewhere in the output

    Returns
    -------
    version_str : str
        A version string of the form "program_name version 2.3.0" or
        "Unrecognized program_name version".
    """
    # Run the command to get the version, split and clean the output
    text = command.run(command_line)
    lines = text.split('\n')
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]

    # Look for an output line with the word "version"
    for line in lines:
        lowerline = line.lower()
        if "version" in lowerline:
            lowerline = lowerline.replace(':', ' ')
            tokens = lowerline.split()
            for index, token in enumerate(tokens):
                if token == "version" and len(tokens) > index + 1:
                    return program_name + " version " + tokens[index + 1]

    # if only one line and only one token, assume it is the version identifier
    if len(lines) == 1:
        tokens = lines[0].split()
        if len(tokens) == 1:
            return program_name + " version " + tokens[0]

    return "Unrecognized " + program_name + " version"


def is_directory_writeable(path):
    """Returns true if the specified directory is writeable.

    Parameters
    ----------
    path : str
        Directory path to create.

    Returns
    -------
    writeable : bool
        True if the path is a writeable directory
    """
    return os.path.isdir(path) and os.access(path, os.W_OK | os.X_OK)


def mkdir_p(path):
    """Python equivalent of bash mkdir -p.

    Parameters
    ----------
    path : str
        Directory path to create.

    Raises
    ------
    OSError if the directory does not already exist and cannot be created
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def remove_file(filename):
    """Remove a file without complaints if the file does not exist.

    Parameters
    ----------
    path : str
        File path to create.

    Raises
    ------
    OSError if an error occurs (other than the file not existing)
    """
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT: # No such file or directory
            raise # Re-raise exception if a different error occurred


def add_file_suffix(path, suffix, enable=True):
    """Insert a suffix at the end of a file name, but before the file extension.

    Parameters
    ----------
    path : str
        File path with a file extension.
    suffix : str
        Suffix to add to the file name.
    enable : bool, optional defaults to True
        If not True, the path is returned unchanged.

    Returns
    -------
    path : str
        Modified file path.

    Examples
    -------
    # disabled
    >>> add_file_suffix("aaa/bbb/ccc", ".suffix", enable=False)
    'aaa/bbb/ccc'

    # No extension
    >>> add_file_suffix("aaa/bbb/ccc", ".suffix", enable=True)
    'aaa/bbb/ccc.suffix'

    # Extension
    >>> add_file_suffix("aaa/bbb/ccc.txt", ".suffix", enable=True)
    'aaa/bbb/ccc.suffix.txt'
    """
    if enable:
        path_without_extension, extension = os.path.splitext(path)
        path = path_without_extension + suffix + extension
    return path


def read_properties(prop_file_path, recognize_vars=False):
    """Read a file of name=value pairs and load them into a dictionary.

    Name and value must be separated by =.
    Spaces are stripped around both name and value.
    Quotes around value are removed.
    Comment lines starting with # are ignored.

    Parameters
    ----------
    prop_file_path : str
        Path to the property file.
    recognize_vars : bool, optional
        If True, values containing dollar prefixed tokens are expanded to
        previously read property values or environment variable values if
        possible.

    Returns
    -------
    properties : dict
        Dictionary mapping property names to str values.

    Examples
    --------
    # Setup tests
    >>> from tempfile import NamedTemporaryFile

    # Create a property file
    >>> f = NamedTemporaryFile(delete=False, mode='w')
    >>> filepath = f.name
    >>> print(" # Comment", file=f)         # space before comment
    >>> print("Aaa=bbb", file=f)            # no spaces around =, no quotes around value
    >>> print(" Bbb = '2' ", file=f)        # leading space, spaces around =, single quoted value
    >>> print("Ccc = \\"33\\"", file=f)     # double quotes around value
    >>> print("Ddd=", file=f)               # no value
    >>> print("Eee='5\\"'", file=f)         # double quote within single-quoted value
    >>> print("Fff=\\"6'\\"", file=f)       # single quote within double-quoted value
    >>> print("Ggg=$Aaa", file=f)           # ignore other parameter when recognize_vars=False
    >>> f.close()

    # Read the properties
    >>> read_properties(filepath)
    OrderedDict([('Aaa', 'bbb'), ('Bbb', '2'), ('Ccc', '33'), ('Ddd', ''), ('Eee', '5"'), ('Fff', "6'"), ('Ggg', '$Aaa')])
    >>> os.unlink(filepath)

    # Create a property file
    >>> os.environ["DUMMY_ENV_VAR"] = "D ummy"
    >>> f = NamedTemporaryFile(delete=False, mode='w')
    >>> filepath = f.name
    >>> print("Aaa = bbb", file=f)                          # define Aaa property
    >>> print("Bbb = $Aaa ", file=f)                        # substitute property
    >>> print("Ccc=$DUMMY_ENV_VAR", file=f)                 # substitute environment variable
    >>> print("Ddd=$DUMMY_ENV_VAR $more", file=f)           # substitute environment variable, ignore something that cannot be resolved
    >>> print("Eee=$DUMMY_ENV_VAR$Aaa $Bbb", file=f)        # substitute two properties and one environment variable
    >>> f.close()

    # Read the properties
    >>> read_properties(filepath, recognize_vars=True)
    OrderedDict([('Aaa', 'bbb'), ('Bbb', 'bbb'), ('Ccc', 'D ummy'), ('Ddd', 'D ummy $more'), ('Eee', 'D ummybbb bbb')])
    >>> os.unlink(filepath)
    """
    properties = OrderedDict()
    with open(prop_file_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            assign_op_idx = line.find('=')
            if assign_op_idx < 0:
                continue
            key = line[:assign_op_idx]
            value = line[assign_op_idx + 1:]
            key = key.strip()
            value = value.strip()
            if value.startswith('"') and value.endswith('"'):
                value = value.strip('"')
            elif value.startswith("'") and value.endswith("'"):
                value = value.strip("'")

            # expand known tokens
            if recognize_vars:
                old_value = None
                token_regex = "[$][a-zA-Z_]+[a-zA-Z0-9_]*"
                while value != old_value:
                    old_value = value
                    match = re.search(token_regex, value)
                    if match:
                        token = match.group(0)[1:] # drop the leading $
                        if token in properties:
                            token_value = properties[token]
                            value = re.sub(token_regex, token_value, value, count=1)
                        elif token in os.environ:
                            token_value = os.environ[token]
                            value = re.sub(token_regex, token_value, value, count=1)

            properties[key] = value

    return properties


def sample_id_from_dir(dir_path):
    """Return the sample id of a sample directory as defined by the cfsan snp pipeline.

    This assumes the files for the sample are placed in a directory whose basename is the sample_id.

    Parameters
    ----------
    dir_path : str
        Path to the directory.

    Returns
    -------
    sample_id : str
        The cfsan snp pipeline defines the sample id as the basename of the parent directory.

    Examples
    --------
    >>> sample_id_from_dir("") == os.path.basename(os.getcwd())
    True

    >>> sample_id_from_dir(".") == os.path.basename(os.getcwd())
    True

    >>> sample_id_from_dir("aaa")
    'aaa'

    >>> sample_id_from_dir("aaa/bbb")
    'bbb'
    """
    sample_dir = os.path.abspath(dir_path)
    sample_id = os.path.basename(sample_dir)
    return sample_id


def sample_id_from_file(file_path):
    """Return the sample id of a file as defined by the cfsan snp pipeline.

    This assumes the file has been placed in a directory whose basename is the sample_id.

    Parameters
    ----------
    file_path : str
        Path to the file.  This could be a fastq, sam, bam. pileup, vcf, etc.

    Returns
    -------
    sample_id : str
        The cfsan snp pipeline defines the sample id as the basename of the parent directory.

    Examples
    --------
    >>> sample_id_from_file("x.fastq") == os.path.basename(os.getcwd())
    True

    >>> sample_id_from_file("aaa/x.fastq")
    'aaa'

    >>> sample_id_from_file("aaa/bbb/x.fastq")
    'bbb'
    """
    sample_abs_path = os.path.abspath(file_path)
    sample_dir = os.path.dirname(sample_abs_path)
    sample_id = os.path.basename(sample_dir)
    return sample_id


def log_error(message):
    """Write an error message to the error log if enabled.
    """
    # Log to error file if it is defined, which happens automatically if running run_snp_pipeline.sh.
    # The errorOutputFile may also be defined manually if running scripts without run_snp_pipeline.sh.
    error_output_file = os.environ.get("errorOutputFile")
    if error_output_file:
        with open(error_output_file, "a") as err_log:
            print(message, file=err_log)


def report_error(message):
    """Send an error message to the error log and to stderr.
    """
    log_error(message)

    # send the detail error message to stderr -- this will put the error
    # message in the process-specific log file.
    sys.stdout.flush() # make sure stdout is flushed before printing the error
    if message:
        print(message, file=sys.stderr)


def fatal_error(message):
    """Log a fatal error message to the error summary file and to stderr and then exit.
    """
    report_error(message)
    sys.exit(1)


def sample_warning(message):
    """Log a warning to the error summary file with special formatting and also print the warning to stderr.
    Do not exit regardless of the StopOnSampleError setting.

    Parameters
    ----------
    message : str
        Warning message text.
    """
    error_output_file = os.environ.get("errorOutputFile")
    if error_output_file:
        with open(error_output_file, "a") as err_log:
            print("%s warning:" % program_name_with_command(), file=err_log)
            print(message, file=err_log)
            print("================================================================================", file=err_log)

    # Also send the detail error message to stderr -- this will put the error message in the
    # process-specific log file.
    sys.stdout.flush() # make sure stdout is flushed before printing the error
    if message:
        print(message, file=sys.stderr)


def global_error(message):
    """
    Log a fatal error to the error summary file and exit with error code 100
    to cause Sun Grid Engine to also detect the error.

    This method always stops pipeline execution, it does not care about the
    StopOnSampleError flag.

    Args:
        message : str
            Error message
    """
    # Log the event to the error log
    error_output_file = os.environ.get("errorOutputFile")
    if error_output_file:
        with open(error_output_file, "a") as err_log:
            print("%s failed." % program_name_with_command(), file=err_log)
            if message:
                print(message, file=err_log)
            print("=" * 80, file=err_log)

    # send the detail error message to stderr -- this will put the error
    # message in the process-specific log file.
    sys.stdout.flush() # make sure stdout is flushed before printing the error
    if message:
        print(message, file=sys.stderr)

    # Exit 100 does two things:
    # 1. Sun Grid Engine will stop execution of dependent jobs
    # 2. run_snp_pipeline.sh will know this error has already been reported
    sys.exit(100)


def sample_error(message, continue_possible=False):
    """
    Log an error to the error summary file and conditionally exit with error
    code 100 to cause Sun Grid Engine to also detect the error.

    The StopOnSampleError and continue_possible flags control the
    pipeline exit / continuation behavior.  Possible behaviors are:
    - Stop this step and all subsequent steps of the pipeline if
      StopOnSampleError is true or unset
    - Stop execution of this step, but continue subsequent steps if
      StopOnSampleError is false and continue_possible is false
    - Allow this step to continue if
      StopOnSampleError is false and continue_possible is true

    Args:
        message : str
            Error message
        continue_possible : boolean
            Indicates if it is possible to continue execution.  Setting this
            flag true may allow the code to continue without exiting if
            configured to do so.
    """
    stop_on_error_env = os.environ.get("StopOnSampleError")
    stop_on_error = stop_on_error_env is None or stop_on_error_env == "true"

    # Log the event to the error log
    error_output_file = os.environ.get("errorOutputFile")
    if error_output_file:
        with open(error_output_file, "a") as err_log:
            if stop_on_error or not continue_possible:
                print("%s failed." % program_name_with_command(), file=err_log)
            else:
                print("%s" % program_name_with_command(), file=err_log)
            print(message, file=err_log)
            print("=" * 80, file=err_log)

    # send the detail error message to stderr -- this will put the error
    # message in the process-specific log file.
    sys.stdout.flush() # make sure stdout is flushed before printing the error
    print(message, file=sys.stderr)

    # Exit 100 does two things:
    # 1. Sun Grid Engine will stop execution of dependent jobs
    # 2. run_snp_pipeline.sh will know this error has already been reported
    if stop_on_error:
        # run_snp_pipeline.sh will know this error has already been reported
        sys.exit(100)
    else:
        # run_snp_pipeline.sh will know this error has already been reported,
        # but it should not stop execution
        if not continue_possible:
            sys.exit(98)


def handle_global_exception(exc_type, exc_value, exc_traceback):
    """
    This function replaces the default python unhandled exception handler.
    It Logs the error and returns error code 100 to cause Sun Grid Engine to
    also detect the error.
    """
    external_program_command = None
    if exc_type == subprocess.CalledProcessError:
        external_program_command = exc_value.cmd

    # Report the exception in the error log if configuired
    error_output_file = os.environ.get("errorOutputFile")
    if error_output_file:
        with open(error_output_file, "a") as err_log:
            trace_entries = traceback.extract_tb(exc_traceback)
            file_name, line_number, function_name, code_text = trace_entries[-1]
            exc_type_name = exc_type.__name__

            print("Error detected while running %s." % program_name_with_command(), file=err_log)
            print("", file=err_log)
            print("The command line was:", file=err_log)
            print("    %s" % command_line_short(), file=err_log)
            print("", file=err_log)

            if external_program_command:
                print("The error occured while running:", file=err_log)
                print("    %s" % external_program_command, file=err_log)
            else:
                print("%s exception in function %s at line %d in file %s" % (exc_type_name, function_name, line_number, file_name), file=err_log)
                print("    %s" % code_text, file=err_log)
            print("=" * 80, file=err_log)

    # Report the exception in the usual way to stderr
    sys.stdout.flush() # make sure stdout is flushed before printing the trace
    if external_program_command:
        print("Error occured while running:", file=sys.stderr)
        print("    %s" % external_program_command, file=sys.stderr)
    else:
        traceback.print_exception(exc_type, exc_value, exc_traceback)

    # Exit 100 does two things:
    # 1. Sun Grid Engine will stop execution of dependent jobs
    # 2. run_snp_pipeline.sh will know this error has already been reported
    sys.exit(100)


def handle_sample_exception(exc_type, exc_value, exc_traceback):
    """
    This function replaces the default python unhandled exception handler.
    It Logs the error and returns error code 100 to cause Sun Grid Engine to
    also detect the error.
    """
    external_program_command = None
    if exc_type == subprocess.CalledProcessError:
        external_program_command = exc_value.cmd

    # Report the exception in the error log if configuired
    error_output_file = os.environ.get("errorOutputFile")
    if error_output_file:
        with open(error_output_file, "a") as err_log:
            trace_entries = traceback.extract_tb(exc_traceback)
            file_name, line_number, function_name, code_text = trace_entries[-1]
            exc_type_name = exc_type.__name__

            print("Error detected while running %s." % program_name_with_command(), file=err_log)
            print("", file=err_log)
            print("The command line was:", file=err_log)
            print("    %s" % command_line_short(), file=err_log)
            print("", file=err_log)

            if external_program_command:
                print("The error occured while running:", file=err_log)
                print("    %s" % external_program_command, file=err_log)
            else:
                print("%s exception in function %s at line %d in file %s" % (exc_type_name, function_name, line_number, file_name), file=err_log)
                print("    %s" % code_text, file=err_log)
            print("=" * 80, file=err_log)

    # Report the exception in the usual way to stderr
    sys.stdout.flush() # make sure stdout is flushed before printing the trace
    if external_program_command:
        print("Error occured while running:", file=sys.stderr)
        print("    %s" % external_program_command, file=sys.stderr)
    else:
        traceback.print_exception(exc_type, exc_value, exc_traceback)

    # Exit 100 does two things:
    # 1. Sun Grid Engine will stop execution of dependent jobs
    # 2. run_snp_pipeline.sh will know this error has already been reported
    stop_on_error_env = os.environ.get("StopOnSampleError")
    stop_on_error = stop_on_error_env is None or stop_on_error_env == "true"
    if stop_on_error:
        # run_snp_pipeline.sh will know this error has already been reported
        sys.exit(100)
    else:
        # run_snp_pipeline.sh will know this error has already been reported,
        # but it should not stop execution
        sys.exit(98)


def which(executable):
    """Search the PATH for the specified executable file.

    Parameters
    ----------
    executable : str
        Name of executable

    Returns
    -------
    path : str, or None
        Path to the executable or None if not found on the PATH.
        If the excutable is available at more than one location,
        Only the first path is returned.
    """
    path = os.environ.get('PATH', "")

    for p in path.split(os.pathsep):
        p = os.path.join(p, executable)
        if os.access(p, os.X_OK): # Can the file be executed?
            return p

    return None


def verify_existing_input_files(error_prefix, file_list, error_handler=None, continue_possible=False):
    """Verify each file in a list of files exists.  It does
    not matter whether the file is empty.
    Missing files are reported in the verbose log.

    Parameters
    ----------
    error_prefix : str
        First part of error message to be logged per bad file
    file_list : list
        List of relative or absolute paths to files
    error_handler : str, optional
        Allowed values are "global", "sample", or None.
        When set to "global", any error are handled by the global_error() function which will
        cause program exit.
        When set to "sample", any error are handled by the sample_error() function which might
        cause program exit when configured to do so.
        When set to None, errors are logged with no possibility of exiting the program
    continue_possible : boolean, optional
        Only used when error_handler is "sample".  Indicates if it is possible
        to continue execution.  Setting this flag true may allow the code
        to continue without exiting if configured to do so.

    Returns
    -------
    count : int
        Number of missing files
    """
    if error_handler not in ["global", "sample", None]:
        raise ValueError("Invalid error_handler: %s" % repr(error_handler))

    err_messages = []
    for file_path in file_list:

        if not os.path.isfile(file_path):
            err_messages.append("%s %s does not exist." % (error_prefix, file_path))
            continue

    if len(err_messages) > 0:
        err_message = '\n'.join(err_messages)
        if error_handler == "global":
            global_error(err_message)
        elif error_handler == "sample":
            sample_error(err_message, continue_possible=continue_possible)
        else:
            report_error(err_message)

    return len(err_messages)


def verify_non_empty_input_files(error_prefix, file_list, error_handler=None, continue_possible=False, empty_ok=False):
    """Verify each file in a list of files exists and is non-empty.
    Missing or empty files are reported in the verbose log.

    Parameters
    ----------
    error_prefix : str
        First part of error message to be logged per bad file
    file_list : list
        List of relative or absolute paths to files
    error_handler : str, optional
        Allowed values are "global", "sample", or None.
        When set to "global", any error are handled by the global_error() function which will
        cause program exit.
        When set to "sample", any error are handled by the sample_error() function which might
        cause program exit when configured to do so.
        When set to None, errors are logged with no possibility of exiting the program
    continue_possible : boolean, optional
        Only used when error_handler is "sample".  Indicates if it is possible
        to continue execution.  Setting this flag true may allow the code
        to continue without exiting if configured to do so.
    empty_ok : bool, optional, defaults to False
        Flag to allow empty files in special cases.

    Returns
    -------
    count : int
        Number of missing or empty files
    """
    if error_handler not in ["global", "sample", None]:
        raise ValueError("Invalid error_handler: %s" % repr(error_handler))

    err_messages = []
    for file_path in file_list:

        if not os.path.isfile(file_path):
            err_messages.append("%s %s does not exist." % (error_prefix, file_path))
            continue
        if not empty_ok and os.path.getsize(file_path) == 0:
            err_messages.append("%s %s is empty." % (error_prefix, file_path))
            continue

    if len(err_messages) > 0:
        err_message = '\n'.join(err_messages)
        if error_handler == "global":
            global_error(err_message)
        elif error_handler == "sample":
            sample_error(err_message, continue_possible=continue_possible)
        else:
            report_error(err_message)

    return len(err_messages)


def verify_non_empty_directory(error_prefix, directory, error_handler=None, continue_possible=False):
    """Verify a directory exists and is non-empty.
    Missing or empty files are reported in the verbose log.

    Parameters
    ----------
    error_prefix : str
        First part of error message to be logged
    directory : list
        List of relative or absolute paths to files
    error_handler : str, optional
        Allowed values are "global", "sample", or None.
        When set to "global", any error are handled by the global_error() function which will
        cause program exit.
        When set to "sample", any error are handled by the sample_error() function which might
        cause program exit when configured to do so.
        When set to None, errors are logged with no possibility of exiting the program
    continue_possible : boolean, optional
        Only used when error_handler is "sample".  Indicates if it is possible
        to continue execution.  Setting this flag true may allow the code
        to continue without exiting if configured to do so.

    Returns
    -------
    okay : bool
        If the function does not exit the program, it returns true if directory exists and is not empty
    """
    if error_handler not in ["global", "sample", None]:
        raise ValueError("Invalid error_handler: %s" % repr(error_handler))

    error_prefix = error_prefix + ' ' if error_prefix else ''
    message = None
    if not os.path.exists(directory):
        message = error_prefix + directory + " does not exist."
    elif not os.path.isdir(directory):
        message = error_prefix + directory + " is not a directory."
    elif len(os.listdir(directory)) == 0:
        message = error_prefix + directory + " is empty."
    else:
        return True

    if error_handler == "global":
        global_error(message)
    elif error_handler == "sample":
        sample_error(message, continue_possible=continue_possible)
    else:
        report_error(message)
    return False


def global_error_on_missing_file(file_path, program):
    """Generate a global error if a specified file is missing or empty after
    running a named program.

    Parameters
    ----------
    file_path : str
        Path to the file to check
    program : str
        Name of the program that should have created the file.

    Returns
    -------
    None
        If the file is missing or empty, this function does not return, the program exits.
    """
    if not os.path.isfile(file_path):
        global_error("Error: %s does not exist after running %s." % (file_path, program))
    if os.path.getsize(file_path) == 0:
        global_error("Error: %s is empty after running %s." % (file_path, program))


def sample_error_on_missing_file(file_path, program, empty_ok=False):
    """Generate a sample error if a specified file is missing or empty after
    running a named program.

    Parameters
    ----------
    file_path : str
        Path to the file to check
    program : str
        Name of the program that should have created the file.
    empty_ok : bool, optional, defaults to False
        Flag to allow empty files in special cases.

    Returns
    -------
    None
        If the file is missing or (empty and not empty_ok), this function does not return, the program exits.
    """
    if not os.path.isfile(file_path):
        sample_error("Error: %s does not exist after running %s." % (file_path, program))
    if not empty_ok and os.path.getsize(file_path) == 0:
        sample_error("Error: %s is empty after running %s." % (file_path, program))


def sample_error_on_file_contains(file_path, text, program):
    """Generate a sample error if a specified file contains a specified string.

    Parameters
    ----------
    file_path : str
        Path to the file to check
    text : str
        The text to serach for in the file
    program : str
        Name of the program that should have created the file.

    Returns
    -------
    None
        If the file is missing or empty, this function does not return, the program exits.
    """
    with open(file_path) as f:
        for line in f:
            if text in line:
                sample_error("Error: %s contains unexpected text: '%s' after running %s." % (file_path, text, program))


def target_needs_rebuild(source_files, target_file):
    """Determine if a target file needs a fresh rebuild, i.e. the target does
    not exist or its modification time is older than any of its source files.

    Parameters
    ----------
    source_files : list of str
        Relative or absolute path to a list of files.
    target_file : str
        Relative or absolute path to target file.
    """
    if not os.path.isfile(target_file):
        return True

    if os.path.getsize(target_file) == 0:
        return True

    target_timestamp = os.stat(target_file).st_mtime

    for source_file in source_files:
        # A non-existing source file should neither force a rebuild, nor prevent a rebuild.
        # You should error-check the existence of the source files before calling this function.
        #
        # An empty source file should force a rebuild if it is newer than the target, just like
        # a regular non-empty source file.
        if not os.path.isfile(source_file):
            continue

        source_timestamp = os.stat(source_file).st_mtime
        if source_timestamp > target_timestamp:
            return True

    return False


def targets_needs_rebuild(source_files, target_files):
    """Determine if target files need a fresh rebuild, i.e. any one of target files does
    not exist or its modification time is older than any of its source files.

    Args:
        source_files : relative or absolute path to a list of files
        target_files : relative or absolute path to a list of file
    """

    if len(target_files) == 0:
        return True
    else:
        if not os.path.isfile(target_files[0]):
            return True

        oldest_timestamp = os.stat(target_files[0]).st_mtime

        if os.path.getsize(target_files) == 0:
            return True

        for target_file in target_files:
            if (not os.path.isfile(target_file)) or (os.path.getsize(target_file) == 0):
                return True

            target_timestamp = os.stat(target_file).st_mtime
            if oldest_timestamp > target_timestamp:
                oldest_timestamp = target_timestamp

        for source_file in source_files:
            # A non-existing source file should neither force a rebuild, nor prevent a rebuild.
            # You should error-check the existence of the source files before calling this function.
            #
            # An empty source file should force a rebuild if it is newer than the target, just like
            # a regular non-empty source file.
            if not os.path.isfile(source_file):
                continue

            source_timestamp = os.stat(source_file).st_mtime
            if source_timestamp > oldest_timestamp:
                return True

    return False


def write_list_of_snps(file_path, snp_dict):
    """Write out list of snps for all samples to a single file.

    Args:
        file_path : path to snplist file to be written
        snp_dict  : dictionary with key = tuple(CHROM, POS) -> value = list[sampleName1, sampleName2, ..., sampleNameN]

    Returns:
        Nothing
    """

    with open(file_path, "w") as snp_list_file_object:
        for key in sorted(snp_dict.keys()):
            sample_list = snp_dict[key]
            snp_list_file_object.write("%s\t%d\t%d\t%s\n" % (key[0], key[1], len(sample_list), "\t".join(sample_list)))


def read_snp_position_list(snp_list_file_path):
    """Read list of snp positions across all samples from the snplist.txt.

    Args:
        snp_list_file_path : path to snplist file to be written

    Returns:
        snp_list  : sorted list of tuple(str(CHROM), int(POS))
    """

    snp_list = list()
    with open(snp_list_file_path, "r") as snp_list_file_object:
        for line in snp_list_file_object:
            chrom, pos = line.split()[0:2]
            snp_list.append((chrom, int(pos)))
    return snp_list


def write_reference_snp_file(reference_file_path, snp_list_file_path,
                             snp_reference_file_path):
    """Write out the snp fasta file for the reference.fasta using the snp
    position file ( snplist.txt).
    """
    #TODO finish documentation
    #TODO actual code is more general than stated. Fix this.

    with open(snp_list_file_path, "r") as snp_list_file:
        position_list = [line.split()[0:2] for line in snp_list_file]
    match_dict = SeqIO.to_dict(SeqIO.parse(reference_file_path, "fasta"))

    with open(snp_reference_file_path, "w") as snp_reference_file_object:
        for ordered_id in sorted(match_dict.keys()):
            ref_str = ""
            for chrom_id, pos in position_list:
                if chrom_id == ordered_id:
                    ref_str += match_dict[ordered_id][int(pos) - 1].upper()
            record = SeqRecord(Seq(ref_str), id=ordered_id, description="")
            SeqIO.write([record], snp_reference_file_object, "fasta")


def convert_vcf_file_to_snp_set(vcf_file_path):
    """convert vcf files to a set of SNPs.

    Args:
        vcf_file_path : relative or absolute path to the sample VCF file

    Returns:
        snp_set  : set of (CHROM, POS) tuples

    """

    snp_set = set()

    with open(vcf_file_path, 'r') as vcf_file_object:
        vcf_reader = vcf.Reader(vcf_file_object)
        for vcf_data_line in vcf_reader:
            key = (vcf_data_line.CHROM, vcf_data_line.POS)
            snp_set.add(key)

    return snp_set


def calculate_sequence_distance(seq1, seq2, case_insensitive=True):
    """Calulate the number of nucleotide differences between two sequences.

    The sequences must be the same length.

    Args:
        seq1 : DNA string 1
        seq2 : DNA string 2
        case_insensitive : optional flag for case insensitive compare, defaults to True

    Returns:
        int number of differences
    """
    if case_insensitive:
        allowed_bases = frozenset(['A', 'C', 'G', 'T'])
        seq1 = seq1.upper()
        seq2 = seq2.upper()
    else:
        allowed_bases = frozenset(['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'])

    mismatches = 0
    for pos in range(len(seq1)):
        base1 = seq1[pos]
        base2 = seq2[pos]
        if base1 not in allowed_bases:
            continue
        if base2 not in allowed_bases:
            continue
        if base1 != base2:
            mismatches += 1
    return mismatches


def merge_regions(regions):
    """Coalesce regions.

    Scans a sorted list of region starting and ending positions looking
    for the outer-most start and end positions to coalesce overlapping
    and contained regions into a smaller list of larger regions.

    Parameters
    ----------
    regions : list of tuples
        List of (start, end) position integers.

    Returns
    -------
    regions : list of tuples
        List of merged (start, end) position integers.

    Examples
    --------
    >>> # Empty list
    >>> merge_regions([])
    []

    >>> # Only one region
    >>> merge_regions([(10,20)])
    [(10, 20)]

    >>> # Discard contained region at left
    >>> merge_regions([(10,20), (10,15)])
    [(10, 20)]

    >>> # Discard contained region at right
    >>> merge_regions([(10,20), (15,20)])
    [(10, 20)]

    >>> # Discard contained region exact match
    >>> merge_regions([(10,20), (10,20)])
    [(10, 20)]

    >>> # Discard contained region fully contained
    >>> merge_regions([(10,20), (11,19)])
    [(10, 20)]

    >>> # Extend region by overlap right
    >>> merge_regions([(10,20), (15,25)])
    [(10, 25)]

    >>> # Extend region by overlap left
    >>> merge_regions([(10,20), (5,15)])
    [(5, 20)]

    >>> # Extend immediately adjacent region by extension
    >>> merge_regions([(10,20), (21,30)])
    [(10, 30)]

    >>> # No overlap
    >>> merge_regions([(40,50), (25,30)])
    [(25, 30), (40, 50)]

    >>> # Single position region : discard contained region
    >>> merge_regions([(40,50), (40,40)])
    [(40, 50)]

    >>> # Single position region : discard contained region
    >>> merge_regions([(40,50), (50,50)])
    [(40, 50)]

    >>> # Single position region : discard contained region
    >>> merge_regions([(40,50), (41,41)])
    [(40, 50)]

    >>> # Single position region : discard contained region
    >>> merge_regions([(40,50), (49,49)])
    [(40, 50)]

    >>> # Single position region : extend immediately adjacent region by extension
    >>> merge_regions([(10,10), (11,21)])
    [(10, 21)]

    >>> # Single position region : extend immediately adjacent region by extension
    >>> merge_regions([(10,20), (21,21)])
    [(10, 21)]

    >>> # Single position region : merge two immediately adjacent single-position regions
    >>> merge_regions([(20,20), (21,21)])
    [(20, 21)]

    >>> # Single position region : no overlap
    >>> merge_regions([(40,50), (60,60)])
    [(40, 50), (60, 60)]

    >>> # Single position region : no overlap
    >>> merge_regions([(40,40), (50,60)])
    [(40, 40), (50, 60)]

    >>> # Single position region : no overlap
    >>> merge_regions([(40,40), (50,50)])
    [(40, 40), (50, 50)]
    """
    if len(regions) == 0:
        return regions
    regions = sorted(regions)
    merged_regions = list()
    merged_regions.append(regions[0])
    for region in regions[1:]:
        last_merged_region = merged_regions[-1]
        last_merged_region_start, last_merged_region_end = last_merged_region
        region_start, region_end = region
        if region_start >= last_merged_region_start and region_end <= last_merged_region_end:
            pass # discard region contained in the last region
        elif region_start <= (last_merged_region_end + 1) and region_end > last_merged_region_end:
            merged_regions[-1] = (last_merged_region_start, region_end) # extend last region by overlapping or adjacent region
        else:
            merged_regions.append(region) # add non-overlapping region to sorted list
    return merged_regions


def in_region(pos, regions):
    """Find whether a position is included in a region.

    Parameters
    ----------
    pos : int
        DNA base position.
    regions : list of tuples
        List of (start, end) position integers.

    Returns
    -------
    bool
        True if the position is within an of the regions, False otherwise.

    Examples
    --------
    # Empty list
    >>> in_region(1, [])
    False

    # In list
    >>> in_region(10, [(3, 5), (9, 12)])
    True

    # Not in list
    >>> in_region(10, [(3, 5), (11, 12)])
    False
    """
    for region in regions:
        if (pos >= region[0]) and (pos <= region[1]):
            return True

    return False


def compute_num_processes_and_threads(max_cpu_cores, threads_per_process):
    """Compute the number of allowed processes and threads given the maximum allowed
    number of CPUs and requested number of threads per process.

    Parameters
    ----------
    max_cpu_cores : int or None
        The maximum allowed number of CPU cores to consume by all instances of the process.
        If set the None, it implies no limit to the number of CPU cores that can be used.
    threads_per_process : int
        The user-requested number of threads to use.

    Returns
    -------
    num_processes : int or None
        The computed maximum number of allowed concurrent processes, or None to allow unlimited.
    threads_per_process : int
        The number of threads per process instance which might be less than requested if the
        number of allowed cpus is less than the requested number of threads.

    Examples
    --------
    # max CPU not set
    >>> compute_num_processes_and_threads(None, 8)
    (None, 8)

    # max CPU set, not multiple of threads
    >>> compute_num_processes_and_threads(20, 8)
    (2, 8)

    # max CPU set, multiple of threads
    >>> compute_num_processes_and_threads(24, 8)
    (3, 8)

    # max CPU less than desired threads
    >>> compute_num_processes_and_threads(2, 8)
    (1, 2)
    """
    if max_cpu_cores is None:
        num_processes = None
    elif max_cpu_cores >= threads_per_process:
        num_processes = int(max_cpu_cores / threads_per_process)
    else:
        num_processes = 1
        threads_per_process = max_cpu_cores
    return num_processes, threads_per_process


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
    threads_option : str, or list of str
        The exact spelling of a command line option to set the number of threads, for example "-n".
        This parameter can also be a list of command line options when there is more than one way to
        specify the number of threads, for example, ["-nt", "--num_threads"].
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

    # Two different ways to configure the number of threads
    >>> os.environ["RealignerTargetCreator_ExtraParams"] = "-nt 7"
    >>> configure_process_threads("RealignerTargetCreator_ExtraParams", ["-nt", "--num_threads"], 8, 24)
    (3, 7)
    >>> os.environ["RealignerTargetCreator_ExtraParams"] = "--num_threads 7"
    >>> configure_process_threads("RealignerTargetCreator_ExtraParams", ["-nt", "--num_threads"], 8, 24)
    (3, 7)
    """
    threads_options = [threads_option] if type(threads_option) is str else threads_option
    match = False
    for threads_option in threads_options:
        regex_str = threads_option + "[ \t]*([0-9]+)"
        extra_params = os.environ.get(extra_params_env_var, "")
        match = re.search(regex_str, extra_params)
        if match:
            break # stop looking for thread control options in the command line as soon as we find the first one
    if match:
        configured_threads_per_process = int(match.group(1))
        threads_per_process = configured_threads_per_process
    else:
        threads_per_process = default_threads_per_process

    max_processes, threads_per_process = compute_num_processes_and_threads(max_cpu_cores, threads_per_process)

    threads_option += ' ' + str(threads_per_process)
    if match and threads_per_process != configured_threads_per_process:
        extra_params = re.sub(regex_str, threads_option, extra_params)
        os.environ[extra_params_env_var] = extra_params
    elif not match:
        if extra_params:
            extra_params += ' '
        os.environ[extra_params_env_var] = extra_params + threads_option

    return max_processes, threads_per_process


def find_path_in_path_list(search_item, env_var, case_sensitive=False):
    """Search a colon-separated environment variable for a specified string.
    Return the path containing the string.

    Parameters
    ----------
    search_item : str
        The string to search for in the path.  This string should have enough
        characters to uniquely identify the path, but it does not have to be the
        entire path to the file.  Any unique substring will work.  Typically, the
        basename of the file is enough.
    env_var : str
        Name of the environment variable containing the path with the search item.
    case_sensitive : bool, optional, defaults to False
        When false, the search will be case insensitive, so any combination of
        uppercase/lowercase search item is fine.

    Returns
    -------
    path : str
        The path to the jar file as specified in the environment variable if found.
        Returns None if not found.

    Examples
    --------
    >>> os.environ["TESTPATH"] = "/a/a/picard.jar:/b/b/gatk.jar:/c/c/VarScan.jar"
    >>> find_path_in_path_list("missing", "TESTPATH") is None
    True
    >>> find_path_in_path_list("Picard", "TESTPATH") # first, case insensitive
    '/a/a/picard.jar'
    >>> find_path_in_path_list("Gatk", "TESTPATH") # middle, case insensitive
    '/b/b/gatk.jar'
    >>> find_path_in_path_list("varscan", "TESTPATH") # last, case insensitive
    '/c/c/VarScan.jar'
    >>> find_path_in_path_list("picard", "TESTPATH", case_sensitive=True) # first, case sensitive
    '/a/a/picard.jar'
    >>> find_path_in_path_list("gatk", "TESTPATH", case_sensitive=True) # middle, case sensitive
    '/b/b/gatk.jar'
    >>> find_path_in_path_list("VarScan", "TESTPATH", case_sensitive=True) # last, case sensitive
    '/c/c/VarScan.jar'
    """
    paths = os.environ.get(env_var, "")
    if not case_sensitive:
        search_item = search_item.lower()

    paths = paths.split(':')
    for path in paths:
        path_to_search = path.lower() if not case_sensitive else path
        if search_item in path_to_search:
            return path

    return None

