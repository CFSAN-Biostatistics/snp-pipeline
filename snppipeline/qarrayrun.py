#!/usr/bin/env python

"""
This script executes a single slot of an array job on an HPC compute node.
It is intended to be used with Sun Grid Engine or Torque job schedulers.
It assumes every instance of the job array runs the same command but with
different arguments.  This script performs the work of looking-up the
arguments in a text file and substituting those arguments into the command
to be executed.

Parameters
----------
The script takes 3 arguments:

1. Name of the environment variable that contains the sub-task number.
   You should use SGE_TASK_ID for grid engine.
   You should use PBS_ARRAYID for torque.

2. Name of the file containing the arguments for each sub-task with one line
   per sub-task.  This script will extract the arguments for this sub-task
   at the line number identified by the environment variable above.  The
   line is parsed and substituted into the command, replacing the parameter placeholders
   with the actual arguments.

3. The remainder of the command line is the command to be executed with parameter
   placeholders of the form {1}, {2}, {3} ...

Examples
--------
# Sort some txt files, writing the sorted output to new files
ls *.txt > files.txt
echo 'qarrayrun.py SGE_TASK_ID files.txt sort -o sorted.{1} {1}' | qsub -t 1-$(cat files.txt | wc -l) -cwd -j y -V -o log

# Your input file might have multiple columns, use {2} for the 2nd column
# Sort the largest files first
ls *.txt | xargs -n 1 wc -c | sort -n -r > files.txt
echo 'qarrayrun.py SGE_TASK_ID files.txt sort -o sorted.{2} {2}' | qsub -t 1-$(cat files.txt | wc -l) -cwd -j y -V -o log

# Use the --shell option and quote your pipeline when you need shell redirection
# Remove blanks before sorting files
ls *.txt > files.txt
echo 'qarrayrun.py --shell SGE_TASK_ID files.txt "cat {1} | tr -d [:blank:] | sort > sorted.{1}"' | qsub -t 1-$(cat files.txt | wc -l) -cwd -j y -V -o log

"""

from __future__ import print_function

import argparse
import os
import re
import shlex
import subprocess
import sys
import textwrap


def get_file_line(file_name, line_num):
    """Returns the line at the specified line number from the specified file.

    Parameters
    ----------
    file_name : str
        Path to the file
    line_num : int
        Line number to extract from the file

    Returns
    -------
    line : str
        Line at line_num in file or None if the line_num is larger than the number of lines in the file

    Raises
    ------
    ValueError if line_num is less than 1
    """
    if line_num <= 0:
        raise ValueError("line_num must be greater than zero")

    line = ""
    line_counter = 0
    with open(file_name) as f:
        for line in f:
            line_counter += 1
            if line_counter == line_num:
                return line
    return None


def substitute_arguments(command_line, arguments):
    """Replace the parameter placeholders in a command line with actual arguments.

    Parameters
    ----------
    command_line : str
        Command line with parameter placeholders like {1}, {2}, {3}
    arguments : list of str
        List of arguments numbered 0,1,2 ....
        arguments[0] corresponds to placeholder {1}.

    Returns
    -------
    command_line : str
        Command line with actual arguments ready for execution

    Examples
    --------
    >>> substitute_arguments("cmd {0}/{1}/{2} -- {3}{4}", ["aa", "bb", "cc"])
    'cmd /aa/bb -- cc'
    """
    args = [""] # put an empty string at index 0
    args.extend(arguments)

    # Get a list of all the parameter numbers appearing in the command line
    param_nums = re.findall("{([0-9]+)}", command_line)
    param_nums = [int(param_num) for param_num in param_nums]

    # Replace the parameters with actual arguments
    for param_num in param_nums:
        placeholder = "{%s}" % param_num
        if param_num == 0 or param_num >= len(args):
            command_line = command_line.replace(placeholder, "")
        else:
            command_line = command_line.replace(placeholder, args[param_num])

    return command_line


def run(argv):
    """Run the qarrayrun program with the passed command line arguments.

    Parameters
    ----------
    argv : list of str
        List of command line arguments.  Usually sys.argv[1:].

    Returns
    -------
    The return code of the executed command is passed to sys.exit()

    Examples
    --------
    # Setup tests
    >>> os.environ["SGE_TASK_ID"] = "2"
    >>> from tempfile import NamedTemporaryFile
    >>> f = NamedTemporaryFile(delete=False, mode='w')
    >>> arg_file = f.name
    >>> print("A B C", file=f)
    >>> print("Argument1 Argument2 Argument3", file=f)
    >>> f.close()
    >>> f = NamedTemporaryFile(delete=False, mode='w')
    >>> out_file = f.name

    # Write the arguments in reverse order to out_file
    >>> cmd = "python -c '"
    >>> cmd += 'f = open("%s", "w");' % out_file
    >>> cmd += 'f.write("{3} {2} {1}"); f.close()'
    >>> cmd += "'"
    >>> run(["SGE_TASK_ID", arg_file, cmd])
    0

    # Read the file just created to verify reverse order
    >>> f = open(out_file)
    >>> s = f.read(); f.close();
    >>> s
    'Argument3 Argument2 Argument1'

    # Verify non-zero exit code is returned
    >>> run(["--shell", "SGE_TASK_ID", arg_file, "exit 100"])
    100

    # Clean up temp files
    >>> os.unlink(arg_file)
    >>> os.unlink(out_file)
    """
    description = textwrap.dedent("""
        Executes a single slot of an array job on an HPC computational node. This is
        intended to be used with Sun Grid Engine or Torque job schedulers when every
        instance of the job array runs the same command but with different arguments.
        This script performs the work of looking-up the arguments in a text file and
        substituting those arguments into the command to be executed.""")

    epilog = textwrap.dedent("""
        Examples
        --------
        # Sort some txt files, writing the sorted output to new files
        ls *.txt > files.txt
        echo 'qarrayrun SGE_TASK_ID files.txt sort -o sorted.{1} {1}' | qsub -t 1-$(cat files.txt | wc -l) -cwd -j y -V -o log

        # Your input file might have multiple columns, use {2} for the 2nd column
        # Sort the largest files first
        ls *.txt | xargs -n 1 wc -c | sort -n -r > files.txt
        echo 'qarrayrun SGE_TASK_ID files.txt sort -o sorted.{2} {2}' | qsub -t 1-$(cat files.txt | wc -l) -cwd -j y -V -o log

        # Use the --shell option and quote your pipeline when you need shell redirection
        # Remove blanks before sorting files
        ls *.txt > files.txt
        echo 'qarrayrun --shell SGE_TASK_ID files.txt "cat {1} | tr -d [:blank:] | sort > sorted.{1}"' | qsub -t 1-$(cat files.txt | wc -l) -cwd -j y -V -o log
        """)

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=formatter_class)

    parser.add_argument(dest="subtask_var", type=str, metavar="NAME", help="""Name of the environment variable that contains the sub-task number.
                                                                              You should use SGE_TASK_ID with Grid Engine and PBS_ARRAYID with Torque.""")
    parser.add_argument(dest="array_file", type=str, help="""Name of the file containing the arguments for each sub-task with one line
                                                             per sub-task.  This script will extract the arguments for this sub-task
                                                             at the line number identified by the sub-task environment variable
                                                             (SGE_TASK_ID or PBS_ARRAYID).  The line is parsed and substituted
                                                             into the command, replacing the parameter placeholders with the actual
                                                             arguments.""")
    parser.add_argument(dest="command", help="""The remainder of the command line is the command to be executed with parameter
                                                placeholders of the form {1}, {2}, {3} ...""")
    parser.add_argument("--shell", dest="shell", action="store_true", help="Run the command through the shell.")

    args, remainder = parser.parse_known_args(argv)
    command = [args.command] + remainder

    # Which sub-task number am I?
    subtask_num = os.environ.get(args.subtask_var)
    if not subtask_num:
        print("Error: the %s environment variable is not defined." % args.subtask_var, file=sys.stderr)
        exit(1)
    subtask_num = int(subtask_num)

    # Read and parse the substitution arguments from the input file
    line = get_file_line(args.array_file, subtask_num)
    if not line:
        exit(0) # Silently ignore attempts to process beyond the end of the file
    arguments = line.split()

    # Build the command with substituted arguments
    command_line = ' '.join(command)
    command_line = substitute_arguments(command_line, arguments)

    # Execute the command
    if args.shell:
        return_code = subprocess.call(command_line, shell=True)
    else:
        command_split = shlex.split(command_line)
        return_code = subprocess.call(command_split, shell=False)
    return return_code


def main():
    """This is the main function which is magically turned into an executable
    qarrayrun script by the setuptools entry_points.  See setup.py.

    To run this function as a script, first install the package:
        $ python setup.py develop
        or
        $ pip install --user snp-pipeline


    Parameters
    ----------
    This function must not take any parameters

    Returns
    -------
    The return code of the executed command is passed to sys.exit()
    """
    return_code = run(sys.argv[1:])
    sys.exit(return_code)


if __name__ == "__main__":
    main()
