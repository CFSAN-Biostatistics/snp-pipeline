"""
This module contains functions to run external commands/programs.
"""

from __future__ import print_function
from __future__ import absolute_import
import subprocess
import sys


def run(command, outfile=None):
    """Run a specified command as a separate process.

    Stdout can either be written to a file or captured and returned in a string.
    The stderr file handle is inherited from the calling Python program.

    Parameters
    ----------
    command : str
        Process to run with command line options and arguments.
    outfile : str, optional
        Path to stdout output file.  If not specified, the stdout is captured and returned in a string.
        You can also specify sys.stdout to inherit the stdout from the calling Python program.

    Raises
    ------
    CalledProcessError on non-zero exitcode

    Returns
    -------
    CalledProcessError exception is raised on error
    Stdout str is returned on sucess if outfile is not specified.
    Nothing is returned on success if outfile is specified.

    Examples
    --------
    # Verify Verify stdout written to a file and stderr written to stderr file handle
    >>> import os
    >>> from tempfile import NamedTemporaryFile
    >>> fout = NamedTemporaryFile(delete=False, mode='w'); fout.close()
    >>> ferr = NamedTemporaryFile(delete=False, mode='w'); ferr.close()
    >>> run("echo Hello written to file; echo Error written to file 2>%s 1>&2" % ferr.name, fout.name)
    >>> f = open(fout.name); out = f.read(); f.close(); os.unlink(fout.name)
    >>> type(out) == type("aaa")
    True
    >>> print(out.strip())
    Hello written to file
    >>> f = open(ferr.name); err = f.read(); f.close(); os.unlink(ferr.name)
    >>> print(err.strip())
    Error written to file

    # Verify stdout captured in a string and stderr written to stderr file handle
    >>> ferr = NamedTemporaryFile(delete=False, mode='w'); ferr.close()
    >>> out = run("echo Hello returned; echo Error written to file 2>%s 1>&2" % ferr.name)
    >>> type(out) == type("aaa")
    True
    >>> print(out.strip())
    Hello returned
    >>> f = open(ferr.name); err = f.read(); f.close(); os.unlink(ferr.name)
    >>> print(err.strip())
    Error written to file
    """
    sys.stdout.flush() # flush stdout to keep the unbuffered stderr in chronological order with stdout

    if outfile is None:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = proc.communicate()
        exitcode = proc.returncode
        if sys.version_info > (3,):
            stdout = stdout.decode("utf-8")  # Python 3 stdout is bytes, not str
        return stdout
    elif outfile == sys.stdout:
        subprocess.check_call(command, shell=True)
    else:
        with open(outfile, "wb") as out:
            subprocess.check_call(command, stdout=out, shell=True)
