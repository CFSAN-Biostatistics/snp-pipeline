"""
This module provides an abstraction layer to run jobs on high performance computers
using torque, grid, or locally with xargs.
"""

from __future__ import print_function
from __future__ import absolute_import

import os
import psutil
import subprocess
import sys


class JobRunnerException(Exception):
    """Raised for fatal JobRunner errors"""

class JobRunner(object):
    def __init__(self, hpc_type, strip_job_array_suffix=True):
        """Initialize an hpc job runner object.

        Parameters
        ----------
        hpc_type : str
            Type of job runner.  Possible values are "grid", "torque", and "local".
        strip_job_array_suffix : bool, optional defaults to True
            When true, the dot and array suffix in the job id is removed before returning the job id.

        Examples
        --------
        >>> runner = JobRunner("foobar")
        Traceback (most recent call last):
        ValueError: hpc_type must be one of: "grid", "torque", "local"
        """
        hpc_type = hpc_type.lower()
        if hpc_type not in ["grid", "torque", "local"]:
            raise ValueError('hpc_type must be one of: "grid", "torque", "local"')

        self.hpc_type = hpc_type
        self.strip_job_array_suffix = strip_job_array_suffix

        if hpc_type == 'grid':
            self.subtask_env_var_name = "SGE_TASK_ID"
        elif hpc_type == 'torque':
            self.subtask_env_var_name = "PBS_ARRAYID"


    def _make_qsub_command(self, job_name, log_file, wait_for=[], wait_for_array=[], slot_dependency=False, threads=1, parallel_environment=None, num_tasks=None, max_processes=None):
        """Create the qsub command line to run a job on a computing cluster.

        Parameters
        ----------
        job_name : str
            Job name that will appear in the job scheduler queue.
        log_file : str
            Path to the combined stdout / stderr log file.
        wait_for : str or list of str, optional defaults to empty list
            Single job id or list of jobs ids to wait for before beginning execution.
        wait_for_array : str or list of str, optional defaults to empty list
            Single array job id or list of array jobs ids to wait for before beginning execution.
        slot_dependency : bool, optional defaults to False
            Ignored for all schedulers but grid engine.
            If true, the sub-tasks of the array job being submitted will be dependent on the
            completion of the corresponding sub-tasks of the jobs in the wait_for_array.  Has
            no effect on the dependencies of non-array jobs.
        threads : int, optional defaults to 1
            Number of CPU threads consumed by the job.
        parallel_environment : str, optional defaults to None
            Name of the grid engine parallel execution environment.  This must be specified when
            consuming more than one thread on grid engine.  Ununsed for any other job scheduler.
        num_tasks : int, optional defaults to None
            When specified, the job becomes an array job with num_tasks sub-tasks.
        max_processes : int, optional defaults to None
            If not None, it sets the maximium number of concurrent processes for an array job.


        Returns
        -------
        qsub_command_line : str
            Qsub command line for Grid or torque.

        Examples
        --------
        >>> runner = JobRunner("local")
        >>> runner._make_qsub_command("JobName", "log")
        Traceback (most recent call last):
        ValueError: _make_qsub_command() does not support hpc type local

        # grid
        >>> runner = JobRunner("grid")
        >>> os.environ["GridEngine_QsubExtraParams"] = "-q service.q"
        >>> runner._make_qsub_command("JobName", "log", "777", "888", threads=8, parallel_environment="mpi")
        'qsub -terse -V -j y -cwd -N JobName -o log -hold_jid 777,888 -pe mpi 8 -q service.q'

        >>> runner._make_qsub_command("JobName", "log", ["666", "777"], ["888", "999"], threads=8, parallel_environment="mpi")
        'qsub -terse -V -j y -cwd -N JobName -o log -hold_jid 666,777,888,999 -pe mpi 8 -q service.q'

        >>> runner._make_qsub_command("JobName", "log", ["666", "777"], ["888", "999"], slot_dependency=True, threads=8, parallel_environment="mpi")
        'qsub -terse -V -j y -cwd -N JobName -o log -hold_jid 666,777 -hold_jid_ad 888,999 -pe mpi 8 -q service.q'

        # array job
        >>> runner._make_qsub_command("JobName", "log", ["666", "777"], ["888", "999"], slot_dependency=True, threads=8, parallel_environment="mpi", num_tasks=99, max_processes=2)
        'qsub -terse -t 1-99 -V -j y -cwd -N JobName -o log -hold_jid 666,777 -hold_jid_ad 888,999 -tc 2 -pe mpi 8 -q service.q'

        # torque
        >>> runner = JobRunner("torque")
        >>> os.environ["Torque_QsubExtraParams"] = "-q short.q"
        >>> cmd = runner._make_qsub_command("JobName", "log", "777", "888", threads=8)
        >>> cmd == "qsub -V -j oe -d %s -N JobName -o log -W depend=afterok:777,afterokarray:888 -l nodes=1:ppn=8 -q short.q" % os.getcwd()
        True

        >>> cmd = runner._make_qsub_command("JobName", "log", ["666", "777"], ["888", "999"], threads=8)
        >>> cmd == "qsub -V -j oe -d %s -N JobName -o log -W depend=afterok:666:777,afterokarray:888:999 -l nodes=1:ppn=8 -q short.q" % os.getcwd()
        True

        # array job
        >>> cmd = runner._make_qsub_command("JobName", "log", ["666", "777"], ["888", "999"], slot_dependency=True, threads=8, num_tasks=44, max_processes=2)
        >>> cmd == "qsub -t 1-44%%2 -V -j oe -d %s -N JobName -o log -W depend=afterok:666:777,afterokarray:888:999 -l nodes=1:ppn=8 -q short.q" % os.getcwd()
        True
        """
        if self.hpc_type not in ["grid", "torque"]:
            raise ValueError("_make_qsub_command() does not support hpc type %s" % self.hpc_type)

        if self.hpc_type == "grid":
            array_option = " -t 1-" + str(num_tasks) if num_tasks else ""
            qsub_command_line = "qsub -terse" + array_option + " -V -j y -cwd -N " + job_name + " -o " + log_file

            if isinstance(wait_for, str):
                wait_for = [wait_for]
            if isinstance(wait_for_array, str):
                wait_for_array = [wait_for_array]
            if not slot_dependency:
                wait_for.extend(wait_for_array) # combine lists
            if len(wait_for) > 0:
                qsub_command_line += " -hold_jid " + ','.join(wait_for)
            if slot_dependency and len(wait_for_array) > 0:
                qsub_command_line += " -hold_jid_ad " + ','.join(wait_for_array)

            if max_processes:
                qsub_command_line += " -tc " + str(max_processes)
            if threads > 1:
                if not parallel_environment:
                    raise ValueError("You must use a parallel environment when consuming more than one thread on grid engine")
                qsub_command_line += " -pe " + parallel_environment + ' ' + str(threads)

            qsub_extra_params = os.environ.get("GridEngine_QsubExtraParams") or ""
            qsub_command_line += ' ' + qsub_extra_params
            return qsub_command_line

        if self.hpc_type == "torque":
            max_processes_option = "%%%i" % max_processes if max_processes else ""
            array_option = " -t 1-" + str(num_tasks) + max_processes_option if num_tasks else ""
            qsub_command_line = "qsub" + array_option + " -V -j oe -d " + os.getcwd() + " -N " + job_name + " -o " + log_file

            if isinstance(wait_for, str):
                wait_for = [wait_for]
            if isinstance(wait_for_array, str):
                wait_for_array = [wait_for_array]
            dependencies = []
            if len(wait_for) > 0:
                dependencies.append("afterok:" + ':'.join(wait_for))
            if len(wait_for_array) > 0:
                dependencies.append("afterokarray:" + ':'.join(wait_for_array))
            if len(dependencies) > 0:
                qsub_command_line += " -W depend=" + ','.join(dependencies)

            if threads > 1:
                qsub_command_line += " -l nodes=1:ppn=" + str(threads)

            qsub_extra_params = os.environ.get("Torque_QsubExtraParams") or ""
            qsub_command_line += ' ' + qsub_extra_params
            return qsub_command_line


    def run(self, command_line, job_name, log_file, wait_for=[], wait_for_array=[], threads=1, parallel_environment=None):
        """Run a non-array job.  Stderr is redirected (joined) to stdout.

        Parameters
        ----------
        command_line : str
            Command with all arguments to be executed.
        job_name : str
            Job name that will appear in the job scheduler queue.
        log_file : str
            Path to the combined stdout / stderr log file.
        wait_for : str or list of str, optional defaults to empty list
            Single job id or list of jobs ids to wait for before beginning execution.
            Ignored when running locally.
        wait_for_array : str or list of str, optional defaults to empty list
            Single array job id or list of array jobs ids to wait for before beginning execution.
            Ignored when running locally.
        threads : int, optional defaults to 1
            Number of CPU threads consumed by the job, unused when running locally.
        parallel_environment : str, optional defaults to None
            Name of the grid engine parallel execution environment.  This must be specified when
            consuming more than one thread on grid engine.  Ununsed for any other job scheduler.

        Returns
        -------
        job_id : str
            Grid or torque job id.  Returns '0' in local mode.

        Raises
        ------
        CalledProcessError on non-zero exitcode

        Examples
        --------
        # Normal case - verify job id is '0', stdout and stderr written to log file
        >>> from tempfile import NamedTemporaryFile
        >>> fout = NamedTemporaryFile(delete=False, mode='w'); fout.close()
        >>> runner = JobRunner("local")
        >>> # Parenthesis are needed when the command line contains multiple commands separated by semicolon
        >>> job_id = runner.run("(echo text to stdout; echo text to stderr 1>&2)", "JobName", fout.name)
        >>> type(job_id) == type("this is a string")
        True
        >>> job_id
        '0'
        >>> f = open(fout.name); out = f.read(); f.close(); os.unlink(fout.name)
        >>> print(out.strip())
        text to stdout
        text to stderr

        # Error case, external program returns non-zero
        >>> job_id = runner.run("exit 100", "JobName", "")
        Traceback (most recent call last):
        CalledProcessError: Command 'set -o pipefail; exit 100 2>&1 | tee ' returned non-zero exit status 100
        """
        if self.hpc_type == "local":
            command_line = "set -o pipefail; " + command_line + " 2>&1 | tee " + log_file

            # flush stdout to keep the unbuffered stderr in chronological order with stdout
            sys.stdout.flush()

            # Run command. Wait for command to complete. If the return code was zero then return, otherwise raise CalledProcessError
            subprocess.check_call(command_line, shell=True, executable="bash")
            return '0'

        else: # grid or torque
            qsub_command_line = self._make_qsub_command(job_name, log_file, wait_for, wait_for_array, threads, parallel_environment)

            # Run command and return its stdout output as a byte string.
            # If the return code was non-zero it raises a CalledProcessError.
            shell_command_line = "echo " + command_line + " | " + qsub_command_line
            job_id = subprocess.check_output(shell_command_line, shell=True)
            return job_id.strip()


    def run_array(self, command_line, job_name, log_file, array_file, num_tasks=None, max_processes=None, wait_for=[], wait_for_array=[], slot_dependency=False, threads=1, parallel_environment=None):
        """Run an array of sub-tasks with the work of each task defined by a single line in the specified array_file.

        Parameters
        ----------
        command_line : str
            Command to be executed with parameter placeholders of the form {1}, {2}, {3} ...
        job_name : str
            Job name that will appear in the job scheduler queue.
        log_file : str
            Path to the combined stdout / stderr log file.  The sub-task number will be automatically appended.
        array_file : str
            Name of the file containing the arguments for each sub-task with one line
            per sub-task.  The arguments for each sub-task are found at the line number
            corresponding to the sub-task number.  The line is parsed and substituted
            into the command, replacing the parameter placeholders with the actual arguments.
        num_tasks : int, optional defaults to None
            Defines the number of subtasks in the job array.  If not specified, the array_file must exist
            and the number of tasks will be equal to the number of lines in the file.  Use this option
            when the array_file does not pre-exist and is created by a process that has not run yet.
        max_processes : int, optional defaults to None
            If None, the number of concurrent processes is limited to available CPU on an HPC
            and limited to the number of CPU cores when run locally.
            If not None, it sets the maximium number of concurrent processes for the array job.
            This works locally with xargs, and with grid and torque.
        wait_for : str or list of str, optional defaults to empty list
            Single job id or list of jobs ids to wait for before beginning execution.
            Ignored when running locally.
        wait_for_array : str or list of str, optional defaults to empty list
            Single array job id or list of array jobs ids to wait for before beginning execution.
            Ignored when running locally.
        slot_dependency : bool, optional defaults to False
            Ignored for all schedulers but grid engine.
            If true, the sub-tasks of the array job being submitted will be dependent on the
            completion of the corresponding sub-tasks of the jobs in the wait_for_array.  Has
            no effect on the dependencies of non-array jobs.
        threads : int, optional defaults to 1
            Number of CPU threads consumed by each sub-task of the job, unused when running locally.
        parallel_environment : str, optional defaults to None
            Name of the grid engine parallel execution environment.
            Ununsed for any other job scheduler.
        """
        # Determine the number of array slots
        if not num_tasks:
            if not os.path.isfile(array_file):
                raise JobRunnerException("The file %s does not exist.\nCannot start array job %s." % (array_file, job_name))
            if os.path.getsize(array_file) == 0:
                raise JobRunnerException("The file %s is empty.\nCannot start array job %s." % (array_file, job_name))

            with open(array_file) as f:
                num_tasks = sum(1 for line in f)

        if self.hpc_type == "grid":
            log_file += "-\$TASK_ID" # TODO: verify this works properly

        if self.hpc_type == "local":
            # Change parameter placeholder into bash variables ready to feed to bash through xargs
            for param_num in range(1, 10):
                placeholder = '{' + str(param_num) + '}'
                command_line = command_line.replace(placeholder, '$' + str(param_num))

            # Use all CPU cores, if no limit requested
            if max_processes is None:
                max_processes = psutil.cpu_count()

            # Number the tasks with nl to get the task number into the log file suffix.
            # Allow up to 9 parameters per command.
            command_line = "nl " + array_file + " | xargs -P " + str(max_processes) + " -n 9 -L 1 bash -c + 'set -o pipefail; " + command_line + " 2>&1 | tee " + log_file + "-$0'"

            # Run command. Wait for command to complete
            subprocess.check_call(command_line, shell=True, executable="bash") # If the return code is non-zero it raises a CalledProcessError
            return '0'

        else: # grid or torque
            qsub_command_line = self._make_qsub_command(job_name, log_file, wait_for, wait_for_array, slot_dependency, threads, parallel_environment, num_tasks, max_processes)

            shell_command_line = "qarrayrun " + self.subtask_env_var_name + ' ' + array_file + ' ' + command_line + " | " + qsub_command_line
            job_id = subprocess.check_output(shell_command_line, shell=True) # If the return code is non-zero it raises a CalledProcessError
            if self.strip_job_array_suffix:
                dot_idx = job_id.find('.')
                if dot_idx > 0:
                    job_id = job_id[0: dot_idx]
            return job_id
