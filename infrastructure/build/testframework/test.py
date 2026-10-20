#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
from abc import ABC, abstractmethod
import collections
import errno
import math
import os
from pathlib import Path
from subprocess import run
import sys
from tempfile import NamedTemporaryFile
from textwrap import dedent
import time
from typing import List


##############################################################################
class AbstractTest(ABC):
    """
    Base functionality of a test. This class is responsible for actually
    running a test and handing the result on to the handler method.
    """

    def __init__(self, executable):
        """
        Constructor.

        parameter executable - Executable command in list form.
        """
        self._executable = executable

    @abstractmethod
    def test(self, returncode, out, err):
        """
        Examines the result of running the test for correctness.

        parameter returncode - OS level result of running the executable.
        parameter out        - String holding standard out from executable.
        parameter err        - String holding standard error from executable.

        Beware, liable to memory exhaustion in the face of large amounts of
        output.

        Throw TestFailed object on finding a mistake.
        """
        pass

    def performTest(self):
        """
        Runs the executable and passes results to handler method.
        """
        result = run([str(arg) for arg in self._executable],
                     capture_output=True, encoding="utf8")
        self.post_execution(result.returncode)
        return self.test(
            result.returncode,
            self.filterOut(result.stdout),
            self.filterErr(result.stderr)
        )

    def post_execution(self, code):
        """
        Perform any management tasks after the executable has been run.

        The default implementation does nothing. Only override if you have
        something to do.
        """
        pass

    def filterOut(self, out):
        """
        Processes the standard output string.

        The default implementation simply returns the string unprocessed.
        Only override if you need to filter.
        """
        return out

    def filterErr(self, err):
        """
        Processes the standard error string.

        The default implementation simply returns the string unprocessed.
        Only override if you need to filter.
        """
        return err


##############################################################################
class Test(AbstractTest):
    """
    Base for serial tests.
    """

    def __init__(self, command=sys.argv[1]):
        if type(command) is not list:
            command = [command]
        super().__init__(command)


##############################################################################
class MpiTest(AbstractTest):
    """
    Base for parallel tests.
    """

    _mpiexec_broken = None

    @staticmethod
    def set_mpiexec_broken():
        MpiTest._mpiexec_broken = True

    def __init__(self, command=sys.argv[1], processes=4):
        self._processes = processes

        if type(command) is not list:
            command = [command]

        commandString = " ".join(command)
        commandName = os.path.basename(command[0])
        self._startTag = "Start {name}".format(name=commandName)
        self._doneTag = "Done {name}".format(name=commandName)

        # If environment variable doesn't exist, variable will be set to
        # None. This will cause the argument to gain its default value,
        # which is the correct behaviour.
        #
        tmp_dir = os.environ.get("SCRATCH")
        with NamedTemporaryFile(mode="wt", delete=False, dir=tmp_dir,
                                prefix="run-", suffix=".sh") as launch_script:
            launch_script.write(
                dedent(
                    f"""
                    #!/bin/sh
                    echo {self._startTag}
                    echo {self._startTag} >&2
                    {commandString}
                    result=$?
                    echo {self._doneTag}
                    echo {self._doneTag} >&2
                    sync
                    exit $result
                    """
                ).lstrip()
            )

            self._script_file = Path(launch_script.name)
            self._script_file.chmod(0o700)

        if MpiTest._mpiexec_broken:
            mpi_launcher = "mpirun"
        else:
            mpi_launcher = "mpiexec"

        mpiCommand = [
            mpi_launcher,
            "-n",
            str(self._processes),
            self._script_file,
        ]
        super().__init__(mpiCommand)

    def __del__(self):
        self._script_file.unlink()

    def __rejectWaffle(self, log: str) -> List[str]:
        """
        Strips out MPI cruft from a log file.
        """
        newLog = []
        state = "spinup"
        processesRunning = 0
        for line in log.splitlines():
            if state == "spinup":
                if line == self._startTag:
                    processesRunning += 1
                    if processesRunning == self._processes:
                        state = "spindown"
            elif state == "spindown":
                if line == self._doneTag:
                    processesRunning -= 1
                    if processesRunning == 0:
                        state = "done"
                else:  # line does not start with 'Done '
                    newLog.append(line)

        # If no markers were found, it probably means something catastrophic
        # so we return the log unchanged.
        #
        if state == "spinup":
            return log
        else:
            return "\n".join(newLog)

    def filterOut(self, out):
        """
        Strips MPI cruft from standard out.
        """
        return self.__rejectWaffle(out)

    def filterErr(self, err):
        """
        Strip MPI cruft from standard out.
        """
        return self.__rejectWaffle(err)


##############################################################################
class LFRicLoggingTest(MpiTest):
    """
    Base for LFRicLogging parallel tests.
    """

    def __init__(
        self, command=sys.argv[1], name="log_mod_error_test.Log", processes=4
    ):
        """
        Constructor.

        Arguments:

          name - String - Name given to LFRicLogging and therefore the one
          which appears in log file names.
        """
        super().__init__(command, processes)
        self._application_name = name
        self._LFRicLoggingLog = collections.defaultdict(lambda: None)

    def getLFRicLoggingLog(self, process=0):
        """
        Returns the log file generated by LFRicLogging on a particular MPI
        process.

        Arguments:

          process - Integer - The process for which a log file is desired.
                              Defaults to zero.
        """
        return self._LFRicLoggingLog[process]

    def performTest(self):
        """
        Removes any old log files and runs the executable.
        """
        # Remove any existing log files
        for filename in os.listdir("."):
            if filename.startswith("PET"):
                try:
                    os.remove(filename)
                except IOError as ex:
                    # If the file doesn't exist continue otherwise raise error
                    if ex.errno is errno.ENOENT:
                        print(f"{ex}, continuing")
                    else:
                        raise

        return super().performTest()

    def post_execution(self, return_code):
        """
        Caches log files for future retrieval.

        When running on a single processor all log messages are written to
        stdout and no log file is generated. Therefore, only attempt to cache
        the log file if the number processes is greater than one.

        Some filing systems (Lustre for instance) suffer a lag between a file
        being closed and it becoming visible for reading. To try and mitigate
        this a number of attempts will be made to open the file with a sleep
        period between them. Only after these attempts are exhausted will
        the file be reported missing.
        """
        if self._processes > 1:
            for number in range(0, self._processes):
                width = int(math.floor(math.log10(self._processes))) + 1
                filenameFormat = "PET{{number:0{width}d}}.{{name}}"
                filename = filenameFormat.format(width=width).format(
                    number=number, name=self._application_name
                )
                number_retries = 10
                delay_seconds = 60
                for attempt in range(1, number_retries + 1):
                    try:
                        with open(filename, "rt") as handle:
                            self._LFRicLoggingLog[number] = handle.read()
                    except IOError:
                        if attempt < number_retries:
                            message = (
                                '{timestamp}: File "{filename}" not found, '
                                + "waiting {delay} seconds..."
                            )
                            print(
                                message.format(
                                    timestamp=time.clock(),
                                    filename=filename,
                                    delay=delay_seconds,
                                ),
                                file=sys.stderr,
                            )
                            time.sleep(delay_seconds)
                        else:
                            raise
