#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,    Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
"""
The Fortran logging module terminates on error. This cannot be tested by the
unit testing framework as it terminates the unit tests as well.
"""

import datetime
import re

from testframework import LFRicLoggingTest, MpiTest, TestEngine, TestFailed


###############################################################################
class LogModErrorSerialTest(MpiTest):  # pylint: disable=too-few-public-methods
    """
    Tests logging in serial scenarios.
    """

    def __init__(self):
        super().__init__(processes=1)

        self._minimum_timestamp = datetime.datetime.utcnow()

    def test(self, returncode: int, out: str, err: str):
        """
        Tests that loggin an error ends execution.
        """
        expected_level = "ERROR"
        expected_message = " An error was logged."

        if returncode == 0:  # pylint: disable=no-else-raise
            message = "Logging an error did not cause termination to end"
            raise TestFailed(message)
        elif returncode == 127:
            raise TestFailed("Test executable not found")
        elif returncode > 128:
            raise TestFailed("Execution fault such as segmentation fault")

        if out != "":
            message = (
                "Expected no output on standard out:\n"
                + f"Standard out: {out}"
            )
            raise TestFailed(message)

        try:
            timestamp_string, level, report = err.split(":", 2)
            timestamp_without_timezone = timestamp_string[:-5]

            timestamp = datetime.datetime.strptime(
                timestamp_without_timezone, "%Y%m%d%H%M%S.%f"
            )
        except Exception as ex:
            message = f"Unable to get timestamp from message: {err}"
            raise TestFailed(message) from ex

        if timestamp < self._minimum_timestamp:
            message = (
                f"Expected a timestamp after {self._minimum_timestamp}"
                f" but read {timestamp}"
            )
            raise TestFailed(message)

        if level != expected_level:
            message = 'Expected "{}" but read "{}"'
            raise TestFailed(message.format(expected_level, level))

        # We only check the first line as compilers tend to print the return
        # code as well. This will remain true until we can use Fortran 2008 and
        # "stop error".
        #
        first, _, _ = report.partition("\n")
        if first != expected_message:
            message = 'Expected "{}" but read "{}"'
            raise TestFailed(message.format(expected_message, first))

        message = "Logging an error caused exit as expected with code {code}"
        return message.format(code=returncode)


##############################################################################
class LogModErrorParallelTest(LFRicLoggingTest):
    """
    Tests logging in MPI parallel scenarios.
    """

    # pylint: disable=too-few-public-methods
    def __init__(self):
        super().__init__(processes=2)

        self._minimum_timestamp = datetime.datetime.now(datetime.timezone.utc)
        line_pattern_string = (
            r"(\d{4})(\d\d)(\d\d)(\d\d)(\d\d)(\d\d)"
            r"\.(\d{3})([+-])(\d\d)(\d\d):P(\d+):\s*(\w+)"
            r":\s+(.+)"
        )
        self._line_pattern = re.compile(line_pattern_string)

    def test(self, returncode: int, out: str, err: str):
        """
        Tests that logging an error terminates execution when run in parallel.
        """
        # pylint: disable=too-many-locals
        expected_level = "ERROR"
        expected_message = "An error was logged."

        if returncode == 0:
            message = "Logging an error did not cause termination to end"
            raise TestFailed(message)

        if out != "":
            message = (
                "Expected no output on standard out:\n"
                + f"Standard out: {out}"
            )
            raise TestFailed(message)

        # ToDo: Ideally we would test for stderr output here but whether some
        # is generated or not is dependent on whether the compiler in use
        # supports generating a backtrace for warnings. Rather than trying to
        # solve that problem we'll leave it as an exercise for later.
        #
        # In light of this we will just ignore the "err" argument we are
        # passed as part of the LFRicLoggingTest interface.

        pet_log = self.getLFRicLoggingLog()
        pet_log = "\n".join(pet_log.splitlines())

        match = self._line_pattern.match(pet_log)
        if match:
            try:
                tzsign = -1 if match.group(8) == "-" else 1
                tzhours = int(match.group(9))
                tzmins = int(match.group(10))
                timezone = datetime.timezone(
                    tzsign * datetime.timedelta(hours=tzhours, minutes=tzmins)
                )
                timestamp = datetime.datetime(
                    int(match.group(1)),  # Year
                    int(match.group(2)),  # Month
                    int(match.group(3)),  # Day
                    int(match.group(4)),  # Hour
                    int(match.group(5)),  # Minute
                    int(match.group(6)),  # Second
                    int(match.group(7)) * 1000,  # Microseconds
                    timezone,  # Timezone
                )
            except Exception as ex:
                message = f"Bad timestamp format: {pet_log}"
                raise TestFailed(message) from ex
            process = int(match.group(11))
            level = match.group(12)
            report = match.group(13)
        else:
            message = f"Unexpected log message: {pet_log}"
            raise TestFailed(message)

        if timestamp < self._minimum_timestamp:
            message = (
                f"Expected a timestamp after {self._minimum_timestamp}"
                f" but read {timestamp}"
            )
            raise TestFailed(message)

        if process < 0:
            message = "Process number went negative"
            raise TestFailed(message)

        if level != expected_level:
            message = "Expected '{}' but read '{}'"
            raise TestFailed(message.format(expected_level, level))

        # We only check the first line as compilers tend to print the return
        # code as well. This will remain true until we can use Fortran 2008 and
        # "stop error".
        #
        first, _, _ = report.partition("\n")
        if first != expected_message:
            message = 'Expected "{}" but read "{}"'
            raise TestFailed(message.format(expected_message, first))

        message = "Logging an error caused exit as expected with code {code}"
        return message.format(code=returncode)


##############################################################################
if __name__ == "__main__":
    TestEngine.run(LogModErrorSerialTest())
    TestEngine.run(LogModErrorParallelTest())
