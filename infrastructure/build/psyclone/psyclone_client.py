#!/usr/bin/env python3
# -----------------------------------------------------------------------------
#  (C) Crown copyright 2026 Met Office. All rights reserved.
#  The file LICENCE, distributed with this code, contains details of the terms
#  under which the code may be used.
# -----------------------------------------------------------------------------
"""
Thin client that submits a single PSyclone job to the persistent
:mod:`psyclone_server` and blocks for the result.

It is a drop-in replacement for the ``psyclone`` command as invoked by
``psyclone_psykal.mk``: everything after the program name is forwarded verbatim
as PSyclone's argument list. The client

* auto-starts the server (once) under a lockfile if it is not running;
* writes the job payload to a per-job directory and sends a *tiny* enqueue
  message (the job-directory path) over the shared request FIFO;
* waits on a per-job response FIFO, replays PSyclone's captured stdout/stderr
  and exits with PSyclone's return code.

If anything goes wrong (the server cannot be started, a timeout occurs, or the
protocol is violated) the client transparently falls back to executing the
real ``psyclone`` binary directly so that a build is never broken by the
optimisation.
"""

import json
import os
import subprocess
import sys
import tempfile
import time


# Environment variable names shared with the server.
ENV_SERVER_DIR = "PSYCLONE_SERVER_DIR"
ENV_WORKERS = "PSYCLONE_WORKERS"
ENV_DISABLE = "PSYCLONE_SERVER_DISABLE"

# How long to wait for the server to become ready, and for a job result.
SERVER_START_TIMEOUT = 60.0
JOB_TIMEOUT = 900.0


def _server_dir():
    """
    Return the directory used for server coordination files, creating it if
    required. Defaults to ``<WORKING_DIR>/.psyclone-server`` and falls back to
    a location under the system temporary directory.

    :rtype: str
    """
    base = os.environ.get(ENV_SERVER_DIR)
    if not base:
        working = os.environ.get("WORKING_DIR", "working")
        base = os.path.join(working, ".psyclone-server")
    os.makedirs(base, exist_ok=True)
    return base


def _process_alive(pid):
    """Return True if a process with the given pid exists."""
    try:
        os.kill(pid, 0)
    except (OSError, ProcessLookupError):
        return False
    return True


def _server_ready(server_dir):
    """Return True if a healthy server appears to be running."""
    ready_file = os.path.join(server_dir, "server.ready")
    pid_file = os.path.join(server_dir, "server.pid")
    if not (os.path.exists(ready_file) and os.path.exists(pid_file)):
        return False
    try:
        with open(pid_file, encoding="utf8") as handle:
            pid = int(handle.read().strip())
    except (OSError, ValueError):
        return False
    return _process_alive(pid)


def _reap_stale(server_dir):
    """Remove coordination files left behind by a crashed server."""
    for name in ("server.ready", "server.pid", "request.fifo"):
        try:
            os.remove(os.path.join(server_dir, name))
        except OSError:
            pass


def _start_server(server_dir):
    """
    Ensure a server is running, starting one under a lockfile if necessary.

    :param str server_dir: server coordination directory.
    :returns: True if a ready server is available, False otherwise.
    :rtype: bool
    """
    import fcntl

    if _server_ready(server_dir):
        return True

    lock_path = os.path.join(server_dir, "server.lock")
    with open(lock_path, "w", encoding="utf8") as lock_file:
        fcntl.flock(lock_file, fcntl.LOCK_EX)
        # Re-check now that we hold the lock: another make process may have
        # started the server while we were waiting.
        if _server_ready(server_dir):
            return True

        # A pid file without a live process means a previous server crashed.
        _reap_stale(server_dir)

        env = dict(os.environ)
        env[ENV_SERVER_DIR] = server_dir

        pid_file = os.path.join(server_dir, "server.pid")
        # Launch the server fully detached so it outlives this make recipe.
        server_module = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "psyclone_server.py")
        process = subprocess.Popen(
            [sys.executable, server_module],
            env=env, stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            start_new_session=True)
        with open(pid_file, "w", encoding="utf8") as handle:
            handle.write(str(process.pid))

        deadline = time.monotonic() + SERVER_START_TIMEOUT
        while time.monotonic() < deadline:
            if os.path.exists(os.path.join(server_dir, "server.ready")):
                return True
            if process.poll() is not None:
                # Server died during start-up.
                _reap_stale(server_dir)
                return False
            time.sleep(0.05)
        return False


def _submit(server_dir, argv):
    """
    Submit a job to the running server and return its result.

    :param str server_dir: server coordination directory.
    :param list argv: PSyclone command-line arguments.
    :returns: (returncode, stdout, stderr).
    :rtype: tuple[int, str, str]
    :raises RuntimeError: if the job cannot be completed by the server.
    """
    job_dir = tempfile.mkdtemp(prefix="job-", dir=server_dir)
    response_fifo = os.path.join(job_dir, "response")
    os.mkfifo(response_fifo)

    request = {
        "argv": argv,
        "cwd": os.getcwd(),
        # PYTHONPATH entries are needed so the server can import optimisation
        # recipes and their psyclone_tools helper (mirrors the makefile's
        # PYTHONPATH=$(LFRIC_BUILD)/psyclone:$$PYTHONPATH).
        "sys_path": [p for p in
                     os.environ.get("PYTHONPATH", "").split(os.pathsep) if p],
    }
    with open(os.path.join(job_dir, "request.json"), "w",
              encoding="utf8") as handle:
        json.dump(request, handle)

    # Tiny, atomic enqueue message: just the job directory path.
    message = (job_dir + "\n").encode("utf8")
    request_fifo = os.path.join(server_dir, "request.fifo")
    fifo_fd = os.open(request_fifo, os.O_WRONLY)
    try:
        os.write(fifo_fd, message)
    finally:
        os.close(fifo_fd)

    # Block until a worker opens the response FIFO and writes the result.
    result = _read_response(response_fifo)
    _cleanup(job_dir)
    return result["returncode"], result["stdout"], result["stderr"]


def _read_response(response_fifo):
    """Read and decode the JSON result from a per-job response FIFO."""
    deadline = time.monotonic() + JOB_TIMEOUT
    # Opening a FIFO for reading blocks until a writer appears; guard the whole
    # read with an alarm so a dead worker cannot hang the build forever.
    import signal

    def _timeout(_signum, _frame):
        raise RuntimeError("timed out waiting for PSyclone server response")

    previous = signal.signal(signal.SIGALRM, _timeout)
    try:
        remaining = max(1.0, deadline - time.monotonic())
        signal.setitimer(signal.ITIMER_REAL, remaining)
        with open(response_fifo, encoding="utf8") as handle:
            payload = handle.read()
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)

    if not payload:
        raise RuntimeError("empty response from PSyclone server")
    return json.loads(payload)


def _cleanup(job_dir):
    """Remove a completed job directory."""
    try:
        for name in os.listdir(job_dir):
            os.remove(os.path.join(job_dir, name))
        os.rmdir(job_dir)
    except OSError:
        pass


def _fallback(argv):
    """Run the real psyclone binary directly, returning its exit code."""
    return subprocess.call(["psyclone"] + argv)


def main(argv=None):
    """
    Entry point.

    :param list argv: PSyclone arguments (defaults to ``sys.argv[1:]``).
    :returns: the process exit code.
    :rtype: int
    """
    if argv is None:
        argv = sys.argv[1:]

    # Allow the server to be disabled entirely (e.g. for debugging).
    if os.environ.get(ENV_DISABLE):
        return _fallback(argv)

    try:
        server_dir = _server_dir()
        if not _start_server(server_dir):
            return _fallback(argv)
        returncode, out, err = _submit(server_dir, argv)
    except Exception:  # pylint: disable=broad-except
        # Any failure in the fast path must not break the build.
        return _fallback(argv)

    if out:
        sys.stdout.write(out)
    if err:
        sys.stderr.write(err)
    return returncode


if __name__ == "__main__":
    sys.exit(main())

