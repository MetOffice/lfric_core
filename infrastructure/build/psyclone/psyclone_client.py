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

Server lifetime is pinned to the build: the *owning* make process is resolved
here (see :func:`_resolve_owner`) and handed to the server, which exits as soon
as that process disappears. The coordination directory is keyed on the same
process, so one server - and only one - serves a whole top-level ``make``,
including all of its recursive sub-makes, and is cleaned up when make finishes.

If anything goes wrong (the server cannot be started, a timeout occurs, or the
protocol is violated) the client transparently falls back to executing the
real ``psyclone`` binary directly so that a build is never broken by the
optimisation.
"""

import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import time

try:
    import psyclone_procs
except ImportError:  # pragma: no cover - executed from an unusual sys.path
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import psyclone_procs


# Environment variable names shared with the server.
ENV_SERVER_DIR = "PSYCLONE_SERVER_DIR"
ENV_WORKERS = "PSYCLONE_WORKERS"
ENV_DISABLE = "PSYCLONE_SERVER_DISABLE"
ENV_OWNER_PID = "PSYCLONE_OWNER_PID"
ENV_OWNER_START = "PSYCLONE_OWNER_START"
ENV_TRANSIENT = "PSYCLONE_SERVER_DIR_TRANSIENT"

# Prefix of the auto-created (per build) server directories.
SERVER_DIR_PREFIX = "psyclone-server"

# How long to wait for the server to become ready, and for a job result.
SERVER_START_TIMEOUT = 60.0
JOB_TIMEOUT = 900.0


def _resolve_owner():
    """
    Identify the build process the server's lifetime should be pinned to.

    Preference order:

    1. ``PSYCLONE_OWNER_PID`` if it names a live process. The LFRic build sets
       this from ``lfric.mk`` to the pid of the *top-level* make, so every
       recursive sub-make agrees on the same owner.
    2. the outermost make process in this client's ancestry - a robust
       fallback when the makefiles have not been updated, or when PSyclone is
       driven by some other make-based build.
    3. the process group leader, then this client's parent, for invocations
       that have nothing to do with make.

    Watching make itself matters: under a non-interactive shell (rose-stem and
    cylc job scripts) make is *not* the process group leader, so a server that
    watched the process group would outlive a manually killed build.

    :returns: (pid, start_time) of the owner, or (None, "") if none was found.
    :rtype: tuple[int or None, str]
    """
    candidate = os.environ.get(ENV_OWNER_PID)
    if candidate:
        try:
            pid = int(candidate)
        except ValueError:
            pid = None
        if pid and psyclone_procs.process_alive(pid):
            started = (os.environ.get(ENV_OWNER_START)
                       or psyclone_procs.start_time(pid))
            return pid, started

    pid = psyclone_procs.outermost_make_pid()
    if pid is None:
        for fallback in (os.getpgrp(), os.getppid()):
            if fallback > 1 and psyclone_procs.process_alive(fallback):
                pid = fallback
                break
    if pid is None:
        return None, ""
    return pid, psyclone_procs.start_time(pid)


def _usable_dir(path):
    """
    Create ``path`` (private to this user) and confirm it is safe to use.

    The auto-generated location is predictable, so refuse to use anything that
    is not a real directory belonging to us.

    :param str path: candidate directory.
    :rtype: bool
    """
    try:
        os.makedirs(path, mode=0o700, exist_ok=True)
        info = os.lstat(path)
    except OSError:
        return False
    if not stat.S_ISDIR(info.st_mode):
        return False
    return info.st_uid == os.getuid()


def _server_dir(owner_pid):
    """
    Return the directory used for server coordination files, creating it if
    required.

    ``PSYCLONE_SERVER_DIR`` wins if set. Otherwise the directory is keyed on
    the owning make process and placed under the system temporary directory:
    one server per build (rather than one per ``WORKING_DIR``), on a local
    filesystem where FIFOs are reliable, and removed when the build ends.
    ``<WORKING_DIR>/.psyclone-server`` remains as a last resort.

    :param owner_pid: pid of the owning build process.
    :type owner_pid: int or None
    :returns: (directory, transient) where ``transient`` marks a directory the
        server should delete when it shuts down.
    :rtype: tuple[str, bool]
    :raises RuntimeError: if no usable directory could be created.
    """
    explicit = os.environ.get(ENV_SERVER_DIR)
    if explicit:
        if not _usable_dir(explicit):
            raise RuntimeError(f"unusable server directory: {explicit}")
        return explicit, False

    if owner_pid:
        shared = os.path.join(tempfile.gettempdir(),
                              f"{SERVER_DIR_PREFIX}-{os.getuid()}-{owner_pid}")
        if _usable_dir(shared):
            return shared, True

    working = os.environ.get("WORKING_DIR", "working")
    fallback = os.path.join(working, ".psyclone-server")
    if _usable_dir(fallback):
        return fallback, False
    raise RuntimeError("no usable PSyclone server directory")


def _sweep_stale_dirs(keep):
    """
    Remove server directories belonging to builds that have finished.

    A server normally deletes its own directory, so this only tidies up after
    one that was killed outright (SIGKILL) along with its build.

    :param str keep: the directory in use by this client.
    """
    root = tempfile.gettempdir()
    marker = f"{SERVER_DIR_PREFIX}-{os.getuid()}-"
    try:
        names = os.listdir(root)
    except OSError:
        return
    keep = os.path.realpath(keep)
    for name in names:
        if not name.startswith(marker):
            continue
        path = os.path.join(root, name)
        if os.path.realpath(path) == keep:
            continue
        try:
            owner = int(name[len(marker):])
        except ValueError:
            continue
        if psyclone_procs.process_alive(owner):
            continue
        try:
            if os.lstat(path).st_uid != os.getuid():
                continue
        except OSError:
            continue
        if _process_alive(_read_pid(os.path.join(path, "server.pid"))):
            continue
        shutil.rmtree(path, ignore_errors=True)


def _process_alive(pid):
    """Return True if a process with the given pid exists."""
    return psyclone_procs.process_alive(pid)


def _read_pid(path):
    """
    Read a pid from a file.

    :param str path: file containing a decimal pid.
    :returns: the pid, or None if it could not be read.
    :rtype: int or None
    """
    try:
        with open(path, encoding="utf8") as handle:
            return int(handle.read().strip())
    except (OSError, ValueError):
        return None


def _server_ready(server_dir):
    """Return True if a healthy server appears to be running."""
    ready_file = os.path.join(server_dir, "server.ready")
    pid_file = os.path.join(server_dir, "server.pid")
    if not (os.path.exists(ready_file) and os.path.exists(pid_file)):
        return False
    return _process_alive(_read_pid(pid_file))


def _reap_stale(server_dir):
    """Remove coordination files left behind by a crashed server."""
    for name in ("server.ready", "server.pid", "request.fifo"):
        try:
            os.remove(os.path.join(server_dir, name))
        except OSError:
            pass


def _start_server(server_dir, owner, transient):
    """
    Ensure a server is running, starting one under a lockfile if necessary.

    :param str server_dir: server coordination directory.
    :param tuple owner: (pid, start_time) of the owning build process.
    :param bool transient: True if the server should remove ``server_dir``
        when it shuts down.
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
        if transient:
            # Tidy up after builds that were killed before their server could
            # remove its own directory.
            _sweep_stale_dirs(server_dir)

        env = dict(os.environ)
        env[ENV_SERVER_DIR] = server_dir
        env[ENV_TRANSIENT] = "1" if transient else "0"
        # Tie the server's lifetime to the owning make process (not to the
        # process group, which under a non-interactive shell is led by the
        # calling script and so outlives make). The server polls this pid and
        # exits as soon as the build finishes or is killed, by SIGINT, SIGTERM
        # or SIGKILL, rather than lingering until its idle timeout. The start
        # time guards against the pid being reused by an unrelated process.
        owner_pid, owner_start = owner
        if owner_pid:
            env[ENV_OWNER_PID] = str(owner_pid)
            env[ENV_OWNER_START] = owner_start or ""
        else:
            env.pop(ENV_OWNER_PID, None)
            env.pop(ENV_OWNER_START, None)

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
        owner = _resolve_owner()
        server_dir, transient = _server_dir(owner[0])
        if not _start_server(server_dir, owner, transient):
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

