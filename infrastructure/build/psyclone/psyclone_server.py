# -----------------------------------------------------------------------------
#  (C) Crown copyright 2026 Met Office. All rights reserved.
#  The file LICENCE, distributed with this code, contains details of the terms
#  under which the code may be used.
# -----------------------------------------------------------------------------
"""
Persistent PSyclone server for the LFRic build system.

The server imports PSyclone once and pre-forks a pool of worker processes
(sized to the build parallelism, ``-j``). Each worker keeps PSyclone resident
in memory and executes transformation jobs *in process* by calling
:func:`psyclone.generator.main`, avoiding the repeated cost of loading the
Python interpreter and PSyclone libraries from disk for every algorithm file.

Jobs are submitted by :mod:`psyclone_client` over a shared POSIX FIFO. To keep
submission atomic (a write of no more than ``PIPE_BUF`` bytes is guaranteed not
to be interleaved) the FIFO only ever carries a *tiny enqueue message*: the
path of a per-job directory. The full request payload (argv, cwd, environment)
is written to a regular file in that directory and the result is returned to
the client through a per-job response FIFO.

Protocol (per job)
------------------
* client writes ``<job-dir>\\n`` to the shared request FIFO (atomic);
* dispatcher reads the line and places the job on a shared work queue;
* a worker reads ``<job-dir>/request.json``, runs PSyclone and writes the
  result to the ``<job-dir>/response`` FIFO as ``result.json`` content;
* client reads the response FIFO, replays the captured output and exits with
  the returned code.

The server shuts itself down after ``PSYCLONE_SERVER_IDLE_TIMEOUT`` seconds
without work, so no explicit stop step is required from the build.

State isolation
---------------
Between jobs a worker resets the PSyclone global state that would otherwise
leak from one algorithm file to the next:

* ``Config._instance`` so that a per-job ``--config`` file is honoured;
* ``LFRicConstants.HAS_BEEN_INITIALISED`` so the constants are rebuilt against
  the fresh configuration;
* any modules imported from *user* locations (optimisation recipes loaded via
  ``-s`` and their ``psyclone_tools`` helpers) are purged from
  :data:`sys.modules`, because PSyclone imports recipes by bare basename and
  different components ship different ``global.py``/``<alg>.py`` files that
  would otherwise collide in the module cache.
"""

import json
import os
import signal
import sys
import time
import traceback
from multiprocessing import Process, Queue
from queue import Empty


# Sentinel placed on the work queue to ask a worker to exit.
_SHUTDOWN = None

# Number of seconds of inactivity after which the server exits.
DEFAULT_IDLE_TIMEOUT = 300.0


def _installed_prefixes():
    """
    Return the set of realpath prefixes that are considered "installed", i.e.
    modules imported from these locations are safe to cache across jobs.

    :returns: realpath prefixes for the standard library and site-packages.
    :rtype: tuple[str, ...]
    """
    import psyclone

    prefixes = {
        os.path.realpath(sys.prefix),
        os.path.realpath(sys.base_prefix),
        # site-packages directory containing the psyclone package.
        os.path.realpath(os.path.dirname(os.path.dirname(psyclone.__file__))),
    }
    return tuple(prefixes)


class _Worker:
    """
    A single pre-forked worker. Holds PSyclone resident and executes jobs
    pulled from the shared queue, resetting global state between them.

    :param queue: shared queue delivering job-directory paths.
    :type queue: multiprocessing.Queue
    """

    def __init__(self, queue):
        self._queue = queue
        # Import PSyclone up front so the cost is paid once per worker.
        import psyclone.generator  # noqa: F401  (imported for side effects)

        self._installed_prefixes = _installed_prefixes()
        # Snapshot the module cache once PSyclone is fully imported so we can
        # detect (and purge) anything imported by individual jobs.
        self._baseline_modules = frozenset(sys.modules)

    # -- state isolation ----------------------------------------------------
    def _reset_state(self):
        """Reset PSyclone global state so the next job starts clean."""
        # Reset the configuration singleton so a per-job --config is re-read.
        try:
            from psyclone.configuration import Config
            Config._instance = None
        except Exception:  # pylint: disable=broad-except
            pass

        # Force the LFRic constants to be rebuilt against the fresh config.
        try:
            from psyclone.domain.lfric import LFRicConstants
            LFRicConstants.HAS_BEEN_INITIALISED = False
        except Exception:  # pylint: disable=broad-except
            pass

        # Purge modules imported from user locations (optimisation recipes and
        # their same-directory helpers). PSyclone loads recipes by bare
        # basename, so different components would otherwise collide here.
        for name in list(set(sys.modules) - self._baseline_modules):
            module = sys.modules.get(name)
            filename = getattr(module, "__file__", None)
            if filename is None:
                # Namespace or builtin module - safe to leave in place.
                continue
            realname = os.path.realpath(filename)
            if not realname.startswith(self._installed_prefixes):
                del sys.modules[name]

    # -- job execution ------------------------------------------------------
    def _run_job(self, job_dir):
        """
        Execute a single PSyclone job described by ``job_dir``.

        :param str job_dir: directory containing ``request.json`` and the
            ``response`` FIFO.
        """
        from io import StringIO

        with open(os.path.join(job_dir, "request.json"),
                  encoding="utf8") as handle:
            request = json.load(handle)

        argv = request["argv"]
        cwd = request["cwd"]
        extra_paths = request.get("sys_path", [])

        # Preserve state that the job is allowed to mutate.
        saved_cwd = os.getcwd()
        saved_argv = sys.argv
        saved_stdout, saved_stderr = sys.stdout, sys.stderr
        saved_sys_path = list(sys.path)

        out, err = StringIO(), StringIO()
        returncode = 0
        try:
            os.chdir(cwd)
            # Make the recipe search locations (PYTHONPATH from the client,
            # e.g. LFRIC_BUILD/psyclone for psyclone_tools) importable.
            for path in reversed(extra_paths):
                if path and path not in sys.path:
                    sys.path.insert(0, path)
            sys.argv = ["psyclone"] + argv
            sys.stdout, sys.stderr = out, err
            from psyclone.generator import main as psyclone_main
            try:
                psyclone_main(argv)
            except SystemExit as exit_error:
                code = exit_error.code
                returncode = 0 if code is None else (
                    code if isinstance(code, int) else 1)
        except Exception:  # pylint: disable=broad-except
            traceback.print_exc(file=err)
            returncode = 1
        finally:
            sys.stdout, sys.stderr = saved_stdout, saved_stderr
            sys.argv = saved_argv
            sys.path[:] = saved_sys_path
            os.chdir(saved_cwd)
            self._reset_state()

        self._respond(job_dir, returncode, out.getvalue(), err.getvalue())

    @staticmethod
    def _respond(job_dir, returncode, stdout, stderr):
        """Write the job result to the per-job response FIFO."""
        payload = json.dumps(
            {"returncode": returncode, "stdout": stdout, "stderr": stderr})
        response_fifo = os.path.join(job_dir, "response")
        try:
            # Opening for write blocks until the client opens for read.
            with open(response_fifo, "w", encoding="utf8") as handle:
                handle.write(payload)
        except OSError:
            # Client gave up (e.g. fell back to a direct invocation); nothing
            # more we can usefully do.
            pass

    # -- main loop ----------------------------------------------------------
    def serve(self):
        """Pull jobs from the queue until asked to shut down."""
        # Workers should ignore SIGINT; the parent handles shutdown.
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        while True:
            job_dir = self._queue.get()
            if job_dir is _SHUTDOWN:
                break
            try:
                self._run_job(job_dir)
            except Exception:  # pylint: disable=broad-except
                # A failure here must never take the worker down.
                traceback.print_exc(file=sys.stderr)


def _worker_entry(queue):
    """Entry point for a worker process."""
    _Worker(queue).serve()


class PsycloneServer:
    """
    The dispatcher process: owns the request FIFO, the worker pool and the
    idle-timeout based shutdown.

    :param str server_dir: directory holding the request FIFO and job dirs.
    :param int workers: number of worker processes to pre-fork.
    :param float idle_timeout: seconds of inactivity before shutting down.
    """

    def __init__(self, server_dir, workers, idle_timeout):
        self._server_dir = server_dir
        self._workers = max(1, int(workers))
        self._idle_timeout = float(idle_timeout)
        self._request_fifo = os.path.join(server_dir, "request.fifo")
        self._pid_file = os.path.join(server_dir, "server.pid")
        self._ready_file = os.path.join(server_dir, "server.ready")
        self._queue = Queue()
        self._pool = []

    def _prefork(self):
        """Start the worker pool."""
        for _ in range(self._workers):
            process = Process(target=_worker_entry, args=(self._queue,))
            process.daemon = True
            process.start()
            self._pool.append(process)

    def _create_request_fifo(self):
        """Create the shared request FIFO if it does not already exist."""
        try:
            os.mkfifo(self._request_fifo)
        except FileExistsError:
            pass

    def _signal_ready(self):
        """Announce that the server is accepting work."""
        with open(self._ready_file, "w", encoding="utf8") as handle:
            handle.write(str(os.getpid()))

    def _shutdown(self):
        """Ask all workers to stop and clean up server files."""
        for _ in self._pool:
            self._queue.put(_SHUTDOWN)
        for process in self._pool:
            process.join(timeout=10)
            if process.is_alive():
                process.terminate()
        for path in (self._ready_file, self._pid_file, self._request_fifo):
            try:
                os.remove(path)
            except OSError:
                pass

    def serve(self):
        """Run the dispatch loop until the idle timeout expires."""
        self._create_request_fifo()
        self._prefork()

        # Open the read end first (non-blocking) so that opening the write end
        # does not deadlock waiting for a reader. The write handle is then kept
        # open for the lifetime of the server so the read end never sees EOF
        # between clients; reads block (with a timeout) rather than spinning.
        reader_fd = os.open(self._request_fifo, os.O_RDONLY | os.O_NONBLOCK)
        writer_fd = os.open(self._request_fifo, os.O_WRONLY)

        # Only advertise readiness once we can actually accept work.
        self._signal_ready()
        try:
            self._dispatch_loop(reader_fd)
        finally:
            os.close(reader_fd)
            os.close(writer_fd)
            self._shutdown()

    def _dispatch_loop(self, reader_fd):
        """Read enqueue messages and hand jobs to workers."""
        import select

        buffer = b""
        last_activity = time.monotonic()
        while True:
            ready, _, _ = select.select([reader_fd], [], [], 1.0)
            if ready:
                chunk = os.read(reader_fd, 65536)
                if chunk:
                    buffer += chunk
                    while b"\n" in buffer:
                        line, buffer = buffer.split(b"\n", 1)
                        job_dir = line.decode("utf8").strip()
                        if job_dir == "__STOP__":
                            return
                        if job_dir:
                            self._queue.put(job_dir)
                            last_activity = time.monotonic()
            if time.monotonic() - last_activity > self._idle_timeout:
                return


def main():
    """Command-line entry point used by the auto-start logic in the client."""
    server_dir = os.environ["PSYCLONE_SERVER_DIR"]
    workers = os.environ.get("PSYCLONE_WORKERS") or os.cpu_count() or 1
    idle_timeout = os.environ.get(
        "PSYCLONE_SERVER_IDLE_TIMEOUT", DEFAULT_IDLE_TIMEOUT)
    server = PsycloneServer(server_dir, workers, idle_timeout)
    server.serve()


if __name__ == "__main__":
    main()

