# -----------------------------------------------------------------------------
#  (C) Crown copyright 2026 Met Office. All rights reserved.
#  The file LICENCE, distributed with this code, contains details of the terms
#  under which the code may be used.
# -----------------------------------------------------------------------------
"""
Process-inspection helpers shared by :mod:`psyclone_client` and
:mod:`psyclone_server`.

The persistent PSyclone server is deliberately detached from the build (it must
outlive the individual make recipe that happens to start it) so it cannot rely
on the usual parent/child signalling to know when to stop. Instead its lifetime
is *pinned to the make process* that owns the build: the server polls that pid
and exits as soon as it disappears, whether the build ended normally or was
killed manually with SIGINT, SIGTERM or SIGKILL.

Identifying "the make process" is done here:

* :func:`outermost_make_pid` walks the process ancestry via ``/proc`` and
  returns the *outermost* make in the chain, so that recursive sub-makes (the
  LFRic build runs PSyclone from ``$(MAKE) -f psyclone_psykal.mk``) all resolve
  to the same top-level build;
* :func:`owner_alive` reports whether that process is still running, treating a
  zombie as dead and using the process start time to defend against pid reuse.

Note that watching the *process group* is not sufficient: when make is run from
a non-interactive shell script (as it is by rose-stem/cylc) there is no job
control, so make is not the group leader and the group outlives it.
"""

import os


# Executable names that count as "make" when scanning the process ancestry.
MAKE_NAMES = frozenset(("make", "gmake", "gnumake", "remake", "dmake",
                        "colormake", "bmake"))

# Whether procfs is available; when it is, a missing entry means the process
# really has gone rather than that we cannot tell.
HAVE_PROCFS = os.path.isdir("/proc/self")

# Guard against pathological (or looping) ancestries.
_MAX_ANCESTRY = 128


def read_proc_stat(pid):
    """
    Read selected fields from ``/proc/<pid>/stat``.

    :param int pid: process to inspect.
    :returns: (comm, state, ppid, start_time) or None if unavailable.
    :rtype: tuple[str, str, int, str] or None
    """
    try:
        with open(f"/proc/{int(pid)}/stat", encoding="utf8",
                  errors="replace") as handle:
            data = handle.read()
    except (OSError, ValueError, TypeError):
        return None
    # Format: "pid (comm) state ppid ...". The comm field is bracketed but may
    # itself contain spaces and brackets, so anchor on the first " (" and the
    # last ") ".
    try:
        _, _, remainder = data.partition(" (")
        comm, _, tail = remainder.rpartition(") ")
        fields = tail.split()
        # tail starts at field 3 (state), so ppid is [1] and starttime (field
        # 22) is [19].
        return comm, fields[0], int(fields[1]), fields[19]
    except (IndexError, ValueError):
        return None


def process_alive(pid):
    """
    Return True if a process with the given pid currently exists.

    :param pid: process identifier.
    :type pid: int or None
    :rtype: bool
    """
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        # Exists, but is owned by somebody else.
        return True
    except (OSError, ValueError, TypeError):
        return False
    return True


def start_time(pid):
    """
    Return the start time of a process, used to detect pid reuse.

    :param pid: process identifier.
    :type pid: int or None
    :returns: the raw procfs start-time field, or "" if unknown.
    :rtype: str
    """
    stat = read_proc_stat(pid) if pid else None
    return stat[3] if stat else ""


def is_make(comm):
    """
    Return True if a process name looks like a make program.

    :param str comm: the executable name reported by procfs.
    :rtype: bool
    """
    name = os.path.basename(comm or "").strip()
    return name in MAKE_NAMES or name.endswith("make")


def ancestry(pid=None, stat_reader=None):
    """
    Walk the process tree upwards from ``pid``.

    :param pid: process to start from (defaults to this process).
    :type pid: int or None
    :param stat_reader: injectable replacement for :func:`read_proc_stat`
        (used by the tests).
    :type stat_reader: callable or None
    :returns: [(pid, comm), ...] ordered from ``pid`` outwards.
    :rtype: list[tuple[int, str]]
    """
    reader = stat_reader or read_proc_stat
    current = os.getpid() if pid is None else int(pid)
    chain = []
    seen = set()
    while current > 1 and current not in seen and len(chain) < _MAX_ANCESTRY:
        seen.add(current)
        stat = reader(current)
        if stat is None:
            break
        comm, _state, ppid, _start = stat
        chain.append((current, comm))
        current = ppid
    return chain


def outermost_make_pid(pid=None, stat_reader=None):
    """
    Find the top-level make process owning this build.

    The *outermost* make is chosen so that the many recursive sub-makes used by
    the LFRic build all agree on a single owner - and therefore share a single
    PSyclone server which lives exactly as long as the build does.

    :param pid: process to start the search from (defaults to this process).
    :type pid: int or None
    :param stat_reader: injectable replacement for :func:`read_proc_stat`.
    :type stat_reader: callable or None
    :returns: pid of the outermost make ancestor, or None if there is none.
    :rtype: int or None
    """
    makes = [process for process, comm
             in ancestry(pid=pid, stat_reader=stat_reader) if is_make(comm)]
    return makes[-1] if makes else None


def owner_alive(pid, expected_start=None, stat_reader=None):
    """
    Return True if the owning build process is still running.

    A process that has become a zombie (killed but not yet reaped) counts as
    dead so that the server does not linger during that window, as does a pid
    whose start time no longer matches - that pid has been reused by an
    unrelated process.

    :param pid: the owner's process identifier; None means "no owner", in
        which case the caller relies solely on its idle timeout.
    :type pid: int or None
    :param expected_start: start time recorded when the owner was resolved.
    :type expected_start: str or None
    :param stat_reader: injectable replacement for :func:`read_proc_stat`.
    :type stat_reader: callable or None
    :rtype: bool
    """
    if not pid:
        return True
    reader = stat_reader or read_proc_stat
    stat = reader(pid)
    if stat is not None:
        _comm, state, _ppid, started = stat
        if state == "Z":
            return False
        if expected_start and started and started != expected_start:
            return False
        return True
    if HAVE_PROCFS:
        # procfs works, so an absent entry means the process has gone.
        return False
    # No procfs (non-Linux): fall back to a signal probe.
    return process_alive(pid)

