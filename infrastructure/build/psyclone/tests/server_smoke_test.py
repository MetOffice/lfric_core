# -----------------------------------------------------------------------------
#  (C) Crown copyright 2026 Met Office. All rights reserved.
#  The file LICENCE, distributed with this code, contains details of the terms
#  under which the code may be used.
# -----------------------------------------------------------------------------
"""
Smoke test for the persistent PSyclone server (psyclone_server.py) and its
client (psyclone_client.py).

It verifies the properties that the in-process, state-reset design depends on:

1. *Fidelity* - PSy layer produced via the server is byte-for-byte identical
   to that produced by invoking the ``psyclone`` binary directly.
2. *State isolation* - two consecutive jobs handled by the same long-lived
   server, each using a differently-behaving optimisation recipe that happens
   to share the module basename ``global``, do not contaminate one another
   (PSyclone loads recipes by bare basename, so a naive cache would reuse the
   first recipe for the second job).
3. *Lifetime* - the server is pinned to the owning make process: it exits when
   that process does, even when the surrounding process group survives (as it
   does when make is run from a shell script, e.g. by rose-stem).

The test can be run directly (``python server_smoke_test.py``) or under pytest.
"""

import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time


HERE = os.path.dirname(os.path.abspath(__file__))
# infrastructure/build/psyclone/tests -> repository root (develop).
REPO_ROOT = os.path.realpath(os.path.join(HERE, "..", "..", "..", ".."))
PSYCLONE_DIR = os.path.realpath(os.path.join(HERE, ".."))
CLIENT = os.path.join(PSYCLONE_DIR, "psyclone_client.py")
CONFIG = os.path.join(REPO_ROOT, "etc", "psyclone.cfg")

sys.path.insert(0, PSYCLONE_DIR)
import psyclone_procs  # noqa: E402  (needs PSYCLONE_DIR on the path)


ALGORITHM = """\
module smoke_alg_mod
  use constants_mod, only: r_def
  use field_mod,     only: field_type

  implicit none

contains

  subroutine smoke_alg( one )
    implicit none
    type(field_type), intent(inout) :: one
    call invoke(setval_c(one, 1.0_r_def))
  end subroutine smoke_alg

end module smoke_alg_mod
"""

# A recipe that duplicates the loop body once (two copies).
RECIPE_DUPLICATE_ONCE = '''\
"""Duplicate the first loop body once."""
from psyclone.psyGen import InvokeSchedule


def trans(psyir):
    for subroutine in psyir.walk(InvokeSchedule):
        loop = subroutine.loops()[0]
        loop.loop_body.addchild(loop.loop_body[0].copy())
    return psyir
'''

# A recipe with the *same* module basename but different behaviour: duplicates
# the loop body twice (three copies). Used to prove the two jobs are isolated.
RECIPE_DUPLICATE_TWICE = '''\
"""Duplicate the first loop body twice."""
from psyclone.psyGen import InvokeSchedule


def trans(psyir):
    for subroutine in psyir.walk(InvokeSchedule):
        loop = subroutine.loops()[0]
        loop.loop_body.addchild(loop.loop_body[0].copy())
        loop.loop_body.addchild(loop.loop_body[0].copy())
    return psyir
'''


def _run(command, workspace, use_server):
    """
    Run PSyclone either directly or through the server client.

    :param list command: PSyclone arguments (after the program name).
    :param str workspace: directory to run in.
    :param bool use_server: True to use the client, False for direct psyclone.
    :returns: CompletedProcess.
    """
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [PSYCLONE_DIR, env.get("PYTHONPATH", "")])
    if use_server:
        env["PSYCLONE_SERVER_DIR"] = os.path.join(workspace, ".server")
        program = [sys.executable, CLIENT]
    else:
        env["PSYCLONE_SERVER_DISABLE"] = "1"
        program = ["psyclone"]
    return subprocess.run(program + command, cwd=workspace, env=env,
                          capture_output=True, text=True, check=False)


def _generate(workspace, recipe_text, tag, use_server):
    """
    Generate a PSy layer for the shared algorithm using ``recipe_text``.

    :returns: (alg_source, psy_source) as strings.
    """
    recipe_dir = os.path.join(workspace, tag)
    os.makedirs(recipe_dir, exist_ok=True)
    # Deliberately share the basename "global" across recipes to exercise the
    # module-cache isolation.
    recipe = os.path.join(recipe_dir, "global.py")
    with open(recipe, "w", encoding="utf8") as handle:
        handle.write(recipe_text)

    alg_in = os.path.join(workspace, f"{tag}_alg.x90")
    with open(alg_in, "w", encoding="utf8") as handle:
        handle.write(ALGORITHM)

    alg_out = os.path.join(workspace, f"{tag}_alg.f90")
    psy_out = os.path.join(workspace, f"{tag}_psy.f90")
    kern_dir = os.path.join(workspace, "kernel")
    os.makedirs(kern_dir, exist_ok=True)

    result = _run(
        ["-api", "lfric", "-l", "all", "--config", CONFIG,
         "-s", recipe, "-okern", kern_dir,
         "-oalg", alg_out, "-opsy", psy_out, alg_in],
        workspace, use_server)
    assert result.returncode == 0, (
        f"PSyclone failed ({'server' if use_server else 'direct'}):\n"
        f"{result.stdout}\n{result.stderr}")

    with open(alg_out, encoding="utf8") as handle:
        alg = handle.read()
    with open(psy_out, encoding="utf8") as handle:
        psy = handle.read()
    return alg, psy


def test_server_matches_direct_and_isolates_state():
    """Server output matches direct psyclone and jobs stay isolated."""
    workspace = tempfile.mkdtemp(prefix="psyclone-smoke-")
    try:
        # Reference outputs from direct psyclone invocations.
        ref_once_alg, ref_once_psy = _generate(
            workspace, RECIPE_DUPLICATE_ONCE, "ref_once", use_server=False)
        ref_twice_alg, ref_twice_psy = _generate(
            workspace, RECIPE_DUPLICATE_TWICE, "ref_twice", use_server=False)

        # Two consecutive jobs through the SAME server, recipes sharing the
        # basename "global" but behaving differently.
        srv_once_alg, srv_once_psy = _generate(
            workspace, RECIPE_DUPLICATE_ONCE, "srv_once", use_server=True)
        srv_twice_alg, srv_twice_psy = _generate(
            workspace, RECIPE_DUPLICATE_TWICE, "srv_twice", use_server=True)

        # 1. Fidelity: server output is byte-identical to direct output.
        assert srv_once_psy == ref_once_psy, "server PSy differs from direct"
        assert srv_once_alg == ref_once_alg, "server alg differs from direct"

        # 2. State isolation: the second job used its own recipe, not a cached
        # "global" module from the first job.
        assert srv_twice_psy == ref_twice_psy, "second job was contaminated"
        assert srv_once_psy != srv_twice_psy, "recipes should differ"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_server_shuts_down_when_owner_dies():
    """The server exits promptly (not after its idle timeout) once the owning
    build process is killed, so a manually interrupted build leaves no
    orphaned servers behind."""
    server_module = os.path.join(PSYCLONE_DIR, "psyclone_server.py")
    workspace = tempfile.mkdtemp(prefix="psyclone-owner-")
    try:
        # A stand-in for the top-level make process.
        owner = subprocess.Popen(["sleep", "300"])
        env = dict(os.environ,
                   PSYCLONE_SERVER_DIR=workspace,
                   PSYCLONE_WORKERS="2",
                   # A long idle timeout so that only owner-death can stop the
                   # server quickly.
                   PSYCLONE_SERVER_IDLE_TIMEOUT="300",
                   PSYCLONE_OWNER_PID=str(owner.pid))
        server = subprocess.Popen([sys.executable, server_module], env=env)
        try:
            ready = os.path.join(workspace, "server.ready")
            for _ in range(120):
                if os.path.exists(ready):
                    break
                time.sleep(0.25)
            assert os.path.exists(ready), "server never became ready"

            # Kill the owner; the server should follow within a few seconds.
            owner.kill()
            owner.wait()
            for _ in range(100):
                if server.poll() is not None:
                    break
                time.sleep(0.1)
            assert server.poll() is not None, (
                "server did not shut down after its owner was killed")
            assert not os.path.exists(ready), "server did not clean up"
        finally:
            if server.poll() is None:
                server.kill()
            if owner.poll() is None:
                owner.kill()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


# A process tree of the shape produced by a rose-stem/cylc build:
#   bash (job script, and process group leader)
#     make            <- top-level build; the owner we want
#       sh
#         make        <- recursive sub-make running psyclone_psykal.mk
#           sh
#             python  <- the client
# Values are (comm, state, ppid, start_time) as returned by read_proc_stat.
FAKE_TREE = {
    611: ("python3", "R", 610, "9006"),
    610: ("sh", "S", 609, "9005"),
    609: ("make", "S", 608, "9004"),
    608: ("sh", "S", 607, "9003"),
    607: ("make", "S", 606, "9002"),
    606: ("bash", "S", 1, "9001"),
}


def test_owner_is_the_outermost_make_process():
    """The owner resolves to the top-level make, not a sub-make and not the
    surrounding shell (which is the process group leader in a script)."""
    reader = FAKE_TREE.get

    assert psyclone_procs.outermost_make_pid(pid=611, stat_reader=reader) == 607
    # Every level of recursive make agrees on the same owner, so the whole
    # build shares one server.
    assert psyclone_procs.outermost_make_pid(pid=609, stat_reader=reader) == 607
    # Nothing is claimed when make is not involved at all.
    assert psyclone_procs.outermost_make_pid(
        pid=606, stat_reader=reader) is None


def test_owner_alive_detects_death_zombies_and_pid_reuse():
    """Liveness must be robust: a reaped pid, an unreaped zombie and a reused
    pid all mean the build has gone."""
    reader = FAKE_TREE.get

    assert psyclone_procs.owner_alive(607, "9002", stat_reader=reader)
    # Gone from the process table.
    assert not psyclone_procs.owner_alive(999, "9002", stat_reader=reader)
    # Killed but not yet reaped by its parent.
    zombie = {**FAKE_TREE, 607: ("make", "Z", 606, "9002")}
    assert not psyclone_procs.owner_alive(607, "9002", stat_reader=zombie.get)
    # The pid has been recycled by an unrelated process.
    reused = {**FAKE_TREE, 607: ("sleep", "S", 1, "9999")}
    assert not psyclone_procs.owner_alive(607, "9002", stat_reader=reused.get)
    # No owner at all: the caller falls back to its idle timeout.
    assert psyclone_procs.owner_alive(None, stat_reader=reader)


# Helper run from a make recipe: starts a server exactly as the client would
# and reports the owner it resolved plus the directory it used.
START_HELPER = '''\
import sys
sys.path.insert(0, {psyclone_dir!r})
import psyclone_client as client

owner = client._resolve_owner()
server_dir, transient = client._server_dir(owner[0])
started = client._start_server(server_dir, owner, transient)
with open(sys.argv[1], "w") as handle:
    handle.write("%s\\n%s\\n%s\\n" % (started, server_dir, owner[0]))
'''

# The make process is deliberately detached from the test's own ancestry (the
# inner subshell exits, so make is reparented) and left in a process group led
# by a shell which outlives it. That is exactly the situation created by a
# rose-stem job script, and the reason watching the process group is not good
# enough.
BUILD_SCRIPT = '''\
( ( "$1" -s -C "$2" all > "$2/make.log" 2>&1 & ) & ) ; sleep 600
'''

BUILD_MAKEFILE = '''\
all:
	@echo $$PPID > {workspace}/make.pid
	@{python} {helper} {workspace}/server.info
	@sleep 600
'''


def _wait_for(path, timeout=120.0):
    """Wait for a file to appear, returning True if it did."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if os.path.exists(path):
            return True
        time.sleep(0.1)
    return False


def test_server_lifetime_is_pinned_to_the_make_process():
    """A server started by a build dies with its make process, even though the
    process group that make belongs to is still alive."""
    make = shutil.which("make") or shutil.which("gmake")
    if make is None:
        print("skipping: no make available")
        return

    workspace = tempfile.mkdtemp(prefix="psyclone-make-")
    wrapper = None
    server_dir = None
    try:
        helper = os.path.join(workspace, "start_server.py")
        with open(helper, "w", encoding="utf8") as handle:
            handle.write(START_HELPER.format(psyclone_dir=PSYCLONE_DIR))
        with open(os.path.join(workspace, "Makefile"), "w",
                  encoding="utf8") as handle:
            handle.write(BUILD_MAKEFILE.format(
                workspace=workspace, python=sys.executable, helper=helper))

        # A private temporary directory keeps the auto-created server
        # directory (and its clean-up) inside the workspace.
        env = dict(os.environ, TMPDIR=workspace)
        env.pop("PSYCLONE_SERVER_DIR", None)
        env.pop("PSYCLONE_OWNER_PID", None)
        wrapper = subprocess.Popen(
            ["sh", "-c", BUILD_SCRIPT, "build", make, workspace],
            env=env, start_new_session=True)

        info = os.path.join(workspace, "server.info")
        assert _wait_for(os.path.join(workspace, "make.pid")), \
            "make never ran the recipe"
        assert _wait_for(info), "the build never started a server"
        with open(info, encoding="utf8") as handle:
            started, server_dir, owner_pid = handle.read().split("\n")[:3]
        with open(os.path.join(workspace, "make.pid"), encoding="utf8") as pid:
            make_pid = int(pid.read().strip())

        assert started == "True", "the server did not become ready"
        # The owner is make itself, not the process group leader.
        assert int(owner_pid) == make_pid, (
            f"server pinned to {owner_pid}, expected make ({make_pid})")
        assert int(owner_pid) != os.getpgid(make_pid), (
            "test did not reproduce a process group led by something other "
            "than make")

        with open(os.path.join(server_dir, "server.pid"),
                  encoding="utf8") as handle:
            server_pid = int(handle.read().strip())
        assert psyclone_procs.process_alive(server_pid)

        # Kill make outright, as an impatient user would.
        os.kill(make_pid, signal.SIGKILL)

        deadline = time.monotonic() + 30.0
        while time.monotonic() < deadline:
            if not psyclone_procs.process_alive(server_pid):
                break
            time.sleep(0.1)
        assert not psyclone_procs.process_alive(server_pid), (
            "the server outlived the make process that owned it")
        # The process group is still alive, proving the server followed make
        # rather than the group.
        assert psyclone_procs.process_alive(wrapper.pid)
        assert not os.path.exists(server_dir), (
            "the server did not remove its working directory")
    finally:
        if wrapper is not None:
            try:
                os.killpg(os.getpgid(wrapper.pid), signal.SIGKILL)
            except OSError:
                pass
            wrapper.wait(timeout=10)
        shutil.rmtree(workspace, ignore_errors=True)


if __name__ == "__main__":
    test_owner_is_the_outermost_make_process()
    test_owner_alive_detects_death_zombies_and_pid_reuse()
    test_server_lifetime_is_pinned_to_the_make_process()
    test_server_matches_direct_and_isolates_state()
    test_server_shuts_down_when_owner_dies()
    print("PSyclone server smoke test passed")

