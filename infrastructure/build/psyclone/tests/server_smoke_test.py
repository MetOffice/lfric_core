# -----------------------------------------------------------------------------
#  (C) Crown copyright 2026 Met Office. All rights reserved.
#  The file LICENCE, distributed with this code, contains details of the terms
#  under which the code may be used.
# -----------------------------------------------------------------------------
"""
Smoke test for the persistent PSyclone server (psyclone_server.py) and its
client (psyclone_client.py).

It verifies two properties that the in-process, state-reset design depends on:

1. *Fidelity* - PSy layer produced via the server is byte-for-byte identical
   to that produced by invoking the ``psyclone`` binary directly.
2. *State isolation* - two consecutive jobs handled by the same long-lived
   server, each using a differently-behaving optimisation recipe that happens
   to share the module basename ``global``, do not contaminate one another
   (PSyclone loads recipes by bare basename, so a naive cache would reuse the
   first recipe for the second job).

The test can be run directly (``python server_smoke_test.py``) or under pytest.
"""

import os
import shutil
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


if __name__ == "__main__":
    test_server_matches_direct_and_isolates_state()
    test_server_shuts_down_when_owner_dies()
    print("PSyclone server smoke test passed")

