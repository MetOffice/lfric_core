##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT
# which you should have received as part of this distribution
##############################################################################
"""
Tests the nf-config wrapper.
"""

from pytest_subprocess.fake_process import FakeProcess

from fab.tools.category import Category
from nf_config import NfConfig


def call_list(fake_process: FakeProcess) -> list[list[str]]:
    """
    Converts FakeProcess calls to strings.

    :returns: List of argument strings per call.
    """
    result: list[list[str]] = []
    for call in fake_process.calls:
        result.append([str(arg) for arg in call])
    return result


def test_constructor() -> None:
    """
    Tests default constructor.
    """
    nfc = NfConfig()
    assert nfc.category == Category.MISC
    assert nfc.name == "nf-config"
    assert nfc.exec_name == "nf-config"
    assert nfc.get_flags() == []


def test_nf_config_check_available(fake_process: FakeProcess) -> None:
    """
    Tests availability functionality.
    """
    fake_process.register(['nf-config', '--version'],
                          returncode=0,
                          stdout="netCDF-Fortran 4.6.1")

    nfc = NfConfig()
    assert nfc.check_available()
    assert call_list(fake_process) == [["nf-config", "--version"]]


def test_nf_config_check_unavailable(fake_process: FakeProcess) -> None:
    """
    Tests availability failure.
    """
    fake_process.register(['nf-config', '--version'],
                          returncode=127,
                          stderr="command 'nf-config' not found")
    nfc = NfConfig()
    assert not nfc.check_available()
    assert call_list(fake_process) == [["nf-config", "--version"]]


def test_nf_config_compiler_flags(fake_process: FakeProcess) -> None:
    """
    Tests getting the compiler flags.
    """
    fake_process.register(['nf-config', '--fflags'],
                          returncode=0,
                          stdout="-I /somewhere")
    nfc = NfConfig()
    assert nfc.get_compiler_flags() == ["-I", "/somewhere"]
    assert call_list(fake_process) == [["nf-config", "--fflags"]]


def test_nf_config_linker_flags(fake_process: FakeProcess) -> None:
    """
    Tests availability failure.
    """
    fake_process.register(['nf-config', '--flibs'],
                          returncode=0,
                          stdout="-L /somewhere -lsomewhat")
    nfc = NfConfig()
    assert nfc.get_linker_flags() == ["-L", "/somewhere", "-lsomewhat"]
    assert call_list(fake_process) == [["nf-config", "--flibs"]]
