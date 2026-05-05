##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Author: J. Lyu, Bureau of Meteorology
# Author: J. Henrichs, Bureau of Meteorology

"""
Tests the LFRicBase class.
"""

from pathlib import Path
import os
import sys
import argparse
import inspect
from unittest import mock
from typing import List, Optional

import pytest

from fab.api import (ArtefactSet, BuildConfig, Category, ToolRepository,
                     Linker)
from fab.tools.compiler import CCompiler, FortranCompiler

from lfric_base import LFRicBase


class MockSiteConfig:
    """
    Creates a mock site config class.
    """
    def __init__(self) -> None:
        self.args: Optional[argparse.Namespace] = None

    def get_valid_profiles(self) -> List[str]:
        """
        :return: list of valid compilation profiles.
        """
        return ["default-profile"]

    def update_toolbox(self, build_config: BuildConfig) -> None:
        """
        Dummy function where the tool box could be modified
        """

    def handle_command_line_options(self, args: argparse.Namespace) -> None:
        """
        Simple function to handle command line options.
        """
        self.args = args

    def get_path_flags(self, _build_config: BuildConfig) -> List[str]:
        """
        :returns: list of path-specific flags.
        """
        return []


@pytest.fixture(name="stub_fortran_compiler", scope='function')
def stub_fortran_compiler_init() -> FortranCompiler:
    """
    Provides a minimal Fortran compiler.
    """
    compiler = FortranCompiler('some Fortran compiler', 'sfc', 'stub',
                               r'([\d.]+)', openmp_flag='-omp',
                               module_folder_flag='-mods')
    return compiler


@pytest.fixture(name="stub_c_compiler", scope='function')
def stub_c_compiler_init() -> CCompiler:
    """
    Provides a minimal C compiler.
    """
    compiler = CCompiler("some C compiler", "scc", "stub",
                         version_regex=r"([\d.]+)", openmp_flag='-omp')
    return compiler


@pytest.fixture(name="stub_linker", scope='function')
def stub_linker_init(stub_c_compiler) -> Linker:
    """
    Provides a minimal linker.
    """
    linker = Linker(stub_c_compiler, None, 'sln')
    return linker


@pytest.fixture(scope="function", autouse=True)
def setup_site_specific_config_environment(tmp_path):
    """
    This sets up the environment for the mocked site_specific config class
    (MockSiteConfig) to be used by tests of LFRicBase class methods without
    errors. This fixture is automatically executed for any test in this file.
    """
    # Creates mock module with __file__ attribute
    mock_site_module = mock.MagicMock()
    mock_site_module.Config = MockSiteConfig
    mock_site_module.__file__ = str(tmp_path / "site_specific" /
                                    "default" / "config.py")

    # Mocks site-specific imports
    sys.modules['site_specific'] = mock.MagicMock()
    sys.modules['site_specific.default'] = mock.MagicMock()
    sys.modules['site_specific.default.config'] = mock_site_module

    # Clears environment variables
    with mock.patch.dict(os.environ, clear=True):
        yield

    # Cleanups
    for module in ['site_specific', 'site_specific.default',
                   'site_specific.default.config']:
        if module in sys.modules:
            del sys.modules[module]


@pytest.fixture(scope="function", autouse=True)
def setup_tool_repository(stub_fortran_compiler, stub_c_compiler,
                          stub_linker):
    '''
    This sets up a ToolRepository that allows the LFRicBase class
    to proceed without raising errors. This fixture is automatically
    executed for any test in this file.
    '''
    # pylint: disable=protected-access
    # Make sure we always get a new ToolRepository to not be affected by
    # other tests:
    ToolRepository._singleton = None

    # Remove all compiler and linker, so we get results independent
    # of the software available on the platform this test is running
    tr = ToolRepository()
    for category in [Category.C_COMPILER, Category.FORTRAN_COMPILER,
                     Category.LINKER]:
        tr[category] = []

    # Add compilers and linkers, and mark them all as available,
    # as well as supporting MPI and OpenMP
    for tool in [stub_c_compiler, stub_fortran_compiler, stub_linker]:
        tool._mpi = True
        tool._openmp_flag = "-some-openmp-flag"
        tool._is_available = True
        tool._version = (1, 2, 3)
        tr.add_tool(tool)

    # Remove environment variables that could affect tests
    with mock.patch.dict(os.environ, clear=True):
        yield

    # Reset tool repository for other tests
    ToolRepository._singleton = None


def test_constructor(monkeypatch) -> None:
    '''
    Tests constructor.
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])
    lfric_base = LFRicBase(name="test_name")

    # Check root symbol defaults to name if not specified
    assert lfric_base.root_symbol == ["test_name"]

    # Check root symbol can be specified
    lfric_base = LFRicBase(name="test_name",
                           root_symbol="root1")
    assert lfric_base.root_symbol == ["root1"]

    # Check root symbol list
    lfric_base = LFRicBase(name="test_name",
                           root_symbol=["root1", "root2"])
    assert lfric_base.root_symbol == ["root1", "root2"]


def test_get_directory(monkeypatch, tmp_path) -> None:
    '''
    Tests the correct setup of lfric_core_root and lfric_apps_root.
    '''

    # Create mock directory structure
    mock_core = tmp_path / "core"
    mock_core.mkdir(parents=True)

    # Create mock LFRic base file location
    mock_base_dir = mock_core / "lfric_build"
    mock_base_dir.mkdir(parents=True)
    mock_base_file = mock_base_dir / "lfric_base.py"
    mock_base_file.write_text("", encoding='utf-8')

    # Mock __file__ attribute
    monkeypatch.setattr('lfric_base.__file__', str(mock_base_file))

    mock_apps = tmp_path / "apps"
    mock_apps.mkdir()
    deps_file = mock_apps / "dependencies.sh"
    deps_file.write_text("", encoding='utf-8')

    mock_caller = mock_apps / "some_app" / "build.py"
    mock_caller.parent.mkdir(parents=True)
    mock_caller.write_text("", encoding='utf-8')

    # Create mock frame objects with proper structure
    def create_frame(filename):
        frame = mock.Mock()
        frame.f_globals = {'__file__': filename}
        return frame

    def create_frame_info(filename):
        return (create_frame(filename), filename, None, None, None, None,
                None, None)

    # Mock inspect.stack() to return our test callers with proper
    # frame info structure
    mock_stack = [
        create_frame_info(str(mock_base_file)),  # First call in base dir
        create_frame_info(str(mock_caller))      # Second call in apps
    ]
    monkeypatch.setattr('inspect.stack', lambda: mock_stack)
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])

    lfric_base = LFRicBase(name="test")

    # Verify core root is set correctly
    assert lfric_base.lfric_core_root == mock_core


def test_precision_definition_without_default(monkeypatch) -> None:
    '''
    Tests specification of precision if no default precision is
    specified on the command line (--precision-default). Tests all
    other ways a precision can be specified: default command line,
    explicit command line, environment variable, and the per
    R_*PRECISION default.
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py",
                                      "--precision_other", "32"])
    monkeypatch.setattr(os, 'environ', {"R_BL_PRECISION": "64"})

    lfric_base = LFRicBase(name="test")
    lfric_base.define_preprocessor_flags_step()
    flags = lfric_base.preprocess_flags_common

    # Explicitly set on command line:
    assert '-DRDEF_PRECISION=32' in flags
    # Original default of this precision
    assert '-DR_SOLVER_PRECISION=32' in flags
    # Original default of this precision
    assert '-DR_TRAN_PRECISION=64' in flags
    # From environment variable
    assert '-DR_BL_PRECISION=64' in flags


@pytest.mark.parametrize('no_xios', [True, False])
@pytest.mark.parametrize('mpi', [True, False])
def test_preprocessor_flags(monkeypatch, no_xios, mpi) -> None:
    """
    Tests setting of preprocessor flags, and also that we get
    the expected defaults for the precision variables.
    """
    argv = ["fab_script", "--no-openmp"]
    if no_xios:
        argv.append("--no-xios")
    if not mpi:
        argv.append("--no-mpi")
    monkeypatch.setattr(sys, "argv", argv)

    # Mark the compiler to have MPI or not, depending on what is needed
    tr = ToolRepository()
    fc = tr.get_tool(Category.FORTRAN_COMPILER, "sfc")
    monkeypatch.setattr(fc, "_mpi", mpi)

    lfric_base = LFRicBase(name="test")
    lfric_base.define_preprocessor_flags_step()

    expected_flags = [
        '-DRDEF_PRECISION=64',
        '-DR_SOLVER_PRECISION=32',
        '-DR_TRAN_PRECISION=64',
        '-DR_BL_PRECISION=64'
    ]
    if not no_xios:
        expected_flags.append("-DUSE_XIOS")
    if not mpi:
        expected_flags.append("-DNO_MPI")
    assert set(lfric_base.preprocess_flags_common) == set(expected_flags)


def test_setup_site_specific_location(monkeypatch) -> None:
    '''
    Tests site specific path setup for LFRicBase.
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])
    lfric_base = LFRicBase(name="test")

    old_path = sys.path.copy()
    lfric_base.setup_site_specific_location()

    # Check paths added correctly
    base_dir = Path(inspect.getfile(LFRicBase)).parent
    assert str(base_dir) in sys.path
    assert str(base_dir / "site_specific") in sys.path

    # Restore path
    sys.path = old_path


def test_get_linker_flags(monkeypatch) -> None:
    '''
    Tests linker flags include required libraries.
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])

    lfric_base = LFRicBase(name="test")
    flags = lfric_base.get_linker_flags()

    expected_libs = ['yaxt', 'xios', 'netcdf', 'hdf5']
    assert set(flags) == set(expected_libs)


def test_grab_files_step(monkeypatch) -> None:
    '''
    Tests grabbing required source files
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])

    # Create mock objects
    mock_grab = mock.MagicMock()
    mock_core = Path("/mock/core")

    # Setup mocks
    monkeypatch.setattr('lfric_base.grab_folder', mock_grab)

    lfric_base = LFRicBase(name="test")
    monkeypatch.setattr(lfric_base, '_lfric_core_root', mock_core)

    # Call method under test
    lfric_base.grab_files_step()

    # Verify grab_folder called for all required directories
    expected_calls = [
        # Source directories
        mock.call(lfric_base.config,
                  src=mock_core/'infrastructure'/'source',
                  dst_label=''),
        mock.call(lfric_base.config,
                  src=mock_core/'components'/'driver'/'source',
                  dst_label=''),
        mock.call(lfric_base.config,
                  src=mock_core/'components'/'inventory'/'source',
                  dst_label=''),
        mock.call(lfric_base.config,
                  src=mock_core/'components'/'science'/'source',
                  dst_label=''),
        mock.call(lfric_base.config,
                  src=mock_core/'components'/'lfric-xios'/'source',
                  dst_label=''),
        # PSyclone config directory
        mock.call(lfric_base.config,
                  src=mock_core/'etc',
                  dst_label='psyclone_config')
    ]

    # Check both number of calls and call arguments
    assert mock_grab.call_count == len(expected_calls)
    mock_grab.assert_has_calls(expected_calls, any_order=True)


def test_find_source_files_step(monkeypatch) -> None:
    '''
    Tests finding and filtering source files
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])

    # Create mocks
    with (mock.patch('lfric_base.FabBase.find_source_files_step') as find_step,
          mock.patch('lfric_base.LFRicBase.templaterator_step') as temp_step,
          mock.patch('lfric_base.LFRicBase.configurator_step') as conf_step,
          mock.patch('lfric_base.Exclude') as mock_exclude):
        lfric_base = LFRicBase(name="test")
        lfric_base.find_source_files_step()

        # Verify exclusion filter added and super called
        mock_exclude.assert_called_once_with('unit-test', '/test/')
        find_step.assert_called_once()
        # Verify configurator and templaterator called
        conf_step.assert_called_once()
        temp_step.assert_called_once_with(lfric_base.config)


def test_configurator_step(monkeypatch) -> None:
    '''
    Tests the configurator setup and execution.
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])

    # Create mock objects
    mock_config = mock.MagicMock()
    mock_meta = mock.MagicMock(return_value="rose_meta.conf")

    # Set up mocks using monkeypatch
    monkeypatch.setattr('lfric_base.configurator', mock_config)

    lfric_base = LFRicBase(name="test")
    monkeypatch.setattr(lfric_base, 'get_rose_meta', mock_meta)

    with pytest.warns(match="_metric_send_conn not set, cannot send metrics"):
        lfric_base.configurator_step()

    # Verify configurator called with correct arguments
    mock_config.assert_called_once_with(
        lfric_base.config,
        lfric_core_source=lfric_base.lfric_core_root,
        rose_meta_conf="rose_meta.conf",
        include_paths=[],
    )


def test_templaterator_step(monkeypatch, tmp_path) -> None:
    '''
    Tests the templaterator step processes template files correctly.
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])

    # Create mock template file
    source_path = tmp_path / "source"
    source_path.mkdir(parents=True)
    template_file = source_path / "field.t90"
    template_file.write_text("template content", encoding='utf-8')

    # Create mock templaterator
    mock_templaterator = mock.MagicMock()
    mock_templaterator_instance = mock.MagicMock()
    mock_templaterator.return_value = mock_templaterator_instance
    monkeypatch.setattr('lfric_base.Templaterator', mock_templaterator)

    # Mock input_to_output_fpath
    mock_output_path = tmp_path / "build_output"
    mock_output_path.mkdir(parents=True)

    # Create mock config with proper artefact store
    mock_artefact_store = mock.MagicMock()
    mock_artefact_store.__getitem__.return_value = set()

    config = mock.MagicMock()
    config.artefact_store = mock_artefact_store
    config.build_output = mock_output_path
    config.source_root = source_path

    # Mock SuffixFilter to return our template file
    mock_filter = mock.MagicMock()
    mock_filter.return_value = {template_file}
    monkeypatch.setattr('lfric_base.SuffixFilter', lambda *args: mock_filter)

    # Create LFRicBase instance
    lfric_base = LFRicBase(name="test")
    monkeypatch.setattr(lfric_base, '_lfric_core_root', tmp_path)

    # Run templaterator step
    with pytest.warns(match="_metric_send_conn not set, cannot send metrics"):
        lfric_base.templaterator_step(config)

    # Verify templaterator initialization
    mock_templaterator.assert_called_once_with(tmp_path / "infrastructure" /
                                               "build" / "tools" /
                                               "Templaterator")

    # Verify template processing
    expected_calls = []
    templates = [
        {"kind": "real32", "type": "real"},
        {"kind": "real64", "type": "real"},
        {"kind": "int32", "type": "integer"}
    ]
    for template in templates:
        out_file = mock_output_path / f"field_{template['kind']}_mod.f90"
        out_file = mock_output_path / f"field_{template['kind']}_mod.f90"
        expected_calls.append(
            mock.call(template_file, out_file, key_values=template)
        )

    assert mock_templaterator_instance.process.call_count == 3
    mock_templaterator_instance.process.assert_has_calls(expected_calls)

    # Verify artefact store add calls
    expected_add_calls = []
    for template in templates:
        out_file = mock_output_path / f"field_{template['kind']}_mod.f90"
        expected_add_calls.append(
            mock.call(ArtefactSet.FORTRAN_COMPILER_FILES, out_file)
        )

    assert mock_artefact_store.add.call_count == 3
    mock_artefact_store.add.assert_has_calls(expected_add_calls)

    # Test empty template files case
    mock_filter.return_value = set()

    with pytest.warns(match="_metric_send_conn not set, cannot send metrics"):
        lfric_base.templaterator_step(config)
    # Call count should remain the same since no new files processed
    assert mock_templaterator_instance.process.call_count == 3


def test_get_rose_meta(monkeypatch) -> None:
    '''
    Tests getting rose meta configuration
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])

    lfric_base = LFRicBase(name="test")
    assert lfric_base.get_rose_meta() is None


def test_analyse_step(monkeypatch) -> None:
    '''Tests analysis step configuration and execution'''

    # Test case 1: No ignore_dependencies argument specified
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])

    # Create mocks
    mock_analyse = mock.MagicMock()
    mock_preprocess = mock.MagicMock()
    mock_psyclone = mock.MagicMock()

    # Setup mocks
    monkeypatch.setattr('fab.fab_base.fab_base.FabBase.analyse_step',
                        mock_analyse)

    lfric_base = LFRicBase(name="test")

    # Mock instance methods
    monkeypatch.setattr(lfric_base, 'preprocess_x90_step', mock_preprocess)
    monkeypatch.setattr(lfric_base, 'psyclone_step', mock_psyclone)

    # The PSyclone step will modify sys.path (to allow import of
    # psyclone_tools by PSyclone scripts). Make sure sys.path is unchanged:
    old_sys_path = sys.path[:]
    # Call analyse_step (which calls PSyclone)
    lfric_base.analyse_step()
    assert sys.path == old_sys_path

    # Verify method calls
    mock_preprocess.assert_called_once()
    mock_psyclone.assert_called_once()

    # Verify analyse called with correct default ignore_dependencies
    expected_ignore = ['netcdf', 'mpi', 'mpi_f08', 'yaxt',
                       'xios', 'icontext', 'mod_wait']
    mock_analyse.assert_called_once_with(
        ignore_dependencies=expected_ignore,
        find_programs=False
    )

    # Test case 2: Custom ignore_dependencies arguments specified
    custom_ignore = ['custom_dep1', 'custom_dep2']
    mock_analyse.reset_mock()
    mock_preprocess.reset_mock()
    mock_psyclone.reset_mock()

    lfric_base = LFRicBase(name="test")
    monkeypatch.setattr(lfric_base, 'preprocess_x90_step', mock_preprocess)
    monkeypatch.setattr(lfric_base, 'psyclone_step', mock_psyclone)

    # Call analyse_step
    lfric_base.analyse_step(ignore_dependencies=custom_ignore)

    # Verify methods still called
    mock_preprocess.assert_called_once()
    mock_psyclone.assert_called_once()

    # Verify analyse called with custom_ignore added to ignore list
    expected_ignore = ['custom_dep1', 'custom_dep2', 'netcdf', 'mpi',
                       'mpi_f08', 'yaxt', 'xios', 'icontext',
                       'mod_wait']
    mock_analyse.assert_called_once_with(
        ignore_dependencies=expected_ignore,
        find_programs=False
    )


def test_preprocess_x90_step(monkeypatch) -> None:
    '''
    Tests preprocessing of X90 files.
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])

    mock_preproc = mock.MagicMock()
    monkeypatch.setattr('lfric_base.preprocess_x90', mock_preproc)

    lfric_base = LFRicBase(name="test")
    lfric_base.add_preprocessor_flags(["-flag1", "-flag2"])
    lfric_base.preprocess_x90_step()

    mock_preproc.assert_called_once_with(
        lfric_base.config,
        common_flags=["-flag1", "-flag2"]
    )


def test_psyclone_step(monkeypatch) -> None:
    '''
    Tests the PSyclone step.
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])

    # Create mock objects
    mock_psy = mock.MagicMock()
    mock_psyclone_config = "/mock/psyclone.cfg"

    # Set up monkeypatch for module level import
    monkeypatch.setattr('lfric_base.psyclone', mock_psy)

    lfric_base = LFRicBase(name="test")

    # Patch instance methods. Return a copy to avoid that
    # PSyclone modified these lists in the lambdas when it modifies the list
    monkeypatch.setattr(lfric_base, 'get_psyclone_config',
                        lambda: mock_psyclone_config)

    # Call method under test
    lfric_base.psyclone_step(additional_parameters=["-additional"])

    # Verify psyclone called with correct arguments
    mock_psy.assert_called_once_with(
        lfric_base.config,
        kernel_roots=[(lfric_base.config.build_output / "kernel")],
        transformation_script=lfric_base.get_transformation_script,
        api="lfric",
        cli_args=(["--config", mock_psyclone_config, "-additional"]),
        ignore_dependencies=None
    )


def test_get_psyclone_config(monkeypatch) -> None:
    '''
    Tests getting PSyclone config.
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])

    lfric_base = LFRicBase(name="test")
    config_args = lfric_base.get_psyclone_config()

    assert config_args == str(lfric_base.config.source_root /
                              'psyclone_config/psyclone.cfg')


def test_get_transformation_script(monkeypatch, tmp_path) -> None:
    '''
    Tests finding PSyclone transformation scripts.
    '''
    monkeypatch.setattr(sys, "argv", ["lfric_base.py"])

    # Create LFRicBase instance with mocked site/platform
    lfric_base = LFRicBase(name="test")

    # Create mock config
    config = mock.MagicMock()
    config.source_root = tmp_path
    config.build_output = tmp_path / "build"
    config.build_output.mkdir()

    # Create x90 test source file
    source_path = tmp_path / "some/path"
    source_path.mkdir(parents=True)
    test_file = source_path / "file.x90"
    test_file.touch()

    # Test case 1: x90 file not in source or build directories
    outside_file = tmp_path.parent / "outside.x90"
    assert lfric_base.get_transformation_script(outside_file, config) is None

    # Test case 2: No optimisation directory, no transformation script
    assert lfric_base.get_transformation_script(test_file, config) is None

    # Test case 3: No PSykal but optimisation directory
    optimisation_folder_path = (tmp_path / "optimisation" / "default-default" /
                                "psykal")
    global_script = optimisation_folder_path / "global.py"
    global_script.parent.mkdir(parents=True)
    global_script.touch()

    # No file-specific transformation script, use global script
    other_file = tmp_path / "other/path/test.x90"
    other_file.parent.mkdir(parents=True)
    other_file.touch()
    assert (lfric_base.get_transformation_script(other_file, config) ==
            global_script)

    # Test case 4: Psykal directory exists
    psykal_path = tmp_path / "optimisation/default-default/psykal"

    # Create specific transformation script in psykal dir
    specific_script = psykal_path / "some/path/file.py"
    specific_script.parent.mkdir(parents=True)
    specific_script.touch()

    # Use specific script in psykal directory
    assert lfric_base.get_transformation_script(test_file, config) == \
        specific_script
