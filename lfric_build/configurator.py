##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Author J. Henrichs, Bureau of Meteorology
# Author J. Lyu, Bureau of Meteorology

"""
This file defines the configurator script sequence for LFRic.
"""

import logging
from pathlib import Path
from typing import cast, Optional

from fab.api import BuildConfig, find_source_files, Category
from fab.tools.shell import Shell

from rose_picker import RosePicker

logger = logging.getLogger('fab')


def configurator(config: BuildConfig,
                 lfric_core_source: Path,
                 rose_meta_conf: Path,
                 include_paths: Optional[list[Path]] = None,
                 config_dir: Optional[Path] = None) -> None:
    """
    This method implements the LFRic configurator tool.

    :param config: the Fab build config instance
    :param lfric_core_source: the path to the LFRic core directory
    :param rose_meta_conf: the path to the rose-meta configuration file
    :param include_paths: additional include paths (each path will be added,
        as well as the path with /'rose-meta')
    :param config_dir: the directory for the generated configuration files
    """

    tools = lfric_core_source / 'infrastructure' / 'build' / 'tools'
    config_dir = config_dir or config.build_output / 'configuration'
    config_dir.mkdir(parents=True, exist_ok=True)

    # rose picker
    # -----------
    # creates rose-meta.json and config_namelists.txt in
    # gungho/build
    logger.info('rose_picker')

    include_dirs = [lfric_core_source / 'rose-meta']
    if include_paths:
        for path in include_paths:
            include_dirs.append(path / 'rose-meta')

    rose_picker = RosePicker()
    rose_picker.execute(rose_meta_conf, config_dir,
                        include_paths=include_dirs)
    rose_meta = config_dir / 'rose-meta.json'

    shell = config.tool_box.get_tool(Category.SHELL)
    shell = cast(Shell, shell)

    # build_config_loaders
    # --------------------
    # builds a bunch of f90s from the json
    logger.info('GenerateNamelistLoader')
    shell.exec(f"{tools / 'GenerateNamelistLoader'} -verbose {rose_meta} "
               f"-directory {config_dir}")

    # create configuration_mod.f90 in source root
    # -------------------------------------------
    logger.info('GenerateConfigLoader')
    with open(config_dir / 'config_namelists.txt', encoding="utf8") as f_in:
        names = [name.strip() for name in f_in.readlines()]

    shell.exec(f"{tools / 'GenerateConfigLoader'} "
               f"{' '.join(names)} "
               f"-o {config_dir}")

    logger.info('GenerateExtendedNamelistType')
    shell.exec(f"{tools / 'GenerateExtendedNamelistType'} {rose_meta} "
               f"-directory {config_dir}")

    duplicates: list[str] = []
    with open(config_dir / 'duplicate_namelists.txt', encoding="utf8") as f_in:
        for name in f_in.readlines():
            duplicates.extend(["-duplicate", name.strip()])

    logger.info('GenerateConfigType')
    shell.exec(f"{tools / 'GenerateConfigType'} "
               f"{' '.join(names)} "
               f"{' '.join(duplicates)} "
               f"-o {config_dir}")

    # create feign_config_mod.f90 in source root
    # ------------------------------------------
    logger.info('GenerateFeigns')
    feign_config_mod_fpath = config_dir / 'feign_config_mod.f90'
    shell.exec(f"{tools / 'GenerateFeigns'} {rose_meta} "
               f"-output {feign_config_mod_fpath}")

    find_source_files(config, source_root=config_dir)
