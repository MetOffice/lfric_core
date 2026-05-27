#!/usr/bin/env python3

##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT
# which you should have received as part of this distribution
##############################################################################

'''
This module reads in a psyclone_info.yaml file.
'''
from pathlib import Path
from typing import Optional, Union
import yaml

from fab.api import BuildConfig
from fab.fab_base.fab_base import FabBase


class PsycloneInfo:

    FILE_SPECIFIC = "file_specific"
    EXCLUDE = "exclude"

    def __init__(self, name: str, fab_base: FabBase) -> None:
        self._fab_base = fab_base
        # This will be initialised/updated each time when reading an info file.
        self._opt_path = Path()
        self._script_dir: str = ""
        self._name: str = name
        self._comment: str = ""
        self._api: str = ""
        self._artefacts: str = ""
        self._rules: list[tuple[str, list[str]]] = []

    def __str__(self) -> str:
        return self._name

    @property
    def name(self) -> str:
        """
        :returns: the name of this phas.
        """
        return self._name

    @property
    def comment(self) -> str:
        """
        :returns: the comment for this phase.
        """
        return self._comment

    @property
    def api(self) -> str:
        """
        :returns: the PSyclone command line options.
        """
        return self._api

    @property
    def artefacts(self) -> str:
        """
        :returns: the artefacts to apply this info to.
        """
        return self._artefacts

    @property
    def opt_path(self) -> Path:
        return self._opt_path

    def update(self, info: dict[str, str]) -> None:
        """
        Update this PSyclone information with data taken from a
        PSyclone info yaml file. This function is used to initially
        read in a first specification, or update a specification based
        on an additional file being read in later.

        :param info: the yaml information taken from a PSyclone info file.
        """
        for rule in info:
            if rule == "comment":
                self._comment = info["comment"]
            elif rule == "api":
                self._api = info["api"]
            elif rule == "artefacts":
                self._artefacts = info["artefacts"]
            elif rule == "script_dir":
                self._script_dir = info["script_dir"]
            else:
                self._read_rule(rule, info[rule])

        # Store the potentially updated optimisation root path, i.e. the site-
        # and platform-specific location, followed by a script dir (typically
        # transmute or psykal). This path is used in a few places.
        self._opt_path = (self._fab_base.config.source_root / "optimisation" /
                          f"{self._fab_base.site}-{self._fab_base.platform}" /
                          self._script_dir)

    def _read_rule(self, rule: str, file_list: str) -> None:
        """
        """
        # Support '*', which is a reserved character in yaml and needs to
        # be escaped or quoted.
        if file_list in ["\\*", "'*'", '"*"']:
            file_list = "*"
        self._rules.append((rule, file_list.split()))

    def view(self) -> str:
        s = f"""{self._name}:
comment: {self.comment}
api: {self.api}
artefacts: {self.artefacts}
script_dir: {self._script_dir}
rules: {self._rules}
"""
        return s

    def file_specific_script(self, fpath: Path) -> Optional[Path]:
        relative_path = None
        # The source file might be either in build_output (e.g. a preprocessed
        # .X90 file), or still in source (.x90 file). Check if the file
        # is in one of the two sub-trees, and use the relative path to
        # check if there is a file-specific optimisation script
        for base_path in [self._fab_base.config.source_root,
                          self._fab_base.config.build_output]:
            try:
                relative_path = fpath.relative_to(base_path)
            except ValueError:
                # The file is not under the `base_path` - keep on checking
                pass

        if relative_path:
            # The file was under either source or build. Check if there
            # is a file-specific optimisation script:
            local_transformation_script = (self.opt_path /
                                           (relative_path.with_suffix('.py')))
            if local_transformation_script.exists():
                return local_transformation_script
        return None

    def get_script(self, file: Path, config: BuildConfig) -> Optional[Path]:
        # Search starting from the end, so last rule wins
        file_str = str(file)

        for rule, file_list in self._rules[::-1]:
            for pattern in file_list:
                if pattern not in file_str and pattern != "*":
                    continue

                # Now the pattern matches. Check which rule is used
                # (note that file_specific might fall through in case that
                # there is no file-specific script)
                if rule == PsycloneInfo.FILE_SPECIFIC:
                    script = self.file_specific_script(file)
                    if script:
                        return script
                    if pattern == "*":
                        # Fall through, i.e. check for other rules
                        continue

                    # Now we have an explicit request for a file-specific
                    # script a file, but that script does not exist.
                    raise FileNotFoundError(
                        f"Cannot find explicitly requested script '{script}'.")

                elif rule == PsycloneInfo.EXCLUDE:
                    # Exclude pattern matches:
                    return None

                else:
                    opt_script = self.opt_path / rule
                    if not opt_script.exists():
                        raise FileExistsError(f"Cannot find script "
                                              f"'{opt_script}'.")
                    return opt_script

        return None


class PsycloneConfig:

    def __init__(self, fab_base: FabBase) -> None:
        self._fab_base = fab_base
        self._all_phases: list[str] = []
        self._psyclone_info: dict[str, PsycloneInfo] = {}

    @property
    def all_phases(self):
        return self._all_phases

    def get_info(self, phase: str) -> PsycloneInfo:
        return self._psyclone_info[phase]

    def view(self):
        s = f"""Phases: {" ".join(self._all_phases)}\n\n"""
        for phase in self._all_phases:
            s += f"{self._psyclone_info[phase].view()}\n"
        return s

    def read(self, filename: Union[str, Path]) -> None:

        with open(filename, "r", encoding="utf8") as stream:
            dependencies = yaml.safe_load(stream)

        # First take phases (if available)
        if dependencies.get("phases", None):
            self._all_phases = dependencies["phases"]

        for key in dependencies:
            if key == "phases":
                # Already handled
                continue
            if key not in self._psyclone_info:
                self._psyclone_info[key] = PsycloneInfo(key, self._fab_base)

            self._psyclone_info[key].update(dependencies[key])
