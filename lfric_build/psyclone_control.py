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


class PsycloneInfo:
    """
    This class stores the set of rules and settings for one specific PSyclone
    phase.

    :param name: the name of this phase.
    :param fab_base: the application script derived from FabBase. Required to
        get access to the build config for paths, and the selected site and
        platform.
    """

    # Define a new 'script paths' to indicate special results. Use : in the
    # name, which is typically not a valid name (at the start)
    RESULT_EXCLUDE = Path("::RESULT_EXCLUDE")
    RESULT_NO_SCRIPT = Path("::RESULT_NO_SCRIPT")

    FILE_SPECIFIC = "file_specific"
    EXCLUDE = "exclude"
    NO_SCRIPT = "no_script"

    def __init__(self, name: str,
                 base_paths: list[Path],
                 script_root: Path) -> None:
        self._base_paths = base_paths
        self._script_root = script_root
        # This will be initialised/updated each time when reading an info file.
        self._opt_path = Path()
        self._relative_script_dir: str = ""
        self._name: str = name
        self._comment: str = ""
        self._api: str = ""
        self._artefacts: str = ""
        self._rules: list[tuple[str, list[str]]] = []

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
        """
        :return: the optimisation root as absolute path (including site-
            and platform-specific settings, and subdirectory, e.g. psykal
            or transmute).
        """
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
                self._relative_script_dir = info["script_dir"]
            else:
                self._read_rule(rule, info[rule])

        # Store the potentially updated optimisation root path, i.e. the site-
        # and platform-specific location, followed by a script dir (typically
        # transmute or psykal). This path is used in a few places.
        self._opt_path = self._script_root / self._relative_script_dir

    def _read_rule(self, rule: str, file_list: str) -> None:
        """
        Parses a single rule from the yaml file. It especially handles
        various way a '*' can be specified in a yaml file.

        :param rule: the name of the rule (typically the script name, or
            special term like `exclude` or `file_specific`).
        :param file_list: the list of files to which to apply the rule to.
        """
        # Support '*', which is a reserved character in yaml and needs to
        # be escaped or quoted.
        if file_list in ["\\*", "'*'", '"*"']:
            file_list = "*"
        self._rules.append((rule, file_list.split()))

    def view(self) -> str:
        """
        :returns: a string representation of this phase in yaml format.
        """
        s = f"""{self._name}:
comment: {self.comment}
api: {self.api}
artefacts: {self.artefacts}
script_dir: {self._relative_script_dir}
rules: {self._rules}
"""
        return s

    def file_specific_script(self, fpath: Path) -> Optional[Path]:
        """
        Searches for a file-specific optimisation script. It will search
        both under the source and the build directories of the project
        directory.

        :param fpath: the file path of the Fortran file.
        :returns: the path of the file-specific optimisation script, or
            None if no such file exists.
        """
        relative_path = None
        # The source file might be either in build_output (e.g. a preprocessed
        # .X90 file), or still in source (.x90 file). Check if the file
        # is in one of the two sub-trees, and use the relative path to
        # check if there is a file-specific optimisation script
        for base_path in self._base_paths:
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

    def get_script(self, fpath: Path) -> Path:
        """
        This method returns the script to be used for a given filename, or
        None if no rule applies (or an explicit exclude rule applies)

        This function will also provided to Fab's PSyclone step, and as such
        it will receive the config object, even though it is not used.

        :param fpath: the Fortran source file for which to find a
            transformation script.

        :returns: the path to the transformation script, or None if no rule
            applies (or an exclude rule applies).
        """

        file_str = str(fpath)
        # Search starting from the end, so last rule wins
        for rule, file_list in self._rules[::-1]:
            for pattern in file_list:
                if pattern not in file_str and pattern != "*":
                    continue

                # Now the pattern matches. Check which rule is used
                # (note that file_specific might fall through in case that
                # there is no file-specific script)
                if rule == PsycloneInfo.FILE_SPECIFIC:
                    script = self.file_specific_script(fpath)
                    if script:
                        return script
                    if pattern == "*":
                        # Fall through, i.e. check for other rules
                        continue

                    # Now we have an explicit request for a file-specific
                    # script, but that script does not exist.
                    raise FileNotFoundError(
                        f"Cannot find explicitly requested script '{script}'.")

                if rule == PsycloneInfo.EXCLUDE:
                    # Exclude pattern matches:
                    return PsycloneInfo.RESULT_EXCLUDE

                if rule == PsycloneInfo.NO_SCRIPT:
                    # Exclude pattern matches:
                    return PsycloneInfo.RESULT_NO_SCRIPT

                opt_script = self.opt_path / rule
                if not opt_script.exists():
                    raise FileNotFoundError(f"Cannot find script "
                                            f"'{opt_script}'.")
                return opt_script

        return PsycloneInfo.RESULT_EXCLUDE


class PsycloneControl:
    """
    This class stores the information from psyclone_info.yaml file(s). Several
    files can be read, and latter information will extend the rules from
    previous files, and replace the phases executed.

    Details of each phase will be stored in PsycloneInfo instances.

    :param base_paths:
    :param fab_base: The FabBase derived application script. This is required
        to get site, platform and config information when searching for
        PSyclone scripts to be executed.
    """

    def __init__(self,
                 script_root: Path,
                 base_paths: list[Path]) -> None:
        # Keep a copy in case that the user modifies the list later
        self._base_paths = base_paths[:]
        self._script_root = script_root
        self._all_phases: list[str] = []
        self._psyclone_info: dict[str, PsycloneInfo] = {}

    @property
    def all_phases(self) -> list[str]:
        """
        :returns: the list of all PSyclone phases to execute.
        """
        return self._all_phases

    def get_info(self, phase: str) -> PsycloneInfo:
        """
        Returns the PSyclone information for the specified phase.

        :param phase: the name of the phase.

        :returns: the PSyclone Information for the specified phase.
        """
        return self._psyclone_info[phase]

    def view(self) -> str:
        """
        This returns a string representation of the combined read yaml files.
        This is useful to be logged to show the actual details used.

        :returns: the string represenation in yaml format of this PSyclone
            control instance.
        """
        s = f"""Phases: {" ".join(self._all_phases)}\n\n"""
        for phase in self._all_phases:
            s += f"{self._psyclone_info[phase].view()}\n"
        return s

    def read(self, filename: Union[str, Path]) -> None:
        """
        Reads a yaml file, and extends the potentially existing information.
        Any phases specified in the new read yaml file will replace the
        phases to be executed (i.e. will overwrite what was previously
        specified). Any new rule sets will be added as new PsycloneInfo
        instance. Rules for an exiting phase will be appended to the
        existing information. The precedence handling means that any later
        rule will overwrite any previous rule.

        :param filename: the filename to read.
        """

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
                self._psyclone_info[key] = PsycloneInfo(key,
                                                        self._base_paths,
                                                        self._script_root)

            self._psyclone_info[key].update(dependencies[key])
