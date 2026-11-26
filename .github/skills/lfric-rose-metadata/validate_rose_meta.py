#!/usr/bin/env python3
"""Run Rose's metadata checks on one or more metadata files or directories."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Iterable, List, Optional, Sequence


def _repository_root() -> Path:
    script_path = Path(__file__).resolve()
    for parent in script_path.parents:
        if (parent / ".git").exists():
            return parent
    return script_path.parents[3]


def _metadata_files(paths: Sequence[Path]) -> List[Path]:
    files = []
    seen = set()

    def add(path: Path) -> None:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            files.append(resolved)

    for path in paths:
        if not path.exists():
            raise ValueError("Path does not exist: {}".format(path))
        if path.is_file():
            if path.name != "rose-meta.conf":
                raise ValueError("Expected a rose-meta.conf file: {}".format(path))
            add(path)
            continue

        direct_metadata = path / "rose-meta.conf"
        if direct_metadata.is_file():
            add(direct_metadata)
            continue

        head_files = []
        for directory, subdirectories, filenames in os.walk(path):
            subdirectories[:] = [
                name for name in subdirectories
                if name not in {".git", "__pycache__", "working"}
            ]
            if Path(directory).name == "HEAD" and "rose-meta.conf" in filenames:
                head_files.append(Path(directory) / "rose-meta.conf")
        for metadata_file in head_files:
            add(metadata_file)

    if not files:
        raise ValueError(
            "No HEAD/rose-meta.conf files found; specify a version directory "
            "or rose-meta.conf file to check another version"
        )
    return files


def _repository_metadata_roots(repository_root: Path) -> List[Path]:
    roots = []
    for directory, subdirectories, _ in os.walk(repository_root):
        subdirectories[:] = [
            name for name in subdirectories
            if name not in {".git", "__pycache__", "working"}
        ]
        if Path(directory).name == "rose-meta":
            roots.append(Path(directory))
            subdirectories[:] = []
    return roots


def _metadata_roots(
    files: Iterable[Path],
    extra_roots: Sequence[Path],
    repository_root: Path,
) -> List[str]:
    roots = []
    seen = set()

    def add(path: Path) -> None:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            roots.append(str(resolved))

    for root in _repository_metadata_roots(repository_root):
        add(root)

    for metadata_file in files:
        for parent in metadata_file.parents:
            if parent.name == "rose-meta":
                add(parent)
                break

    for root in extra_roots:
        if not root.is_dir():
            raise ValueError("Metadata search root is not a directory: {}".format(root))
        add(root)

    for root in os.environ.get("ROSE_META_PATH", "").split(os.pathsep):
        if root:
            add(Path(root))

    return roots


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Validate Rose metadata using the installed `rose metadata-check` "
            "command. A metadata tree defaults to HEAD; supply a specific "
            "rose-meta.conf file or version directory to check another version."
        )
    )
    parser.add_argument(
        "paths",
        nargs="+",
        type=Path,
        help="Metadata files, version directories, or trees (trees select HEAD)",
    )
    parser.add_argument(
        "--meta-path",
        action="append",
        default=[],
        type=Path,
        metavar="DIR",
        help="Additional Rose metadata search root; may be supplied more than once",
    )
    parser.add_argument(
        "--source-root",
        default=_repository_root(),
        type=Path,
        metavar="DIR",
        help="Repository root to search for metadata (default: this checker's repository)",
    )
    parser.add_argument(
        "--rose",
        default="rose",
        help="Rose command to run (default: rose)",
    )
    args = parser.parse_args(argv)

    try:
        metadata_files = _metadata_files(args.paths)
        if not args.source_root.is_dir():
            raise ValueError("Source root is not a directory: {}".format(args.source_root))
        metadata_roots = _metadata_roots(
            metadata_files, args.meta_path, args.source_root
        )
    except ValueError as error:
        parser.error(str(error))

    rose_command = shutil.which(args.rose)
    if rose_command is None:
        parser.error("Rose command not found: {}".format(args.rose))

    environment = os.environ.copy()
    if metadata_roots:
        environment["ROSE_META_PATH"] = os.pathsep.join(metadata_roots)

    print("Repository source root: {}".format(args.source_root.resolve()))
    failures = 0
    for metadata_file in metadata_files:
        print("Checking {}".format(metadata_file))
        result = subprocess.run(
            [
                rose_command,
                "metadata-check",
                "--property=fail-if",
                "--property=warn-if",
                "-C",
                str(metadata_file.parent),
            ],
            env=environment,
            check=False,
        )
        if result.returncode:
            print(
                "FAIL: {} (rose metadata-check exit {})".format(
                    metadata_file, result.returncode
                ),
                file=sys.stderr,
            )
            failures += 1
        else:
            print("PASS: {}".format(metadata_file))

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
