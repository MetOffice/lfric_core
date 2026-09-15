# ############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT
# which you should have received as part of this distribution
# ############################################################################

"""
Unit tests for the psyclone_control module.
"""

from pathlib import Path
import pytest
import yaml

from psyclone_control import PsycloneInfo, PsycloneControl


def test_psyclone_info_properties():
    """
    Test initial properties and getters of PsycloneInfo.
    """
    base_paths = [Path("/base1"), Path("/base2")]
    script_root = Path("/scripts")

    info = PsycloneInfo(
        name="test_phase", base_paths=base_paths, script_root=script_root
    )

    assert info.name == "test_phase"
    assert info.comment == ""
    assert info.api == ""
    assert info.opt_path == Path()


@pytest.mark.parametrize("escaped_wildcard", ["\\*", "'*'", '"*"'])
def test_psyclone_info_wildcard_handling(escaped_wildcard):
    """
    Verify different yaml wildcard string variants parse as expected.
    """
    info = PsycloneInfo("phase", [], Path())
    # Indirectly hit _read_rule via update
    yaml_dict = {"some_script.py": escaped_wildcard}
    info.update(yaml_dict)

    assert info._rules == [("some_script.py", ["*"])]


def test_psyclone_control_read_and_to_yaml(tmp_path):
    """
    Test full workflow: reading configuration files and generating YAML output.
    """
    yaml_content_1 = """
phases:
  - dsl
dsl:
  comment: "PSyclone DSL Phase"
  api: lfric
  script_dir: psykal
  global.py: \\*
  file_specific: \\*
"""
    yaml_file_1 = tmp_path / "psyclone_info.yaml"
    yaml_file_1.write_text(yaml_content_1, encoding="utf-8")

    # Second YAML payload tests appending/overwriting phase specifications
    yaml_content_2 = """
phases:
  - dsl
  - secondary
secondary:
  comment: "Secondary phase"
  script_dir: alternative
"""
    yaml_file_2 = tmp_path / "psyclone_info_override.yaml"
    yaml_file_2.write_text(yaml_content_2, encoding="utf-8")

    base_paths = [tmp_path / "src", tmp_path / "build"]
    script_root = tmp_path / "scripts"

    pc = PsycloneControl(script_root=script_root, base_paths=base_paths)

    # Read first file
    pc.read(yaml_file_1)
    assert pc.all_phases == ["dsl"]

    info_dsl = pc.get_info("dsl")
    assert info_dsl.comment == "PSyclone DSL Phase"
    assert info_dsl.api == "lfric"
    assert info_dsl.opt_path == script_root / "psykal"

    # Read second file to test incremental overrides
    pc.read(yaml_file_2)
    assert pc.all_phases == ["dsl", "secondary"]
    assert pc.get_info("secondary").comment == "Secondary phase"

    # Test YAML text generation output matches structures
    yaml_out = pc.to_yaml()
    assert f"#       {yaml_file_1.resolve()}" in yaml_out
    assert f"#       {yaml_file_2.resolve()}" in yaml_out

    parsed_out = yaml.safe_load(yaml_out)
    assert parsed_out["phases"] == ["dsl", "secondary"]
    assert parsed_out["dsl"]["api"] == "lfric"


def test_file_specific_script_resolution(tmp_path):
    """
    Validate looking up file-specific scripts inside source tree trees.
    """
    base_src = tmp_path / "src"
    base_build = tmp_path / "build"
    script_root = tmp_path / "scripts"

    # Define an active script path destination directory
    opt_dir = script_root / "psykal"
    opt_dir.mkdir(parents=True)

    info = PsycloneInfo(
        name="dsl", base_paths=[base_src, base_build], script_root=script_root
    )
    info.update({"script_dir": "psykal"})

    # Case 1: Target file is out of any base path boundary
    external_file = tmp_path / "outside" / "some_mod.x90"
    assert info.file_specific_script(external_file) is None

    # Case 2: Target inside src directory, but no companion exists
    src_file = base_src / "kernel" / "some_mod.x90"
    assert info.file_specific_script(src_file) is None

    # Case 3: Script target exists matching the relative structure
    expected_script = opt_dir / "kernel" / "some_mod.py"
    expected_script.parent.mkdir(parents=True, exist_ok=True)
    expected_script.touch()

    assert info.file_specific_script(src_file) == expected_script


def test_get_script_matching_logic(tmp_path):
    """
    Verify filtering behavior, fallback hierarchies, and exception conditions.
    """
    base_src = tmp_path / "src"
    script_root = tmp_path / "scripts"
    opt_dir = script_root / "psykal"
    opt_dir.mkdir(parents=True)

    info = PsycloneInfo(
        name="dsl", base_paths=[base_src], script_root=script_root
    )

    # Setup rules: non-matching pattern, NO_SCRIPT rule, then explicit scripts
    yaml_config = {
        "script_dir": "psykal",
        "exclude": "ignored_module.x90",
        "no_script": "skipped_module.x90",
        "missing_script.py": "broken_module.x90",
        "valid_script.py": "good_module.x90",
        "file_specific": "\\*",
    }
    info.update(yaml_config)

    # 1. Test standard pattern mismatch fallback
    assert (
        info.get_script(Path("unmatched_file.x90"))
        == PsycloneInfo.RESULT_EXCLUDE
    )

    # 2. Test explicit EXCLUDE rules matching
    assert (
        info.get_script(Path("ignored_module.x90"))
        == PsycloneInfo.RESULT_EXCLUDE
    )

    # 3. Test explicit NO_SCRIPT matching
    assert (
        info.get_script(Path("skipped_module.x90"))
        == PsycloneInfo.RESULT_NO_SCRIPT
    )

    # 4. Test explicit rule script missing physically on storage disk
    with pytest.raises(
        FileNotFoundError, match="Cannot find script '.*missing_script.py'"
    ):
        info.get_script(Path("broken_module.x90"))

    # 5. Test explicit rule script that exists successfully
    valid_script_path = opt_dir / "valid_script.py"
    valid_script_path.touch()
    assert info.get_script(Path("good_module.x90")) == valid_script_path

    # 6. Test file_specific script with non-existent explicit definition
    info_strict = PsycloneInfo(
        name="dsl", base_paths=[base_src], script_root=script_root
    )
    info_strict.update(
        {"script_dir": "psykal", "file_specific": "explicit_custom.x90"}
    )
    with pytest.raises(
        FileNotFoundError, match="Cannot find explicitly requested script"
    ):
        info_strict.get_script(base_src / "explicit_custom.x90")

    # 7. File_specific if the script exists:
    valid_script_path = opt_dir / "file_specific.py"
    valid_script_path.touch()
    assert info.get_script(base_src / "file_specific.x90") == valid_script_path
