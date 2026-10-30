###############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
###############################################################################
# Some of the content of this file has been produced with the assistance of
# Met Office Github Copilot Enterprise.

"""Utility module for parsing and validating field output from Fortran programs.
Provides helpers for comparing expected field values (per-rank basis) against
actual output from Fortran test executables. Supports creating expected data
easily and validates that each field appears exactly once in the output.

The Fortran side (test_helpers_mod) prints one line per rank in the form:
    field_name (rank_id) = value1 value2 ... valueN
with values in ES format at full precision, so exact round-tripping of 64-bit
values is possible. 32-bit fields need a looser tolerance (around 1e-6
relative to the magnitude of the values).

Usage:
    # Create expected field data easily
    expected = uniform_field_data(6, [1.0]*27)
    # Or with different values per rank
    expected = per_rank_field_data({0: [1.0]*27, 1: [2.0]*27})
    # Parse and compare in one step
    extract_field_and_compare('field_name', expected, stdout, stderr)
"""
from math import isclose, isfinite
from re import compile as re_compile, escape as re_escape
from typing import Dict, List, Mapping, Sequence
from testframework import TestFailed


def _to_float_list(values: Sequence[float], context: str) -> List[float]:
    """Validate that values is a list/tuple of real numbers and copy it.
    Args:
        values: Candidate sequence of numbers
        context: Description used in error messages
    Returns:
        New list of floats
    Raises:
        ValueError: If values is not a list/tuple of int/float
    """
    # Only accept a list or tuple. Strings, sets, generators and numpy arrays
    # are rejected so mistakes are caught early rather than in the comparison.
    if not isinstance(values, (list, tuple)):
        raise ValueError(f"{context} must be a list or tuple")
    result: List[float] = []
    for idx, value in enumerate(values):
        # bool is a subclass of int in Python, so True/False would otherwise
        # be accepted as 1/0. Reject them, along with anything non-numeric
        # (e.g. the string '1.0'), reporting the offending index.
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(
                f"{context} index {idx} must be a real number, got {value!r}")
        # Store as float so all later comparisons work on the same type.
        result.append(float(value))
    # A new list is returned, so the caller's data is never shared or mutated.
    return result


def _is_rank(value: object) -> bool:
    """Return True if value is a valid (non-negative, non-bool) integer."""
    # Exclude bool explicitly because isinstance(True, int) is True.
    return (isinstance(value, int) and not isinstance(value, bool)
            and value >= 0)


def uniform_field_data(num_ranks: int,
                       values: Sequence[float]) -> Dict[int, List[float]]:
    """Create expected field data with the same values on every rank.
    Example:
        uniform_field_data(6, [1.0] * 27)
    Args:
        num_ranks: Number of MPI ranks (must be at least 1)
        values: Field values expected on each rank
    Returns:
        Dictionary mapping rank ID to list of field values
    Raises:
        ValueError: If num_ranks is not a positive integer or values is not
                    a list or tuple of real numbers
    """
    # At least one rank is required: zero ranks would produce empty expected
    # data, which would make any comparison pass without checking anything.
    if not _is_rank(num_ranks) or num_ranks < 1:
        raise ValueError("num_ranks must be a positive integer")
    # Validate and convert the values once, before copying them per rank.
    checked = _to_float_list(values, "values")
    # Give each rank its own copy of the list, so changing one rank's
    # expected values cannot affect any other rank.
    return {rank: list(checked) for rank in range(num_ranks)}


def per_rank_field_data(
        data: Mapping[int, Sequence[float]]) -> Dict[int, List[float]]:
    """Create expected field data with different values on each rank.
    Example:
        per_rank_field_data({0: [1.0, 2.0], 1: [3.0, 4.0]})
    Args:
        data: Non-empty mapping from rank ID to the values on that rank
    Returns:
        New dictionary mapping rank ID to list of field values
    Raises:
        ValueError: If data is not a non-empty mapping, a key is not a
                    non-negative integer rank, or a value is not a list or
                    tuple of real numbers
    """
    # Accept any mapping type (dict, OrderedDict, ...), but nothing else.
    if not isinstance(data, Mapping):
        raise ValueError("data must be a mapping of rank ID to values")
    # Empty expected data would make the comparison pass vacuously.
    if not data:
        raise ValueError("data must contain at least one rank")
    # Build a new dictionary rather than returning the caller's object, so
    # later changes to the input cannot change the expected data.
    result: Dict[int, List[float]] = {}
    for rank, values in data.items():
        # Rank IDs must be real MPI rank numbers: integers 0, 1, 2, ...
        if not _is_rank(rank):
            raise ValueError(
                f"Rank keys must be non-negative integers, got {rank!r}")
        # Validate and copy this rank's values (rank is included in errors).
        result[rank] = _to_float_list(values, f"Values for rank {rank}")
    return result


def parse_field_output(output_text: str,
                       field_name: str) -> Dict[int, List[float]]:
    """Parse field output from Fortran program output.
    Extracts field values for each MPI rank from lines formatted as:
    "field_name (rank_id) = value1 value2 value3 ..."
    The field name must match exactly: it may not be preceded by a word
    character, so 'after_caching' will not match 'field_after_caching'.
    Args:
        output_text: Complete output from the Fortran executable
        field_name: Name of the field to extract (e.g., 'field_cache_field')
    Returns:
        Dictionary mapping rank ID to list of field values (empty if the
        field does not appear)
    Raises:
        ValueError: If the field appears more than once for the same rank,
                    or its values cannot be parsed as numbers
    """
    # Build the pattern for one line of field output:
    #   (?<!\w)            - field name must not follow a letter, digit or _,
    #                        so it cannot match the end of a longer name
    #   re_escape(...)     - treat every character in the name literally,
    #                        even regex symbols such as '.', '+', '(' or '['
    #   \s*\((\d+)\)       - the rank number in brackets (captured as group 1)
    #   \s*=\s*(.*)        - '=' followed by the values (captured as group 2)
    pattern = re_compile(
        rf'(?<!\w){re_escape(field_name)}\s*\((\d+)\)\s*=\s*(.*)')
    result: Dict[int, List[float]] = {}
    # Check every line of the output; lines that don't match are ignored
    # (logging, other fields, etc.).
    for line in output_text.splitlines():
        match = pattern.search(line)
        if not match:
            continue
        rank = int(match.group(1))
        values_str = match.group(2).strip()
        # Each field should be printed once per rank. A second line for the
        # same rank means duplicate output, which could hide a real
        # problem, so treat it as an error rather than overwrite the data.
        if rank in result:
            raise ValueError(
                f"Field '{field_name}' appears multiple times for rank "
                f"{rank}. This indicates duplicate output or an ambiguous "
                f"field name.")
        # Values are separated by spaces. float() also accepts 'NaN' and
        # 'Infinity', which are checked later in assert_field_output. Anything
        # else that isn't a number (e.g. Fortran's '*******' when a value is
        # too wide for its format) raises an error here.
        try:
            result[rank] = [float(v) for v in values_str.split()]
        except ValueError as err:
            # Re-raise with the field and rank included; 'from err' keeps the
            # original error in the traceback.
            raise ValueError(
                f"Could not parse field values for '{field_name}' rank "
                f"{rank}. Values string: '{values_str}'. Error: {err}"
            ) from err
    # An empty result means the field never appeared; the caller decides
    # whether that is an error (assert_field_output reports missing ranks).
    return result


def assert_field_output(field_name: str,
                        expected: Mapping[int, Sequence[float]],
                        actual: Mapping[int, Sequence[float]],
                        stdout: str,
                        stderr: str,
                        tolerance: float = 1e-10) -> None:
    """Assert that actual field output matches expected values.
    Values are compared with an absolute tolerance. NaN never matches
    anything; infinities only match an identical infinity.
    Args:
        field_name: Name of the field being tested (for error messages)
        expected: Expected values per rank (must not be empty)
        actual: Actual values per rank (from parse_field_output)
        stdout: Standard output from test (for error context)
        stderr: Standard error from test (for error context)
        tolerance: Absolute comparison tolerance (default 1e-10)
    Raises:
        TestFailed: If validation fails, with detailed message and context
    """
    # With no expected data every check below would pass without testing
    # anything, so treat it as a failure instead.
    if not expected:
        raise TestFailed(
            f"No expected data supplied for field '{field_name}'",
            stdout=stdout, stderr=stderr)
    # Step 1: check that the output contains exactly the expected ranks.
    expected_ranks = set(expected.keys())
    actual_ranks = set(actual.keys())
    # Ranks that should have printed the field but didn't. This also covers
    # the field being missing from the output altogether.
    missing_ranks = expected_ranks - actual_ranks
    if missing_ranks:
        message = (f"Field '{field_name}' missing data for ranks: "
                   f"{sorted(missing_ranks)}. "
                   f"Expected ranks: {sorted(expected_ranks)}, "
                   f"Actual ranks: {sorted(actual_ranks)}")
        raise TestFailed(message, stdout=stdout, stderr=stderr)
    # Ranks that printed the field but were not expected, e.g. the test was
    # run with more MPI processes than the expected data describes.
    extra_ranks = actual_ranks - expected_ranks
    if extra_ranks:
        message = (f"Field '{field_name}' has unexpected data for ranks: "
                   f"{sorted(extra_ranks)}. "
                   f"Expected ranks: {sorted(expected_ranks)}, "
                   f"Actual ranks: {sorted(actual_ranks)}")
        raise TestFailed(message, stdout=stdout, stderr=stderr)
    # Step 2: compare values rank by rank. Sorting the ranks means the
    # lowest failing rank is always the one reported.
    for rank in sorted(expected_ranks):
        expected_values = expected[rank]
        actual_values = actual[rank]
        # Check the counts first: zip() below stops at the shorter list, so
        # a missing or extra value would otherwise go unnoticed.
        if len(expected_values) != len(actual_values):
            message = (f"Field '{field_name}' rank {rank}: "
                       f"Expected {len(expected_values)} values but got "
                       f"{len(actual_values)}. "
                       f"Expected: {list(expected_values)}, "
                       f"Actual: {list(actual_values)}")
            raise TestFailed(message, stdout=stdout, stderr=stderr)
        # Compare each value in order, reporting the first mismatch.
        for idx, (expected_val, actual_val) in enumerate(
                zip(expected_values, actual_values)):
            # NaN always fails (NaN != NaN, even if NaN was expected), and an
            # infinity fails unless the same infinity was expected. This check
            # is separate so the message clearly says a non-finite value was
            # produced, rather than giving a confusing difference such as
            # 'nan' or 'inf'.
            if not isfinite(actual_val) and actual_val != expected_val:
                message = (f"Field '{field_name}' rank {rank} index {idx}: "
                           f"Expected {expected_val} but got non-finite "
                           f"value {actual_val}")
                raise TestFailed(message, stdout=stdout, stderr=stderr)
            # Absolute-tolerance comparison (rel_tol=0 turns off the relative
            # check). Unlike 'abs(a - b) > tol', isclose() returns False for
            # NaN, so a NaN can never pass by accident.
            if not isclose(actual_val, expected_val,
                           rel_tol=0.0, abs_tol=tolerance):
                message = (f"Field '{field_name}' rank {rank} index {idx}: "
                           f"Expected {expected_val} but got {actual_val} "
                           f"(difference: {abs(expected_val - actual_val)}, "
                           f"tolerance: {tolerance})")
                raise TestFailed(message, stdout=stdout, stderr=stderr)


def extract_field_and_compare(field_name: str,
                              expected_data: Mapping[int, Sequence[float]],
                              stdout: str,
                              stderr: str,
                              tolerance: float = 1e-10) -> None:
    """Parse a field from stdout and compare it with expected values.
    Equivalent to parse_field_output followed by assert_field_output, except
    that parse errors are reported as TestFailed with output context.
    Args:
        field_name: Name of the field to extract and validate
        expected_data: Expected values per rank
        stdout: Standard output from the test executable (parsed for data)
        stderr: Standard error from the test executable (for error context)
        tolerance: Absolute comparison tolerance
    Raises:
        TestFailed: If parsing or validation fails
    """
    # Parse the field from stdout. Parse problems (duplicates, values that
    # aren't numbers) are re-raised as TestFailed so the test framework
    # reports them as a test failure, with the program output attached,
    # rather than as a crash of the test script.
    try:
        actual = parse_field_output(stdout, field_name)
    except ValueError as err:
        raise TestFailed(str(err), stdout=stdout, stderr=stderr) from err
    # Compare the parsed values against what was expected.
    assert_field_output(field_name, expected_data, actual, stdout, stderr,
                        tolerance)
