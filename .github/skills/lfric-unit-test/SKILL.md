---
name: lfric-unit-test
description: "Use when creating, modifying, reviewing, or debugging LFRic pFUnit unit tests (*.pf files), test fixtures, assertions, or kernel unit tests."
---

<!--
-----------------------------------------------------------------------------
Some of the content of this file has been produced with the assistance of
Met Office Github Copilot Enterprise.
-----------------------------------------------------------------------------
-->

# LFRic Unit Tests

## Scope and guidance

- Apply `.github/instructions/lfric-fortran.instructions.md` to test source, including copyright, attribution, declaration grouping, indentation, and continuation rules. Tests are excluded from Doxygen requirements: use ordinary comments for test purpose and non-obvious setup, not new Doxygen blocks or AI documentation-review notes. Preserve existing documentation.
- Use `documentation/source/how_to_use_it/testing/unit_testing.rst` as the repository testing guide. For canned kernel inputs, consult `documentation/source/how_to_use_it/testing/unit_test_canned_data_routines.rst` and the actual helper interfaces before calling them.
- For tests of kernels or their PSyclone interfaces, also load `.github/skills/lfric-kernel/SKILL.md`. Kernel-source-specific formatting rules do not automatically apply to test procedure headers.

## Locate the test

- Start from the procedure and behaviour under test. Inspect its implementation and the nearest existing test before editing; identify one observable expectation and the narrowest check that could disprove it.
- Reuse an existing test module and helpers where suitable. Otherwise, place the test under the owning project's `unit-test` directory, mirroring the relative `source` path and neighbouring module and filename conventions, typically `<module>_test.pf`.
- Prefer `.pf`. Use `.PF` only when preprocessing before pFUnit is required, such as for templated source, and follow the owning project's build conventions. Do not edit generated Fortran test files or generated drivers.

## Test structure

- Prefer standalone `@test` procedures. Give each test a descriptive name tied to one behaviour or failure mode, and make test entry points accessible to the pFUnit-generated driver using the neighbouring module's visibility conventions.
- Align each decorator, including `@test`, `@before`, and `@after`, with its associated procedure header. Module-contained procedure headers and matching terminators remain in column 1 under the general Fortran rules.
- Keep each pFUnit directive on one physical line; do not continue `@assert...` or other `@` directives with `&`. If necessary, compute an intermediate expected value or logical condition in ordinary Fortran before the assertion. Prefer concise directives; do not break their syntax merely to meet line-length guidance.
- Follow the local pFUnit import and supported assertion conventions. Use explicit `only:` imports for ordinary modules; do not copy broad production-module imports from older tests. If the framework's generated code needs additional public symbols or imports, check the installed pFUnit interface and a comparable working test rather than guessing.
- Prefer standalone `@before` and `@after` fixtures. Derive a `TestCase` only when fixture data must be held on the case object; use `MPITestCase` only when the behaviour requires MPI. Do not introduce a fixture class solely to run setup and cleanup procedures.
- Initialise every input and fixture explicitly. Ensure tests can run independently and in any order. Release owned allocations, reset modified global state, and finalise initialised services using the appropriate cleanup path, including after assertion failures. Do not deallocate borrowed pointer targets.

## Inputs and assertions

- Test the smallest useful public behaviour rather than an entire workflow. Keep setup dependencies minimal and avoid unrelated refactoring or new helper abstractions for a single straightforward case.
- Use small, deterministic inputs with distinct values that expose indexing, ordering, sign, and accumulation errors. Avoid only-zero or only-uniform inputs when they could hide defects. Reuse established deterministic random-data helpers only when needed.
- Derive expected results independently: prefer explicit constants or simple calculations from known constants, not a copy of the implementation under test. For regressions, add a case that would fail with the original defect.
- Put expected values before actual values in `@assertEqual`. For floating-point comparisons, choose a justified tolerance of the appropriate kind and scale, following the local pFUnit API. Do not relax tolerances merely to make a test pass; use exact comparisons only when exact equality is part of the contract.
- Assert relevant values, shapes, bounds, state changes, and unchanged inputs or untouched regions. Prefer assertions that report useful expected and actual values; add concise messages where a failure would otherwise be ambiguous.
- Include representative normal cases, valid boundary cases, and supported error paths. Do not invoke an unconditional abort or `stop` in the shared test process to test a failure; use an existing isolated-process harness if the project provides one, otherwise report that coverage limitation.

## Kernel tests

- Call the kernel code procedure directly with primitive scalars and arrays matching its complete generated interface. Do not construct fields, function spaces, meshes, or a PSy layer solely to prepare kernel inputs unless those infrastructure objects are themselves under test.
- Reuse suitable canned-data helpers for dofmaps, basis functions, quadrature, and sizes. Verify dimensions, kinds, function-space assumptions, and layer bounds against the code interface; do not invent helper names or argument lists.
- Where relevant, exercise multiple dofs and layers, nontrivial dofmaps, and nonzero initial output values. Check that read-only inputs and unmapped output entries remain unchanged.
- Distinguish cell-wise contributions from once-per-dof updates. For increment kernels on shared dofs, construct expected results for the actual number of contributions, not an assumed global field addition.

## MPI and resources

- Keep tests serial unless the behaviour depends on communication, rank ownership, or distributed data. Follow neighbouring tests for `@test` rank parameters and MPI fixtures rather than imposing one rank layout on all tests.
- Use the test framework's communicator and initialise and finalise LFRic MPI or halo services consistently with existing fixtures. Do not assume `MPI_COMM_WORLD` is the communicator assigned to the test.
- Keep collective participation consistent across ranks, including assertion and cleanup paths; avoid rank-local early returns that leave other ranks waiting in collectives.
- Reuse project test-data locations and isolated temporary-file conventions. Do not rely on an arbitrary working directory or leave generated files that affect later tests.

## Validation

- Default to non-compiling checks: inspect decorators, visibility, declarations, assertion signatures, input initialisation, expected results, and resource lifetimes. A plain Fortran parser may not accept pFUnit directives; use the configured pFUnit preprocessing path when available and authorised, without treating preprocessing as proof that tests pass.
- Do not automatically compile, build, or run unit-test binaries. When the user explicitly requests executable validation, load `module load lfric`, inspect the owning project's Makefile and test options, and select the narrowest supported test command. Use the existing `unit-tests` target where appropriate; do not invent filtering flags or install tools from the internet.
- Local MPI tests must use at most four ranks. Submit tests requiring more than four ranks as SLURM batch jobs using established site conventions.
- Report checks actually performed, results, and any unverified behaviour. Keep unrelated failing tests outside the edit scope and distinguish static review, preprocessing, compilation, and executed tests.