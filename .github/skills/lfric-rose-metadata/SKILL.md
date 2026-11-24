---
name: lfric-rose-metadata
description: "Use when creating, modifying, reviewing, or debugging LFRic Rose metadata (rose-meta.conf), namelist variables, metadata descriptions and help, Configurator extensions, imports, or metadata version upgrades."
---

<!--
-----------------------------------------------------------------------------
Some of the content of this file has been produced with the assistance of
Met Office Github Copilot Enterprise.
-----------------------------------------------------------------------------
-->

# LFRic Rose Metadata

## Scope and references

- Rose metadata describes configuration settings, GUI presentation, and validation. LFRic extends it to generate Fortran namelist readers, configuration types, accessors, enumerations, and test configuration helpers. Keep metadata definitions distinct from runtime configuration values and from PSyclone kernel metadata.
- Consult `documentation/source/how_it_works/build_system/configuration_meta.rst` for LFRic extensions and `documentation/source/how_to_use_it/configuration/using_configuration.rst` for the generated configuration API. Consult the Rose documentation at https://metomi.github.io/rose/doc/html/ for standard metadata properties.
- For the extraction and generation pipeline, consult `documentation/source/how_it_works/build_system/picker.rst` and `documentation/source/how_it_works/build_system/configurator.rst`. For configuration-dependent tests, consult the "Feigning Configuration" section of `documentation/source/how_to_use_it/testing/unit_testing.rst`. Treat documented filenames, command options, and legacy examples as guidance to verify against the current tools, not templates to copy uncritically.
- Check the installed Rose and Configurator versions and nearby working metadata when syntax is uncertain. Repository documentation and examples can differ, especially for expressions and ignored settings; do not assume a documented extension is accepted by every tool version.
- When changing generated API consumers, apply `.github/instructions/lfric-fortran.instructions.md` to Fortran edits and load the unit-test or kernel skill when relevant. Do not edit generated configuration source or generated metadata JSON directly.

## Locate the owning metadata

- Start from the requested setting, consuming configuration accessor, or owning namelist. Read its existing metadata section and one nearby example before editing. Establish the datatype, units, meaning, valid values, and whether the value is required; do not invent these from its name.
- Metadata normally lives under `<owner>/rose-meta/lfric-<name>/HEAD/rose-meta.conf`, alongside released version directories. Follow the actual ownership and import chain, not an assumed location. Top-level `rose-meta` entries may expose metadata held elsewhere; avoid creating a second independent definition.
- Follow `import=` entries and confirm the application selects the modified metadata version. Change the owning definition instead of introducing an accidental override. Preserve unrelated imports, property ordering, comments, and formatting.
- Modify development metadata under `HEAD` unless the task explicitly targets a release. Do not rewrite historical version directories merely to keep them identical. For compatibility changes, inspect neighbouring upgrade macros and versioning conventions before deciding what migration is required.

## Syntax and LFRic extensions

- Use `[namelist:<name>]` for a namelist and `[namelist:<name>=<member>]` for a member. Prefer lowercase identifiers consistently across metadata, runtime values, and generated accessor usage; do not rename existing settings as unrelated cleanup.
- Use standard Rose properties such as `type`, `compulsory`, `description`, `help`, `values`, `value-titles`, `range`, `fail-if`, `warn-if`, and `trigger` according to their supported syntax. Follow local `title`, `ns`, and `sort-key` conventions for GUI placement. Do not add properties without a clear purpose.
- The long-term project policy is `compulsory=true` for input metadata entries, except derived inputs that Rose must not emit into namelists. Apply this policy to new or changed inputs while assessing compatibility; do not make unrelated existing inputs compulsory as cleanup. A derived member still needs metadata so the Configurator can calculate it at runtime.
- Normally give every namelist member a `type`. For LFRic enumerations supported directly in Fortran, use `values` and `!enumeration` without a standard `type`: Rose overrides the type for enumerated choices and may report a metadata validation error if it is present. For readable enumerated choices backed by legacy integer interfaces without LFRic enumeration support, use `!type=integer` as required by the Configurator rather than `type=integer`; verify the mapping against the consuming code.
- Preserve meaningful exclamation marks. LFRic uses `!`-prefixed properties as generator extensions even though Rose treats them as comments. An ignored section such as `[!namelist:...]` has a different role; do not remove its marker as comment cleanup. Use ordinary `#` comments for prose following local conventions.
- `!kind` is optional and defaults to `default`. The documented default maps to `i_def`, `r_def`, or `l_def` according to datatype, not universally to `r_def`; check supported alternatives against the generator and `constants_mod`. Use an explicit property when needed for clarity or a non-default kind. The Fortran preference for omitting `kind=` in datatype declarations does not apply to metadata `!kind=` properties.
- `!string_length` is optional and defaults to `default`, mapping to `str_def`; use `!string_length=filename` for `str_max_filename` where appropriate. Character length is not a character kind. Preserve the quoting required for character values in Rose and Fortran namelists.
- Prefer meaningful enumerations over magic integer choices. For LFRic enumerated choices, use `values` and `!enumeration=true`, and provide matching `value-titles` where possible so GUI labels avoid raw quotes, underscores, and unexplained abbreviations. Keep titles in the same order as values. Adopting enumerations may require changes in the consuming Fortran or its interfaces; verify generated constant names and use them rather than guessed numeric codes or strings.
- Distinguish selectable GUI choices from an LFRic enumeration: `values` and `value-titles` do not themselves require conversion to generated integer constants. For unrestricted character input, use `type=character` without `!enumeration`. For restricted character choices, use `values`, matching `value-titles`, and `!type=character` rather than a standard `type` property, retaining character data in the generated interface. Specify `!string_length` where needed, using its documented default otherwise.
- Omit `!enumeration` when enumeration is false; do not add redundant `!enumeration=false`. Only add `!enumeration=true` when the consuming Fortran actually requires an LFRic integer enumeration. Do not silently change a character setting or logical switch into an integer enumeration merely to provide point-and-click choices; update the complete suggested metadata block to reflect the requested representation.
- Do not promise radio buttons solely from `values` and `value-titles`: they define the selectable options, while the actual radio-button or dropdown presentation depends on the Rose editor and its widget support. Verify any explicit widget property or class against the installed version rather than inventing a widget identifier.
- For arrays, distinguish Rose `length` validation from generated Fortran `!bounds`. If a dynamic array's size depends on another namelist input, use `length=:` with `!bounds` referencing that size and suitable `fail-if` or `warn-if` checks for the actual length. Dynamic logical arrays require `!bounds`; do not use `length=:` and `type=logical` alone. Check fixed versus variable extents and the generated accessor interface.
- Every namelist definition with `duplicate=true` must also contain a non-empty `!instance_key_member=<member>` property in that same namelist section. This is mandatory, not optional; flag a missing or empty instance key during writing and review. Check that the named member is defined in the namelist metadata, each runtime instance supplies a unique key value, and consumers select the intended instance.
- `!instance_key_member` must always name a member whose namelist input datatype is string/character. Do not use integer, real, or logical members as instance keys; flag non-character keys during writing and review. Use `type=character` for an unrestricted string or the appropriate generator-visible character type for restricted string choices. Verify the effective input datatype rather than inferring it from a GUI label or a generated enumeration accessor.
- Documented non-default kind choices include integer `short`, `medium`, and `long`; real `native`, `single`, `double`, and `second`; and logical `native`. These map to the corresponding `i_*`, `r_*`, or `l_*` constants. Do not infer that a choice supported for one datatype is valid for another; check generator support and the exported constant before changing precision.
- When a value can be derived wholly from namelist inputs and source-code parameters, represent it as a derived input using `!expression` rather than an independently configurable value. Verify exact expression syntax, references, constants, and kind suffixes against a working example and the generator. Provide `description` and `help` even though derived inputs are not visible in the Rose configuration editor and users may need to inspect the metadata directly.

Example of restricted character choices, not an LFRic integer enumeration:

```ini
[namelist:bert=sausage]
compulsory=true
description=Sausage variety
!type=character
!string_length=default
sort-key=Panel-A02
values='small', 'jumbo', 'battered', 'curry'
value-titles=Small, Jumbo, Battered, Curry
```

## Extraction and generated interfaces

- `rose_picker` converts application metadata into JSON and a namelist-name list consumed by the Configurator. Supply each required metadata search root with a separate `-include_dirs` option; these roots contain `rose-meta` packages, not just individual imported files. Use `-directory` to isolate output, and verify exact options with local `-help`. Existing files use `import=` even though older picker documentation refers to an "include" directive; preserve the syntax accepted by the current metadata loader.
- `GenerateNamelistLoader` generates per-namelist reading and MPI-broadcast procedures. `GenerateConfigLoader` coordinates loading and storing namelist objects in the application's configuration object. Keep their recognised namelist set consistent; adding metadata without including it in the application's resolved set will not expose it automatically.
- `GenerateExtendedNamelistType` generates the namelist-specific `<name>_nml_type` and read-only accessors. `GenerateConfigType` generates the application-specific container. Duplicate-namelist flags must be propagated wherever the owning generation pipeline requires them; use that pipeline rather than running a partial hand-written sequence with inconsistent options.
- `GenerateFeigns` generates configuration setup routines for tests. These are not a production configuration-loading interface. Do not replace application loading with feign calls or bypass the generated MPI distribution with rank-local configuration assignments.

## Loading and consuming configuration

- Initialise each `config_type` object and load it using the generated `read_configuration` interface, following the owning driver. The driver configuration component can orchestrate loading and checking that required namelists are present; do not assume GUI validation alone ensures application-level completeness.
- Access ordinary configuration values through the relevant configuration object, typically `modeldb%config` in applications, using generated namelist accessors such as `config%base_mesh%mesh_name()`. Keep multiple loaded configurations independent; do not introduce module-global state that silently selects one configuration for all consumers.
- Import named enumeration parameters from the appropriate `<namelist>_config_mod` when needed, but do not import mutable namelist values directly from that module. The documentation discourages this legacy access because it cannot distinguish multiple configuration objects or repeated namelist instances. Do not assume a documented migration deadline proves that old access has already been removed in the checked-out version.
- Enumeration settings are represented by `integer(i_def)` choices with generated member-prefixed constants, for example `geometry_planar`. An unset choice is not automatically one of the valid enumerators; retain appropriate requiredness checks or invalid-choice handling instead of assuming the first option is a default.
- A repeated namelist is a collection, not a single ordinary namelist component. Use its generated iterator or supported lookup API and identify the intended instance from the configured key/profile. Verify generated iterator type names, `get_profile_name` or other key accessors, ownership, and lifetimes against current generated code. Do not copy direct single-instance access into duplicate-namelist consumers or assume every instance-key representation is itself an enumeration integer.

## Configuration-dependent tests

- Explicitly initialise every configuration value reached by the unit under test, including values read by helper procedures. Do not rely on defaults, a previous test, or partial setup to initialise omitted configuration.
- Where configuration feigning is needed, use generated `feign_config_mod` routines with named arguments. Verify the current signatures and how the routines populate the configuration used by the test; older documentation describes module-global interfaces that may differ from the current configuration-object API.
- Reset configuration through the matching fixture cleanup and avoid leaking state between tests. Load `.github/skills/lfric-unit-test/SKILL.md` for fixture and assertion rules; do not add configuration fixtures to primitive-array kernel tests that do not need them.

## Relevance and validation constraints

- Use `trigger` wherever possible so users see only settings relevant to their configuration. Trigger-ignored settings are not emitted into generated namelist input files. Keep variables that trigger a namelist outside that namelist, and put each trigger condition on the metadata entry supplying that condition, not on unrelated entries.
- When the same condition controls every member of a namelist, target the whole namelist in the controlling entry's `trigger` rather than listing every member separately. If the whole namelist can be trigger-ignored so that no active instances remain, every Rose file-output `source` entry referencing it must mark that source as optional using parentheses, for example `(namelist:bert)`. Apply this to application configurations and their applicable optional configurations or templates; preserve any existing instance selector inside the parentheses for repeated namelists. This is file-output optionality, not a reason to set individual members to `compulsory=false`, and it does not replace checking whether application code requires the namelist.
- Conditions targeting the same setting combine with logical AND. For example, a geometry choice and a profile choice can each trigger the same specialised namelist; both conditions must hold. Check combined conditions for accidentally unreachable inputs.
- Use `fail-if` for invalid combinations that triggers cannot express, such as array-length consistency, choices whose validity depends on another input, or consistency between settings that must remain active. Use `warn-if` for advisory conditions, not hard errors.
- Add suitable `range` bounds even when they seem obvious: they provide immediate GUI feedback. Use `range` for sensible permitted bounds, not a narrower recommendation; put recommended settings in `help`.

For a `bert` namelist that can be wholly trigger-ignored, a file-output entry in `rose-app.conf` can be:

```ini
[file:configuration.nml]
mode=auto
source=(namelist:bert)
```

The parentheses allow `bert` to be absent without making this file source a missing-required-source error. When active, the namelist is included normally. In a source list, wrap only those sources that may be absent, leaving required sources unparenthesised.

## GUI organisation

- Use `ns` to organise an intuitive GUI hierarchy and group related inputs. A setting can appear in a different panel from its namelist, but do not use this to disguise poor namelist ownership; place settings in the appropriate namelist first.
- Follow `sort-key=Section-A00`-style keys for namespace sections in the GUI tree and `sort-key=Panel-A00`-style keys for variables within a panel, adjusting indices to neighbouring entries.
- Use `title` as the user-facing panel or tree label, independently of `ns`. Mutually exclusive namespaces presenting alternative inputs for the same conceptual section should normally share a title and section sort position.
- Put common context, terminology, and acronyms in namespace or panel help where that avoids repetition. Consider a separate advanced panel for expert-only settings, but keep its explanations understandable to non-specialists.

## Descriptions and help

- Write `description` as a short, clear summary displayed beneath the code variable's name in the GUI. State what the setting represents or does; do not repeat the identifier, restate its datatype, or use vague wording such as "amount" without context. Include units or dimensions in square brackets, for example `[kg]`, and scope or reference frame where needed. For logical settings, use positive logic: describe what happens when true. Keep it concise; if necessary, wrap over two lines using Rose continuation syntax to avoid excessive editor-panel whitespace. Put derivations, caveats, and references in `help`, not `description`.
- Use British English. Explain domain-specific abbreviations where needed and avoid implementation details that are irrelevant to configuration users. The Doxygen-specific wording and alignment rules do not automatically apply to Rose prose.
- Write `help` so a general user can decide whether the setting matters to their case, make a sensible choice when it does, and identify a more detailed reference or expert when necessary. Do not assume knowledge of the science or require users to read papers merely to discover a setting is irrelevant.
- Give brief context, expand acronyms and jargon, explain units and expressions, describe each enumeration choice where applicable, and provide justified recommended settings. Avoid repeating `description` or restating the datatype. Do not use LaTeX, which Rose help does not render, or add hyperlinks to other pages that create link-maintenance overhead; use descriptive references or expert contacts where needed.
- Format `help` with readable whitespace. Follow Rose continuation syntax: subsequent value lines are indented and begin with `=`; use a continuation containing only `=` for a paragraph break. Move shared explanations into panel or namespace help where appropriate.
- Do not replace triggers with prose lists of valid and invalid configurations. Encode relevance with `trigger`, permitted values with `range` or `values`, and remaining consistency requirements with `fail-if` or `warn-if`; explain necessary context without duplicating those lists in help. Do not claim a runtime default, unit conversion, or optional behaviour unless the application implements it. GUI defaults and metadata requiredness do not by themselves establish runtime fallback values.

Example for a required positive integer in an existing solver namelist:

```ini
[namelist:solver=max_iterations]
compulsory=true
description=Maximum number of solver iterations
help=Upper limit on the number of iterations performed by the solver.
    =The solve may finish earlier when its convergence criterion is met.
!kind=default
range=1:
type=integer
```

Use this example only where the solver contract actually has those semantics.

## Adding or changing a setting

- Confirm whether the namelist already exists, add or update the member definition, and keep metadata imports and application selection consistent. Ask for unresolved datatype, units, requiredness, or compatibility requirements rather than guessing.
- Update affected runtime configurations, examples, and test fixtures with suitable values, particularly when adding a compulsory member. Do not assume changing `description`, `help`, or `compulsory` inserts a value into existing configurations.
- For renamed, removed, newly required, or semantically changed settings, assess version upgrades and migration macros. Preserve compatibility where required and test the upgrade using existing project procedures; do not invent a macro name or version transition.
- Use the normal configuration-generation path before relying on new accessors or enumeration constants. Inspect `lfric_build/configurator.py` and the owning build setup for how `rose_picker` resolves imports and supplies metadata to the generators. Read configuration through the generated API; for kernels, pass values through algorithm-layer arguments rather than passing a configuration object into the kernel.

## Validation

- Check section identifiers, duplicate definitions, import resolution, kinds, bounds, enumeration choices, units, prose, required runtime values, and affected consumers. Keep validation scoped to the changed metadata and its direct dependencies.
- Use installed Rose metadata and configuration validators and the repository's non-compiling extraction or generation tools where available. Inspect their local help and established commands first; do not invent flags or treat generic INI parsing as validation of Rose imports, ignored states, or LFRic extensions.
- Check `fail-if` and `warn-if` explicitly: the Rose editor may require selecting its `Metadata-Check fail-if,warn-if` check rather than presenting these checks automatically. Use `validate_rose_meta` where available and verify its supported invocation and coverage. Derived inputs are not covered by the Rose validation script, so independently inspect or test their generated expressions and dependencies; passing ordinary metadata validation is not sufficient for them.
- Use the local environment, loading `module load lfric` or the appropriate installed tooling module as needed. Do not fetch untrusted tools from the internet. Run generators into an isolated output directory rather than overwriting source or existing build artefacts unnecessarily.
- Do not automatically compile or build Fortran. Validate extraction and generated interfaces without compilation where possible; compile or execute tests only when explicitly requested. For pFUnit changes, follow the unit-test skill.
- Report checks actually performed and distinguish metadata syntax validation, import resolution, generated-code inspection, migration tests, and executable testing. State unavailable tools and unverified requirements explicitly.