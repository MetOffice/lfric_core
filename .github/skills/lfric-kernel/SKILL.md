---
name: lfric-kernel
description: "Use when creating, modifying, or reviewing LFRic Fortran kernel modules, kernel metadata, or PSyclone kernel interfaces."
---

# LFRic Kernel Writing

## Rule precedence

Apply `.github/instructions/lfric-fortran.instructions.md` as the general Fortran and LFRic baseline. Within repository guidance, an explicit kernel-specific rule in this skill takes precedence over a conflicting general Fortran rule for kernel code. All non-conflicting general rules still apply. This precedence does not override higher-priority instructions or an explicit user request.

## Reference-first workflow

- Before scanning application or component kernel sources, read
  `examples/README.md` and select the closest curated template or reference.
- Treat the skill examples as a fast starting point, not as proof that an
  interface is valid for every PSyclone version or metadata combination.
- If an example does not cover the required feature, inspect only the narrowest
  comparable local kernel and its algorithm or generated PSy-layer call. Use a
  broad source search only when those references leave a specific question
  unanswered.
- For DOF and domain kernels, use the family-specific material rather than
  deriving an interface from a cell-column template. Treat entries labelled as
  interface sketches or references as non-drop-in, and confirm the exact family
  rules in the configured LFRic API.

## Kernel interface

- Place kernel modules under the owning application's or component's `source/kernel` directory, following nearby naming and extension conventions.
- Choose the kernel family and iteration entity before defining the interface. General-purpose cell-column kernels, DOF/owned-DOF kernels, domain kernels, CMA kernels, and inter-grid kernels have different legal arguments and execution rules; check the configured LFRic API documentation before using a specialised family.
- Declare metadata in a public type extending `kernel_type`, and bind the code procedure with `nopass`. Set `operates_on` to the intended iteration entity; do not assume all user kernels operate on cell-columns.
- Use concise Doxygen `@brief` descriptions for kernel modules, metadata types, and code procedures. Add an optional `@details` entry when it provides useful context beyond the summary; when both are present, align their description text as specified by the general Fortran instructions.
- Add an inline comment to every `meta_args` element with the corresponding code dummy name. For continued elements, place it after the continuation marker, for example: `arg_type(...), & ! temperature`.
- Keep `meta_args` entries in algorithm-argument order. There is one entry per algorithm-level scalar, scalar array, field, field vector, or operator; a field vector is one metadata entry but expands to multiple field dummies in the code procedure. Match metadata datatype, access, function space, and stencil/operator properties to the actual arguments and code.
- Use only supported argument-kind and datatype combinations for the selected kernel family. For general-purpose user kernels, scalar and scalar-array arguments are read-only; at least one valid data argument must be modified, and a cell-column kernel needs a field, field vector, or operator to define its domain.
- Pass runtime source data through kernel metadata and code arguments. Do not access runtime values through use-associated module variables; imported variables must be named constants declared with `parameter`.
- Read configuration through accessors in the algorithm layer and pass required values as scalar kernel arguments. Do not pass `config_type` into the kernel.
- Index fields through supplied dofmaps and dimensions, and iterate only over valid dofs and layers. Do not assume one dof per cell or access unavailable stencil entries.
- Keep writes within metadata-declared writable arguments; avoid hidden state changes and I/O in kernels.
- Structure kernel computation to minimise branching, especially repeated `if` tests inside dof and layer loops. Where semantics permit, move loop-invariant decisions outside loops and resolve configuration-dependent choices in the algorithm layer by selecting the appropriate kernel. Prefer simple, equivalent expressions or suitable intrinsics when they remove unnecessary branches without changing numerical behaviour. Retain branches needed for correctness, bounds safety, or exceptional cases; do not replace protective tests with expressions that may evaluate invalid operations or out-of-bounds accesses. Do not assume Fortran logical expressions short-circuit or that `merge` evaluates only its selected argument. Avoid excessive code duplication merely to remove a branch.
- Format counted loops using the general Fortran form `do index=start, end`: put no spaces around `=`, use one space after each comma, and write simple integer upper-bound decrements without spaces around `-` (for example, `nlayers-1`). Indent each nested construct by two spaces and align each `end do` with its opener.

## PSyclone metadata and execution

See the [LFRic API documentation](https://psyclone.readthedocs.io/en/latest/user_guide/lfric.html) and [PSyclone repository](https://github.com/stfc/PSyclone) for the full rules. Use the version configured by this repository as the authority: the `latest` documentation and `main` branch may describe a different API than the installed PSyclone.

- For a general-purpose cell-column kernel, PSyclone determines loop continuity from modified arguments. If multiple modified fields require different iteration spaces, the most demanding space controls the loop. Generic `ANY_SPACE` metadata is conservatively treated as continuous when its actual space is not known; use a discontinuous generic space only when that is guaranteed by the kernel contract. Operator writes and reads have halo limits that also affect legal iteration, so check the API rules before changing operator access.
- Choose field access from both data flow and parallel semantics. `GH_READ` is read-only. `GH_WRITE` requires the first access to be a write; on shared dofs, all iterations must write the same value. `GH_INC` denotes contributions to shared values and requires the appropriate synchronisation or colouring. Use `GH_READINC` when the kernel reads the existing value before incrementing it. `GH_READWRITE` is for supported non-shared updates, not a generic replacement for increment access. PSyclone's allowed modes depend on argument kind and function space.
- `meta_funcs` requests basis and/or differential-basis data for function spaces used by the arguments. Request only what the computation consumes. Keep `gh_shape`, any evaluator targets, quadrature metadata, and algorithm-layer quadrature objects consistent; these determine the shape and availability of generated basis and weight arguments. Use `meta_reference_element` or `meta_mesh` only when the kernel needs those properties and the configured API supports them.
- Specialised kernel constraints are not interchangeable: DOF kernels operate on one dof and have restricted argument/function-space rules; domain kernels receive domain-level data; CMA kernels are cell-column-only and restrict stencils and field vectors; inter-grid kernels require mesh metadata and fields on both meshes. Confirm the exact restrictions in the configured API version, especially before composing these kernels in one invoke.

## Deriving and checking the code interface

- Keep the algorithm, metadata, PSy layer, and code procedure roles distinct. Metadata drives PSyclone's generated call and data-management logic; the code procedure contains only the local computation.
- The code-procedure dummy list includes both user data and generated support arguments. Depending on the metadata and kernel family, support arguments include cell/layer or dof indices, dof counts and maps, stencil sizes and maps, operator extents, basis tables, and quadrature counts or weights. Do not add these generated arguments to `meta_args`.
- Format the bound kernel code procedure's argument list (typically `<name>_code`) across multiple lines to make its correspondence with the kernel type's `meta_args` array clear. Put each dummy representing a metadata argument on its own line, following metadata order; if a metadata entry expands to multiple data dummies, put each expanded data dummy on its own line. Generated support arguments may share a line, for example `ndf, undf, map`. Preserve the exact PSyclone interface order, including support arguments interleaved with data arguments; do not reorder arguments merely to group them. Align continuation arguments under the first argument after `(` and align trailing `&` markers as specified below.
- Format dummy declarations using the general Fortran rules: declare one data entity per statement, keep dummy and local declarations in separate groups, and align `intent`, `dimension` where present, and `::` columns within each contiguous dummy-declaration block. Pad between declaration fields as needed, leave fields empty when an attribute does not apply, and keep exactly one space after `::`. Keep local loop-counter declarations in their own block with the required `! Local loop counters` comment.
- In kernel source, when a logical statement has multiple trailing continuation markers, align every trailing `&` in one shared column within that statement. This applies to procedure headers, metadata constructors, calls, and expressions, and overrides the general Fortran option of using one space before each marker. Set the shared column to the earliest column that fits every marker: find the rightmost code character on marker-bearing lines (excluding existing padding and comments), leave one separator space, then place `&` in the next column. Pad shorter lines only as far as needed to reach that column; do not align markers to an arbitrary distant column such as column 79. Change only padding immediately before the markers; preserve operand alignment, all non-whitespace tokens, and comments after the markers. Retain required separators, including exactly one space in outer array-constructor `, &` separators; if these constraints conflict, realign the constructor elements rather than adding padding after the comma. Verify that marker columns match and that the shared column is the earliest valid one, as well as checking line lengths across the complete statement. For a statement with only one trailing marker, use exactly one space before it.

```fortran
subroutine example_copy_code(result, &
                             source, &
                             ndf, undf, map)
```

- Derive the exact dummy order, rank, extent, kind, and intent using the configured API documentation and a comparable generated PSy-layer call or kernel with the same metadata features. Check field-vector expansion and argument-specific map/basis needs; do not copy an interface across different iteration entities, stencils, operator kinds, quadrature shapes, or PSyclone versions.
- Validate metadata with the repository's PSyclone parsing or generation path when available. Resolve metadata/interface errors before compiling the Fortran kernel; use static and non-compiling validation by default.

## Validation

Use static diagnostics and non-compiling checks by default. Do not compile or build automatically; do so only when explicitly requested by the user. Do not invoke editor diagnostics or validation commands that trigger compilation (including compiler-backed `get_errors`) unless the user explicitly requests a compile or build. Prefer metadata parsing, PSyclone generation checks that do not compile Fortran, and source-level formatting checks. Apply the general Fortran style, spacing, indentation, documentation, and line-length guidance from the linked instruction file unless this skill explicitly states a kernel-specific exception.