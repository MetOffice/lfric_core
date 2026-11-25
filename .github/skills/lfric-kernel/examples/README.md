# LFRic kernel examples and templates

Read this catalogue before scanning application or component kernel sources.
All examples and templates are contained in this directory; this catalogue
does not link to repository source files that may not be present in a checkout.

These are curated starting points, not substitutes for the LFRic API supported
by the active PSyclone installation. Derive a complete code-procedure interface
from matching metadata and a comparable generated PSy-layer call; do not copy
support arguments across kernel families.

## Catalogue

- `cell-column.md`: minimal field copy on one discontinuous function space.
- `dof.md`: scalar field values for one degree of freedom.
- `domain.md`: segment-based domain metadata and interface shape.
- `scalar-arguments.md`: thresholding kernel with scalar real arguments.
- `integer-arguments.md`: kernel with integer field and scalar arguments.
- `logical-arguments.md`: kernel with logical scalar flags controlling branching.
- `mixed-datatypes.md`: kernel combining real and integer fields on same space.
- `cma.md`: columnwise-operator metadata and interface considerations.
- `inter-mesh.md`: coarse/fine field metadata.
- `basis-quadrature.md`: basis, mesh, reference-element, and face quadrature
  metadata.

The cell-column and DOF files contain starter code. The other files contain
focused metadata and interface examples rather than drop-in kernels because
their generated arguments depend strongly on the exact metadata and API
version. No example in this directory has been validated against one known
PSyclone version.
