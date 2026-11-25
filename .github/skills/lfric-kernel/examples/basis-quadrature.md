# Basis, mesh, reference-element, and quadrature example

This metadata pattern combines:

| Metadata declaration | Example values |
| --- | --- |
| `meta_funcs` | `W2: GH_BASIS`; `W3: GH_BASIS` |
| `meta_mesh` | `adjacent_face` |
| `meta_reference_element` | `normals_to_horizontal_faces`; `outward_normals_to_faces` |
| `gh_shape` | `GH_QUADRATURE_face` |

This is a metadata example, not a complete kernel. The algorithm layer must
pass a compatible face-quadrature object; the code procedure receives
generated basis tables, normals, and mesh data in metadata-dependent order.
Keep `gh_shape`, the quadrature object, evaluator targets, requested basis data,
and code-procedure dummies consistent.

