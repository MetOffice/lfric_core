# Columnwise operator (CMA) example

This metadata example has three operator arguments, three function-space pairs,
and a cell-column iteration entity.

| Code dummy | Access | Function spaces |
| --- | --- | --- |
| `columnwise_matrix_A` | `GH_READ` | `ANY_SPACE_1` to `ANY_SPACE_2` |
| `columnwise_matrix_B` | `GH_READ` | `ANY_SPACE_2` to `ANY_SPACE_3` |
| `columnwise_matrix_C` | `GH_READWRITE` | `ANY_SPACE_1` to `ANY_SPACE_3` |

The metadata entries are `GH_COLUMNWISE_OPERATOR` arguments, and the iteration
entity is `CELL_COLUMN`. This is not a complete CMA code-procedure template.
Support arguments include operator-specific extents and banded-matrix
parameters. Derive their exact order, rank, and intent from the metadata and
generated call.

