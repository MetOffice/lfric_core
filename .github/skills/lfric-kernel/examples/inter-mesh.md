# Coarse/fine inter-mesh example

This example has three field arguments spanning coarse and fine meshes:

- `fine_field`: `GH_WRITE`, fine mesh, `ANY_DISCONTINUOUS_SPACE_2`
- `coarse_field`: `GH_READ`, coarse mesh, `ANY_SPACE_2`
- `weights`: `GH_READ`, fine mesh, `ANY_DISCONTINUOUS_SPACE_2`

The metadata entries use the corresponding mesh tags:

```fortran
arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2, &
         mesh_arg=GH_FINE) ! fine_field
arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_SPACE_2, &
         mesh_arg=GH_COARSE) ! coarse_field
arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2, &
         mesh_arg=GH_FINE) ! weights
```

This is not a general inter-mesh interface. Prolongation may require additional
scalar arguments, cell maps, and per-mesh dofmaps. Multidata transfers can also
change field extents and weights. Match every mesh tag, function space, and
generated support argument to the exact algorithm and metadata.

