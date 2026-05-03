# Mesh Recovery Strategies

The Nikolaisen2022 particle STLs are reconstructed surfaces, and many of them
contain holes, duplicated/non-manifold edges, overlapping facets, or other
triangle topology that exact CAD-style surface meshing does not like. For the
micromagnetic ensemble use case, the package now defaults to a pragmatic
fallback sequence: keep the rough particle shape, produce a closed tetrahedral
volume, and record which strategy was used.

## Default Order

`process_nikolaisen_merrill_meshes` and `process_nikolaisen_size_bin` try these
strategies in order:

1. `gmsh`: use Gmsh Delaunay tetrahedralization on the original STL surface.
2. `pymeshfix_gmsh`: repair the STL with MeshFix, then retry Gmsh.
3. `voxel`: voxelize the rough particle volume and split voxels into
   tetrahedra.
4. `hull_delaunay`: mesh the STL point cloud's Delaunay hull as the final,
   most permissive fallback.

The report CSV includes `mesh_strategy`, `attempted_mesh_strategies`, and
`strategy_errors` so later filtering can distinguish exact meshes from repaired
or brute-force approximations.

## Brute-Force Voxel Fallback

`tetrahedralize_stl_with_voxels` is the workhorse fallback for malformed
surfaces. It rasterizes the particle into voxels, then converts each voxel to
tetrahedra. This tends to remove small holes, non-manifold edge ambiguity, and
surface self-intersection trouble at the cost of stair-stepped geometry.

The requested default spacing is still derived from the 9 nm target edge length.
For very large bounding boxes, spacing is coarsened by
`DEFAULT_BRUTE_FORCE_MAX_CELLS_PER_AXIS` to avoid runaway file sizes. The
realized edge lengths are still reported after conversion, so downstream
Merrill.jl screening can account for the actual mesh size.

## Useful Functions

```python
from stl2fem.conversion import (
    tetrahedralize_repaired_stl_for_merrill,
    tetrahedralize_bruteforce_stl_for_merrill,
    tetrahedralize_stl_with_voxels,
    tetrahedralize_stl_with_delaunay_hull,
)
from stl2fem.workflow import process_nikolaisen_merrill_meshes
```

For a single difficult STL:

```python
tetrahedralize_bruteforce_stl_for_merrill(
    "particle.stl",
    "particle_native.msh",
    "particle_meters.msh",
    input_unit="um",
    target_edge_length_m=9e-9,
    strategy="voxel",
)
```

For the Nikolaisen2022 Plag batch, the default command already uses the full
fallback sequence:

```bash
stl2fem process-nikolaisen-merrill --overwrite
```

To force the bluntest path for every mesh:

```bash
stl2fem process-nikolaisen-merrill \
  --overwrite \
  --mesh-strategy voxel \
  --mesh-strategy hull_delaunay
```

This is appropriate when exact surface fidelity matters less than getting a
complete ensemble of roughly shape-preserving tetrahedral volume meshes.
