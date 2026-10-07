# Gergov et al. (2025) Hekla and Vesuvius magnetite meshes for merrill.jl

## Source

Gergov, H., Muxworthy, A. R., Williams, W., & Cowan, A. (2025). Magnetic recording fidelity of basalts through
3D nanotomography. *Geochemistry, Geophysics, Geosystems*. https://doi.org/10.1029/2024GC011776

- **Data:** Supplementary dataset, Zenodo, https://doi.org/10.5281/zenodo.11369780. Licence CC-BY-4.0.
- **What was used:** `Individual Meshes Hekla.7z` and `Individual Meshes Vesuvius.7z`.
  - These are FIB-SEM slice-and-view reconstructions of (titano)magnetite in the 1991 Hekla and 1944 Vesuvius
    basalts.
  - There is one smoothed, CUBIT-meshed tetrahedral grain per file, in Patran neutral format, in micrometres.
  - Their MD5 checksums matched the Zenodo record on download (2026-10-06).
- **Also used:** `HeklaGrainMetrics.csv` and `VesuviusGrainMetrics.csv`, for the published per-grain metrics.

## Conversion

`scripts/convert_gergov2025.py` (reader: `src/stl2fem/patran.py`).

- **Mesh preserved exactly.** Every node, every tetrahedron and their order are kept. The only change is scaling
  coordinates from µm to m, the stl2fem convention.
- **Output format:** binary Gmsh 2.2 with one physical volume named `particle`. This is the same layout as the
  Nikolaisen2022 meshes, and it loads with merrill.jl `load_tetrahedral_mesh` and with the PINT compact-mesh loader.
- **Checks:** a sample was read back with both gmsh 4.15 and meshio, and node and element counts matched.

```bash
python scripts/convert_gergov2025.py --workers 3      # convert all meshes
python scripts/convert_gergov2025.py --metadata-only  # re-join published metrics only
```

The meshes are written to `data/Gergov2025/merrill_msh/<ID>.msh` (1.94 GB, gitignored like every mesh set in this
repository). `inventory.csv` here pins every output by SHA-256, with the source `.pat` SHA-256 and quality metrics.

## Inventory summary (976 meshes)

**Grains**

| | Value |
|---|---|
| Meshes | 805 Hekla (H###), 171 Vesuvius (V###) |
| Equivalent sphere diameter | 26–624 nm, median 76 nm |
| Grains below 100 nm | 644 |
| Nodes | 418–311,904, median 3,929 |
| Meshes with 5,000 nodes or fewer | 547 |
| Meshes above 100,000 nodes | 15 |

**Mesh quality**

- **Median edge length:**
  - about 5 nm for grains below 300 nm;
  - about 9 nm for grains of 300 nm and above;
  - never above 9.1 nm. The longest single edge is 18 nm.
- **Integrity:** every mesh is a single connected piece with a closed, manifold boundary. There are no inverted,
  degenerate or orphaned elements, and no non-tetrahedral elements. Property ID is 1 throughout.
- **Element shape:** measured by a mean-ratio score (1 = regular tetrahedron).
  - Median of each mesh's 1st percentile: 0.55.
  - 31 meshes contain a sliver below 0.1, and 18 below 0.05 (`shape_min` column).

**Published metrics** (`published_*` columns)

- 948 meshes join a published row.
- **ID fix:** one block of Hekla rows is labelled by STL file name (`828.stl`). Those rows are mapped to the
  matching `H828`-style ID.
- **No published row:** 28 meshes. Their IDs are listed by `--metadata-only`.
- **ESD agreement:**
  - The meshed equivalent sphere diameter agrees with the published one to a median of −0.2%.
  - 26 meshes differ by more than 5%, up to 15%. The largest is H009: 496 nm meshed vs 560 nm published.
  - The published values may refer to a different stage of the surface processing. The `esd_rel_diff_vs_published`
    column records the difference.
- **Published ground states:** 565 SD, 315 SV, 49 BSV, 20 MV.

## Using the meshes

The coordinates are in metres, in the FIB sample frame. Load them the same way as the Nikolaisen2022 `.msh`
files. As with the published study, the micromagnetic material and temperature are set in the solver, not in
the mesh.
