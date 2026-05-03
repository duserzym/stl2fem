# Nikolaisen2022 Examples

These notebooks use `data/Nikolaisen2022/Plag Binary meshes/*.stl` as the
production example input folder. The package detects the real STL encoding,
because the dataset's binary folders can contain ASCII STL files.

The conversion notebooks default to the production target:

- phase: Plag only
- metadata filter: `EVSD < 1 um`
- output folder: `data/Nikolaisen2022_merrill_msh`
- retained mesh files: Merrill-ready meter-scale `.msh` files only
- temporary native-unit diagnostic meshes are deleted after quality metrics are
  recorded
- per-mesh timeout: package default, with failures recorded in the report CSV
  so one difficult STL does not block the rest of a bin
- recovery strategy: original Gmsh, MeshFix+Gmsh, voxel tetrahedra, then
  Delaunay hull fallback when rough shape preservation is preferable to losing
  a particle from the ensemble

Run order:

1. `00_basic_functions_demo.ipynb`
2. `01_edge_size_realization_demo.ipynb`
3. `Nikolaisen2022_00_inventory_and_size_bins.ipynb`
4. `Nikolaisen2022_01_bin_00_smallest.ipynb`
5. `Nikolaisen2022_02_bin_01_small.ipynb`
6. `Nikolaisen2022_03_bin_02_large.ipynb`
7. `Nikolaisen2022_04_bin_03_largest.ipynb`

The processing notebooks are split by sorted STL file size to keep PyVista
display output manageable. They print/report every mesh in their bin, but the
default display limit is intentionally small so executed notebook files do not
become enormous. Increase `DISPLAY_LIMIT` inside a notebook for deeper visual
inspection.

The inventory notebook also estimates volume-equivalent grain size in nm for
each STL and overlays reference cube sizes with 10, 50, 100, and 200 nm edge
lengths.

STL files do not reliably store units. The examples treat the Nikolaisen2022
coordinates as micrometers by dataset convention and write Merrill-ready meshes
in meters. The default physical target tetrahedron edge length is `9e-9 m`,
which is `0.009` native units for micrometer-scale coordinates. Each conversion
report prints realized tetrahedron edge-length statistics in both native units
and meters because Gmsh can realize a different size distribution than
requested.

The recovery strategy used for each particle is recorded in the report CSV.
See `../MESH_RECOVERY_STRATEGIES.md` for the source functions and tradeoffs.
