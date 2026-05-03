# Nikolaisen2022 Examples

These notebooks use `data/Nikolaisen2022/* Binary meshes/*.stl` as the example
input folders. The package detects the real STL encoding, because the dataset's
binary folders can contain ASCII STL files.

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

STL files do not reliably store units. The examples treat the Nikolaisen2022
coordinates as micrometers by dataset convention and write Merrill-ready meshes
in meters. The default physical target tetrahedron edge length is `9e-9 m`,
which is `0.009` native units for micrometer-scale coordinates. Each conversion
report prints realized tetrahedron edge-length statistics in both native units
and meters because Gmsh can realize a different size distribution than
requested.
