# stl2fem

Python-first tools for converting STL particle surfaces into tetrahedral FEM
meshes, screening mesh quality, and estimating downstream Merrill.jl
micromagnetic modeling cost.

The package is being developed around the Nikolaisen2022 magnetite inclusion
dataset, but the core functions are generic STL-to-`.msh` helpers built on
PyVista/VTK and Gmsh.

## Units

STL files do not reliably store physical units. They store triangle vertex
coordinates as plain numbers, so a coordinate like `22.88` only means `22.88`
in the coordinate system chosen by the file creator.

For Merrill.jl compatibility, `stl2fem` treats units explicitly:

- declare the input STL coordinate unit with `--input-unit` or
  `input_unit=`;
- define the requested tetrahedron edge length in meters with
  `--target-edge-length-m` or `target_edge_length_m=`;
- mesh in the STL's native coordinate scale for numerical convenience;
- write a separate Merrill-ready `.msh` with coordinates scaled to meters.

The Nikolaisen2022 metadata and coordinate magnitudes are micrometer-scale, so
the examples use `input_unit="um"`. The default target edge length is
`9e-9 m`, i.e. 9 nm for magnetite:

```text
9e-9 m / 1e-6 m per micrometer = 0.009 native STL units
```

Every conversion report prints both the requested physical target and realized
tetrahedron edge lengths, because Gmsh may not produce exactly the requested
edge size.

## Install

```bash
python -m pip install -e ".[examples]"
```

The optional `examples` extra installs notebook/display dependencies. If you
only want the CLI and source package, use:

```bash
python -m pip install -e .
```

## Quick Start

Inventory the Nikolaisen2022 STL files and split them into size bins:

```bash
stl2fem nikolaisen-inventory \
  --dataset-root data/Nikolaisen2022 \
  --out processed/Nikolaisen2022/04_quality_reports/nikolaisen_binary_stl_inventory.csv
```

Process one size bin by checking surface quality, generating Gmsh Delaunay
tetrahedral `.msh` files, checking tetrahedron quality, and estimating Merrill.jl
memory:

```bash
stl2fem process-nikolaisen-bin \
  --dataset-root data/Nikolaisen2022 \
  --output-root processed/Nikolaisen2022 \
  --bin-index 0
```

Generic single-file conversion:

```bash
stl2fem convert-stl input.stl output_native.msh \
  --merrill-msh output_meters.msh \
  --input-unit um \
  --target-edge-length-m 9e-9 \
  --algorithm delaunay
```

By default, conversion targets a physical tetrahedron edge length of `9e-9 m`.
You can override it per run:

```bash
stl2fem convert-stl input.stl output_native.msh --target-edge-length-m 2e-8
```

The reported output includes the realized tetrahedron edge lengths
`edge_length_min`, `edge_length_median`, `edge_length_p95`,
`edge_length_max`, and their `_m` meter-scaled counterparts.

Generated meshes and reports are written under `processed/Nikolaisen2022/`.

## Examples

Notebook examples live in `examples/`:

- `00_basic_functions_demo.ipynb`
- `01_edge_size_realization_demo.ipynb`
- `Nikolaisen2022_00_inventory_and_size_bins.ipynb`
- `Nikolaisen2022_01_bin_00_smallest.ipynb`
- `Nikolaisen2022_02_bin_01_small.ipynb`
- `Nikolaisen2022_03_bin_02_large.ipynb`
- `Nikolaisen2022_04_bin_03_largest.ipynb`

The notebooks are split by sorted STL file size so PyVista rendering output can
be kept manageable. They report every mesh in each bin, but cap display output
with `DISPLAY_LIMIT` by default.

Important dataset note: the Nikolaisen2022 folders named `Binary meshes` can
contain ASCII STL files. `stl2fem` detects the actual file encoding instead of
trusting the folder label.

## Development

See `ENVIRONMENT.md` for reproducible `venv` and `conda` setup commands.
