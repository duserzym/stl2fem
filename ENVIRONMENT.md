# Environment

The tested local environment uses Python 3.10 and the direct package versions in
`requirements-dev.txt`.

## venv

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e ".[examples,dev]"
```

For exact direct dependency pins:

```bash
python -m pip install -r requirements-dev.txt
python -m pip install -e .
```

## conda/mamba

```bash
mamba env create -f environment.yml
mamba activate stl2fem
```

## Notes

- `pyvista` pulls in VTK, which is the largest dependency.
- STL files do not reliably store physical units. Declare the input unit with
  `--input-unit` or `--input-scale-to-meters`.
- The default physical target tetrahedron edge length is `9e-9 m`. For
  micrometer-scale STLs, that is `0.009` native coordinate units. Use a larger
  `--target-edge-length-m` for quick demos and pilot runs.
- Generated meshes and reports are written to `processed/`, which is ignored by
  git.
