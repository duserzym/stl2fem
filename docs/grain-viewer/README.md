# Magnetite grain viewer

A static three.js page for reviewing every plagioclase- (PLAG001–246) and
orthopyroxene-hosted (OPX001–082) magnetite grain of the Nikolaisen et al.
(2022) FIB-SEM dataset. For each grain it shows:

- **Original**: the published per-particle STL;
- **Processed**: the boundary surface of the tetrahedral mesh that the
  micromagnetic (MERRILL) campaigns use.

The page has four view modes: original, processed, overlay and side by side.
The camera is orthographic, so the scale bar is exact anywhere in the view.
The page also shows automatic geometry/topology checks, provenance hashes and a
per-grain review log. Review notes are kept in the reviewer's browser and are
exported or imported as CSV.

## Rebuilding the data

```bash
python scripts/build_grain_viewer.py
```

The script needs only numpy. It reads:

- the raw STLs from `data/Nikolaisen2022/`;
- the meshes from the `data/Nikolaisen2022*_merrill_msh*` directories;
- the per-grain mesh selection from the PINT_with_reversal inventories:
  - `--batches`, default `../PINT_campaign_opx/.../production_batches`;
  - files `nikolaisen2022_mesh_inventory_245_surface_fill.csv` and
    `nikolaisen2022_opx_mesh_inventory_82_9nm_v2.csv`.

It verifies every mesh against the SHA-256 in its inventory. It then writes
`data/index.json` plus one quantized, gzipped surface per grain and variant.
Use `--only PLAG124 OPX048` to refresh individual grains.

Most published PLAG STLs (about 199 of 246) are open, non-manifold and not
consistently orientable. Their enclosed volume is therefore undefined, and the
index records it as `null`. For those grains the processed mesh is checked
against the published volume only.

The large surface-fill mesh sets are gitignored. The committed `data/`
bundle is therefore the only copy of the processed surfaces in the repository.

## Local preview

```bash
python -m http.server --directory docs 8765
```

Then open <http://localhost:8765/grain-viewer/>.

URLs take the form `#PLAG124` or `#OPX048/split` (views: `raw`, `proc`,
`overlay`, `split`). Use these to link reviewers to a specific grain.

three.js r170 is vendored under `vendor/` (MIT, see `vendor/three.LICENSE.txt`).
