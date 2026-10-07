# Magnetite grain viewer

A static three.js page for reviewing FIB-SEM magnetite grains against the tetrahedral meshes that stl2fem makes
for micromagnetic modeling. It is published at <https://duserzym.github.io/stl2fem/grain-viewer/>.

| Dataset | Grains | What is shown |
|---|---|---|
| Nikolaisen et al. (2022), plagioclase- and orthopyroxene-hosted | PLAG001–246, OPX001–082 | **Original**: the published per-particle STL.<br>**Processed**: the boundary of the tetrahedral mesh that the PINT campaigns use. |
| Gergov et al. (2025), Hekla 1991 and Vesuvius 1944 basalts | 805 H###, 171 V### | The published smoothed CUBIT mesh, converted exactly on branch `Gergov2025_stl2msh`. These grains were published as meshes, so there is no separate original surface. |

**View modes:** original, processed, overlay and side by side.

**Scale bar:** the camera is orthographic, so the scale bar is exact anywhere in the view.

**Review aids:** the page also shows automatic geometry and topology checks, and provenance hashes.

**Review log:** each grain has a review log. Notes are kept in the reviewer's browser and are exported or
imported as CSV.

## Rebuilding the data

```bash
python scripts/build_grain_viewer.py                      # both datasets
python scripts/build_grain_viewer.py --datasets gergov    # one dataset, merged into the existing index
python scripts/build_grain_viewer.py --only PLAG124 H009  # individual grains
```

The script needs only numpy.

**Nikolaisen inputs:**

- the raw STLs, from `data/Nikolaisen2022/`;
- the meshes, from the `data/Nikolaisen2022*_merrill_msh*` directories;
- the per-grain mesh selection, from the PINT_with_reversal inventories:
  - `--batches`, default `../PINT_campaign_opx/.../production_batches`;
  - files `nikolaisen2022_mesh_inventory_245_surface_fill.csv` and `nikolaisen2022_opx_mesh_inventory_82_9nm_v2.csv`.

**Gergov inputs:** a checkout of branch `Gergov2025_stl2msh` with converted meshes.

- Pass it with `--gergov-root`; the default is `../stl2fem_gergov2025`.
- The script reads its `reports/Gergov2025/inventory.csv` and `data/Gergov2025/merrill_msh/`.

**Outputs:**

- Every mesh is verified against the SHA-256 in its inventory.
- The script writes `data/index.json` plus one quantized, gzipped surface per grain and variant.
- The data bundle is about 92 MB. It is gitignored on the source branches.

**Volume checks:** most published PLAG STLs (about 199 of 246) are open, non-manifold and not consistently
orientable.

- Their enclosed volume is therefore undefined, and the index records it as `null`.
- For those grains, and for all Gergov grains, the mesh is checked against the published volume only.

## Publishing

```powershell
./scripts/publish_site.ps1
```

This replaces the `gh-pages` branch with a single commit holding `docs/` and the current data. GitHub Pages serves
that branch, so rebuilds never accumulate in the repository history.

## Local preview

```bash
python -m http.server --directory docs 8765
```

Then open <http://localhost:8765/grain-viewer/>.

URLs take the form `#PLAG124`, `#OPX048/split` or `#H009`. The views are `raw`, `proc`, `overlay` and `split`.

three.js r170 is vendored under `vendor/` (MIT, see `vendor/three.LICENSE.txt`).
