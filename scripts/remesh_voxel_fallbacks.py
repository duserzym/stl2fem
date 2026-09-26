"""Re-mesh every particle that fell back to the voxel staircase with surface_fill.

Only voxel-fallback particles are touched; particles that the gmsh or MeshFix+gmsh strategies meshed keep their
existing files byte for byte, because those meshes are hash-pinned by downstream campaigns and re-running gmsh on a
different machine does not reproduce them exactly. New meshes go to a separate output root, with a report that
records, per particle, the old and new strategy, node and tetrahedron counts, the mesh volume against the STL's
enclosed volume, whether MeshFix was needed, and the SHA-256 of the metre-scaled mesh.

Particles are processed smallest first, and a particle that already has an output is skipped, so the run can be
interrupted and resumed.

Usage: python scripts/remesh_voxel_fallbacks.py <phase: PLAG|OPX> <existing quality report CSV> <output root>
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import time
import traceback
from pathlib import Path

import meshio
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from stl2fem.conversion import tetrahedralize_bruteforce_stl_for_merrill  # noqa: E402


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    phase, report_path, out_root = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    out_root.mkdir(parents=True, exist_ok=True)
    rows = [r for r in csv.DictReader(open(report_path, newline="")) if r["mesh_strategy"] == "voxel"]
    rows.sort(key=lambda r: int(float(r["n_nodes"] or 0)))
    out_csv = out_root / "_reports" / "surface_fill_remesh.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if out_csv.exists():
        done = {r["particle_id"] for r in csv.DictReader(open(out_csv, newline=""))}
    fields = ["particle_id", "phase", "source_path", "old_strategy", "old_n_nodes", "new_strategy", "n_nodes",
              "n_tets", "mesh_volume_m3", "stl_volume_m3", "volume_relative_error", "repaired", "first_error",
              "edge_median_nm", "edge_p90_nm", "surface_max_edge_nm", "interior_pinned",
              "merrill_msh_path", "sha256", "seconds", "status", "error"]
    new_file = not out_csv.exists()
    with open(out_csv, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        if new_file:
            w.writeheader()
        for r in rows:
            pid = r["particle_id"]
            if pid in done:
                continue
            stl = ROOT / r["source_path"]
            native = out_root / "_native" / f"{pid}.msh"
            meter = out_root / f"{pid}.msh"
            t0 = time.time()
            rec = {"particle_id": pid, "phase": phase, "source_path": r["source_path"], "old_strategy": "voxel",
                   "old_n_nodes": r["n_nodes"], "new_strategy": "", "n_nodes": 0, "n_tets": 0,
                   "mesh_volume_m3": "", "stl_volume_m3": "", "volume_relative_error": "", "repaired": "",
                   "first_error": "", "edge_median_nm": "", "edge_p90_nm": "", "surface_max_edge_nm": "",
                   "interior_pinned": "", "merrill_msh_path": "", "sha256": "", "seconds": 0, "status": "", "error": ""}
            # Each particle runs in its own interpreter: MeshFix and VTK can crash natively on degenerate surfaces
            # (a 29-point open surface took the whole run down without a Python traceback), and a crash must cost
            # only that particle.
            try:
                done_proc = subprocess.run([sys.executable, __file__, "--one", str(stl), str(native), str(meter)],
                                           capture_output=True, text=True, timeout=3600)
                result = json.loads(done_proc.stdout.strip().splitlines()[-1]) if done_proc.returncode == 0 else None
                if result is None:
                    tail = (done_proc.stderr or "").strip().splitlines()[-1:] or [f"exit code {done_proc.returncode}"]
                    rec.update(status="failed", error=f"exit {done_proc.returncode}: {tail[0][:300]}")
                else:
                    rec.update(result, merrill_msh_path=str(meter), sha256=sha256(meter), status="ok")
            except subprocess.TimeoutExpired:
                rec.update(status="failed", error="timeout after 3600 s")
            rec["seconds"] = round(time.time() - t0, 1)
            w.writerow(rec)
            fh.flush()
            print(f"{pid} {rec['status']} nodes {rec['n_nodes']} ({rec['seconds']} s)", flush=True)


def one(stl: Path, native: Path, meter: Path) -> None:
    _, _, units, info = tetrahedralize_bruteforce_stl_for_merrill(
        stl, native, meter, input_unit="um", target_edge_length_m=9e-9, strategy="surface_fill", overwrite=True)
    m = meshio.read(meter)
    tets = m.cells_dict["tetra"]
    p = m.points
    a, b, c, d = (p[tets[:, i]] for i in range(4))
    vol = float(np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))).sum() / 6.0)
    stl_vol = float(info["surface_fill_surface_volume_native"]) * units.input_scale_to_meters ** 3
    if not (np.isfinite(stl_vol) and stl_vol > 0):
        raise RuntimeError("surface encloses no volume")
    e = np.vstack([tets[:, [i, j]] for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))])
    e = np.unique(np.sort(e, axis=1), axis=0)
    lengths = np.linalg.norm(p[e[:, 0]] - p[e[:, 1]], axis=1) * 1e9
    print(json.dumps({"new_strategy": "surface_fill", "n_nodes": int(len(np.unique(tets))), "n_tets": int(len(tets)),
                      "mesh_volume_m3": vol, "stl_volume_m3": stl_vol, "volume_relative_error": vol / stl_vol - 1.0,
                      "repaired": bool(info["surface_fill_repaired"]),
                      "first_error": str(info["surface_fill_first_error"])[:200],
                      "edge_median_nm": float(np.median(lengths)), "edge_p90_nm": float(np.percentile(lengths, 90)),
                      "surface_max_edge_nm": float(info["surface_fill_surface_max_edge_native"]) * 1000.0,
                      "interior_pinned": bool(info["surface_fill_interior_pinned"])}))


if __name__ == "__main__":
    if len(sys.argv) == 5 and sys.argv[1] == "--one":
        one(Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]))
    else:
        main()
