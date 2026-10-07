"""Convert the Gergov et al. (2025) Hekla/Vesuvius magnetite meshes to merrill.jl-ready Gmsh 2.2 files.

Source: Gergov, Muxworthy, Williams & Cowan (2025), "Magnetic recording fidelity of basalts through 3D
nanotomography", G-cubed, doi:10.1029/2024GC011776; data doi:10.5281/zenodo.11369780 (CC-BY-4.0).
The archives ``Individual Meshes Hekla.7z`` / ``Individual Meshes Vesuvius.7z`` hold one smoothed CUBIT
tetrahedral mesh per grain as a Patran neutral file in micrometres.

The tetrahedral meshes are kept exactly (same nodes, same elements, same order); only the coordinates are
scaled from micrometres to metres, the stl2fem convention, and written as binary Gmsh 2.2 with one physical
volume named ``particle``. Quality metrics are computed from the volume mesh and joined with the published
per-grain metrics.

Usage (from the repository root, stl2fem venv):
    python scripts/convert_gergov2025.py [--data data/Gergov2025] [--workers 2] [--only H009 V002]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from stl2fem.patran import read_patran_neutral  # noqa: E402

UM_TO_M = 1e-6
METRIC_FILES = {"Hekla": "HeklaGrainMetrics.csv", "Vesuvius": "VesuviusGrainMetrics.csv"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_msh22_binary(path: Path, points_m: np.ndarray, tets: np.ndarray) -> None:
    """Binary Gmsh 2.2, node tags 1..n in input order, one tetra block, physical/geometrical tag 1."""
    n, m = len(points_m), len(tets)
    nodes = np.empty(n, dtype=[("tag", "<i4"), ("xyz", "<f8", 3)])
    nodes["tag"] = np.arange(1, n + 1)
    nodes["xyz"] = points_m
    elems = np.empty((m, 7), dtype="<i4")
    elems[:, 0] = np.arange(1, m + 1)
    elems[:, 1:3] = 1  # physical, geometrical
    elems[:, 3:] = tets + 1
    tmp = path.with_suffix(".msh.part")
    with tmp.open("wb") as fh:
        fh.write(b"$MeshFormat\n2.2 1 8\n")
        fh.write(np.array([1], dtype="<i4").tobytes())
        fh.write(b"\n$EndMeshFormat\n$PhysicalNames\n1\n3 1 \"particle\"\n$EndPhysicalNames\n")
        fh.write(f"$Nodes\n{n}\n".encode())
        fh.write(nodes.tobytes())
        fh.write(f"\n$EndNodes\n$Elements\n{m}\n".encode())
        fh.write(np.array([4, m, 2], dtype="<i4").tobytes())  # type tetra, count, 2 tags
        fh.write(elems.tobytes())
        fh.write(b"\n$EndElements\n")
    tmp.replace(path)


def quality(coords_nm: np.ndarray, tets: np.ndarray) -> dict:
    a, b, c, d = (coords_nm[tets[:, k]] for k in range(4))
    signed = np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a)) / 6.0
    vol = np.abs(signed)
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    sq = np.stack([((coords_nm[tets[:, i]] - coords_nm[tets[:, j]]) ** 2).sum(1) for i, j in pairs], axis=1)
    # mean-ratio style shape measure: 1 for a regular tetrahedron, -> 0 for slivers
    shape = 6.0 * np.sqrt(2.0) * vol / np.maximum(sq.mean(1), 1e-300) ** 1.5
    edges = np.unique(np.sort(tets[:, pairs].reshape(-1, 2), axis=1), axis=0)
    elen = np.linalg.norm(coords_nm[edges[:, 0]] - coords_nm[edges[:, 1]], axis=1)
    faces = np.sort(tets[:, [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]].reshape(-1, 3), axis=1)
    ufaces, fcount = np.unique(faces, axis=0, return_counts=True)
    bfaces = ufaces[fcount == 1]
    bedges = np.sort(bfaces[:, [[0, 1], [1, 2], [0, 2]]].reshape(-1, 2), axis=1)
    _, ecount = np.unique(bedges, axis=0, return_counts=True)
    # connected components of the tetra graph (shared nodes)
    label = np.arange(len(coords_nm))
    links = np.concatenate([tets[:, [0, k]] for k in (1, 2, 3)])
    while True:
        low = np.minimum(label[links[:, 0]], label[links[:, 1]])
        new = label.copy()
        np.minimum.at(new, links[:, 0], low)
        np.minimum.at(new, links[:, 1], low)
        new = new[new]
        if np.array_equal(new, label):
            break
        label = new
    used = np.unique(tets)
    volume = float(vol.sum())
    return {
        "n_nodes_used": int(len(used)),
        "n_tets": int(len(tets)),
        "volume_um3": volume * 1e-9,
        "esd_nm": float((6.0 * volume / np.pi) ** (1.0 / 3.0)),
        "n_components": int(len(np.unique(label[used]))),
        "n_negative_tets": int((signed < 0).sum()),
        "n_degenerate_tets": int((vol <= 1e-12 * np.median(vol)).sum()),
        "shape_min": float(shape.min()),
        "shape_p01": float(np.percentile(shape, 1)),
        "shape_median": float(np.median(shape)),
        "edge_min_nm": float(elen.min()),
        "edge_median_nm": float(np.median(elen)),
        "edge_p95_nm": float(np.percentile(elen, 95)),
        "edge_max_nm": float(elen.max()),
        "n_boundary_triangles": int(len(bfaces)),
        "boundary_closed_manifold": bool((ecount == 2).all()),
        **{f"extent_{ax}_nm": float(e) for ax, e in zip("xyz", np.ptp(coords_nm[used], axis=0))},
    }


def convert_one(job: tuple[str, str, str, str]) -> dict:
    grain, locality, src, dst = job
    src, dst = Path(src), Path(dst)
    mesh = read_patran_neutral(src)
    pids = sorted(set(mesh.property_ids.tolist()))
    orphan = len(mesh.coords) - len(np.unique(mesh.tets))
    write_msh22_binary(dst, mesh.coords * UM_TO_M, mesh.tets)
    row = {
        "grain_id": grain,
        "locality": locality,
        "source_pat": src.name,
        "source_sha256": sha256(src),
        "msh_file": dst.name,
        "msh_sha256": sha256(dst),
        "msh_bytes": dst.stat().st_size,
        "n_nodes": len(mesh.coords),
        "n_orphan_nodes": int(orphan),
        "property_ids": ";".join(map(str, pids)),
        "skipped_non_tet_elements": mesh.skipped_elements,
    }
    row.update(quality(mesh.coords * 1e3, mesh.tets))
    return row


def published_metrics(zenodo: Path) -> dict[str, dict]:
    """Published per-grain metrics keyed by mesh ID.

    A block of Hekla rows is labelled by STL file name ("828.stl") instead of mesh ID ("H828"); those are mapped
    to the Hekla ID with the same number.
    """
    out = {}
    for locality, name in METRIC_FILES.items():
        with (zenodo / name).open(newline="", encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                r = {k.strip(): (v.strip() if isinstance(v, str) else v) for k, v in r.items()}
                key = r["Mesh ID"]
                if key.endswith(".stl") and key[:-4].isdigit():
                    key = f"{locality[0]}{int(key[:-4]):03d}"
                out[key] = r
    return out


PUBLISHED_FIELDS = {
    "published_esd_nm": "Equivalent sphere diameter (nm)",
    "published_volume_um3": "Volume (um3)",
    "published_flinn_ratio": "Flinn Ratio",
    "published_ground_state": "Ground State",
}


def join_published(row: dict, meta: dict[str, dict]) -> dict:
    pub = meta.get(row["grain_id"], {})
    for col, src in PUBLISHED_FIELDS.items():
        row[col] = pub.get(src, "")
    row["published_lem_states"] = ";".join(
        s for s in (pub.get("LEM State 1"), pub.get("LEM State 2"), pub.get("LEM State 3")) if s)
    try:
        row["esd_rel_diff_vs_published"] = float(row["esd_nm"]) / float(row["published_esd_nm"]) - 1.0
    except (TypeError, ValueError, ZeroDivisionError):
        row["esd_rel_diff_vs_published"] = ""
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    root = Path(__file__).resolve().parents[1]
    ap.add_argument("--data", type=Path, default=root / "data" / "Gergov2025")
    ap.add_argument("--reports", type=Path, default=root / "reports" / "Gergov2025")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--metadata-only", action="store_true",
                    help="re-join the published metrics into an existing inventory without reconverting")
    args = ap.parse_args()

    if args.metadata_only:
        inv = args.reports / "inventory.csv"
        with inv.open(newline="") as fh:
            rows = list(csv.DictReader(fh))
        meta = published_metrics(args.data / "_zenodo")
        rows = [join_published(r, meta) for r in rows]
        with inv.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        missing = [r["grain_id"] for r in rows if not r["published_esd_nm"]]
        print(f"re-joined published metrics for {len(rows)} meshes; {len(missing)} without a published row: {missing}")
        return

    src_dir = args.data / "_extracted"
    out_dir = args.data / "merrill_msh"
    out_dir.mkdir(parents=True, exist_ok=True)
    args.reports.mkdir(parents=True, exist_ok=True)
    jobs = []
    for pat in sorted(src_dir.glob("*_mesh.pat")):
        grain = pat.name.split("_")[0]
        if args.only and grain not in args.only:
            continue
        locality = {"H": "Hekla", "V": "Vesuvius"}[grain[0]]
        jobs.append((grain, locality, str(pat), str(out_dir / f"{grain}.msh")))
    # largest first so the long files do not straggle at the end
    jobs.sort(key=lambda j: -Path(j[2]).stat().st_size)
    meta = published_metrics(args.data / "_zenodo")

    rows, failures = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for job, fut in [(j, pool.submit(convert_one, j)) for j in jobs]:
            try:
                row = fut.result()
            except Exception as e:  # keep going; report at the end
                failures.append((job[0], repr(e)))
                print(f"{job[0]}: FAILED {e!r}", flush=True)
                continue
            rows.append(join_published(row, meta))
            print(f"{row['grain_id']}: {row['n_nodes']} nodes, {row['n_tets']} tets, "
                  f"edge median {row['edge_median_nm']:.1f} nm", flush=True)

    rows.sort(key=lambda r: r["grain_id"])
    if args.only and (args.reports / "inventory.csv").exists():
        with (args.reports / "inventory.csv").open(newline="") as fh:
            old = {r["grain_id"]: r for r in csv.DictReader(fh)}
        old.update({r["grain_id"]: r for r in rows})
        rows = [old[k] for k in sorted(old)]
    with (args.reports / "inventory.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    if failures:
        with (args.reports / "failures.csv").open("w", newline="") as fh:
            csv.writer(fh).writerows([("grain_id", "error"), *failures])
    print(f"converted {len(rows)} meshes, {len(failures)} failures")


if __name__ == "__main__":
    main()
