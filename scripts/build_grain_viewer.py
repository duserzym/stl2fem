"""Build the data bundle for the static magnetite-grain review viewer.

For every PLAG- and OPX-hosted magnetite grain in the Nikolaisen et al. (2022)
FIB-SEM dataset this writes two compact surface meshes into
``docs/grain-viewer/data``:

* ``raw/<ID>.bin.gz``  - the published per-particle STL (vertices welded);
* ``proc/<ID>.bin.gz`` - the boundary surface of the tetrahedral Gmsh 2.2 mesh
  currently used for micromagnetic FEM modeling, as selected by the mesh
  inventories (PLAG: 245 surface-fill inventory, OPX: 82 9 nm v2 inventory).

Both meshes are expressed in nanometres relative to the centre of the raw STL
bounding box, so the two geometries stay registered in the sample frame.
``index.json`` carries the per-grain bookkeeping metadata and provenance.

Binary layout (little endian, then gzip):
    0  char[4]   magic "GVM1"
    4  uint32    n_vertices
    8  uint32    n_triangles
    12 uint32    flags (bit 0: uint16 indices, else uint32)
    16 float32x3 quantization origin (nm)
    28 float32x3 quantization step (nm)
    40 uint16    positions [n_vertices * 3], padded to 4 bytes
       uint16/32 triangle indices [n_triangles * 3]

The mesh inventories and the 100-grain cohort list belong to the downstream
PINT_with_reversal campaigns; by default they are read from sibling checkouts.
Only numpy is required.  Example:
    python scripts/build_grain_viewer.py
    python scripts/build_grain_viewer.py --only PLAG124 OPX048
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import struct
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BATCHES = REPO_ROOT.parent / "PINT_campaign_opx" / "code" / "reversal_paleointensity_probability" / "production_batches"
DEFAULT_COHORT = REPO_ROOT.parent / "PINT_with_reversal" / "cohort" / "index.csv"
DEFAULT_OUT = REPO_ROOT / "docs" / "grain-viewer" / "data"

HOSTS = {
    "PLAG": {
        "label": "Plagioclase",
        "stl_dir": "Plag Binary meshes",
        "stl_name": "{id}-binary.stl",
        "stl_meta": "Plag_stl_B16.csv",
        "hyst_meta": "Plag_B16.csv",
        "inventory": "nikolaisen2022_mesh_inventory_245_surface_fill.csv",
    },
    "OPX": {
        "label": "Orthopyroxene",
        "stl_dir": "OPX Binary meshes",
        "stl_name": "{id}-binary.stl",
        "stl_meta": "OPX_stl_B16.csv",
        "hyst_meta": "OPX_B16.csv",
        "inventory": "nikolaisen2022_opx_mesh_inventory_82_9nm_v2.csv",
    },
}


# ----------------------------------------------------------------------------
# Readers
# ----------------------------------------------------------------------------

def read_stl(path: Path) -> np.ndarray:
    """Return an (n, 3, 3) float64 triangle soup; binary or ASCII STL."""
    data = path.read_bytes()
    if len(data) >= 84:
        n = struct.unpack_from("<I", data, 80)[0]
        if len(data) == 84 + 50 * n:
            rec = np.dtype([("normal", "<f4", 3), ("v", "<f4", (3, 3)), ("attr", "<u2")])
            return np.frombuffer(data, dtype=rec, count=n, offset=84)["v"].astype(np.float64)
    # ASCII (some files in the "Binary" folders are ASCII).
    coords = [
        [float(t) for t in line.split()[1:4]]
        for line in data.decode("ascii", errors="replace").splitlines()
        if line.strip().startswith("vertex")
    ]
    return np.asarray(coords, dtype=np.float64).reshape(-1, 3, 3)


def weld(soup: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    flat = soup.reshape(-1, 3)
    verts, inverse = np.unique(flat, axis=0, return_inverse=True)
    faces = inverse.reshape(-1, 3)
    keep = (faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2]) & (faces[:, 0] != faces[:, 2])
    return verts, faces[keep]


def read_gmsh22(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised reader for binary (or ASCII) Gmsh 2.2 meshes -> nodes, tets."""
    data = path.read_bytes()

    def next_line(pos: int) -> int:
        return data.index(b"\n", pos) + 1

    fmt_pos = next_line(data.index(b"$MeshFormat"))
    version, file_type, _ = data[fmt_pos:next_line(fmt_pos)].split()
    if not version.startswith(b"2."):
        raise ValueError(f"{path}: unsupported Gmsh version {version!r}")
    binary = file_type == b"1"

    pos = next_line(data.index(b"$Nodes"))
    end = next_line(pos)
    n_nodes = int(data[pos:end])
    pos = end
    if binary:
        rec = np.dtype([("tag", "<i4"), ("xyz", "<f8", 3)])
        block = np.frombuffer(data, dtype=rec, count=n_nodes, offset=pos)
        tags, xyz = block["tag"].astype(np.int64), block["xyz"].astype(np.float64)
    else:
        rows = np.loadtxt(data[pos:data.index(b"$EndNodes")].decode().splitlines(), ndmin=2)
        tags, xyz = rows[:, 0].astype(np.int64), rows[:, 1:4]
    lookup = np.full(tags.max() + 1, -1, dtype=np.int64)
    lookup[tags] = np.arange(n_nodes)

    pos = next_line(data.index(b"$Elements"))
    end = next_line(pos)
    n_elem = int(data[pos:end])
    pos = end
    nodes_per = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5, 8: 3, 9: 6, 11: 10, 15: 1}
    tets = []
    if binary:
        seen = 0
        while seen < n_elem:
            etype, count, ntags = struct.unpack_from("<iii", data, pos)
            pos += 12
            width = 1 + ntags + nodes_per[etype]
            block = np.frombuffer(data, dtype="<i4", count=count * width, offset=pos).reshape(count, width)
            pos += 4 * count * width
            if etype == 4:
                tets.append(block[:, -4:])
            seen += count
    else:
        for line in data[pos:data.index(b"$EndElements")].decode().splitlines():
            parts = [int(p) for p in line.split()]
            if len(parts) > 2 and parts[1] == 4:
                tets.append(np.asarray([parts[-4:]]))
    if not tets:
        raise ValueError(f"{path}: no tetrahedra")
    return xyz, lookup[np.concatenate(tets).astype(np.int64)]


# ----------------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------------

def boundary_surface(xyz: np.ndarray, tets: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Outward-oriented boundary triangles of a tetrahedral mesh."""
    local = np.array([[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 3, 1], [1, 2, 3, 0]])
    cand = tets[:, local].reshape(-1, 4)  # three face nodes + opposite node
    key = np.sort(cand[:, :3], axis=1)
    _, inverse, counts = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    faces = cand[counts[inverse.ravel()] == 1]
    a, b, c, opp = (xyz[faces[:, i]] for i in range(4))
    outward = np.einsum("ij,ij->i", np.cross(b - a, c - a), (a + b + c) / 3.0 - opp) >= 0
    tri = faces[:, :3].copy()
    tri[~outward] = tri[~outward][:, [0, 2, 1]]
    used, remap = np.unique(tri, return_inverse=True)
    return xyz[used], remap.reshape(-1, 3)


def signed_volume(verts: np.ndarray, faces: np.ndarray) -> float:
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def surface_area(verts: np.ndarray, faces: np.ndarray) -> float:
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum())


def edge_lengths(verts: np.ndarray, cells: np.ndarray) -> np.ndarray:
    k = cells.shape[1]
    pairs = np.array([(i, j) for i in range(k) for j in range(i + 1, k)])
    edges = np.sort(cells[:, pairs].reshape(-1, 2), axis=1)
    edges = np.unique(edges, axis=0)
    return np.linalg.norm(verts[edges[:, 0]] - verts[edges[:, 1]], axis=1)


def topology(n_vertices: int, faces: np.ndarray) -> dict:
    """Connected components and edge-manifold/orientation counts of a triangle surface."""
    directed = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    _, ucount = np.unique(np.sort(directed, axis=1), axis=0, return_counts=True)
    _, dcount = np.unique(directed, axis=0, return_counts=True)
    label = np.arange(n_vertices)
    while True:  # label propagation to the minimum vertex index per component
        low = np.minimum(label[directed[:, 0]], label[directed[:, 1]])
        new = label.copy()
        np.minimum.at(new, directed[:, 0], low)
        np.minimum.at(new, directed[:, 1], low)
        new = new[new]
        if np.array_equal(new, label):
            break
        label = new
    return {
        "components": int(len(np.unique(label[faces[:, 0]]))),
        "boundary_edges": int((ucount == 1).sum()),
        "non_manifold_edges": int((ucount > 2).sum()),
        "misoriented_edges": int((dcount > 1).sum()),
    }


def tet_volumes(xyz: np.ndarray, tets: np.ndarray) -> np.ndarray:
    a, b, c, d = (xyz[tets[:, i]] for i in range(4))
    return np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))) / 6.0


def encode(verts_nm: np.ndarray, faces: np.ndarray) -> bytes:
    lo = verts_nm.min(axis=0)
    span = np.maximum(verts_nm.max(axis=0) - lo, 1e-9)
    step = span / 65535.0
    q = np.rint((verts_nm - lo) / step).astype("<u2")
    small = len(verts_nm) < 65536
    head = b"GVM1" + struct.pack("<III", len(verts_nm), len(faces), 1 if small else 0)
    head += np.asarray(lo, "<f4").tobytes() + np.asarray(step, "<f4").tobytes()
    pos = q.tobytes()
    pos += b"\0" * (-len(pos) % 4)
    idx = faces.astype("<u2" if small else "<u4").tobytes()
    return gzip.compress(head + pos + idx, compresslevel=9, mtime=0)


def stats(values: np.ndarray) -> dict:
    return {"min": float(values.min()), "median": float(np.median(values)), "max": float(values.max())}


# ----------------------------------------------------------------------------
# Metadata
# ----------------------------------------------------------------------------

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_rows(path: Path, key: str) -> dict[str, dict]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        return {row[key].strip(): row for row in csv.DictReader(fh)}


def to_float(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def relocate(mesh_path: str, stl2fem: Path) -> Path:
    """Inventory paths are absolute; re-anchor them on the stl2fem data root."""
    p = Path(mesh_path)
    return stl2fem / "data" / p.parent.name / p.name


def git_head(path: Path) -> str | None:
    """Resolve HEAD without needing git on PATH."""
    git = path / ".git"
    try:
        if git.is_file():  # worktree: "gitdir: <path>"
            git = Path(git.read_text().split(":", 1)[1].strip())
        head = (git / "HEAD").read_text().strip()
        if not head.startswith("ref:"):
            return head
        ref = head.split(" ", 1)[1]
        common = Path((git / "commondir").read_text().strip()) if (git / "commondir").exists() else git
        common = common if common.is_absolute() else (git / common).resolve()
        for base in (git, common):
            if (base / ref).exists():
                return (base / ref).read_text().strip()
        for line in (common / "packed-refs").read_text().splitlines():
            if line.endswith(" " + ref):
                return line.split()[0]
    except OSError:
        pass
    return None


# ----------------------------------------------------------------------------

def build(args: argparse.Namespace) -> None:
    out = args.out
    (out / "raw").mkdir(parents=True, exist_ok=True)
    (out / "proc").mkdir(parents=True, exist_ok=True)
    data_root = REPO_ROOT / "data" / "Nikolaisen2022"
    cohort = read_rows(args.cohort, "grain_id") if args.cohort.exists() else {}

    grains, inventories = [], {}
    for host, cfg in HOSTS.items():
        geo = read_rows(data_root / cfg["stl_meta"], "Filename")
        hyst = read_rows(data_root / cfg["hyst_meta"], "Filename")
        inv_path = args.batches / cfg["inventory"]
        inventory = read_rows(inv_path, "grain_id")
        inventories[host] = {"file": cfg["inventory"], "sha256": sha256(inv_path), "rows": len(inventory)}

        ids = sorted(geo)
        if args.only:
            ids = [g for g in ids if g in args.only]
        for gid in ids:
            print(f"{gid} ...", end=" ", flush=True)
            meta = geo[gid]
            entry = {
                "id": gid,
                "host": host,
                "published": {
                    "volume_um3": to_float(meta.get("Volume") or meta.get("Volume (µm³)")),
                    "evsd_um": to_float(meta.get("EVSD (mu)")),
                    "oblate_sphere_prolate": to_float(meta.get("Ob_Sph_Pro")),
                },
                "cohort100": gid in cohort,
            }
            if gid in cohort:
                entry["cohort_status"] = cohort[gid].get("status")
            h = hyst.get(gid)
            if h and (h.get("Domain_Structure") or "").strip():
                entry["published"]["domain_state"] = h["Domain_Structure"].strip().replace("_", " ")

            # --- raw STL ------------------------------------------------------
            stl_path = data_root / cfg["stl_dir"] / cfg["stl_name"].format(id=gid)
            verts_um, faces = weld(read_stl(stl_path))
            verts_nm = verts_um * 1e3
            centre = 0.5 * (verts_nm.min(axis=0) + verts_nm.max(axis=0))
            topo = topology(len(verts_nm), faces)
            # Most published PLAG surfaces are open, non-manifold and not consistently orientable;
            # their enclosed volume is undefined, so it is only reported for closed manifolds.
            closed = not (topo["boundary_edges"] or topo["non_manifold_edges"] or topo["misoriented_edges"])
            vol = signed_volume(verts_nm, faces)
            if closed and vol < 0:
                faces = faces[:, [0, 2, 1]]
            vol = abs(vol) if closed else None
            blob = encode(verts_nm - centre, faces)
            (out / "raw" / f"{gid}.bin.gz").write_bytes(blob)
            entry["centre_nm"] = [round(float(c), 3) for c in centre]
            entry["raw"] = {
                "file": f"raw/{gid}.bin.gz",
                "source": f"Nikolaisen2022/{cfg['stl_dir']}/{stl_path.name}",
                "source_sha256": sha256(stl_path),
                "n_vertices": int(len(verts_nm)),
                "n_triangles": int(len(faces)),
                "volume_nm3": vol,
                "area_nm2": surface_area(verts_nm, faces),
                "extent_nm": [float(e) for e in np.ptp(verts_nm, axis=0)],
                "edge_nm": stats(edge_lengths(verts_nm, faces)),
                "closed_manifold": closed,
                **topo,
                "bytes": len(blob),
            }
            pub_vol = entry["published"]["volume_um3"]
            if pub_vol and closed:
                entry["raw"]["volume_rel_diff_vs_published"] = (vol - pub_vol * 1e9) / (pub_vol * 1e9)

            # --- processed FEM mesh --------------------------------------------
            row = inventory.get(gid)
            strategy = (row or {}).get("mesh_strategy") or "unknown"
            if row is None:
                entry["proc"] = None
                entry["proc_note"] = "No processed mesh in the current inventory."
            else:
                msh = relocate(row["mesh_path"], REPO_ROOT)
                digest = sha256(msh)
                xyz_m, tets = read_gmsh22(msh)
                xyz_nm = xyz_m * 1e9
                sverts, sfaces = boundary_surface(xyz_nm, tets)
                blob = encode(sverts - centre, sfaces)
                (out / "proc" / f"{gid}.bin.gz").write_bytes(blob)
                pvol = float(tet_volumes(xyz_nm, tets).sum())
                entry["proc"] = {
                    "file": f"proc/{gid}.bin.gz",
                    "source": f"{msh.parent.name}/{msh.name}",
                    "strategy": strategy,
                    "sha256": digest,
                    "sha256_matches_inventory": digest == row.get("sha256"),
                    "n_nodes": int(len(xyz_nm)),
                    "n_tets": int(len(tets)),
                    "n_surface_vertices": int(len(sverts)),
                    "n_surface_triangles": int(len(sfaces)),
                    "volume_nm3": pvol,
                    "area_nm2": surface_area(sverts, sfaces),
                    "volume_rel_diff_vs_raw": (pvol - vol) / vol if closed else None,
                    "extent_nm": [float(e) for e in np.ptp(xyz_nm, axis=0)],
                    "tet_edge_nm": stats(edge_lengths(xyz_nm, tets)),
                    "components": topology(len(sverts), sfaces)["components"],
                    "bytes": len(blob),
                }
                if pub_vol:
                    entry["proc"]["volume_rel_diff_vs_published"] = (pvol - pub_vol * 1e9) / (pub_vol * 1e9)
                note = {
                    "voxel": "Voxel-fallback mesh (stair-stepped); surface-fill remesh failed.",
                    "voxel_not_surface_filled": "Voxel-fallback mesh (stair-stepped); surface-fill remesh not applied.",
                    "surface_fill_9nm_sliver_refused": "9 nm surface-fill remesh refused by the sliver gate.",
                }.get(strategy)
                if note:
                    entry["proc_note"] = note
            grains.append(entry)
            print("ok")

    index = {
        "schema": 1,
        "generated": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
        "units": "nm; coordinates relative to the raw-STL bounding-box centre in the FIB sample frame",
        "provenance": {
            "dataset": "Nikolaisen et al. (2022), silicate-hosted magnetite FIB-SEM STL meshes",
            "stl2fem_commit": git_head(REPO_ROOT),
            "inventory_repo_commit": git_head(args.batches.parents[2]),
            "inventories": inventories,
            "cohort100_source": "PINT_with_reversal/cohort/index.csv" if cohort else None,
            "builder": "scripts/build_grain_viewer.py",
        },
        "grains": grains,
    }
    if args.only and (out / "index.json").exists():  # partial rebuild: merge
        old = json.loads((out / "index.json").read_text())
        keep = {g["id"]: g for g in old["grains"]}
        keep.update({g["id"]: g for g in grains})
        index["grains"] = sorted(keep.values(), key=lambda g: (g["host"] != "PLAG", g["id"]))
    (out / "index.json").write_text(json.dumps(index, separators=(",", ":")), encoding="utf-8")
    total = sum(f.stat().st_size for f in out.rglob("*") if f.is_file())
    print(f"wrote {len(index['grains'])} grains, {total / 2**20:.1f} MiB -> {out}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--batches", type=Path, default=DEFAULT_BATCHES, help="directory holding the mesh inventory CSVs")
    p.add_argument("--cohort", type=Path, default=DEFAULT_COHORT, help="PINT_with_reversal cohort/index.csv (optional)")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--only", nargs="*", help="rebuild only these grain IDs (merged into an existing index.json)")
    build(p.parse_args())


if __name__ == "__main__":
    sys.exit(main())
