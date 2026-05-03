"""Mesh quality inspection using PyVista/VTK."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np


def _require_pyvista():
    try:
        import pyvista as pv
    except ImportError as exc:
        raise ImportError(
            "PyVista is required for mesh quality inspection. "
            "Install with `pip install -e .[examples]` or `pip install pyvista`."
        ) from exc
    return pv


def _stats(prefix: str, values: Iterable[float]) -> dict[str, float]:
    arr = np.asarray(list(values), dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {
            f"{prefix}_min": np.nan,
            f"{prefix}_median": np.nan,
            f"{prefix}_p95": np.nan,
            f"{prefix}_max": np.nan,
        }
    return {
        f"{prefix}_min": float(np.min(arr)),
        f"{prefix}_median": float(np.median(arr)),
        f"{prefix}_p95": float(np.percentile(arr, 95)),
        f"{prefix}_max": float(np.max(arr)),
    }


def _quality_stats(mesh, measures: tuple[str, ...]) -> dict[str, float]:
    rows: dict[str, float] = {}
    try:
        quality = mesh.cell_quality(list(measures))
    except Exception:
        return rows

    for measure in measures:
        if measure in quality.cell_data:
            rows.update(_stats(measure, quality.cell_data[measure]))
    return rows


def _tetra_edge_lengths(mesh) -> np.ndarray:
    """Return all tetrahedron edge lengths, counting shared edges per cell."""

    if mesh.n_cells == 0:
        return np.asarray([], dtype=float)

    cells = np.asarray(mesh.cells)
    if cells.size % 5 != 0:
        return np.asarray([], dtype=float)

    tets = cells.reshape(-1, 5)
    if not np.all(tets[:, 0] == 4):
        return np.asarray([], dtype=float)

    ids = tets[:, 1:5]
    points = np.asarray(mesh.points)
    edge_pairs = np.asarray(
        [
            [0, 1],
            [0, 2],
            [0, 3],
            [1, 2],
            [1, 3],
            [2, 3],
        ]
    )
    starts = points[ids[:, edge_pairs[:, 0]]]
    ends = points[ids[:, edge_pairs[:, 1]]]
    return np.linalg.norm(ends - starts, axis=2).ravel()


def load_surface(path: str | Path):
    pv = _require_pyvista()
    mesh = pv.read(str(path))
    if not isinstance(mesh, pv.PolyData):
        mesh = mesh.extract_surface()
    return mesh.triangulate().clean()


def inspect_surface_stl(path: str | Path) -> dict[str, object]:
    """Inspect an STL surface for topology and triangle quality."""

    path = Path(path)
    surface = load_surface(path)
    bounds = surface.bounds
    diagonal = np.linalg.norm(
        [
            bounds[1] - bounds[0],
            bounds[3] - bounds[2],
            bounds[5] - bounds[4],
        ]
    )

    boundary_edges = surface.extract_feature_edges(
        boundary_edges=True,
        non_manifold_edges=False,
        feature_edges=False,
        manifold_edges=False,
    )
    non_manifold_edges = surface.extract_feature_edges(
        boundary_edges=False,
        non_manifold_edges=True,
        feature_edges=False,
        manifold_edges=False,
    )

    connected = surface.connectivity()
    region_ids = connected.cell_data.get("RegionId", [])
    n_components = int(len(np.unique(region_ids))) if len(region_ids) else 1

    row: dict[str, object] = {
        "path": str(path),
        "file_size_bytes": path.stat().st_size,
        "n_points_surface": int(surface.n_points),
        "n_triangles_surface": int(surface.n_cells),
        "surface_area": float(surface.area),
        "surface_volume": float(surface.volume),
        "bbox_xmin": float(bounds[0]),
        "bbox_xmax": float(bounds[1]),
        "bbox_ymin": float(bounds[2]),
        "bbox_ymax": float(bounds[3]),
        "bbox_zmin": float(bounds[4]),
        "bbox_zmax": float(bounds[5]),
        "bbox_diagonal": float(diagonal),
        "boundary_edge_count": int(boundary_edges.n_cells),
        "non_manifold_edge_count": int(non_manifold_edges.n_cells),
        "connected_components": n_components,
        "is_watertight": bool(boundary_edges.n_cells == 0),
    }
    row.update(_quality_stats(surface, ("area", "min_angle", "max_angle", "aspect_ratio")))
    return row


def load_volume_mesh(path: str | Path):
    pv = _require_pyvista()
    path = Path(path)
    try:
        mesh = pv.read(str(path))
        if hasattr(mesh, "celltypes"):
            tet_mask = mesh.celltypes == pv.CellType.TETRA
            if tet_mask.any():
                return mesh.extract_cells(tet_mask)
        return mesh
    except Exception:
        try:
            import meshio
        except ImportError as exc:
            raise ImportError(
                "Reading this volume mesh requires meshio. Install with `pip install meshio`."
            ) from exc

        meshio_mesh = meshio.read(path)
        tets = meshio_mesh.cells_dict.get("tetra")
        if tets is None:
            raise ValueError(f"No tetrahedral cells found in {path}")
        return pv.UnstructuredGrid({pv.CellType.TETRA: tets}, meshio_mesh.points)


def inspect_volume_mesh(path: str | Path) -> dict[str, object]:
    """Inspect a tetrahedral volume mesh."""

    path = Path(path)
    mesh = load_volume_mesh(path)
    bounds = mesh.bounds
    row: dict[str, object] = {
        "msh_path": str(path),
        "msh_size_bytes": path.stat().st_size,
        "msh_size_mib": path.stat().st_size / 1024**2,
        "n_nodes": int(mesh.n_points),
        "n_tets": int(mesh.n_cells),
        "bbox_xmin": float(bounds[0]),
        "bbox_xmax": float(bounds[1]),
        "bbox_ymin": float(bounds[2]),
        "bbox_ymax": float(bounds[3]),
        "bbox_zmin": float(bounds[4]),
        "bbox_zmax": float(bounds[5]),
    }
    row.update(
        _quality_stats(
            mesh,
            (
                "volume",
                "scaled_jacobian",
                "radius_ratio",
                "aspect_ratio",
                "min_angle",
                "max_angle",
            ),
        )
    )
    row.update(_stats("edge_length", _tetra_edge_lengths(mesh)))
    return row
