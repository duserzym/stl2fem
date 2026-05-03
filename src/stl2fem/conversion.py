"""STL to tetrahedral mesh conversion."""

from __future__ import annotations

from pathlib import Path
import math

import numpy as np

from .quality import load_surface
from .repair import repair_surface_with_pymeshfix
from .units import DEFAULT_TARGET_EDGE_LENGTH_M, make_unit_context


GMSH_3D_ALGORITHMS = {
    "delaunay": 1,
    "frontal": 4,
    "hxt": 10,
}

DEFAULT_TARGET_EDGE_LENGTH = 0.009
DEFAULT_BRUTE_FORCE_MAX_CELLS_PER_AXIS = 220


def _require_meshio():
    try:
        import meshio
    except ImportError as exc:
        raise ImportError("meshio is required to write Gmsh meshes.") from exc
    return meshio


def _require_pyvista():
    try:
        import pyvista as pv
    except ImportError as exc:
        raise ImportError(
            "PyVista is required for brute-force fallback meshing. "
            "Install with `pip install -e .[examples]` or `pip install pyvista`."
        ) from exc
    return pv


def tetrahedralize_stl_with_gmsh(
    stl_path: str | Path,
    msh_path: str | Path,
    *,
    algorithm: str = "delaunay",
    target_edge_length: float | None = DEFAULT_TARGET_EDGE_LENGTH,
    characteristic_length: float | None = None,
    classify_angle_degrees: float = 40.0,
    msh_version: str = "2.2",
    optimize: bool = True,
    optimize_method: str = "Gmsh",
    overwrite: bool = False,
    verbose: bool = False,
) -> Path:
    """Convert an STL surface to a tetrahedral Gmsh ``.msh`` volume mesh.

    ``target_edge_length`` is expressed in the same coordinate units as the STL.
    For the Nikolaisen2022 micrometer-scale particle STLs, the default ``0.009``
    corresponds to 9 nm, close to the magnetite exchange length. Gmsh may not
    realize this target exactly; inspect the converted mesh edge-length stats.
    """

    stl_path = Path(stl_path)
    msh_path = Path(msh_path)
    if msh_path.exists() and not overwrite:
        return msh_path

    algorithm_key = algorithm.lower()
    if algorithm_key not in GMSH_3D_ALGORITHMS:
        known = ", ".join(sorted(GMSH_3D_ALGORITHMS))
        raise ValueError(f"Unknown Gmsh 3D algorithm {algorithm!r}; expected {known}")
    if characteristic_length is not None:
        target_edge_length = characteristic_length

    try:
        import gmsh
    except ImportError as exc:
        raise ImportError(
            "The gmsh Python package is required for tetrahedralization. "
            "Install with `pip install gmsh` or `pip install -e .`."
        ) from exc

    msh_path.parent.mkdir(parents=True, exist_ok=True)
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", GMSH_3D_ALGORITHMS[algorithm_key])
        gmsh.option.setNumber("Mesh.MshFileVersion", float(msh_version))
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        if target_edge_length is not None:
            gmsh.option.setNumber("Mesh.MeshSizeMax", target_edge_length)
            gmsh.option.setNumber("Mesh.MeshSizeMin", target_edge_length)
            gmsh.option.setNumber("Mesh.CharacteristicLengthMax", target_edge_length)
            gmsh.option.setNumber("Mesh.CharacteristicLengthMin", target_edge_length)

        gmsh.model.add(stl_path.stem)
        gmsh.merge(str(stl_path))
        gmsh.model.mesh.classifySurfaces(
            math.radians(classify_angle_degrees),
            True,
            True,
            math.pi,
        )
        gmsh.model.mesh.createGeometry()

        surface_tags = [tag for dim, tag in gmsh.model.getEntities(2)]
        if not surface_tags:
            raise RuntimeError(f"Gmsh could not classify any surfaces from {stl_path}")

        surface_loop = gmsh.model.geo.addSurfaceLoop(surface_tags)
        volume = gmsh.model.geo.addVolume([surface_loop])
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(3, [volume], 1)
        gmsh.model.setPhysicalName(3, 1, "particle")

        gmsh.model.mesh.generate(3)
        if optimize:
            gmsh.model.mesh.optimize(optimize_method)
        gmsh.write(str(msh_path))
    finally:
        gmsh.finalize()

    return msh_path


def write_tetrahedral_grid_as_msh(
    grid,
    output_mesh_path: str | Path,
    *,
    coordinate_scale: float = 1.0,
    overwrite: bool = False,
) -> Path:
    """Write a PyVista tetrahedral grid as a Gmsh 2.2 ``.msh`` file.

    This is used by brute-force fallbacks that build a tetrahedral grid without
    going through Gmsh's STL surface parametrization. The output contains one
    physical/geometrical volume group named by integer tag ``1``.
    """

    pv = _require_pyvista()
    meshio = _require_meshio()

    output_mesh_path = Path(output_mesh_path)
    if output_mesh_path.exists() and not overwrite:
        return output_mesh_path

    if not hasattr(grid, "celltypes"):
        raise ValueError("Expected an unstructured grid with tetrahedral cells.")

    tet_mask = grid.celltypes == pv.CellType.TETRA
    if not tet_mask.any():
        raise ValueError("No tetrahedral cells were generated.")

    tet_grid = grid.extract_cells(tet_mask)
    cells = np.asarray(tet_grid.cells)
    if cells.size % 5 != 0:
        raise ValueError("Unexpected tetrahedral cell array layout.")
    tets = cells.reshape(-1, 5)
    if not np.all(tets[:, 0] == 4):
        raise ValueError("Unexpected non-tetrahedral cells in extracted grid.")

    tetra_cells = np.asarray(tets[:, 1:5], dtype=np.int64)
    points = np.asarray(tet_grid.points, dtype=float) * coordinate_scale
    tags = np.ones(len(tetra_cells), dtype=np.int32)
    mesh = meshio.Mesh(
        points=points,
        cells=[("tetra", tetra_cells)],
        cell_data={"gmsh:physical": [tags], "gmsh:geometrical": [tags.copy()]},
    )

    output_mesh_path.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(output_mesh_path, mesh, file_format="gmsh22")
    return output_mesh_path


def _capped_voxel_spacing(
    surface,
    target_edge_length: float,
    max_cells_per_axis: int | None,
) -> float:
    bounds = surface.bounds
    lengths = np.asarray(
        [bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]],
        dtype=float,
    )
    spacing = float(target_edge_length)
    if max_cells_per_axis is not None and max_cells_per_axis > 0:
        spacing = max(spacing, float(np.max(lengths) / max_cells_per_axis))
    if not np.isfinite(spacing) or spacing <= 0:
        raise ValueError("Could not determine a positive voxel spacing.")
    return spacing


def tetrahedralize_stl_with_voxels(
    stl_path: str | Path,
    msh_path: str | Path,
    *,
    target_edge_length: float = DEFAULT_TARGET_EDGE_LENGTH,
    max_cells_per_axis: int | None = DEFAULT_BRUTE_FORCE_MAX_CELLS_PER_AXIS,
    overwrite: bool = False,
) -> tuple[Path, dict[str, object]]:
    """Approximate an STL volume by voxels and split voxels into tetrahedra.

    This fallback is intentionally blunt. It ignores many surface-topology
    defects that stop Gmsh, keeps the rough particle envelope, and produces a
    tetrahedral mesh suitable for ensemble-level screening. If the requested
    edge size would create too many voxels along the longest bounding-box axis,
    spacing is increased to respect ``max_cells_per_axis``.
    """

    stl_path = Path(stl_path)
    msh_path = Path(msh_path)
    if msh_path.exists() and not overwrite:
        return msh_path, {"mesh_strategy": "voxel", "voxel_spacing_native": np.nan}

    surface = load_surface(stl_path)
    spacing = _capped_voxel_spacing(surface, target_edge_length, max_cells_per_axis)
    voxel_grid = surface.voxelize(spacing=spacing)
    tetra_grid = voxel_grid.triangulate()
    write_tetrahedral_grid_as_msh(tetra_grid, msh_path, overwrite=True)

    return msh_path, {
        "mesh_strategy": "voxel",
        "voxel_spacing_native": spacing,
        "voxel_cells": int(voxel_grid.n_cells),
        "voxel_points": int(voxel_grid.n_points),
    }


def tetrahedralize_stl_with_delaunay_hull(
    stl_path: str | Path,
    msh_path: str | Path,
    *,
    overwrite: bool = False,
) -> tuple[Path, dict[str, object]]:
    """Create a tetrahedral mesh from the STL point cloud's Delaunay hull.

    This is the final, most permissive fallback. It preserves the broad
    bounding geometry but can smooth concavities because it meshes the point
    cloud's 3D Delaunay hull rather than the exact surface.
    """

    stl_path = Path(stl_path)
    msh_path = Path(msh_path)
    if msh_path.exists() and not overwrite:
        return msh_path, {"mesh_strategy": "hull_delaunay"}

    surface = load_surface(stl_path)
    grid = surface.delaunay_3d(alpha=0.0, tol=0.001, offset=2.5)
    write_tetrahedral_grid_as_msh(grid, msh_path, overwrite=True)
    return msh_path, {"mesh_strategy": "hull_delaunay"}


def write_meter_scaled_mesh(
    input_mesh_path: str | Path,
    output_mesh_path: str | Path,
    *,
    input_scale_to_meters: float,
    overwrite: bool = False,
) -> Path:
    """Write a copy of a mesh with point coordinates scaled to meters."""

    input_mesh_path = Path(input_mesh_path)
    output_mesh_path = Path(output_mesh_path)
    if output_mesh_path.exists() and not overwrite:
        return output_mesh_path

    meshio = _require_meshio()

    mesh = meshio.read(input_mesh_path)
    mesh.points = np.asarray(mesh.points, dtype=float) * input_scale_to_meters
    output_mesh_path.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(output_mesh_path, mesh, file_format="gmsh22")
    return output_mesh_path


def tetrahedralize_stl_for_merrill(
    stl_path: str | Path,
    native_msh_path: str | Path,
    meter_msh_path: str | Path | None = None,
    *,
    input_unit: str = "um",
    input_scale_to_meters: float | None = None,
    target_edge_length_m: float = DEFAULT_TARGET_EDGE_LENGTH_M,
    target_edge_length_native: float | None = None,
    algorithm: str = "delaunay",
    overwrite: bool = False,
) -> tuple[Path, Path | None, object]:
    """Mesh an STL in native units and optionally export a meter-scale copy.

    STL files do not store units, so callers must declare ``input_unit`` or
    ``input_scale_to_meters``. The native Gmsh mesh is useful for diagnostics;
    the meter-scaled mesh is the one intended for Merrill.jl.
    """

    units = make_unit_context(
        input_unit=input_unit,
        input_scale_to_meters=input_scale_to_meters,
        target_edge_length_m=target_edge_length_m,
        target_edge_length_native=target_edge_length_native,
    )
    native_path = tetrahedralize_stl_with_gmsh(
        stl_path,
        native_msh_path,
        algorithm=algorithm,
        target_edge_length=units.target_edge_length_native,
        overwrite=overwrite,
    )
    meter_path = None
    if meter_msh_path is not None:
        meter_path = write_meter_scaled_mesh(
            native_path,
            meter_msh_path,
            input_scale_to_meters=units.input_scale_to_meters,
            overwrite=overwrite,
        )
    return native_path, meter_path, units


def tetrahedralize_repaired_stl_for_merrill(
    stl_path: str | Path,
    native_msh_path: str | Path,
    repair_stl_path: str | Path,
    meter_msh_path: str | Path | None = None,
    *,
    input_unit: str = "um",
    input_scale_to_meters: float | None = None,
    target_edge_length_m: float = DEFAULT_TARGET_EDGE_LENGTH_M,
    target_edge_length_native: float | None = None,
    algorithm: str = "delaunay",
    overwrite: bool = False,
) -> tuple[Path, Path | None, object, dict[str, object]]:
    """Repair an STL with MeshFix, then run the standard Gmsh conversion."""

    repair_surface_with_pymeshfix(stl_path, repair_stl_path, overwrite=overwrite)
    native_path, meter_path, units = tetrahedralize_stl_for_merrill(
        repair_stl_path,
        native_msh_path,
        meter_msh_path,
        input_unit=input_unit,
        input_scale_to_meters=input_scale_to_meters,
        target_edge_length_m=target_edge_length_m,
        target_edge_length_native=target_edge_length_native,
        algorithm=algorithm,
        overwrite=overwrite,
    )
    return native_path, meter_path, units, {"mesh_strategy": "pymeshfix_gmsh"}


def tetrahedralize_bruteforce_stl_for_merrill(
    stl_path: str | Path,
    native_msh_path: str | Path,
    meter_msh_path: str | Path | None = None,
    *,
    input_unit: str = "um",
    input_scale_to_meters: float | None = None,
    target_edge_length_m: float = DEFAULT_TARGET_EDGE_LENGTH_M,
    target_edge_length_native: float | None = None,
    strategy: str = "voxel",
    max_cells_per_axis: int | None = DEFAULT_BRUTE_FORCE_MAX_CELLS_PER_AXIS,
    overwrite: bool = False,
) -> tuple[Path, Path | None, object, dict[str, object]]:
    """Run a permissive non-Gmsh fallback and optionally export meters.

    Supported strategies are ``"voxel"`` and ``"hull_delaunay"``.
    """

    units = make_unit_context(
        input_unit=input_unit,
        input_scale_to_meters=input_scale_to_meters,
        target_edge_length_m=target_edge_length_m,
        target_edge_length_native=target_edge_length_native,
    )

    strategy_key = strategy.lower()
    if strategy_key == "voxel":
        native_path, info = tetrahedralize_stl_with_voxels(
            stl_path,
            native_msh_path,
            target_edge_length=units.target_edge_length_native,
            max_cells_per_axis=max_cells_per_axis,
            overwrite=overwrite,
        )
    elif strategy_key == "hull_delaunay":
        native_path, info = tetrahedralize_stl_with_delaunay_hull(
            stl_path,
            native_msh_path,
            overwrite=overwrite,
        )
    else:
        raise ValueError("strategy must be 'voxel' or 'hull_delaunay'")

    meter_path = None
    if meter_msh_path is not None:
        write_meter_scaled_mesh(
            native_path,
            meter_msh_path,
            input_scale_to_meters=units.input_scale_to_meters,
            overwrite=overwrite,
        )
        meter_path = Path(meter_msh_path)

    return native_path, meter_path, units, info
