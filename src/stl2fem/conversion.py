"""STL to tetrahedral mesh conversion."""

from __future__ import annotations

from pathlib import Path
import math

import numpy as np

from .units import DEFAULT_TARGET_EDGE_LENGTH_M, make_unit_context


GMSH_3D_ALGORITHMS = {
    "delaunay": 1,
    "frontal": 4,
    "hxt": 10,
}

DEFAULT_TARGET_EDGE_LENGTH = 0.009


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

    try:
        import meshio
    except ImportError as exc:
        raise ImportError("meshio is required to write meter-scaled meshes.") from exc

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
