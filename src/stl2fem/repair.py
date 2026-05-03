"""Surface repair helpers."""

from __future__ import annotations

from pathlib import Path

from .quality import inspect_surface_stl, load_surface


def repair_surface_with_pymeshfix(
    input_path: str | Path,
    output_path: str | Path,
    *,
    overwrite: bool = False,
    verbose: bool = False,
) -> dict[str, object]:
    """Repair a triangle surface with PyMeshFix and write a cleaned STL.

    MeshFix assumes the input represents one solid object. Use this for small
    defects, and manually review large holes or self-intersections before
    accepting the repaired geometry for micromagnetic modeling.
    """

    input_path = Path(input_path)
    output_path = Path(output_path)
    if output_path.exists() and not overwrite:
        return inspect_surface_stl(output_path)

    try:
        import pymeshfix
    except ImportError as exc:
        raise ImportError(
            "PyMeshFix is required for surface repair. Install with `pip install pymeshfix`."
        ) from exc

    surface = load_surface(input_path)
    faces = surface.faces.reshape(-1, 4)[:, 1:]
    meshfix = pymeshfix.MeshFix(surface.points, faces)
    try:
        meshfix.repair(verbose=verbose)
    except TypeError:
        meshfix.repair()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    meshfix.mesh.save(output_path)
    return inspect_surface_stl(output_path)
