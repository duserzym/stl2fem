"""PyVista display helpers."""

from __future__ import annotations

from pathlib import Path

from .quality import load_volume_mesh


def display_volume_mesh(
    msh_path: str | Path,
    *,
    scalars: str | None = None,
    show_edges: bool = True,
    opacity: float = 0.35,
    window_size: tuple[int, int] = (900, 650),
    jupyter_backend: str = "trame",
):
    """Display a converted tetrahedral mesh in a notebook or PyVista window."""

    import pyvista as pv

    mesh = load_volume_mesh(msh_path)
    plotter = pv.Plotter(window_size=window_size, notebook=True)
    plotter.add_mesh(mesh, scalars=scalars, show_edges=show_edges, opacity=opacity)
    plotter.add_axes()
    plotter.show_grid()
    return plotter.show(jupyter_backend=jupyter_backend)
