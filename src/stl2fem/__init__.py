"""STL-to-FEM meshing helpers for micromagnetic modeling."""

from .conversion import tetrahedralize_stl_for_merrill, tetrahedralize_stl_with_gmsh
from .datasets import assign_size_bins, detect_stl_format, nikolaisen_inventory
from .memory import estimate_merrill_memory
from .quality import inspect_surface_stl, inspect_volume_mesh
from .repair import repair_surface_with_pymeshfix
from .units import DEFAULT_TARGET_EDGE_LENGTH_M, make_unit_context, scale_for_unit
from .workflow import process_nikolaisen_size_bin

__all__ = [
    "assign_size_bins",
    "DEFAULT_TARGET_EDGE_LENGTH_M",
    "detect_stl_format",
    "estimate_merrill_memory",
    "inspect_surface_stl",
    "inspect_volume_mesh",
    "nikolaisen_inventory",
    "process_nikolaisen_size_bin",
    "repair_surface_with_pymeshfix",
    "make_unit_context",
    "scale_for_unit",
    "tetrahedralize_stl_for_merrill",
    "tetrahedralize_stl_with_gmsh",
]
