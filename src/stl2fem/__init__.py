"""STL-to-FEM meshing helpers for micromagnetic modeling."""

from .conversion import (
    DEFAULT_BRUTE_FORCE_MAX_CELLS_PER_AXIS,
    tetrahedralize_bruteforce_stl_for_merrill,
    tetrahedralize_repaired_stl_for_merrill,
    tetrahedralize_stl_for_merrill,
    tetrahedralize_stl_with_delaunay_hull,
    tetrahedralize_stl_with_gmsh,
    tetrahedralize_stl_with_voxels,
    write_tetrahedral_grid_as_msh,
)
from .datasets import (
    add_nikolaisen_stl_metadata,
    assign_size_bins,
    detect_stl_format,
    nikolaisen_inventory,
    nikolaisen_stl_metadata,
)
from .memory import estimate_merrill_memory
from .quality import estimate_grain_size_stl, inspect_surface_stl, inspect_volume_mesh
from .repair import repair_surface_with_pymeshfix
from .units import DEFAULT_TARGET_EDGE_LENGTH_M, make_unit_context, scale_for_unit
from .workflow import (
    DEFAULT_MESH_TIMEOUT_SECONDS,
    DEFAULT_MESH_STRATEGIES,
    process_nikolaisen_merrill_meshes,
    process_nikolaisen_size_bin,
)

__all__ = [
    "assign_size_bins",
    "add_nikolaisen_stl_metadata",
    "DEFAULT_BRUTE_FORCE_MAX_CELLS_PER_AXIS",
    "DEFAULT_TARGET_EDGE_LENGTH_M",
    "DEFAULT_MESH_TIMEOUT_SECONDS",
    "DEFAULT_MESH_STRATEGIES",
    "detect_stl_format",
    "estimate_merrill_memory",
    "estimate_grain_size_stl",
    "inspect_surface_stl",
    "inspect_volume_mesh",
    "nikolaisen_inventory",
    "nikolaisen_stl_metadata",
    "process_nikolaisen_size_bin",
    "process_nikolaisen_merrill_meshes",
    "repair_surface_with_pymeshfix",
    "make_unit_context",
    "scale_for_unit",
    "tetrahedralize_bruteforce_stl_for_merrill",
    "tetrahedralize_repaired_stl_for_merrill",
    "tetrahedralize_stl_for_merrill",
    "tetrahedralize_stl_with_delaunay_hull",
    "tetrahedralize_stl_with_gmsh",
    "tetrahedralize_stl_with_voxels",
    "write_tetrahedral_grid_as_msh",
]
