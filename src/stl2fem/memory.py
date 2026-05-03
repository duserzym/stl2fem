"""Order-of-magnitude memory estimates for downstream micromagnetic modeling."""

from __future__ import annotations


def human_bytes(num_bytes: float) -> str:
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    value = float(num_bytes)
    for unit in units:
        if abs(value) < 1024 or unit == units[-1]:
            return f"{value:.2f} {unit}"
        value /= 1024
    return f"{value:.2f} TiB"


def estimate_merrill_memory(
    n_nodes: int,
    n_tets: int,
    *,
    bytes_per_float: int = 8,
    bytes_per_index: int = 8,
    vector_work_arrays: int = 12,
    scalar_work_arrays: int = 12,
    nnz_per_node: int = 80,
    sparse_operator_count: int = 3,
    safety_factor: float = 1.5,
) -> dict[str, float | str | int]:
    """Estimate memory for screening Merrill.jl runs.

    This is intentionally conservative and should be calibrated against actual
    Merrill.jl runs. It estimates mesh storage, magnetization/work arrays, and
    sparse FEM-like operators; it does not claim to model every allocation in
    Merrill.jl.
    """

    n_nodes = int(n_nodes)
    n_tets = int(n_tets)

    mesh_bytes = (
        n_nodes * 3 * bytes_per_float + n_tets * 4 * bytes_per_index
    )
    vector_state_bytes = n_nodes * 3 * bytes_per_float * vector_work_arrays
    scalar_state_bytes = n_nodes * bytes_per_float * scalar_work_arrays
    tet_workspace_bytes = n_tets * bytes_per_float * 8
    sparse_operator_bytes = (
        n_nodes
        * nnz_per_node
        * sparse_operator_count
        * (bytes_per_float + bytes_per_index)
    )
    estimated_bytes = safety_factor * (
        mesh_bytes
        + vector_state_bytes
        + scalar_state_bytes
        + tet_workspace_bytes
        + sparse_operator_bytes
    )

    return {
        "n_nodes": n_nodes,
        "n_tets": n_tets,
        "mesh_storage_bytes": mesh_bytes,
        "vector_state_bytes": vector_state_bytes,
        "scalar_state_bytes": scalar_state_bytes,
        "tet_workspace_bytes": tet_workspace_bytes,
        "sparse_operator_bytes": sparse_operator_bytes,
        "estimated_memory_bytes": estimated_bytes,
        "estimated_memory_mib": estimated_bytes / 1024**2,
        "estimated_memory_gib": estimated_bytes / 1024**3,
        "estimated_memory_human": human_bytes(estimated_bytes),
        "estimate_note": "Screening estimate; calibrate with real Merrill.jl runs.",
    }

