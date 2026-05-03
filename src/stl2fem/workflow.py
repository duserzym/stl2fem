"""End-to-end workflows used by notebooks and the CLI."""

from __future__ import annotations

from pathlib import Path
import traceback

import pandas as pd

from .conversion import tetrahedralize_stl_for_merrill
from .datasets import assign_size_bins, nikolaisen_inventory
from .memory import estimate_merrill_memory
from .quality import inspect_surface_stl, inspect_volume_mesh
from .units import DEFAULT_TARGET_EDGE_LENGTH_M, add_meter_scaled_columns, make_unit_context


def process_nikolaisen_size_bin(
    dataset_root: str | Path = "data/Nikolaisen2022",
    output_root: str | Path = "processed/Nikolaisen2022",
    *,
    bin_index: int = 0,
    n_bins: int = 4,
    phases: tuple[str, ...] = ("OPX", "PLAG"),
    input_unit: str = "um",
    input_scale_to_meters: float | None = None,
    target_edge_length_m: float = DEFAULT_TARGET_EDGE_LENGTH_M,
    target_edge_length: float | None = None,
    characteristic_length: float | None = None,
    overwrite: bool = False,
    max_meshes: int | None = None,
    continue_on_error: bool = True,
) -> pd.DataFrame:
    """Quality-check, convert, and estimate memory for one file-size bin."""

    if characteristic_length is not None:
        target_edge_length = characteristic_length
    units = make_unit_context(
        input_unit=input_unit,
        input_scale_to_meters=input_scale_to_meters,
        target_edge_length_m=target_edge_length_m,
        target_edge_length_native=target_edge_length,
    )

    output_root = Path(output_root)
    reports_dir = output_root / "04_quality_reports"
    msh_dir = output_root / "03_tet_msh"
    merrill_dir = output_root / "05_merrill_ready"
    reports_dir.mkdir(parents=True, exist_ok=True)
    msh_dir.mkdir(parents=True, exist_ok=True)
    merrill_dir.mkdir(parents=True, exist_ok=True)

    inventory = assign_size_bins(
        nikolaisen_inventory(dataset_root, phases=phases),
        n_bins=n_bins,
    )
    subset = inventory[inventory["size_bin_index"] == bin_index].copy()
    if max_meshes is not None:
        subset = subset.head(max_meshes).copy()

    rows: list[dict[str, object]] = []
    for _, record in subset.iterrows():
        particle_id = record["particle_id"]
        stl_path = Path(record["source_path"])
        msh_path = msh_dir / f"{particle_id}.msh"
        merrill_msh_path = merrill_dir / f"{particle_id}.msh"
        row = record.to_dict()
        row["input_unit"] = units.input_unit
        row["input_scale_to_meters"] = units.input_scale_to_meters
        row["target_edge_length_native"] = units.target_edge_length_native
        row["target_edge_length_m"] = units.target_edge_length_m

        try:
            row.update(
                {
                    f"surface_{key}": value
                    for key, value in inspect_surface_stl(stl_path).items()
                }
            )
            tetrahedralize_stl_for_merrill(
                stl_path,
                msh_path,
                merrill_msh_path,
                algorithm="delaunay",
                input_unit=units.input_unit,
                input_scale_to_meters=units.input_scale_to_meters,
                target_edge_length_m=units.target_edge_length_m,
                overwrite=overwrite,
            )
            row["msh_path"] = str(msh_path)
            row["merrill_msh_path"] = str(merrill_msh_path)
            row["merrill_msh_size_bytes"] = merrill_msh_path.stat().st_size
            row["merrill_msh_size_mib"] = merrill_msh_path.stat().st_size / 1024**2
            volume_quality = add_meter_scaled_columns(
                inspect_volume_mesh(msh_path),
                input_scale_to_meters=units.input_scale_to_meters,
            )
            row.update(volume_quality)
            row.update(estimate_merrill_memory(row["n_nodes"], row["n_tets"]))
            row["status"] = "ok"
        except Exception as exc:
            row["status"] = "failed"
            row["error_type"] = exc.__class__.__name__
            row["error"] = str(exc)
            row["traceback"] = traceback.format_exc()
            if not continue_on_error:
                raise
        rows.append(row)

        partial = pd.DataFrame(rows)
        partial.to_csv(reports_dir / f"quality_memory_bin_{bin_index:02d}.csv", index=False)

    return pd.DataFrame(rows)
