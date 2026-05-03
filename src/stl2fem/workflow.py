"""End-to-end workflows used by notebooks and the CLI."""

from __future__ import annotations

import multiprocessing as mp
from pathlib import Path
from queue import Empty
import traceback

import pandas as pd

from .conversion import (
    DEFAULT_BRUTE_FORCE_MAX_CELLS_PER_AXIS,
    tetrahedralize_bruteforce_stl_for_merrill,
    tetrahedralize_repaired_stl_for_merrill,
    tetrahedralize_stl_for_merrill,
)
from .datasets import add_nikolaisen_stl_metadata, assign_size_bins, nikolaisen_inventory
from .memory import estimate_merrill_memory
from .quality import inspect_surface_stl, inspect_volume_mesh
from .units import DEFAULT_TARGET_EDGE_LENGTH_M, add_meter_scaled_columns, make_unit_context


DEFAULT_NIKOLAISEN_DATASET_ROOT = "data/Nikolaisen2022"
DEFAULT_NIKOLAISEN_MERRILL_OUTPUT_ROOT = "data/Nikolaisen2022_merrill_msh"
DEFAULT_NIKOLAISEN_PHASES = ("PLAG",)
DEFAULT_MAX_EVSD_UM = 1.0
DEFAULT_MESH_TIMEOUT_SECONDS = 60.0
DEFAULT_MESH_STRATEGIES = ("gmsh", "pymeshfix_gmsh", "voxel", "hull_delaunay")


def _cleanup_mesh_outputs(*paths: str | Path | None) -> None:
    for path in paths:
        if path is None:
            continue
        mesh_path = Path(path)
        if mesh_path.exists():
            mesh_path.unlink()


def _convert_and_inspect_mesh(payload: dict[str, object]) -> dict[str, object]:
    stl_path = Path(payload["stl_path"])
    msh_path = Path(payload["msh_path"])
    merrill_msh_path = Path(payload["merrill_msh_path"])
    strategy = str(payload.get("mesh_strategy", "gmsh"))
    info: dict[str, object] = {"mesh_strategy": strategy}

    if strategy == "gmsh":
        tetrahedralize_stl_for_merrill(
            stl_path,
            msh_path,
            merrill_msh_path,
            algorithm="delaunay",
            input_unit=str(payload["input_unit"]),
            input_scale_to_meters=float(payload["input_scale_to_meters"]),
            target_edge_length_m=float(payload["target_edge_length_m"]),
            overwrite=bool(payload["overwrite"]),
        )
    elif strategy == "pymeshfix_gmsh":
        _, _, _, info = tetrahedralize_repaired_stl_for_merrill(
            stl_path,
            msh_path,
            Path(payload["repair_stl_path"]),
            merrill_msh_path,
            algorithm="delaunay",
            input_unit=str(payload["input_unit"]),
            input_scale_to_meters=float(payload["input_scale_to_meters"]),
            target_edge_length_m=float(payload["target_edge_length_m"]),
            overwrite=bool(payload["overwrite"]),
        )
    elif strategy in {"voxel", "hull_delaunay"}:
        _, _, _, info = tetrahedralize_bruteforce_stl_for_merrill(
            stl_path,
            msh_path,
            merrill_msh_path,
            input_unit=str(payload["input_unit"]),
            input_scale_to_meters=float(payload["input_scale_to_meters"]),
            target_edge_length_m=float(payload["target_edge_length_m"]),
            strategy=strategy,
            max_cells_per_axis=int(payload["bruteforce_max_cells_per_axis"]),
            overwrite=bool(payload["overwrite"]),
        )
    else:
        raise ValueError(f"Unknown mesh strategy {strategy!r}")

    result: dict[str, object] = {
        "msh_path": str(msh_path),
        "merrill_msh_path": str(merrill_msh_path),
        "merrill_msh_size_bytes": merrill_msh_path.stat().st_size,
        "merrill_msh_size_mib": merrill_msh_path.stat().st_size / 1024**2,
    }
    result.update(info)
    volume_quality = add_meter_scaled_columns(
        inspect_volume_mesh(msh_path),
        input_scale_to_meters=float(payload["input_scale_to_meters"]),
    )
    result.update(volume_quality)
    result.update(estimate_merrill_memory(result["n_nodes"], result["n_tets"]))
    result["status"] = "ok"

    if not bool(payload["save_native_mesh"]):
        _cleanup_mesh_outputs(msh_path)
        result["msh_path"] = None
    return result


def _convert_and_inspect_mesh_worker(
    payload: dict[str, object],
    queue: mp.Queue,
) -> None:
    try:
        queue.put(("ok", _convert_and_inspect_mesh(payload)))
    except BaseException as exc:
        queue.put(
            (
                "failed",
                {
                    "status": "failed",
                    "error_type": exc.__class__.__name__,
                    "error": str(exc),
                    "traceback": traceback.format_exc(),
                },
            )
        )


def _run_conversion_with_timeout(
    payload: dict[str, object],
    mesh_timeout_seconds: float | None,
) -> dict[str, object]:
    if mesh_timeout_seconds is None:
        return _convert_and_inspect_mesh(payload)

    ctx = mp.get_context("spawn")
    queue = ctx.Queue()
    process = ctx.Process(target=_convert_and_inspect_mesh_worker, args=(payload, queue))
    process.start()
    process.join(mesh_timeout_seconds)

    if process.is_alive():
        process.terminate()
        process.join(5)
        if process.is_alive():
            process.kill()
            process.join(5)
        if not bool(payload["save_native_mesh"]):
            _cleanup_mesh_outputs(payload["msh_path"], payload["merrill_msh_path"])
        return {
            "status": "failed",
            "error_type": "TimeoutError",
            "error": (
                f"Mesh conversion exceeded {mesh_timeout_seconds:g} seconds. "
                "Increase mesh_timeout_seconds to retry this STL."
            ),
            "traceback": "",
        }

    try:
        _, result = queue.get_nowait()
    except Empty:
        result = {
            "status": "failed",
            "error_type": "RuntimeError",
            "error": f"Mesh worker exited with code {process.exitcode} without a result.",
            "traceback": "",
        }

    if result.get("status") != "ok" and not bool(payload["save_native_mesh"]):
        _cleanup_mesh_outputs(payload["msh_path"], payload["merrill_msh_path"])
    return result


def _run_conversion_strategies(
    payload: dict[str, object],
    strategies: tuple[str, ...],
    mesh_timeout_seconds: float | None,
) -> dict[str, object]:
    errors: list[str] = []
    last_result: dict[str, object] | None = None
    for strategy in strategies:
        strategy_payload = {**payload, "mesh_strategy": strategy}
        _cleanup_mesh_outputs(payload["msh_path"], payload["merrill_msh_path"])
        result = _run_conversion_with_timeout(strategy_payload, mesh_timeout_seconds)
        last_result = result
        if result.get("status") == "ok":
            result["attempted_mesh_strategies"] = ";".join(strategies[: strategies.index(strategy) + 1])
            result["strategy_errors"] = " | ".join(errors)
            return result

        errors.append(
            f"{strategy}: {result.get('error_type', 'Error')}: {result.get('error', '')}"
        )

    if last_result is None:
        last_result = {
            "status": "failed",
            "error_type": "RuntimeError",
            "error": "No mesh strategies were attempted.",
            "traceback": "",
        }
    last_result["attempted_mesh_strategies"] = ";".join(strategies)
    last_result["strategy_errors"] = " | ".join(errors)
    return last_result


def process_nikolaisen_size_bin(
    dataset_root: str | Path = DEFAULT_NIKOLAISEN_DATASET_ROOT,
    output_root: str | Path = DEFAULT_NIKOLAISEN_MERRILL_OUTPUT_ROOT,
    *,
    bin_index: int = 0,
    n_bins: int = 4,
    phases: tuple[str, ...] = DEFAULT_NIKOLAISEN_PHASES,
    input_unit: str = "um",
    input_scale_to_meters: float | None = None,
    target_edge_length_m: float = DEFAULT_TARGET_EDGE_LENGTH_M,
    target_edge_length: float | None = None,
    characteristic_length: float | None = None,
    max_evsd_um: float | None = DEFAULT_MAX_EVSD_UM,
    save_native_mesh: bool = False,
    flat_merrill_output: bool = True,
    overwrite: bool = False,
    max_meshes: int | None = None,
    mesh_strategies: tuple[str, ...] = DEFAULT_MESH_STRATEGIES,
    bruteforce_max_cells_per_axis: int = DEFAULT_BRUTE_FORCE_MAX_CELLS_PER_AXIS,
    mesh_timeout_seconds: float | None = DEFAULT_MESH_TIMEOUT_SECONDS,
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
    reports_dir = output_root / "_reports"
    msh_dir = output_root / "_native_tmp"
    repair_dir = output_root / "_repair_tmp"
    merrill_dir = output_root if flat_merrill_output else output_root / "05_merrill_ready"
    reports_dir.mkdir(parents=True, exist_ok=True)
    msh_dir.mkdir(parents=True, exist_ok=True)
    repair_dir.mkdir(parents=True, exist_ok=True)
    merrill_dir.mkdir(parents=True, exist_ok=True)

    inventory = add_nikolaisen_stl_metadata(
        nikolaisen_inventory(dataset_root, phases=phases),
        dataset_root=dataset_root,
    )
    if max_evsd_um is not None:
        inventory = inventory[inventory["metadata_evsd_um"] < max_evsd_um].copy()
    inventory = assign_size_bins(inventory, n_bins=n_bins)
    subset = inventory[inventory["size_bin_index"] == bin_index].copy()
    if max_meshes is not None:
        subset = subset.head(max_meshes).copy()

    rows: list[dict[str, object]] = []
    for _, record in subset.iterrows():
        particle_id = record["particle_id"]
        stl_path = Path(record["source_path"])
        msh_path = msh_dir / f"{particle_id}.msh"
        repair_stl_path = repair_dir / f"{particle_id}_pymeshfix.stl"
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
            row.update(
                _run_conversion_strategies(
                    {
                        "stl_path": str(stl_path),
                        "msh_path": str(msh_path),
                        "repair_stl_path": str(repair_stl_path),
                        "merrill_msh_path": str(merrill_msh_path),
                        "input_unit": units.input_unit,
                        "input_scale_to_meters": units.input_scale_to_meters,
                        "target_edge_length_m": units.target_edge_length_m,
                        "overwrite": overwrite,
                        "save_native_mesh": save_native_mesh,
                        "bruteforce_max_cells_per_axis": bruteforce_max_cells_per_axis,
                    },
                    strategies=mesh_strategies,
                    mesh_timeout_seconds=mesh_timeout_seconds,
                )
            )
            if row.get("status") != "ok" and not continue_on_error:
                raise RuntimeError(str(row.get("error", "mesh conversion failed")))
        except Exception as exc:
            row["status"] = "failed"
            row["error_type"] = exc.__class__.__name__
            row["error"] = str(exc)
            row["traceback"] = traceback.format_exc()
            if not save_native_mesh:
                _cleanup_mesh_outputs(msh_path, merrill_msh_path)
            if not continue_on_error:
                raise
        rows.append(row)

        partial = pd.DataFrame(rows)
        partial.to_csv(reports_dir / f"quality_memory_bin_{bin_index:02d}.csv", index=False)

    if not save_native_mesh and msh_dir.exists():
        for native_mesh in msh_dir.glob("*.msh"):
            native_mesh.unlink()
        try:
            msh_dir.rmdir()
        except OSError:
            pass
    if repair_dir.exists():
        for repaired_stl in repair_dir.glob("*.stl"):
            repaired_stl.unlink()
        try:
            repair_dir.rmdir()
        except OSError:
            pass

    return pd.DataFrame(rows)


def process_nikolaisen_merrill_meshes(
    dataset_root: str | Path = DEFAULT_NIKOLAISEN_DATASET_ROOT,
    output_root: str | Path = DEFAULT_NIKOLAISEN_MERRILL_OUTPUT_ROOT,
    *,
    phases: tuple[str, ...] = DEFAULT_NIKOLAISEN_PHASES,
    input_unit: str = "um",
    input_scale_to_meters: float | None = None,
    target_edge_length_m: float = DEFAULT_TARGET_EDGE_LENGTH_M,
    target_edge_length: float | None = None,
    max_evsd_um: float | None = DEFAULT_MAX_EVSD_UM,
    save_native_mesh: bool = False,
    overwrite: bool = False,
    max_meshes: int | None = None,
    mesh_strategies: tuple[str, ...] = DEFAULT_MESH_STRATEGIES,
    bruteforce_max_cells_per_axis: int = DEFAULT_BRUTE_FORCE_MAX_CELLS_PER_AXIS,
    mesh_timeout_seconds: float | None = DEFAULT_MESH_TIMEOUT_SECONDS,
    continue_on_error: bool = True,
) -> pd.DataFrame:
    """Process all selected Nikolaisen particles into Merrill-ready meter meshes.

    The defaults are intentionally the production target requested for this
    repository: Plag binary STL inputs, EVSD below 1 micrometer, 9 nm target edge
    length, and only Merrill-ready meter-scale ``.msh`` files retained.
    """

    frames = []
    remaining = max_meshes
    for bin_index in range(4):
        frame = process_nikolaisen_size_bin(
            dataset_root,
            output_root,
            bin_index=bin_index,
            n_bins=4,
            phases=phases,
            input_unit=input_unit,
            input_scale_to_meters=input_scale_to_meters,
            target_edge_length_m=target_edge_length_m,
            target_edge_length=target_edge_length,
            max_evsd_um=max_evsd_um,
            save_native_mesh=save_native_mesh,
            flat_merrill_output=True,
            overwrite=overwrite,
            max_meshes=remaining,
            mesh_strategies=mesh_strategies,
            bruteforce_max_cells_per_axis=bruteforce_max_cells_per_axis,
            mesh_timeout_seconds=mesh_timeout_seconds,
            continue_on_error=continue_on_error,
        )
        frames.append(frame)
        if remaining is not None:
            remaining -= len(frame)
            if remaining <= 0:
                break

    if not frames:
        return pd.DataFrame()

    result = pd.concat(frames, ignore_index=True)
    reports_dir = Path(output_root) / "_reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(reports_dir / "quality_memory_all.csv", index=False)
    return result
