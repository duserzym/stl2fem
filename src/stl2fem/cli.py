"""Command line interface."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .conversion import tetrahedralize_stl_for_merrill
from .datasets import assign_size_bins, nikolaisen_inventory
from .memory import estimate_merrill_memory
from .quality import inspect_surface_stl, inspect_volume_mesh
from .units import DEFAULT_TARGET_EDGE_LENGTH_M, add_meter_scaled_columns
from .workflow import (
    DEFAULT_MESH_TIMEOUT_SECONDS,
    DEFAULT_MESH_STRATEGIES,
    DEFAULT_MAX_EVSD_UM,
    DEFAULT_NIKOLAISEN_MERRILL_OUTPUT_ROOT,
    DEFAULT_NIKOLAISEN_PHASES,
    process_nikolaisen_merrill_meshes,
    process_nikolaisen_size_bin,
)


def _float_format(value: float) -> str:
    return f"{value:.12g}"


def _print_series(values: dict[str, object]) -> None:
    print(pd.Series(values).to_string(float_format=_float_format))


def _print_frame(frame: pd.DataFrame) -> None:
    print(frame.to_string(index=False, float_format=_float_format))


def _mesh_timeout(value: float | None) -> float | None:
    if value is not None and value <= 0:
        return None
    return value


def _inventory(args: argparse.Namespace) -> None:
    frame = assign_size_bins(
        nikolaisen_inventory(args.dataset_root, phases=tuple(args.phase or ["OPX", "PLAG"])),
        n_bins=args.bins,
    )
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(out, index=False)
    _print_frame(frame)


def _process_bin(args: argparse.Namespace) -> None:
    frame = process_nikolaisen_size_bin(
        args.dataset_root,
        args.output_root,
        bin_index=args.bin_index,
        n_bins=args.bins,
        phases=tuple(args.phase or DEFAULT_NIKOLAISEN_PHASES),
        input_unit=args.input_unit,
        input_scale_to_meters=args.input_scale_to_meters,
        target_edge_length_m=args.target_edge_length_m,
        target_edge_length=args.target_edge_length_native,
        max_evsd_um=args.max_evsd_um,
        save_native_mesh=args.save_native_mesh,
        overwrite=args.overwrite,
        max_meshes=args.max_meshes,
        mesh_strategies=tuple(args.mesh_strategy or DEFAULT_MESH_STRATEGIES),
        bruteforce_max_cells_per_axis=args.bruteforce_max_cells_per_axis,
        mesh_timeout_seconds=_mesh_timeout(args.mesh_timeout_seconds),
        continue_on_error=not args.stop_on_error,
    )
    _print_frame(frame)


def _process_merrill(args: argparse.Namespace) -> None:
    frame = process_nikolaisen_merrill_meshes(
        args.dataset_root,
        args.output_root,
        phases=tuple(args.phase or DEFAULT_NIKOLAISEN_PHASES),
        input_unit=args.input_unit,
        input_scale_to_meters=args.input_scale_to_meters,
        target_edge_length_m=args.target_edge_length_m,
        target_edge_length=args.target_edge_length_native,
        max_evsd_um=args.max_evsd_um,
        save_native_mesh=args.save_native_mesh,
        overwrite=args.overwrite,
        max_meshes=args.max_meshes,
        mesh_strategies=tuple(args.mesh_strategy or DEFAULT_MESH_STRATEGIES),
        bruteforce_max_cells_per_axis=args.bruteforce_max_cells_per_axis,
        mesh_timeout_seconds=_mesh_timeout(args.mesh_timeout_seconds),
        continue_on_error=not args.stop_on_error,
    )
    _print_frame(frame)


def _inspect_surface(args: argparse.Namespace) -> None:
    _print_series(inspect_surface_stl(args.stl))


def _convert_stl(args: argparse.Namespace) -> None:
    native_path, merrill_path, units = tetrahedralize_stl_for_merrill(
        args.stl,
        args.native_msh,
        args.merrill_msh,
        algorithm=args.algorithm,
        input_unit=args.input_unit,
        input_scale_to_meters=args.input_scale_to_meters,
        target_edge_length_m=args.target_edge_length_m,
        target_edge_length_native=args.target_edge_length_native,
        overwrite=args.overwrite,
    )
    volume = add_meter_scaled_columns(
        inspect_volume_mesh(native_path),
        input_scale_to_meters=units.input_scale_to_meters,
    )
    memory = estimate_merrill_memory(volume["n_nodes"], volume["n_tets"])
    unit_report = {
        "input_unit": units.input_unit,
        "input_scale_to_meters": units.input_scale_to_meters,
        "target_edge_length_native": units.target_edge_length_native,
        "target_edge_length_m": units.target_edge_length_m,
        "native_msh_path": str(native_path),
        "merrill_msh_path": str(merrill_path) if merrill_path else None,
    }
    _print_series({**unit_report, **volume, **memory})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="stl2fem")
    subparsers = parser.add_subparsers(required=True)

    inventory = subparsers.add_parser("nikolaisen-inventory")
    inventory.add_argument("--dataset-root", default="data/Nikolaisen2022")
    inventory.add_argument("--phase", action="append")
    inventory.add_argument("--bins", type=int, default=4)
    inventory.add_argument("--out")
    inventory.set_defaults(func=_inventory)

    process_bin = subparsers.add_parser("process-nikolaisen-bin")
    process_bin.add_argument("--dataset-root", default="data/Nikolaisen2022")
    process_bin.add_argument("--output-root", default=DEFAULT_NIKOLAISEN_MERRILL_OUTPUT_ROOT)
    process_bin.add_argument("--phase", action="append")
    process_bin.add_argument("--bins", type=int, default=4)
    process_bin.add_argument("--bin-index", type=int, required=True)
    process_bin.add_argument("--input-unit", default="um")
    process_bin.add_argument("--input-scale-to-meters", type=float)
    process_bin.add_argument(
        "--target-edge-length-m",
        type=float,
        default=DEFAULT_TARGET_EDGE_LENGTH_M,
        help="Physical target tetrahedron edge length in meters. Default: 9e-9.",
    )
    process_bin.add_argument(
        "--target-edge-length-native",
        "--target-edge-length",
        "--characteristic-length",
        dest="target_edge_length_native",
        type=float,
        default=None,
        help=(
            "Override target edge length in native STL coordinate units. "
            "If omitted, it is derived from --target-edge-length-m and input units."
        ),
    )
    process_bin.add_argument("--max-evsd-um", type=float, default=DEFAULT_MAX_EVSD_UM)
    process_bin.add_argument("--save-native-mesh", action="store_true")
    process_bin.add_argument("--overwrite", action="store_true")
    process_bin.add_argument("--max-meshes", type=int)
    process_bin.add_argument(
        "--mesh-strategy",
        action="append",
        choices=["gmsh", "pymeshfix_gmsh", "voxel", "hull_delaunay"],
        help=(
            "Meshing strategy to try, in order. Repeat to override the default "
            "gmsh -> pymeshfix_gmsh -> voxel -> hull_delaunay sequence."
        ),
    )
    process_bin.add_argument(
        "--bruteforce-max-cells-per-axis",
        type=int,
        default=220,
        help="Coarsen voxel fallback spacing if needed to keep the longest axis under this count.",
    )
    process_bin.add_argument(
        "--mesh-timeout-seconds",
        type=float,
        default=DEFAULT_MESH_TIMEOUT_SECONDS,
        help=(
            "Maximum seconds to spend on one STL before recording a timeout. "
            "Use 0 to disable the timeout."
        ),
    )
    process_bin.add_argument("--stop-on-error", action="store_true")
    process_bin.set_defaults(func=_process_bin)

    process_merrill = subparsers.add_parser("process-nikolaisen-merrill")
    process_merrill.add_argument("--dataset-root", default="data/Nikolaisen2022")
    process_merrill.add_argument("--output-root", default=DEFAULT_NIKOLAISEN_MERRILL_OUTPUT_ROOT)
    process_merrill.add_argument("--phase", action="append")
    process_merrill.add_argument("--input-unit", default="um")
    process_merrill.add_argument("--input-scale-to-meters", type=float)
    process_merrill.add_argument(
        "--target-edge-length-m",
        type=float,
        default=DEFAULT_TARGET_EDGE_LENGTH_M,
        help="Physical target tetrahedron edge length in meters. Default: 9e-9.",
    )
    process_merrill.add_argument(
        "--target-edge-length-native",
        "--target-edge-length",
        "--characteristic-length",
        dest="target_edge_length_native",
        type=float,
        default=None,
        help=(
            "Override target edge length in native STL coordinate units. "
            "If omitted, it is derived from --target-edge-length-m and input units."
        ),
    )
    process_merrill.add_argument("--max-evsd-um", type=float, default=DEFAULT_MAX_EVSD_UM)
    process_merrill.add_argument("--save-native-mesh", action="store_true")
    process_merrill.add_argument("--overwrite", action="store_true")
    process_merrill.add_argument("--max-meshes", type=int)
    process_merrill.add_argument(
        "--mesh-strategy",
        action="append",
        choices=["gmsh", "pymeshfix_gmsh", "voxel", "hull_delaunay"],
        help=(
            "Meshing strategy to try, in order. Repeat to override the default "
            "gmsh -> pymeshfix_gmsh -> voxel -> hull_delaunay sequence."
        ),
    )
    process_merrill.add_argument(
        "--bruteforce-max-cells-per-axis",
        type=int,
        default=220,
        help="Coarsen voxel fallback spacing if needed to keep the longest axis under this count.",
    )
    process_merrill.add_argument(
        "--mesh-timeout-seconds",
        type=float,
        default=DEFAULT_MESH_TIMEOUT_SECONDS,
        help=(
            "Maximum seconds to spend on one STL before recording a timeout. "
            "Use 0 to disable the timeout."
        ),
    )
    process_merrill.add_argument("--stop-on-error", action="store_true")
    process_merrill.set_defaults(func=_process_merrill)

    inspect_surface = subparsers.add_parser("inspect-surface")
    inspect_surface.add_argument("stl")
    inspect_surface.set_defaults(func=_inspect_surface)

    convert_stl = subparsers.add_parser("convert-stl")
    convert_stl.add_argument("stl")
    convert_stl.add_argument("native_msh")
    convert_stl.add_argument("--merrill-msh")
    convert_stl.add_argument("--algorithm", default="delaunay")
    convert_stl.add_argument("--input-unit", default="um")
    convert_stl.add_argument("--input-scale-to-meters", type=float)
    convert_stl.add_argument(
        "--target-edge-length-m",
        type=float,
        default=DEFAULT_TARGET_EDGE_LENGTH_M,
        help="Physical target tetrahedron edge length in meters. Default: 9e-9.",
    )
    convert_stl.add_argument(
        "--target-edge-length-native",
        "--target-edge-length",
        "--characteristic-length",
        dest="target_edge_length_native",
        type=float,
        default=None,
        help=(
            "Override target edge length in native STL coordinate units. "
            "If omitted, it is derived from --target-edge-length-m and input units."
        ),
    )
    convert_stl.add_argument("--overwrite", action="store_true")
    convert_stl.set_defaults(func=_convert_stl)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
