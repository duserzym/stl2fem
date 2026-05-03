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
from .workflow import process_nikolaisen_size_bin


def _float_format(value: float) -> str:
    return f"{value:.12g}"


def _print_series(values: dict[str, object]) -> None:
    print(pd.Series(values).to_string(float_format=_float_format))


def _print_frame(frame: pd.DataFrame) -> None:
    print(frame.to_string(index=False, float_format=_float_format))


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
        phases=tuple(args.phase or ["OPX", "PLAG"]),
        input_unit=args.input_unit,
        input_scale_to_meters=args.input_scale_to_meters,
        target_edge_length_m=args.target_edge_length_m,
        target_edge_length=args.target_edge_length_native,
        overwrite=args.overwrite,
        max_meshes=args.max_meshes,
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
    process_bin.add_argument("--output-root", default="processed/Nikolaisen2022")
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
    process_bin.add_argument("--overwrite", action="store_true")
    process_bin.add_argument("--max-meshes", type=int)
    process_bin.add_argument("--stop-on-error", action="store_true")
    process_bin.set_defaults(func=_process_bin)

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
