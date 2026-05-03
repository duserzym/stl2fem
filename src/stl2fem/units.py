"""Physical unit helpers for STL and Merrill.jl mesh workflows.

STL coordinates do not carry standardized units. Callers must declare how to
interpret the input coordinate numbers before we can make Merrill.jl-compatible
meter-scale meshes.
"""

from __future__ import annotations

from dataclasses import dataclass


UNIT_TO_METERS = {
    "m": 1.0,
    "meter": 1.0,
    "meters": 1.0,
    "mm": 1e-3,
    "millimeter": 1e-3,
    "millimeters": 1e-3,
    "um": 1e-6,
    "micrometer": 1e-6,
    "micrometers": 1e-6,
    "micron": 1e-6,
    "microns": 1e-6,
    "nm": 1e-9,
    "nanometer": 1e-9,
    "nanometers": 1e-9,
}

DEFAULT_INPUT_UNIT = "um"
DEFAULT_INPUT_SCALE_TO_METERS = UNIT_TO_METERS[DEFAULT_INPUT_UNIT]
DEFAULT_TARGET_EDGE_LENGTH_M = 9e-9


@dataclass(frozen=True)
class UnitContext:
    input_unit: str
    input_scale_to_meters: float
    target_edge_length_m: float

    @property
    def target_edge_length_native(self) -> float:
        return self.target_edge_length_m / self.input_scale_to_meters


def scale_for_unit(unit: str) -> float:
    """Return the number of meters represented by one input coordinate unit."""

    key = unit.strip().lower().replace("µ", "u").replace("μ", "u")
    try:
        return UNIT_TO_METERS[key]
    except KeyError as exc:
        known = ", ".join(sorted(UNIT_TO_METERS))
        raise ValueError(f"Unknown unit {unit!r}; expected one of {known}") from exc


def make_unit_context(
    *,
    input_unit: str = DEFAULT_INPUT_UNIT,
    input_scale_to_meters: float | None = None,
    target_edge_length_m: float = DEFAULT_TARGET_EDGE_LENGTH_M,
    target_edge_length_native: float | None = None,
) -> UnitContext:
    """Create a consistent unit context for meshing and reporting."""

    scale = scale_for_unit(input_unit) if input_scale_to_meters is None else input_scale_to_meters
    if scale <= 0:
        raise ValueError("input_scale_to_meters must be positive")

    target_m = target_edge_length_m
    if target_edge_length_native is not None:
        target_m = target_edge_length_native * scale
    if target_m <= 0:
        raise ValueError("target edge length must be positive")

    return UnitContext(
        input_unit=input_unit,
        input_scale_to_meters=float(scale),
        target_edge_length_m=float(target_m),
    )


def add_meter_scaled_columns(
    row: dict[str, object],
    *,
    input_scale_to_meters: float,
    prefixes: tuple[str, ...] = ("edge_length", "bbox"),
) -> dict[str, object]:
    """Add meter-scaled length columns to a quality/report row."""

    result = dict(row)
    if "edge_length" in prefixes:
        for key in ("edge_length_min", "edge_length_median", "edge_length_p95", "edge_length_max"):
            if key in row:
                result[f"{key}_m"] = float(row[key]) * input_scale_to_meters
    if "bbox" in prefixes:
        for axis in ("x", "y", "z"):
            for suffix in ("min", "max"):
                key = f"bbox_{axis}{suffix}"
                if key in row:
                    result[f"{key}_m"] = float(row[key]) * input_scale_to_meters
    return result

