"""Dataset discovery helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import struct

import pandas as pd


@dataclass(frozen=True)
class NikolaisenPhase:
    name: str
    folder: str
    prefix: str


NIKOLAISEN_PHASES = {
    "OPX": NikolaisenPhase("OPX", "OPX Binary meshes", "OPX"),
    "PLAG": NikolaisenPhase("PLAG", "Plag Binary meshes", "PLAG"),
}

SIZE_BIN_LABELS = ("00_smallest", "01_small", "02_large", "03_largest")
NIKOLAISEN_METADATA_FILES = {
    "OPX": "OPX_stl_B16.csv",
    "PLAG": "Plag_stl_B16.csv",
}


def detect_stl_format(path: str | Path) -> str:
    """Detect whether an STL is binary or ASCII from file content.

    The Nikolaisen2022 "Binary meshes" folders contain files that may be ASCII
    STL despite the folder name, so examples should use this rather than folder
    labels.
    """

    path = Path(path)
    size = path.stat().st_size
    with path.open("rb") as handle:
        head = handle.read(512)

    if size >= 84:
        with path.open("rb") as handle:
            handle.seek(80)
            raw_count = handle.read(4)
        if len(raw_count) == 4:
            triangle_count = struct.unpack("<I", raw_count)[0]
            if size == 84 + 50 * triangle_count:
                return "binary"

    stripped = head.lstrip().lower()
    if stripped.startswith(b"solid") and (b"facet" in stripped or b"\n" in head):
        return "ascii"

    return "unknown"


def particle_id_from_path(path: str | Path) -> str:
    stem = Path(path).stem
    stem = re.sub(r"[-_](binary|ascii)$", "", stem, flags=re.IGNORECASE)
    return stem.upper()


def nikolaisen_inventory(
    dataset_root: str | Path = "data/Nikolaisen2022",
    phases: tuple[str, ...] = ("OPX", "PLAG"),
) -> pd.DataFrame:
    """Inventory Nikolaisen2022 individual STL files from the binary folders."""

    dataset_root = Path(dataset_root)
    rows: list[dict[str, object]] = []

    for phase_name in phases:
        phase_key = phase_name.upper()
        if phase_key not in NIKOLAISEN_PHASES:
            known = ", ".join(NIKOLAISEN_PHASES)
            raise ValueError(f"Unknown phase {phase_name!r}; expected one of {known}")

        phase = NIKOLAISEN_PHASES[phase_key]
        folder = dataset_root / phase.folder
        for path in sorted(folder.glob("*.stl")):
            size_bytes = path.stat().st_size
            rows.append(
                {
                    "particle_id": particle_id_from_path(path),
                    "phase": phase.name,
                    "source_path": str(path),
                    "source_folder": phase.folder,
                    "stl_size_bytes": size_bytes,
                    "stl_size_mib": size_bytes / 1024**2,
                    "stl_format": detect_stl_format(path),
                }
            )

    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    return frame.sort_values(["stl_size_bytes", "particle_id"]).reset_index(drop=True)


def nikolaisen_stl_metadata(
    dataset_root: str | Path = "data/Nikolaisen2022",
    phases: tuple[str, ...] = ("OPX", "PLAG"),
) -> pd.DataFrame:
    """Load Nikolaisen STL metadata with normalized particle IDs."""

    dataset_root = Path(dataset_root)
    rows: list[pd.DataFrame] = []
    for phase_name in phases:
        phase_key = phase_name.upper()
        metadata_name = NIKOLAISEN_METADATA_FILES.get(phase_key)
        if metadata_name is None:
            known = ", ".join(NIKOLAISEN_METADATA_FILES)
            raise ValueError(f"Unknown phase {phase_name!r}; expected one of {known}")

        path = dataset_root / metadata_name
        frame = pd.read_csv(path)
        frame = frame.rename(
            columns={
                "Filename": "particle_id",
                "Volume": "metadata_volume_um3",
                "EVSD (mu)": "metadata_evsd_um",
            }
        )
        frame["particle_id"] = frame["particle_id"].astype(str).str.upper()
        frame["phase"] = phase_key
        rows.append(frame[["particle_id", "phase", "metadata_volume_um3", "metadata_evsd_um"]])

    if not rows:
        return pd.DataFrame(
            columns=["particle_id", "phase", "metadata_volume_um3", "metadata_evsd_um"]
        )
    return pd.concat(rows, ignore_index=True)


def add_nikolaisen_stl_metadata(
    inventory: pd.DataFrame,
    dataset_root: str | Path = "data/Nikolaisen2022",
) -> pd.DataFrame:
    """Attach Nikolaisen volume and EVSD metadata to an inventory table."""

    if inventory.empty:
        return inventory.copy()

    phases = tuple(sorted(inventory["phase"].dropna().unique()))
    metadata = nikolaisen_stl_metadata(dataset_root, phases=phases)
    return inventory.merge(metadata, on=["particle_id", "phase"], how="left")


def assign_size_bins(df: pd.DataFrame, n_bins: int = 4) -> pd.DataFrame:
    """Sort meshes by STL file size and assign balanced notebook bins."""

    if n_bins < 1:
        raise ValueError("n_bins must be >= 1")
    if df.empty:
        result = df.copy()
        result["size_rank"] = []
        result["size_bin"] = []
        return result

    labels = list(SIZE_BIN_LABELS)
    if n_bins != len(labels):
        labels = [f"{index:02d}_bin" for index in range(n_bins)]

    result = df.sort_values(["stl_size_bytes", "particle_id"]).reset_index(drop=True)
    count = len(result)
    result["size_rank"] = range(count)
    result["size_bin_index"] = [
        min(n_bins - 1, int(rank * n_bins / count)) for rank in range(count)
    ]
    result["size_bin"] = result["size_bin_index"].map(lambda index: labels[index])
    return result
