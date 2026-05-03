from pathlib import Path

from stl2fem.datasets import (
    add_nikolaisen_stl_metadata,
    assign_size_bins,
    detect_stl_format,
    nikolaisen_inventory,
)


def test_detect_nikolaisen_stl_format():
    path = Path("data/Nikolaisen2022/Plag Binary meshes/PLAG246-binary.stl")
    if not path.exists():
        return
    assert detect_stl_format(path) in {"ascii", "binary"}


def test_nikolaisen_inventory_bins_if_data_available():
    root = Path("data/Nikolaisen2022")
    if not root.exists():
        return

    inventory = assign_size_bins(nikolaisen_inventory(root), n_bins=4)

    assert len(inventory) == 328
    assert set(inventory["phase"]) == {"OPX", "PLAG"}
    assert inventory["size_bin"].nunique() == 4
    assert inventory["stl_size_bytes"].is_monotonic_increasing


def test_nikolaisen_metadata_filter_for_plag_if_data_available():
    root = Path("data/Nikolaisen2022")
    if not root.exists():
        return

    inventory = add_nikolaisen_stl_metadata(
        nikolaisen_inventory(root, phases=("PLAG",)),
        dataset_root=root,
    )

    assert len(inventory) == 246
    assert (inventory["metadata_evsd_um"] < 1.0).sum() == 245
