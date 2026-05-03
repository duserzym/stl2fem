import numpy as np
import pyvista as pv

from stl2fem.conversion import write_meter_scaled_mesh
from stl2fem.conversion import tetrahedralize_bruteforce_stl_for_merrill
from stl2fem.quality import estimate_grain_size_stl, inspect_volume_mesh
from stl2fem.units import make_unit_context, scale_for_unit


def test_unit_context_converts_9nm_to_micrometer_native_units():
    units = make_unit_context(input_unit="um", target_edge_length_m=9e-9)

    assert units.input_scale_to_meters == 1e-6
    assert np.isclose(units.target_edge_length_native, 0.009)


def test_unit_aliases():
    assert scale_for_unit("µm") == 1e-6
    assert scale_for_unit("micron") == 1e-6
    assert scale_for_unit("nm") == 1e-9


def test_write_meter_scaled_mesh(tmp_path):
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    grid = pv.UnstructuredGrid({pv.CellType.TETRA: np.array([[0, 1, 2, 3]])}, points)
    native = tmp_path / "native.vtu"
    scaled = tmp_path / "scaled.msh"
    grid.save(native)

    write_meter_scaled_mesh(native, scaled, input_scale_to_meters=1e-6)
    stats = inspect_volume_mesh(scaled)

    assert np.isclose(stats["edge_length_min"], 1e-6)
    assert np.isclose(stats["edge_length_max"], np.sqrt(2.0) * 1e-6)


def test_estimate_grain_size_for_cube_stl(tmp_path):
    cube = pv.Cube(bounds=(0, 0.1, 0, 0.1, 0, 0.1)).triangulate()
    path = tmp_path / "cube.stl"
    cube.save(path)

    stats = estimate_grain_size_stl(path, input_unit="um")

    assert np.isclose(stats["grain_equivalent_cube_edge_nm"], 100.0)
    assert np.isclose(stats["grain_bbox_max_nm"], 100.0)
    assert np.isclose(stats["grain_volume_nm3"], 1_000_000.0)


def test_voxel_bruteforce_fallback_writes_meter_msh(tmp_path):
    cube = pv.Cube(bounds=(0, 0.027, 0, 0.027, 0, 0.027)).triangulate()
    stl_path = tmp_path / "rough_cube.stl"
    native_path = tmp_path / "rough_cube_native.msh"
    meter_path = tmp_path / "rough_cube_meter.msh"
    cube.save(stl_path)

    _, _, _, info = tetrahedralize_bruteforce_stl_for_merrill(
        stl_path,
        native_path,
        meter_path,
        input_unit="um",
        target_edge_length_m=9e-9,
        strategy="voxel",
        overwrite=True,
    )
    stats = inspect_volume_mesh(meter_path)

    assert info["mesh_strategy"] == "voxel"
    assert stats["n_tets"] > 0
    assert stats["edge_length_median"] < 2e-8
