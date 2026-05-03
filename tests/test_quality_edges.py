import numpy as np
import pyvista as pv

from stl2fem.quality import inspect_volume_mesh


def test_inspect_volume_mesh_reports_realized_edge_lengths(tmp_path):
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    grid = pv.UnstructuredGrid({pv.CellType.TETRA: np.array([[0, 1, 2, 3]])}, points)
    path = tmp_path / "one_tet.vtu"
    grid.save(path)

    stats = inspect_volume_mesh(path)

    assert stats["n_nodes"] == 4
    assert stats["n_tets"] == 1
    assert np.isclose(stats["edge_length_min"], 1.0)
    assert np.isclose(stats["edge_length_max"], np.sqrt(2.0))
