from pathlib import Path

import numpy as np

from stl2fem.patran import read_patran_neutral


def _node(nid, x, y, z):
    return (f"01{nid:8d}{0:8d}{2:8d}{0:8d}{0:8d}{0:8d}{0:8d}{0:8d}\n"
            f"{x:16.8e}{y:16.8e}{z:16.8e}\n"
            f"0G       6       0       0  000000\n")


def _tet(eid, nodes, pid=1):
    return (f"02{eid:8d}{5:8d}{2:8d}{0:8d}{0:8d}{0:8d}{0:8d}{0:8d}\n"
            f"{4:8d}{0:8d}{pid:8d}{0:8d}{0.0:16.8e}{0.0:16.8e}{0.0:16.8e}\n"
            + "".join(f"{n:8d}" for n in nodes) + "\n")


def test_read_patran_neutral_tetrahedra(tmp_path: Path):
    text = (f"25{0:8d}{0:8d}{1:8d}{0:8d}{0:8d}{0:8d}{0:8d}{0:8d}\nunit test\n"
            + _node(10, 0, 0, 0) + _node(11, 1, 0, 0) + _node(12, 0, 1, 0) + _node(13, 0, 0, 1)
            + _node(14, 1, 1, 1)
            + _tet(1, (10, 11, 12, 13)) + _tet(2, (11, 12, 13, 14), pid=2)
            + f"99{0:8d}{0:8d}{1:8d}{0:8d}{0:8d}{0:8d}{0:8d}{0:8d}\n")
    path = tmp_path / "two_tets.pat"
    path.write_text(text)
    mesh = read_patran_neutral(path)
    assert mesh.title == "unit test"
    assert mesh.node_ids.tolist() == [10, 11, 12, 13, 14]
    np.testing.assert_allclose(mesh.coords[4], [1.0, 1.0, 1.0])
    assert mesh.tets.tolist() == [[0, 1, 2, 3], [1, 2, 3, 4]]
    assert mesh.property_ids.tolist() == [1, 2]
    assert mesh.skipped_elements == 0
