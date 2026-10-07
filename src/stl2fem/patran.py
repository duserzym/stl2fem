"""Reader for Patran neutral (``.pat``) tetrahedral meshes, the MERRILL legacy mesh format.

Patran neutral files are packets of fixed-width cards. Each header card is ``IT ID IV KC N1..N5`` (I2, 8 x I8),
where ``IT`` is the packet type and ``KC`` the number of data cards that follow. This reader uses:

* packet 01 (node): card 1 holds X, Y, Z (3E16.9);
* packet 02 (element): ``IV`` is the shape code (5 = tetrahedron); card 1 is ``NODES CONFIG PID CEID ...``, card 2
  the node IDs (I8 fields).

Other packets (title, summary, properties, materials) are skipped by their card counts.
"""

from __future__ import annotations

from array import array
from dataclasses import dataclass
from pathlib import Path

import numpy as np

TET_SHAPE = 5


@dataclass
class PatranMesh:
    node_ids: np.ndarray  # (n,) int64, as written in the file
    coords: np.ndarray  # (n, 3) float64, file units
    element_ids: np.ndarray  # (m,) int64
    tets: np.ndarray  # (m, 4) int64, zero-based indices into coords
    property_ids: np.ndarray  # (m,) int64 (PID; MERRILL uses it as the material/block index)
    title: str
    skipped_elements: int  # non-tetrahedral elements ignored


def _ints(card: str) -> list[int]:
    return [int(card[k:k + 8]) for k in range(0, len(card.rstrip("\r\n")), 8) if card[k:k + 8].strip()]


def read_patran_neutral(path: str | Path) -> PatranMesh:
    path = Path(path)
    node_ids, xyz = array("q"), array("d")
    elem_ids, conn, pids = array("q"), array("q"), array("q")
    title, skipped = "", 0
    with path.open("r", encoding="ascii", errors="replace") as fh:
        for header in fh:
            if len(header.strip()) == 0:
                continue
            it = int(header[0:2])
            ident = int(header[2:10])
            iv = int(header[10:18])
            kc = int(header[18:26])
            if it == 99:
                break
            cards = [next(fh) for _ in range(kc)]
            if it == 1:
                c = cards[0]
                node_ids.append(ident)
                xyz.extend((float(c[0:16]), float(c[16:32]), float(c[32:48])))
            elif it == 2:
                head = _ints(cards[0][:32])
                nnodes, pid = head[0], head[2]
                if iv != TET_SHAPE or nnodes != 4:
                    skipped += 1
                    continue
                nodes = []
                for c in cards[1:]:
                    nodes.extend(_ints(c))
                elem_ids.append(ident)
                conn.extend(nodes[:4])
                pids.append(pid)
            elif it == 25 and cards:
                title = cards[0].strip()
    nid = np.frombuffer(node_ids, dtype=np.int64).copy()
    coords = np.frombuffer(xyz, dtype=np.float64).reshape(-1, 3).copy()
    raw = np.frombuffer(conn, dtype=np.int64).reshape(-1, 4)
    lookup = np.full(int(nid.max()) + 1, -1, dtype=np.int64)
    lookup[nid] = np.arange(len(nid))
    tets = lookup[raw]
    if (tets < 0).any():
        raise ValueError(f"{path}: elements reference undefined nodes")
    return PatranMesh(nid, coords, np.frombuffer(elem_ids, dtype=np.int64).copy(), tets,
                      np.frombuffer(pids, dtype=np.int64).copy(), title, skipped)
