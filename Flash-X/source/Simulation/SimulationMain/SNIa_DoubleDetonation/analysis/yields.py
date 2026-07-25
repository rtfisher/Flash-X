#!/usr/bin/env python3
"""
Integrate species masses from a Flash-X (PARAMESH) 2D cylindrical checkpoint.

Sums over leaf blocks only (node type == 1).  In 2D cylindrical (r, z) geometry
the cell volume is dV = 2*pi*r_c*dr*dz, with r_c the cell-centre radius, so the
radial weighting must be applied per cell rather than per block.

Usage:  yields.py <checkpoint> [<checkpoint> ...]
"""

import sys

import h5py
import numpy as np

MSUN = 1.98892e33

SPECIES = [
    "he4 ", "c12 ", "o16 ", "ne20", "mg24", "si28",
    "s32 ", "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
]


def cell_volumes(f):
    """Return per-cell volumes for leaf blocks, shape (nleaf, nz, ny, nx)."""
    ntype = f["node type"][:]
    leaf = ntype == 1
    bbox = f["bounding box"][:][leaf]          # (nleaf, 3, 2)
    shape = f["dens"].shape                     # (nblocks, nz, ny, nx)
    nz, ny, nx = shape[1], shape[2], shape[3]

    rlo, rhi = bbox[:, 0, 0], bbox[:, 0, 1]
    zlo, zhi = bbox[:, 1, 0], bbox[:, 1, 1]

    dr = (rhi - rlo) / nx
    dz = (zhi - zlo) / ny

    # cell-centre radii, shape (nleaf, nx)
    i = np.arange(nx) + 0.5
    r_c = rlo[:, None] + dr[:, None] * i[None, :]

    # dV = 2*pi*r*dr*dz  -> broadcast to (nleaf, nz, ny, nx)
    dV = 2.0 * np.pi * r_c * dr[:, None]        # (nleaf, nx)
    dV = dV[:, None, None, :] * dz[:, None, None, None]
    dV = np.broadcast_to(dV, (bbox.shape[0], nz, ny, nx))
    return leaf, dV


def analyse(path):
    with h5py.File(path, "r") as f:
        leaf, dV = cell_volumes(f)
        dens = f["dens"][:][leaf]
        cell_mass = dens * dV

        total = cell_mass.sum()
        out = {"file": path, "total_Msun": total / MSUN, "nleaf": int(leaf.sum())}

        for s in SPECIES:
            if s in f:
                out[s.strip()] = (f[s][:][leaf] * cell_mass).sum() / MSUN

        # simulation time, if recorded
        try:
            rr = f["real scalars"][:]
            for name, val in rr:
                if name.decode().strip() == "time":
                    out["time"] = val
        except Exception:
            pass
    return out


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    rows = [analyse(p) for p in sys.argv[1:]]

    hdr = f"{'FILE':<34} {'t(s)':>7} {'M_tot':>9} {'ni56':>10} {'si28':>10} {'c12':>9} {'o16':>9}"
    print(hdr)
    for r in rows:
        print(
            f"{r['file'].split('/')[-1]:<34} "
            f"{r.get('time', float('nan')):>7.4f} "
            f"{r['total_Msun']:>9.5f} "
            f"{r.get('ni56', 0.0):>10.3e} "
            f"{r.get('si28', 0.0):>10.3e} "
            f"{r.get('c12', 0.0):>9.5f} "
            f"{r.get('o16', 0.0):>9.5f}"
        )


if __name__ == "__main__":
    main()
