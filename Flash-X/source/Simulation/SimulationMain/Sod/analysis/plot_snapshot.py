#!/usr/bin/env python3
"""
plot_snapshot.py -- make a publication-style slice plot of one output file.

WHAT THIS DOES
    Renders a two-dimensional slice of a chosen variable.  In a 2D run the
    "slice" is really the whole computational domain -- yt just needs to be
    told which coordinate to look along.  Which one that is depends on the
    geometry of the run, so this script works it out for you:

        2D cylindrical (r, z)   -- as in the supernova -- look along 'theta'
        2D cartesian  (x, y)    -- as in Sod and Sedov -- look along 'z'

    That means the same command works on all three of the workshop's problems.

HOW TO RUN
    python3 plot_snapshot.py <file> [field] [zoom]

    python3 plot_snapshot.py sod_hdf5_plt_cnt_0001 dens
    python3 plot_snapshot.py sedov_hdf5_plt_cnt_0005 dens
    python3 plot_snapshot.py snd_hdf5_plt_cnt_0005 enuc
    python3 plot_snapshot.py snd_hdf5_plt_cnt_0005 temp 4

    field  : dens, temp, pres, enuc, velx, vely   (plotfiles)
             ...plus all 13 species               (checkpoints only)
    zoom   : optional zoom factor, default 1 (whole domain)

NOTE ON FIELD NAMES
    Flash-X pads short species names to four characters, so the fields are
    literally 'c12 ', 'o16 ', 'he4 ', 's32 ' -- WITH a trailing space.  Four
    character names such as 'ni56' and 'si28' have none.  This trips everyone
    up at least once.

OUTPUT
    snapshot_<field>_<basename>.png
"""

import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import yt

warnings.filterwarnings("ignore")
yt.set_log_level(50)

# Sensible defaults so beginners get a readable picture first time.
PRESETS = {
    "dens": dict(log=True,  cmap="viridis", label=r"$\rho$ [g cm$^{-3}$]"),
    "temp": dict(log=True,  cmap="hot",     label="T [K]"),
    "pres": dict(log=True,  cmap="plasma",  label=r"P [erg cm$^{-3}$]"),
    "enuc": dict(log=True,  cmap="inferno", label=r"$\dot{e}_{\rm nuc}$ [erg g$^{-1}$ s$^{-1}$]"),
    "ni56": dict(log=False, cmap="magma",   label=r"X($^{56}$Ni)"),
    "si28": dict(log=False, cmap="magma",   label=r"X($^{28}$Si)"),
    "c12 ": dict(log=False, cmap="Blues",   label=r"X($^{12}$C)"),
    "o16 ": dict(log=False, cmap="Greens",  label=r"X($^{16}$O)"),
}


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    path = sys.argv[1]
    field = sys.argv[2] if len(sys.argv) > 2 else "temp"
    zoom = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0

    # Flash-X pads species names; be forgiving about the trailing space.
    ds = yt.load(path)
    available = [f[1] for f in ds.field_list]
    if field not in available:
        if field + " " in available:
            field = field + " "
        else:
            print(f"Field '{field}' not in this file.")
            print("Available:", sorted(available))
            print("\n(Species live in checkpoints -- snd_hdf5_chk_* -- not plotfiles.)")
            sys.exit(1)

    opts = PRESETS.get(field, dict(log=True, cmap="viridis", label=field))

    print(f"file  : {path}")
    print(f"time  : {float(ds.current_time):.5f} s")
    print(f"field : '{field}'")

    # Which coordinate do we look ALONG?  In 2D the third axis is the one the
    # simulation does not resolve, and its name depends on the geometry:
    # cylindrical (r, z) calls it 'theta'; cartesian (x, y) calls it 'z'.
    #
    # We do not guess.  The geometry is a runtime parameter you set in
    # flash.par, Flash-X writes it into every output file, and yt reads it
    # back as ds.geometry -- so the file tells us what it is.  Sod and Sedov
    # are cartesian, the supernova is cylindrical, and this serves all three.
    normal = {"cylindrical": "theta", "polar": "z", "spherical": "phi"}.get(
        str(ds.geometry), "z")
    print(f"geom  : {ds.geometry} -- slicing along '{normal}'")

    # origin='native' keeps the axes in physical coordinates rather than
    # recentring them on the middle of the domain.
    p = yt.SlicePlot(ds, normal, ("flash", field), origin="native")

    # Choosing colour limits needs a little care.  'enuc' is zero wherever
    # nothing is burning and can even go NEGATIVE where neutrino losses exceed
    # nuclear heating, so a naive log scale raises "Invalid vmin or vmax".
    # Set explicit limits from the positive part of the data instead.
    use_log = opts["log"]
    if use_log:
        import numpy as np
        vals = np.array(ds.all_data()["flash", field])
        pos = vals[vals > 0.0]
        if pos.size == 0:
            use_log = False                      # nothing positive; go linear
            print("  (no positive values -- switching to a linear scale)")
        else:
            vmax = float(pos.max())
            vmin = max(vmax * 1.0e-8, float(pos.min()))
            p.set_zlim(("flash", field), vmin, vmax)
            print(f"  colour range: [{vmin:.3e}, {vmax:.3e}]")

    p.set_log(("flash", field), use_log)
    p.set_cmap(("flash", field), opts["cmap"])
    p.set_colorbar_label(("flash", field), opts["label"])
    # Pick a sane length unit.  The supernova domain is thousands of km; the
    # Sod tube is about 1 cm.  Forcing km on Sod would label every axis 1e-5.
    width_cm = float(ds.domain_width.in_units("cm").max())
    p.set_axes_unit("km" if width_cm > 1.0e6 else "cm")
    p.annotate_timestamp(corner="upper_left", draw_inset_box=True)
    if zoom != 1.0:
        p.zoom(zoom)

    out = f"snapshot_{field.strip()}_{path.split('/')[-1]}.png"
    p.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
