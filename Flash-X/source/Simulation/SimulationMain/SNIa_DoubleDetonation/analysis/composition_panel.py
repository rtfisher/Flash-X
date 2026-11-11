#!/usr/bin/env python3
"""
composition_panel.py -- four-panel map of the nucleosynthesis.

WHAT THIS DOES
    Renders X(C12), X(O16), X(Si28) and X(Ni56) side by side from a checkpoint,
    revealing the layered "onion" structure the explosion produces.

WHY CHECKPOINTS AND NOT PLOTFILES
    Plotfiles contain only the variables named by plot_var_1..N in the
    parameter file -- typically dens, temp, pres, enuc and the velocities.
    Checkpoints contain EVERY variable, including all 13 species.  So all
    composition analysis must use snd_hdf5_chk_* files.

WHAT TO LOOK FOR
    Ni56 sits in the hottest, densest, most completely burned material at the
    centre.  Si28 forms a shell further out where burning stalled short of
    nuclear statistical equilibrium.  Unburned C and O remain outside that.
    This stratification is the origin of the characteristic spectral evolution
    of a Type Ia supernova: you see silicon early, from the outer layers, and
    iron-peak elements later as the ejecta thin and the interior is revealed.

HOW TO RUN
    python3 composition_panel.py snd_hdf5_chk_0001

OUTPUT
    composition_panel.png
"""

import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import yt

warnings.filterwarnings("ignore")
yt.set_log_level(50)

# (field name as stored, pretty label, colour map)
PANELS = [
    ("c12 ", r"X($^{12}$C)  unburned fuel",       "Blues"),
    ("o16 ", r"X($^{16}$O)  unburned fuel",       "Greens"),
    ("si28", r"X($^{28}$Si) intermediate mass",   "Oranges"),
    ("ni56", r"X($^{56}$Ni) iron peak / NSE",     "magma"),
]


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)

    path = sys.argv[1]
    if "chk" not in path:
        print("WARNING: that does not look like a checkpoint.")
        print("Species are only stored in snd_hdf5_chk_* files.\n")

    ds = yt.load(path)
    print(f"file : {path}")
    print(f"time : {float(ds.current_time):.5f} s")

    have = [f[1] for f in ds.field_list]
    fields = [(f, lab, cm) for (f, lab, cm) in PANELS if f in have]
    if not fields:
        print("No species fields found -- is this a plotfile rather than a checkpoint?")
        print("Available:", sorted(have))
        sys.exit(1)

    # yt can build a multi-panel figure directly from a list of fields.
    p = yt.SlicePlot(ds, "theta", [("flash", f) for f, _, _ in fields],
                     origin="native")
    for f, lab, cm in fields:
        p.set_log(("flash", f), False)          # mass fractions: linear
        p.set_zlim(("flash", f), 0.0, 1.0)      # common scale, 0 to 1
        p.set_cmap(("flash", f), cm)
        p.set_colorbar_label(("flash", f), lab)
    p.set_axes_unit("km")
    p.annotate_timestamp(corner="upper_left", draw_inset_box=True)

    out = "composition_panel.png"
    p.save(out)
    print(f"wrote {out}")

    # A quick numerical summary to accompany the picture.
    ad = ds.all_data()
    print("\npeak mass fractions in the domain:")
    for f, lab, _ in fields:
        print(f"  {f.strip():>5} : {float(ad.max(('flash', f))):.4f}")


if __name__ == "__main__":
    main()
