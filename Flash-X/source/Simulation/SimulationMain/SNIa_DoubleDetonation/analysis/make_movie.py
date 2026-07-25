#!/usr/bin/env python3
"""
make_movie.py -- turn a sequence of plotfiles into an animation.

WHAT THIS DOES
    Renders one frame per plotfile with a FIXED colour scale (so the eye sees
    real change, not autoscaling), then stitches them with ffmpeg if available.

WHY enuc RATHER THAN temp
    Try both.  'temp' is dominated by the hot ash left behind the front, which
    swamps the front itself.  'enuc' is the *rate* of nuclear energy release,
    so only actively burning material lights up and the detonation front
    appears as a thin travelling shell.  Visualising a rate rather than a
    state is often the difference between seeing the physics and not.

HOW TO RUN
    python3 make_movie.py                 # enuc, all plotfiles here
    python3 make_movie.py temp
    python3 make_movie.py enuc 1e15 1e22  # explicit colour limits

OUTPUT
    frames/frame_XXXX.png  and (if ffmpeg is present) movie_<field>.mp4
"""

import glob
import os
import subprocess
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import yt

warnings.filterwarnings("ignore")
yt.set_log_level(50)

CMAP = {"enuc": "inferno", "temp": "hot", "dens": "viridis", "pres": "plasma"}


def main():
    field = sys.argv[1] if len(sys.argv) > 1 else "enuc"
    vmin = float(sys.argv[2]) if len(sys.argv) > 3 else None
    vmax = float(sys.argv[3]) if len(sys.argv) > 3 else None

    files = sorted(glob.glob("snd_hdf5_plt_cnt_*"))
    if not files:
        print("No plotfiles (snd_hdf5_plt_cnt_*) in this directory.")
        print("Run the simulation first, or cd into its run directory.")
        sys.exit(1)
    print(f"found {len(files)} plotfiles")

    # Fixed colour limits matter.  With autoscaling every frame renormalises
    # and a growing detonation can look completely static.  If the user did
    # not supply limits, take them from the LAST frame, which is the most
    # developed.
    if vmin is None:
        ds_last = yt.load(files[-1])
        ad = ds_last.all_data()
        fmax = float(ad.max(("flash", field)))
        vmax = fmax
        vmin = max(fmax * 1e-8, 1e-30)     # eight decades is usually readable
        print(f"auto colour range from final frame: [{vmin:.3e}, {vmax:.3e}]")

    os.makedirs("frames", exist_ok=True)

    skipped = []
    for i, f in enumerate(files):
        # One unreadable plotfile should cost you one frame, not the movie.
        try:
            ds = yt.load(f)
            p = yt.SlicePlot(ds, "theta", ("flash", field), origin="native")
            p.set_log(("flash", field), True)
            p.set_zlim(("flash", field), vmin, vmax)       # <- the fixed scale
            p.set_cmap(("flash", field), CMAP.get(field, "viridis"))
            p.set_axes_unit("km")
            p.annotate_timestamp(corner="upper_left", draw_inset_box=True)
            p.save(f"frames/frame_{i:04d}.png")
        except Exception as exc:
            skipped.append((f, type(exc).__name__))
            print(f"  frame {i:4d}/{len(files)-1}   SKIPPED ({type(exc).__name__})")
            continue
        print(f"  frame {i:4d}/{len(files)-1}   t = {float(ds.current_time):.4f} s")

    if skipped:
        print(f"\n  {len(skipped)} file(s) could not be read and were left out "
              "of the movie:")
        for f, why in skipped:
            print(f"    {f}  ({why})")

    # ---- stitch, if ffmpeg exists ----------------------------------------
    out = f"movie_{field}.mp4"
    try:
        subprocess.run(
            ["ffmpeg", "-y", "-framerate", "10",
             "-pattern_type", "glob", "-i", "frames/frame_*.png",
             "-c:v", "libx264", "-pix_fmt", "yuv420p", "-vf",
             "pad=ceil(iw/2)*2:ceil(ih/2)*2", out],
            check=True, capture_output=True)
        print(f"\nwrote {out}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("\nffmpeg not available (or failed) -- frames/ still written.")
        print("You can animate them elsewhere, or view them in order.")


if __name__ == "__main__":
    main()
