#!/usr/bin/env python3
"""
detonation_speed.py -- measure how fast the burning front travels.

WHAT THIS DOES
    For each plotfile, finds the radius at which nuclear energy generation
    peaks -- that is the detonation front -- and fits a straight line to
    radius versus time.  The slope is the front speed.

WHY THIS IS THE BEST EXERCISE IN THE SET
    It is a measurement extracted from a simulation and compared against an
    independent theoretical prediction.  For carbon-oxygen material the
    Chapman-Jouguet detonation speed is of order 1.1e9 cm/s (11,000 km/s).
    If your measurement lands near that, the simulation is doing something
    recognisably right.  If it does not, the interesting question is why.

    Expect disagreement at coarse resolution.  A real detonation reaction zone
    is centimetres thick; our cells are tens of kilometres.  What the code
    actually propagates is a numerically broadened front whose speed is set by
    how fast energy release can push the shock, not by the true reaction-zone
    physics.  Quantifying that discrepancy IS the exercise.

HOW TO RUN
    cd into a run directory, then
        python3 detonation_speed.py

OUTPUT
    detonation_speed.png  and a printed table
"""

import glob
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt

warnings.filterwarnings("ignore")
yt.set_log_level(50)

KM = 1.0e5
CJ_SPEED_CGS = 1.1e9      # approximate C/O Chapman-Jouguet speed, cm/s


def front_radius(ds):
    """Radius of peak nuclear energy generation, in cm.

    Find the cell with the largest enuc and report its spherical radius
    sqrt(r^2 + z^2).  In 2D cylindrical geometry yt names the coordinates
    ('index','r') and ('index','z') -- NOT 'x' and 'y', which is the natural
    guess and gives a YTFieldNotFound.
    """
    ad = ds.all_data()
    enuc = np.array(ad["flash", "enuc"])
    if enuc.max() <= 0.0:
        return np.nan
    k = int(np.argmax(enuc))
    r = float(np.array(ad["index", "r"])[k])
    z = float(np.array(ad["index", "z"])[k])
    return np.sqrt(r * r + z * z)


def main():
    files = sorted(glob.glob("snd_hdf5_plt_cnt_*"))
    if not files:
        print("No plotfiles here.  cd into a run directory first.")
        return

    times, radii = [], []
    print(f"{'t [s]':>10} {'R_front [km]':>14}")
    for f in files:
        ds = yt.load(f)
        t = float(ds.current_time)
        R = front_radius(ds)
        if np.isfinite(R):
            times.append(t); radii.append(R)
            print(f"{t:10.5f} {R / KM:14.1f}")

    times = np.array(times); radii = np.array(radii)
    if len(times) < 3:
        print("\nNot enough frames with active burning to fit a speed.")
        return

    # Fit only the part where the front is genuinely propagating: skip the
    # first couple of frames, which are still the matchhead switching on.
    m = slice(2, None)
    slope, intercept = np.polyfit(times[m], radii[m], 1)

    print("\n--- result ---")
    print(f"measured front speed : {slope:.3e} cm/s  =  {slope/KM:.0f} km/s")
    print(f"C/O Chapman-Jouguet  : {CJ_SPEED_CGS:.3e} cm/s  =  {CJ_SPEED_CGS/KM:.0f} km/s")
    print(f"ratio measured/CJ    : {slope / CJ_SPEED_CGS:.3f}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(times, radii / KM, "o", label="peak $\\dot{e}_{nuc}$ location")
    ax.plot(times[m], (slope * times[m] + intercept) / KM, "r-",
            label=f"fit: {slope/KM:.0f} km/s")
    ax.plot(times, (CJ_SPEED_CGS * times + intercept) / KM, "k--",
            label=f"CJ speed: {CJ_SPEED_CGS/KM:.0f} km/s")
    ax.set_xlabel("time [s]"); ax.set_ylabel("front radius [km]")
    ax.set_title("Detonation front propagation")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig("detonation_speed.png", dpi=130)
    print("\nwrote detonation_speed.png")


if __name__ == "__main__":
    main()
