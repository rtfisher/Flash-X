#!/usr/bin/env python3
"""
sedov_radius.py -- measure the blast-wave radius over time, and check it
                   against the Sedov-Taylor prediction.

WHAT THIS DOES
    This is the verification lab.  You are not going to trust the code --
    you are going to check it.

    For each output file it finds the shock: the blast wave sweeps gas into
    a thin dense shell, so the densest cell in the box IS the shock front.
    Its distance from the explosion is the shock radius, R.  Do that for
    every file and you have R at a series of times.

    Then it compares with theory -- and the theory is something you can
    derive yourself, with no calculus, in four lines:

      1. All the blast has is the energy you put in, split between motion
         and heat:                       E = kinetic + internal

      2. The blast is SELF-SIMILAR -- a bigger one is just a scaled-up copy
         of a smaller one -- so that split never changes:
                                         internal is proportional to kinetic
         which means we can just track    E ~ M v^2
         (drop the 1/2; a factor of two cannot change a power law)

      3. Put in what M and v actually are.  The shock has gone a distance R
         in a time t, and has swept up everything behind it:
                                         v ~ R/t        M ~ rho R^d
         where d is the number of dimensions it expands into.

      4. Substitute:   E ~ rho R^d (R/t)^2 = rho R^(d+2) / t^2

         and rearrange:      R  =  xi * ( E t^2 / rho ) ^ (1/(d+2))

    with xi a number close to 1 that the full solution pins down.  So the
    dimension is not a detail -- it IS the answer:

        in 2D (this run):   M ~ rho R^2,  so  R ~ (E t^2/rho)^(1/4) ~ t^(1/2)
        in 3D:              M ~ rho R^3,  so  R ~ (E t^2/rho)^(1/5) ~ t^(2/5)

    The script measures the power from your data and tells you what it got.
    That measured exponent, landing on the predicted one, is the verification.

HOW TO RUN
    python3 sedov_radius.py                 # finds the sedov files here
    python3 sedov_radius.py sedov_hdf5_chk_*

OUTPUT
    sedov_radius.png   -- your measured points against the theory curve
    sedov_radius.dat   -- the (t, R) table, if you want to plot it yourself

NOTE
    The t=0 file has no shock yet, so it is skipped automatically.
"""

import glob
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt

warnings.filterwarnings("ignore")
yt.set_log_level(50)


def shock_radius(path):
    """Return (time, radius) for one output file.

    The shock is where the gas piles up, so we simply look for the cell with
    the highest density and measure how far it sits from the blast centre.
    """
    ds = yt.load(path)
    ad = ds.all_data()

    dens = np.array(ad["flash", "dens"])
    # Cell centres, in the same units as the domain.
    x = np.array(ad["index", "x"])
    y = np.array(ad["index", "y"])

    # The explosion went off at the centre of the box.
    cx, cy = (float(c) for c in ds.domain_center[:2])

    r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    return float(ds.current_time), float(r[np.argmax(dens)])


def main():
    files = sys.argv[1:]
    if not files:
        # Checkpoints carry every variable and are written on a clean cadence,
        # so they are the easiest thing to sweep.
        files = sorted(glob.glob("sedov_hdf5_chk_*"))
    if not files:
        print("No sedov output files found here.")
        print("Run the simulation first, or pass filenames on the command line.")
        sys.exit(1)

    print(f"{'time':>12}  {'radius':>12}")
    print(f"{'-'*12}  {'-'*12}")

    times, radii = [], []
    for f in files:
        t, r = shock_radius(f)
        if t <= 0.0:
            continue                      # t=0 has no shock yet
        times.append(t)
        radii.append(r)
        print(f"{t:12.5f}  {r:12.5f}")

    if len(times) < 2:
        print("\nNeed at least two output files with t > 0 to fit a curve.")
        sys.exit(1)

    t = np.array(times)
    R = np.array(radii)
    np.savetxt("sedov_radius.dat", np.column_stack([t, R]),
               header="time  radius", fmt="%.6e")

    # Fit a power law by fitting a straight line to the logarithms:
    # if R = A t^n then log R = log A + n log t, so the SLOPE is the power.
    n, logA = np.polyfit(np.log(t), np.log(R), 1)
    A = np.exp(logA)

    # This run is 2D, so the blast expands into 2 dimensions.
    d = 2
    expected = 2.0 / (d + 2)

    print()
    print()
    print("  Why does the slope give the power?  If  R = A t^n  then taking logs")
    print("  gives  log R = log A + n log t  -- a straight line in log-log, with")
    print("  slope n.  So fitting a LINE to the logs measures the POWER directly.")
    print()
    print(f"measured   R ~ t^{n:.3f}")
    print(f"predicted  R ~ t^{expected:.3f}   (Sedov-Taylor in {d}D)")
    print(f"difference {abs(n - expected) / expected * 100:.1f}%")

    # The similarity constant xi, recovered from your own data.  E and rho
    # are both 1 in this setup, which keeps the arithmetic honest and simple.
    E, rho = 1.0, 1.0
    xi = A / (E / rho) ** (1.0 / (d + 2))
    print(f"similarity constant xi = {xi:.3f}   (theory says close to 1)")

    # Two panels. LEFT: the measurement as you took it, against theory.
    # RIGHT: the same data on log-log, where a power law becomes a STRAIGHT
    # LINE whose slope is the power -- which is how the exponent was actually
    # obtained, and is worth seeing rather than taking on trust.
    fig, (ax, axl) = plt.subplots(1, 2, figsize=(11.0, 4.6))

    ax.plot(t, R, "o", ms=9, color="#E0A030", label="your measurement", zorder=3)
    tt = np.linspace(t.min() * 0.9, t.max() * 1.05, 200)
    ax.plot(tt, xi * (E * tt ** 2 / rho) ** (1.0 / (d + 2)), "-", lw=2,
            color="#3C5A6E", label=rf"Sedov-Taylor:  $R \propto t^{{{expected:.2f}}}$")
    ax.set_xlabel("time")
    ax.set_ylabel("shock radius")
    ax.set_title("as measured")
    ax.legend(frameon=False)

    axl.loglog(t, R, "o", ms=9, color="#E0A030", label="your measurement", zorder=3)
    axl.loglog(tt, A * tt ** n, "-", lw=2, color="#B03A2E",
               label=rf"straight-line fit:  slope $n = {n:.3f}$")
    axl.loglog(tt, xi * (E * tt ** 2 / rho) ** (1.0 / (d + 2)), "--", lw=1.6,
               color="#3C5A6E", label=rf"theory slope {expected:.3f}")
    axl.set_xlabel("log time")
    axl.set_ylabel("log shock radius")
    axl.set_title("log-log: a power law is a straight line")
    axl.legend(frameon=False, fontsize=9)

    fig.suptitle(f"measured  $R \\propto t^{{{n:.3f}}}$   vs  theory  $t^{{{expected:.3f}}}$")
    fig.tight_layout()
    fig.savefig("sedov_radius.png", dpi=140)
    print("\nwrote sedov_radius.png and sedov_radius.dat")


if __name__ == "__main__":
    main()
