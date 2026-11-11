#!/usr/bin/env python3
"""
sod_exact.py -- plot the Sod shock tube against its EXACT answer.

WHAT THIS DOES
    The Sod problem is one of the few in fluid dynamics whose answer is known
    exactly, on paper, forever.  So we can do something you cannot usually do:
    put the simulation and the truth on the same axes and look.

    Three features appear, and the plot labels them for you:

        rarefaction -- a smooth fan where the gas expands and cools
        contact     -- a jump in density only; pressure and velocity glide through
        shock       -- a jump in everything, the fastest-moving feature

    Where the points sit on the line, the code is right.  Where they smear
    across a jump, that is the grid: a shock is infinitely thin in theory and
    a couple of cells wide on any grid, which is not an error but a limit of
    resolution.  Rerun at higher resolution and the smear narrows.

HOW TO RUN
    python3 sod_exact.py                      # finds the last plotfile here
    python3 sod_exact.py sod_hdf5_chk_0001
    python3 sod_exact.py sod_hdf5_chk_0001 pres

    field : dens (default), pres, velx  -- or 'all' for a three-panel figure

OUTPUT
    sod_exact_<field>.png
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

# Fallbacks only.  The real values are read out of your own output file, so
# if you change the initial condition in flash.par the comparison follows.
RHO_L, P_L, U_L = 1.0, 1.0, 0.0
RHO_R, P_R, U_R = 0.125, 0.1, 0.0
GAMMA = 1.4
X0 = 0.5                       # position of the initial diaphragm


def read_setup(ds):
    """Read the initial condition back out of the output file.

    Flash-X stores every runtime parameter in its output, so we do not have
    to assume what was in flash.par -- we can just look.
    """
    p = ds.parameters
    g = {}
    for name, key, default in (("rho_l", "sim_rholeft", RHO_L),
                               ("rho_r", "sim_rhoright", RHO_R),
                               ("p_l", "sim_pleft", P_L),
                               ("p_r", "sim_pright", P_R),
                               ("u_l", "sim_uleft", U_L),
                               ("u_r", "sim_uright", U_R),
                               ("gamma", "gamma", GAMMA),
                               ("posn", "sim_posn", X0),
                               ("xangle", "sim_xangle", 0.0)):
        g[name] = float(p.get(key, default))
    return g


def exact_riemann(x, t, s):
    """The analytic Sod solution at time t, evaluated at positions x.

    `s` is the setup dictionary from read_setup(), so this works for whatever
    initial condition your run actually used.

    This is the textbook Riemann solution.  The only thing that needs solving
    numerically is the pressure in the star region between the rarefaction and
    the shock; everything else follows from it in closed form.
    """
    g = s["gamma"]
    RHO_L, P_L, U_L = s["rho_l"], s["p_l"], s["u_l"]
    RHO_R, P_R, U_R = s["rho_r"], s["p_r"], s["u_r"]
    X0 = s["posn"]
    cL = np.sqrt(g * P_L / RHO_L)          # sound speed, left
    cR = np.sqrt(g * P_R / RHO_R)          # sound speed, right

    # --- pressure in the star region, by bisection on the velocity mismatch.
    def u_left(p):                          # velocity behind a left-moving wave
        if p > P_L:                         # shock
            A, B = 2.0 / ((g + 1) * RHO_L), (g - 1) / (g + 1) * P_L
            return U_L - (p - P_L) * np.sqrt(A / (p + B))
        # rarefaction
        return U_L - 2 * cL / (g - 1) * ((p / P_L) ** ((g - 1) / (2 * g)) - 1)

    def u_right(p):                         # velocity behind a right-moving wave
        if p > P_R:                         # shock
            A, B = 2.0 / ((g + 1) * RHO_R), (g - 1) / (g + 1) * P_R
            return U_R + (p - P_R) * np.sqrt(A / (p + B))
        return U_R + 2 * cR / (g - 1) * ((p / P_R) ** ((g - 1) / (2 * g)) - 1)

    lo, hi = 1.0e-8, 10.0 * max(P_L, P_R)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if u_left(mid) - u_right(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    p_star = 0.5 * (lo + hi)
    u_star = 0.5 * (u_left(p_star) + u_right(p_star))

    # --- densities either side of the contact
    # Left side: a rarefaction here, so the gas follows an isentrope.
    rho_star_L = RHO_L * (p_star / P_L) ** (1.0 / g)
    # Right side: a shock, so use the Rankine-Hugoniot density jump.
    rr = p_star / P_R
    rho_star_R = RHO_R * (rr + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * rr + 1)

    # --- wave speeds
    c_star_L = np.sqrt(g * p_star / rho_star_L)
    s_head = U_L - cL                       # rarefaction head, moving left
    s_tail = u_star - c_star_L              # rarefaction tail
    s_shock = U_R + cR * np.sqrt((g + 1) / (2 * g) * rr + (g - 1) / (2 * g))

    rho = np.empty_like(x, dtype=float)
    pre = np.empty_like(x, dtype=float)
    vel = np.empty_like(x, dtype=float)

    xi = (x - X0) / t                       # self-similar coordinate
    for i, s in enumerate(xi):
        if s < s_head:                                   # undisturbed left
            rho[i], pre[i], vel[i] = RHO_L, P_L, U_L
        elif s < s_tail:                                 # inside the fan
            vel[i] = 2.0 / (g + 1) * (cL + (g - 1) / 2 * U_L + s)
            c = cL - (g - 1) / 2 * (vel[i] - U_L)
            rho[i] = RHO_L * (c / cL) ** (2.0 / (g - 1))
            pre[i] = P_L * (c / cL) ** (2 * g / (g - 1))
        elif s < u_star:                                 # left of the contact
            rho[i], pre[i], vel[i] = rho_star_L, p_star, u_star
        elif s < s_shock:                                # right of the contact
            rho[i], pre[i], vel[i] = rho_star_R, p_star, u_star
        else:                                            # undisturbed right
            rho[i], pre[i], vel[i] = RHO_R, P_R, U_R

    return dict(dens=rho, pres=pre, velx=vel,
                waves=dict(rarefaction=(X0 + s_head * t, X0 + s_tail * t),
                           contact=X0 + u_star * t,
                           shock=X0 + s_shock * t))


def simulation_profile(path):
    """Pull a 1D cut along the tube out of the (2D) Flash-X output."""
    ds = yt.load(path)
    t = float(ds.current_time)
    s = read_setup(ds)

    # The tube runs along x. A ray straight down the middle gives the profile.
    mid_y = float(ds.domain_center[1])
    ray = ds.ray((float(ds.domain_left_edge[0]), mid_y, 0.0),
                 (float(ds.domain_right_edge[0]), mid_y, 0.0))

    x = np.array(ray["index", "x"])
    order = np.argsort(x)
    prof = {"x": x[order]}
    for f in ("dens", "pres", "velx"):
        prof[f] = np.array(ray["flash", f])[order]
    return t, prof, s


LABELS = {"dens": r"density  $\rho$", "pres": "pressure  $P$",
          "velx": "velocity  $u$"}


def panel(ax, prof, exact, field, annotate):
    ax.plot(prof["x"], prof[field], "o", ms=3.5, color="#E0A030",
            label="simulation", zorder=3)
    ax.plot(prof["x"], exact[field], "-", lw=1.8, color="#3C5A6E",
            label="exact answer")
    ax.set_ylabel(LABELS[field])
    if annotate:
        w = exact["waves"]
        for xpos, name in ((0.5 * sum(w["rarefaction"]), "rarefaction"),
                           (w["contact"], "contact"),
                           (w["shock"], "shock")):
            ax.axvline(xpos, color="#B8B2A7", lw=0.8, ls=":", zorder=1)
            ax.annotate(name, xy=(xpos, ax.get_ylim()[1]), xytext=(0, -12),
                        textcoords="offset points", ha="center", fontsize=8,
                        color="#666666", rotation=90, va="top")


def main():
    args = [a for a in sys.argv[1:]]
    field = "dens"
    if args and args[-1] in ("dens", "pres", "velx", "all"):
        field = args.pop()

    path = args[0] if args else None
    if path is None:
        cands = sorted(glob.glob("sod_hdf5_chk_*") + glob.glob("sod_hdf5_plt_cnt_*"))
        if not cands:
            print("No sod output files found here. Run the simulation first.")
            sys.exit(1)
        path = cands[-1]

    t, prof, s = simulation_profile(path)
    print(f"file  : {path}")
    print(f"time  : {t:.5f} s")
    print(f"setup : rho {s['rho_l']:g}/{s['rho_r']:g}  "
          f"P {s['p_l']:g}/{s['p_r']:g}  gamma {s['gamma']:g}")

    if t <= 0.0:
        print("\nThis file is at t = 0, so nothing has happened yet.")
        print("Use a later output -- the exact solution is undefined at t = 0.")
        sys.exit(1)

    # Flash-X's own Sod test tilts the diaphragm (sim_xangle = 45) to exercise
    # the multidimensional solver.  That is a fine code test but a poor first
    # picture: a horizontal cut then crosses the wave at an angle, so the
    # profile is stretched and the velocity is foreshortened, and the exact
    # 1D answer will NOT lie on it.  Use the workshop parfile instead.
    if abs(s["xangle"]) > 1.0e-6:
        print(f"\n*** This run has sim_xangle = {s['xangle']:g} degrees, so the shock")
        print("    tube is tilted and a straight cut along x does not follow it.")
        print("    Rerun with the workshop parfile (sim_xangle = 0, sim_posn = 0.5),")
        print("    which is the standard textbook Sod problem.")
        sys.exit(1)

    exact = exact_riemann(prof["x"], t, s)
    w = exact["waves"]
    print(f"shock at x = {w['shock']:.4f}   contact at x = {w['contact']:.4f}")

    fields = ["dens", "pres", "velx"] if field == "all" else [field]
    fig, axes = plt.subplots(len(fields), 1, sharex=True,
                             figsize=(7.0, 3.0 * len(fields) + 0.6))
    axes = np.atleast_1d(axes)
    for ax, f in zip(axes, fields):
        panel(ax, prof, exact, f, annotate=(f == fields[0]))
    axes[0].legend(frameon=False, loc="upper right")
    axes[-1].set_xlabel("x")
    axes[0].set_title(f"Sod shock tube at t = {t:.3f}")
    fig.tight_layout()

    out = f"sod_exact_{field}.png"
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
