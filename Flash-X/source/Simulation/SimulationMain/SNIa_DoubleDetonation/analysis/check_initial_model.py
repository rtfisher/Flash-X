#!/usr/bin/env python3
"""
check_initial_model.py -- verify the simulation's initial condition.

WHAT THIS DOES
    Loads the t = 0 checkpoint, extracts density, temperature and composition
    along the symmetry axis, and overplots the one-dimensional profile file
    that the simulation was built from.  If the two disagree, everything
    downstream is meaningless.

WHY IT MATTERS
    A simulation inherits every error in its initial condition.  Several
    profiles shipped with SNIa_DoubleDetonation are missing their oxygen
    column, and a star initialised from one of those is nearly pure carbon --
    with nothing in the log to warn you.  This script is how you catch that.

HOW TO RUN
    python3 check_initial_model.py snd_hdf5_chk_0000 Profiles/<name>.dat

    e.g.
    python3 check_initial_model.py \
        snd_hdf5_chk_0000 Profiles/wd_profile_rhoc_415e5_M_106_correct.dat

OUTPUT
    initial_model_check.png
"""

import sys
import warnings

import matplotlib
matplotlib.use("Agg")            # render to file, no interactive window needed
import matplotlib.pyplot as plt
import numpy as np
import yt

warnings.filterwarnings("ignore")
yt.set_log_level(50)             # quiet yt down; raise to 20 for chatter

KM = 1.0e5                       # centimetres in a kilometre


def read_profile(path):
    """Read a Flash-X 1D WD profile file.

    Format:  line 1 = '# radius density temperature <species...>'
             line 2 = number of points
             then that many rows of numbers, then optional trailing comments.
    """
    lines = open(path).read().splitlines()
    names = lines[0].lstrip("#").split()
    npts = int(lines[1].split()[0])
    data = np.array([[float(v) for v in row.split()]
                     for row in lines[2:2 + npts]])
    return names, data


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    chk_path, prof_path = sys.argv[1], sys.argv[2]

    # ---- the 1D input profile -------------------------------------------
    names, prof = read_profile(prof_path)
    print(f"profile columns : {names}")
    r_prof = prof[:, 0] / KM
    rho_prof = prof[:, 1]
    T_prof = prof[:, 2]

    # ---- the simulation's t=0 state -------------------------------------
    ds = yt.load(chk_path)
    print(f"checkpoint time : {float(ds.current_time):.4e} s")

    # A "ray" is a 1D line through the dataset.  We shoot one along +z at
    # r = 0, i.e. straight up the symmetry axis.  In Flash-X 2D cylindrical
    # the axes are (r, z, theta).
    ray = ds.ray(ds.domain_left_edge, [0.0, float(ds.domain_right_edge[1]), 0.0])
    order = np.argsort(ray["index", "z"])
    z_sim = np.array(ray["index", "z"][order]) / KM
    rho_sim = np.array(ray["flash", "dens"][order])
    T_sim = np.array(ray["flash", "temp"][order])

    # ---- plot -------------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

    ax[0].semilogy(r_prof, rho_prof, "k-", lw=2, label="1D profile (input)")
    ax[0].semilogy(z_sim, rho_sim, "r--", lw=1.5, label="simulation, t=0")
    ax[0].set_xlabel("radius [km]"); ax[0].set_ylabel(r"$\rho$ [g cm$^{-3}$]")
    ax[0].set_title("Density"); ax[0].legend(); ax[0].set_xlim(0, r_prof.max() * 1.2)

    ax[1].semilogy(r_prof, T_prof, "k-", lw=2, label="1D profile (input)")
    ax[1].semilogy(z_sim, T_sim, "r--", lw=1.5, label="simulation, t=0")
    ax[1].set_xlabel("radius [km]"); ax[1].set_ylabel("T [K]")
    ax[1].set_title("Temperature  (note the hot matchhead at r=0!)")
    ax[1].legend(); ax[1].set_xlim(0, r_prof.max() * 1.2)

    # composition: plot every species column the profile declares
    for j, nm in enumerate(names[3:]):
        ax[2].plot(r_prof, prof[:, 3 + j], lw=2, label=f"{nm} (input)")
    # and the same species as the simulation actually initialised them
    for nm in ("c12 ", "o16 ", "he4 "):
        if ("flash", nm) in ds.field_list:
            ax[2].plot(z_sim, np.array(ray["flash", nm][order]), "--", lw=1.2,
                       label=f"{nm.strip()} (sim)")
    ax[2].set_xlabel("radius [km]"); ax[2].set_ylabel("mass fraction")
    ax[2].set_title("Composition"); ax[2].legend(fontsize=8)
    ax[2].set_xlim(0, r_prof.max() * 1.2); ax[2].set_ylim(-0.05, 1.05)

    fig.tight_layout()
    fig.savefig("initial_model_check.png", dpi=130)
    print("wrote initial_model_check.png")

    # ---- an automatic sanity check ---------------------------------------
    print("\n--- sanity check ---")
    if ("flash", "o16 ") in ds.field_list:
        xo16 = float(np.array(ray["flash", "o16 "][order])[len(z_sim) // 3])
        print(f"X(O16) one third of the way out : {xo16:.4f}")
        if xo16 < 1.0e-3:
            print("  *** WARNING: oxygen is essentially absent.")
            print("  *** Your profile file probably has no o16 column, so the")
            print("  *** star has been renormalised to nearly pure carbon.")
        else:
            print("  looks like a genuine C/O white dwarf.")


if __name__ == "__main__":
    main()
