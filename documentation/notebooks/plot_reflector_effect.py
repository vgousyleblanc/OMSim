"""
Plots the effect of the gelpad side-wall reflector (--pom_reflector) on the full
P-OM module's effective area, in water, as a function of theta.

Reads the raw .dat output of `OMSim_effective_area --detector_type 8` runs (with and
without --pom_reflector) from reflector_water_scan_data/, named
norefl_t<theta>.dat / refl_t<theta>.dat.

To regenerate the underlying data (from the build/ directory):
    for theta in 0 15 30 45 60 75 90 105 120 135 150 165 180; do
        ./OMSim_effective_area --detector_type 8 --environment 3 -n 25000 \
            -t $theta -f 45 -l 400 -o reflector_water_scan_data/norefl_t${theta}
        ./OMSim_effective_area --detector_type 8 --pom_reflector --environment 3 -n 25000 \
            -t $theta -f 45 -l 400 -o reflector_water_scan_data/refl_t${theta}
    done
(environment 3 = water, matching P-OM's real deployment medium)
"""

import os
import numpy as np
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(__file__), "reflector_water_scan_data")
THETAS = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]


def read_ea(path):
    """Returns (EA_Total, EA_Total_error) from the last two tab-separated columns.

    Row layout: phi, theta, wavelength, 16x per-PMT weighted hits, total_hits,
    EA_Total, EA_Total_error - i.e. EA_Total is the *second-to-last* field, not
    third-to-last (that's total_hits, a raw summed hit count, not an area).
    """
    row = open(path).read().strip().split("\n")[-1].split("\t")
    return float(row[-2]), float(row[-1])


def main():
    ea_n, err_n, ea_r, err_r = [], [], [], []
    for theta in THETAS:
        n, ne = read_ea(os.path.join(DATA_DIR, f"norefl_t{theta}.dat"))
        r, re = read_ea(os.path.join(DATA_DIR, f"refl_t{theta}.dat"))
        ea_n.append(n); err_n.append(ne)
        ea_r.append(r); err_r.append(re)

    thetas = np.array(THETAS)
    ea_n, err_n = np.array(ea_n), np.array(err_n)
    ea_r, err_r = np.array(ea_r), np.array(err_r)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8, 8), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]}
    )

    ax1.errorbar(thetas, ea_n, yerr=err_n, fmt="o-", color="#1f77b4",
                 label="No reflector", capsize=3, markersize=6)
    ax1.errorbar(thetas, ea_r, yerr=err_r, fmt="s-", color="#d62728",
                 label="With reflector (all 16 gelpads)", capsize=3, markersize=6)
    ax1.set_ylabel("Effective Area (cm²)")
    ax1.set_title("P-OM full module, water, λ=400nm, φ=45°\n"
                   "Effect of a gelpad side-wall reflector")
    ax1.legend()
    ax1.grid(alpha=0.3)

    ratio = ea_r / ea_n
    ratio_err = ratio * np.sqrt((err_r / ea_r) ** 2 + (err_n / ea_n) ** 2)
    ax2.errorbar(thetas, ratio, yerr=ratio_err, fmt="D-", color="#2ca02c",
                 capsize=3, markersize=5)
    ax2.axhline(1.0, color="k", linestyle="--", linewidth=1, alpha=0.6)
    ax2.set_xlabel("Theta (degrees)")
    ax2.set_ylabel("Ratio (reflector / no reflector)")
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    outpath = os.path.join(os.path.dirname(__file__), "reflector_water_plot.png")
    plt.savefig(outpath, dpi=150)
    print("Saved to", outpath)


if __name__ == "__main__":
    main()
