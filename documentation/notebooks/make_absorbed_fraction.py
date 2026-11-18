"""Absorbed fraction of the photocathode vs wavelength, as OMSim's thin-film photocathode actually absorbs it, for the
PMT response's QE weight (weight = measured QE / absorbed fraction, jAbsorbedFractionFileName in the PMT .dat).

The file has to match the simulated film: if it says the photocathode absorbs more than the film does, every hit is
weighted too low and the simulated QE falls short of the measured one. The 4-inch file (LOM_Hamamatsu_4inch_intrinsic_QE)
gives 0.733 at 400 nm, while OMSim's film absorbs 0.613, so 4-inch hits were weighted ~16 % too low.

Input: one run per wavelength of the bare PMT in air, normal incidence, 10 mm beam on the window centre, 100k photons,
with OMSIM_HIT_DUMP (one line per photon absorbed in the photocathode):

    OMSIM_HIT_DUMP=<dir>/m<model>_l<nm>.txt ./OMSim_effective_area --detector_type 1 --pmt_model <model> \
        --environment 0 -n 100000 -t 0 -f 0 -r 10 -l <nm> --no_header -o <dir>/o_m<model>_l<nm>

    python make_absorbed_fraction.py <dir>
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = 100000
WAVELENGTHS = np.arange(250, 796, 5)
QE_DIR = "../../common/data/PMTs/measurement_matching_data/QE"
PMTS = {2: ("4-inch", "LOM_Hamamatsu_4inch", "POM_Hamamatsu_4inch_intrinsic_QE_OMSim.dat"),
        0: ("3-inch R15458", "mDOM_Hamamatsu_R15458", None)}
OLD = {2: "LOM_Hamamatsu_4inch_intrinsic_QE.dat", 0: "mDOM_Hamamatsu_R15458_CT_intrinsic_QE.dat"}
MEAN_QE = {2: "LOM_Hamamatsu_4inch_mean_QE.dat", 0: "mDOM_Hamamatsu_R15458_mean_QE.dat"}

scan_dir = sys.argv[1]
fig, (ax, axq) = plt.subplots(1, 2, figsize=(14, 5))
for model, (label, _, out_name) in PMTS.items():
    counts = np.array([sum(1 for _ in open(f)) if os.path.exists(f := f"{scan_dir}/m{model}_l{l}.txt") else 0
                       for l in WAVELENGTHS])
    p = counts / N
    err = np.sqrt(p * (1 - p) / N)
    old = np.loadtxt(f"{QE_DIR}/{OLD[model]}", comments="#")
    qe = np.loadtxt(f"{QE_DIR}/{MEAN_QE[model]}", comments="#")
    i400 = np.searchsorted(WAVELENGTHS, 400)
    print(f"{label}: absorbed fraction at 400 nm, OMSim film {p[i400]:.3f}, old file {np.interp(400, old[:, 0], old[:, 1]):.3f}")
    if out_name:
        path = f"{QE_DIR}/{out_name}"
        np.savetxt(path, np.c_[WAVELENGTHS, p, err], delimiter="\t",
                   header="Wavelength(nm) \t QE \t error\nAbsorbed fraction of OMSim's thin-film photocathode "
                          "(Surf_Generic_Photocathode_20nm), bare PMT in air, normal incidence, window centre, "
                          f"{N} photons per wavelength (make_absorbed_fraction.py)")
        print(f"  wrote {os.path.normpath(path)}")
    line, = ax.plot(WAVELENGTHS, p, "-", lw=2, label=f"{label}: OMSim film")
    ax.plot(old[:, 0], old[:, 1], "--", color=line.get_color(), lw=1.5, label=f"{label}: old file")
    # simulated QE at the window centre = absorbed x (measured QE / absorbed-fraction file)
    q = np.interp(WAVELENGTHS, qe[:, 0], qe[:, 1])
    with np.errstate(divide="ignore", invalid="ignore"):
        sim_old = np.where(np.interp(WAVELENGTHS, old[:, 0], old[:, 1]) > 0, p * q / np.interp(WAVELENGTHS, old[:, 0], old[:, 1]), 0)
    axq.plot(WAVELENGTHS, q, "-", color=line.get_color(), lw=2, label=f"{label}: measured QE (= simulated with the OMSim file)")
    axq.plot(WAVELENGTHS, sim_old, ":", color=line.get_color(), lw=2, label=f"{label}: simulated QE with the old file")
ax.set_xlabel("wavelength [nm]")
ax.set_ylabel("absorbed fraction")
ax.set_title("Photocathode absorbed fraction (air, normal incidence)", loc="left")
ax.legend(frameon=False, fontsize=9)
ax.grid(alpha=0.3)
axq.set_xlabel("wavelength [nm]")
axq.set_ylabel("QE at the window centre")
axq.set_title("Simulated vs measured QE", loc="left")
axq.legend(frameon=False, fontsize=9)
axq.grid(alpha=0.3)
fig.patch.set_facecolor("white")
fig.tight_layout()
fig.savefig("absorbed_fraction.png", dpi=110)
print("saved absorbed_fraction.png")
