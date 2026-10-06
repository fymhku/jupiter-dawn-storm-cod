# -*- coding: utf-8 -*-
"""
Figures for the Supplemental Material on the quiet-period event only (2017-05-10).

    S_quiet_bandpass.png   band-passed |B_perp| with the active period shaded
    S_quiet_psd_Bperp.png  B_perp PSD alone (no comparison, no band overlay)
    S_quiet_psd_Bpar.png   B_parallel PSD alone

Nothing from the 2017-05-16 dawn storm is included.
"""
import os
import sys

WS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(WS)
sys.path.insert(0, os.path.join(ROOT, "juno_ulf_analysis"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from juno_ulf import bandpass as bp
from juno_ulf import components, cwt, data, psd

MAG_DIR = os.environ.get("JUNO_MAG_DIR", os.path.join(ROOT, "data", "mag"))
OUT = os.path.join(os.environ.get("JUNO_ULF_OUT", os.path.join(ROOT, "output")), "quiet", "supplementary")
os.makedirs(OUT, exist_ok=True)

DATE = "2017_05_10"
EVENT_ID = "2017-05-10"
WIN = ("2017-05-10 08:00:00", "2017-05-10 10:00:00")
ACTIVE = ("2017-05-10 08:44:00", "2017-05-10 08:56:00")
BAND = (2e-3, 5e-3)

CWT_CFG = dict(wavelet="cmor1.5-1.0", n_scales=200,
               freq_range=(1e-3, 50e-3), ion_mass_amu=16.0)
PSD_FMIN, PSD_FMAX = 1e-3, 50e-3

Br, Bth, Bph = data.load_mag_window(DATE, WIN[0], WIN[1], MAG_DIR)
B_total = np.sqrt(Br ** 2 + Bth ** 2 + Bph ** 2)
B_par, B_perp = components.parallel_perp_by_window(Br, Bth, Bph, 20 * 60)

# ---------------------------------------------------------------- band pass
lo, hi = BAND
sig = bp.bandpass(B_perp.to_numpy(), lo, hi, fs=1.0, order=4)
env = bp.envelope(sig, smooth_sec=60)
i0, i1 = bp.detect_last_packet(env, 0.30)
tw = B_perp.index
print("detected last packet %s - %s (peak envelope %.3f nT)"
      % (tw[i0].strftime("%H:%M:%S"), tw[i1].strftime("%H:%M:%S"), env.max()))

fig, ax = plt.subplots(figsize=(9, 3.4))
ax.plot(tw, sig, color="#1f77b4", lw=1.3, label="band pass")
ax.axvspan(pd.Timestamp(ACTIVE[0]), pd.Timestamp(ACTIVE[1]),
           color="#fff59d", alpha=0.45, lw=0, zorder=0, label="active period")
ax.set_ylabel(r"$|B_\perp|$ (nT)")
ax.set_xlabel("UTC")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right", fontsize=9)
ax.set_title(r"$B_\perp$ band pass %g-%g mHz, %s   %s-%s UTC"
             % (lo * 1e3, hi * 1e3, EVENT_ID, WIN[0][11:], WIN[1][11:]), fontsize=11)
fig.tight_layout()
p = os.path.join(OUT, "S_quiet_bandpass.png")
fig.savefig(p, dpi=200)
plt.close(fig)
print("saved", p)

# ---------------------------------------------------------------- PSDs
for series, label, fn in ((B_perp, r"$B_\perp$", "S_quiet_psd_Bperp.png"),
                          (B_par, r"$B_{||}$", "S_quiet_psd_Bpar.png")):
    coef, amp, freqs, t_sec, t_utc, f_ion = cwt.compute_cwt(
        series, B_total, **CWT_CFG)
    tf, energy, _ = psd.frequency_energy_spectrum(
        freqs, coef, PSD_FMIN, PSD_FMAX, l1_normalize=False)
    fig = psd.plot_energy_curve(
        tf, energy,
        title="PSD - %s\n%s-%s UTC   %s" % (EVENT_ID, WIN[0][11:], WIN[1][11:], label),
        data_units="nT", l1_normalize=False, figsize=(9, 4.0))
    ax = fig.axes[0]
    ax.set_xlim(PSD_FMIN, PSD_FMAX)
    fig.tight_layout()
    p = os.path.join(OUT, fn)
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("saved", p)
