# -*- coding: utf-8 -*-
"""Regenerate 2017-05-16 perpendicular band-pass figures with the widened
active period (09:17:38-09:46:52 UTC) highlighted, overwriting the
consolidated images."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.signal import butter, sosfiltfilt

EVENT_DATE = '2017_05_16'
T_START = '2017-05-16 07:52:00'
T_STOP  = '2017-05-16 10:03:00'
SRC = (r"D:\Juno_Data\Juno_data_txt\mag_tvar\MAG_JSS_1s\juno_B1s_JSSRTP"
       rf"\JUNO_MAG_JSS_1s_{EVENT_DATE}.txt")
OUT_DIR = r"D:\summer research\pic\bandpass_analysis"
ACT0, ACT1 = '2017-05-16 09:17:38', '2017-05-16 09:46:52'   # widened active period
ORDER = 4
BLOCK_MIN = 20

# (lo, hi, label, outfile)
SERIES = [(6e-3, 9e-3, '6-9', '2017_05_16_perp_bandpass_6-9mHz.png'),
          (2e-3, 5e-3, '2-5', '2017_05_16_perp_bandpass_2-5mHz.png')]

raw = np.genfromtxt(SRC, dtype=None, names=True, encoding="utf-8")
names = raw.dtype.names
def col(key):
    for n in names:
        base = n.replace("(", "").replace(")", "").replace("nT", "").replace("NT", "").strip()
        if base == key:
            return n
    raise KeyError(key)
t = np.array([x.replace("/", " ") for x in raw[col("UT_Time")]]).astype("datetime64[s]")
B = np.column_stack([raw[col("Br")].astype(float),
                     raw[col("Btheta")].astype(float),
                     raw[col("Bphi")].astype(float)])
start = np.datetime64(T_START.replace(" ", "T"))
stop  = np.datetime64(T_STOP.replace(" ", "T"))
mask = (t >= start) & (t <= stop)
Bw, tw = B[mask], t[mask]
N = Bw.shape[0]
sec = (tw - tw[0]).astype("int64")

block_len_s = BLOCK_MIN * 60
block_idx = np.clip(sec // block_len_s, 0, int(np.ceil(sec[-1] / block_len_s)) - 1)
pperp = np.zeros(N)
for b in range(int(block_idx.max()) + 1):
    sel = block_idx == b
    Bmean = Bw[sel].mean(axis=0)
    nrm = np.linalg.norm(Bmean)
    if nrm == 0:
        continue
    bhat = Bmean / nrm
    b_par = Bw[sel] @ bhat
    pperp[sel] = np.sqrt(np.clip(np.sum(Bw[sel] ** 2, axis=1) - b_par ** 2, 0, None))

pa = int(np.argmin(np.abs(tw.astype("int64") - np.datetime64(ACT0.replace(" ", "T")).astype("int64"))))
pb = int(np.argmin(np.abs(tw.astype("int64") - np.datetime64(ACT1.replace(" ", "T")).astype("int64"))))
X = tw.astype("int64") / 86400.0

for lo, hi, lbl, out in SERIES:
    sos = butter(ORDER, [lo, hi], btype="bandpass", fs=1.0, output="sos")
    sig = sosfiltfilt(sos, pperp)
    fig, ax = plt.subplots(figsize=(12, 5.5))
    ax.plot(X, sig, color="#1f77b4", lw=1.5, label="band pass")
    ax.axvspan(X[pa], X[pb], color="#fff59d", alpha=0.45, lw=0, zorder=0,
               label="active period")
    ax.set_ylabel("|B_perp| (nT)")
    ax.set_xlabel("2017-05-16 UTC")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    ax.set_title(rf"$B_\perp$ band pass {lo*1e3:g}-{hi*1e3:g} mHz, "
                 rf"{T_START[11:]}-{T_STOP[11:]} UTC", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, out), dpi=150)
    plt.close(fig)
    print("Saved", os.path.join(OUT_DIR, out))
