# -*- coding: utf-8 -*-
"""
Juno MAG (JSS-1s, RTP) ULF analysis -- PERPENDICULAR component only.

Per the request:
  * keep only the perpendicular (transverse) field,
  * no raw field curves in the images,
  * show how the 6-9 mHz wave evolves (time series + dynamic spectrogram).

Method: per 20-min block, the mean field vector defines the background
direction b_hat; the perpendicular field is the leftover B - (B.b_hat)b_hat.
That transverse vector is expressed in an orthonormal basis (e1,e2) and the
two components are band-pass filtered 6-9 mHz (zero-phase Butterworth).
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.signal import butter, sosfiltfilt, spectrogram

SRC = (r"D:\Juno_Data\Juno_data_txt\mag_tvar\MAG_JSS_1s\juno_B1s_JSSRTP"
       r"\JUNO_MAG_JSS_1s_2017_03_17.txt")
OUT_DIR = r"D:\summer research\pic\bandpass_2017_03_17_perp"
os.makedirs(OUT_DIR, exist_ok=True)

FS = 1.0
T_START = "2017-03-17 08:56:30"
T_STOP  = "2017-03-17 10:56:30"
LO, HI = 6e-3, 9e-3
ORDER = 4
BLOCK_MIN = 20

# --- load ---------------------------------------------------------------
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
Bw = B[mask]
tw = t[mask]
N = Bw.shape[0]
sec = (tw - tw[0]).astype("int64")

# --- 20-min block field-aligned decomposition ---------------------------
block_len_s = BLOCK_MIN * 60
block_idx = np.clip(sec // block_len_s, 0, int(np.ceil(sec[-1] / block_len_s)) - 1)
n_blocks = int(block_idx.max()) + 1

z1 = np.zeros(N)      # transverse comp 1
z2 = np.zeros(N)      # transverse comp 2
for b in range(n_blocks):
    sel = block_idx == b
    Bmean = Bw[sel].mean(axis=0)
    nrm = np.linalg.norm(Bmean)
    if nrm == 0:
        continue
    bhat = Bmean / nrm
    axes = np.eye(3)
    ref = axes[int(np.argmin(np.abs(bhat @ axes)))]
    e1 = np.cross(bhat, ref); e1 /= np.linalg.norm(e1)
    e2 = np.cross(bhat, e1)
    z1[sel] = Bw[sel] @ e1
    z2[sel] = Bw[sel] @ e2

# --- band-pass 6-9 mHz (zero-phase) -------------------------------------
sos = butter(ORDER, [LO, HI], btype="bandpass", fs=FS, output="sos")
z1_f = sosfiltfilt(sos, z1)
z2_f = sosfiltfilt(sos, z2)
amp_f = np.sqrt(z1_f ** 2 + z2_f ** 2)      # transverse wave envelope

print("z1 band-passed rms %.4f nT; z2 rms %.4f nT"
      % (np.std(z1_f), np.std(z2_f)))
print("transverse envelope max %.4f nT" % amp_f.max())

X = tw.astype("int64") / 86400.0

# --- Fig 1: band-passed perpendicular components + envelope -------------
fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
axes[0].plot(X, z1_f, color="#1f77b4", lw=1.2)
axes[0].set_ylabel("perp comp 1 (nT)")
axes[0].grid(True, alpha=0.3)
axes[1].plot(X, z2_f, color="#d62728", lw=1.2)
axes[1].set_ylabel("perp comp 2 (nT)")
axes[1].grid(True, alpha=0.3)
axes[2].plot(X, amp_f, color="#2ca02c", lw=1.3)
axes[2].set_ylabel("transverse amp (nT)")
axes[2].set_xlabel("2017-03-17 UT")
axes[2].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axes[2].grid(True, alpha=0.3)
axes[0].set_title("Perpendicular field, band-pass 6-9 mHz, 08:56:30-10:56:30 UT",
                  fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "Fig01_perp_bandpass_timeseries.png"), dpi=150)
plt.close(fig)

# --- Fig 2: dynamic spectrogram of perpendicular field -------------------
# Use the transverse amplitude evolution across the band.  Spectrogram of the
# perpendicular magnitude (background removed) shows where power sits in
# time-frequency; the 6-9 mHz band is outlined.
perp_mag = np.sqrt(z1 ** 2 + z2 ** 2)     # non-negative perpendicular amplitude
nperseg = 600                             # 10-min window -> ~1.7 mHz resolution
f, tspec, S = spectrogram(perp_mag - perp_mag.mean(), fs=FS,
                          nperseg=nperseg, noverlap=nperseg // 2,
                          detrend="linear", scaling="density")

fig, ax = plt.subplots(figsize=(12, 5.5))
# tspec is in seconds -> convert to matplotlib date units (days since epoch)
Xspec = tw[0].astype("int64") / 86400.0 + tspec / 86400.0
im = ax.pcolormesh(Xspec, f * 1e3, np.log10(S + 1e-30),
                   shading="auto", cmap="viridis")
ax.axhspan(LO * 1e3, HI * 1e3, color="none", ec="w", lw=1.6, ls="--")
ax.text(0.02, 0.93, "6-9 mHz band", transform=ax.transAxes,
        color="w", fontsize=8, va="top")
ax.set_yscale("linear")
ax.set_ylim(0, 30)
ax.set_xlabel("2017-03-17 UT")
ax.set_ylabel("Frequency (mHz)")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.set_title("Dynamic spectrogram of perpendicular field, 6-9 mHz band",
             fontsize=11)
cb = fig.colorbar(im, ax=ax)
cb.set_label("log10 PSD (nT²/Hz)")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "Fig02_perp_spectrogram.png"), dpi=150)
plt.close(fig)

print("Saved figures to:", OUT_DIR)
