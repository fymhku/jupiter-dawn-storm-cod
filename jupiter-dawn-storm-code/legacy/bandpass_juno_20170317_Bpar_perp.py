# -*- coding: utf-8 -*-
"""
Juno MAG (JSS-1s, RTP) ULF analysis: parallel / perpendicular decomposition
then band-pass 6-9 mHz.

Interpretation
--------------
The 2-h window (2017-03-17 08:56:30 -> 10:56:30 UT) is split into six 20-minute
blocks.  For each block the *mean field vector* <B> is taken as the background
(parallel) direction b_hat = <B>/|<B>|.  The field is then decomposed field-
aligned:

    B_||(t)      = B(t) . b_hat                 (parallel, along background)
    B_perp(t)    = B(t) - B_||(t) * b_hat       (leftover / perpendicular field)

Both B_|| and B_perp are band-pass filtered 6-9 mHz.  A transverse hodogram
(measured in an orthonormal basis e1,e2 perpendicular to b_hat) is also produced.

Raw-field figures from the previous run are NOT regenerated here; they are
removed separately.  This script saves only parallel/perpendicular results.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.signal import butter, sosfiltfilt, welch

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
SRC = (r"D:\Juno_Data\Juno_data_txt\mag_tvar\MAG_JSS_1s\juno_B1s_JSSRTP"
       r"\JUNO_MAG_JSS_1s_2017_03_17.txt")
OUT_DIR = r"D:\summer research\pic\bandpass_2017_03_17_Bpar_perp"
os.makedirs(OUT_DIR, exist_ok=True)

FS = 1.0
T_START = "2017-03-17 08:56:30"
T_STOP  = "2017-03-17 10:56:30"
LO, HI = 6e-3, 9e-3
ORDER = 4
BLOCK_MIN = 20                 # length of each background-averaging block

# ----------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------
print("Loading", SRC)
raw = np.genfromtxt(SRC, dtype=None, names=True, encoding="utf-8")
names = raw.dtype.names

def col(key):
    for n in names:
        base = n.replace("(", "").replace(")", "").replace("nT", "").replace("NT", "").strip()
        if base == key:
            return n
    raise KeyError(key)

t = np.array([x.replace("/", " ") for x in raw[col("UT_Time")]]).astype("datetime64[s]")
br = raw[col("Br")].astype(float)
bt = raw[col("Btheta")].astype(float)
bp = raw[col("Bphi")].astype(float)
B = np.column_stack([br, bt, bp])          # (N, 3)

start = np.datetime64(T_START.replace(" ", "T"))
stop  = np.datetime64(T_STOP.replace(" ", "T"))
mask = (t >= start) & (t <= stop)
Bw = B[mask]
tw = t[mask]
N = Bw.shape[0]
print("Points in window:", N)

# Seconds since window start (for block indexing)
t0 = tw[0]
sec = (tw - t0).astype("int64")

# ----------------------------------------------------------------------
# 20-minute block field-aligned decomposition
# ----------------------------------------------------------------------
block_len_s = BLOCK_MIN * 60
block_idx = np.clip(sec // block_len_s, 0, int(np.ceil(sec[-1] / block_len_s)) - 1)
n_blocks = int(block_idx.max()) + 1
print("Number of 20-min blocks:", n_blocks)

Bpar = np.zeros(N)            # parallel component
Bperp = np.zeros((N, 3))      # perpendicular vector
Bperp_mag = np.zeros(N)
z1 = np.zeros(N)              # perpendicular in-plane component 1
z2 = np.zeros(N)              # perpendicular in-plane component 2
bhat_store = np.zeros((n_blocks, 3))

for b in range(n_blocks):
    sel = block_idx == b
    Bmean = Bw[sel].mean(axis=0)
    norm = np.linalg.norm(Bmean)
    if norm == 0:
        continue
    bhat = Bmean / norm
    bhat_store[b] = bhat

    # reference axis least aligned with bhat -> build transverse basis
    axes = np.eye(3)
    dots = np.abs(bhat @ axes)
    ref = axes[int(np.argmin(dots))]
    e1 = np.cross(bhat, ref); e1 /= np.linalg.norm(e1)
    e2 = np.cross(bhat, e1)

    b_par = Bw[sel] @ bhat
    b_perp = Bw[sel] - np.outer(b_par, bhat)
    Bpar[sel] = b_par
    Bperp[sel] = b_perp
    Bperp_mag[sel] = np.linalg.norm(b_perp, axis=1)
    z1[sel] = Bw[sel] @ e1
    z2[sel] = Bw[sel] @ e2

# ----------------------------------------------------------------------
# Band-pass 6-9 mHz (zero-phase Butterworth)
# ----------------------------------------------------------------------
sos = butter(ORDER, [LO, HI], btype="bandpass", fs=FS, output="sos")

def bp(x):
    return sosfiltfilt(sos, x)

Bpar_f = bp(Bpar)
Bperp_mag_f = bp(Bperp_mag)
z1_f = bp(z1)
z2_f = bp(z2)

for nm, sig in [("B_||", Bpar_f), ("|B_perp|", Bperp_mag_f)]:
    print(f"{nm}: band-passed rms {np.std(sig):.4f} nT, "
          f"range [{sig.min():.3f}, {sig.max():.3f}] nT")

# x-axis in matplotlib date units (days since 1970-01-01)
X = tw.astype("int64") / 86400.0

# ----------------------------------------------------------------------
# Figure 1: parallel & perpendicular time series (raw + band-passed)
# ----------------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
ax = axes[0]
ax.plot(X, Bpar, color="0.6", lw=0.9, label="B_|| (raw, along background)")
ax.plot(X, Bpar_f, color="#d62728", lw=1.3, label=f"B_|| band-pass {LO*1e3}-{HI*1e3} mHz")
ax.set_ylabel("B_|| (nT)")
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right", fontsize=8)

ax = axes[1]
ax.plot(X, Bperp_mag, color="0.6", lw=0.9, label="|B_perp| (raw)")
ax.plot(X, Bperp_mag_f, color="#1f77b4", lw=1.3, label=f"|B_perp| band-pass {LO*1e3}-{HI*1e3} mHz")
ax.set_ylabel("|B_perp| (nT)")
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right", fontsize=8)

# draw 20-min block boundaries and label background direction
bound = [t0 + np.timedelta64(b * block_len_s, "s") for b in range(1, n_blocks)]
for bd in bound:
    axes[0].axvline(bd.astype("int64") / 86400.0, color="k", ls=":", lw=0.7, alpha=0.6)
    axes[1].axvline(bd.astype("int64") / 86400.0, color="k", ls=":", lw=0.7, alpha=0.6)
axes[0].set_title("Juno MAG JSS-1s 6-9 mHz band-pass, 08:56:30-10:56:30 UT\n"
                  "per 20-min block background (field-aligned) decomposition",
                  fontsize=10)
axes[-1].set_xlabel("2017-03-17 UT")
axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "Fig01_Bpar_perp_timeseries.png"), dpi=150)
plt.close(fig)

# ----------------------------------------------------------------------
# Figure 2: PSD of band-passed parallel & perpendicular
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 5.5))
nperseg = 2048
for sig, clr, lbl in [(Bpar_f, "#d62728", "B_||"),
                      (Bperp_mag_f, "#1f77b4", "|B_perp|")]:
    f, p = welch(sig, fs=FS, nperseg=nperseg, detrend="linear")
    ax.semilogy(f * 1e3, p, color=clr, lw=1.4, label=f"band-passed {lbl}")
ax.axvspan(LO * 1e3, HI * 1e3, color="orange", alpha=0.2, label="6-9 mHz band")
ax.set_xlabel("Frequency (mHz)")
ax.set_ylabel("PSD (nT²/Hz)")
ax.set_xlim(0, 50)
ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="upper right", fontsize=8)
ax.set_title("Welch PSD of band-passed parallel / perpendicular components")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "Fig02_Bpar_perp_psd.png"), dpi=150)
plt.close(fig)

# ----------------------------------------------------------------------
# Figure 3: transverse hodogram of band-passed perpendicular field
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 7))
sc = ax.scatter(z1_f, z2_f, c=block_idx, cmap="viridis", s=6)
ax.axhline(0, color="0.7", lw=0.8)
ax.axvline(0, color="0.7", lw=0.8)
ax.set_xlabel("perp component 1 (nT, band-passed)")
ax.set_ylabel("perp component 2 (nT, band-passed)")
ax.set_aspect("equal")
cb = fig.colorbar(sc, ax=ax, ticks=range(n_blocks))
cb.set_label("20-min block")
ax.set_title("Transverse hodogram (band-passed perpendicular field)")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "Fig03_hodogram.png"), dpi=150)
plt.close(fig)

print("Saved figures to:", OUT_DIR)
