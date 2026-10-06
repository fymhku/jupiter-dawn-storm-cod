# -*- coding: utf-8 -*-
"""
Juno MAG (JSS-1s, RTP) -- PERPENDICULAR field band-pass analysis.

Same workflow as the 2017-03-17 event, applied to a given window:
  1) per 20-min block, mean field vector -> background (parallel) direction;
  2) B_perp = B - (B.b_hat) b_hat ; magnitude = sqrt(|B|^2 - (B.b_hat)^2);
  3) band-pass 6-9 mHz (zero-phase) the perpendicular magnitude;
  4) detect the LAST wave packet (envelope) = "active period", highlight yellow;
  5) optional dynamic spectrogram.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.signal import butter, sosfiltfilt, hilbert, spectrogram

# ---------------------- event configuration ------------------------------
EVENT_DATE = '2017_05_16'
T_START = '2017-05-16 07:52:00'
T_STOP  = '2017-05-16 10:03:00'
OUT_DIR = r"D:\summer research\pic\bandpass_2017_05_16_perp_2_5mHz"
SRC = (r"D:\Juno_Data\Juno_data_txt\mag_tvar\MAG_JSS_1s\juno_B1s_JSSRTP"
       rf"\JUNO_MAG_JSS_1s_{EVENT_DATE}.txt")

FS = 1.0
LO, HI = 2e-3, 5e-3
ORDER = 4
BLOCK_MIN = 20
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------- load ----------------------------------------------
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
print(f"Event {EVENT_DATE} {T_START} -> {T_STOP}: {N} points")

# ---------------------- perpendicular magnitude ---------------------------
block_len_s = BLOCK_MIN * 60
block_idx = np.clip(sec // block_len_s, 0, int(np.ceil(sec[-1] / block_len_s)) - 1)
n_blocks = int(block_idx.max()) + 1
pperp = np.zeros(N)
for b in range(n_blocks):
    sel = block_idx == b
    Bmean = Bw[sel].mean(axis=0)
    nrm = np.linalg.norm(Bmean)
    if nrm == 0:
        continue
    bhat = Bmean / nrm
    b_par = Bw[sel] @ bhat
    pperp[sel] = np.sqrt(np.clip(np.sum(Bw[sel] ** 2, axis=1) - b_par ** 2, 0, None))

# ---------------------- band-pass 6-9 mHz --------------------------------
sos = butter(ORDER, [LO, HI], btype="bandpass", fs=FS, output="sos")
pperp_f = sosfiltfilt(sos, pperp)
print("|B_perp| raw range [%.3f, %.3f] nT; band-passed rms %.4f nT"
      % (pperp.min(), pperp.max(), np.std(pperp_f)))

# ---------------------- detect LAST wave packet (active period) ----------
env = np.abs(hilbert(pperp_f))
env_s = np.convolve(env, np.ones(60) / 60.0, mode="same")
thresh = 0.30 * env_s.max()
above = env_s > thresh
runs, i = [], 0
while i < N:
    if above[i]:
        j = i
        while j < N and above[j]:
            j += 1
        runs.append((i, j - 1))
        i = j
    else:
        i += 1
pa, pb = runs[-1]
def hms(ts):
    return str(np.datetime_as_string(ts, unit="s"))[11:19]
t_pkt = (tw[pa], tw[pb])
tstr = f"{hms(t_pkt[0])} - {hms(t_pkt[1])} UTC"
print("LAST wave packet (active period):", tstr, " peak env %.4f nT"
      % env_s[pa:pb + 1].max())

X = tw.astype("int64") / 86400.0

# ---------------------- Fig 1: band-passed perpendicular -----------------
fig, ax = plt.subplots(figsize=(12, 5.5))
ax.plot(X, pperp_f, color="#1f77b4", lw=1.5, label="band pass")
ax.axvspan(t_pkt[0].astype("int64") / 86400.0, t_pkt[1].astype("int64") / 86400.0,
           color="#fff59d", alpha=0.45, lw=0, zorder=0, label="active period")
ax.set_ylabel("|B_perp| (nT)")
ax.set_xlabel("2017-05-16 UTC")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right", fontsize=9)
ax.set_title(rf"$B_\perp$ band pass {LO*1e3:g}-{HI*1e3:g} mHz, "
             rf"{T_START[11:]}-{T_STOP[11:]} UTC", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "Fig01_perp_magnitude.png"), dpi=150)
plt.close(fig)

# ---------------------- Fig 2: dynamic spectrogram -----------------------
nperseg = 600
f, tspec, S = spectrogram(pperp - pperp.mean(), fs=FS, nperseg=nperseg,
                          noverlap=nperseg // 2, detrend="linear",
                          scaling="density")
Xspec = tw[0].astype("int64") / 86400.0 + tspec / 86400.0
fig, ax = plt.subplots(figsize=(12, 5.5))
im = ax.pcolormesh(Xspec, f * 1e3, np.log10(S + 1e-30), shading="auto", cmap="viridis")
ax.axhspan(LO * 1e3, HI * 1e3, color="none", ec="w", lw=1.6, ls="--")
ax.text(0.02, 0.93, f"{LO*1e3:g}-{HI*1e3:g} mHz band", transform=ax.transAxes,
        color="w", fontsize=8, va="top")
ax.set_ylim(0, 30)
ax.set_xlabel("2017-05-16 UTC")
ax.set_ylabel("Frequency (mHz)")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.set_title(f"Dynamic spectrogram of perpendicular field, "
             f"{LO*1e3:g}-{HI*1e3:g} mHz band", fontsize=11)
fig.colorbar(im, ax=ax).set_label("log10 PSD (nT²/Hz)")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "Fig02_perp_spectrogram.png"), dpi=150)
plt.close(fig)

print("Saved figures to:", OUT_DIR)
