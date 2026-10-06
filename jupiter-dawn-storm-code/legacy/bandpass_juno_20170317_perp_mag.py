# -*- coding: utf-8 -*-
"""
Juno MAG (JSS-1s, RTP) -- PERPENDICULAR field as a single scalar.

Per request: do NOT split into in-plane directions.  Everything that is not
the parallel (background) direction counts as perpendicular.  Concretely, for
each 20-min block the mean field vector gives the background direction b_hat;
the perpendicular magnetic field is the leftover

        B_perp(t) = B(t) - (B(t) . b_hat) b_hat

and its magnitude |B_perp| = sqrt(|B|^2 - (B.b_hat)^2) is the single
perpendicular series.  It is band-pass filtered 6-9 mHz (zero-phase) and
plotted over the window.  (The light-gray curve is the raw |B_perp| for
comparison only; the large background field is already removed.)
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.signal import butter, sosfiltfilt, hilbert

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

# --- perpendicular magnitude (everything except parallel) ---------------
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
    b_perp_mag = np.sqrt(np.clip(np.sum(Bw[sel]**2, axis=1) - b_par**2, 0, None))
    pperp[sel] = b_perp_mag

# --- band-pass 6-9 mHz ---------------------------------------------------
sos = butter(ORDER, [LO, HI], btype="bandpass", fs=FS, output="sos")
pperp_f = sosfiltfilt(sos, pperp)
print("|B_perp| raw range [%.3f, %.3f] nT; band-passed rms %.4f nT"
      % (pperp.min(), pperp.max(), np.std(pperp_f)))

# --- detect LAST wave packet from the amplitude envelope ------------------
env = np.abs(hilbert(pperp_f))
env_s = np.convolve(env, np.ones(60) / 60.0, mode="same")   # 1-min smooth
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
pa, pb = runs[-1]                                   # last packet
t_pkt = (tw[pa], tw[pb])
def hms(ts):
    return str(np.datetime_as_string(ts, unit="s"))[11:19]
tstr = f"{hms(t_pkt[0])} - {hms(t_pkt[1])} UTC"
print("LAST wave packet:", tstr, " peak env %.4f nT" % env_s[pa:pb + 1].max())

X = tw.astype("int64") / 86400.0

fig, ax = plt.subplots(figsize=(12, 5.5))
ax.plot(X, pperp_f, color="#1f77b4", lw=1.5, label="band pass")
# highlight the last wave packet in light yellow
ax.axvspan(t_pkt[0].astype("int64") / 86400.0, t_pkt[1].astype("int64") / 86400.0,
           color="#fff59d", alpha=0.45, lw=0, zorder=0,
           label="active period")
ax.set_ylabel("|B_perp| (nT)")
ax.set_xlabel("2017-03-17 UTC")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right", fontsize=9)
ax.set_title(rf"$B_\perp$ band pass 6-9 mHz, {T_START[11:]}-{T_STOP[11:]} UTC",
             fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "Fig01_perp_magnitude.png"), dpi=150)
plt.close(fig)

print("Saved:", OUT_DIR)
