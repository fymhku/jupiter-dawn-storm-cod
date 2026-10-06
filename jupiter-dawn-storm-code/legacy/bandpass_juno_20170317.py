# -*- coding: utf-8 -*-
"""
Juno MAG (JSS 1s, RTP frame) band-pass filtering for a ULF/Pc5 interval.

Time window : 2017-03-17 08:56:30 - 10:56:30 UT  (2 h, 1-s cadence)
Band        : 6e-3 - 9e-3 Hz  (period ~ 111-167 s)

Outputs three figures into:
    D:\\summer research\\pic\\bandpass_2017_03_17\\
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless backend
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.signal import butter, sosfiltfilt, welch

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
SRC = (r"D:\Juno_Data\Juno_data_txt\mag_tvar\MAG_JSS_1s\juno_B1s_JSSRTP"
       r"\JUNO_MAG_JSS_1s_2017_03_17.txt")
OUT_DIR = r"D:\summer research\pic\bandpass_2017_03_17"
os.makedirs(OUT_DIR, exist_ok=True)

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
FS = 1.0                       # Hz  (1-s cadence)
T_START = "2017-03-17 08:56:30"
T_STOP  = "2017-03-17 10:56:30"
LO, HI = 6e-3, 9e-3            # Hz
ORDER = 4
COMPONENTS = ["Br", "Btheta", "Bphi"]
COLORS = {"Br": "#d62728", "Btheta": "#1f77b4", "Bphi": "#2ca02c"}

# ----------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------
print("Loading", SRC)
# Column names are on the first line; timestamps look like 2017-03-17/00:00:00
raw = np.genfromtxt(SRC, dtype=None, names=True, encoding="utf-8")
names = raw.dtype.names

# Header names discard punctuation, e.g. "BrnT". Strip trailing unit/number and
# match the base name; also tolerate "(nT)" style just in case.
def col(key):
    for n in names:
        base = n.replace("(", "").replace(")", "").replace("nT", "").replace("NT", "").strip()
        if base == key:
            return n
    raise KeyError(key)

time_name = col("UT_Time")
tstr = np.array([x.replace("/", " ") for x in raw[time_name]])
t = tstr.astype("datetime64[s]")

br = raw[col("Br")].astype(float)
bt = raw[col("Btheta")].astype(float)
bp = raw[col("Bphi")].astype(float)
data = {"Br": br, "Btheta": bt, "Bphi": bp}
print("Headers:", list(names))

# ----------------------------------------------------------------------
# Select time window
# ----------------------------------------------------------------------
start = np.datetime64(T_START.replace(" ", "T").replace(":", ":"))
stop  = np.datetime64(T_STOP.replace(" ", "T").replace(":", ":"))

# Build a clean datetime64 mask
t64 = t
mask = (t64 >= start) & (t64 <= stop)
print("Points in window:", int(mask.sum()))
print("Window span (min):", (t64[mask][-1] - t64[mask][0]) / np.timedelta64(1, "s") / 60.0)

# x-axis in Matplotlib "date" units (days since 1970-01-01) so we can format
# tick labels as HH:MM UT.
X = t64[mask].astype("int64") / 86400.0   # datetime64[s] -> seconds -> days

# ----------------------------------------------------------------------
# Band-pass filter each component (zero-phase Butterworth)
# ----------------------------------------------------------------------
sos = butter(ORDER, [LO, HI], btype="bandpass", fs=FS, output="sos")
filtered = {}
for c in COMPONENTS:
    x = data[c][mask]
    filtered[c] = sosfiltfilt(sos, x)
    print(f"{c}: range raw [{np.nanmin(x):.3f},{np.nanmax(x):.3f}] nT, "
          f"bandpassed rms {np.std(filtered[c]):.4f} nT")

# Optional total field magnitude (context only, unfiltered)
bmag = np.sqrt(br**2 + bt**2 + bp**2)[mask]

# ----------------------------------------------------------------------
# Figure 1: raw + band-passed time series, one panel per component
# ----------------------------------------------------------------------
fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
for ax, c in zip(axes, COMPONENTS):
    x = data[c][mask]
    ax.plot(X, x, color="lightgray", lw=0.8, label="raw")
    ax.plot(X, filtered[c], color=COLORS[c], lw=1.1, label=f"band-pass {LO*1e3}-{HI*1e3} mHz")
    ax.set_ylabel(f"{c} (nT)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
ax = axes[-1]
ax.set_xlabel("Time on 2017-03-17 (UT)")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.set_title("Juno MAG JSS-1s (RTP): band-pass 6-9 mHz, 08:56:30-10:56:30 UT",
             fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "Fig01_bandpass_timeseries.png"), dpi=150)
plt.close(fig)

# ----------------------------------------------------------------------
# Figure 2: overlaid raw vs band-passed (zoomed, first coherent cycle)
# ----------------------------------------------------------------------
fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
for ax, c in zip(axes, COMPONENTS):
    ax.plot(X, data[c][mask], color="gray", lw=0.9, alpha=0.6, label="raw")
    ax.plot(X, filtered[c], color=COLORS[c], lw=1.2, label="band-pass")
    ax.set_ylabel(f"{c} (nT)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
axes[-1].set_xlabel("2017-03-17 UT")
axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
fig.suptitle("Band-passed 6-9 mHz components (full window)", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "Fig02_components_overlay.png"), dpi=150)
plt.close(fig)

# ----------------------------------------------------------------------
# Figure 3: Welch PSD of raw vs band-passed, with the band shaded
# ----------------------------------------------------------------------
fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
for ax, c in zip(axes, COMPONENTS):
    x = data[c][mask]
    f1, p1 = welch(x, fs=FS, nperseg=2048, detrend="linear")
    f2, p2 = welch(filtered[c], fs=FS, nperseg=2048, detrend="linear")
    ax.semilogy(f1 * 1e3, p1, color="gray", lw=1, label="raw PSD")
    ax.semilogy(f2 * 1e3, p2, color=COLORS[c], lw=1.4, label="band-passed PSD")
    ax.axvspan(LO * 1e3, HI * 1e3, color="orange", alpha=0.2, label="6-9 mHz band")
    ax.set_ylabel(f"{c} (nT²/Hz)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
axes[-1].set_xlabel("Frequency (mHz)")
axes[-1].set_xlim(0, 50)
fig.suptitle("Welch PSD / band-pass 6-9 mHz", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "Fig03_psd.png"), dpi=150)
plt.close(fig)

print("Saved figures to:", OUT_DIR)
