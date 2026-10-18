"""Recompute the supplied S1 Butterworth calculation on the bundled 2017 field window."""
from pathlib import Path
import os
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
SOURCE_ROOT = Path(__file__).resolve().parents[1]
PACKAGE = Path(os.environ.get("JUPITER_RUN_ROOT", str(SOURCE_ROOT/"output"))).resolve()
HERE = PACKAGE/"figures"
HERE.mkdir(parents=True,exist_ok=True)
frame=pd.read_csv(PACKAGE/"data/main_timeseries.csv")
frame.index=pd.DatetimeIndex(pd.to_datetime(frame.UTC)).as_unit("ns")
br,bt,bf=[frame[c] for c in ["Br_nT","Btheta_nT","Bphi_nT"]]
b=np.column_stack([br.to_numpy(),bt.to_numpy(),bf.to_numpy()])
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["DejaVu Sans","Arial"],
"font.size":9,"axes.labelsize":9,"axes.titlesize":9,"xtick.labelsize":9,"ytick.labelsize":9,
"legend.fontsize":9,"axes.linewidth":.65,"lines.linewidth":.9,"pdf.fonttype":42,
"svg.fonttype":"none","axes.spines.top":False,"axes.spines.right":False})
seconds = np.arange(len(br))
block_idx = np.clip(seconds//1200, 0, int(np.ceil(seconds[-1]/1200))-1)
original_perp = np.zeros(len(br))
for block in range(int(block_idx.max())+1):
    select = block_idx == block
    mean = b[select].mean(axis=0)
    bhat = mean / np.linalg.norm(mean)
    par = b[select] @ bhat
    original_perp[select] = np.sqrt(np.clip((b[select]**2).sum(axis=1)-par**2, 0, None))
sos = butter(4, [2e-3,5e-3], btype="bandpass", fs=1., output="sos")
filtered = sosfiltfilt(sos, original_perp)
active = [pd.Timestamp("2017-05-16 09:17:40"),pd.Timestamp("2017-05-16 09:46:44")]
fig, ax = plt.subplots(figsize=(160/25.4,60/25.4))
fig.subplots_adjust(left=0.12, right=.985, bottom=.22, top=.75)
ax.plot(br.index, filtered, color="#0072B2", lw=.9, label="Band-pass signal")
ax.axvspan(*active, color="#fff59d", alpha=.45, lw=0, zorder=0, label="Active period")
ax.set_ylabel("Band-passed B⊥ (nT)")
ax.set_xlabel("Time (UT), 2017 May 16")
ax.set_xlim(br.index[0], br.index[-1])
ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0,15,30,45]))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.grid(True, alpha=.18, lw=.5)
ax.legend(loc="lower right", bbox_to_anchor=(1,1.02), ncol=2,
          frameon=False, handlelength=1.6, borderaxespad=0, columnspacing=1.)
fig.text(.12,.92,"2–5 mHz band-pass; 07:52:00–10:03:00 UT",fontsize=10,ha="left")
fig.savefig(HERE/"FigS1_bandpass_20170516.png",dpi=600)
fig.savefig(HERE/"FigS1_bandpass_20170516.pdf",dpi=600)
fig.savefig(HERE/"FigS1_bandpass_20170516.svg",dpi=600)
plt.close(fig)
pd.DataFrame({"UTC":br.index.astype(str), "Bperp_nT":original_perp,
              "bandpassed_Bperp_nT":filtered}).to_csv(HERE/"FigS1_source_data.csv",index=False)

print("Saved FigS1 and its recomputed numerical series")
