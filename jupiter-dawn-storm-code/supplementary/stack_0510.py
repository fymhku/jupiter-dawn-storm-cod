# -*- coding: utf-8 -*-
"""Assemble the 2017-05-10 five-panel figure used in the SM.

(a) Magnetic-field components without time shading;
(b,c) CWT spectrograms of B_perp and B_parallel, 1-50 mHz, with f_ci curves;
(d,e) time-averaged wavelet power spectra without percentile-band overlays.
Panel labels are rendered as upright bold text. Numerical settings and
panel assembly are retained from the paper-workspace script.
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

from juno_ulf import components, cwt, data, psd
import textlabel

MAG_DIR = os.environ.get("JUNO_MAG_DIR", os.path.join(ROOT, "data", "mag"))
DATE = "2017_05_10"
EVENT_ID = "2017-05-10"
WIN = ("2017-05-10 08:00:00", "2017-05-10 10:00:00")
CROSSING = ("2017-05-10 08:44:00", "2017-05-10 08:56:00")
ACTIVE = CROSSING
OUT = os.path.join(os.environ.get("JUNO_ULF_OUT", os.path.join(ROOT, "output")), "quiet", "figures")
os.makedirs(OUT, exist_ok=True)

CWT_CFG = dict(wavelet="cmor1.5-1.0", n_scales=200,
               freq_range=(1e-3, 50e-3), ion_mass_amu=16.0)
PSD_FMIN, PSD_FMAX = 1e-3, 50e-3
PERCENTILES = (25, 75)

DPI = 150
FIG_W_IN = 12.0
LABEL_GAP = 10
PANEL_GAP = 22
SECTION_GAP = 40
TARGET_W = 1800                 # resample every panel to this width

fs, ph, ah, bh = textlabel.calibrate(dpi=400)
print("tag: %.2f pt @400dpi -> parenthesis %d px (upright bold)" % (fs, ph))
TAGS = [textlabel.label_bitmap("(%s)" % c, fs, dpi=400) for c in "abcde"]


def rgba_of(fig):
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba())


def resize_w(rgba, width=TARGET_W):
    """bilinear resample to `width`, preserving the aspect ratio"""
    h, w = rgba.shape[:2]
    if w == width:
        return rgba
    nh = max(1, int(round(h * width / w)))
    yi = np.linspace(0, h - 1, nh)
    xi = np.linspace(0, w - 1, width)
    y0 = np.floor(yi).astype(int)
    y1 = np.clip(y0 + 1, 0, h - 1)
    x0 = np.floor(xi).astype(int)
    x1 = np.clip(x0 + 1, 0, w - 1)
    wy = (yi - y0)[:, None, None]
    wx = (xi - x0)[None, :, None]
    a = rgba[y0][:, x0] * (1 - wy) * (1 - wx) + rgba[y1][:, x0] * wy * (1 - wx) \
        + rgba[y0][:, x1] * (1 - wy) * wx + rgba[y1][:, x1] * wy * wx
    return np.clip(a, 0, 255).astype(np.uint8)


# ---------------------------------------------------------------- data
Br, Bth, Bph = data.load_mag_window(DATE, WIN[0], WIN[1], MAG_DIR)
B_total = np.sqrt(Br ** 2 + Bth ** 2 + Bph ** 2)
B_par, B_perp = components.parallel_perp_by_window(Br, Bth, Bph, 20 * 60)

order = []          # (tag index, rgba)

# ---------------------------------------------------------------- (a) field
# clean field panel: no shaded interval, so it carries only the measurements
fig, ax = plt.subplots(figsize=(FIG_W_IN, 4.2), dpi=DPI)
ax.plot(Br.index, Br.values, lw=0.9, label=r"$B_r$")
ax.plot(Bth.index, Bth.values, lw=0.9, label=r"$B_\theta$")
ax.plot(Bph.index, Bph.values, lw=0.9, label=r"$B_\phi$")
ax.set_ylabel("B (nT)")
ax.set_xlabel("UTC")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.grid(True, alpha=0.3)
ax.legend(loc="upper left", fontsize=10, ncol=3, framealpha=0.9)
ax.set_title("%s   %s-%s UTC" % (EVENT_ID, WIN[0][11:], WIN[1][11:]), fontsize=12)
fig.tight_layout()
order.append((0, resize_w(rgba_of(fig))))
plt.close(fig)

# ---------------------------------------------------- (b)(c) CWT spectrograms
psd_panels = []
for idx, (series, label, comp) in enumerate(
        ((B_perp, r"$B_\perp$", "B_perp"),
         (B_par, r"$B_{||}$", "B_parallel"))):
    coef, amp, freqs, t_sec, t_utc, f_ion = cwt.compute_cwt(
        series, B_total, **CWT_CFG)

    fig = cwt.plot_spectrogram(
        amp, freqs, t_sec, t_utc, f_ion, comp,
        title="%s   %s   %s-%s UTC" % (EVENT_ID, label, WIN[0][11:], WIN[1][11:]),
        figsize=(FIG_W_IN, 5.0))
    ax = fig.axes[0]
    # the package already draws the f_ci curve and its legend -- keep both,
    # only widen the frequency axis to 1-50 mHz
    ax.set_ylim(freqs[0], freqs[-1])
    fig.tight_layout()
    order.append((1 + idx, resize_w(rgba_of(fig))))
    plt.close(fig)

    # PSD panel prepared here but appended after both spectrograms
    tf, energy, _ = psd.frequency_energy_spectrum(
        freqs, coef, PSD_FMIN, PSD_FMAX, l1_normalize=False)
    fig = psd.plot_energy_curve(
        tf, energy,
        title="PSD - %s\n%s-%s UTC   %s"
              % (EVENT_ID, WIN[0][11:], WIN[1][11:], label),
        data_units="nT", l1_normalize=False, figsize=(FIG_W_IN, 5.0))
    ax = fig.axes[0]
    ax.set_xlim(PSD_FMIN, PSD_FMAX)
    # no ion-cyclotron band on the PSD plots
    fig.tight_layout()
    psd_panels.append(resize_w(rgba_of(fig)))
    plt.close(fig)

order.append((3, psd_panels[0]))     # (d)
order.append((4, psd_panels[1]))     # (e)

# ---------------------------------------------------------------- stack
W = TARGET_W
H = 0
for n, (ti, r) in enumerate(order):
    H += LABEL_GAP + TAGS[ti].shape[0] + r.shape[0]
    if n < len(order) - 1:
        H += SECTION_GAP if n == 2 else PANEL_GAP

canvas = np.full((H, W, 4), 255, dtype=np.uint8)


def blit(rgba, x0, y0):
    a = rgba[:, :, 3:4].astype(np.float32) / 255.0
    rgb = rgba[:, :, :3].astype(np.float32)
    h, w = a.shape[:2]
    dst = canvas[y0:y0 + h, x0:x0 + w, :3].astype(np.float32)
    canvas[y0:y0 + h, x0:x0 + w, :3] = np.clip(
        rgb * a + dst * (1 - a), 0, 255).astype(np.uint8)


y = 0
for n, (ti, r) in enumerate(order):
    tag = TAGS[ti]
    y += LABEL_GAP
    blit(tag, 8, y)
    y += tag.shape[0]
    canvas[y:y + r.shape[0], :r.shape[1], :] = r
    y += r.shape[0]
    if n < len(order) - 1:
        y += SECTION_GAP if n == 2 else PANEL_GAP

out = os.path.join(OUT, "2017_05_10_field_cwt_stacked.png")
plt.imsave(out, canvas)
print("\npanel order: %s" % " ".join(
    "(%s)" % "abcde"[ti] for ti, _ in order))
print("wrote %s  (%dx%d, %.2f MB)"
      % (out, W, H, os.path.getsize(out) / 1048576))
